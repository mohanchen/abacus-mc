#include "dftu_nao_op.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_base/parallel_reduce.h"

#include "dftu_nao_adj.h"
#include "dftu_nao_fs_r.h"
#include "dftu_nao_ijr.h"
#include "dftu_nao_pots.h"

template <typename TK, typename TR>
hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::DFTU(HS_Matrix_K<TK>* hsk_in,
                                                 const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                                 hamilt::HContainer<TR>* hR_in,
                                                 const UnitCell& ucell_in,
                                                 const Grid_Driver* GridD_in,
                                                 const TwoCenterIntegrator* intor,
                                                 const std::vector<double>& orb_cutoff,
                                                 Plus_U_Base* p_dftu,
                                                 const int nspin_in,
                                                 const double onsite_radius,
                                                 const elecstate::DensityMatrix<TK, double>* dm_in)
    : hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in), intor_(intor), orb_cutoff_(orb_cutoff)
{
    ModuleBase::timer::start("DFTU", "DFTU");
    this->cal_type = calculation_type::lcao_dftu;
    this->ucell = &ucell_in;
    this->dftu = p_dftu;
    this->dm_ = dm_in;

    assert(this->ucell != nullptr);
    assert(this->dm_ != nullptr);

    // build the adjacent-atom lists for all Hubbard atoms. The size of HR
    // will not change in DFTU, because the DFT+U correction only touches
    // atom pairs already covered by the Nonlocal operator.
    this->adjs_all = DFTU_LCAO::build_adjacent_atoms(this->ucell, this->dftu, GridD_in, this->orb_cutoff_, onsite_radius);

    // set nspin
    this->nspin = nspin_in;
    ModuleBase::timer::end("DFTU", "DFTU");
}

// get the read-only real-space density matrix of target spin from the solver-owned DensityMatrix
template <typename TK, typename TR>
const hamilt::HContainer<double>* hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::get_dmr(int ispin) const
{
    assert(ispin >= 0);
    // a not-yet-calculated DMR means the first SCF iteration before the first diagonalization
    if (this->dm_ == nullptr || !this->dm_->is_dmr_ready())
    {
        return nullptr;
    }
    return this->dm_->get_DMR_pointer(ispin + 1);
}

// contributeHR()
/**
 * @brief Contribute DFT+U Hamiltonian to real-space HR matrix
 * 
 * @details This function handles different scenarios based on:
 * 1. Whether occ_mat (occupation matrix) is read from file (is_occmat_ready)
 * 2. Spin configuration (nspin=1, 2, or 4)
 * 3. SCF iteration stage (first vs subsequent iterations)
 * 
 * Case 1: Occ_mat NOT ready (!is_occmat_ready)
 *   - First electronic iteration: calculates occupation matrix from density matrix (DMR)
 *     * Uses get_dmr(current_spin) to get real-space density matrix
 *     * Accumulates contributions from all atom pairs via DFTU_LCAO::cal_occ_ijr()
 *     * Performs MPI reduction to sum occ across processes
 *     * Stores result via set_occ_mat_flat() for use in pot_onsite calculation
 *     * For nspin=1: occ is scaled by 0.5 (since only one spin channel computed)
 *   - Subsequent iterations: occ_mat is computed fresh each iteration from updated DMR
 * 
 * Case 2: Occ_mat IS ready (is_occmat_ready, i.e., read from dm_onsite.txt file)
 *   - First electronic iteration: uses pre-read occ_mat directly without DMR calculation
 *     * Skips DMR-based occ calculation entirely
 *     * Reads occ_mat from stored data via get_occ_mat()
 *     * Different indexing for nspin=4 vs nspin=1/2 (see below)
 *   - After first iteration: set_occmat_stale() is called to force recomputation
 * 
 * Spin configurations:
 *   nspin=1 (non-spin-polarized):
 *     - Single spin channel, occ computed once
 *     - Energy correction doubled at end (set_double_energy)
 *     - current_spin always 0
 *   
 *   nspin=2 (collinear spin-polarized):
 *     - Two separate spin channels (spin-up: 0, spin-down: 1)
 *     - current_spin toggles between 0 and 1 across iterations
 *     - set_occmat_stale() called when current_spin == 1 (last spin)
 *     - HR accumulated separately for each spin
 *   
 *   nspin=4 (non-collinear/SOC):
 *     - Single 4x4 Pauli matrix representation per atom
 *     - occ has 4*(2l+1)^2 elements (spin_fold=4)
 *     - get_occ_mat uses spin=0, ipol indices for Pauli blocks
 *     - set_occmat_stale() always called (current_spin check always true)
 *     - No current_spin toggling (all spins handled simultaneously)
 * 
 * @warning THREAD SAFETY: DFTU_LCAO::cal_hr_ijr() updates shared HR matrix entries.
 *          Different iat0 may contribute to same HR(iat1, iat2, R), requiring
 *          critical section protection for multithreaded correctness.
 *          TODO: Consider refactoring to atom_row_list pattern (see nonlocal.cpp)
 *          for better parallel performance instead of critical section.
 */
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("DFTU", "contributeHR");
    // Early exit: DMR not available AND occ_mat not yet initialized
    const bool dmr_null = (this->get_dmr(0) == nullptr);
    const bool occ_mat_not_init = !this->dftu->is_occmat_ready();

    if (dmr_null && occ_mat_not_init)
    {
        return;
    }
    if (this->current_spin == 0)
    {
        this->dftu->set_energy(0.0);
    }
    ModuleBase::timer::start("DFTU", "contributeHR");

    const Parallel_Orbitals* pv = this->hR->get_atom_pair(0).get_paraV();
    if (this->nlm_tot.empty())
    {
        // lazily precompute the <phi|alpha^I> overlap table; it only depends
        // on the structure and is reused across SCF iterations
        this->nlm_tot = DFTU_LCAO::cal_nlm_all(*this->ucell, *this->dftu, *this->intor_, this->adjs_all, *pv);
    }

    // loop over all Hubbard-projector center atoms (iat0)
    int atom_index = 0;
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        int T0 = 0;
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!this->dftu->has_l_channel(T0))
        {
            continue;
        }
        const int target_L = this->dftu->get_l_channel(T0);
        const int tlp1 = 2 * target_L + 1;
        AdjacentAtomInfo& adjs = this->adjs_all[atom_index++];

        const int spin_fold = (this->nspin == 4) ? 4 : 1;
        std::vector<double> occ(tlp1 * tlp1 * spin_fold, 0.0);

        // compute or load occupation matrix
        if (!this->dftu->is_occmat_ready())
        {
            this->compute_occ_from_dmr(iat0, target_L, adjs, pv, occ);
        }
        else
        {
            DFTU_LCAO::load_occ_from_file(*this->dftu, iat0, target_L, this->nspin, this->current_spin, occ);
        }

        // compute Hubbard potential and energy
        const double u_value = this->dftu->get_u_current(T0);
        std::vector<double> pot_onsite_tmp(occ.size());
        double u_energy = this->dftu->get_energy();
        DFTU_LCAO::cal_pot_onsite(occ, tlp1, u_value, pot_onsite_tmp.data(), u_energy);
        this->dftu->set_energy(u_energy);

        std::vector<TR> pot_onsite(occ.size());
        DFTU_LCAO::transfer_pot_onsite(pot_onsite_tmp, pot_onsite);

        // accumulate HR contributions from neighbor pairs
        this->accumulate_HR_for_iat0(iat0, adjs, pv, pot_onsite);
    }

    // post-processing: energy doubling for nspin=1
    if (this->nspin == 1)
    {
        this->dftu->set_double_energy();
    }
    // mark occ_mat stale for next iteration
    if (this->current_spin == this->nspin - 1 || this->nspin == 4)
    {
        this->dftu->set_occmat_stale();
    }
    // toggle spin channel for nspin=2
    if (this->nspin == 2)
    {
        this->current_spin = 1 - this->current_spin;
    }

    ModuleBase::timer::end("DFTU", "contributeHR");
}

// compute_occ_from_dmr: BRANCH 1 of contributeHR
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::compute_occ_from_dmr(
    int iat0,
    int target_L,
    const AdjacentAtomInfo& adjs,
    const Parallel_Orbitals* pv,
    std::vector<double>& occ)
{
    const hamilt::HContainer<double>* dmR_current = this->get_dmr(this->current_spin);
    for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
    {
        const int T1 = adjs.ntype[ad1];
        const int I1 = adjs.natom[ad1];
        const int iat1 = ucell->itia2iat(T1, I1);
        const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
        const std::unordered_map<int, std::vector<double>>& nlm1 = nlm_tot[iat0][ad1];
        for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
        {
            const int T2 = adjs.ntype[ad2];
            const int I2 = adjs.natom[ad2];
            const int iat2 = ucell->itia2iat(T2, I2);
            const std::unordered_map<int, std::vector<double>>& nlm2 = nlm_tot[iat0][ad2];
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
            ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                              R_index2[1] - R_index1[1],
                                              R_index2[2] - R_index1[2]);
            const hamilt::BaseMatrix<double>* tmp
                = dmR_current->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
            if (tmp != nullptr)
            {
                DFTU_LCAO::cal_occ_ijr(iat1,
                                       iat2,
                                       ucell->get_npol(),
                                       *pv,
                                       nlm1,
                                       nlm2,
                                       tmp->get_pointer(),
                                       occ);
            }
        }
    }
    Parallel_Reduce::reduce_all(occ.data(), occ.size());
    if (this->nspin == 1)
    {
        for (double& v : occ) { v *= 0.5; }
    }
    this->dftu->occmat().set_flat(iat0, target_L, this->current_spin, occ);
}

// accumulate_HR_for_iat0: Step 5 of contributeHR
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::accumulate_HR_for_iat0(
    int iat0,
    const AdjacentAtomInfo& adjs,
    const Parallel_Orbitals* pv,
    const std::vector<TR>& pot_onsite)
{
    for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
    {
        const int T1 = adjs.ntype[ad1];
        const int I1 = adjs.natom[ad1];
        const int iat1 = ucell->itia2iat(T1, I1);
        const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
        const std::unordered_map<int, std::vector<double>>& nlm1 = nlm_tot[iat0][ad1];
        for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
        {
            const int T2 = adjs.ntype[ad2];
            const int I2 = adjs.natom[ad2];
            const int iat2 = ucell->itia2iat(T2, I2);
            const std::unordered_map<int, std::vector<double>>& nlm2 = nlm_tot[iat0][ad2];
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
            ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                              R_index2[1] - R_index1[1],
                                              R_index2[2] - R_index1[2]);
            hamilt::BaseMatrix<TR>* tmp = this->hR->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
            if (tmp != nullptr)
            {
#ifdef _OPENMP
#pragma omp critical(dftu_hr_update)
#endif
                {
                    DFTU_LCAO::cal_hr_ijr<TR>(iat1,
                                              iat2,
                                              ucell->get_npol(),
                                              *pv,
                                              nlm1,
                                              nlm2,
                                              pot_onsite,
                                              tmp->get_pointer());
                }
            }
        }
    }
}

template class hamilt::DFTU<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
