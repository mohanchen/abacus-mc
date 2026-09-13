#ifndef DFTU_NAO_OP_H
#define DFTU_NAO_OP_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

#include <unordered_map>

namespace elecstate
{
template <typename TK, typename TR>
class DensityMatrix;
} // namespace elecstate

namespace hamilt
{

/// The DFTU class template inherits from class T
/// it is used to calculate the non-local pseudopotential of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK, TR> or OperatorPW<TK>
template <class T>
class DFTU : public T
{
};

/// DFTU class template specialization for OperatorLCAO<TK, TR> base class.
/// Adds the DFT+U on-site correction to the real-space Hamiltonian, which is
/// then folded to k-space by the OperatorLCAO machinery:
///   HR(mu,nu;I,J,R) = <phi_{mu,I,0}|chi_m> pot_onsite(m,m') <chi_m'|phi_{nu,J,R}>
///   HK = sum_R e^{ikR} HR
/// where chi_m are the Hubbard projectors of the correlated shell.
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class DFTU<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    DFTU<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                               const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                               hamilt::HContainer<TR>* hR_in,
                               const UnitCell& ucell_in,
                               const Grid_Driver* gridD_in,
                               const TwoCenterIntegrator* intor,
                               const std::vector<double>& orb_cutoff,
                               Plus_U_Base* p_dftu,
                               const int nspin_in,
                               const double onsite_radius,
                               const elecstate::DensityMatrix<TK, double>* dm_in);
    ~DFTU<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() calculates the HR matrix
     * <phi_{\mu, 0}|chi_m> pot_onsite(m,m') <chi_m'|phi_{\nu, R}>
     */
    virtual void contributeHR() override;

    /**
     * @brief get the real-space density matrix of target spin from the solver-owned DensityMatrix
     * @param ispin spin index (0 based): 0 for nspin=1/4, 0/1 for nspin=2
     * @return read-only DMR pointer, or nullptr when DMR has not been calculated yet
     */
    const hamilt::HContainer<double>* get_dmr(int ispin) const;

    /// calculate force and stress for DFT+U
    void cal_force_stress(const bool cal_force,
                          const bool cal_stress,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& stress);

    // Getters for free functions in dftu_nao_fs_r/dftu_nao_for_r/dftu_nao_str_r
    const UnitCell* get_ucell() const { return ucell; }
    Plus_U_Base* get_dftu() const { return dftu; }
    const TwoCenterIntegrator* get_intor() const { return intor_; }
    int get_nspin() const { return nspin; }
    const std::vector<AdjacentAtomInfo>& get_adjs_all() const { return adjs_all; }

    /// transfer pot_onsite format from pauli matrix to normal for non-collinear spin case
    void transfer_pot_onsite(std::vector<double>& pot_onsite_tmp, std::vector<TR>& pot_onsite);

  private:
    const UnitCell* ucell = nullptr;

    Plus_U_Base* dftu = nullptr;

    /// @brief solver-owned density matrix providing DMR; lifetime covers each ionic step
    const elecstate::DensityMatrix<TK, double>* dm_ = nullptr;

    const TwoCenterIntegrator* intor_ = nullptr;

    std::vector<double> orb_cutoff_;

    /// @brief the number of spin components, 1 for no-spin, 2 for collinear spin case and 4 for non-collinear spin case
    int nspin = 0;

    /**
     * @brief build the adjacent-atom lists for all Hubbard atoms and save
     *        them into this->adjs_all. The size of HR will not change in
     *        DFTU, because the DFT+U correction only touches atom pairs
     *        already covered by the Nonlocal operator.
     */
    void initialize_HR(const Grid_Driver* gridD_in, const double onsite_radius);

    /**
     * @brief calculate the <phi|alpha^I> overlap values and save them in this->nlm_tot
     * it will be reused in the calculation of calculate_HR()
     */
    void cal_nlm_all(const Parallel_Orbitals* pv);

    /// @brief occupation matrix of one Hubbard atom (iat0) from the DMR:
    ///        occ(m,m') = sum_R DMR(I,J,R) * <phi_0|chi_m(I)> * <chi_m'(J)|phi_R>
    void compute_occ_from_dmr(int iat0,
                              int target_L,
                              const AdjacentAtomInfo& adjs,
                              const Parallel_Orbitals* pv,
                              std::vector<double>& occ);

    /// @brief load the occupation matrix of one Hubbard atom (iat0) from a
    ///        pre-read occ_mat file
    void load_occ_from_file(int iat0,
                            int target_L,
                            std::vector<double>& occ);

    /// @brief accumulate the HR contributions of one Hubbard atom (iat0)
    ///        from the precomputed pot_onsite
    void accumulate_HR_for_iat0(int iat0,
                                const AdjacentAtomInfo& adjs,
                                const Parallel_Orbitals* pv,
                                const std::vector<TR>& pot_onsite);

    std::vector<AdjacentAtomInfo> adjs_all;
    /// @brief if the nlm_tot is calculated
    bool precal_nlm_done = false;
    /// @brief the overlap values for all [atoms][nerghbors][orb_index(iw) in NAOs][m of target_l in Projectors]
    std::vector<std::vector<std::unordered_map<int, std::vector<double>>>> nlm_tot;
};

} // namespace hamilt
#endif
