#ifndef DFTU_NAO_OP_H
#define DFTU_NAO_OP_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_lcao/module_dftu/dftu_nao_adj.h"

#include <unordered_map>
#include <vector>

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
    DFTU(HS_Matrix_K<TK>* hsk_in,
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
    ~DFTU() = default;

    /**
     * @brief contributeHR() calculates the HR matrix
     * <phi_{\mu, 0}|chi_m> pot_onsite(m,m') <chi_m'|phi_{\nu, R}>
     */
    void contributeHR() override;

    /**
     * @brief get the real-space density matrix of target spin from the solver-owned DensityMatrix
     * @param ispin spin index (0 based): 0 for nspin=1/4, 0/1 for nspin=2
     * @return read-only DMR pointer, or nullptr when DMR has not been calculated yet
     */
    const hamilt::HContainer<double>* get_dmr(int ispin) const;

  private:
    const UnitCell* ucell = nullptr;

    Plus_U_Base* dftu = nullptr;

    /// @brief solver-owned density matrix providing DMR; lifetime covers each ionic step
    const elecstate::DensityMatrix<TK, double>* dm_ = nullptr;

    const TwoCenterIntegrator* intor_ = nullptr;

    std::vector<double> orb_cutoff_;

    /// @brief the number of spin components, 1 for no-spin, 2 for collinear spin case and 4 for non-collinear spin case
    int nspin = 0;

    /// @brief occupation matrix of one Hubbard atom (iat0) from the DMR:
    ///        occ(m,m') = sum_R DMR(I,J,R) * <phi_0|chi_m(I)> * <chi_m'(J)|phi_R>
    void compute_occ_from_dmr(int iat0,
                              int target_L,
                              const AdjacentAtomInfo& adjs,
                              const Parallel_Orbitals* pv,
                              std::vector<double>& occ);

    /// @brief accumulate the HR contributions of one Hubbard atom (iat0)
    ///        from the precomputed pot_onsite
    void accumulate_HR_for_iat0(int iat0,
                                const AdjacentAtomInfo& adjs,
                                const Parallel_Orbitals* pv,
                                const std::vector<TR>& pot_onsite);

    std::vector<AdjacentAtomInfo> adjs_all;
    /// @brief cached <phi|alpha^I> overlap values; empty until first contributeHR() call
    DFTU_LCAO::NlmTot nlm_tot;
};

} // namespace hamilt
#endif
