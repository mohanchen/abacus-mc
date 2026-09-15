#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include "force_lcao.h"
#include "source_base/global_function.h"
#include "source_base/matrix.h"
#include "source_pw/module_pwdft/force_pw.h"
#include "source_pw/module_pwdft/stress_func.h"
#include "source_pw/module_pwdft/stru_fac.h"
#include "source_io/module_parameter/input_conv.h"
#include "source_psi/psi.h"
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri_interface.h"
#endif
#include "force_stress_arrays.h"
#include "source_lcao/setup_exx.h" // for exx, mohan add 20251008
#include "source_lcao/setup_deepks.h" // for deepks, mohan add 20251010
#include "source_lcao/setup_dm.h" // mohan add 2025-11-03
#include "source_pw/module_pwdft/dftu_base.h" // mohan add 2025-11-07
#include "source_hamilt/hamilt.h"
#include "source_hamilt/module_xc/exx_info.h"

namespace vdw
{
struct VdwResult;
}

class TwoCenterBundle;

// Force/stress component matrices assembled by getForceStress. Grouping them
// into a struct lets the assembly/print helpers take one reference instead of
// ~19 individual matrix arguments. Members are default-constructed and only
// created (allocated) when the corresponding term is active.
struct LCAOForceParts
{
    ModuleBase::matrix foverlap;
    ModuleBase::matrix ftvnl_dphi;
    ModuleBase::matrix fvnl_dbeta;
    ModuleBase::matrix fvl_dphi;
    ModuleBase::matrix fvl_dvl;
    ModuleBase::matrix fewalds;
    ModuleBase::matrix fcc;
    ModuleBase::matrix fscc;
    ModuleBase::matrix fvnl_dalpha; // deepks
    ModuleBase::matrix fpothybrid;
    ModuleBase::matrix force_u;
    ModuleBase::matrix force_dspin;
    ModuleBase::matrix force_exx;
    ModuleBase::matrix force_vdw;
    ModuleBase::matrix fefield;
    ModuleBase::matrix fefield_tddft;
    ModuleBase::matrix fgate;
    ModuleBase::matrix fsol;
};

struct LCAOStressParts
{
    ModuleBase::matrix soverlap;
    ModuleBase::matrix stvnl_dphi;
    ModuleBase::matrix svnl_dbeta;
    ModuleBase::matrix svl_dphi;
    ModuleBase::matrix sigmadvl;
    ModuleBase::matrix sigmahar;
    ModuleBase::matrix sigmaewa;
    ModuleBase::matrix sigmacc;
    ModuleBase::matrix sigmaxc;
    ModuleBase::matrix svnl_dalpha; // deepks
    ModuleBase::matrix stress_u;
    ModuleBase::matrix stress_dspin;
    ModuleBase::matrix stress_exx;
    ModuleBase::matrix stress_vdw;
};


template <typename T>
class Force_Stress_LCAO
{
    // mohan add 2021-02-09
    friend class md;
    friend void Input_Conv::Convert();
    friend class ions;

  public:
    Force_Stress_LCAO(Record_adj& ra, const int nat_in);
    ~Force_Stress_LCAO();

    void getForceStress(UnitCell& ucell,
                        const vdw::VdwResult* vdw_result,
                        const bool isforce,
                        const bool isstress,
                        const bool istestf,
                        const bool istests,
                        const Grid_Driver& gd,
                        Parallel_Orbitals& pv,
                        const elecstate::ElecState* pelec,
                        LCAO_domain::Setup_DM<T> &dmat, // mohan add 2025-11-03
                        const psi::Psi<T>* psi,
                        const TwoCenterBundle& two_center_bundle,
                        const LCAO_Orbitals& orb,
                        ModuleBase::matrix& fcs,
                        ModuleBase::matrix& scs,
                        const pseudopot_cell_vl& locpp,
                        const Structure_Factor& sf,
                        const K_Vectors& kv,
                        ModulePW::PW_Basis* rhopw,
						surchem& solvent,
						Plus_U_Base &dftu, // mohan add 2025-11-07
                        Setup_DeePKS<T> &deepks,
                        Exx_NAO<T> &exx_nao,
                        ModuleSymmetry::Symmetry* symm,
                        const Exx_Info& exx_info,
                        const int td_stype = 0,
                        hamilt::Hamilt<T>* p_hamilt = nullptr);

  private:
    int nat;
    Record_adj* RA = nullptr;
    Force_LCAO<T> flk;
    Stress_Func<double> sc_pw;

    // Operator-based force/stress terms: kinetic, overlap, nonlocal,
    // rt-TDDFT hybrid gauge, local-potential Pulay term, and DeltaSpin.
    void cal_operator_fs(const UnitCell& ucell,
                         const Grid_Driver& gd,
                         Parallel_Orbitals& pv,
                         const elecstate::ElecState* pelec,
                         LCAO_domain::Setup_DM<T>& dmat,
                         const psi::Psi<T>* psi,
                         const TwoCenterBundle& two_center_bundle,
                         const LCAO_Orbitals& orb,
                         const K_Vectors& kv,
                         const bool isforce,
                         const bool isstress,
                         const int td_stype,
                         hamilt::Hamilt<T>* p_hamilt,
                         LCAOForceParts& parts,
                         LCAOStressParts& sparts);

    // DeePKS correction force/stress (only active under __MLALGO).
    void cal_deepks_fs(const UnitCell& ucell,
                       const Grid_Driver& gd,
                       const LCAO_Orbitals& orb,
                       const K_Vectors& kv,
                       const bool isforce,
                       const bool isstress,
                       Setup_DeePKS<T>& deepks,
                       LCAOForceParts& parts,
                       LCAOStressParts& sparts);

    // EXX force/stress (only active under __EXX).
    void cal_exx_fs(const UnitCell& ucell,
                    const bool isforce,
                    const bool isstress,
                    const Exx_Info& exx_info,
                    Exx_NAO<T>& exx_nao,
                    LCAOForceParts& parts,
                    LCAOStressParts& sparts);

    // vdW force/stress and external-field forces: E-field, rt-TDDFT E-field,
    // gate field and the implicit solvation model.
    void cal_vdw_and_fields_fs(const vdw::VdwResult* vdw_result,
                               UnitCell& ucell,
                               surchem& solvent,
                               ModulePW::PW_Basis* rhopw,
                               const pseudopot_cell_vl& locpp,
                               const bool isforce,
                               const bool isstress,
                               LCAOForceParts& parts,
                               LCAOStressParts& sparts);

    // DFT+U force/stress.
    void cal_dftu_fs(UnitCell& ucell,
                     const Grid_Driver& gd,
                     Parallel_Orbitals& pv,
                     const LCAO_Orbitals& orb,
                     const K_Vectors& kv,
                     LCAO_domain::Setup_DM<T>& dmat,
                     Plus_U_Base& dftu,
                     const bool isforce,
                     const bool isstress,
                     LCAOForceParts& parts,
                     LCAOStressParts& sparts);

    void forceSymmetry(const UnitCell& ucell, ModuleBase::matrix& fcs, ModuleSymmetry::Symmetry* symm);

    // Sum the computed force parts into fcs, apply symmetry and the net-force
    // (drift) correction, then print the per-term and total forces.
    void assemble_and_print_force(const UnitCell& ucell,
                                  const bool istestf,
                                  const vdw::VdwResult* vdw_result,
                                  const Exx_Info& exx_info,
                                  ModuleSymmetry::Symmetry* symm,
                                  Setup_DeePKS<T>& deepks,
                                  const LCAOForceParts& parts,
                                  ModuleBase::matrix& fcs);

    // Sum the computed stress parts into scs, symmetrize, subtract the external
    // pressure and print the per-term and total stresses.
    void assemble_and_print_stress(const UnitCell& ucell,
                                   const bool istests,
                                   const vdw::VdwResult* vdw_result,
                                   const Exx_Info& exx_info,
                                   ModuleSymmetry::Symmetry* symm,
                                   Setup_DeePKS<T>& deepks,
                                   const LCAOStressParts& parts,
                                   ModuleBase::matrix& scs);

    void calForcePwPart(UnitCell& ucell,
                        ModuleBase::matrix& fvl_dvl,
                        ModuleBase::matrix& fewalds,
                        ModuleBase::matrix& fcc,
                        ModuleBase::matrix& fscc,
                        const double& etxc,
                        const ModuleBase::matrix& vnew,
                        const bool vnew_exist,
                        const Charge* const chr,
                        ModulePW::PW_Basis* rhopw,
                        const pseudopot_cell_vl& locpp,
                        const Structure_Factor& sf);

    void calStressPwPart(UnitCell& ucell,
                         ModuleBase::matrix& sigmadvl,
                         ModuleBase::matrix& sigmahar,
                         ModuleBase::matrix& sigmaewa,
                         ModuleBase::matrix& sigmacc,
                         ModuleBase::matrix& sigmaxc,
                         const double& etxc,
                         const Charge* const chr,
                         ModulePW::PW_Basis* rhopw,
                         const pseudopot_cell_vl& locpp,
                         const Structure_Factor& sf);

    static double force_invalid_threshold_ev;
};

template <typename T>
double Force_Stress_LCAO<T>::force_invalid_threshold_ev = 0.00;

// only for DFT+U, mohan add 2025-11-04
template <typename T>
void assign_dmk_ptr(
    elecstate::DensityMatrix<T,double>* dm,
    std::vector<std::vector<double>>*& dmk_d,
    std::vector<std::vector<std::complex<double>>>*& dmk_c,
    bool gamma_only_local
);

#endif
