#include "force_stress_lcao.h"

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_base/parallel_reduce.h"
#include "source_pw/module_pwdft/dftu_base.h" //Quxin add for DFT+U on 20201029
#include "source_lcao/module_dftu/dftu_nao_fs_k.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"
// new
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_pot/h_tddft_pw.h"       // Taoni add 2025-02-20
#include "source_estate/module_pot/efield.h"           // liuyu add 2022-05-18
#include "source_estate/module_pot/gatefield.h"        // liuyu add 2022-09-13
#include "source_hamilt/module_surchem/surchem.h" //sunml add 2022-08-10
#include "source_hamilt/module_vdw/vdw.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __MLALGO
#include "source_lcao/module_deepks/lcao_deepks.h"    //caoyu add for deepks 2021-06-03
#include "source_lcao/module_deepks/lcao_deepks_io.h" // mohan add 2024-07-22
#include "source_lcao/module_deepks/deepks_force.h"
#endif
#include "source_lcao/module_dftu/dftu_nao_adj.h"
#include "source_lcao/module_dftu/dftu_nao_fs_r.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/overlap.h"
#include "source_lcao/module_operator_lcao/td_pot_hybrid.h"
#include "source_lcao/pulay_fs.h"
#include "source_lcao/module_rt/force_rt_overlap.h"


// mohan add 2025-11-04
template <>
void assign_dmk_ptr<double>(
    elecstate::DensityMatrix<double,double>* dm,
    std::vector<std::vector<double>>*& dmk_d,
    std::vector<std::vector<std::complex<double>>>*& dmk_c,
    bool gamma_only_local
) {
    auto& dmk_tmp = dm->get_DMK_vector();
    dmk_d = &dmk_tmp;
    dmk_c = nullptr;
}

template <>
void assign_dmk_ptr<std::complex<double>>(
    elecstate::DensityMatrix<std::complex<double>,double>* dm,
    std::vector<std::vector<double>>*& dmk_d,
    std::vector<std::vector<std::complex<double>>>*& dmk_c,
    bool gamma_only_local
) {
    auto& dmk_tmp = dm->get_DMK_vector();
    dmk_c = &dmk_tmp;
    dmk_d = nullptr;
}



template <typename T>
Force_Stress_LCAO<T>::Force_Stress_LCAO(Record_adj& ra, const int nat_in) : RA(&ra), nat(nat_in)
{
}
template <typename T>
Force_Stress_LCAO<T>::~Force_Stress_LCAO()
{
}
template <typename T>
void Force_Stress_LCAO<T>::getForceStress(UnitCell& ucell,
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
                                          Setup_DeePKS<T>& deepks,
                                          Exx_NAO<T> &exx_nao,
                                          ModuleSymmetry::Symmetry* symm,
                                          const Exx_Info& exx_info,
                                          const int td_stype,
                                          hamilt::Hamilt<T>* p_hamilt)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "getForceStress");
    ModuleBase::timer::start("Force_Stress_LCAO", "getForceStress");

    if (!isforce && !isstress)
    {
        ModuleBase::timer::end("Force_Stress_LCAO", "getForceStress");
        return;
    }

    const int nat = ucell.nat;

    // NOTE: ForceStressArrays is no longer needed as we use operator-based force calculation
    // ForceStressArrays fsr; // removed - no longer needed

    // total force : ModuleBase::matrix fcs;

    // part of total force / stress, grouped so the assembly helpers can take a
    // single container reference.
    LCAOForceParts parts;
    LCAOStressParts sparts;

    parts.fvl_dphi.create(nat, 3); // must do it now, update it later, noted by zhengdy

    if (isforce)
    {
        fcs.create(nat, 3);
        parts.foverlap.create(nat, 3); // overlap force
        parts.ftvnl_dphi.create(nat, 3); // pulay force of NAO
        parts.fvnl_dbeta.create(nat, 3); // pulay force of non-local projectors
        parts.fvl_dvl.create(nat, 3); // force from local potentials
        parts.fewalds.create(nat, 3); // Ewald force
        parts.fcc.create(nat, 3); // force due to core correction
        parts.fscc.create(nat, 3); // force due to self-consistent field
        parts.fvnl_dalpha.create(nat, 3); // deepks
        parts.fpothybrid.create(nat, 3); // pulay force for hybrid gauge rt-tddft

        // calculate basic terms in Force, same method with PW base
        this->calForcePwPart(ucell, parts.fvl_dvl, parts.fewalds, parts.fcc, parts.fscc, pelec->f_en.etxc,
              pelec->vnew, pelec->vnew_exist, pelec->charge, rhopw, locpp, sf);
    }

    // total stress : ModuleBase::matrix scs

    //! stress
    if (isstress)
    {
        scs.create(3, 3);
        sparts.sigmacc.create(3, 3);
        sparts.sigmadvl.create(3, 3);
        sparts.sigmaewa.create(3, 3);
        sparts.sigmaxc.create(3, 3);
        sparts.sigmahar.create(3, 3);

        sparts.soverlap.create(3, 3);
        sparts.stvnl_dphi.create(3, 3);
        sparts.svnl_dbeta.create(3, 3);
        sparts.svl_dphi.create(3, 3);
        sparts.svnl_dalpha.create(3, 3);

        // calculate basic terms in Stress, similar method with PW base
        this->calStressPwPart(ucell, sparts.sigmadvl, sparts.sigmahar, sparts.sigmaewa, sparts.sigmacc,
          sparts.sigmaxc, pelec->f_en.etxc, pelec->charge, rhopw, locpp, sf);
    }
    // Calculate operator-based force/stress terms (kinetic, overlap,
    // nonlocal, rt-TDDFT hybrid gauge, local Pulay term and DeltaSpin).
    this->cal_operator_fs(ucell, gd, pv, pelec, dmat, psi, two_center_bundle,
                          orb, kv, isforce, isstress, td_stype, p_hamilt, parts, sparts);

    // MPI reduction for forces
    if (isforce)
    {
        Parallel_Reduce::reduce_pool(parts.fvl_dphi.c, parts.fvl_dphi.nr * parts.fvl_dphi.nc);
    }

    // MPI reduction for stresses
    if (isstress)
    {
        Parallel_Reduce::reduce_pool(sparts.svl_dphi.c, sparts.svl_dphi.nr * sparts.svl_dphi.nc);
    }

    // Handle DeePKS forces if enabled
    this->cal_deepks_fs(ucell, gd, orb, kv, isforce, isstress, deepks, parts, sparts);

    // vdW force/stress and external-field forces
    this->cal_vdw_and_fields_fs(vdw_result, ucell, solvent, rhopw, locpp,
                                isforce, isstress, parts, sparts);

    // DFT+U force/stress
    this->cal_dftu_fs(ucell, gd, pv, orb, kv, dmat, dftu, isforce, isstress, parts, sparts);


    // NOTE: finish_ftable is no longer needed as we don't use ForceStressArrays for overlap/kinetic
    // if (!PARAM.globalv.gamma_only_local)
    // {
    //     this->flk.finish_ftable(fsr);
    // }

    // EXX force/stress
    this->cal_exx_fs(ucell, isforce, isstress, exx_info, exx_nao, parts, sparts);
    //--------------------------------
    // begin calculate and output force
    //--------------------------------
    if (isforce)
    {
        this->assemble_and_print_force(ucell, istestf, vdw_result, exx_info, symm, deepks, parts, fcs);
    } // end of force calculation
    //---------------------------------
    // begin calculate and output stress
    //---------------------------------
    if (isstress)
    {
        this->assemble_and_print_stress(ucell, istests, vdw_result, exx_info, symm, deepks, sparts, scs);
    } // end of stress calculation

    ModuleBase::timer::end("Force_Stress_LCAO", "getForceStress");
    return;
}

// Operator-based force/stress terms: kinetic, overlap, nonlocal,
// rt-TDDFT hybrid gauge, local-potential Pulay term, and DeltaSpin.
template <typename T>
void Force_Stress_LCAO<T>::cal_operator_fs(const UnitCell& ucell,
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
                                             LCAOStressParts& sparts)
{

    // Calculate forces and stresses using new operator-based methods
    // Step 1: Calculate Energy Density Matrix (EDM) for overlap force
    // EDM = Σ_k w_k * ε_k * |ψ_k><ψ_k|
    elecstate::DensityMatrix<T, double> edm = flk.cal_edm(pelec, *psi, *dmat.dm, kv, pv,
                                                           PARAM.inp.nspin, PARAM.inp.nbands, ucell, *this->RA);

    // Step 2: Handle different spin cases
    if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2)
    {
        // For nspin=1 or nspin=2, use double precision
        // Switch to spin channel 1 for DMR access
        if (PARAM.inp.nspin == 2)
        {
            dmat.dm->switch_dmr(1);
            edm.switch_dmr(1);
        }

        const hamilt::HContainer<double>* dmR = dmat.dm->get_DMR_pointer(1);
        const hamilt::HContainer<double>* edmR = edm.get_DMR_pointer(1);

        // Calculate kinetic force/stress (uses DM)
        if (PARAM.inp.t_in_h)
        {
            hamilt::EKinetic<hamilt::OperatorLCAO<T, double>> tmp_ekinetic(
                nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
                two_center_bundle.kinetic_orb.get());
            tmp_ekinetic.cal_force_stress(isforce, isstress, dmR, parts.ftvnl_dphi, sparts.stvnl_dphi);
        }

        // Calculate overlap force/stress (uses EDM)
        hamilt::Overlap<hamilt::OperatorLCAO<T, double>> tmp_overlap(
            nullptr, kv.kvec_d, nullptr, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb.get());
        if(td_stype != 2)
        {
            tmp_overlap.cal_force_stress(isforce, isstress, edmR, parts.foverlap, sparts.soverlap);
        }

        // Calculate nonlocal force/stress (uses DM)
        hamilt::Nonlocal<hamilt::OperatorLCAO<T, double>> tmp_nonlocal(
            nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb_beta.get());
        tmp_nonlocal.cal_force_stress(isforce, isstress, dmR, parts.fvnl_dbeta, sparts.svnl_dbeta);
        
        if(td_stype == 2)
        {
            hamilt::TD_pot_hybrid<hamilt::OperatorLCAO<T, double>> tmp_hybrid(
                nullptr, &kv, nullptr, nullptr, orb, &ucell, orb.cutoffs(), &gd, nullptr);
            tmp_hybrid.cal_force_stress(isforce, dmR, parts.fpothybrid);

            cal_foverlap_rt(parts.foverlap, dmat, p_hamilt, kv, pv, ucell);
        }

        // Switch back to spin channel 0
        if (PARAM.inp.nspin == 2)
        {
            dmat.dm->switch_dmr(0);
            edm.switch_dmr(0);
        }

        // Calculate local potential force/stress (vl_dphi)
        // This uses grid integration, not operator-based method
        flk.ParaV = dmat.dm->get_paraV_pointer();
        PulayForceStress::cal_pulay_fs(parts.fvl_dphi, sparts.svl_dphi, *dmat.dm, ucell, pelec->pot,
                                       isforce, isstress, false /*reset dm to gint*/);
    }
    else if (PARAM.inp.nspin == 4)
    {

        // Calculate kinetic force/stress (uses DM)
        if (PARAM.inp.t_in_h)
        {
            hamilt::EKinetic<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_ekinetic(
                nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
                two_center_bundle.kinetic_orb.get());
            tmp_ekinetic.cal_force_stress(isforce, isstress, dmat.dm->get_DMR_pointer(1), parts.ftvnl_dphi, sparts.stvnl_dphi);
        }

        // Calculate overlap force/stress (uses EDM)
        hamilt::Overlap<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_overlap(
            nullptr, kv.kvec_d, nullptr, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb.get());
        tmp_overlap.cal_force_stress(isforce, isstress, edm.get_DMR_pointer(1), parts.foverlap, sparts.soverlap);

        // For nspin=4 (non-collinear), need complex DMR
        // Create temporary complex DMR for DM
        hamilt::HContainer<std::complex<double>> tmp_dmr(dmat.dm->get_DMR_pointer(1)->get_paraV());
        std::vector<int> ijrs = dmat.dm->get_DMR_pointer(1)->get_ijr_info();
        tmp_dmr.insert_ijrs(&ijrs);
        tmp_dmr.allocate();
        dmat.dm->cal_DMR_full(&tmp_dmr);
        // Calculate nonlocal force/stress (uses DM)
        hamilt::Nonlocal<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>> tmp_nonlocal(
            nullptr, kv.kvec_d, nullptr, &ucell, orb.cutoffs(), &gd,
            two_center_bundle.overlap_orb_beta.get());
        tmp_nonlocal.cal_force_stress(isforce, isstress, &tmp_dmr, parts.fvnl_dbeta, sparts.svnl_dbeta);

        // Calculate local potential force/stress (vl_dphi)
        flk.ParaV = dmat.dm->get_paraV_pointer();
        PulayForceStress::cal_pulay_fs(parts.fvl_dphi, sparts.svl_dphi, *dmat.dm, ucell, pelec->pot,
                                       isforce, isstress, false /*reset dm to gint*/);
    }

    // atomic force and stress for DeltaSpin
    if (PARAM.inp.sc_mag_switch)
    {
        if (isforce)
        {
            parts.force_dspin.create(ucell.nat, 3);
        }
        if (isstress)
        {
            sparts.stress_dspin.create(3, 3);
        }

        hamilt::DeltaSpin<hamilt::OperatorLCAO<T, double>> tmp_dspin(nullptr,
                                                                     kv.kvec_d,
                                                                     nullptr,
                                                                     ucell,
                                                                     &gd,
                                                                     two_center_bundle.overlap_orb_onsite.get(),
                                                                     orb.cutoffs());

        if (PARAM.inp.nspin == 2)
        {
            dmat.dm->switch_dmr(2);
        }
        const hamilt::HContainer<double>* dmr = dmat.dm->get_DMR_pointer(1);
        tmp_dspin.cal_force_stress(isforce, isstress, dmr, parts.force_dspin, sparts.stress_dspin);
        if (PARAM.inp.nspin == 2)
        {
            dmat.dm->switch_dmr(0);
        }
    }
}

// DeePKS correction force/stress (only active under __MLALGO).
template <typename T>
void Force_Stress_LCAO<T>::cal_deepks_fs(const UnitCell& ucell,
                                           const Grid_Driver& gd,
                                           const LCAO_Orbitals& orb,
                                           const K_Vectors& kv,
                                           const bool isforce,
                                           const bool isstress,
                                           Setup_DeePKS<T>& deepks,
                                           LCAOForceParts& parts,
                                           LCAOStressParts& sparts)
{
    // Handle DeePKS forces if enabled
#ifdef __MLALGO
    if (PARAM.inp.deepks_scf)
    {
        const int nks = (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2) ? 1 : kv.get_nks();
        if (PARAM.globalv.gamma_only_local)
        {
            DeePKS_domain::cal_f_delta<double>(
                ucell,
                orb,
                gd,
                *flk.ParaV,
                nks,
                deepks.ld.deepks_param,
                kv.kvec_d,
                deepks.ld.phialpha,
                parts.fvnl_dalpha,
                isstress,
                sparts.svnl_dalpha,
                deepks.ld.dm_r,
                deepks.ld.gedm,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.dm_r_mag : nullptr,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.gedm_mag : nullptr);
        }
        else
        {
            DeePKS_domain::cal_f_delta<std::complex<double>>(
                ucell,
                orb,
                gd,
                *flk.ParaV,
                nks,
                deepks.ld.deepks_param,
                kv.kvec_d,
                deepks.ld.phialpha,
                parts.fvnl_dalpha,
                isstress,
                sparts.svnl_dalpha,
                deepks.ld.dm_r,
                deepks.ld.gedm,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.dm_r_mag : nullptr,
                (PARAM.inp.nspin == 2 && !PARAM.inp.deepks_equiv) ? deepks.ld.gedm_mag : nullptr);
        }

        if (isforce)
        {
            Parallel_Reduce::reduce_pool(parts.fvnl_dalpha.c, parts.fvnl_dalpha.nr * parts.fvnl_dalpha.nc);
        }
        if (isstress)
        {
            Parallel_Reduce::reduce_pool(sparts.svnl_dalpha.c, sparts.svnl_dalpha.nr * sparts.svnl_dalpha.nc);
        }
    }
#endif
}

// EXX force/stress (only active under __EXX).
template <typename T>
void Force_Stress_LCAO<T>::cal_exx_fs(const UnitCell& ucell,
                                        const bool isforce,
                                        const bool isstress,
                                        const Exx_Info& exx_info,
                                        Exx_NAO<T>& exx_nao,
                                        LCAOForceParts& parts,
                                        LCAOStressParts& sparts)
{
#ifdef __EXX
    bool cal_exx = exx_info.info_global.cal_exx;
    bool real_number = exx_info.info_ri.real_number;
    double hybrid_alpha = exx_info.info_global.hybrid_alpha;

    if (cal_exx)
    {
        if (isforce)
        {
            if (real_number)
            {
                exx_nao.exd->cal_exx_force(ucell.nat);
                parts.force_exx = hybrid_alpha * exx_nao.exd->get_force();
            }
            else
            {
                exx_nao.exc->cal_exx_force(ucell.nat);
                parts.force_exx = hybrid_alpha * exx_nao.exc->get_force();
            }
        }
        if (isstress)
        {
            if (real_number)
            {
                exx_nao.exd->cal_exx_stress(ucell.omega, ucell.lat0);
                sparts.stress_exx = hybrid_alpha * exx_nao.exd->get_stress();
            }
            else
            {
                exx_nao.exc->cal_exx_stress(ucell.omega, ucell.lat0);
                sparts.stress_exx = hybrid_alpha * exx_nao.exc->get_stress();
            }
        }
    }
#endif
}

// vdW force/stress and external-field forces: E-field, rt-TDDFT E-field,
// gate field and the implicit solvation model.
template <typename T>
void Force_Stress_LCAO<T>::cal_vdw_and_fields_fs(const vdw::VdwResult* vdw_result,
                                                   UnitCell& ucell,
                                                   surchem& solvent,
                                                   ModulePW::PW_Basis* rhopw,
                                                   const pseudopot_cell_vl& locpp,
                                                   const bool isforce,
                                                   const bool isstress,
                                                   LCAOForceParts& parts,
                                                   LCAOStressParts& sparts)
{
    //! forces and stress from vdw
    //  Peize Lin add 2014-04-04, update 2021-03-09
    //  jiyy add 2019-05-18, update 2021-05-02
    if (vdw_result != nullptr)
    {
        if (isforce)
        {
            if (!vdw_result->has_force || vdw_result->force.size() != static_cast<std::size_t>(ucell.nat))
            {
                ModuleBase::WARNING_QUIT("Force_Stress_LCAO::getForceStress",
                                         "The cached vdW force is unavailable or has an invalid size.");
            }
            parts.force_vdw.create(ucell.nat, 3);
            for (int iat = 0; iat < nat; ++iat)
            {
                parts.force_vdw(iat, 0) = vdw_result->force[iat].x;
                parts.force_vdw(iat, 1) = vdw_result->force[iat].y;
                parts.force_vdw(iat, 2) = vdw_result->force[iat].z;
            }
        }
        if (isstress)
        {
            if (!vdw_result->has_stress)
            {
                ModuleBase::WARNING_QUIT("Force_Stress_LCAO::getForceStress",
                                         "The cached vdW stress is unavailable.");
            }
            sparts.stress_vdw = vdw_result->stress.to_matrix();
        }
    }

    //! forces from E-field
    if (PARAM.inp.efield_flag && isforce)
    {
        parts.fefield.create(ucell.nat, 3);
        elecstate::Efield::compute_force(ucell, parts.fefield);
    }

    //! atomic forces from E-field of rt-TDDFT
    if (PARAM.inp.esolver_type == "tddft" && isforce)
    {
        parts.fefield_tddft.create(ucell.nat, 3);
        elecstate::H_TDDFT_pw::compute_force(ucell, parts.fefield_tddft);
    }

    //! atomic forces from gate field
    if (PARAM.inp.gate_flag && isforce)
    {
        parts.fgate.create(ucell.nat, 3);
        elecstate::Gatefield::compute_force(ucell, parts.fgate);
    }

    //! atomic forces from implicit solvation model
    if (PARAM.inp.imp_sol && isforce)
    {
        parts.fsol.create(ucell.nat, 3);
        solvent.cal_force_sol(ucell, rhopw, locpp.vloc, PARAM.inp.nspin, parts.fsol);
    }
}

// DFT+U force/stress.
template <typename T>
void Force_Stress_LCAO<T>::cal_dftu_fs(UnitCell& ucell,
                                         const Grid_Driver& gd,
                                         Parallel_Orbitals& pv,
                                         const LCAO_Orbitals& orb,
                                         const K_Vectors& kv,
                                         LCAO_domain::Setup_DM<T>& dmat,
                                         Plus_U_Base& dftu,
                                         const bool isforce,
                                         const bool isstress,
                                         LCAOForceParts& parts,
                                         LCAOStressParts& sparts)
{
    //! atomic forces from DFT+U (Quxin version)

    if (PARAM.inp.dft_plus_u) // Quxin add for DFT+U on 20201029
    {
        if (isforce)
        {
            parts.force_u.create(ucell.nat, 3);
        }
        if (isstress)
        {
            sparts.stress_u.create(3, 3);
        }
        if (PARAM.inp.dft_plus_u == 2)
        {
            // The legacy dft_plus_u==2 force/stress path is currently broken.
            //
            // Background: Plus_U::force_stress relies on ForceStressArrays
            // members DSloc_x/y/z (gamma_only) or DSloc_Rx/Ry/Rz (multik)
            // and DH_r being pre-allocated and filled with dS/dR data by the
            // main force flow (formerly ForceLcaoGamma::ftable). The DFT+U
            // step 2 refactor (commit 70c54c9d5a, 2026-01-23) removed the
            // main-flow ForceStressArrays because the operator-based force
            // calculation no longer needs it, but the legacy dft_plus_u==2
            // path still depends on it. The local fsr_dftu below is declared
            // without allocating those arrays, so any call into
            // cal_force_gamma / cal_stress_gamma / folding_matrix_k would
            // pass nullptr to pdgemm_ and crash with SIGSEGV.
            //
            // Until the legacy path is restored or re-implemented, we
            // explicitly reject dft_plus_u==2 with cal_force or cal_stress
            // enabled. SCF-only runs (no force/stress) are unaffected
            // because the energy is computed in cal_energy_correction,
            // which does not touch DSloc arrays. Use dft_plus_u=1 for
            // force/stress calculations.
            if (isforce || isstress)
            {
                ModuleBase::WARNING_QUIT("Force_Stress_LCAO::getForceStress",
                    "dft_plus_u==2 with cal_force or cal_stress is currently broken; "
                    "please use dft_plus_u=1 instead. See notes in source/source_lcao/force_stress_lcao.cpp.");
            }
            ForceStressArrays fsr_dftu;
            std::vector<std::vector<double>>* dmk_d = nullptr;
            std::vector<std::vector<std::complex<double>>>* dmk_c = nullptr;
            assign_dmk_ptr<T>(dmat.dm, dmk_d, dmk_c, PARAM.globalv.gamma_only_local);
            DFTU_LCAO::DftuFsEnv dftu_fs_env(dftu, ucell, gd, pv, fsr_dftu,
                                             orb.cutoffs(), PARAM.inp.ks_solver);
            DFTU_LCAO::force_stress(dftu_fs_env, isforce, isstress,
                                    dmk_d, dmk_c, parts.force_u, sparts.stress_u, kv,
                                    PARAM.globalv.gamma_only_local);
        }
        else
        {
            // Build DFT+U force/stress inputs directly without constructing a
            // full DFTU operator (hsk/hR are irrelevant for this path).
            auto adjs_all = DFTU_LCAO::build_adjacent_atoms(
                &ucell, &dftu, &gd, orb.cutoffs(), PARAM.inp.onsite_radius);

            // The DensityMatrix holds nspin_dm = (nspin==2 ? 2 : 1) real-space DMR
            // channels: nspin=4 (non-collinear) packs all four Pauli components
            // into a single complex DMR, so only one channel exists (cf. setup_dm.cpp
            // and the is0 = nspin==2 ? is : 0 indexing in cal_for/str_IJR_nao_r).
            const int nspin_dm = (PARAM.inp.nspin == 2) ? 2 : 1;
            std::vector<const hamilt::HContainer<double>*> dmR_tmp(nspin_dm, nullptr);
            for (int is = 0; is < nspin_dm; ++is)
            {
                dmR_tmp[is] = dmat.dm->get_DMR_pointer(is + 1);
            }

            DFTU_LCAO::cal_fs_nao_r(&ucell, &dftu,
                                    two_center_bundle.overlap_orb_onsite.get(),
                                    PARAM.inp.nspin,
                                    adjs_all, dmR_tmp,
                                    isforce, isstress, parts.force_u, sparts.stress_u);
        }
    }
}

// local pseudopotential, ewald, core correction, scc terms in force
template <typename T>
void Force_Stress_LCAO<T>::calForcePwPart(UnitCell& ucell,
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
                                          const Structure_Factor& sf)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "calForcePwPart");
#ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Forces<double, base_device::DEVICE_GPU> f_pw(nat);
        f_pw.cal_force_loc(ucell, fvl_dvl, rhopw, locpp.vloc, chr);
        f_pw.cal_force_ew(ucell, fewalds, rhopw, &sf);
        f_pw.cal_force_cc(fcc, rhopw, chr, locpp.numeric, ucell);
        f_pw.cal_force_scc(fscc, rhopw, vnew, vnew_exist, locpp.numeric, ucell);
    }
    else
#endif
    {
        Forces<double, base_device::DEVICE_CPU> f_pw(nat);
        f_pw.cal_force_loc(ucell, fvl_dvl, rhopw, locpp.vloc, chr);
        f_pw.cal_force_ew(ucell, fewalds, rhopw, &sf);
        f_pw.cal_force_cc(fcc, rhopw, chr, locpp.numeric, ucell);
        f_pw.cal_force_scc(fscc, rhopw, vnew, vnew_exist, locpp.numeric, ucell);
    }

    return;
}

// vlocal, hartree, ewald, core correction, exchange-correlation terms in stress
template <typename T>
void Force_Stress_LCAO<T>::calStressPwPart(UnitCell& ucell,
                                           ModuleBase::matrix& sigmadvl,
                                           ModuleBase::matrix& sigmahar,
                                           ModuleBase::matrix& sigmaewa,
                                           ModuleBase::matrix& sigmacc,
                                           ModuleBase::matrix& sigmaxc,
                                           const double& etxc,
                                           const Charge* const chr,
                                           ModulePW::PW_Basis* rhopw,
                                           const pseudopot_cell_vl& locpp,
                                           const Structure_Factor& sf)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "calStressPwPart");

    // local pseudopotential stress:
    sc_pw.stress_loc(ucell, sigmadvl, rhopw, locpp.vloc, &sf, 0, chr);

    // hartree term
    sc_pw.stress_har(ucell, sigmahar, rhopw, 0, chr);

    // ewald stress: use plane wave only.
    sc_pw.stress_ewa(ucell, sigmaewa, rhopw, 0); // remain problem

    // stress due to core correlation.
    sc_pw.stress_cc(sigmacc, rhopw, ucell, &sf, 0, locpp.numeric, chr);

    // stress due to self-consistent charge.
    for (int i = 0; i < 3; i++)
    {
        sigmaxc(i, i) = -etxc / ucell.omega;
    }
    // Exchange-correlation for PBE
    sc_pw.stress_gga(ucell, sigmaxc, rhopw, chr);

    return;
}

#include "source_base/mathzone.h"
template <typename T>
void Force_Stress_LCAO<T>::assemble_and_print_force(const UnitCell& ucell,
                                                   const bool istestf,
                                                   const vdw::VdwResult* vdw_result,
                                                   const Exx_Info& exx_info,
                                                   ModuleSymmetry::Symmetry* symm,
                                                   Setup_DeePKS<T>& deepks,
                                                   const LCAOForceParts& parts,
                                                   ModuleBase::matrix& fcs)
{
    const int nat = ucell.nat;
    //---------------------------------
    // sum all parts of force!
    //---------------------------------
    ModuleBase::Vector3<double> net_force = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++)
    {
        for (int iat = 0; iat < nat; iat++)
        {
            fcs(iat, i) += parts.foverlap(iat, i) + parts.ftvnl_dphi(iat, i) + parts.fvnl_dbeta(iat, i) + parts.fvl_dphi(iat, i)
                           + parts.fvl_dvl(iat, i) // derivative of local potential force (pw)
                           + parts.fewalds(iat, i) // ewald force (pw)
                           + parts.fcc(iat, i)     // nonlinear core correction force (pw)
                           + parts.fscc(iat, i)    // self consistent corretion force (pw)
                           + parts.fpothybrid(iat, i); // pulay force for hybrid gauge rt-tddft

            // Force contribution from DFT+U, Quxin add on 20201029
            if (PARAM.inp.dft_plus_u)
            {
                fcs(iat, i) += parts.force_u(iat, i);
            }
            if (PARAM.inp.sc_mag_switch)
            {
                fcs(iat, i) += parts.force_dspin(iat, i);
            }
#ifdef __EXX
            // Force contribution from exx
            if (exx_info.info_global.cal_exx)
            {
                fcs(iat, i) += parts.force_exx(iat, i);
            }
#endif
            // VDW force of vdwd2 or vdwd3
            if (vdw_result != nullptr)
            {
                fcs(iat, i) += parts.force_vdw(iat, i);
            }
            // E-field force
            if (PARAM.inp.efield_flag)
            {
                fcs(iat, i) += parts.fefield(iat, i);
            }
            // E-field force of tddft
            if (PARAM.inp.esolver_type == "tddft")
            {
                fcs(iat, i) += parts.fefield_tddft(iat, i);
            }
            // Gate field force
            if (PARAM.inp.gate_flag)
            {
                fcs(iat, i) += parts.fgate(iat, i);
            }
            // implicit solvation model
            if (PARAM.inp.imp_sol)
            {
                fcs(iat, i) += parts.fsol(iat, i);
            }
#ifdef __MLALGO
            // mohan add 2021-08-04
            if (PARAM.inp.deepks_scf)
            {
                fcs(iat, i) += parts.fvnl_dalpha(iat, i);
            }
#endif
        }
    }

    if (PARAM.inp.gate_flag || PARAM.inp.efield_flag)
    {
        GlobalV::ofs_running << "Atomic forces are not shifted if gate_flag or efield_flag == true!" << std::endl;
    }

    // pengfei 2016-12-20
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        this->forceSymmetry(ucell, fcs, symm);
    }

    // The net force should be evaluated AFTER the symmetrization.
    // With symmetry switched on, the forces assembled above are built from IBZ-reduced
    // quantities and only become physical after the symmetrization, forceSymmetry(). 
    // Force symmetrization is linear, so it commutes with the removal of a
    // uniform shift: the resulting fcs is identical to the previous ordering.
    for (int i = 0; i < 3; i++)
    {
        double sum = 0.0;

        for (int iat = 0; iat < nat; iat++)
        {
            // sum total force for correction
            sum += fcs(iat, i);
        }
        net_force[i]=sum;
        if (!(PARAM.inp.gate_flag || PARAM.inp.efield_flag))
        {
            for (int iat = 0; iat < nat; ++iat)
            {
                fcs(iat, i) -= sum / nat;
            }
        }
    }

    // compute forces using the DeePKS model
    deepks.write_forces(fcs, parts.fvnl_dalpha, PARAM.inp);

    // print Rydberg force or not
    bool ry = false;
    if (istestf)
    {
        // test
        // ModuleBase::matrix fvlocal;
        // fvlocal.create(nat,3);
        ModuleBase::matrix ftvnl;
        ftvnl.create(nat, 3);
        for (int iat = 0; iat < nat; iat++)
        {
            for (int i = 0; i < 3; i++)
            {
                // fvlocal(iat,i) = parts.fvl_dphi(iat,i) + parts.fvl_dvl(iat,i);
                ftvnl(iat, i) = parts.ftvnl_dphi(iat, i) + parts.fvnl_dbeta(iat, i);
            }
        }

        GlobalV::ofs_running << "\n PARTS OF FORCE: " << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
        GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(8) << std::endl;
        //-----------------------------
        // regular force terms test.
        //-----------------------------
        // this->print_force("OVERLAP    FORCE",parts.foverlap,1,ry);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "OVERLAP    FORCE", parts.foverlap, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "TVNL_DPHI  force",parts.ftvnl_dphi,false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "VNL_DBETA  force",parts.fvnl_dbeta,false);
        // this->print_force("T_VNL      FORCE",ftvnl,1,ry);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "T_VNL      FORCE", ftvnl, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "VL_dPHI    FORCE", parts.fvl_dphi, false);
        // this->print_force("VL_dPHI    FORCE",parts.fvl_dphi,1,ry);
        // this->print_force("VL_dVL     FORCE",parts.fvl_dvl,1,ry);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "VL_dVL     FORCE", parts.fvl_dvl, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "EWALD      FORCE", parts.fewalds, false);
        // this->print_force("VLOCAL     FORCE",fvlocal,PARAM.inp.test_force);
        // this->print_force("EWALD      FORCE",parts.fewalds,1,ry);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "NLCC       FORCE", parts.fcc, false);
        ModuleIO::print_force(GlobalV::ofs_running, ucell, "SCC        FORCE", parts.fscc, false);
        // this->print_force("NLCC       FORCE",parts.fcc,1,ry);
        // this->print_force("SCC        FORCE",parts.fscc,1,ry);
        //-------------------------------
        // put extra force here for test!
        //-------------------------------
        if (PARAM.inp.efield_flag)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "EFIELD     FORCE", parts.fefield, false);
            // this->print_force("EFIELD     FORCE",parts.fefield,1,ry);
        }
        if (PARAM.inp.esolver_type == "tddft")
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "EFIELD_TDDFT     FORCE", parts.fefield_tddft, false);
            // this->print_force("EFIELD_TDDFT     FORCE",parts.fefield_tddft,1,ry);
        }
        if (PARAM.inp.gate_flag)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "GATEFIELD     FORCE", parts.fgate, false);
            // this->print_force("GATEFIELD     FORCE",parts.fgate,1,ry);
        }
        if (PARAM.inp.imp_sol)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "IMP_SOL     FORCE", parts.fsol, false);
            // this->print_force("IMP_SOL     FORCE",parts.fsol,1,ry);
        }
        if (vdw_result != nullptr)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "VDW        FORCE", parts.force_vdw, false);
            // this->print_force("VDW        FORCE",parts.force_vdw,1,ry);
        }
        if (PARAM.inp.dft_plus_u)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "DFT+U      FORCE", parts.force_u, false);
        }
        if (PARAM.inp.sc_mag_switch)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "DeltaSpin  FORCE", parts.force_dspin, false);
        }
#ifdef __MLALGO
        // caoyu add 2021-06-03
        if (PARAM.inp.deepks_scf)
        {
            ModuleIO::print_force(GlobalV::ofs_running, ucell, "DeePKS     FORCE", parts.fvnl_dalpha, true);
        }
#endif
    }

    GlobalV::ofs_running << std::setiosflags(std::ios::left);

    // this->printforce_total(ry, istestf, fcs);
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "TOTAL-FORCE (eV/Angstrom)", fcs, false);
    net_force*= ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Net force vector (eV/Ang)", net_force.x, net_force.y, net_force.z);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Total drift (ev/Ang)", net_force.norm());
    if (istestf)
    {
        GlobalV::ofs_running << "\n FORCE INVALID TABLE." << std::endl;
        GlobalV::ofs_running << " " << std::setw(8) << "atom" << std::setw(5) << "x" << std::setw(5) << "y"
                             << std::setw(5) << "z" << std::endl;
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            GlobalV::ofs_running << " " << std::setw(8) << iat;
            for (int i = 0; i < 3; i++)
            {
                if (std::abs(fcs(iat, i) * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A)
                    < Force_Stress_LCAO::force_invalid_threshold_ev)
                {
                    fcs(iat, i) = 0.0;
                    GlobalV::ofs_running << std::setw(5) << "1";
                }
                else
                {
                    GlobalV::ofs_running << std::setw(5) << "0";
                }
            }
            GlobalV::ofs_running << std::endl;
        }
    }
}

template <typename T>
void Force_Stress_LCAO<T>::assemble_and_print_stress(const UnitCell& ucell,
                                                     const bool istests,
                                                     const vdw::VdwResult* vdw_result,
                                                     const Exx_Info& exx_info,
                                                     ModuleSymmetry::Symmetry* symm,
                                                     Setup_DeePKS<T>& deepks,
                                                     const LCAOStressParts& sparts,
                                                     ModuleBase::matrix& scs)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            scs(i, j) += sparts.soverlap(i, j) + sparts.stvnl_dphi(i, j) + sparts.svnl_dbeta(i, j) + sparts.svl_dphi(i, j)
                         + sparts.sigmadvl(i, j)  // derivative of local potential stress (pw)
                         + sparts.sigmaewa(i, j)  // ewald stress (pw)
                         + sparts.sigmacc(i, j)   // nonlinear core correction stress (pw)
                         + sparts.sigmaxc(i, j)   // exchange corretion stress
                         + sparts.sigmahar(i, j); // hartree stress

            // VDW stress from linpz and jiyy
            if (vdw_result != nullptr)
            {
                scs(i, j) += sparts.stress_vdw(i, j);
            }
            // DFT plus U stress from qux
            if (PARAM.inp.dft_plus_u)
            {
                scs(i, j) += sparts.stress_u(i, j);
            }
            if (PARAM.inp.sc_mag_switch)
            {
                scs(i, j) += sparts.stress_dspin(i, j);
            }
#ifdef __EXX
            // Stress contribution from exx
            if (exx_info.info_global.cal_exx)
            {
                scs(i, j) += sparts.stress_exx(i, j);
            }
#endif
#ifdef __MLALGO
            if (PARAM.inp.deepks_scf)
            {
                scs(i, j) += sparts.svnl_dalpha(i, j);
            }
#endif
        }
    }
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        symm->symmetrize_mat3(scs, ucell.lat);
    } // end symmetry

    deepks.write_stress(scs, sparts.svnl_dalpha, ucell.omega, PARAM.inp);

    // print Rydberg stress or not
    bool ry = false;

    // test stress each terms if needed
    if (istests)
    {
        // test
        ModuleBase::matrix svlocal;
        svlocal.create(3, 3);
        ModuleBase::matrix stvnl;
        stvnl.create(3, 3);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                svlocal(i, j) = sparts.svl_dphi(i, j) + sparts.sigmadvl(i, j);
                stvnl(i, j) = sparts.stvnl_dphi(i, j) + sparts.svnl_dbeta(i, j);
            }
        }

        const bool screen = PARAM.inp.test_stress;

        GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
        GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(10) << std::endl;
        ModuleIO::print_stress("OVERLAP  STRESS", sparts.soverlap, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("T        STRESS", sparts.stvnl_dphi, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("VNL      STRESS", sparts.svnl_dbeta, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("T_VNL    STRESS", stvnl, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("VL_dPHI  STRESS", sparts.svl_dphi, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("VL_dVL   STRESS", sparts.sigmadvl, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("HAR      STRESS", sparts.sigmahar, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("EWALD    STRESS", sparts.sigmaewa, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("cc       STRESS", sparts.sigmacc, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("XC       STRESS", sparts.sigmaxc, screen, ry, GlobalV::ofs_running);
        if (vdw_result != nullptr)
        {
            ModuleIO::print_stress("VDW      STRESS", sparts.stress_vdw, screen, ry, GlobalV::ofs_running);
        }
        if (PARAM.inp.dft_plus_u)
        {
            ModuleIO::print_stress("DFTU     STRESS", sparts.stress_u, screen, ry, GlobalV::ofs_running);
        }
        if (PARAM.inp.sc_mag_switch)
        {
            ModuleIO::print_stress("DeltaSpin  STRESS", sparts.stress_dspin, screen, ry, GlobalV::ofs_running);
        }
#ifdef __EXX
        if (exx_info.info_global.cal_exx)
        {
            ModuleIO::print_stress("EXX      STRESS", sparts.stress_exx, screen, ry, GlobalV::ofs_running);
        }
#endif
        ModuleIO::print_stress("TOTAL    STRESS", scs, screen, ry, GlobalV::ofs_running);
    } // end of test

    GlobalV::ofs_running << std::setiosflags(std::ios::left);

    // print total stress
    bool screen_normal = true;
    ModuleIO::print_stress("TOTAL-STRESS", scs, screen_normal, ry, GlobalV::ofs_running);

    double unit_transform = 0.0;
    unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {PARAM.inp.press1, PARAM.inp.press2, PARAM.inp.press3};

    for (int i = 0; i < 3; i++)
    {
        scs(i, i) -= external_stress[i] / unit_transform;
    }
}

// do symmetry for total force
template <typename T>
void Force_Stress_LCAO<T>::forceSymmetry(const UnitCell& ucell, ModuleBase::matrix& fcs, ModuleSymmetry::Symmetry* symm)
{
    double d1, d2, d3;
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        ModuleBase::Mathzone::Cartesian_to_Direct(fcs(iat, 0), fcs(iat, 1), fcs(iat, 2),
          ucell.a1.x, ucell.a1.y, ucell.a1.z, ucell.a2.x, ucell.a2.y, ucell.a2.z,
          ucell.a3.x, ucell.a3.y, ucell.a3.z, d1, d2, d3);

        fcs(iat, 0) = d1;
        fcs(iat, 1) = d2;
        fcs(iat, 2) = d3;
    }
    symm->symmetrize_vec3_nat(fcs.c);
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        ModuleBase::Mathzone::Direct_to_Cartesian(fcs(iat, 0), fcs(iat, 1), fcs(iat, 2),
          ucell.a1.x, ucell.a1.y, ucell.a1.z, ucell.a2.x, ucell.a2.y, ucell.a2.z,
          ucell.a3.x, ucell.a3.y, ucell.a3.z, d1, d2, d3);

        fcs(iat, 0) = d1;
        fcs(iat, 1) = d2;
        fcs(iat, 2) = d3;
    }
    return;
}

template class Force_Stress_LCAO<double>;
template class Force_Stress_LCAO<std::complex<double>>;
