#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_esolver/esolver_ks_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_io/berryphase.h"
#include "module_io/get_pchg_lcao.h"
#include "module_io/get_wf_lcao.h"
#include "module_io/io_npz.h"
#include "module_io/to_wannier90_lcao.h"
#include "module_io/to_wannier90_lcao_in_pw.h"
#include "module_io/write_HS_R.h"
#include "module_parameter/parameter.h"
#include "module_elecstate/elecstate_tools.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_base/formatter.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_io/cube_io.h"
#include "module_io/read_wfc_nao.h"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_wfc_nao.h"
#ifdef __EXX
#include "module_io/restart_exx_csr.h"
#endif

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_scf");

    //! 1) call before_scf() of ESolver_KS.
    ESolver_KS<TK>::before_scf(ucell, istep);

    //! 2) find search radius
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                   PARAM.inp.out_level,
                                                   orb_.get_rcutmax_Phi(),
                                                   ucell.infoNL.get_rcutmax_Beta(),
                                                   PARAM.globalv.gamma_only_local);

    //! 3) use search_radius to search adj atoms
    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         this->gd,
                         ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    //! 4) initialize NAO basis set 
    double dr_uniform = 0.001;
    std::vector<double> rcuts;
    std::vector<std::vector<double>> psi_u;
    std::vector<std::vector<double>> dpsi_u;
    std::vector<std::vector<double>> d2psi_u;

    Gint_Tools::init_orb(dr_uniform, rcuts, ucell, orb_, psi_u, dpsi_u, d2psi_u);

    //! 5) set periodic boundary conditions
    this->GridT.set_pbc_grid(this->pw_rho->nx,
                             this->pw_rho->ny,
                             this->pw_rho->nz,
                             this->pw_big->bx,
                             this->pw_big->by,
                             this->pw_big->bz,
                             this->pw_big->nbx,
                             this->pw_big->nby,
                             this->pw_big->nbz,
                             this->pw_big->nbxx,
                             this->pw_big->nbzp_start,
                             this->pw_big->nbzp,
                             this->pw_rho->ny,
                             this->pw_rho->nplane,
                             this->pw_rho->startz_current,
                             ucell,
                             this->gd,
                             dr_uniform,
                             rcuts,
                             psi_u,
                             dpsi_u,
                             d2psi_u,
                             PARAM.inp.nstream);

    //! 6) prepare grid integral
#ifdef __NEW_GINT
    auto gint_info = std::make_shared<ModuleGint::GintInfo>(
        this->pw_big->nbx,
        this->pw_big->nby,
        this->pw_big->nbz,
        this->pw_rho->nx,
        this->pw_rho->ny,
        this->pw_rho->nz,
        0,
        0,
        this->pw_big->nbzp_start,
        this->pw_big->nbx,
        this->pw_big->nby,
        this->pw_big->nbzp,
        orb_.Phi,
        ucell,
        this->gd);
    ModuleGint::Gint::init_gint_info(gint_info);
#endif

    psi_u.clear();
    psi_u.shrink_to_fit();
    dpsi_u.clear();
    dpsi_u.shrink_to_fit();
    d2psi_u.clear();
    d2psi_u.shrink_to_fit();

    // 7) For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    this->RA.for_2d(ucell, this->gd, this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());


    // 8) initialize electronic wave function psi
    // this code belongs to before_all_runners? mohan add 2025-03-10
    if (this->psi == nullptr)
    {
        int nsk = 0;
        int ncol = 0;
        if (PARAM.globalv.gamma_only_local)
        {
            nsk = PARAM.inp.nspin;
            ncol = this->pv.ncol_bands;
            if (PARAM.inp.ks_solver == "genelpa" || PARAM.inp.ks_solver == "elpa" || PARAM.inp.ks_solver == "lapack"
                || PARAM.inp.ks_solver == "pexsi" || PARAM.inp.ks_solver == "cusolver"
                || PARAM.inp.ks_solver == "cusolvermp")
            {
                ncol = this->pv.ncol;
            }
        }
        else
        {
            nsk = this->kv.get_nks();
#ifdef __MPI
            ncol = this->pv.ncol_bands;
#else
            ncol = PARAM.inp.nbands;
#endif
        }
        this->psi = new psi::Psi<TK>(nsk, ncol, this->pv.nrow, this->kv.ngk, true);
    }

    // 9) read psi from file
    // this code belongs to before_all_runners? mohan add 2025-03-10
    if (istep == 0 && PARAM.inp.init_wfc == "file")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir, this->pv, *(this->psi), this->pelec))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "read electronic wave functions failed");
        }
    }


    // 10) after ions move, prepare grid in Gint
    LCAO_domain::grid_prepare(this->GridT, this->GG, this->GK, ucell, orb_, *this->pw_rho, *this->pw_big);


    // 11) initialize the Hamiltonian operators
    // if atom moves, then delete old pointer and add a new one
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr)
    {
        elecstate::DensityMatrix<TK, double>* DM = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();

        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(
            PARAM.globalv.gamma_only_local ? &(this->GG) : nullptr,
            PARAM.globalv.gamma_only_local ? nullptr : &(this->GK),
            ucell,
            this->gd,
            &this->pv,
            this->pelec->pot,
            this->kv,
            two_center_bundle_,
            orb_,
            DM
#ifdef __DEEPKS
            ,
            &this->ld
#endif
#ifdef __EXX
            ,
            istep,
            GlobalC::exx_info.info_ri.real_number ? &this->exd->two_level_step : &this->exc->two_level_step,
            GlobalC::exx_info.info_ri.real_number ? &exx_lri_double->Hexxs : nullptr,
            GlobalC::exx_info.info_ri.real_number ? nullptr : &exx_lri_complex->Hexxs
#endif
        );
    }




#ifdef __DEEPKS
    // 12) for each ionic step, the overlap <phi|alpha> must be rebuilt
    // since it depends on ionic positions
    if (PARAM.globalv.deepks_setorb)
    {
        const Parallel_Orbitals* pv = &this->pv;
        // allocate <phi(0)|alpha(R)>, phialpha is different every ion step, so it is allocated here
        DeePKS_domain::allocate_phialpha(PARAM.inp.cal_force, ucell, orb_, this->gd, pv, this->ld.phialpha);
        // build and save <phi(0)|alpha(R)> at beginning
        DeePKS_domain::build_phialpha(PARAM.inp.cal_force,
                                      ucell,
                                      orb_,
                                      this->gd,
                                      pv,
                                      *(two_center_bundle_.overlap_orb_alpha),
                                      this->ld.phialpha);

        if (PARAM.inp.deepks_out_unittest)
        {
            DeePKS_domain::check_phialpha(PARAM.inp.cal_force,
                                          ucell,
                                          orb_,
                                          this->gd,
                                          pv,
                                          this->ld.phialpha,
                                          GlobalV::MY_RANK);
        }
    }
#endif

    // 13) prepare sc calculation
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.init_sc(PARAM.inp.sc_thr,
                   PARAM.inp.nsc,
                   PARAM.inp.nsc_min,
                   PARAM.inp.alpha_trial,
                   PARAM.inp.sccut,
                   PARAM.inp.sc_drop_thr,
                   ucell,
                   &(this->pv),
                   PARAM.inp.nspin,
                   this->kv,
                   this->p_hamilt,
                   this->psi,
                   this->pelec);
    }

    // 14) set xc type before the first cal of xc in pelec->init_scf
    // Peize Lin add 2016-12-03
#ifdef __EXX
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exd->exx_beforescf(istep, this->kv, *this->p_chgmix, ucell, orb_);
        }
        else
        {
            this->exc->exx_beforescf(istep, this->kv, *this->p_chgmix, ucell, orb_);
        }
    }
#endif

    // 15) init_scf, should be before_scf? mohan add 2025-03-10
    this->pelec->init_scf(istep, ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, ucell.symm);


    // 16) initalize DMR
    // DMR should be same size with Hamiltonian(R)
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
        ->get_DM()
        ->init_DMR(*(dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt)->getHR()));


    // 17) two cases are considered:
    // 1. DMK in DensityMatrix is not empty (istep > 0), then DMR is initialized by DMK
    // 2. DMK in DensityMatrix is empty (istep == 0), then DMR is initialized by zeros
    if (istep > 0)
    {
        dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->cal_DMR();
    }

    if (PARAM.inp.dm_to_rho)
    {
        // file name of DM
        std::string zipname = "output_DM0.npz";
        elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
      
        // read DM from file
        ModuleIO::read_mat_npz(&(this->pv), ucell, zipname, *(dm->get_DMR_pointer(1)));

        // if nspin=2, need extra reading
        if (PARAM.inp.nspin == 2)
        {
            zipname = "output_DM1.npz";
            ModuleIO::read_mat_npz(&(this->pv), ucell, zipname, *(dm->get_DMR_pointer(2)));
        }

        elecstate::calculate_weights(this->pelec->ekb,
                                     this->pelec->wg,
                                     this->pelec->klist,
                                     this->pelec->eferm,
                                     this->pelec->f_en,
                                     this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
      
        this->pelec->psiToRho(*this->psi);

        int nspin0 = PARAM.inp.nspin == 2 ? 2 : 1;

        for (int is = 0; is < nspin0; is++)
        {
            std::string fn = PARAM.globalv.global_out_dir + "/SPIN" + std::to_string(is + 1) + "_CHG.cube";

            // write electron density
            ModuleIO::write_vdata_palgrid(this->Pgrid,
                                          this->chr.rho[is],
                                          is,
                                          PARAM.inp.nspin,
                                          istep,
                                          fn,
                                          this->pelec->eferm.get_efval(is),
                                          &(ucell),
                                          3,
                                          1);
        }

        // why we need to return here? mohan add 2025-03-10
        ModuleBase::timer::tick("ESolver_KS_LCAO", "before_scf");
        return;
    }

    // 18) the electron charge density should be symmetrized,
    // here is the initialization
    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, this->chr, this->pw_rho, ucell.symm);
    }

    // 19) why we need to set this sentence? mohan add 2025-03-10
    this->p_hamilt->non_first_scf = istep;


    // 20) update of RDMFT, added by jghan 
    if (PARAM.inp.rdmft == true)
    {
        // necessary operation of these parameters have be done with p_esolver->Init() in source/driver_run.cpp
        rdmft_solver.update_ion(ucell,
                                *(this->pw_rho),
                                this->locpp.vloc,
                                this->sf.strucFac);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "before_scf");
    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
