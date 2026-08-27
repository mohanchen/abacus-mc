#include "vnl_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/clebsch_gordan_coeff.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/output.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_device/device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/kernels/vnl_op.h"

#include "source_base/parallel_comm.h" // use POOL_WORLD


pseudopot_cell_vnl::pseudopot_cell_vnl()
{
    this->use_gpu_ = (PARAM.inp.device == "gpu");
}

pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
    delete[] indv_ijkb0;
}

void pseudopot_cell_vnl::release_memory()
{
    if (this->nhm <= 0 || memory_released) {
        return;
}
    if (this->use_gpu_)
    {
        delmem_sd_op()(this->s_deeq);
        delmem_sd_op()(this->s_nhtol);
        delmem_sd_op()(this->s_nhtolm);
        delmem_sd_op()(this->s_indv);
        delmem_sd_op()(this->s_tab);
        delmem_sd_op()(this->s_qq_nt);
        delmem_cd_op()(this->c_deeq_nc);
        delmem_cd_op()(this->c_vkb);
        delmem_cd_op()(this->c_qq_so);
        delmem_zd_op()(this->z_deeq_nc);
        delmem_zd_op()(this->z_qq_so);
        delmem_dd_op()(this->d_deeq);
        delmem_zd_op()(this->z_vkb);
        delmem_dd_op()(this->d_tab);
        delmem_dd_op()(this->d_indv);
        delmem_dd_op()(this->d_nhtol);
        delmem_dd_op()(this->d_nhtolm);
        delmem_dd_op()(this->d_qq_nt);
    }
    else
    {
        delmem_sh_op()(this->s_deeq);
        delmem_sh_op()(this->s_nhtol);
        delmem_sh_op()(this->s_nhtolm);
        delmem_sh_op()(this->s_indv);
        delmem_sh_op()(this->s_tab);
        delmem_sh_op()(this->s_qq_nt);
        delmem_ch_op()(this->c_deeq_nc);
        delmem_ch_op()(this->c_vkb);
        delmem_ch_op()(this->c_qq_so);
#ifdef __DSP
        if (this->z_vkb != nullptr)
        {
            base_device::memory::delete_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>()(this->z_vkb);
            this->z_vkb = nullptr;
        }
#endif
        // There's no need to delete double precision pointers while in a CPU environment.
    }
    memory_released = true;
}

//-----------------------------------
// setup lmaxkb, nhm, nkb, lmaxq
// allocate vkb, PARAM.globalv.nqx, tab, tab_at
//-----------------------------------
void pseudopot_cell_vnl::init(const UnitCell& ucell,
                              Structure_Factor* psf_in,
                              const ModulePW::PW_Basis_K* wfc_basis,
                              const bool allocate_vkb)
{
    const int ntype = ucell.ntype;
    ModuleBase::TITLE("pseudopot_cell_vnl", "init");
    ModuleBase::timer::start("ppcell_vnl", "init");

    GlobalV::ofs_running << "\n SETUP NONLOCAL PSEUDOPOTENTIALS IN PLANE WAVE BASIS" << std::endl;

    int it = 0;
    this->wfcpw = wfc_basis;
    this->psf = psf_in;
    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : lmaxkb(max angular momentum,(see pseudo_h))
    //----------------------------------------------------------
    this->lmaxkb = -1;
    for (it = 0; it < ntype; it++)
    {
        GlobalV::ofs_running << " " << ucell.atoms[it].label << " non-local projectors:" << std::endl;
        for (int ibeta = 0; ibeta < ucell.atoms[it].ncpp.nbeta; ibeta++)
        {
            GlobalV::ofs_running << " projector " << ibeta + 1 << " L=" << ucell.atoms[it].ncpp.lll[ibeta]
                                 << std::endl;
            this->lmaxkb = std::max(this->lmaxkb, ucell.atoms[it].ncpp.lll[ibeta]);
        }
    }

    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : nhm(max number of different beta functions per atom)
    // NAME : nbetam(max number of beta functions)
    // NAME : nwfcm(max number of wavefunctions)
    //----------------------------------------------------------
    this->nhm = 0;
    this->nbetam = 0;
    int nwfcm = 0;
    for (it = 0; it < ntype; it++)
    {
        this->nhm = std::max(nhm, ucell.atoms[it].ncpp.nh);
        this->nbetam = std::max(nbetam, ucell.atoms[it].ncpp.nbeta);
        nwfcm = std::max(nwfcm, ucell.atoms[it].ncpp.nchi);
    }

    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : nkb(total number of beta functions, with struct.fact.)
    //----------------------------------------------------------
    this->nkb = 0;
    for (it = 0; it < ntype; it++)
    {
        this->nkb += ucell.atoms[it].ncpp.nh * ucell.atoms[it].na;
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "TOTAL NUMBER OF NONLOCAL PROJECTORS", nkb);

    if (this->nhm > 0)
    {
        this->indv.create(ntype, this->nhm);
        this->nhtol.create(ntype, this->nhm);
        this->nhtolm.create(ntype, this->nhm);
        this->nhtoj.create(ntype, this->nhm);
        this->deeq.create(PARAM.inp.nspin, ucell.nat, this->nhm, this->nhm);
        this->deeq_nc.create(PARAM.inp.nspin, ucell.nat, this->nhm, this->nhm);
        this->qq_nt.create(ntype, this->nhm, this->nhm);
        this->qq_so.create(ntype, 4, this->nhm, this->nhm);
        if (this->use_gpu_)
        {
            if (PARAM.globalv.has_float_data)
            {
                resmem_sd_op()(s_deeq, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
                resmem_sd_op()(s_nhtol, ntype * this->nhm);
                resmem_sd_op()(s_nhtolm, ntype * this->nhm);
                resmem_sd_op()(s_indv, ntype * this->nhm);
                resmem_sd_op()(s_qq_nt, ntype * this->nhm * this->nhm);
                resmem_cd_op()(c_deeq_nc, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
                resmem_cd_op()(c_qq_so, ntype * 4 * this->nhm * this->nhm);
            }
            if (PARAM.globalv.has_double_data)
            {
                resmem_zd_op()(z_deeq_nc, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
                resmem_zd_op()(z_qq_so, ntype * 4 * this->nhm * this->nhm);
            }
            resmem_dd_op()(d_deeq, PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm);
            resmem_dd_op()(d_indv, ntype * this->nhm);
            resmem_dd_op()(d_nhtol, ntype * this->nhm);
            resmem_dd_op()(d_nhtolm, ntype * this->nhm);
            resmem_dd_op()(d_qq_nt, ntype * this->nhm * this->nhm);
        }
        else
        {
            if (PARAM.globalv.has_float_data)
            {
                resmem_sh_op()(s_deeq,
                               PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm,
                               "VNL::s_deeq");
                resmem_sh_op()(s_nhtol, ntype * this->nhm, "VNL::s_nhtol");
                resmem_sh_op()(s_nhtolm, ntype * this->nhm, "VNL::s_nhtolm");
                resmem_sh_op()(s_indv, ntype * this->nhm, "VNL::s_indv");
                resmem_sh_op()(s_qq_nt, ntype * this->nhm * this->nhm, "VNL::s_qq_nt");
                resmem_ch_op()(c_deeq_nc,
                               PARAM.inp.nspin * ucell.nat * this->nhm * this->nhm,
                               "VNL::c_deeq_nc");
                resmem_ch_op()(c_qq_so, ntype * 4 * this->nhm * this->nhm, "VNL::c_qq_so");
            }
            if (PARAM.globalv.has_double_data)
            {
                this->z_deeq_nc = this->deeq_nc.ptr;
                this->z_qq_so = this->qq_so.ptr;
            }
            this->d_deeq = this->deeq.ptr;
            this->d_indv = this->indv.c;
            this->d_nhtol = this->nhtol.c;
            this->d_nhtolm = this->nhtolm.c;
            this->d_qq_nt = this->qq_nt.ptr;
            // There's no need to delete double precision pointers while in a CPU environment.
        }
        this->dvan.create(ntype, this->nhm, this->nhm);
        this->dvan_so.create(PARAM.inp.nspin, ntype, this->nhm, this->nhm);

        this->ijtoh.create(ntype, this->nhm, this->nhm);
        this->qq_at.create(ucell.nat, this->nhm, this->nhm);
    }
    else
    {
        GlobalV::ofs_running << "\n nhm = 0, not allocate some matrix.";
    }

    // nqxq = ((sqrt(gcutm)+sqrt(xqq[1]*xqq[1]+xqq[2]*xqq[2]+xqq[3]*xqq[3])/
    // dq+4)*cell_factor;
    this->lmaxq = 2 * this->lmaxkb + 1;
    int npwx = this->wfcpw->npwk_max;
    this->vkbnc = npwx;
    if (nkb > 0 && allocate_vkb)
    {
        if (!this->use_gpu_)
    {
        vkb.create(nkb, npwx);
        ModuleBase::Memory::record("VNL::vkb", nkb * npwx * sizeof(std::complex<double>));
        }
        // GPU path: vkb ComplexMatrix is not allocated.
        // Column dimension is stored in vkbnc for gemm/gemv leading dimension.
        // Actual GPU buffers (c_vkb/z_vkb) are allocated below.
    }

    // this->nqx = 10000;		// calculted in allocate_nlpot.f90
    // PARAM.globalv.nqxq = static_cast<int>(((sqrt(INPUT.ecutrho) + qnorm) / PARAM.globalv.dq + 4.0) * cell_factor);

    // mohan update 2021-02-22
    // liuyu update 2023-09-28
    if (nbetam > 0)
    {
        const int nbrx_nc = 2 * nbetam;
        // nbetam: max number of beta functions
        if (PARAM.inp.nspin != 4)
        {
            this->tab.create(ntype, nbetam, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab", ntype * nbetam * PARAM.globalv.nqx * sizeof(double));
        }
        else
        {
            this->tab.create(ntype, nbrx_nc, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab", ntype * nbrx_nc * PARAM.globalv.nqx * sizeof(double));
        }

        if (lmaxq > 0)
        {
            this->qrad.create(ntype, lmaxq, nbetam * (nbetam + 1) / 2, PARAM.globalv.nqxq);
        }
    }

    // mohan update 2021-02-22
    // liuyu update 2023-09-28
    if (nwfcm > 0)
    {
        int nchix_nc = 2 * nwfcm;
        // nwfcm : max number of atomic wavefunctions per atom
        if (PARAM.inp.nspin != 4)
        {
            this->tab_at.create(ntype, nwfcm, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab_at", ntype * nwfcm * PARAM.globalv.nqx * sizeof(double));
        }
        else
        {
            this->tab_at.create(ntype, nchix_nc, PARAM.globalv.nqx);
            ModuleBase::Memory::record("VNL::tab_at", ntype * nchix_nc * PARAM.globalv.nqx * sizeof(double));
        }
    }
    if (this->use_gpu_)
    {
        if (PARAM.globalv.has_float_data)
        {
            resmem_sd_op()(s_tab, this->tab.getSize());
            resmem_cd_op()(c_vkb, nkb * npwx);
        }
        resmem_zd_op()(z_vkb, nkb * npwx);
        resmem_dd_op()(d_tab, this->tab.getSize());
    }
    else
    {
        if (PARAM.globalv.has_float_data)
        {
            resmem_sh_op()(s_tab, this->tab.getSize());
            resmem_ch_op()(c_vkb, nkb * npwx);
        }
#ifdef __DSP
        base_device::memory::resize_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>()
        (this->z_vkb, this->vkb.size, "VNL::z_vkb");
        // memcpy(this->z_vkb,this->vkb.c,this->vkb.size*16);
#else
        this->z_vkb = this->vkb.c;
#endif
        this->d_tab = this->tab.ptr;
        // There's no need to delete double precision pointers while in a CPU environment.
    }

    ModuleBase::timer::end("ppcell_vnl", "init");
    return;
}

void pseudopot_cell_vnl::init_vnl(UnitCell& cell, const ModulePW::PW_Basis* rho_basis)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "init_vnl");
    ModuleBase::timer::start("ppcell_vnl", "init_vnl");

    this->omega_old = cell.omega;

    // from init_us_1
    //    a) For each non vanderbilt pseudopotential it computes the D and
    //       the betar in the same form of the Vanderbilt pseudopotential.
    //    b) It computes the indices indv which establish the correspondence
    //       nh <-> beta in the atom
    //    c) It computes the indices nhtol which establish the correspondence
    //       nh <-> angular momentum of the beta function
    //    d) It computes the indices nhtolm which establish the correspondence
    //       nh <-> combined (l,m) index for the beta function.

    // the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
    // but in some versions of the PP files lmax is not set to the maximum
    // l of the beta functions but includes the l of the local potential
    for (int it = 0; it < cell.ntype; it++)
    {
        if (cell.atoms[it].ncpp.tvanp)
        {
            cell.atoms[it].ncpp.nqlc = std::min(cell.atoms[it].ncpp.nqlc, lmaxq);
            if (cell.atoms[it].ncpp.nqlc < 0) {
                cell.atoms[it].ncpp.nqlc = 0;
}
        }
    }

    // In the spin-orbit case we need the unitary matrix u which rotates the
    // real spherical harmonics and yields the complex ones.
    soc.fcoef.create(cell.ntype, this->nhm, this->nhm);
    if (PARAM.inp.lspinorb)
    {
        soc.rot_ylm(this->lmaxkb);
    }

    // For each pseudopotential we initialize the indices nhtol, nhtolm,
    // nhtoj, indv, and if the pseudopotential is of KB type we initialize
    // the atomic D terms
    this->dvan.zero_out();
    this->dvan_so.zero_out(); // added by zhengdy-soc
    delete[] indv_ijkb0;
    this->indv_ijkb0 = new int[cell.nat];
    int ijkb0 = 0;
    for (int it = 0; it < cell.ntype; it++)
    {
        int BetaIndex = 0;
        const int Nprojectors = cell.atoms[it].ncpp.nh;
        for (int ib = 0; ib < cell.atoms[it].ncpp.nbeta; ib++)
        {
            const int l = cell.atoms[it].ncpp.lll[ib];
            const double j = cell.atoms[it].ncpp.jjj[ib];
            for (int m = 0; m < 2 * l + 1; m++)
            {
                this->nhtol(it, BetaIndex) = l;
                this->nhtolm(it, BetaIndex) = l * l + m;
                this->nhtoj(it, BetaIndex) = j;
                this->indv(it, BetaIndex) = ib;
                ++BetaIndex;
            }
        }

        // ijtoh map augmentation channel indexes ih and jh to composite
        // "triangular" index ijh
        for (int ih1 = 0; ih1 < nhm; ih1++)
        {
            for (int ih2 = ih1; ih2 < nhm; ih2++)
            {
                this->ijtoh(it, ih1, ih2) = -1;
                this->ijtoh(it, ih2, ih1) = -1;
            }
        }
        int ijv = 0;
        for (int ih1 = 0; ih1 < Nprojectors; ih1++)
        {
            for (int ih2 = ih1; ih2 < Nprojectors; ih2++)
            {
                ijv++; // liuyu I'm not sure from 0 or 1
                this->ijtoh(it, ih1, ih2) = ijv;
                this->ijtoh(it, ih2, ih1) = ijv;
            }
        }

        // ijkb0 points to the first beta "in the solid" for atom ia
        // i.e. ijkb0,.. ijkb0+nh(ityp(ia))-1 are the nh beta functions of
        // atom ia in the global list of beta functions (ijkb0=0 for ia=0)
        for (int ia = 0; ia < cell.nat; ia++)
        {
            if (it == cell.iat2it[ia])
            {
                this->indv_ijkb0[ia] = ijkb0;
                ijkb0 += Nprojectors;
            }
        }

        // From now on the only difference between KB and US pseudopotentials
        // is in the presence of the q and Q functions.
        // Here we initialize the D of the solid
        if (cell.atoms[it].ncpp.has_so)
        {
            for (int ip = 0; ip < Nprojectors; ip++)
            {
                const int l1 = this->nhtol(it, ip);
                const double j1 = this->nhtoj(it, ip);
                const int m1 = this->nhtolm(it, ip) - l1 * l1;
                // const int v1 = static_cast<int>( indv(it, ip ) );
                for (int ip2 = 0; ip2 < Nprojectors; ip2++)
                {
                    const int l2 = this->nhtol(it, ip2);
                    const double j2 = this->nhtoj(it, ip2);
                    const int m2 = this->nhtolm(it, ip2) - l2 * l2;
                    // const int v2 = static_cast<int>( indv(it, ip2 ) );
                    if (l1 == l2 && fabs(j1 - j2) < 1e-7)
                    {
                        for (int is1 = 0; is1 < 2; is1++)
                        {
                            for (int is2 = 0; is2 < 2; is2++)
                            {
                                soc.set_fcoef(l1, l2, is1, is2, m1, m2, j1, j2, it, ip, ip2);
                            }
                        }
                    }
                }
            }
            //
            //   and calculate the bare coefficients
            //
            for (int ip = 0; ip < Nprojectors; ++ip)
            {
                const int ir = static_cast<int>(indv(it, ip));
                for (int ip2 = 0; ip2 < Nprojectors; ++ip2)
                {
                    const int is = static_cast<int>(indv(it, ip2));
                    int ijs = 0;
                    for (int is1 = 0; is1 < 2; ++is1)
                    {
                        for (int is2 = 0; is2 < 2; ++is2)
                        {
                            this->dvan_so(ijs, it, ip, ip2)
                                = cell.atoms[it].ncpp.dion(ir, is) * soc.fcoef(it, is1, is2, ip, ip2);
                            ++ijs;
                            if (ir != is) {
                                soc.fcoef(it, is1, is2, ip, ip2) = std::complex<double>(0.0, 0.0);
}
                        }
                    }
                }
            }
        }
        else {
            for (int ip = 0; ip < Nprojectors; ip++)
            {
                for (int ip2 = 0; ip2 < Nprojectors; ip2++)
                {
                    if (this->nhtol(it, ip) == nhtol(it, ip2) && this->nhtolm(it, ip) == nhtolm(it, ip2))
                    {
                        const int ir = static_cast<int>(indv(it, ip));
                        const int is = static_cast<int>(indv(it, ip2));
                        if (PARAM.inp.lspinorb)
                        {
                            this->dvan_so(0, it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                            this->dvan_so(3, it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                        }
                        else
                        {
                            this->dvan(it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
                        }
                    }
                }
            }
}
    }

    // e) It computes the coefficients c_{LM}^{nm} which relates the
    //    spherical harmonics in the Q expansion
    // f) It computes the interpolation table "qrad" for Q(G)
    // g) It computes the qq terms which define the S matrix.

    // compute Clebsch-Gordan coefficients
    if (PARAM.globalv.use_uspp)
    {
        ModuleBase::Clebsch_Gordan::clebsch_gordan(lmaxkb + 1, ap, lpx, lpl);
    }

    // here for the US types we compute the Fourier transform of the Q functions.
    if (lmaxq > 0)
    {
        this->compute_qrad(cell);
    }

    // compute the qq coefficients by integrating the Q.
    // The qq are the g=0 components of Q
    if (rho_basis->ig_gge0 >= 0)
    {
        ModuleBase::matrix ylmk0(lmaxq * lmaxq, 1);
        const double qnorm = 0.0; // only G=0 term
        std::complex<double> qgm(0.0, 0.0);
        ModuleBase::YlmReal::Ylm_Real(lmaxq * lmaxq, 1, &(rho_basis->gcar[rho_basis->ig_gge0]), ylmk0);
        for (int it = 0; it < cell.ntype; it++)
        {
            Atom_pseudo* upf = &cell.atoms[it].ncpp;
            if (upf->tvanp)
            {
                if (upf->has_so)
                {
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = 0; jh < upf->nh; jh++)
                        {
                            this->radial_fft_q(1, ih, jh, it, &qnorm, ylmk0, &qgm);
                            this->qq_nt(it, ih, jh) = cell.omega * qgm.real();
                            for (int kh = 0; kh < upf->nh; kh++)
                            {
                                for (int lh = 0; lh < upf->nh; lh++)
                                {
                                    int ijs = 0;
                                    for (int is1 = 0; is1 < 2; ++is1)
                                    {
                                        for (int is2 = 0; is2 < 2; ++is2)
                                        {
                                            for (int is = 0; is < 2; is++)
                                            {
                                                this->qq_so(it, ijs, kh, lh) += cell.omega * qgm.real()
                                                                                * soc.fcoef(it, is1, is, kh, ih)
                                                                                * soc.fcoef(it, is, is2, jh, lh);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = ih; jh < upf->nh; jh++)
                        {
                            this->radial_fft_q(1, ih, jh, it, &qnorm, ylmk0, &qgm);
                            if (PARAM.inp.lspinorb)
                            {
                                this->qq_so(it, 0, ih, jh) = cell.omega * qgm.real();
                                this->qq_so(it, 0, jh, ih) = this->qq_so(it, 0, ih, jh);
                                this->qq_so(it, 3, ih, jh) = this->qq_so(it, 0, ih, jh);
                                this->qq_so(it, 3, jh, ih) = this->qq_so(it, 0, ih, jh);
                            }
                            this->qq_nt(it, ih, jh) = cell.omega * qgm.real();
                            this->qq_nt(it, jh, ih) = cell.omega * qgm.real();
                        }
                    }
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(this->qq_nt.ptr, this->qq_nt.getSize());
    Parallel_Reduce::reduce_pool(this->qq_so.ptr, this->qq_so.getSize());
#endif

    // set the atomic specific qq_at matrices
    for (int ia = 0; ia < cell.nat; ia++)
    {
        int it = cell.iat2it[ia];
        for (int ih = 0; ih < nhm; ih++)
        {
            for (int jh = ih; jh < nhm; jh++)
            {
                this->qq_at(ia, ih, jh) = qq_nt(it, ih, jh);
                this->qq_at(ia, jh, ih) = qq_nt(it, jh, ih);
            }
        }
    }

    // h) It fills the interpolation table for the beta functions
    /**********************************************************
    // He Lixin: this block is used for non-local potential
    // fill the interpolation table tab
    ************************************************************/

    const double pref = ModuleBase::FOUR_PI / sqrt(cell.omega);
    this->tab.zero_out();
    GlobalV::ofs_running << "\n Init Non-Local PseudoPotential table : ";
    for (int it = 0; it < cell.ntype; it++)
    {
        const int nbeta = cell.atoms[it].ncpp.nbeta;
        int kkbeta = cell.atoms[it].ncpp.kkbeta;

        // mohan modify 2008-3-31
        // mohan add kkbeta>0 2009-2-27
        if ((kkbeta % 2 == 0) && kkbeta > 0)
        {
            kkbeta--;
        }

        double* jl = new double[kkbeta];
        double* aux = new double[kkbeta];
        for (int ib = 0; ib < nbeta; ib++)
        {
            const int l = cell.atoms[it].ncpp.lll[ib];
            for (int iq = 0; iq < PARAM.globalv.nqx; iq++)
            {
                const double q = iq * PARAM.globalv.dq;
                ModuleBase::Sphbes::Spherical_Bessel(kkbeta, cell.atoms[it].ncpp.r.data(), q, l, jl);
                for (int ir = 0; ir < kkbeta; ir++)
                {
		            aux[ir] = cell.atoms[it].ncpp.betar(ib, ir) * jl[ir] * cell.atoms[it].ncpp.r[ir];
                }
                double vqint=0.0;
                ModuleBase::Integral::Simpson_Integral(kkbeta, aux, cell.atoms[it].ncpp.rab.data(), vqint);
                this->tab(it, ib, iq) = vqint * pref;
            }
        }
        delete[] aux;
        delete[] jl;
    }
    if (this->use_gpu_)
    {
        if (PARAM.globalv.has_float_data)
        {
            castmem_d2s_h2d_op()(this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2d_op()(this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2d_op()(this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2d_op()(this->s_tab, this->tab.ptr, this->tab.getSize());
            if (this->nhm > 0)
            {
                castmem_d2s_h2d_op()(this->s_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
                castmem_z2c_h2d_op()(this->c_qq_so, this->qq_so.ptr, this->qq_so.getSize());
            }
        }
        if (PARAM.globalv.has_double_data && this->nhm > 0)
        {
            syncmem_z2z_h2d_op()(this->z_qq_so, this->qq_so.ptr, this->qq_so.getSize());
        }
        // Even when the single precision flag is enabled,
        // these variables are utilized in the Force/Stress calculation as well.
        // modified by denghuilu at 2023-05-15
        syncmem_d2d_h2d_op()(this->d_indv, this->indv.c, this->indv.nr * this->indv.nc);
        syncmem_d2d_h2d_op()(this->d_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
        syncmem_d2d_h2d_op()(this->d_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
        syncmem_d2d_h2d_op()(this->d_tab, this->tab.ptr, this->tab.getSize());
        if (this->nhm > 0)
        {
            syncmem_d2d_h2d_op()(this->d_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
        }
    }
    else
    {
        if (PARAM.globalv.has_float_data)
        {
            castmem_d2s_h2h_op()(this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2h_op()(this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2h_op()(this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2h_op()(this->s_tab, this->tab.ptr, this->tab.getSize());
            if (this->nhm > 0)
            {
                castmem_d2s_h2h_op()(this->s_qq_nt, this->qq_nt.ptr, this->qq_nt.getSize());
                castmem_z2c_h2h_op()(this->c_qq_so, this->qq_so.ptr, this->qq_so.getSize());
            }
        }
        // There's no need to synchronize double precision pointers while in a CPU environment.
    }
    ModuleBase::timer::end("ppcell_vnl", "init_vnl");
    GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done." << std::endl;
    return;
}

#ifdef __LCAO
std::complex<double> pseudopot_cell_vnl::Cal_C(int alpha, int lu, int mu, int L, int M) // pengfei Li  2018-3-23
{
    std::complex<double> cf;
    if (alpha == 0)
    {
        cf = -sqrt(4 * ModuleBase::PI / 3) * CG(lu, mu, 1, 1, L, M);
    }
    else if (alpha == 1)
    {
        cf = -sqrt(4 * ModuleBase::PI / 3) * CG(lu, mu, 1, 2, L, M);
    }
    else if (alpha == 2)
    {
        cf = sqrt(4 * ModuleBase::PI / 3) * CG(lu, mu, 1, 0, L, M);
    }
    else
    {
        ModuleBase::WARNING_QUIT("pseudopot_cell_vnl_alpha", "alpha must be 0~2");
    }

    return cf;
}

double pseudopot_cell_vnl::CG(int l1, int m1, int l2, int m2, int L, int M) // pengfei Li 2018-3-23
{
    int dim = L * L + M;
    int dim1 = l1 * l1 + m1;
    int dim2 = l2 * l2 + m2;

    // double A = MGT.Gaunt_Coefficients(dim1, dim2, dim);

    return MGT.Gaunt_Coefficients(dim1, dim2, dim);
}

// void pseudopot_cell_vnl::getvnl_alpha(const int &ik)           // pengfei Li  2018-3-23
// {
// 	if(PARAM.inp.test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl_alpha");
// 	ModuleBase::timer::start("pp_cell_vnl","getvnl_alpha");

// 	if(lmaxkb < 0)
// 	{
// 		return;
// 	}

// 	const int npw = this->wfcpw->npwk[ik];
// 	int ig, ia, nb, ih, lu, mu;

// 	double *vq = new double[npw];
// 	const int x1= (lmaxkb + 2)*(lmaxkb + 2);

// 	ModuleBase::matrix ylm(x1, npw);
// 	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
// 	for (ig = 0;ig < npw;ig++)
// 	{
// 		gk[ig] = this->wfcpw->getgpluskcar(ik,ig);
// 	}

// 	vkb1_alpha = new std::complex<double>**[3];
// 	for(int i=0; i<3; i++)
// 	{
// 		vkb1_alpha[i] = new std::complex<double>*[nhm];
// 		for(int j=0; j<nhm; j++)
// 		{
// 			vkb1_alpha[i][j] = new std::complex<double>[npw];
// 		}
// 	}

// 	vkb_alpha = new std::complex<double>**[3];
// 	for(int i=0; i<3; i++)
// 	{
// 		vkb_alpha[i] = new std::complex<double>*[nkb];
// 		for(int j=0; j<nkb; j++)
// 		{
// 			vkb_alpha[i][j] = new std::complex<double>[this->wfcpw->npwk_max];
// 		}
// 	}

// 	ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

// 	MGT.init_Gaunt_CH( lmaxkb + 2 );
// 	MGT.init_Gaunt( lmaxkb + 2 );

// 	int jkb = 0;
// 	for(int it = 0;it < ucell.ntype;it++)
// 	{
// 		if(PARAM.inp.test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
// 		// calculate beta in G-space using an interpolation table
// 		const int nbeta = ucell.atoms[it].ncpp.nbeta;
// 		const int nh = ucell.atoms[it].ncpp.nh;

// 		if(PARAM.inp.test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

// 		for(int i=0; i<3; i++)
// 			for(int j=0; j<nhm; j++)
// 			{
// 				ModuleBase::GlobalFunc::ZEROS(vkb1_alpha[i][j], npw);
// 			}

// 		for (ih = 0;ih < nh; ih++)
// 		{
// 			lu = static_cast<int>( nhtol(it, ih));
// 			mu = static_cast<int>( nhtolm(it, ih)) - lu * lu;
// 			nb = static_cast<int>( indv(it, ih));

// 			for (int L= std::abs(lu - 1); L<= (lu + 1); L++)
// 			{
// 				for (ig = 0;ig < npw;ig++)
// 				{
// 					const double gnorm = gk[ig].norm() * ucell.tpiba;
// 					vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
// 							this->tab_alpha, it, nb, L, PARAM.globalv.nqx, PARAM.globalv.dq, gnorm);

// 					for (int M=0; M<2*L+1; M++)
// 					{
// 						int lm = L*L + M;
// 						for (int alpha=0; alpha<3; alpha++)
// 						{
// 							std::complex<double> c = Cal_C(alpha,lu, mu, L, M);
// 							/*if(alpha == 0)
// 							{
// 								std::cout<<"lu= "<<lu<<"  mu= "<<mu<<"  L= "<<L<<"  M= "<<M<<" alpha = "<<alpha<<"
// "<<c<<std::endl;
// 							}*/
// 							vkb1_alpha[alpha][ih][ig] += c * vq[ig] * ylm(lm, ig) * pow( ModuleBase::NEG_IMAG_UNIT, L);
// 						}
// 					}
// 				}
// 			}
// 		} // end nbeta

// 		for (ia=0; ia<ucell.atoms[it].na; ia++)
// 		{
// 			std::complex<double> *sk = this->psf->get_sk(ik, it, ia,this->wfcpw);
// 			for (ih = 0;ih < nh;ih++)
// 			{
// 				for (ig = 0;ig < npw;ig++)
// 				{
// 					for(int alpha=0; alpha<3; alpha++)
// 					{
// 						vkb_alpha[alpha][jkb][ig] = vkb1_alpha[alpha][ih][ig] * sk [ig];
// 					}
// 				}
// 				++jkb;
// 			} // end ih
// 			delete [] sk;
// 		} // end ia
// 	} // enddo

// 	delete [] gk;
// 	delete [] vq;
// 	ModuleBase::timer::end("pp_cell_vnl","getvnl_alpha");
// 	return;
// }
#endif

void pseudopot_cell_vnl::init_vnl_alpha(const UnitCell& ucell) // pengfei Li 2018-3-23
{
    if (PARAM.inp.test_pp) {
        ModuleBase::TITLE("pseudopot_cell_vnl", "init_vnl_alpha");
}
    ModuleBase::timer::start("ppcell_vnl", "init_vnl_alpha");

    for (int it = 0; it < ucell.ntype; it++)
    {
        int BetaIndex = 0;
        // const int Nprojectors = ucell.atoms[it].nh;
        for (int ib = 0; ib < ucell.atoms[it].ncpp.nbeta; ib++)
        {
            const int l = ucell.atoms[it].ncpp.lll[ib];
            for (int m = 0; m < 2 * l + 1; m++)
            {
                this->nhtol(it, BetaIndex) = l;
                this->nhtolm(it, BetaIndex) = l * l + m;
                this->indv(it, BetaIndex) = ib;
                ++BetaIndex;
            }
        }
    }

    // max number of beta functions
    const int nbrx = 10;

    const double pref = ModuleBase::FOUR_PI / sqrt(ucell.omega);
    this->tab_alpha.create(ucell.ntype, nbrx, lmaxkb + 2, PARAM.globalv.nqx);
    this->tab_alpha.zero_out();
    GlobalV::ofs_running << "\n Init Non-Local PseudoPotential table( including L index) : ";
    for (int it = 0; it < ucell.ntype; it++)
    {
        const int nbeta = ucell.atoms[it].ncpp.nbeta;
        int kkbeta = ucell.atoms[it].ncpp.kkbeta;

        // mohan modify 2008-3-31
        // mohan add kkbeta>0 2009-2-27
        if ((kkbeta % 2 == 0) && kkbeta > 0)
        {
            kkbeta--;
        }

        double* jl = new double[kkbeta];
        double* aux = new double[kkbeta];

        for (int ib = 0; ib < nbeta; ib++)
        {
            for (int L = 0; L <= lmaxkb + 1; L++)
            {
                for (int iq = 0; iq < PARAM.globalv.nqx; iq++)
                {
                    const double q = iq * PARAM.globalv.dq;
                    ModuleBase::Sphbes::Spherical_Bessel(kkbeta, ucell.atoms[it].ncpp.r.data(), q, L, jl);

                    for (int ir = 0; ir < kkbeta; ir++)
                    {
                        aux[ir] = ucell.atoms[it].ncpp.betar(ib, ir) * jl[ir]
                                  * ucell.atoms[it].ncpp.r[ir] * ucell.atoms[it].ncpp.r[ir];
                    }
                    double vqint = 0.0;
                    ModuleBase::Integral::Simpson_Integral(kkbeta, aux, ucell.atoms[it].ncpp.rab.data(), vqint);
                    this->tab_alpha(it, ib, L, iq) = vqint * pref;
                }
            }
        }
        delete[] aux;
        delete[] jl;
    }
    ModuleBase::timer::end("ppcell_vnl", "init_vnl_alpha");
    GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done(including L)." << std::endl;
    return;
}

void pseudopot_cell_vnl::print_vnl(std::ofstream& ofs)
{
    output::printr3_d(ofs, " tab : ", tab);
}

// ----------------------------------------------------------------------
// scale the non-local pseudopotential tables
void pseudopot_cell_vnl::rescale_vnl(const double& omega_in)
{
    const double ratio = this->omega_old / omega_in;
    const double sqrt_ratio = std::sqrt(ratio);
    this->omega_old = omega_in;

    for (int i = 0; i < this->tab.getSize(); i++)
    {
        this->tab.ptr[i] *= sqrt_ratio;
    }
    for (int i = 0; i < this->tab_at.getSize(); i++)
    {
        this->tab_at.ptr[i] *= sqrt_ratio;
    }
    for (int i = 0; i < this->qrad.getSize(); i++)
    {
        this->qrad.ptr[i] *= ratio;
    }
}

template <>
float* pseudopot_cell_vnl::get_nhtol_data() const
{
    return this->s_nhtol;
}
template <>
double* pseudopot_cell_vnl::get_nhtol_data() const
{
    return this->d_nhtol;
}

template <>
float* pseudopot_cell_vnl::get_nhtolm_data() const
{
    return this->s_nhtolm;
}
template <>
double* pseudopot_cell_vnl::get_nhtolm_data() const
{
    return this->d_nhtolm;
}

template <>
float* pseudopot_cell_vnl::get_indv_data() const
{
    return this->s_indv;
}
template <>
double* pseudopot_cell_vnl::get_indv_data() const
{
    return this->d_indv;
}

template <>
float* pseudopot_cell_vnl::get_tab_data() const
{
    return this->s_tab;
}
template <>
double* pseudopot_cell_vnl::get_tab_data() const
{
    return this->d_tab;
}

template <>
float* pseudopot_cell_vnl::get_deeq_data() const
{
    return this->s_deeq;
}
template <>
double* pseudopot_cell_vnl::get_deeq_data() const
{
    return this->d_deeq;
}

template <>
float* pseudopot_cell_vnl::get_qq_nt_data() const
{
    return this->s_qq_nt;
}
template <>
double* pseudopot_cell_vnl::get_qq_nt_data() const
{
    return this->d_qq_nt;
}

template <>
std::complex<float>* pseudopot_cell_vnl::get_vkb_data() const
{
    return this->c_vkb;
}
template <>
std::complex<double>* pseudopot_cell_vnl::get_vkb_data() const
{
    return this->z_vkb;
}

template <>
std::complex<float>* pseudopot_cell_vnl::get_deeq_nc_data() const
{
    return this->c_deeq_nc;
}
template <>
std::complex<double>* pseudopot_cell_vnl::get_deeq_nc_data() const
{
    return this->z_deeq_nc;
}

template <>
std::complex<float>* pseudopot_cell_vnl::get_qq_so_data() const
{
    return this->c_qq_so;
}
template <>
std::complex<double>* pseudopot_cell_vnl::get_qq_so_data() const
{
    return this->z_qq_so;
}
