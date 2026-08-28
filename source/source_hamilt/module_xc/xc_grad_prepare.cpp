#include "xc_functional.h"
#include "source_base/timer.h"
#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_io/module_parameter/parameter.h"
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_types.h>
#include <source_hamilt/module_xc/kernels/xc_functional_op.h>

#ifdef __LIBXC
#include "libxc_abacus.h"
#ifdef __EXX
#include "source_hamilt/module_xc/exx_info.h"
#endif
#endif

void XC_Functional::gradcorr_prepare_rho(
    const Charge* const chr,
    ModulePW::PW_Basis* rhopw,
    const UnitCell* ucell,
    const int nspin,
    const int nspin0,
    const double fac,
    const bool need_laplacian,
    const bool is_stress,
    const bool domag,
    const bool domag_z,
    ModuleBase::matrix& v,
    double*& rhotmp1,
    double*& rhotmp2,
    std::complex<double>*& rhogsum1,
    std::complex<double>*& rhogsum2,
    ModuleBase::Vector3<double>*& gdr1,
    ModuleBase::Vector3<double>*& gdr2,
    ModuleBase::Vector3<double>*& h1,
    ModuleBase::Vector3<double>*& h2,
    double*& neg,
    double**& vsave,
    double**& vgg,
    std::vector<double>& lapl1,
    std::vector<double>& lapl2,
    std::vector<double>& vlapl_arr1,
    std::vector<double>& vlapl_arr2)
{
    // doing FFT to get rho in G space: rhog1
    rhopw->real2recip(chr->rho[0], chr->rhog[0]);
    if(nspin==2)
    {
        rhopw->real2recip(chr->rho[1], chr->rhog[1]);
    }
    rhopw->real2recip(chr->rho_core, chr->rhog_core);

    // sum up (rho_core+rho) for each spin in real space
    // and reciprocal space.

    // for spin unpolarized case,
    // calculate the gradient of (rho_core+rho) in reciprocal space.
    rhotmp1 = new double[rhopw->nrxx];
    rhogsum1 = new std::complex<double>[rhopw->npw];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ir=0; ir<rhopw->nrxx; ir++)
    {
        rhotmp1[ir] = chr->rho[0][ir] + fac * chr->rho_core[ir];
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ig=0; ig<rhopw->npw; ig++)
    {
        rhogsum1[ig] = chr->rhog[0][ig] + fac * chr->rhog_core[ig];
    }

    gdr1 = new ModuleBase::Vector3<double>[rhopw->nrxx];
    if(!is_stress)
    {
        h1 = new ModuleBase::Vector3<double>[rhopw->nrxx];
    }

    XC_Functional::grad_rho( rhogsum1 , gdr1, rhopw, ucell->tpiba);

    if(need_laplacian)
    {
        lapl1.resize(rhopw->nrxx);
        XC_Functional::laplacian_rho(rhogsum1, lapl1.data(), rhopw, ucell->tpiba);
        if(use_libxc)
        {
            vlapl_arr1.resize(rhopw->nrxx, 0.0);
        }
    }

    // for spin polarized case;
    // calculate the gradient of (rho_core+rho) in reciprocal space.
    if(nspin==2)
    {
        rhotmp2 = new double[rhopw->nrxx];
        rhogsum2 = new std::complex<double>[rhopw->npw];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            rhotmp2[ir] = chr->rho[1][ir] + fac * chr->rho_core[ir];
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rhopw->npw; ig++)
        {
            rhogsum2[ig] = chr->rhog[1][ig] + fac * chr->rhog_core[ig];
        }

        gdr2 = new ModuleBase::Vector3<double>[rhopw->nrxx];
        if(!is_stress)
        {
            h2 = new ModuleBase::Vector3<double>[rhopw->nrxx];
        }

        XC_Functional::grad_rho( rhogsum2 , gdr2, rhopw, ucell->tpiba);

        if(need_laplacian)
        {
            lapl2.resize(rhopw->nrxx);
            XC_Functional::laplacian_rho(rhogsum2, lapl2.data(), rhopw, ucell->tpiba);
            if(use_libxc)
            {
                vlapl_arr2.resize(rhopw->nrxx, 0.0);
            }
        }
    }

    if(nspin == 4&&(domag||domag_z))
    {
        rhotmp2 = new double[rhopw->nrxx];
        rhogsum2 = new std::complex<double>[rhopw->npw];
        neg = new double [rhopw->nrxx];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            rhotmp1[ir] = 0.0;
            rhotmp2[ir] = 0.0;
            neg[ir] = 0.0;
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rhopw->npw; ig++)
        {
            rhogsum1[ig] = 0.0;
            rhogsum2[ig] = 0.0;
        }
        if(!is_stress)
        {
            vsave = new double* [nspin];
            for(int is = 0;is<nspin;is++)
            {
                vsave[is]= new double [rhopw->nrxx];
            }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int is = 0;is<nspin;is++)
            {
                for(int ir =0;ir<rhopw->nrxx;ir++)
                {
                    vsave[is][ir] = v(is,ir);
                    v(is,ir) = 0;
                }
            }
            vgg = new double* [nspin0];
            for(int is = 0;is<nspin0;is++)
            {
                vgg[is] = new double[rhopw->nrxx];
            }
        }
        noncolin_rho(rhotmp1, rhotmp2, neg, chr->rho, rhopw->nrxx, ucell->magnet.ux_, ucell->magnet.lsign_);
        rhopw->real2recip(rhotmp1, rhogsum1);
        rhopw->real2recip(rhotmp2, rhogsum2);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            rhotmp2[ir] += fac * chr->rho_core[ir];
            rhotmp1[ir] += fac * chr->rho_core[ir];
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ig=0; ig<rhopw->npw; ig++)
        {
            rhogsum2[ig] += fac * chr->rhog_core[ig];
            rhogsum1[ig] += fac * chr->rhog_core[ig];
        }

        gdr2 = new ModuleBase::Vector3<double>[rhopw->nrxx];
        if(!is_stress)
        {
            h2 = new ModuleBase::Vector3<double>[rhopw->nrxx];
        }

        XC_Functional::grad_rho( rhogsum1 , gdr1, rhopw, ucell->tpiba);
        XC_Functional::grad_rho( rhogsum2 , gdr2, rhopw, ucell->tpiba);

        if(need_laplacian)
        {
            lapl1.resize(rhopw->nrxx);
            XC_Functional::laplacian_rho(rhogsum1, lapl1.data(), rhopw, ucell->tpiba);
            lapl2.resize(rhopw->nrxx);
            XC_Functional::laplacian_rho(rhogsum2, lapl2.data(), rhopw, ucell->tpiba);
            if(use_libxc)
            {
                vlapl_arr1.resize(rhopw->nrxx, 0.0);
                vlapl_arr2.resize(rhopw->nrxx, 0.0);
            }
        }
    }
}
