// This file contains subroutines realted to gradient calculations
// it contains 5 subroutines:
// 1. gradcorr, which calculates gradient correction
// 2. grad_wfc, which calculates gradient of wavefunction
//		it is used in stress_func_mgga.cpp
// 3. grad_rho, which calculates gradient of density
// 4. grad_dot, which calculates divergence of something
// 5. noncolin_rho, which diagonalizes the spin density matrix
//  and gives the spin up and spin down components of the charge.

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

void XC_Functional::gradcorr(
    double &etxc,
    double &vtxc,
    ModuleBase::matrix &v,
    const Charge* const chr,
    ModulePW::PW_Basis* rhopw,
    const UnitCell *ucell,
    std::vector<double> &stress_gga,
    const bool is_stress,
    const int nspin,
    const bool domag,
    const bool domag_z,
    const double hybrid_alpha_in,
    const double hse_omega_in)
{
    ModuleBase::TITLE("XC_Functional","gradcorr");

    if((func_type == 3 || func_type == 5) && nspin==4)
    {
        ModuleBase::WARNING_QUIT("gradcorr","meta-GGA has not been implemented for nspin = 4 yet");
    }

    if(func_type == 0 || func_type == 1)
    {
        return;
    }

    // No (semi-)local functional at all. `set_xc_type("HF")` sets func_type = 4 but leaves func_id empty, 
    // Without this block, it will read back whatever the previous functional left in the cleared-but-not-freed
    // buffer (PBE, from set_xc_first_loop), producing a phantom XC energy and potential. (issue deepmodeling/abacus-develop#5404)
    if(func_id.empty())
    {
        // `Stress_Func::stress_gga` only guards on func_type (0/1) and unconditionally reads
        // stress_gga[0..8] afterwards, so hand back an explicit zero tensor rather than an
        // untouched (empty) vector.
        if(is_stress)
        {
            stress_gga.assign(9, 0.0);
        }
        return;
    }

    bool igcc_is_lyp = false;
    // func_id may hold a single entry (e.g. PBE0 -> {XC_HYB_GGA_XC_PBEH}), so guard the index.
    if( func_id.size() > 1 && func_id[1] == XC_GGA_C_LYP)
    {
        igcc_is_lyp = true;
    }

    int nspin0 = nspin;
    if(nspin==4)
    {
        nspin0 =1;
    }
    if(nspin==4&&(domag||domag_z))
    {
        nspin0 = 2;
    }

    assert(nspin0>0);
    const double fac = 1.0/ nspin0;

    // Use cached need_laplacian from set_xc_type
    bool need_laplacian = XC_Functional::get_need_laplacian();

    if(is_stress)
    {
        stress_gga.resize(9);
        for(int i=0;i<9;i++)
        {
            stress_gga[i] = 0.0;
        }
    }

    // sum up (rho_core+rho) for each spin in real space
    // and reciprocal space.
    double* rhotmp1 = nullptr;
    double* rhotmp2 = nullptr;
    std::complex<double>* rhogsum1 = nullptr;
    std::complex<double>* rhogsum2 = nullptr;
    ModuleBase::Vector3<double>* gdr1 = nullptr;
    ModuleBase::Vector3<double>* gdr2 = nullptr;
    ModuleBase::Vector3<double>* h1 = nullptr;
    ModuleBase::Vector3<double>* h2 = nullptr;
    double* neg = nullptr;
    double** vsave = nullptr;
    double** vgg = nullptr;
    std::vector<double> lapl1;
    std::vector<double> lapl2;
    std::vector<double> vlapl_arr1;
    std::vector<double> vlapl_arr2;

    gradcorr_prepare_rho(chr, rhopw, ucell, nspin, nspin0, fac,
        need_laplacian, is_stress, domag, domag_z, v, rhotmp1, rhotmp2,
        rhogsum1, rhogsum2, gdr1, gdr2, h1, h2, neg,
        vsave, vgg, lapl1, lapl2, vlapl_arr1, vlapl_arr2);

    double vtxcgc = 0.0;
    double etxcgc = 0.0;

    gradcorr_xc_kernel(chr, rhopw, nspin, nspin0, fac,
        need_laplacian, is_stress, igcc_is_lyp, domag, domag_z,
        hybrid_alpha_in, hse_omega_in,
        rhotmp1, rhotmp2, gdr1, gdr2, lapl1, lapl2, neg,
        vtxcgc, etxcgc, stress_gga, h1, h2,
        vlapl_arr1, vlapl_arr2, v);

    // Add Laplacian stress contribution from meta-GGA functionals
    if(is_stress && use_libxc && !vlapl_arr1.empty())
    {
        const int ng = rhopw->npw;
        const int nrxx = rhopw->nrxx;
        const double tpiba2 = ucell->tpiba * ucell->tpiba;
        std::vector<std::complex<double>> vlapl_g(rhopw->nmaxgr);

        for(int is = 0; is < nspin0; is++)
        {
            double* vlapl_ptr = (is == 0) ? vlapl_arr1.data() : vlapl_arr2.data();
            double* rho_ptr = (is == 0) ? rhotmp1 : rhotmp2;
            if(vlapl_ptr == nullptr || rho_ptr == nullptr) continue;

            std::vector<std::complex<double>> rho_g(rhopw->nmaxgr);
            for(int ir = 0; ir < nrxx; ir++)
                rho_g[ir] = std::complex<double>(rho_ptr[ir], 0.0);
            rhopw->real2recip(rho_g.data(), rho_g.data());

            for(int ir = 0; ir < nrxx; ir++)
                vlapl_g[ir] = std::complex<double>(vlapl_ptr[ir], 0.0);
            rhopw->real2recip(vlapl_g.data(), vlapl_g.data());

            for(int l = 0; l < 3; l++)
            {
                for(int m = 0; m <= l; m++)
                {
                    double sum = 0.0;
                    for(int ig = 0; ig < ng; ig++)
                    {
                        double g_prod = rhopw->gcar[ig][l] * rhopw->gcar[ig][m] * tpiba2;
                        sum += g_prod * (rho_g[ig].real() * vlapl_g[ig].real()
                                       + rho_g[ig].imag() * vlapl_g[ig].imag());
                    }
                    stress_gga[l*3+m] -= 2.0 * static_cast<double>(rhopw->nxyz) * ModuleBase::e2 * sum;
                }
            }
        }
    }

    if(!is_stress)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for(int ir=0; ir<rhopw->nrxx; ir++)
        {
            rhotmp1[ir] -= fac * chr->rho_core[ir];
        }
        if(nspin0==2)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir=0; ir<rhopw->nrxx; ir++)
            {
                rhotmp2[ir] -= fac * chr->rho_core[ir];
            }
        }

        // second term of the gradient correction :
        // sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )

        // dh is in real sapce.
        double* dh = new double[rhopw->nrxx];

        for(int is=0; is<nspin0; is++)
        {
            if(is==0)
            {
                XC_Functional::grad_dot(h1,dh,rhopw,ucell->tpiba);
            }
            if(is==1)
            {
                XC_Functional::grad_dot(h2,dh,rhopw,ucell->tpiba);
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir=0; ir<rhopw->nrxx; ir++)
            {
                v(is, ir) -= dh[ir];
            }

            double sum = 0.0;
            if(is==0)
            {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) schedule(static, 256)
#endif
                for(int ir=0; ir<rhopw->nrxx; ir++)
                {
                    sum += dh[ir] * rhotmp1[ir];
                }
            }
            else if(is==1)
            {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) schedule(static, 256)
#endif
                for(int ir=0; ir<rhopw->nrxx; ir++)
                {
                    sum += dh[ir] * rhotmp2[ir];
                }
            }
            vtxcgc -= sum;
        }

        delete[] dh;

        vtxc += vtxcgc;
        etxc += etxcgc;

        // Add Laplacian contribution from meta-GGA functionals
        // v_xc += nabla^2(vlapl) where vlapl = d(rho*eps_xc)/d(nabla^2 rho)
        if(use_libxc && !vlapl_arr1.empty())
        {
            const int ng = rhopw->npw;
            const int nrxx = rhopw->nrxx;
            const double tpiba2 = ucell->tpiba * ucell->tpiba;
            std::vector<std::complex<double>> vlapl_g(rhopw->nmaxgr);

            for(int is = 0; is < nspin0; is++)
            {
                double* vlapl_ptr = (is == 0) ? vlapl_arr1.data() : vlapl_arr2.data();
                if(vlapl_ptr == nullptr) continue;

                for(int ir = 0; ir < nrxx; ir++)
                    vlapl_g[ir] = std::complex<double>(vlapl_ptr[ir], 0.0);
                rhopw->real2recip(vlapl_g.data(), vlapl_g.data());
                for(int ig = 0; ig < ng; ig++)
                {
                    double g2 = 0.0;
                    for(int i = 0; i < 3; i++)
                        g2 += rhopw->gcar[ig][i] * rhopw->gcar[ig][i];
                    vlapl_g[ig] *= -g2 * tpiba2;
                }
                rhopw->recip2real(vlapl_g.data(), vlapl_g.data());
                for(int ir = 0; ir < nrxx; ir++)
                {
                    v(is, ir) += ModuleBase::e2 * vlapl_g[ir].real();
                }
            }
        }

        if(nspin == 4 && (domag||domag_z))
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
            for(int is=0;is<nspin;is++)
            {
                for(int ir=0;ir<rhopw->nrxx;ir++)
                {
                    if(is<nspin0)
                    {
                        vgg[is][ir] = v(is,ir);
                    }
                    v(is,ir) = vsave[is][ir];
                }
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
            for(int ir=0;ir<rhopw->nrxx;ir++)
            {
                v(0,ir) += 0.5 * (vgg[0][ir] + vgg[1][ir]);
                double amag = sqrt(pow(chr->rho[1][ir],2)+pow(chr->rho[2][ir],2)+pow(chr->rho[3][ir],2));
                if(amag>1e-12)
                {
                    for(int i=1;i<4;i++)
                    {
                        v(i,ir)+= neg[ir] * 0.5 *(vgg[0][ir]-vgg[1][ir])*chr->rho[i][ir]/amag;
                    }
                }
            }
        }
    }
    // deacllocate
    delete[] rhotmp1;
    delete[] rhogsum1;
    delete[] gdr1;
    if(!is_stress)
    {
        delete[] h1;
    }

    if(nspin==2)
    {
        delete[] rhotmp2;
        delete[] rhogsum2;
        delete[] gdr2;
        if(!is_stress)
        {
            delete[] h2;
        }
    }
    if(nspin == 4 && (domag||domag_z))
    {
        delete[] neg;
        if(!is_stress)
        {
            for(int i=0; i<nspin0; i++)
            {
                delete[] vgg[i];
            }
            delete[] vgg;
            for(int i=0; i<nspin; i++)
            {
                delete[] vsave[i];
            }
            delete[] vsave;
            delete[] h2;
        }
        delete[] rhotmp2;
        delete[] rhogsum2;
        delete[] gdr2;
    }

    return;
}

