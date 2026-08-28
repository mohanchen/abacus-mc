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

    gradcorr_assemble_vxc(chr, rhopw, ucell, nspin, nspin0, fac,
        is_stress, domag, domag_z, vtxcgc, etxcgc,
        vtxc, etxc, stress_gga, v,
        rhotmp1, rhotmp2, h1, h2,
        vlapl_arr1, vlapl_arr2, neg, vsave, vgg);
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

