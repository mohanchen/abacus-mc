#include "esolver_ks.h"
#include <iostream>

#include "module_base/timer.h"
#include "module_io/input.h"
#include "module_io/json_output/init_info.h"
#include "module_io/print_info.h"
#include "module_parameter/parameter.h"
//--------------Temporary----------------
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//---------------------------------------
#include "module_io/json_output/output_info.h"

namespace ModuleESolver
{


//------------------------------------------------------------------------------
//! the 6th function of ESolver_KS: print_wfcfft
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename T, typename Device>
void ESolver_KS<T, Device>::print_wfcfft(const Input_para& inp, std::ofstream& ofs)
{
    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " | Setup plane waves of wave functions:                            "
           "   |"
        << std::endl;
    ofs << " | Use the energy cutoff and the lattice vectors to generate the   "
           "   |"
        << std::endl;
    ofs << " | dimensions of FFT grid. The number of FFT grid on each "
           "processor   |"
        << std::endl;
    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space "
           "is   |"
        << std::endl;
    ofs << " | different for charege/potential and wave functions. We also set "
           "   |"
        << std::endl;
    ofs << " | the 'sticks' for the parallel of FFT. The number of plane wave "
           "of  |"
        << std::endl;
    ofs << " | each k-point is 'npwk[ik]' in each processor                    "
           "   |"
        << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<"
        << std::endl;
    ofs << "\n\n\n\n";
    ofs << "\n SETUP PLANE WAVES FOR WAVE FUNCTIONS" << std::endl;

    double ecut = inp.ecutwfc;
    if (std::abs(ecut - this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2) > 1e-6)
    {
        ecut = this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2;
        ofs << "Energy cutoff for wavefunc is incompatible with nx, ny, nz and "
               "it will be reduced!"
            << std::endl;
    }
    ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for wavefunc (unit:Ry)", ecut);
    ModuleBase::GlobalFunc::OUT(ofs,
                                "fft grid for wave functions",
                                this->pw_wfc->nx,
                                this->pw_wfc->ny,
                                this->pw_wfc->nz);
    ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", this->pw_wfc->npwtot);
    ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", this->pw_wfc->nstot);

    ofs << "\n PARALLEL PW FOR WAVE FUNCTIONS" << std::endl;
    ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

    for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
    {
        ofs << " " << std::setw(8) << i + 1 << std::setw(15) << this->pw_wfc->nst_per[i] << std::setw(15)
            << this->pw_wfc->npw_per[i] << std::endl;
    }

    ofs << " --------------- sum -------------------" << std::endl;
    ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_wfc->nstot << std::setw(15)
        << this->pw_wfc->npwtot << std::endl;
    ModuleBase::GlobalFunc::DONE(ofs, "INIT PLANEWAVE");
}

}
