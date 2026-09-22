#include "density_matrix.h"

#include "source_base/libm/libm.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"

namespace module_dm
{

template <>
void DensityMatrix_Tools::func_exp_mul_dmk<double>(
    const std::complex<double> kphase,
    const std::vector<std::complex<double>>& DMK_mat_trans,
    double* target_DMR_mat)
{
    const std::size_t mat_size = DMK_mat_trans.size();
    for(std::size_t i = 0; i < mat_size; i++)
    {
        target_DMR_mat[i]
            += kphase.real() * DMK_mat_trans[i].real()
            - kphase.imag() * DMK_mat_trans[i].imag();
    }
}

template <>
void DensityMatrix_Tools::func_exp_mul_dmk<std::complex<double>>(
    const std::complex<double> kphase,
    const std::vector<std::complex<double>>& DMK_mat_trans,
    std::complex<double>* target_DMR_mat)
{
    BlasConnector::axpy(DMK_mat_trans.size(),
                        kphase,
                        DMK_mat_trans.data(),
                        1,
                        target_DMR_mat,
                        1);
}

template <>
void DensityMatrix_Tools::func_xyz_to_updown<double>(
    const std::complex<double> tmp[4],
    const int icol,
    const int step_trace[4],
    double* target_DMR_mat)
{
    target_DMR_mat[icol + step_trace[0]] = tmp[0].real() + tmp[3].real();  // rho_0 = (rho_upup + rho_downdown).real()
    target_DMR_mat[icol + step_trace[1]] = tmp[1].real() + tmp[2].real();  // rho_x = (rho_updown + rho_downup).real()
    // rho_y: the stored DM block is the complex conjugate of the physical 1-RDM P (cal_dm_psi builds
    // DM_{ab}=sum conj(c_a) c_b = conj(P), so tmp[1]=DM_{ud}=conj(P_{ud})). Extracting m_y from the
    // CONJUGATED block therefore carries the opposite sign of the bare-textbook formula; m_x/m_z read
    // Re() and are conjugation-invariant. Using the bare formula (PR #7664) sign-flips m_y and quenches
    // in-plane non-collinear moments (e.g. Mn3Sn 120-deg AFM); see issue #7831.
    target_DMR_mat[icol + step_trace[2]] = tmp[1].imag() - tmp[2].imag();  // rho_y = Im(P_updown) - Im(P_downup)
    target_DMR_mat[icol + step_trace[3]] = tmp[0].real() - tmp[3].real();  // rho_z = (rho_upup - rho_downdown).real()
}

template <>
void DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(
    const std::complex<double> tmp[4],
    const int icol,
    const int step_trace[4],
    std::complex<double>* target_DMR_mat)
{
    target_DMR_mat[icol + step_trace[0]] = tmp[0] + tmp[3];  // rho_0 = (rho_upup + rho_downdown)
    target_DMR_mat[icol + step_trace[1]] = tmp[1] + tmp[2];  // rho_x = (rho_updown + rho_downup)
    // rho_y sign accounts for the conjugated stored DM block (conj(P)); see the <double> specialization above.
    target_DMR_mat[icol + step_trace[2]]
        = -ModuleBase::IMAG_UNIT * (tmp[1] - tmp[2]);  // rho_y = -i*(rho_updown - rho_downup)
    target_DMR_mat[icol + step_trace[3]] = tmp[0] - tmp[3];  // rho_z = (rho_upup - rho_downdown)
}

} // namespace module_dm
