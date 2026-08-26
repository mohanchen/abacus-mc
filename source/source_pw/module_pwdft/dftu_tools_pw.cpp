#include "source_pw/module_pwdft/dftu_tools_pw.h"

namespace dftu_pw {

void pauli_to_spin_basis(std::complex<double>* vu, int m_size)
{
    const int size = m_size * m_size;
    for (int m1 = 0; m1 < m_size; m1++)
    {
        for (int m2 = 0; m2 < m_size; m2++)
        {
            int index[4];
            index[0] = m1 * m_size + m2;
            index[1] = m1 * m_size + m2 + size;
            index[2] = m1 * m_size + m2 + size * 2;
            index[3] = m1 * m_size + m2 + size * 3;
            std::complex<double> vu_tmp[4];
            for (int i = 0; i < 4; i++)
            {
                vu_tmp[i] = vu[index[i]];
            }
            vu[index[0]] = 0.5 * (vu_tmp[0] + vu_tmp[3]);
            vu[index[3]] = 0.5 * (vu_tmp[0] - vu_tmp[3]);
            vu[index[1]] = 0.5 * (vu_tmp[1] + std::complex<double>(0.0, 1.0) * vu_tmp[2]);
            vu[index[2]] = 0.5 * (vu_tmp[1] - std::complex<double>(0.0, 1.0) * vu_tmp[2]);
        }
    }
}

double compute_vu_spinor(
    std::complex<double>* vu,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size)
{
    double energy_u = 0.0;
    const int m_size2 = m_size * m_size;
    for (int is = 0; is < 4; ++is)
    {
        int start = is * m_size2;
        double diag = (is == 0) ? diag_coeff : 0.0;
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for (int m2 = 0; m2 < m_size; m2++)
            {
                vu[start + m1 * m_size + m2] = u_value *
                    (diag * (m1 == m2) - occ[start + m2 * m_size + m1]);
                energy_u += u_value * weight_eu
                    * occ[start + m2 * m_size + m1]
                    * occ[start + m1 * m_size + m2];
            }
        }
    }
    pauli_to_spin_basis(vu, m_size);
    return energy_u;
}

double compute_vu_scalar(
    std::complex<double>* vu,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size)
{
    double energy_u = 0.0;
    for (int m1 = 0; m1 < m_size; m1++)
    {
        for (int m2 = 0; m2 < m_size; m2++)
        {
            vu[m1 * m_size + m2] = u_value *
                (diag_coeff * (m1 == m2) - occ[m2 * m_size + m1]);
            energy_u += u_value * weight_eu
                * occ[m2 * m_size + m1]
                * occ[m1 * m_size + m2];
        }
    }
    return energy_u;
}

void accumulate_occ_spinor(
    double* occ_mat_out,
    const std::complex<double>* becp,
    int nbands,
    int npol,
    int nkb,
    int begin_ih,
    int m_begin,
    int tlp1,
    const ModuleBase::matrix& wg,
    int ik)
{
    const int tlp1_2 = tlp1 * tlp1;
    for (int ib = 0; ib < nbands; ib++)
    {
        const double weight = wg(ik, ib);
        int ind_m1m2 = 0;
        for (int m1 = 0; m1 < tlp1; m1++)
        {
            const int index_m1 = ib * npol * nkb + begin_ih + m_begin + m1;
            for (int m2 = 0; m2 < tlp1; m2++)
            {
                const int index_m2 = ib * npol * nkb + begin_ih + m_begin + m2;
                std::complex<double> occ[4];
                occ[0] = weight * std::conj(becp[index_m1]) * becp[index_m2];
                occ[1] = weight * std::conj(becp[index_m1]) * becp[index_m2 + nkb];
                occ[2] = weight * std::conj(becp[index_m1 + nkb]) * becp[index_m2];
                occ[3] = weight * std::conj(becp[index_m1 + nkb]) * becp[index_m2 + nkb];
                occ_mat_out[ind_m1m2] += (occ[0] + occ[3]).real();
                occ_mat_out[ind_m1m2 + tlp1_2] += (occ[1] + occ[2]).real();
                occ_mat_out[ind_m1m2 + 2 * tlp1_2] += (occ[1] - occ[2]).imag();
                occ_mat_out[ind_m1m2 + 3 * tlp1_2] += (occ[0] - occ[3]).real();
                ind_m1m2++;
            }
        }
    }
}

void accumulate_occ_scalar(
    double* occ_mat_out,
    const std::complex<double>* becp,
    int nbands,
    int nkb,
    int begin_ih,
    int m_begin,
    int tlp1,
    const ModuleBase::matrix& wg,
    int ik)
{
    for (int ib = 0; ib < nbands; ib++)
    {
        const double weight = wg(ik, ib);
        int ind_m1m2 = 0;
        for (int m1 = 0; m1 < tlp1; m1++)
        {
            const int index_m1 = ib * nkb + begin_ih + m_begin + m1;
            for (int m2 = 0; m2 < tlp1; m2++)
            {
                const int index_m2 = ib * nkb + begin_ih + m_begin + m2;
                occ_mat_out[ind_m1m2] += weight * (std::conj(becp[index_m1]) * becp[index_m2]).real();
                ind_m1m2++;
            }
        }
    }
}

} // namespace dftu_pw
