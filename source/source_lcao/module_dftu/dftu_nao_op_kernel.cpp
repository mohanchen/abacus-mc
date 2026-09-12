/// @file dftu_nao_op_kernel.cpp
/// @brief Per-pair kernel functions for the DFTU LCAO operator, split from
///        dftu_nao_op.cpp to keep the main file under the 500-line limit.
///        These are template member-function definitions that are explicitly
///        instantiated at the end of this translation unit.
#include "dftu_nao_op.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_base/parallel_reduce.h"

// Include the free function implementations for force/stress in real space
#include "dftu_nao_fs_r.h"

// cal_HR_IJR: accumulate HR += <psi_I|beta_m> * pot_onsite(m,m') * <beta_m'|psi_{J,R}>
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(
    const int& iat1,
    const int& iat2,
    const Parallel_Orbitals* pv,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<TR>& pot_onsite,
    TR* data_pointer)
{
    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    auto row_indexes = pv->get_indexes_row(iat1);
    auto col_indexes = pv->get_indexes_col(iat2);
    const int m_size = int(sqrt(pot_onsite.size()) / npol);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = pv->get_ncol_atom(iat2) * is + is2;
        }
    }
    // calculate the local matrix
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is = 0; is < npol * npol; ++is)
            {
                int start = is * m_size * m_size;
                TR nlm_tmp = TR(0);
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for (int m2 = 0; m2 < m_size; m2++)
                    {
                        nlm_tmp += nlm1[m1] * nlm2[m2] * pot_onsite[m1 * m_size + m2 + start];
                    }
                }
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

// cal_occ: compute occ_mm' = sum_R DMR*<phi_0|alpha^I_m'><alpha^I_m'|phi_R>
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_occ(const int& iat1,
                                                         const int& iat2,
                                                         const Parallel_Orbitals* pv,
                                                         const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                                         const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                                         const double* dm_pointer,
                                                         std::vector<double>& occ)
{
    const int npol = this->ucell->get_npol();
    auto row_indexes = pv->get_indexes_row(iat1);
    auto col_indexes = pv->get_indexes_col(iat2);
    const int m_size = int(sqrt(occ.size()) / npol);
    const int m_size2 = m_size * m_size;
#ifdef __DEBUG
    assert(m_size * m_size == occ.size());
#endif
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = pv->get_ncol_atom(iat2) * is + is2;
        }
    }
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is1 = 0; is1 < npol; ++is1)
            {
                for (int is2 = 0; is2 < npol; ++is2)
                {
                    for (int m1 = 0; m1 < m_size; ++m1)
                    {
                        for (int m2 = 0; m2 < m_size; ++m2)
                        {
                            occ[m1 * m_size + m2 + (is1 * npol + is2) * m_size2]
                                += nlm1[m1] * nlm2[m2] * dm_pointer[step_trace[is1 * npol + is2]];
                        }
                    }
                }
            }
            dm_pointer += npol;
        }
        dm_pointer += (npol - 1) * col_indexes.size();
    }
}

// transfer_pot_onsite (generic: identity copy)
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::transfer_pot_onsite(std::vector<double>& pot_onsite_tmp, std::vector<TR>& pot_onsite)
{
#ifdef __DEBUG
    assert(pot_onsite.size() == pot_onsite_tmp.size());
#endif
    for (int i = 0; i < pot_onsite_tmp.size(); i++)
    {
        pot_onsite[i] = pot_onsite_tmp[i];
    }
}

// transfer_pot_onsite (noncollinear specialization: Pauli-to-spinor)
template <>
void hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>::transfer_pot_onsite(
    std::vector<double>& pot_onsite_tmp,
    std::vector<std::complex<double>>& pot_onsite)
{
#ifdef __DEBUG
    assert(pot_onsite.size() == pot_onsite_tmp.size());
#endif

    // Pauli-to-spinor conversion for DFT+U potential:
    // V = V_0*I + V_x*sigma_x + V_y*sigma_y + V_z*sigma_z
    const int m_size = int(sqrt(pot_onsite.size()) / 2);
    const int m_size2 = m_size * m_size;
    pot_onsite.resize(pot_onsite_tmp.size());
    for (int m1 = 0; m1 < m_size; m1++)
    {
        for (int m2 = 0; m2 < m_size; m2++)
        {
            int index[4];
            index[0] = m1 * m_size + m2;
            index[1] = m1 * m_size + m2 + m_size2;
            index[2] = m2 * m_size + m1 + m_size2 * 2;
            index[3] = m2 * m_size + m1 + m_size2 * 3;
            pot_onsite[index[0]] = 0.5 * (pot_onsite_tmp[index[0]] + pot_onsite_tmp[index[3]]);
            pot_onsite[index[3]] = 0.5 * (pot_onsite_tmp[index[0]] - pot_onsite_tmp[index[3]]);
            pot_onsite[index[1]] = 0.5 * (pot_onsite_tmp[index[1]] - std::complex<double>(0.0, 1.0) * pot_onsite_tmp[index[2]]);
            pot_onsite[index[2]] = 0.5 * (pot_onsite_tmp[index[1]] + std::complex<double>(0.0, 1.0) * pot_onsite_tmp[index[2]]);
        }
    }
}

// cal_pot_onsite: pot = U * (1/2*delta - occ), energy = U * 1/2 * occ * occ
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_pot_onsite(const std::vector<double>& occ,
                                                            const int m_size,
                                                            const double u_value,
                                                            double* pot_onsite,
                                                            double& eu)
{
    int spin_fold = occ.size() / m_size / m_size;
    if (spin_fold < 4) {
        for (int is = 0; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0.5 * (m1 == m2) - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.5 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    } else
    {
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for (int m2 = 0; m2 < m_size; m2++)
            {
                pot_onsite[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - occ[m2 * m_size + m1]);
                eu += u_value * 0.25 * occ[m2 * m_size + m1] * occ[m1 * m_size + m2];
            }
        }
        for (int is = 1; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0 - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.25 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    }
}

// explicit template instantiation (matches dftu_nao_op.cpp)
template class hamilt::DFTU<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
