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

// The real-space atom-pair kernels cal_hr_ijr<TR>() and cal_occ_ijr() live
// in dftu_nao_ijr.h as free functions of namespace DFTU_LCAO.

namespace DFTU_LCAO
{

/// On-site potential and Hubbard energy for one correlated shell:
///   pot_onsite(m,m') = U_eff * (0.5 * delta_{m,m'} - occ(m,m'))
///   EU = (U_eff / 2) * sum_{m,m'} occ(m,m') * (delta_{m,m'} - occ(m',m))
void cal_pot_onsite(const std::vector<double>& occ, const int m_size, const double u_value,
                    double* pot_onsite, double& eu)
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

} // namespace DFTU_LCAO

// explicit template instantiation (matches dftu_nao_op.cpp)
template class hamilt::DFTU<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
