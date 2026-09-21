#ifndef CAL_EDM_TDDFT_H
#define CAL_EDM_TDDFT_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_hamilt/hamilt.h"
#include "source_lcao/setup_dm.h"

namespace module_dm
{
void cal_edm_tddft(Parallel_Orbitals& pv,
                   LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                   K_Vectors& kv,
                   hamilt::Hamilt<std::complex<double>>* p_hamilt);

template <typename Device>
void cal_edm_tddft_tensor_lapack(Parallel_Orbitals& pv,
                                 LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                                 K_Vectors& kv,
                                 hamilt::Hamilt<std::complex<double>>* p_hamilt);
} // namespace module_dm
#endif // CAL_EDM_TDDFT_H
