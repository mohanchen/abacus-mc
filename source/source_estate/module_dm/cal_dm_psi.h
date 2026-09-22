#ifndef CAL_DM_PSI_H
#define CAL_DM_PSI_H

#include "density_matrix.h"
#include "source_base/matrix.h"
#include "source_psi/psi.h"

namespace module_dm
{
// for Gamma-Only case where DMK is double
void cal_dm_psi(const Parallel_Orbitals* ParaV,
                const ModuleBase::matrix& wg,
                const psi::Psi<double>& wfc,
                module_dm::DensityMatrix<double, double>& DM);

// for Multi-k case where DMK is std::complex<double>
template <typename TR>
void cal_dm_psi(const Parallel_Orbitals* ParaV,
                const ModuleBase::matrix& wg,
                const psi::Psi<std::complex<double>>& wfc,
                module_dm::DensityMatrix<std::complex<double>, TR>& DM);

#ifdef __MPI
// for Gamma-Only case with MPI
void psi2dm_mpi(const psi::Psi<double>& psi1,
                const psi::Psi<double>& psi2,
                double* dm_out,
                const int* desc_psi,
                const int* desc_dm);

// for multi-k case with MPI
void psi2dm_mpi(const psi::Psi<std::complex<double>>& psi1,
                const psi::Psi<std::complex<double>>& psi2,
                std::complex<double>* dm_out,
                const int* desc_psi,
                const int* desc_dm);

#else
// for Gamma-Only case without MPI
void psi2dm(const psi::Psi<double>& psi1, const psi::Psi<double>& psi2, double* dm_out);

// for multi-k case without MPI
void psi2dm(const psi::Psi<std::complex<double>>& psi1,
            const psi::Psi<std::complex<double>>& psi2,
            std::complex<double>* dm_out);
#endif
} // namespace module_dm
#endif
