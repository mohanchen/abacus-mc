#ifndef DFTU_HAMILT_H
#define DFTU_HAMILT_H

#include <complex>
#include <vector>

class Plus_U_Base;
class Parallel_Orbitals;
class UnitCell;

namespace DFTU_LCAO {

/// DFT+U effective potential in the LCAO basis, with V = pot_onsite (the
/// on-site Hubbard correction potential in the full basis) and S the overlap
/// matrix. Two usage patterns exist below:
///
/// k-space potential: returns the symmetrized U-term potential
///   pot_uterm = (V*S + (V*S)^T) / 2
///
/// real-space Hamiltonian: accumulates the U-term into HR directly
///   HR += (V*SR + SR*V) / 2

/// @brief Compute the LCAO-basis U-term effective potential matrix (complex).
void pot_uterm_complex(Plus_U_Base& dftu,
                       const UnitCell& ucell,
                       const Parallel_Orbitals* pv,
                       const int ik,
                       std::complex<double>* pot_uterm,
                       const std::vector<int>& isk,
                       const std::complex<double>* sk);

/// @brief Compute the LCAO-basis U-term effective potential matrix (real).
void pot_uterm_real(Plus_U_Base& dftu,
                    const UnitCell& ucell,
                    const Parallel_Orbitals* pv,
                    const int ik,
                    double* pot_uterm,
                    const std::vector<int>& isk,
                    const double* sk);

/// @brief Accumulate the DFT+U term into the real-space HR (double).
void pot_uterm_HR_real(const Plus_U_Base& dftu,
                       const UnitCell& ucell,
                       const Parallel_Orbitals* pv,
                       const int ispin,
                       double* SR,
                       double* HR);

/// @brief Accumulate the DFT+U term into the real-space HR (complex).
void pot_uterm_HR_complex(const Plus_U_Base& dftu,
                          const UnitCell& ucell,
                          const Parallel_Orbitals* pv,
                          const int ispin,
                          std::complex<double>* SR,
                          std::complex<double>* HR);

} // namespace DFTU_LCAO

#endif
