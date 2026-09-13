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
/// matrix. Returns the symmetrized k-space U-term potential:
///   pot_uterm = (V*S + (V*S)^T) / 2

/// @brief Compute the LCAO-basis U-term effective potential matrix.
///
/// @tparam T           matrix element type (double or std::complex<double>)
/// @param dftu         Plus_U state
/// @param ucell        unit cell
/// @param pv           parallel orbitals descriptor
/// @param spin         spin channel (isk[ik] from caller)
/// @param pot_uterm    output buffer (length pv->nloc)
/// @param sk           overlap matrix in k-space (length pv->nloc)
template <typename T>
void cal_pot_uterm(Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   const Parallel_Orbitals* pv,
                   const int spin,
                   T* pot_uterm,
                   const T* sk);

// Explicit instantiations
extern template void cal_pot_uterm<double>(Plus_U_Base& dftu,
                                           const UnitCell& ucell,
                                           const Parallel_Orbitals* pv,
                                           const int spin,
                                           double* pot_uterm,
                                           const double* sk);
extern template void cal_pot_uterm<std::complex<double>>(Plus_U_Base& dftu,
                                                         const UnitCell& ucell,
                                                         const Parallel_Orbitals* pv,
                                                         const int spin,
                                                         std::complex<double>* pot_uterm,
                                                         const std::complex<double>* sk);

} // namespace DFTU_LCAO

#endif
