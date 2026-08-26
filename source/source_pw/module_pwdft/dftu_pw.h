#ifndef DFTU_PW_H
#define DFTU_PW_H

#include <complex>

/// Free functions for DFT+U PW basis calculations.
///
/// These functions are pure (no access to Plus_U_Base members) so they can be
/// unit-tested directly by including this header. The member functions in
/// dftu_pw.cpp call them after computing per-atom offsets and fetching the
/// relevant member state (occ_mat, eff_pot_pw, u_current, etc.).
namespace dftu_pw {

/// transform VU from Pauli basis to spin basis (in-place, nspin==4 only).
///
/// vu points to the per-atom VU block of size 4 * m_size * m_size, storing
/// 4 contiguous Pauli blocks (is=0 charge, is=1 sigma_x, is=2 sigma_y,
/// is=3 sigma_z). After the transform, the same memory holds the spin
/// representation:
///   vu[0]        <- 0.5 * (vu_pauli[0] + vu_pauli[3])
///   vu[3*size]   <- 0.5 * (vu_pauli[0] - vu_pauli[3])
///   vu[size]     <- 0.5 * (vu_pauli[1] + i * vu_pauli[2])
///   vu[2*size]   <- 0.5 * (vu_pauli[1] - i * vu_pauli[2])
void pauli_to_spin_basis(std::complex<double>* vu, int m_size);

} // namespace dftu_pw

#endif
