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

/// compute VU and energy contribution for one atom (nspin==4, spinor).
///
/// Writes 4 Pauli blocks of VU into vu (size 4 * m_size * m_size) and
/// returns the energy_u increment. Internally calls pauli_to_spin_basis
/// to convert vu to spin basis in-place.
///
/// vu:  pointer to eff_pot_pw[eff_pot_pw_index[iat]]
/// occ: pointer to occ_mat[iat][target_l][0][0].c (4 Pauli blocks packed)
double compute_vu_spinor(
    std::complex<double>* vu,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size);

/// compute VU and energy contribution for one atom, one spin channel
/// (nspin==1 or nspin==2). Returns the energy_u increment.
///
/// vu:  pointer to the spin channel's VU block (size m_size * m_size)
/// occ: pointer to occ_mat[iat][target_l][0][is].c for this channel
double compute_vu_scalar(
    std::complex<double>* vu,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size);

} // namespace dftu_pw

#endif
