#ifndef CHARGE_MATH_H
#define CHARGE_MATH_H

// Free numerical kernels extracted from the Charge class so that the
// density math (summation, electron-count integration, non-linear core
// correction) can be tested and reused without dragging in Charge's state.
// All grid / geometry inputs are passed explicitly instead of being read
// from Charge members or globals.

namespace charge_math
{

// Sum the (spin-resolved) charge density over the real-space grid and
// convert it to a total charge using the cell volume and grid size.
// rho[is][ir] is the density for spin channel is. nspin0 is the number of
// spin channels to include (2 for nspin==2, 1 otherwise).
double sum_rho(double* const* rho,
               const int nspin0,
               const int nrxx,
               const double omega,
               const int nxyz);

// Integrate a single spin channel rho_in over the grid to obtain the
// electron number, scaled by omega / nxyz. Reduction over the pool is
// performed internally under __MPI.
double cal_rho2ne(const double* rho_in,
                  const int nrxx,
                  const double omega,
                  const int nxyz);

// Non-linear core correction: Fourier transform of the (numeric) core
// charge. gg_uniq / ngg supply the reciprocal grid shells previously read
// from Charge::rhopw.
void non_linear_core_correction(const bool numeric,
                                const double omega,
                                const double tpiba2,
                                const int mesh,
                                const double* r,
                                const double* rab,
                                const double* rhoc,
                                double* rhocg,
                                const double* gg_uniq,
                                const int ngg);

} // namespace charge_math

#endif // CHARGE_MATH_H
