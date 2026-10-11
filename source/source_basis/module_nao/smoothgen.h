#ifndef SMOOTHGEN_H_
#define SMOOTHGEN_H_

#include <vector>

/**
 * @brief Smoothly truncates a radial function at a cutoff radius.
 *
 * Starting with a radial function chi on grid r, this function generates a
 * new radial function alpha on the truncated grid (points below rcut) such
 * that alpha decays smoothly to zero at rcut. The smoothness parameter sigma
 * is optimized by bisection so that the overlap between the derivatives of
 * alpha and chi is minimized, and alpha is normalized so that the overlap
 * integral with chi equals one.
 *
 * @param[in]   nr      number of grid points
 * @param[in]   r       radial grid
 * @param[in]   chi     radial function
 * @param[in]   rcut    cutoff radius
 * @param[out]  alpha   smoothed radial function on the truncated grid
 */
void smoothgen(
    const int nr,
    const double* r,
    const double* chi,
    const double rcut,
    std::vector<double>& alpha);

#endif
