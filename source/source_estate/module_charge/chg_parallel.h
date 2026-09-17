#ifndef CHG_PARALLEL_H
#define CHG_PARALLEL_H

// MPI reductions of charge-density arrays across k-point pools and band
// groups. Stateless free functions extracted from Charge member functions;
// the charge buffers and the parallel grid are supplied by the Charge
// argument. The pool/band conditions are still read from GlobalV and PARAM
// as in the original implementation (migration-neutral).

#ifdef __MPI

class Charge;

namespace module_charge
{

/**
 * @brief Reduce a real-space array across k-point pools and band groups.
 *
 * @param array_rho real-space array [chr.nrxx], reduced in place
 * @param chr charge object supplying the parallel grid and the local grid size
 */
void reduce_diff_pools(double* array_rho, const Charge& chr);

/**
 * @brief Reduce rho across pools; also reduce kin_r for meta-GGA or ELF.
 *
 * @param chr charge object supplying rho/kin_r buffers
 */
void rho_mpi(Charge& chr);

/**
 * @brief Reduce kin_r across pools for meta-GGA or ELF calculations.
 *
 * @param chr charge object supplying kin_r buffers
 */
void kin_r_mpi(Charge& chr);

} // namespace module_charge

#endif

#endif // CHG_PARALLEL_H
