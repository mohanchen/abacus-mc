#ifndef CHG_DMR_H
#define CHG_DMR_H

// Stateless real-space density-matrix (DMR) mixing kernels extracted from
// Charge_Mixing. The mixing history (Mixing_Data) and the Mixing strategy
// object remain owned by Charge_Mixing and are passed in explicitly; the
// functions do not read Charge_Mixing members or PARAM/GlobalV.

#include "chg_mix_cfg.h"

namespace Base_Mixing
{
class Mixing;
class Mixing_Data;
} // namespace Base_Mixing

namespace elecstate
{
template <typename TK, typename TR>
class DensityMatrix;
} // namespace elecstate

namespace module_charge
{

/**
 * @brief Allocate the mixing buffer for the real-space density matrix and
 *        clear its history.
 *
 * The buffer cannot be allocated in Charge_Mixing::set_mixing(): its length
 * nnr (number of non-zero R-matrix elements) is only known after
 * DensityMatrix::init_DMR(), which runs later in beforescf().
 *
 * @param mixing    mixing strategy object, non-null
 * @param mdata     mixing history buffer for DMR, resized and reset in place
 * @param nnr       number of real-space density-matrix elements per spin, > 0
 * @param cfg       mixing config (nspin and scf_thr_type select the path)
 */
void init_mixing_dmr(Base_Mixing::Mixing* mixing,
                     Base_Mixing::Mixing_Data& mdata,
                     const int nnr,
                     const MixingConfig& cfg);

/**
 * @brief Mix the real-space density matrix (LCAO calculations only).
 *
 * For nspin == 1/4 the single spin channel is mixed directly; for nspin == 2
 * the up/down channels are transformed into charge/magnetization channels,
 * mixed with independent betas, and transformed back.
 *
 * @tparam TK        scalar type of the density matrix (double or
 *                   std::complex<double>); the DMR storage itself is real
 * @param dm         density-matrix object supplying DMR and DMR_save
 * @param mixing     mixing strategy object, non-null
 * @param mdata      DMR mixing history buffer
 * @param cfg        mixing config (nspin and the two mixing betas)
 */
template <typename TK>
void mix_dmr(elecstate::DensityMatrix<TK, double>* dm,
             Base_Mixing::Mixing* mixing,
             Base_Mixing::Mixing_Data& mdata,
             const MixingConfig& cfg);

} // namespace module_charge

#endif // CHG_DMR_H
