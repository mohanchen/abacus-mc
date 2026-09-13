#ifndef DFTU_NAO_ADJ_H
#define DFTU_NAO_ADJ_H

/// @file dftu_nao_adj.h
/// @brief Build the adjacent-atom lists for all Hubbard atoms, shared by the
///        DFTU LCAO operator and the real-space force/stress path.

#include <vector>

class UnitCell;
class Plus_U_Base;
class Grid_Driver;
class AdjacentAtomInfo;

namespace DFTU_LCAO
{

/**
 * @brief build the adjacent-atom lists for all Hubbard atoms.
 *
 * For every atom whose type carries a Hubbard-U channel, find its neighbor
 * atoms within orb_cutoff + onsite_radius and filter the AdjacentAtomInfo
 * accordingly.
 *
 * @param ucell         [in] unit cell
 * @param dftu          [in] DFT+U base object (per-type U channels)
 * @param gridD         [in] grid driver for neighbor search
 * @param orb_cutoff    [in] orbital cutoff radius per atom type
 * @param onsite_radius [in] onsite projector cutoff radius
 * @return one AdjacentAtomInfo per Hubbard atom, in atom order
 */
std::vector<AdjacentAtomInfo> build_adjacent_atoms(const UnitCell* ucell,
                                                   Plus_U_Base* dftu,
                                                   const Grid_Driver* gridD,
                                                   const std::vector<double>& orb_cutoff,
                                                   const double onsite_radius);

} // namespace DFTU_LCAO

#endif // DFTU_NAO_ADJ_H
