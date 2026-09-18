#ifndef CHG_INIT_H
#define CHG_INIT_H

#include "source_base/complexmatrix.h"
#include "source_base/parallel_grid.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"

class Charge;
class UnitCell;

namespace module_charge
{

/**
 * @brief Initialize charge density from file, atomic superposition,
 *        restart cache, or wavefunctions, depending on INPUT.init_chg.
 *
 * @param chr [inout] Charge object whose rho/rhog/kin_r buffers are filled.
 * @param ucell [in] unit cell
 * @param pgrid [in] parallel grid descriptor
 * @param strucFac [in] structure factor for atomic-charge superposition
 * @param symm [in] symmetry operations (used by wfc-based init)
 * @param klist [in] k-point list pointer (K_Vectors*), needed only for wfc init
 * @param wfcpw [in] PW_Basis_K pointer, needed only for wfc init
 */
void init_rho(Charge& chr,
              const UnitCell& ucell,
              const Parallel_Grid& pgrid,
              const ModuleBase::ComplexMatrix& strucFac,
              ModuleSymmetry::Symmetry& symm,
              const void* klist,
              const void* wfcpw);

} // namespace module_charge

#endif // CHG_INIT_H
