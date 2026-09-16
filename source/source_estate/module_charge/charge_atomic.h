#ifndef CHARGE_ATOMIC_H
#define CHARGE_ATOMIC_H

#include "source_base/complexmatrix.h"
#include "source_basis/module_pw/pw_basis.h"

class UnitCell;

namespace module_charge
{

// Superposition of atomic charges contained in the array rho_at
// (read from pseudopotential files).
//
// spin_number_need is the number of spin components to be calculated:
//   1 -> total atomic charge density
//   2 -> spin up/down densities assuming uniform atomic polarization
//        equal to start_mag(it)
//   4 -> noncollinear case: total density in component 0, magnetization
//        vector in components 1..3
//
// NB: spin_number_need may differ from nspin (e.g. in update only the
// total charge is needed even in an LSDA calculation).
//
// All grid / basis inputs are passed explicitly via rhopw instead of
// being read from Charge members.
void atomic_rho(const int spin_number_need,
                const double& omega,
                double** rho_in,
                const ModuleBase::ComplexMatrix& strucFac,
                const UnitCell& ucell,
                const ModulePW::PW_Basis* rhopw);

} // namespace module_charge

#endif // CHARGE_ATOMIC_H
