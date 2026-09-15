#ifndef FORCE_STRESS_PW_H
#define FORCE_STRESS_PW_H

#include "source_base/matrix.h"
#include "source_pw/module_pwdft/stress_func.h"

// Free functions for the plane-wave part of the LCAO force/stress and for
// force symmetrization. None of them depends on the electronic template type
// T (they operate purely on double-precision quantities), so they are kept
// outside the Force_Stress_LCAO<T> class template.

class UnitCell;
class Charge;
class pseudopot_cell_vl;
class Structure_Factor;
namespace ModulePW
{
class PW_Basis;
}
namespace ModuleSymmetry
{
class Symmetry;
}

namespace LCAO_domain
{

// vlocal, hartree, ewald, core correction, exchange-correlation stress terms.
void cal_stress_pw(Stress_Func<double>& sc_pw,
                   UnitCell& ucell,
                   ModuleBase::matrix& sigmadvl,
                   ModuleBase::matrix& sigmahar,
                   ModuleBase::matrix& sigmaewa,
                   ModuleBase::matrix& sigmacc,
                   ModuleBase::matrix& sigmaxc,
                   const double& etxc,
                   const Charge* const chr,
                   ModulePW::PW_Basis* rhopw,
                   const pseudopot_cell_vl& locpp,
                   const Structure_Factor& sf);

// Symmetrize the total force: Cartesian -> direct, symmetrize, direct ->
// Cartesian.
void symmetrize_force(const UnitCell& ucell, ModuleBase::matrix& fcs, ModuleSymmetry::Symmetry* symm);

} // namespace LCAO_domain

#endif
