#ifndef FORCE_STRESS_ASSEMBLE_H
#define FORCE_STRESS_ASSEMBLE_H

#include "force_stress_lcao.h"

// Free functions that assemble the per-term LCAO force/stress parts into the
// total force/stress matrices and print the breakdown. They only depend on the
// electronic template type T through the Setup_DeePKS<T> argument, so they are
// function templates with explicit instantiation in the .cpp.

class UnitCell;
namespace ModuleSymmetry
{
class Symmetry;
}
namespace vdw
{
struct VdwResult;
}

namespace LCAO_domain
{

// Sum the computed force parts into fcs, apply symmetry and the net-force
// (drift) correction, then print the per-term and total forces.
// force_threshold is Force_Stress_LCAO<T>::force_invalid_threshold_ev, passed
// in explicitly so this function need not be a class member.
template <typename T>
void assemble_print_force(const UnitCell& ucell,
                          const bool istestf,
                          const vdw::VdwResult* vdw_result,
                          const Exx_Info& exx_info,
                          ModuleSymmetry::Symmetry* symm,
                          Setup_DeePKS<T>& deepks,
                          const LCAOForceParts& parts,
                          const double force_threshold,
                          ModuleBase::matrix& fcs);

// Sum the computed stress parts into scs, symmetrize, subtract the external
// pressure and print the per-term and total stresses.
template <typename T>
void assemble_print_stress(const UnitCell& ucell,
                           const bool istests,
                           const vdw::VdwResult* vdw_result,
                           const Exx_Info& exx_info,
                           ModuleSymmetry::Symmetry* symm,
                           Setup_DeePKS<T>& deepks,
                           const LCAOStressParts& sparts,
                           ModuleBase::matrix& scs);

} // namespace LCAO_domain

#endif
