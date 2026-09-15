#include "force_stress_pw.h"

#include "source_base/mathzone.h"
#include "source_base/timer.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_pw/module_pwdft/stru_fac.h"
#include "source_pw/module_pwdft/vl_pw.h"

namespace LCAO_domain
{

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
                   const Structure_Factor& sf)
{
    ModuleBase::TITLE("Force_Stress_LCAO", "cal_stress_pw");

    // local pseudopotential stress:
    sc_pw.stress_loc(ucell, sigmadvl, rhopw, locpp.vloc, &sf, 0, chr);

    // hartree term
    sc_pw.stress_har(ucell, sigmahar, rhopw, 0, chr);

    // ewald stress: use plane wave only.
    sc_pw.stress_ewa(ucell, sigmaewa, rhopw, 0); // remain problem

    // stress due to core correlation.
    sc_pw.stress_cc(sigmacc, rhopw, ucell, &sf, 0, locpp.numeric, chr);

    // stress due to self-consistent charge.
    for (int i = 0; i < 3; i++)
    {
        sigmaxc(i, i) = -etxc / ucell.omega;
    }
    // Exchange-correlation for PBE
    sc_pw.stress_gga(ucell, sigmaxc, rhopw, chr);

    return;
}

void symmetrize_force(const UnitCell& ucell, ModuleBase::matrix& fcs, ModuleSymmetry::Symmetry* symm)
{
    double d1;
    double d2;
    double d3;
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        ModuleBase::Mathzone::Cartesian_to_Direct(fcs(iat, 0), fcs(iat, 1), fcs(iat, 2),
          ucell.a1.x, ucell.a1.y, ucell.a1.z, ucell.a2.x, ucell.a2.y, ucell.a2.z,
          ucell.a3.x, ucell.a3.y, ucell.a3.z, d1, d2, d3);

        fcs(iat, 0) = d1;
        fcs(iat, 1) = d2;
        fcs(iat, 2) = d3;
    }
    symm->symmetrize_vec3_nat(fcs.c);
    for (int iat = 0; iat < ucell.nat; iat++)
    {
        ModuleBase::Mathzone::Direct_to_Cartesian(fcs(iat, 0), fcs(iat, 1), fcs(iat, 2),
          ucell.a1.x, ucell.a1.y, ucell.a1.z, ucell.a2.x, ucell.a2.y, ucell.a2.z,
          ucell.a3.x, ucell.a3.y, ucell.a3.z, d1, d2, d3);

        fcs(iat, 0) = d1;
        fcs(iat, 1) = d2;
        fcs(iat, 2) = d3;
    }
    return;
}

} // namespace LCAO_domain
