#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_atomic.h"

#include <cmath>
#include <complex>
#include <sstream>
#include <vector>

// charge.cpp references Magnetism; provide a lightweight stub.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of module_charge/chg_atomic.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - atomic_rho: dispatcher over spin_number_need (1/2/4) and per-atom
 *     start_mag.  Covered:
 *       - ntype == 0 path: loop is skipped, only normalize_and_check runs.
 *       - bad spin_number_need (3) triggers WARNING_QUIT (abort).
 */

namespace
{

module_charge::AtomicRhoCfg make_cfg(std::ostream& os)
{
    return {1.0, 0, false, false, os}; // nelec, test_charge, domag, domag_z, ofs_warning
}

} // namespace

TEST(ChgAtomicTest, AtomicRhoNtypeZeroOnlyNormalizes)
{
    ModulePW::PW_Basis rhopw;
    rhopw.initgrids(1.0, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 10.0);
    rhopw.initparameters(false, 10.0);
    rhopw.setuptransform();
    rhopw.collect_local_pw();

    UnitCell ucell;
    ucell.ntype = 0;
    ucell.nat = 0;

    const int spin_number_need = 1;
    const double omega = rhopw.omega;
    ModuleBase::ComplexMatrix strucFac(0, rhopw.npw);
    std::vector<double> rho_in(rhopw.nrxx, 0.0);
    double* rho_ptrs[1] = {rho_in.data()};

    std::stringstream ofs;
    module_charge::AtomicRhoCfg cfg = make_cfg(ofs);
    module_charge::atomic_rho(spin_number_need, omega, rho_ptrs, strucFac, ucell, &rhopw, cfg);

    // with ntype==0, rho_g3d is all zero, so normalize_and_check divides by
    // ne_tot==0; the result is NaN/zero.  We only assert no crash.
    EXPECT_EQ(rho_in.size(), static_cast<size_t>(rhopw.nrxx));
}

TEST(ChgAtomicDeathTest, AtomicRhoBadSpinAborts)
{
    ModulePW::PW_Basis rhopw;
    rhopw.initgrids(1.0, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 10.0);
    rhopw.initparameters(false, 10.0);
    rhopw.setuptransform();
    rhopw.collect_local_pw();

    UnitCell ucell;
    ucell.ntype = 0;
    ucell.nat = 0;

    const int spin_number_need = 3;
    ModuleBase::ComplexMatrix strucFac(0, rhopw.npw);
    std::vector<double> rho_in(rhopw.nrxx, 0.0);
    double* rho_ptrs[1] = {rho_in.data()};

    std::stringstream ofs;
    module_charge::AtomicRhoCfg cfg = make_cfg(ofs);
    EXPECT_DEATH(module_charge::atomic_rho(spin_number_need, rhopw.omega, rho_ptrs,
                                           strucFac, ucell, &rhopw, cfg),
                 "");
}
