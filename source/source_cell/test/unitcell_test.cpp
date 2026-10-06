#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_cell/read_stru.h"
#include "memory"
#include "source_cell/read_stru.h"
#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "prepare_unitcell.h"
#include "source_cell/read_stru.h"
#include <streambuf>
#include <valarray>
#include <vector>


Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of class UnitCell
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor:
 *     - UnitCell() and ~UnitCell()
 *   - Setup:
 *     - setup_from_input(): to set latname, ntype, lmaxmax, init_vel, and lc
 *   - Index
 *     - set_iat2iait(): set index relations in two arrays of Unitcell: iat2it[nat], iat2ia[nat]
 *     - iat2iait(): depends on the above function, but can find both ia & it from iat
 *     - ijat2iaitjajt(): find ia, it, ja, jt from ijat (ijat_max = nat*nat)
 *         which collapses it, ia, jt, ja loop into a single loop
 *   - GetAtomCounts
 *     - get_atom_Counts(): get atomCounts, which is a map from atom type to atom number
 *   - GetOrbitalCounts
 *     - get_orbitalCounts(): get orbitalCounts, which is a map from atom type to orbital number
 *   - CheckDTau
 *     - check_dtau(): move all atomic coordinates into the first unitcell, i.e. in between [0,1)
 *   - CheckTau
 *     - check_tau(): check if any "two atoms are too close"
 *   - ReadAtomSpecies
 *     - read_atom_species(): a successful case
 *   - ReadAtomSpeciesWarning1
 *     - read_atom_species(): unrecognized pseudopotential type.
 *   - ReadAtomSpeciesWarning2
 *     - read_atom_species(): lat0<=0.0
 *   - ReadAtomSpeciesWarning3
 *     - read_atom_species(): do not use LATTICE_PARAMETERS without explicit specification of lattice type
 *   - ReadAtomSpeciesWarning4
 *     - read_atom_species():do not use LATTICE_VECTORS along with explicit specification of lattice type
 *   - ReadAtomSpeciesWarning5
 *     - read_atom_species():latname not supported
 *   - ReadAtomSpeciesLatName
 *     - read_atom_species(): various latname
 *   - ReadAtomPositionsS1
 *     - read_atom_positions(): spin 1 case
 *   - ReadAtomPositionsS2
 *     - read_atom_positions(): spin 2 case
 *   - ReadAtomPositionsS4Noncolin
 *     - read_atom_positions(): spin 4 noncolinear case
 *   - ReadAtomPositionsS4Colin
 *     - read_atom_positions(): spin 4 colinear case
 *   - ReadAtomPositionsC
 *     - read_atom_positions(): Cartesian coordinates
 *   - ReadAtomPositionsCA
 *     - read_atom_positions(): Cartesian_angstrom coordinates
 *   - ReadAtomPositionsCACXY
 *     - read_atom_positions(): Cartesian_angstrom_center_xy coordinates
 *   - ReadAtomPositionsCACXZ
 *     - read_atom_positions(): Cartesian_angstrom_center_xz coordinates
 *   - ReadAtomPositionsCACXYZ
 *     - read_atom_positions(): Cartesian_angstrom_center_xyz coordinates
 *   - ReadAtomPositionsCAU
 *     - read_atom_positions(): Cartesian_au coordinates
 *   - ReadAtomPositionsWarning1
 *     - read_atom_positions(): unknown type of coordinates
 *   - ReadAtomPositionsWarning2
 *     - read_atom_positions(): atomic label inconsistency between ATOM_POSITIONS
 *                              and ATOM_SPECIES
 *   - ReadAtomPositionsWarning3
 *     - read_atom_positions(): warning :  atom number < 0
 *   - ReadAtomPositionsWarning4
 *     - read_atom_positions(): mismatch in atom number for atom type
 *   - ReadAtomPositionsWarning5
 *     - read_atom_positions(): no atoms can move in MD simulations!
 */

class UcellTest : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell{new UnitCell};
    std::string output;
};

using UcellDeathTest = UcellTest;

TEST_F(UcellTest, Constructor)
{
    EXPECT_EQ(ucell->Coordinate, "Direct");
    EXPECT_EQ(ucell->latName, "user_defined_lattice");
    EXPECT_DOUBLE_EQ(ucell->lat0, 0.0);
    EXPECT_DOUBLE_EQ(ucell->lat0_angstrom, 0.0);
    EXPECT_EQ(ucell->ntype, 0);
    EXPECT_EQ(ucell->nat, 0);
    EXPECT_EQ(ucell->namax, 0);
    EXPECT_EQ(ucell->nwmax, 0);
    EXPECT_EQ(ucell->iat2it, nullptr);
    EXPECT_EQ(ucell->iat2ia, nullptr);
    EXPECT_EQ(ucell->iwt2iat, nullptr);
    EXPECT_EQ(ucell->iwt2iw, nullptr);
    EXPECT_DOUBLE_EQ(ucell->tpiba, 0.0);
    EXPECT_DOUBLE_EQ(ucell->tpiba2, 0.0);
    EXPECT_DOUBLE_EQ(ucell->omega, 0.0);
    EXPECT_FALSE(ucell->set_atom_flag);
}

TEST_F(UcellTest, Setup)
{
    std::string latname_in = "bcc";
    int ntype_in = 1;
    int lmaxmax_in = 2;
    bool init_vel_in = false;
    std::vector<std::string> fixed_axes_in = {"None", "volume", "shape", "a", "b", "c", "ab", "ac", "bc", "abc"};
    for (int i = 0; i < fixed_axes_in.size(); ++i)
    {
        ucell->setup_from_input(latname_in, ntype_in, lmaxmax_in, init_vel_in, fixed_axes_in[i]);
        EXPECT_EQ(ucell->latName, latname_in);
        EXPECT_EQ(ucell->ntype, ntype_in);
        EXPECT_EQ(ucell->lmaxmax, lmaxmax_in);
        EXPECT_EQ(ucell->init_vel, init_vel_in);
        if (fixed_axes_in[i] == "None" || fixed_axes_in[i] == "volume" || fixed_axes_in[i] == "shape")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "a")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "b")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "c")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
        else if (fixed_axes_in[i] == "ab")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "ac")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
        else if (fixed_axes_in[i] == "bc")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
        else if (fixed_axes_in[i] == "abc")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
    }
}

TEST_F(UcellTest, Index)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    int iat = 0;
    for (int it = 0; it < utp.natom.size(); ++it)
    {
        for (int ia = 0; ia < utp.natom[it]; ++ia)
        {
            EXPECT_EQ(ucell->iat2it[iat], it);
            EXPECT_EQ(ucell->iat2ia[iat], ia);
            // test iat2iait
            int ia_beg, it_beg;
            ucell->iat2iait(iat, &ia_beg, &it_beg);
            EXPECT_EQ(it_beg, it);
            EXPECT_EQ(ia_beg, ia);
            ++iat;
        }
    }
    // test iat2iait: case of (iat >= nat)
    int ia_beg2;
    int it_beg2;
    long long iat2 = ucell->nat + 1;
    EXPECT_FALSE(ucell->iat2iait(iat2, &ia_beg2, &it_beg2));
    // test ijat2iaitjajt
    int ia_test;
    int it_test;
    int ja_test;
    int jt_test;
    long long ijat = 0;
    for (int it = 0; it < utp.natom.size(); ++it)
    {
        for (int ia = 0; ia < utp.natom[it]; ++ia)
        {
            for (int jt = 0; jt < utp.natom.size(); ++jt)
            {
                for (int ja = 0; ja < utp.natom[jt]; ++ja)
                {
                    ucell->ijat2iaitjajt(ijat, &ia_test, &it_test, &ja_test, &jt_test);
                    EXPECT_EQ(ia_test, ia);
                    EXPECT_EQ(it_test, it);
                    EXPECT_EQ(ja_test, ja);
                    EXPECT_EQ(jt_test, jt);
                    ++ijat;
                }
            }
        }
    }
}

TEST_F(UcellTest, GetAtomCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    std::map<int, int> atomCounts = ucell->get_atom_Counts();
    EXPECT_EQ(atomCounts[0], 1);
    EXPECT_EQ(atomCounts[1], 2);
}

TEST_F(UcellTest, GetOrbitalCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    std::map<int, int> orbitalCounts = ucell->get_orbital_Counts();
    EXPECT_EQ(orbitalCounts[0], 9);
    EXPECT_EQ(orbitalCounts[1], 9);
}

TEST_F(UcellTest, GetLnchiCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    std::map<int, std::map<int, int>> LnchiCounts = ucell->get_lnchi_Counts();
    EXPECT_EQ(LnchiCounts[0][0], 1);
    EXPECT_EQ(LnchiCounts[0][1], 1);
    EXPECT_EQ(LnchiCounts[0][2], 1);
    EXPECT_EQ(LnchiCounts[1][0], 1);
    EXPECT_EQ(LnchiCounts[1][1], 1);
    EXPECT_EQ(LnchiCounts[1][2], 1);
}

TEST_F(UcellTest, CheckDTau)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-CheckDTau"];
    ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    for (int it = 0; it < utp.natom.size(); ++it)
    {
        for (int ia = 0; ia < utp.natom[it]; ++ia)
        {
            EXPECT_GE(ucell->atoms[it].taud[ia].x, 0);
            EXPECT_GE(ucell->atoms[it].taud[ia].y, 0);
            EXPECT_GE(ucell->atoms[it].taud[ia].z, 0);
            EXPECT_LT(ucell->atoms[it].taud[ia].x, 1);
            EXPECT_LT(ucell->atoms[it].taud[ia].y, 1);
            EXPECT_LT(ucell->atoms[it].taud[ia].z, 1);
        }
    }
}

TEST_F(UcellTest, CheckTauFalse)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-CheckTau"];
    ucell = utp.SetUcellInfo();
    GlobalV::ofs_warning.open("checktau_warning");
    unitcell::check_tau(ucell->atoms ,ucell->ntype, ucell->lat0);
    GlobalV::ofs_warning.close();
    std::ifstream ifs;
    ifs.open("checktau_warning");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("two atoms are too close!"));
    ifs.close();
    remove("checktau_warning");
}

TEST_F(UcellTest, CheckTauTrue)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-CheckTau"];
    ucell = utp.SetUcellInfo();
    GlobalV::ofs_warning.open("checktau_warning");
    int atom=0;
    //cause the ucell->lat0 is 0.5,if the type of the check_tau has 
    //an int type,it will set to zero,and it will not pass the unittest
    ucell->lat0=0.5;
    ucell->nat=3;
    for (int it=0;it<ucell->ntype;it++)
    {
        for(int ia=0; ia<ucell->atoms[it].na; ++ia)
        {
            
            for (int i=0;i<3;i++)
            {
                ucell->atoms[it].tau[ia][i]=((atom+i)/(ucell->nat*3.0));
            }
            atom+=3;
        }
    }
    EXPECT_EQ(unitcell::check_tau(ucell->atoms ,ucell->ntype, ucell->lat0),true);
    GlobalV::ofs_warning.close();
}

#ifdef __LCAO
class UcellTestReadStru : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell{new UnitCell};
    std::string output;
      void SetUp() override
    {
        ucell->ntype = 2;
        ucell->pseudo_fn.resize(ucell->ntype);
        ucell->pseudo_type.resize(ucell->ntype);
        ucell->orbital_fn.resize(ucell->ntype);
    }
    void TearDown() override
    {
        ucell->orbital_fn.shrink_to_fit();
    }
};

TEST_F(UcellTestReadStru, ReadAtomSpecies)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running, ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

TEST_F(UcellTestReadStru, ReadAtomSpeciesWarning1)
{
    std::string fn = "./support/STRU_MgO_Warning1";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.txt");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("unrecognized pseudopotential type."));
    ofs_running.close();
    ifa.close();
    //remove("read_atom_species.txt");
}

TEST_F(UcellTestReadStru, ReadLatticeConstantWarning1)
{
    std::string fn = "./support/STRU_MgO_Warning2";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species1.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Lattice constant <= 0.0"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species1.tmp");
}

TEST_F(UcellTestReadStru, ReadLatticeConstantWarning2)
{
    std::string fn = "./support/STRU_MgO_Warning3";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("do not use LATTICE_PARAMETERS without explicit specification of lattice type"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

TEST_F(UcellTestReadStru, ReadLatticeConstantWarning3)
{
    std::string fn = "./support/STRU_MgO_Warning4";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    ucell->latName = "bcc";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,
                testing::HasSubstr("do not use LATTICE_VECTORS along with explicit specification of lattice type"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

TEST_F(UcellTestReadStru, ReadAtomSpeciesLatName)
{
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    std::vector<std::string> latName_in = {"sc",
                                           "fcc",
                                           "bcc",
                                           "hexagonal",
                                           "trigonal",
                                           "st",
                                           "bct",
                                           "so",
                                           "baco",
                                           "fco",
                                           "bco",
                                           "sm",
                                           "bacm",
                                           "triclinic"};
    for (int i = 0; i < latName_in.size(); ++i)
    {
        std::string fn = "./support/STRU_MgO_LatName";
        std::ifstream ifa(fn.c_str());
        std::ofstream ofs_running;
        ofs_running.open("read_atom_species.tmp");
        ucell->latName = latName_in[i];
        EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
        if (ucell->latName == "sc")
        {
            EXPECT_DOUBLE_EQ(ucell->latvec.e11, 1.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e22, 1.0);
            EXPECT_DOUBLE_EQ(ucell->latvec.e33, 1.0);
        }
        ofs_running.close();
        ifa.close();
        remove("read_atom_species.tmp");
    }
}

TEST_F(UcellDeathTest, ReadAtomSpeciesWarning5)
{
    std::string fn = "./support/STRU_MgO_LatName";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_species.tmp");
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    ucell->latName = "arbitrary";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("latname not supported"));
    ofs_running.close();
    ifa.close();
    remove("read_atom_species.tmp");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsS1)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsS2)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 2;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsS4Noncolin)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 4;
    const bool fixed_atoms = false;
    const bool noncolin = true;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsS4Colin)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 4;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsC)
{
    std::string fn = "./support/STRU_MgO_c";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsCA)
{
    std::string fn = "./support/STRU_MgO_ca";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsCACXY)
{
    std::string fn = "./support/STRU_MgO_cacxy";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsCACXZ)
{
    std::string fn = "./support/STRU_MgO_cacxz";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsCACYZ)
{
    std::string fn = "./support/STRU_MgO_cacyz";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsCACXYZ)
{
    std::string fn = "./support/STRU_MgO_cacxyz";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsCAU)
{
    std::string fn = "./support/STRU_MgO_cau";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = true;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsAutosetMag)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    int nspin = 2;
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    for (int it = 0; it < ucell->ntype; it++)
    {
        for (int ia = 0; ia < ucell->atoms[it].na; ia++)
        {
            EXPECT_DOUBLE_EQ(ucell->atoms[it].mag[ia], 1.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].z, 1.0);
        }
    }
    // for nspin == 4
    // Issue #5939: nspin=4 with no mag in STRU no longer autosets (1,1,1);
    // all moments stay zero and a warning is emitted instead.
    nspin = 4;
    testing::internal::CaptureStdout();
    unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0);
    const std::string stdout_output = testing::internal::GetCapturedStdout();
    for (int it = 0; it < ucell->ntype; it++)
    {
        for (int ia = 0; ia < ucell->atoms[it].na; ia++)
        {
            EXPECT_DOUBLE_EQ(ucell->atoms[it].mag[ia], 0.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].x, 0.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].y, 0.0);
            EXPECT_DOUBLE_EQ(ucell->atoms[it].m_loc_[ia].z, 0.0);
        }
    }
    // The zero-moment warning must reach both stdout and the running log.
    EXPECT_NE(stdout_output.find("no initial magnetization is set in STRU"),
              std::string::npos);
    ofs_running.flush();
    std::ifstream ifs_log("read_atom_positions.tmp");
    std::string log_content;
    std::string log_line;
    while (std::getline(ifs_log, log_line))
    {
        log_content += log_line;
    }
    ifs_log.close();
    EXPECT_NE(log_content.find("no initial magnetization is set in STRU"),
              std::string::npos);
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsWarning1)
{
    std::string fn = "./support/STRU_MgO_WarningC1";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("There are several options for you:"));
    EXPECT_THAT(str, testing::HasSubstr("Direct"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_au"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xy"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xz"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_yz"));
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom_center_xyz"));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsWarning2)
{
    std::string fn = "./support/STRU_MgO_WarningC2";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Label read from ATOMIC_POSITIONS is Mo"));
    EXPECT_THAT(str, testing::HasSubstr("Label from ATOMIC_SPECIES is Mg"));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsWarning3)
{
    std::string fn = "./support/STRU_MgO_WarningC3";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_positions.tmp");
    GlobalV::ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, GlobalV::ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    GlobalV::ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("read_atom_positions  warning :  atom number < 0."));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsWarning4)
{
    std::string fn = "./support/STRU_MgO_WarningC4";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    ofs_running.open("read_atom_positions.tmp");
    ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->orbital_fn.resize(ucell->ntype);
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_atom_positions(*ucell, ifa, ofs_running, ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("read_atom_positions, mismatch in atom number for atom type: Mg"));
    ofs_running.close();
    ofs_warning.close();
    ifa.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}

TEST_F(UcellTestReadStru, ReadAtomPositionsWarning5)
{
    std::string fn = "./support/STRU_MgO";
    std::ifstream ifa(fn.c_str());
    std::ofstream ofs_running;
    ofs_running.open("read_atom_positions.tmp");
    GlobalV::ofs_warning.open("read_atom_positions.warn");
    // mandatory preliminaries
    ucell->ntype = 2;
    ucell->atoms = new Atom[ucell->ntype];
    ucell->set_atom_flag = true;
    const std::string basis_type = "lcao";
    const std::string orbital_dir = "";
    const std::string init_wfc = "";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = true;
    const bool rpa = false;
    const int nspin = 1;
    const bool fixed_atoms = true;
    const bool noncolin = false;
    const std::string calculation = "md";
    const std::string esolver_type = "ksdft";
    EXPECT_NO_THROW(unitcell::read_atom_species(ifa, ofs_running, *ucell,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa));
    EXPECT_NO_THROW(unitcell::read_lattice_constant(ifa, ofs_running,ucell->lat));
    EXPECT_DOUBLE_EQ(ucell->latvec.e11, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e22, 4.27957);
    EXPECT_DOUBLE_EQ(ucell->latvec.e33, 4.27957);
    EXPECT_NO_THROW(unitcell::read_atom_positions(*ucell, ifa, ofs_running, GlobalV::ofs_warning, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, fixed_atoms, noncolin,
        calculation, esolver_type, 0));
    ofs_running.close();
    GlobalV::ofs_warning.close();
    ifa.close();
    // check warning file
    std::ifstream ifs_tmp;
    ifs_tmp.open("read_atom_positions.warn");
    std::string str((std::istreambuf_iterator<char>(ifs_tmp)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("read_atoms  warning : no atoms can move in MD simulations!"));
    ifs_tmp.close();
    remove("read_atom_positions.tmp");
    remove("read_atom_positions.warn");
}
#endif
