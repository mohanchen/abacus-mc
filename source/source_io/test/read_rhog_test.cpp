#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_io/module_chgpot/rhog_io.h"
#include "source_base/module_parallel/para_world.h"
#include "source_base/module_parallel/para_tag.h"
#include "source_base/module_parallel/para_bridge.h"
#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif

/**
 * - Tested Functions:
 *  - read_rhog()
 */

class ReadRhogTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis rhopw;
    std::vector<std::vector<std::complex<double>>> rhog_data;
    std::vector<std::complex<double>*> rhog;
    Parallel::ParaWorld pw_world = Parallel::make_pw_world();

    virtual void SetUp()
    {
        rhog_data.resize(1, std::vector<std::complex<double>>(1471));
        rhog.push_back(rhog_data[0].data());
    }
};

// Test the read_rhog function
TEST_F(ReadRhogTest, ReadRhog)
{
    std::string filename = "./support/charge-density.dat";
#ifdef __MPI
    rhopw.initmpi(pw_world.size(), pw_world.rank(), pw_world.comm());
#endif
    rhopw.initgrids(6.5, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 120);
    rhopw.initparameters(false, 120);
    rhopw.setuptransform();
    rhopw.collect_local_pw();

    bool result = ModuleIO::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &GlobalV::ofs_warning)

    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(rhog[0][0].real(), -1.0304462993299456e-05);
    EXPECT_DOUBLE_EQ(rhog[0][0].imag(), -1.2701788626185278e-13);
    EXPECT_DOUBLE_EQ(rhog[0][1].real(), -0.0003875762482855959);
    EXPECT_DOUBLE_EQ(rhog[0][1].imag(), -4.2556814316812048e-12);
    EXPECT_DOUBLE_EQ(rhog[0][1470].real(), -3.5683133614445107e-05);
    EXPECT_DOUBLE_EQ(rhog[0][1470].imag(), 1.6176615686863767e-12);
}

// Test the read_rhog function when the file is not found
TEST_F(ReadRhogTest, NotFoundFile)
{
    std::string filename = "notfound.txt";

    GlobalV::ofs_warning.open("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &GlobalV::ofs_warning)
    GlobalV::ofs_warning.close();

    std::ifstream ifs_running("test_read_rhog.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = " ModuleIO::read_rhog  warning : Can't open file notfound.txt\n";

    EXPECT_FALSE(result);
    EXPECT_EQ(file_content, expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when tgamma_only is inconsistent
TEST_F(ReadRhogTest, InconsistentGammaOnly)
{
    std::string filename = "./support/charge-density.dat";
    rhopw.gamma_only = true;

    GlobalV::ofs_warning.open("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, &rhopw, 2, rhog.data(), pw_world, &GlobalV::ofs_warning);
    GlobalV::ofs_warning.close();

    std::ifstream ifs_running("test_read_rhog.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content
        = " ModuleIO::read_rhog  warning : some planewaves in file are not used\n ModuleIO::read_rhog  warning : some "
          "spin channels in file are missing\n ModuleIO::read_rhog  warning : gamma_only read from file is "
          "inconsistent with INPUT\n";

    EXPECT_FALSE(result);
    EXPECT_EQ(file_content, expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when some planewaves in file are missing
TEST_F(ReadRhogTest, SomePWMissing)
{
    std::string filename = "./support/charge-density.dat";
    rhopw.npwtot = 2000;

    GlobalV::ofs_warning.open("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &GlobalV::ofs_warning)
    GlobalV::ofs_warning.close();

    std::ifstream ifs_running("test_read_rhog.txt");
    std::stringstream ss;
    ss << ifs_running.rdbuf();
    std::string file_content = ss.str();
    ifs_running.close();

    std::string expected_content = " ModuleIO::read_rhog  warning : some planewaves in file are missing\n";

    EXPECT_TRUE(result);
    EXPECT_EQ(file_content, expected_content);
    std::remove("test_read_rhog.txt");
}

int main(int argc, char** argv)
{
#ifdef __MPI
    int nproc = 1;
    int myrank = 0;
    int nproc_in_pool = 1;
    const int kpar = 1;
    int mypool = 0;
    int rank_in_pool = 0;
    setupmpi(argc, argv, nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    finishmpi();
#endif
    return result;
}
