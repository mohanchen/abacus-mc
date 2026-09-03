#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_estate/rhog_io.h"
#include "source_base/module_parallel/para_world.h"
#include "source_base/module_parallel/para_tag.h"
#include "source_base/module_parallel/para_bridge.h"
#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#include "mpi.h"
#endif
#include <fstream>
#include <sstream>

/**
 * - Tested Functions:
 *  - read_rhog()
 *  - write_rhog()
 */

class ReadRhogTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis rhopw;
    std::vector<std::vector<std::complex<double>>> rhog_data;
    std::vector<std::complex<double>*> rhog;
    Parallel::ParaWorld pw_world = Parallel::make_pw_world();
    std::ofstream warning_stream;

    void setup_pw_basis()
    {
#ifdef __MPI
        rhopw.initmpi(pw_world.size(), pw_world.rank(), pw_world.comm());
#endif
        rhopw.initgrids(6.5, ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0), 120);
        rhopw.initparameters(false, 120);
        rhopw.setuptransform();
        rhopw.collect_local_pw();
    }

    void open_warning(const std::string& path)
    {
        warning_stream.open(path);
    }

    void close_warning()
    {
        if (warning_stream.is_open())
        {
            warning_stream.close();
        }
    }

    std::string read_warning_file(const std::string& path)
    {
        std::ifstream ifs(path);
        std::stringstream ss;
        ss << ifs.rdbuf();
        ifs.close();
        return ss.str();
    }

    virtual void SetUp()
    {
        rhog_data.resize(1, std::vector<std::complex<double>>(1471));
        rhog.push_back(rhog_data[0].data());
    }

    virtual void TearDown()
    {
        close_warning();
    }
};

// Test the read_rhog function with normal file
TEST_F(ReadRhogTest, ReadRhog)
{
    std::string filename = "./support/charge-density.dat";
    setup_pw_basis();

    bool result = ModuleIO::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, nullptr);

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

    open_warning("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content = " ModuleIO::read_rhog  warning : Can't open file notfound.txt\n";
    EXPECT_FALSE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when gamma_only is inconsistent
TEST_F(ReadRhogTest, InconsistentGammaOnly)
{
    std::string filename = "./support/charge-density.dat";
    rhopw.gamma_only = true;

    open_warning("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, &rhopw, 2, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content
        = " ModuleIO::read_rhog  warning : some planewaves in file are not used\n ModuleIO::read_rhog  warning : some "
          "spin channels in file are missing\n ModuleIO::read_rhog  warning : gamma_only read from file is "
          "inconsistent with INPUT\n";

    EXPECT_FALSE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove("test_read_rhog.txt");
}

// Test the read_rhog function when some planewaves in file are missing
TEST_F(ReadRhogTest, SomePWMissing)
{
    std::string filename = "./support/charge-density.dat";
    rhopw.npwtot = 2000;

    open_warning("test_read_rhog.txt");
    bool result = ModuleIO::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, &warning_stream);
    close_warning();

    std::string expected_content = " ModuleIO::read_rhog  warning : some planewaves in file are missing\n";
    EXPECT_TRUE(result);
    EXPECT_EQ(read_warning_file("test_read_rhog.txt"), expected_content);
    std::remove("test_read_rhog.txt");
}

// Test read_rhog with os_warning=nullptr (silent mode, must not crash)
TEST_F(ReadRhogTest, OsNullptrSilent)
{
    std::string filename = "notfound.txt";
    bool result = ModuleIO::read_rhog(filename, &rhopw, 1, rhog.data(), pw_world, nullptr);
    EXPECT_FALSE(result);
}

// Test write_rhog round-trip: write then read back, verify data consistency
TEST_F(ReadRhogTest, WriteRoundTrip)
{
    setup_pw_basis();

    // initialize some rhog data
    rhog_data[0].assign(rhopw.npw, std::complex<double>(1.5, 2.5));

    std::string tmpfile = "test_rhog_roundtrip.dat";

    // write
    bool write_result = ModuleIO::write_rhog(
        tmpfile, rhopw.gamma_only, &rhopw, 1,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    // read back into a fresh buffer
    std::vector<std::vector<std::complex<double>>> rhog_read_data(
        1, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    rhog_read.push_back(rhog_read_data[0].data());

    bool read_result = ModuleIO::read_rhog(tmpfile, &rhopw, 1, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

    // compare: within MPI precision tolerance
    int diff_count = 0;
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        if (std::abs(rhog[0][ig] - rhog_read[0][ig]) > 1e-10)
        {
            ++diff_count;
        }
    }
    EXPECT_EQ(diff_count, 0) << diff_count << " planewave values differ after round-trip";

    std::remove(tmpfile.c_str());
}

// Test write_rhog when the output path is not writable
TEST_F(ReadRhogTest, WriteFileFail)
{
    setup_pw_basis();
    rhog_data[0].assign(rhopw.npw, std::complex<double>(1.0, 0.0));

    // try to write to a directory path (not a file) — should fail
    bool result = ModuleIO::write_rhog(
        "/tmp", rhopw.gamma_only, &rhopw, 1,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_FALSE(result);
}

// Test write_rhog with nspin=2, round-trip both channels
TEST_F(ReadRhogTest, WriteRoundTripNspin2)
{
    setup_pw_basis();

    // expand to nspin=2
    rhog_data.resize(2, std::vector<std::complex<double>>(rhopw.npw));
    rhog.clear();
    rhog.push_back(rhog_data[0].data());
    rhog.push_back(rhog_data[1].data());

    // initialize distinct values for each spin channel
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        rhog_data[0][ig] = std::complex<double>(1.0 * ig, 0.1 * ig);
        rhog_data[1][ig] = std::complex<double>(2.0 * ig, 0.2 * ig);
    }

    std::string tmpfile = "test_rhog_roundtrip_nspin2.dat";

    // write nspin=2
    bool write_result = ModuleIO::write_rhog(
        tmpfile, rhopw.gamma_only, &rhopw, 2,
        ModuleBase::Matrix3(-0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.5, 0.0),
        rhog.data(), pw_world, nullptr);
    EXPECT_TRUE(write_result);

    // read back
    std::vector<std::vector<std::complex<double>>> rhog_read_data(
        2, std::vector<std::complex<double>>(rhopw.npw));
    std::vector<std::complex<double>*> rhog_read;
    rhog_read.push_back(rhog_read_data[0].data());
    rhog_read.push_back(rhog_read_data[1].data());

    bool read_result = ModuleIO::read_rhog(tmpfile, &rhopw, 2, rhog_read.data(), pw_world, nullptr);
    EXPECT_TRUE(read_result);

    int diff_count = 0;
    for (int is = 0; is < 2; ++is)
    {
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            if (std::abs(rhog[is][ig] - rhog_read[is][ig]) > 1e-10)
            {
                ++diff_count;
            }
        }
    }
    EXPECT_EQ(diff_count, 0) << diff_count << " planewave values differ after nspin=2 round-trip";

    std::remove(tmpfile.c_str());
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
