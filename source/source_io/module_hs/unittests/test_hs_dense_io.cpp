/**
 * @file test_hs_dense_io.cpp
 * @brief Unit tests for ModuleIO::save_mat
 *
 * save_mat writes a square matrix in text or binary format. When built
 * with __MPI (the project default), save_mat uses MPI collectives and
 * MPI_Barrier. This test provides a main() that calls MPI_Init first
 * (like test_hsk_writer.cpp) so the MPI path is safe on a single
 * rank. With a default-constructed Parallel_2D (all global2local return
 * -1), rank 0 owns no local elements, so the matrix body is all zeros.
 * Tests verify the text header format (step, dimension, gamma flag)
 * and the binary dim header, not the matrix values themselves.
 */
#include <gtest/gtest.h>

#include "source_io/module_hs/hs_dense_io.h"
#include "source_base/parallel_2d.h"
#include "source_base/parallel_comm.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#ifdef __MPI
#include <mpi.h>
#endif

namespace
{
/**
 * @brief Read the entire content of a text file into a string.
 */
std::string read_text_file(const std::string& filename)
{
    std::ifstream ifs(filename);
    if (!ifs.is_open())
    {
        return "";
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

/**
 * @brief Read just the dimension int from a binary matrix file.
 */
int read_binary_dim(const std::string& filename)
{
    FILE* fp = std::fopen(filename.c_str(), "rb");
    if (fp == nullptr)
    {
        return -1;
    }
    int dim = 0;
    std::fread(&dim, sizeof(int), 1, fp);
    std::fclose(fp);
    return dim;
}

/**
 * @brief Create a serial Parallel_2D that owns all rows/cols on rank 0.
 */
Parallel_2D make_serial_pv(int dim)
{
    Parallel_2D pv;
    pv.set_serial(dim, dim);
    return pv;
}
} // namespace

// ---------------------------------------------------------------------------
// Text mode: verify header fields (MPI path, drank=0)
// ---------------------------------------------------------------------------

TEST(HsDenseIoSaveMat, TextHeaderContainsStepAndDimension)
{
    const std::string filename = "test_save_mat_text_header.dat";
    const int dim = 3;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(0, mat.data(), dim, false, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // MPI text path writes header lines starting with '#'
    EXPECT_NE(content.find("ionic step"), std::string::npos);
    EXPECT_NE(content.find("1"), std::string::npos); // istep+1 = 1
    EXPECT_NE(content.find("gamma only"), std::string::npos);
    EXPECT_NE(content.find("rows"), std::string::npos);
    EXPECT_NE(content.find("3"), std::string::npos); // dim
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, TextHeaderStepTwo)
{
    const std::string filename = "test_save_mat_text_header_step2.dat";
    const int dim = 2;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(1, mat.data(), dim, false, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // istep=1 -> "ionic step 2"
    EXPECT_NE(content.find("ionic step 2"), std::string::npos);
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, TextHeaderGammaOnlyDouble)
{
    const std::string filename = "test_save_mat_text_gamma.dat";
    const int dim = 2;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(0, mat.data(), dim, false, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // gamma_only is true for double (std::is_same<T, double>::value)
    EXPECT_NE(content.find("gamma only 1"), std::string::npos);
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, TextHeaderGammaOnlyComplex)
{
    const std::string filename = "test_save_mat_text_gamma_complex.dat";
    const int dim = 2;
    const std::vector<std::complex<double>> mat(dim * dim);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<std::complex<double>>(0, mat.data(), dim, false, 8,
                               false, false, filename, pv, 0, "genelpa", false);
    const std::string content = read_text_file(filename);
    // gamma_only is false for complex<double>
    EXPECT_NE(content.find("gamma only 0"), std::string::npos);
    std::remove(filename.c_str());
}

// ---------------------------------------------------------------------------
// Binary mode: verify dim header (MPI path, drank=0)
// ---------------------------------------------------------------------------

TEST(HsDenseIoSaveMat, BinaryDimHeaderDouble)
{
    const std::string filename = "test_save_mat_bin_dim_double.dat";
    const int dim = 4;
    const std::vector<double> mat(dim * dim, 0.0);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<double>(0, mat.data(), dim, true, 8, false, false,
                               filename, pv, 0, "genelpa", false);
    const int dim_out = read_binary_dim(filename);
    EXPECT_EQ(dim_out, dim);
    std::remove(filename.c_str());
}

TEST(HsDenseIoSaveMat, BinaryDimHeaderComplex)
{
    const std::string filename = "test_save_mat_bin_dim_complex.dat";
    const int dim = 3;
    const std::vector<std::complex<double>> mat(dim * dim);
    const Parallel_2D pv = make_serial_pv(dim);
    ModuleIO::save_mat<std::complex<double>>(0, mat.data(), dim, true, 8,
                               false, false, filename, pv, 0, "genelpa", false);
    const int dim_out = read_binary_dim(filename);
    EXPECT_EQ(dim_out, dim);
    std::remove(filename.c_str());
}

// ---------------------------------------------------------------------------
// Main: initialize MPI before running tests (required by save_mat's MPI path)
// ---------------------------------------------------------------------------

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    DIAG_WORLD = MPI_COMM_WORLD;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
