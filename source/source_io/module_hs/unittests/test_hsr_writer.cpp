/**
 * @file test_hsr_writer.cpp
 * @brief Unit tests for the filename-generation helpers in hsr_writer.cpp
 *
 * These functions build spin-dependent CSR/DAT filenames for H(R), S(R),
 * and dH/dR outputs. They are pure string functions with no physics state;
 * each test verifies one naming rule combination (append/overwrite,
 * step present/absent, csr/dat extension).
 */
#include <gtest/gtest.h>

#include "source_io/module_hs/hsr_writer.h"

#include <string>

namespace
{
/**
 * @brief Helper: compare two filenames with a clear failure message.
 */
void expect_fname_eq(const std::string& expected, const std::string& actual)
{
    EXPECT_EQ(expected, actual) << "Expected \"" << expected << "\", got \"" << actual << "\"";
}
} // namespace

// ---------------------------------------------------------------------------
// hsr_gen_fname: H(R) spin-dependent filename
// ---------------------------------------------------------------------------

TEST(HsrWriterFname, HsrOverwriteWithStep)
{
    // prefix="data-", ispin=0, append=false, istep=0
    // -> "data-1g1_nao.csr"
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, 0);
    expect_fname_eq("data-1g1_nao.csr", fname);
}

TEST(HsrWriterFname, HsrOverwriteWithStepTwo)
{
    // istep=2 -> g3
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 1, false, 2);
    expect_fname_eq("data-2g3_nao.csr", fname);
}

TEST(HsrWriterFname, HsrAppendSuppressesStep)
{
    // append=true -> step number is omitted
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, true, 5);
    expect_fname_eq("data-1_nao.csr", fname);
}

TEST(HsrWriterFname, HsrNegativeStepSuppressesStep)
{
    // istep=-1 -> step number is omitted (no "g" suffix)
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, -1);
    expect_fname_eq("data-1_nao.csr", fname);
}

TEST(HsrWriterFname, HsrDatExtension)
{
    // out_type=2 -> ".dat" extension
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, 0, 2);
    expect_fname_eq("data-1g1_nao.dat", fname);
}

TEST(HsrWriterFname, HsrCsrExtensionDefault)
{
    // out_type=1 (default) -> ".csr" extension
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, 0, 1);
    expect_fname_eq("data-1g1_nao.csr", fname);
}

// ---------------------------------------------------------------------------
// sr_gen_fname: S(R) spin-independent filename
// ---------------------------------------------------------------------------

TEST(HsrWriterFname, SrOverwriteWithStep)
{
    // append=false, istep=0 -> "srg1_nao.csr"
    const std::string fname = ModuleIO::sr_gen_fname(false, 0);
    expect_fname_eq("srg1_nao.csr", fname);
}

TEST(HsrWriterFname, SrAppendSuppressesStep)
{
    // append=true -> "sr_nao.csr" (no step number)
    const std::string fname = ModuleIO::sr_gen_fname(true, 5);
    expect_fname_eq("sr_nao.csr", fname);
}

TEST(HsrWriterFname, SrNegativeStepSuppressesStep)
{
    // istep=-1 -> "sr_nao.csr"
    const std::string fname = ModuleIO::sr_gen_fname(false, -1);
    expect_fname_eq("sr_nao.csr", fname);
}

TEST(HsrWriterFname, SrDatExtension)
{
    // out_type=2 -> ".dat"
    const std::string fname = ModuleIO::sr_gen_fname(false, 0, 2);
    expect_fname_eq("srg1_nao.dat", fname);
}

// ---------------------------------------------------------------------------
// dhr_gen_fname: dH/dR spin-dependent filename
// ---------------------------------------------------------------------------

TEST(HsrWriterFname, DhrOverwriteWithStep)
{
    // prefix="H", ispin=0, append=false, istep=0
    // -> "Hrs1g1_nao.csr"
    const std::string fname = ModuleIO::dhr_gen_fname("H", 0, false, 0);
    expect_fname_eq("Hrs1g1_nao.csr", fname);
}

TEST(HsrWriterFname, DhrAppendSuppressesStep)
{
    // append=true -> step omitted
    const std::string fname = ModuleIO::dhr_gen_fname("H", 1, true, 3);
    expect_fname_eq("Hrs2_nao.csr", fname);
}

TEST(HsrWriterFname, DhrNegativeStepSuppressesStep)
{
    // istep=-1 -> step omitted
    const std::string fname = ModuleIO::dhr_gen_fname("dH", 0, false, -1);
    expect_fname_eq("dHrs1_nao.csr", fname);
}
