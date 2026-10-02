#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../relax_driver.h"
#include "source_base/matrix.h"
#include "source_io/module_parameter/input_parameter.h"

// Test build_stru_header for stru_out / final_out header generation.
// The header logic is the only part of Relax_Driver that does not depend on
// PARAM.globalv or MPI rank, so it can be unit-tested in isolation.

class RelaxDriverHeaderTest : public testing::Test
{
protected:
    Input_para inp;
    ModuleBase::matrix stress;
    const double etot = -10.5; // Ry

    void SetUp() override
    {
        stress.create(3, 3);
        stress(0, 0) = 1.0; stress(0, 1) = 0.1; stress(0, 2) = 0.2;
        stress(1, 0) = 0.1; stress(1, 1) = 1.0; stress(1, 2) = 0.3;
        stress(2, 0) = 0.2; stress(2, 1) = 0.3; stress(2, 2) = 1.0;
    }
};

TEST_F(RelaxDriverHeaderTest, CalForceStressTrue)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(0, etot, stress, inp, false, true);

    EXPECT_THAT(header, testing::HasSubstr("# RELAX STEP 1, Energy:"));
    EXPECT_THAT(header, testing::HasSubstr("# Stress (kbar): "));
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("N/A")));
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("Forces not computed")));
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("NOTE: geometry proposed")));
}

TEST_F(RelaxDriverHeaderTest, CalStressFalse)
{
    inp.cal_force = true;
    inp.cal_stress = false;

    const std::string header = Relax_Driver::build_stru_header(0, etot, stress, inp, false, true);

    // Exactly 3 N/A lines for stress
    int count = 0;
    size_t pos = 0;
    while ((pos = header.find("# Stress (kbar): N/A N/A N/A", pos)) != std::string::npos)
    {
        ++count;
        pos += 1;
    }
    EXPECT_EQ(count, 3);
    EXPECT_THAT(header, testing::HasSubstr("(N/A = not computed, cal_stress=0)"));
    // Force is computed, so no force note
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("Forces not computed")));
}

TEST_F(RelaxDriverHeaderTest, CalForceFalse)
{
    inp.cal_force = false;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(0, etot, stress, inp, false, true);

    EXPECT_THAT(header, testing::HasSubstr("# Forces not computed (cal_force=0); per-atom f fields omitted intentionally"));
    EXPECT_THAT(header, testing::HasSubstr("# Stress (kbar): "));
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("N/A")));
}

TEST_F(RelaxDriverHeaderTest, FinalStep)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(4, etot, stress, inp, true, true);

    EXPECT_THAT(header, testing::HasSubstr("# RELAX STEP 5 (FINAL), Energy:"));
}

TEST_F(RelaxDriverHeaderTest, GeometryNotEvaluated)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(4, etot, stress, inp, true, false);

    // Stress marked N/A because the proposed geometry was not evaluated
    int count = 0;
    size_t pos = 0;
    while ((pos = header.find("# Stress (kbar): N/A N/A N/A", pos)) != std::string::npos)
    {
        ++count;
        pos += 1;
    }
    EXPECT_EQ(count, 3);
    EXPECT_THAT(header, testing::HasSubstr("(N/A = proposed geometry not evaluated)"));
    EXPECT_THAT(header, testing::HasSubstr("# NOTE: geometry proposed by optimizer but not evaluated"));
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("Forces not computed")));
}
