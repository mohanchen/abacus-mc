#include "source_pw/module_proj/projgen.h"

#include "source_base/math_integral.h"

#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace
{

double radial_norm(const std::vector<double>& r, const std::vector<double>& radial)
{
    std::vector<double> dr(radial.size());
    std::vector<double> integrand(radial.size());
    std::adjacent_difference(r.begin(), r.begin() + radial.size(), dr.begin());
    for (std::size_t i = 0; i < radial.size(); ++i)
    {
        integrand[i] = r[i] * r[i] * radial[i] * radial[i];
    }
    return ModuleBase::Integral::simpson(radial.size(), integrand.data(), &dr[1]);
}

} // namespace

TEST(ProjgenTest, GeneratesNormalizedTruncatedProjector)
{
    const int nr = 201;
    const double dr = 0.05;
    const double rcut = 5.0;
    std::vector<double> r(nr);
    std::vector<double> chi(nr);
    for (int i = 0; i < nr; ++i)
    {
        r[i] = i * dr;
        chi[i] = std::exp(-r[i]);
    }

    std::vector<double> alpha;
    projgen(0, nr, r.data(), chi.data(), rcut, 4, alpha);

    ASSERT_EQ(alpha.size(), 101U);
    EXPECT_TRUE(std::all_of(alpha.begin(), alpha.end(), [](const double value) { return std::isfinite(value); }));
    EXPECT_NEAR(radial_norm(r, alpha), 1.0, 1.0e-10);
}

