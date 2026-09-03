#include "gtest/gtest.h"

#include "../para_rgrid_world.h"

TEST(ParaRgridWorldTest, SerialMode)
{
    const Parallel::ParaRgridWorld world(4, 5, 6);
    EXPECT_EQ(world.tag(), "rgrid");
    EXPECT_EQ(world.ncx(), 4);
    EXPECT_EQ(world.ncy(), 5);
    EXPECT_EQ(world.ncz(), 6);
    EXPECT_EQ(world.nczp(), 6);
    EXPECT_EQ(world.nrxx(), 120);
    EXPECT_EQ(world.numz(0), 6);
    EXPECT_EQ(world.startz(0), 0);
}

TEST(ParaRgridWorldTest, WhichproSerial)
{
    const Parallel::ParaRgridWorld world(2, 2, 8);
    for (int iz = 0; iz < 8; ++iz)
    {
        EXPECT_EQ(world.whichpro(iz), 0);
    }
}

TEST(ParaRgridWorldTest, Validity)
{
    const Parallel::ParaRgridWorld world(1, 1, 1);
    EXPECT_TRUE(world.valid());
}
