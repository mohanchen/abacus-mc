#include "gtest/gtest.h"

#include "../para_kmesh_world.h"

// Run with: mpirun -np 4 ./MODULE_BASE_para_kmesh_world_mpi
//
// Covers the band-parallel reduction path (bndpar > 1) that the legacy
// Parallel_Reduce::reduce_double_allpool provided and that was dropped in the
// first ParaKmeshWorld migration, breaking tests/11_PW_GPU/scf_bpcg
// (kpar=1, bndpar=2).

TEST(ParaKmeshWorldMpiTest, ReduceAcrossBandGroupsBndpar2)
{
    int nprocs = 0;
    int myrank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    ASSERT_EQ(nprocs, 4);

    // kpar=1, bndpar=2, 4 ranks: 2 ranks per band group. Each band group
    // holds a partial occupation sum of 14 (28 electrons split in two).
    Parallel::ParaKmeshWorld kmesh(MPI_COMM_WORLD, 1, 0, 0, 1, 2);
    EXPECT_EQ(kmesh.npool(), 2);

    double sumk = 14.0;
    kmesh.reduce_across_pools(sumk);
    EXPECT_DOUBLE_EQ(sumk, 28.0);

    // max/min must span the band groups as well: the two shards see
    // different eigenvalue windows.
    double eup = (myrank < 2) ? 40.0 : 45.0;
    kmesh.reduce_max_across_pools(eup);
    EXPECT_DOUBLE_EQ(eup, 45.0);

    double elw = (myrank < 2) ? -1.0 : -5.0;
    kmesh.reduce_min_across_pools(elw);
    EXPECT_DOUBLE_EQ(elw, -5.0);
}

TEST(ParaKmeshWorldMpiTest, ReduceAcrossKpoolsKpar2)
{
    int nprocs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    // kpar=2, bndpar=1: 2 k-point pools of 2 ranks each.
    const int my_pool = (nprocs > 1) ? 0 : 0; // distribution detail unused here
    Parallel::ParaKmeshWorld kmesh(MPI_COMM_WORLD, 2, my_pool, 4, 1, 1);
    EXPECT_EQ(kmesh.npool(), 2);

    double sumk = 3.5;
    kmesh.reduce_across_pools(sumk);
    EXPECT_DOUBLE_EQ(sumk, 7.0);
}

TEST(ParaKmeshWorldMpiTest, SinglePoolIsNoOp)
{
    int nprocs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    // npool == 1: reduction must be a no-op regardless of the world size.
    Parallel::ParaKmeshWorld kmesh(MPI_COMM_WORLD, 1, 0, 0, 1, 1);
    double sumk = 42.0;
    kmesh.reduce_across_pools(sumk);
    EXPECT_DOUBLE_EQ(sumk, 42.0);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
