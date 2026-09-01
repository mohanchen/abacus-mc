#ifdef __MPI

#include "source_base/parallel_global.h"
#include "source_base/parallel_topology.h"

#include "mpi.h"

#include "gtest/gtest.h"
#include <numeric>
#include <vector>

/************************************************
 *  Unit tests for:
 *   - Parallel_Global::divide_mpi_groups (pure arithmetic)
 *   - ProcessTopology construction and accessors
 ***********************************************/

namespace
{

/// Run divide_mpi_groups for every rank in [0, procs) with the same
/// inputs and return the per-rank results as vectors. Asserts that
/// (1) every rank's (group, rank_in_group) pair is unique inside
/// its group and (2) the group sizes match procs_in_group counts.
struct DivideResult
{
    std::vector<int> procs_in_group; // size = num_groups, per-call value (same for every rank)
    std::vector<int> my_group;       // [0..procs-1] -> group of rank r
    std::vector<int> rank_in_group;  // [0..procs-1] -> in-group rank of r
};

DivideResult divide_all(const int procs, const int num_groups, const bool even = false)
{
    DivideResult out;
    out.my_group.assign(procs, -1);
    out.rank_in_group.assign(procs, -1);
    for (int r = 0; r < procs; ++r)
    {
        int p_in_g = 0;
        int mg = -1;
        int r_in_g = -1;
        Parallel_Global::divide_mpi_groups(procs, num_groups, r, p_in_g, mg, r_in_g, even);
        // divide_mpi_groups mutates procs_in_group on output for the current
        // rank: for ranks in a "plus-one" group it returns base+1, otherwise
        // base. The group-size vector is rank-independent.
        if (static_cast<int>(out.procs_in_group.size()) != num_groups)
        {
            out.procs_in_group.assign(num_groups, 0);
            const int base = procs / num_groups;
            const int extra = procs % num_groups;
            for (int g = 0; g < num_groups; ++g)
            {
                out.procs_in_group[g] = base + (g < extra ? 1 : 0);
            }
        }
        out.my_group[r] = mg;
        out.rank_in_group[r] = r_in_g;
        EXPECT_EQ(p_in_g, out.procs_in_group[mg]);
    }
    return out;
}

void validate_divide(const int procs, const int num_groups, const DivideResult& out)
{
    ASSERT_EQ(static_cast<int>(out.procs_in_group.size()), num_groups);
    ASSERT_EQ(std::accumulate(out.procs_in_group.begin(), out.procs_in_group.end(), 0), procs);

    // per group: collect (rank_in_group -> global rank) mapping
    std::vector<std::vector<int>> seen(num_groups);
    for (int g = 0; g < num_groups; ++g)
    {
        seen[g].assign(out.procs_in_group[g], -1);
    }
    for (int r = 0; r < procs; ++r)
    {
        const int g = out.my_group[r];
        ASSERT_GE(g, 0);
        ASSERT_LT(g, num_groups);
        const int rg = out.rank_in_group[r];
        ASSERT_GE(rg, 0);
        ASSERT_LT(rg, out.procs_in_group[g]);
        EXPECT_EQ(seen[g][rg], -1) << "duplicate in-group rank " << rg << " in group " << g;
        seen[g][rg] = r;
    }
}

} // namespace

// ---- divide_mpi_groups tests (pure arithmetic, callable in MPI_Init'd process) ---

TEST(DivideMpiGroups, EvenDivision)
{
    const int procs = 8;
    const int kpar = 4;
    const auto out = divide_all(procs, kpar, /*even=*/true);
    validate_divide(procs, kpar, out);
    for (int g = 0; g < kpar; ++g)
    {
        EXPECT_EQ(out.procs_in_group[g], 2);
    }
    // First 2 ranks in group 0: r=0 -> group 0, rank 0; r=1 -> group 0, rank 1
    EXPECT_EQ(out.my_group[0], 0);
    EXPECT_EQ(out.rank_in_group[0], 0);
    EXPECT_EQ(out.my_group[1], 0);
    EXPECT_EQ(out.rank_in_group[1], 1);
    // Last 2 ranks -> group 3
    EXPECT_EQ(out.my_group[6], 3);
    EXPECT_EQ(out.rank_in_group[6], 0);
    EXPECT_EQ(out.my_group[7], 3);
    EXPECT_EQ(out.rank_in_group[7], 1);
}

TEST(DivideMpiGroups, UnevenDivision)
{
    // 5 procs / 2 groups -> first group 3, second 2.
    const int procs = 5;
    const int kpar = 2;
    const auto out = divide_all(procs, kpar, /*even=*/false);
    validate_divide(procs, kpar, out);
    EXPECT_EQ(out.procs_in_group[0], 3);
    EXPECT_EQ(out.procs_in_group[1], 2);
    EXPECT_EQ(out.my_group[2], 0);
    EXPECT_EQ(out.rank_in_group[2], 2);
    EXPECT_EQ(out.my_group[3], 1);
    EXPECT_EQ(out.rank_in_group[3], 0);
    EXPECT_EQ(out.my_group[4], 1);
    EXPECT_EQ(out.rank_in_group[4], 1);
}

TEST(DivideMpiGroups, ExactlyOneProcessPerGroup)
{
    const int procs = 6;
    const auto out = divide_all(procs, procs);
    validate_divide(procs, procs, out);
    for (int r = 0; r < procs; ++r)
    {
        EXPECT_EQ(out.my_group[r], r);
        EXPECT_EQ(out.rank_in_group[r], 0);
        EXPECT_EQ(out.procs_in_group[r], 1);
    }
}

TEST(DivideMpiGroups, OneGroup)
{
    const int procs = 5;
    const auto out = divide_all(procs, /*num_groups=*/1);
    validate_divide(procs, /*num_groups=*/1, out);
    EXPECT_EQ(out.procs_in_group[0], 5);
    for (int r = 0; r < procs; ++r)
    {
        EXPECT_EQ(out.my_group[r], 0);
        EXPECT_EQ(out.rank_in_group[r], r);
    }
}

TEST(DivideMpiGroups, UnevenLargePools)
{
    // Typical ABACUS 24 proc / 5 pools -> sizes [5,5,5,5,4].
    // Groups 0..3 span ranks [0..5), [5..10), [10..15), [15..20);
    // last group spans ranks [20..24).
    const int procs = 24;
    const int kpar = 5;
    const auto out = divide_all(procs, kpar);
    validate_divide(procs, kpar, out);
    for (int g = 0; g < kpar; ++g)
    {
        EXPECT_EQ(out.procs_in_group[g], g < 4 ? 5 : 4);
    }
    // Rank 14 -> last rank of group 2 -> rank_in_group 4.
    EXPECT_EQ(out.my_group[14], 2);
    EXPECT_EQ(out.rank_in_group[14], 4);
    // Rank 19 -> last rank of group 3 (size 5) -> rank_in_group 4.
    EXPECT_EQ(out.my_group[19], 3);
    EXPECT_EQ(out.rank_in_group[19], 4);
    // Rank 20 -> first rank of last (size=4) group -> group 4, rank 0.
    EXPECT_EQ(out.my_group[20], 4);
    EXPECT_EQ(out.rank_in_group[20], 0);
    // Last rank (23) -> group 4, rank 3.
    EXPECT_EQ(out.my_group[23], 4);
    EXPECT_EQ(out.rank_in_group[23], 3);
}

// ---- ProcessTopology tests (value semantics, no MPI calls in body) ----------------

TEST(ProcessTopology, DefaultConstructorIsSingleProcess)
{
    const ProcessTopology t;
    EXPECT_EQ(t.world_size(), 1);
    EXPECT_EQ(t.world_rank(), 0);

    EXPECT_EQ(t.kpar(), 1);
    EXPECT_EQ(t.my_pool(), 0);
    EXPECT_EQ(t.rank_in_pool(), 0);
    ASSERT_EQ(static_cast<int>(t.nproc_in_pool().size()), 1);
    EXPECT_EQ(t.nproc_in_pool(0), 1);

    EXPECT_EQ(t.bndpar(), 1);
    EXPECT_EQ(t.my_bndgroup(), 0);
    EXPECT_EQ(t.rank_in_bgroup(), 0);
    EXPECT_EQ(t.nproc_in_bgroup(), 1);

    EXPECT_EQ(t.pool_root_rank(0), 0);
    EXPECT_EQ(t.pool_root_rank(-1), -1);
    EXPECT_EQ(t.pool_root_rank(1), -1);
}

TEST(ProcessTopology, ConstructAndAccessors)
{
    // Procs: 10, KPAR=3 -> pool sizes [4,3,3], BNDPAR=2 -> bgroups 5+5
    const std::vector<int> pool_sizes = {4, 3, 3};
#ifdef __MPI
    const ProcessTopology t(/*world=*/10,
                           /*my_rank=*/6,
                           /*kpar=*/3,
                           /*my_pool=*/2,
                           /*rank_in_pool=*/2,
                           pool_sizes,
                           /*bndpar=*/2,
                           /*my_bndgroup=*/1,
                           /*rank_in_bgroup=*/1,
                           /*nproc_in_bgroup=*/5,
                           /*pool_comm=*/MPI_COMM_SELF,
                           /*kp_world_comm=*/MPI_COMM_NULL,
                           /*int_bgroup_comm=*/MPI_COMM_SELF,
                           /*bp_world_comm=*/MPI_COMM_NULL,
                           /*grid_comm=*/MPI_COMM_NULL,
                           /*diag_comm=*/MPI_COMM_NULL);
#else
    const ProcessTopology t(10, 6, 3, 2, 2, pool_sizes, 2, 1, 1, 5);
#endif

    EXPECT_EQ(t.world_size(), 10);
    EXPECT_EQ(t.world_rank(), 6);

    EXPECT_EQ(t.kpar(), 3);
    EXPECT_EQ(t.my_pool(), 2);
    EXPECT_EQ(t.rank_in_pool(), 2);
    EXPECT_EQ(t.nproc_in_pool(0), 4);
    EXPECT_EQ(t.nproc_in_pool(1), 3);
    EXPECT_EQ(t.nproc_in_pool(2), 3);

    EXPECT_EQ(t.bndpar(), 2);
    EXPECT_EQ(t.my_bndgroup(), 1);
    EXPECT_EQ(t.rank_in_bgroup(), 1);
    EXPECT_EQ(t.nproc_in_bgroup(), 5);

    // pool_root_rank: offsets = [0, 4, 7]
    EXPECT_EQ(t.pool_root_rank(0), 0);
    EXPECT_EQ(t.pool_root_rank(1), 4);
    EXPECT_EQ(t.pool_root_rank(2), 7);
    EXPECT_EQ(t.pool_root_rank(3), -1);

    // Copyable: copy should behave identically.
    ProcessTopology copy = t;
    EXPECT_EQ(copy.world_rank(), t.world_rank());
    EXPECT_EQ(copy.pool_root_rank(2), t.pool_root_rank(2));
    EXPECT_EQ(copy.nproc_in_pool(1), t.nproc_in_pool(1));
}

#endif // __MPI
