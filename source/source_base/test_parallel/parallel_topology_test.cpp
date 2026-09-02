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
    EXPECT_EQ(t.my_band_group(), 0);
    EXPECT_EQ(t.rank_in_band_group(), 0);
    EXPECT_EQ(t.nproc_in_band_group(), 1);

    EXPECT_EQ(t.pool_root_rank(0), 0);
    EXPECT_EQ(t.pool_root_rank(-1), -1);
    EXPECT_EQ(t.pool_root_rank(1), -1);
    EXPECT_EQ(t.band_group_root_rank(0), 0);
    EXPECT_EQ(t.band_group_root_rank(-1), -1);
    EXPECT_EQ(t.band_group_root_rank(1), -1);

#ifdef __MPI
    // Serial fallback: matrix/atom domains are not yet bound.
    EXPECT_EQ(t.matrix_world_comm(), MPI_COMM_NULL);
    EXPECT_EQ(t.atom_world_comm(), MPI_COMM_NULL);
#endif
}

TEST(ProcessTopology, ConstructAndAccessors)
{
    // Real ABACUS layout with 12 world ranks 0..11, KPAR=3, BNDPAR=2.
    // divide_pools constraint (BNDPAR>1): NPROC % (BNDPAR*KPAR) == 0,
    // i.e. 12 % 6 == 0, so pools are even: {4,4,4} covering ranks
    //   pool0: 0..3, pool1: 4..7, pool2: 8..11
    // Each pool is split by BNDPAR=2 (even) into 2 slices of
    // pool_size/bndpar = 2 ranks:
    //   bg0: {0,1} {4,5} {8,9}   -> 6 ranks, root = world rank 0
    //   bg1: {2,3} {6,7} {10,11} -> 6 ranks, root = world rank 2
    // Global band-group size nproc_in_band_group = kpar * 2 = 6 and
    // bndpar * nproc_in_band_group = 12 = world.
    // my_rank = 6 sits in pool1 (start 4) at rank_in_pool 2, which
    // belongs to bg1 (pool1 slice [2,4)) at global band-group rank
    // RANK_IN_BPGROUP = my_pool * 2 + 0 = 2.
    const std::vector<int> pool_sizes = {4, 4, 4};
#ifdef __MPI
    const ProcessTopology t(/*world=*/12,
                           /*my_rank=*/6,
                           /*kpar=*/3,
                           /*my_pool=*/1,
                           /*rank_in_pool=*/2,
                           pool_sizes,
                           /*bndpar=*/2,
                           /*my_band_group=*/1,
                           /*rank_in_band_group=*/2,
                           /*nproc_in_band_group=*/6,
                           /*pw_world=*/MPI_COMM_SELF,
                           /*kmesh_world=*/MPI_COMM_NULL,
                           /*bsame_kdiff=*/MPI_COMM_SELF,
                           /*bdiff_ksame=*/MPI_COMM_NULL,
                           /*rgrid_world=*/MPI_COMM_NULL,
                           /*diag_world=*/MPI_COMM_NULL,
                           /*matrix_world=*/MPI_COMM_SELF,
                           /*atom_world=*/MPI_COMM_WORLD);
#else
    const ProcessTopology t(12, 6, 3, 1, 2, pool_sizes, 2, 1, 2, 6);
#endif

    EXPECT_EQ(t.world_size(), 12);
    EXPECT_EQ(t.world_rank(), 6);

    EXPECT_EQ(t.kpar(), 3);
    EXPECT_EQ(t.my_pool(), 1);
    EXPECT_EQ(t.rank_in_pool(), 2);
    EXPECT_EQ(t.nproc_in_pool(0), 4);
    EXPECT_EQ(t.nproc_in_pool(1), 4);
    EXPECT_EQ(t.nproc_in_pool(2), 4);

    EXPECT_EQ(t.bndpar(), 2);
    EXPECT_EQ(t.my_band_group(), 1);
    EXPECT_EQ(t.rank_in_band_group(), 2);
    EXPECT_EQ(t.nproc_in_band_group(), 6);

    // pool_root_rank prefix offsets -> pool start world rank.
    EXPECT_EQ(t.pool_root_rank(0), 0);
    EXPECT_EQ(t.pool_root_rank(1), 4);
    EXPECT_EQ(t.pool_root_rank(2), 8);
    EXPECT_EQ(t.pool_root_rank(3), -1);

    // band_group_root_rank: first member of each band-group union.
    //   bg0 root = 0 * (4/2) = 0
    //   bg1 root = 1 * (4/2) = 2  (world rank 2, first member of bg1)
    EXPECT_EQ(t.band_group_root_rank(0), 0);
    EXPECT_EQ(t.band_group_root_rank(1), 2);
    EXPECT_EQ(t.band_group_root_rank(2), -1);

#ifdef __MPI
    // The 8 `_world_comm` accessors all round-trip what we injected.
    EXPECT_NE(t.pw_world_comm(), MPI_COMM_NULL); // MPI_COMM_SELF is non-null
    EXPECT_EQ(t.kmesh_world_comm(), MPI_COMM_NULL);
    EXPECT_NE(t.bsame_kdiff_world_comm(), MPI_COMM_NULL);
    EXPECT_EQ(t.bdiff_ksame_world_comm(), MPI_COMM_NULL);
    EXPECT_EQ(t.rgrid_world_comm(), MPI_COMM_NULL);
    EXPECT_EQ(t.diag_world_comm(), MPI_COMM_NULL);
    EXPECT_NE(t.matrix_world_comm(), MPI_COMM_NULL); // caller-bound this time
    EXPECT_EQ(t.atom_world_comm(), MPI_COMM_WORLD);
#endif

    // Copyable: copy should behave identically.
    ProcessTopology copy = t;
    EXPECT_EQ(copy.world_rank(), t.world_rank());
    EXPECT_EQ(copy.pool_root_rank(2), t.pool_root_rank(2));
    EXPECT_EQ(copy.band_group_root_rank(1), t.band_group_root_rank(1));
    EXPECT_EQ(copy.nproc_in_pool(1), t.nproc_in_pool(1));
    EXPECT_EQ(copy.nproc_in_band_group(), t.nproc_in_band_group());
#ifdef __MPI
    EXPECT_EQ(copy.atom_world_comm(), t.atom_world_comm());
    EXPECT_EQ(copy.matrix_world_comm(), t.matrix_world_comm());
#endif
}

// Integration: Parallel_Global::create_topology builds a real topology
// on a live MPI_COMM_WORLD. This test is only meaningful when run via
// mpirun with exactly 4 ranks.
TEST(ParallelGlobalCreateTopology, FourRanksKpar2Bndpar2DiagNp2)
{
    int world_size = 1;
    int world_rank = 0;
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
    if (world_size != 4)
    {
        GTEST_SKIP() << "create_topology integration case requires exactly 4 MPI ranks (got "
                     << world_size << ").";
    }
    const int kpar = 2;
    const int bndpar = 2;
    const int diag_np = 2;
    const int grid_np = 2; // reserved; topology today derives rgrid from diag_np
    const ProcessTopology t = Parallel_Global::create_topology(
        world_size, world_rank, kpar, bndpar, diag_np, grid_np);

    // ---- scalar invariants (same for every rank) -------------------
    EXPECT_EQ(t.world_size(), 4);
    EXPECT_EQ(t.world_rank(), world_rank);
    EXPECT_EQ(t.kpar(), 2);
    EXPECT_EQ(t.bndpar(), 2);
    EXPECT_EQ(t.bndpar() * t.nproc_in_band_group(), t.world_size());

    // 4 ranks, KPAR=2 (even=false): nproc_in_pool = {2, 2}
    ASSERT_EQ(static_cast<int>(t.nproc_in_pool().size()), 2);
    EXPECT_EQ(t.nproc_in_pool(0), 2);
    EXPECT_EQ(t.nproc_in_pool(1), 2);
    EXPECT_EQ(t.pool_root_rank(0), 0);
    EXPECT_EQ(t.pool_root_rank(1), 2);

    // BNDPAR=2 (even=true) in each pool of 2 procs -> each bg gets 1 proc
    // Global union bg size = kpar*1 = 2 = nproc_in_band_group
    EXPECT_EQ(t.nproc_in_band_group(), 2);
    // bg0 = {0, 2} -> root 0; bg1 = {1, 3} -> root = 1 * (2/2) = 1
    EXPECT_EQ(t.band_group_root_rank(0), 0);
    EXPECT_EQ(t.band_group_root_rank(1), 1);

    // ---- per-rank local values -------------------------------------
    // World ranks 0..3, KPAR=2 (even=false) -> pools {0: {0,1}, 1: {2,3}}.
    // BNDPAR=2 even split inside each pool -> for a pool of size 2,
    // each band-group receives exactly 1 process (nprocs_in_group =
    // pool_size / bndpar = 1). Therefore within every pool the
    // bndpar-group rank_in_group is always 0 for whichever
    // band-group the rank belongs to, and NPROC_IN_POOL (aka the
    // band-group slice size per pool) equals 1 globally across all
    // 4 ranks.
    //
    // Per-rank expected (verified by factory stdout against real
    // divide_pools run on R0..R3):
    //   wr | pool | rank_in_pool | band_group | rank_in_band_group
    //    0 |   0  |       0      |     0      | pool*1 + 0 = 0
    //    1 |   0  |       0      |     1      | pool*1 + 0 = 0
    //    2 |   1  |       0      |     0      | pool*1 + 0 = 1
    //    3 |   1  |       0      |     1      | pool*1 + 0 = 1
    const int expected_pool = (world_rank < 2) ? 0 : 1;
    const int expected_band_group = world_rank % 2; // bndpar even slices by world order inside pool
    const int expected_rank_in_pool = 0; // 1 process per (pool, band-group) tile
    const int expected_rk_in_bg = expected_pool * 1 + 0; // formula with nproc_per_bg_in_pool = 1

    EXPECT_EQ(t.my_pool(), expected_pool);
    EXPECT_EQ(t.rank_in_pool(), expected_rank_in_pool);
    EXPECT_EQ(t.my_band_group(), expected_band_group);
    EXPECT_EQ(t.rank_in_band_group(), expected_rk_in_bg);

#ifdef __MPI
    // ---- communicator sizes / memberships --------------------------
    int size = -1;
    int rk = -1;
    // pw_world_comm : (same pool, same band group) intersection of
    // 2x2 = 4 total PW tiles -> singleton size 1 each
    ASSERT_NE(t.pw_world_comm(), MPI_COMM_NULL);
    MPI_Comm_size(t.pw_world_comm(), &size);
    MPI_Comm_rank(t.pw_world_comm(), &rk);
    EXPECT_EQ(size, 1);
    EXPECT_EQ(rk, 0);

    // kmesh_world_comm (KP_WORLD inter_comm bridge): non-null only on
    // ranks where rank_in_pool == 0 (pool 0 rank 0 = wr 0; pool 1 rank 0
    // = wr 2) -> size 2 on those ranks, MPI_COMM_NULL otherwise.
    if (t.rank_in_pool() == 0)
    {
        ASSERT_NE(t.kmesh_world_comm(), MPI_COMM_NULL);
        MPI_Comm_size(t.kmesh_world_comm(), &size);
        MPI_Comm_rank(t.kmesh_world_comm(), &rk);
        EXPECT_EQ(size, kpar);
        EXPECT_EQ(rk, t.my_pool()); // key order in inter_comm = my_group
    }
    else
    {
        EXPECT_EQ(t.kmesh_world_comm(), MPI_COMM_NULL);
    }

    // bsame_kdiff_world_comm (INT_BGROUP): same band group across all pools
    // size = kpar * nproc_per_pool_per_bg = 2 * 1 = 2
    ASSERT_NE(t.bsame_kdiff_world_comm(), MPI_COMM_NULL);
    MPI_Comm_size(t.bsame_kdiff_world_comm(), &size);
    MPI_Comm_rank(t.bsame_kdiff_world_comm(), &rk);
    EXPECT_EQ(size, 2);
    EXPECT_EQ(rk, t.rank_in_band_group()); // matches expected_rk_in_bg

    // bdiff_ksame_world_comm (BP_WORLD duped): same pool, same
    // rank_in_pool (i.e. rank_in_band_group value) pairs. With 2 bg in
    // 1 pool of 2 processes, each BP_WORLD subgroup also has size 2.
    ASSERT_NE(t.bdiff_ksame_world_comm(), MPI_COMM_NULL);
    MPI_Comm_size(t.bdiff_ksame_world_comm(), &size);
    MPI_Comm_rank(t.bdiff_ksame_world_comm(), &rk);
    EXPECT_EQ(size, bndpar);
    EXPECT_EQ(rk, t.my_band_group());

    // rgrid_world_comm (GRID_WORLD): diag_np=2 groups, even=false -> 2+2
    ASSERT_NE(t.rgrid_world_comm(), MPI_COMM_NULL);
    MPI_Comm_size(t.rgrid_world_comm(), &size);
    EXPECT_EQ(size, world_size / diag_np); // 4/2 = 2

    // diag_world_comm (DIAG_WORLD): diag_np=2, divide_mpi_groups over
    // nproc=4 groups=2 even=false -> groups {0,1} {2,3}; color = rank's
    // in-group rank. All 4 ranks take part in MPI_Comm_split with
    // their color so there are diag_np different DIAG_WORLD subgroups
    // each of size diag_np = 2 (colors 0 & 1 both have 2 members:
    // color 0 from {wr 0, 2}; color 1 from {wr 1, 3}).
    ASSERT_NE(t.diag_world_comm(), MPI_COMM_NULL);
    MPI_Comm_size(t.diag_world_comm(), &size);
    MPI_Comm_rank(t.diag_world_comm(), &rk);
    EXPECT_EQ(size, diag_np); // color 0 and 1 each 2 members
    EXPECT_GE(rk, 0);
    EXPECT_LT(rk, size);

    // matrix / atom domains are caller-bound today -> must be null.
    EXPECT_EQ(t.matrix_world_comm(), MPI_COMM_NULL);
    EXPECT_EQ(t.atom_world_comm(), MPI_COMM_NULL);
#endif
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int rc = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return rc;
}
#else

// Top-level test binary used when __MPI is not defined: fall back to
// googletest's default main via link target. This branch is never
// compiled in the normal unit-test build path because
// MODULE_BASE_ProcessTopology links only under the __MPI build
// configuration.
#endif // __MPI
