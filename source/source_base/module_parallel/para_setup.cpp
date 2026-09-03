#include "para_setup.h"

#include <cassert>
#include <iostream>

#ifdef __MPI
#include <mpi.h>
#endif

namespace Parallel
{

void divide_mpi_groups(int nproc, int num_groups, int rank, bool even,
                       int& procs_in_group, int& my_group, int& rank_in_group)
{
    assert(num_groups > 0);
    assert(nproc >= num_groups);

    procs_in_group = nproc / num_groups;
    int extra_procs = nproc % num_groups;

    if (even && extra_procs != 0)
    {
        std::cerr << "Error: " << nproc << " processes not evenly divisible by "
                  << num_groups << " groups." << std::endl;
        assert(false);
    }

    if (rank < extra_procs * (procs_in_group + 1))
    {
        procs_in_group++;
        my_group = rank / procs_in_group;
        rank_in_group = rank % procs_in_group;
    }
    else
    {
        my_group = (rank - extra_procs) / procs_in_group;
        rank_in_group = (rank - extra_procs) % procs_in_group;
    }
}

#ifdef __MPI

namespace {

// Helper: split a communicator into ngroup sub-communicators
struct GroupSplitResult
{
    MPI_Comm group_comm;       // intra-group communicator
    MPI_Comm inter_comm;       // inter-group communicator (MPI_COMM_NULL if single group)
    int ngroups;
    int nprocs_in_group;
    int my_group;
    int rank_in_group;
};

GroupSplitResult split_comm_group(MPI_Comm parent, int ngroup, bool even)
{
    GroupSplitResult res;
    res.group_comm = MPI_COMM_NULL;
    res.inter_comm = MPI_COMM_NULL;
    res.ngroups = ngroup;

    int gsize = 0, grank = 0;
    MPI_Comm_size(parent, &gsize);
    MPI_Comm_rank(parent, &grank);

    divide_mpi_groups(gsize, ngroup, grank, even,
                      res.nprocs_in_group, res.my_group, res.rank_in_group);

    MPI_Comm_split(parent, res.my_group, res.rank_in_group, &res.group_comm);

    if (ngroup > 1)
    {
        MPI_Comm_split(parent, res.rank_in_group, res.my_group, &res.inter_comm);
    }

    return res;
}

} // anonymous namespace

void split_pools(int nproc, int my_rank, int bndpar, int kpar,
                 int& nproc_in_pool, int& rank_in_pool, int& my_pool,
                 int& nproc_in_bndgroup, int& rank_in_bpgroup, int& my_bndgroup,
                 ParaWorld& pw_world, ParaWorld& kmesh_world,
                 ParaWorld& bgroup_int, ParaWorld& bgroup_bp)
{
    if (bndpar > 1 && nproc % (bndpar * kpar) != 0)
    {
        std::cerr << "Error: NPROC (" << nproc
                  << ") must be divisible by BNDPAR*KPAR ("
                  << bndpar * kpar << ")." << std::endl;
        assert(false);
    }

    // k-point parallelization: split WORLD into kpar pools
    GroupSplitResult kpar_res = split_comm_group(MPI_COMM_WORLD, kpar, false);

    // band parallelization: split each pool into bndpar groups
    GroupSplitResult bndpar_res = split_comm_group(kpar_res.group_comm, bndpar, true);

    // Set output indices
    nproc_in_pool = bndpar_res.nprocs_in_group;
    rank_in_pool = bndpar_res.rank_in_group;
    my_pool = kpar_res.my_group;

    // POOL_WORLD
    MPI_Comm pool_comm;
    MPI_Comm_dup(bndpar_res.group_comm, &pool_comm);
    pw_world = ParaWorld::make_mpi(ParaTag::pw, pool_comm);

    // KP_WORLD (inter-pool communicator)
    if (kpar_res.inter_comm != MPI_COMM_NULL)
    {
        MPI_Comm kp_comm;
        MPI_Comm_dup(kpar_res.inter_comm, &kp_comm);
        kmesh_world = ParaWorld::make_mpi(ParaTag::kmesh, kp_comm);
    }
    else
    {
        kmesh_world = ParaWorld::make_mpi(ParaTag::kmesh, MPI_COMM_NULL);
    }

    // Band group communicators
    if (bndpar > 1)
    {
        nproc_in_bndgroup = kpar_res.ngroups * bndpar_res.nprocs_in_group;
        rank_in_bpgroup = kpar_res.my_group * bndpar_res.nprocs_in_group + bndpar_res.rank_in_group;
        my_bndgroup = bndpar_res.my_group;

        MPI_Comm int_bgroup;
        MPI_Comm_split(MPI_COMM_WORLD, my_bndgroup, rank_in_bpgroup, &int_bgroup);
        bgroup_int = ParaWorld::make_mpi(ParaTag::bsame_kdiff, int_bgroup);

        MPI_Comm bp_comm;
        MPI_Comm_dup(bndpar_res.inter_comm, &bp_comm);
        bgroup_bp = ParaWorld::make_mpi(ParaTag::bdiff_ksame, bp_comm);
    }
    else
    {
        nproc_in_bndgroup = nproc;
        rank_in_bpgroup = my_rank;
        my_bndgroup = 0;

        MPI_Comm int_bgroup;
        MPI_Comm_dup(MPI_COMM_WORLD, &int_bgroup);
        bgroup_int = ParaWorld::make_mpi(ParaTag::bsame_kdiff, int_bgroup);

        MPI_Comm bp_comm;
        MPI_Comm_split(MPI_COMM_WORLD, my_rank, 0, &bp_comm);
        bgroup_bp = ParaWorld::make_mpi(ParaTag::bdiff_ksame, bp_comm);
    }
}

ParaWorld split_diag_world(int diag_np, int nproc, int my_rank,
                           int& drank, int& dsize, int& dcolor)
{
    assert(diag_np > 0);

    int procs_in_group = 0, my_group = 0, rank_in_group = 0;
    divide_mpi_groups(nproc, diag_np, my_rank, false,
                      procs_in_group, my_group, rank_in_group);

    MPI_Comm diag_comm;
    MPI_Comm_split(MPI_COMM_WORLD, my_group, rank_in_group, &diag_comm);

    MPI_Comm_rank(diag_comm, &drank);
    MPI_Comm_size(diag_comm, &dsize);
    dcolor = my_group;

    return ParaWorld::make_mpi(ParaTag::diag, diag_comm);
}

ParaWorld split_grid_world(int diag_np, int nproc, int my_rank,
                           int& grank, int& gsize)
{
    assert(diag_np > 0);

    int procs_in_group = 0, my_group = 0, rank_in_group = 0;
    divide_mpi_groups(nproc, diag_np, my_rank, false,
                      procs_in_group, my_group, rank_in_group);

    MPI_Comm grid_comm;
    MPI_Comm_split(MPI_COMM_WORLD, my_group, rank_in_group, &grid_comm);

    MPI_Comm_rank(grid_comm, &grank);
    MPI_Comm_size(grid_comm, &gsize);

    return ParaWorld::make_mpi(ParaTag::rgrid, grid_comm);
}

ParaCollection setup_para_worlds(int nproc, int my_rank, int bndpar, int kpar, int diag_np)
{
    ParaCollection worlds;

    // 1. POOL_WORLD + KP_WORLD + band group comms
    int nproc_in_pool, rank_in_pool, my_pool;
    int nproc_in_bndgroup, rank_in_bpgroup, my_bndgroup;

    ParaWorld pw_world = ParaWorld::make_mpi(ParaTag::pw, MPI_COMM_NULL);
    ParaWorld kmesh_world = ParaWorld::make_mpi(ParaTag::kmesh, MPI_COMM_NULL);
    ParaWorld bgroup_int = ParaWorld::make_mpi(ParaTag::bsame_kdiff, MPI_COMM_NULL);
    ParaWorld bgroup_bp = ParaWorld::make_mpi(ParaTag::bdiff_ksame, MPI_COMM_NULL);

    split_pools(nproc, my_rank, bndpar, kpar,
                nproc_in_pool, rank_in_pool, my_pool,
                nproc_in_bndgroup, rank_in_bpgroup, my_bndgroup,
                pw_world, kmesh_world, bgroup_int, bgroup_bp);

    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::pw, pw_world.comm()));
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::kmesh, kmesh_world.comm()));
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::bsame_kdiff, bgroup_int.comm()));
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::bdiff_ksame, bgroup_bp.comm()));

    // 2. DIAG_WORLD
    int drank, dsize, dcolor;
    ParaWorld diag_world = split_diag_world(diag_np, nproc, my_rank, drank, dsize, dcolor);
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::diag, diag_world.comm()));

    // 3. GRID_WORLD
    int grank, gsize;
    ParaWorld grid_world = split_grid_world(diag_np, nproc, my_rank, grank, gsize);
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::rgrid, grid_world.comm()));

    // 4. matrix domain (serial for now, will get its own split later)
    worlds.add(ParaWorld::make_serial(ParaTag::matrix));

    return worlds;
}

#endif // __MPI

} // namespace Parallel
