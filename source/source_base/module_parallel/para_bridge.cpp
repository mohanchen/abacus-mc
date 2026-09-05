#include "para_bridge.h"
#include "para_tag.h"

#ifdef __MPI
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#endif

namespace Parallel
{
namespace
{
#ifdef __MPI
// Number of band-parallel groups, derived from the pool layout to avoid a
// dependency on the INPUT-parameter module from source_base.
// nproc = (KPAR * bndpar) * NPROC_IN_POOL, so bndpar = NPROC / (KPAR *
// NPROC_IN_POOL). Falls back to 1 if the globals are inconsistent.
int bndpar_from_layout()
{
    const int denom = GlobalV::KPAR * GlobalV::NPROC_IN_POOL;
    if (denom < 1 || GlobalV::NPROC % denom != 0)
    {
        return 1;
    }
    return GlobalV::NPROC / denom;
}
#endif
} // namespace

// Temporary bridge: construct a pw-domain ParaWorld from the old globals.
// Delete this file once ParaCollection is wired into driver initialization.
ParaWorld make_pw_world()
{
#ifdef __MPI
    return ParaWorld::make_mpi(ParaTag::pw, POOL_WORLD);
#else
    return ParaWorld::serial(ParaTag::pw);
#endif
}

// Reduce-only overload: no k-point distribution data.
ParaKmeshWorld make_kmesh_world()
{
#ifdef __MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    const int bndpar = bndpar_from_layout();
    if (mpi_initialized && GlobalV::KPAR * bndpar > 1)
    {
        // Occupation/energy partial sums are distributed across both k-point
        // pools and band groups, so the reduce domain must span all of
        // MPI_COMM_WORLD (KP_WORLD alone only links same-band ranks across
        // k-point pools and would drop the band-parallel contributions).
        // Build from globals but skip distribute_kpoints (nkstot=0).
        return ParaKmeshWorld(MPI_COMM_WORLD, GlobalV::KPAR, GlobalV::MY_POOL,
                             0, 1, bndpar);
    }
#endif
    return ParaKmeshWorld();
}

// Temporary bridge: construct a kmesh-domain ParaKmeshWorld from the old
// globals. Delete this file once ParaCollection is wired into driver init.
ParaKmeshWorld make_kmesh_world(int nkstot, int nspin)
{
#ifdef __MPI
    // Fall back to a serial single-pool domain when MPI is not initialized
    // (e.g. unit tests linked against the MPI-compiled base library) or when
    // there is only one distribution pool, so that no MPI call is made on
    // an unset communicator.
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    const int bndpar = bndpar_from_layout();
    if (mpi_initialized && GlobalV::KPAR * bndpar > 1)
    {
        return ParaKmeshWorld(MPI_COMM_WORLD, GlobalV::KPAR, GlobalV::MY_POOL,
                             nkstot, nspin, bndpar);
    }
#endif
    return ParaKmeshWorld(nkstot, nspin);
}

// Temporary bridge: construct a bgroup-domain ParaBgroupWorld from the old
// globals. Delete this file once ParaCollection is wired into driver init.
ParaBgroupWorld make_bgroup_world()
{
#ifdef __MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    // NPROC_IN_BNDGROUP stays 0 until divide_pools has run, which also
    // guards unit tests that link the MPI base library without a layout.
    if (mpi_initialized && INT_BGROUP != MPI_COMM_NULL && BP_WORLD != MPI_COMM_NULL
        && GlobalV::NPROC_IN_BNDGROUP > 1)
    {
        int nbndgroup = 1;
        MPI_Comm_size(BP_WORLD, &nbndgroup);
        return ParaBgroupWorld(INT_BGROUP, BP_WORLD, nbndgroup);
    }
#endif
    return ParaBgroupWorld();
}

} // namespace Parallel
