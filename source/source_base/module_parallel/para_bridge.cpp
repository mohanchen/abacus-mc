#include "para_bridge.h"
#include "para_tag.h"

#ifdef __MPI
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#endif

namespace Parallel
{

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
    if (mpi_initialized && GlobalV::KPAR > 1 && KP_WORLD != MPI_COMM_NULL)
    {
        // Build from globals but skip distribute_kpoints (nkstot=0).
        return ParaKmeshWorld(KP_WORLD, GlobalV::KPAR, GlobalV::MY_POOL,
                             GlobalV::NPROC, 0, 1);
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
    // (e.g. unit tests linked against the MPI-compiled base library) or
    // when there is only one k-point pool, so that no MPI call is made on
    // an unset communicator.
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized && GlobalV::KPAR > 1 && KP_WORLD != MPI_COMM_NULL)
    {
        return ParaKmeshWorld(KP_WORLD, GlobalV::KPAR, GlobalV::MY_POOL,
                             GlobalV::NPROC, nkstot, nspin);
    }
#endif
    return ParaKmeshWorld(nkstot, nspin);
}

} // namespace Parallel
