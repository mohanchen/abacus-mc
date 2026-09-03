#ifndef PARA_SETUP_H
#define PARA_SETUP_H

#include <vector>

#include "para_world.h"
#include "para_collection.h"
#include "para_tag.h"

namespace Parallel
{

/**
 * @brief Divide nproc processes into num_groups groups.
 *
 * Replaces Parallel_Global::divide_mpi_groups.
 * In serial mode, returns group=0, rank_in_group=0, procs_in_group=1.
 *
 * @param[in] nproc        total number of processes
 * @param[in] num_groups   desired number of groups
 * @param[in] rank         global rank
 * @param[in] even         if true, require even split (assert on failure)
 * @param[out] procs_in_group  processes per group
 * @param[out] my_group    which group this rank belongs to
 * @param[out] rank_in_group   rank within the group
 */
void divide_mpi_groups(int nproc, int num_groups, int rank, bool even,
                       int& procs_in_group, int& my_group, int& rank_in_group);

#ifdef __MPI

/**
 * @brief Split MPI_COMM_WORLD into k-pools and band groups.
 *
 * Replaces Parallel_Global::divide_pools.
 * Returns ParaWorld objects for pw (POOL_WORLD), kmesh (KP_WORLD),
 * and bgroup (INT_BGROUP + BP_WORLD) domains.
 *
 * @param[in] nproc    total number of MPI processes
 * @param[in] my_rank  global rank
 * @param[in] bndpar   number of band groups
 * @param[in] kpar     number of k-pools
 * @param[out] nproc_in_pool    processes per pool
 * @param[out] rank_in_pool     rank within pool
 * @param[out] my_pool          pool index
 * @param[out] nproc_in_bndgroup  processes per band group
 * @param[out] rank_in_bpgroup   rank within band group
 * @param[out] my_bndgroup       band group index
 * @param[out] pw_world     ParaWorld for POOL_WORLD
 * @param[out] kmesh_world  ParaWorld for KP_WORLD
 * @param[out] bgroup_int   ParaWorld for INT_BGROUP
 * @param[out] bgroup_bp   ParaWorld for BP_WORLD
 */
void split_pools(int nproc, int my_rank, int bndpar, int kpar,
                 int& nproc_in_pool, int& rank_in_pool, int& my_pool,
                 int& nproc_in_bndgroup, int& rank_in_bpgroup, int& my_bndgroup,
                 ParaWorld& pw_world, ParaWorld& kmesh_world,
                 ParaWorld& bgroup_int, ParaWorld& bgroup_bp);

/**
 * @brief Split for diagonalization domain.
 *
 * Replaces Parallel_Global::split_diag_world.
 *
 * @param[in] diag_np   number of diag groups
 * @param[in] nproc     total processes
 * @param[in] my_rank   global rank
 * @param[out] drank    rank in diag domain
 * @param[out] dsize    size of diag domain
 * @param[out] dcolor   color (group index)
 * @return  ParaWorld for DIAG_WORLD
 */
ParaWorld split_diag_world(int diag_np, int nproc, int my_rank,
                           int& drank, int& dsize, int& dcolor);

/**
 * @brief Split for real-space grid domain.
 *
 * Replaces Parallel_Global::split_grid_world.
 *
 * @param[in] diag_np   number of grid groups (same param as diag)
 * @param[in] nproc     total processes
 * @param[in] my_rank   global rank
 * @param[out] grank    rank in grid domain
 * @param[out] gsize    size of grid domain
 * @return  ParaWorld for GRID_WORLD
 */
ParaWorld split_grid_world(int diag_np, int nproc, int my_rank,
                           int& grank, int& gsize);

/**
 * @brief Assemble all 8 parallel domains into a ParaCollection.
 *
 * This is the top-level initialization function that driver.cpp will call
 * instead of the old divide_pools + split_diag_world + split_grid_world
 * sequence. Each domain is constructed and added to the collection.
 *
 * @param[in] nproc    total MPI processes
 * @param[in] my_rank  global rank
 * @param[in] bndpar   number of band groups
 * @param[in] kpar     number of k-pools
 * @param[in] diag_np  number of diag groups
 * @return  ParaCollection with all 8 domains
 */
ParaCollection setup_para_worlds(int nproc, int my_rank, int bndpar, int kpar, int diag_np);

#endif // __MPI

} // namespace Parallel

#endif // PARA_SETUP_H
