//==========================================================
// AUTHOR : Fang Wei, Mohan Chen
// DATE : 2008
// LAST UPDATE : 2009-3-23 mohan add GATHER_MINIMUM_DOUBLE
//==========================================================
#ifndef PARALLEL_GLOBAL_H
#define PARALLEL_GLOBAL_H

#include "parallel_comm.h"
#include "parallel_topology.h"

namespace Parallel_Global
{
extern int mpi_number;
extern int omp_number;
//---------------------------
// call at the very first.
//---------------------------

// changed from read_mpi_parameters in 2024-1018
void read_pal_param(int argc, char** argv, int& NPROC, int& NTHREAD_PER_PROC, int& MY_RANK);

/**
 * @brief Build a ProcessTopology snapshot for the given parallel parameters.
 *
 * This is the single factory that knows how to:
 *   * split MPI_COMM_WORLD into KPAR k-pools via divide_group_comm(even=false);
 *   * split each pool into BNDPAR band groups via divide_group_comm(even=true);
 *   * derive INT_BGROUP (bsame_kdiff_world) / BP_WORLD (bdiff_ksame_world);
 *   * split diag_np-based DIAG_WORLD and diag_np-grouped GRID_WORLD.
 *
 * The returned ProcessTopology::pw_world_comm (previously POOL_WORLD) is the
 * smallest PW tile: the intersection of one k-pool and one band-group.
 *
 * matrix_world_comm and atom_world_comm are left as MPI_COMM_NULL in the
 * returned value; callers that know which distributed view is required for
 * a given step (Parallel_2D / Parallel_Orbitals / DomainDecomposition) are
 * expected to fill them in from the appropriate view before passing the
 * topology down.
 *
 * Note: The factory is only available under __MPI. The non-MPI build path
 * uses the default ProcessTopology constructor which already produces the
 * single-process trivial topology.
 *
 * @param[in] world_nproc Size of MPI_COMM_WORLD
 * @param[in] my_rank    Rank in MPI_COMM_WORLD
 * @param[in] kpar       KPAR   from INPUT (k-point parallelism)
 * @param[in] bndpar     BNDPAR from INPUT (band  parallelism)
 * @param[in] diag_np    Number of diag worlds (also serves as group count
 *                       for the real-space grid world: GRID_WORLD groups are
 *                       the contiguous blocks produced by split_grid_world).
 * @param[in] grid_np    Reserved; currently the real-space grid world is
 *                       tied to diag_np via split_grid_world(diag_np, ...).
 */
ProcessTopology create_topology(int world_nproc,
                                int my_rank,
                                int kpar,
                                int bndpar,
                                int diag_np,
                                int grid_np);


/**-------------------------------------------
 * call to split the "diago world"
 * the unit of first proc of each grid group
 * us the diag world.
 * for example, if we have 64 processors,
 * and diago_proc = 4,
 * then we have 4 'grid world', which
 * have 16, 16, 16, 16 processors each,
 * and the head of each 'grid world'
 * leads to the 'diag world', diag
 * is only carried out using those 4 proc.
 */
void split_diag_world(const int& diag_np, const int& nproc, const int& my_rank, int& drank, int& dsize, int& dcolor);
void split_grid_world(const int diag_np, const int& nproc, const int& my_rank, int& grank, int& gsize);

/**
 * @brief An interface function to call "Parallel_Global::divide_pools()"
 *
 */
void init_pools(const int& NPROC,
                const int& MY_RANK,
                const int& BNDPAR,
                const int& KPAR,
                int& NPROC_IN_BNDGROUP,
                int& RANK_IN_BPGROUP,
                int& MY_BNDGROUP,
                int& NPROC_IN_POOL,
                int& RANK_IN_POOL,
                int& MY_POOL);

void divide_pools(const int& NPROC,
                  const int& MY_RANK,
                  const int& BNDPAR,
                  const int& KPAR,
                  int& NPROC_IN_BNDGROUP,
                  int& RANK_IN_BPGROUP,
                  int& MY_BNDGROUP,
                  int& NPROC_IN_POOL,
                  int& RANK_IN_POOL,
                  int& MY_POOL);

/**
 * @brief Divide MPI processes into groups
 * @param[in] procs Number of MPI processes
 * @param[in] num_groups Number of groups
 * @param[in] rank Rank of the process
 * @param[out] procs_in_group Number of processes in each group
 * @param[out] my_group Group number of the process
 * @param[out] rank_in_group Rank of the process in the group
 * @param[in] even If true, require the number of processes in each group is the same
 */
void divide_mpi_groups(const int& procs,
                       const int& num_groups,
                       const int& rank,
                       int& procs_in_group,
                       int& my_group,
                       int& rank_in_group,
                       const bool even = false);

/**
 * @brief Release MPI communicator and resources
 *
 */
#ifdef __MPI
void finalize_mpi();
#endif

} // namespace Parallel_Global

#endif // PARALLEL_GLOBAL_H
