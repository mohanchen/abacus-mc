#ifndef PARALLEL_TOPOLOGY_H
#define PARALLEL_TOPOLOGY_H

#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

/**
 * @brief Immutable description of the ABACUS process topology.
 *
 * Replaces the six raw-global MPI_Comm (POOL_WORLD / KP_WORLD /
 * INT_BGROUP / BP_WORLD / GRID_WORLD / DIAG_WORLD) together with the
 * associated pool / band-group rank & size info from GlobalV.
 * Values captured at construction time; the object is copyable
 * (value semantics) and exposes no mutable workflow switches.
 *
 * Consumers currently reading the global MPI_Comm / GlobalV fields
 * should migrate to taking a `const ProcessTopology&` from their
 * callers. Parallel_Global::create_topology(...) produces one such
 * instance at startup whose communicators are also aliased back to
 * the legacy globals so existing code keeps working during the
 * migration.
 */
class ProcessTopology
{
  public:
    /// Trivially-constructible serial fallback. Produces a single-
    /// process, single-pool topology suitable for non-__MPI builds.
    ProcessTopology();

    /**
     * @brief Full constructor. Intended for Parallel_Global factory.
     *
     * All inputs are plain values; GlobalV is never read inside the
     * constructor. Pass MPI_COMM_NULL for communicators that are
     * unused on this rank (same semantics as legacy divide_pools).
     */
    ProcessTopology(int world_nproc_in,
                    int my_rank_in,
                    int kpar_in,
                    int my_pool_in,
                    int rank_in_pool_in,
                    const std::vector<int>& nproc_in_pool_in,
                    int bndpar_in,
                    int my_bndgroup_in,
                    int rank_in_bgroup_in,
                    int nproc_in_bgroup_in
#ifdef __MPI
                        ,
                    MPI_Comm pool_comm_in,
                    MPI_Comm kp_world_comm_in,
                    MPI_Comm int_bgroup_comm_in,
                    MPI_Comm bp_world_comm_in,
                    MPI_Comm grid_comm_in,
                    MPI_Comm diag_comm_in
#endif
    );

    // world ----------------------------------------------------------
    int world_size() const { return world_nproc_; }
    int world_rank() const { return my_rank_; }

    // k-point pools -------------------------------------------------
    int kpar() const { return kpar_; }
    int my_pool() const { return my_pool_; }
    int rank_in_pool() const { return rank_in_pool_; }
    const std::vector<int>& nproc_in_pool() const { return nproc_in_pool_; }
    int nproc_in_pool(int pool) const { return nproc_in_pool_[pool]; }
    /// World rank of the root (rank 0) of a given pool; -1 on bad index.
    int pool_root_rank(int pool) const;

    // band groups (bndpar) ------------------------------------------
    int bndpar() const { return bndpar_; }
    int my_bndgroup() const { return my_bndgroup_; }
    int rank_in_bgroup() const { return rank_in_bgroup_; }
    int nproc_in_bgroup() const { return nproc_in_bgroup_; }

#ifdef __MPI
    // communicators -------------------------------------------------
    MPI_Comm pool_comm() const { return pool_comm_; }
    MPI_Comm kp_world_comm() const { return kp_world_comm_; }
    MPI_Comm int_bgroup_comm() const { return int_bgroup_comm_; }
    MPI_Comm bp_world_comm() const { return bp_world_comm_; }
    MPI_Comm grid_comm() const { return grid_comm_; }
    MPI_Comm diag_comm() const { return diag_comm_; }
#endif

  private:
    int world_nproc_ = 1;
    int my_rank_ = 0;

    int kpar_ = 1;
    int my_pool_ = 0;
    int rank_in_pool_ = 0;
    std::vector<int> nproc_in_pool_;

    int bndpar_ = 1;
    int my_bndgroup_ = 0;
    int rank_in_bgroup_ = 0;
    int nproc_in_bgroup_ = 1;

#ifdef __MPI
    MPI_Comm pool_comm_ = MPI_COMM_SELF;
    MPI_Comm kp_world_comm_ = MPI_COMM_NULL;
    MPI_Comm int_bgroup_comm_ = MPI_COMM_SELF;
    MPI_Comm bp_world_comm_ = MPI_COMM_NULL;
    MPI_Comm grid_comm_ = MPI_COMM_NULL;
    MPI_Comm diag_comm_ = MPI_COMM_NULL;
#endif
};

#endif // PARALLEL_TOPOLOGY_H
