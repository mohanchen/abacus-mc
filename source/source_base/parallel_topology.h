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
 * associated pool / band-group rank & size info from GlobalV, and
 * also exposes explicit handles for two long-implicit domains:
 *
 *   - matrix_world_comm: communicator backing the 2D block-cyclic
 *     BLACS grid used by Parallel_2D / Parallel_Orbitals to lay out
 *     distributed H/S/projector matrices. Historically callers chose
 *     pw_world / diag_world / MPI_COMM_WORLD ad-hoc; here we give the
 *     concept a first-class name so call sites inject the appropriate
 *     view explicitly in step 2.
 *   - atom_world_comm : communicator backing the 3D Cartesian domain
 *     decomposition used by the MD / neighlist
 *     DomainDecomposition to distribute atoms spatially. Today the
 *     concrete communicator is always MPI_COMM_WORLD (wrapped via the
 *     CommunicationDomain thin shell in parallel_cell.h).
 *
 * All 8 communicator handles share the consistent `_world_comm`
 * suffix. The pair (bdiff_ksame_world, bsame_kdiff_world) is named
 * directly after the relation between band side and k side:
 *
 *   bdiff_ksame  ->  band groups differ, k is the same
 *                    (legacy BP_WORLD, intra-pool cross-band pairs)
 *   bsame_kdiff  ->  band groups are identical, k differs
 *                    (legacy INT_BGROUP, intra-band-group union across
 *                    all k pools)
 *
 * Values are captured at construction time; the object is copyable
 * (value semantics) and exposes no mutable workflow switches.
 * Consumers currently reading the global MPI_Comm / GlobalV fields
 * should migrate to taking a `const ProcessTopology&` from their
 * callers. Parallel_Global::create_topology(...) produces one such
 * instance at startup whose communicators are also aliased back to
 * the legacy globals so existing code keeps working during the
 * migration.
 *
 * matrix_world_comm / atom_world_comm are not derived inside the
 * topology constructor because the correct choice depends on where
 * the view is used (LCAO diag, GK diag, MD step, ...). Callers that
 * need a proper matrix or atom domain are expected to fill in the
 * handle before passing the topology down to distributed modules.
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
     * matrix_world_comm / atom_world_comm default to MPI_COMM_NULL
     * and are filled in later by whichever call site knows which
     * domain view is needed for a given distributed calculation.
     */
    ProcessTopology(int world_nproc_in,
                    int my_rank_in,
                    int kpar_in,
                    int my_pool_in,
                    int rank_in_pool_in,
                    const std::vector<int>& nproc_in_pool_in,
                    int bndpar_in,
                    int my_band_group_in,
                    int rank_in_band_group_in,
                    int nproc_in_band_group_in
#ifdef __MPI
                        ,
                    MPI_Comm pw_world_comm_in,
                    MPI_Comm kmesh_world_comm_in,
                    MPI_Comm bsame_kdiff_world_comm_in,
                    MPI_Comm bdiff_ksame_world_comm_in,
                    MPI_Comm rgrid_world_comm_in,
                    MPI_Comm diag_world_comm_in,
                    MPI_Comm matrix_world_comm_in = MPI_COMM_NULL,
                    MPI_Comm atom_world_comm_in   = MPI_COMM_NULL
#endif
    );

    // world ----------------------------------------------------------
    int world_size() const { return world_nproc_; }
    int world_rank() const { return my_rank_; }

    // k-point pools --------------------------------------------------
    int kpar() const { return kpar_; }
    int my_pool() const { return my_pool_; }
    int rank_in_pool() const { return rank_in_pool_; }
    const std::vector<int>& nproc_in_pool() const { return nproc_in_pool_; }
    int nproc_in_pool(int pool) const { return nproc_in_pool_[pool]; }
    /// World rank of the root (rank 0) of a given pool; -1 on bad index.
    int pool_root_rank(int pool) const;

    // band groups (BNDPAR) ------------------------------------------
    int bndpar() const { return bndpar_; }
    int my_band_group() const { return my_band_group_; }
    int rank_in_band_group() const { return rank_in_band_group_; }
    int nproc_in_band_group() const { return nproc_in_band_group_; }
    /// World rank of the first member (root) of a band-group union.
    /// Band-group slices are stripe-contiguous inside each k pool, and
    /// pool 0 starts at world rank 0, so the root is
    /// band_group * (nproc_in_pool[0] / bndpar). This equals
    /// band_group * nproc_in_band_group() only when kpar == 1.
    /// -1 on invalid index.
    int band_group_root_rank(int band_group) const;

#ifdef __MPI
    // communicators (consistent `_world_comm` suffix) ---------------
    MPI_Comm pw_world_comm()       const { return pw_world_comm_; }
    MPI_Comm kmesh_world_comm()    const { return kmesh_world_comm_; }
    MPI_Comm bsame_kdiff_world_comm() const { return bsame_kdiff_world_comm_; }
    MPI_Comm bdiff_ksame_world_comm() const { return bdiff_ksame_world_comm_; }
    MPI_Comm rgrid_world_comm()    const { return rgrid_world_comm_; }
    MPI_Comm diag_world_comm()     const { return diag_world_comm_; }
    MPI_Comm matrix_world_comm()   const { return matrix_world_comm_; }
    MPI_Comm atom_world_comm()     const { return atom_world_comm_; }
#endif

  private:
    int world_nproc_ = 1;
    int my_rank_ = 0;

    int kpar_ = 1;
    int my_pool_ = 0;
    int rank_in_pool_ = 0;
    std::vector<int> nproc_in_pool_;

    int bndpar_ = 1;
    int my_band_group_ = 0;
    int rank_in_band_group_ = 0;
    int nproc_in_band_group_ = 1;

#ifdef __MPI
    MPI_Comm pw_world_comm_       = MPI_COMM_SELF;
    MPI_Comm kmesh_world_comm_    = MPI_COMM_NULL;
    MPI_Comm bsame_kdiff_world_comm_ = MPI_COMM_SELF;
    MPI_Comm bdiff_ksame_world_comm_ = MPI_COMM_NULL;
    MPI_Comm rgrid_world_comm_    = MPI_COMM_NULL;
    MPI_Comm diag_world_comm_     = MPI_COMM_NULL;
    MPI_Comm matrix_world_comm_   = MPI_COMM_NULL;
    MPI_Comm atom_world_comm_     = MPI_COMM_NULL;
#endif
};

#endif // PARALLEL_TOPOLOGY_H
