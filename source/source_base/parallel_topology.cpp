#include "parallel_topology.h"

#include <cassert>

ProcessTopology::ProcessTopology()
    : world_nproc_(1),
      my_rank_(0),
      kpar_(1),
      my_pool_(0),
      rank_in_pool_(0),
      nproc_in_pool_(1, 1),
      bndpar_(1),
      my_band_group_(0),
      rank_in_band_group_(0),
      nproc_in_band_group_(1)
#ifdef __MPI
      ,
      pw_world_comm_(MPI_COMM_SELF),
      kmesh_world_comm_(MPI_COMM_NULL),
      bsame_kdiff_world_comm_(MPI_COMM_SELF),
      bdiff_ksame_world_comm_(MPI_COMM_NULL),
      rgrid_world_comm_(MPI_COMM_NULL),
      diag_world_comm_(MPI_COMM_NULL),
      matrix_world_comm_(MPI_COMM_NULL),
      atom_world_comm_(MPI_COMM_NULL)
#endif
{
}

ProcessTopology::ProcessTopology(int world_nproc_in,
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
                                 MPI_Comm matrix_world_comm_in,
                                 MPI_Comm atom_world_comm_in
#endif
                                 )
    : world_nproc_(world_nproc_in),
      my_rank_(my_rank_in),
      kpar_(kpar_in),
      my_pool_(my_pool_in),
      rank_in_pool_(rank_in_pool_in),
      nproc_in_pool_(nproc_in_pool_in),
      bndpar_(bndpar_in),
      my_band_group_(my_band_group_in),
      rank_in_band_group_(rank_in_band_group_in),
      nproc_in_band_group_(nproc_in_band_group_in)
#ifdef __MPI
      ,
      pw_world_comm_(pw_world_comm_in),
      kmesh_world_comm_(kmesh_world_comm_in),
      bsame_kdiff_world_comm_(bsame_kdiff_world_comm_in),
      bdiff_ksame_world_comm_(bdiff_ksame_world_comm_in),
      rgrid_world_comm_(rgrid_world_comm_in),
      diag_world_comm_(diag_world_comm_in),
      matrix_world_comm_(matrix_world_comm_in),
      atom_world_comm_(atom_world_comm_in)
#endif
{
    assert(world_nproc_ >= 1);
    assert(my_rank_ >= 0 && my_rank_ < world_nproc_);
    assert(kpar_ >= 1);
    assert(static_cast<int>(nproc_in_pool_.size()) == kpar_);
    int total = 0;
    for (int s : nproc_in_pool_)
    {
        assert(s >= 0);
        total += s;
    }
    assert(total == world_nproc_);
    assert(my_pool_ >= 0 && my_pool_ < kpar_);
    assert(rank_in_pool_ >= 0 && rank_in_pool_ < nproc_in_pool_[my_pool_]);

    assert(bndpar_ >= 1);
    assert(my_band_group_ >= 0 && my_band_group_ < bndpar_);
    assert(nproc_in_band_group_ >= 1);
    assert(bndpar_ * nproc_in_band_group_ == world_nproc_);
    assert(rank_in_band_group_ >= 0 && rank_in_band_group_ < nproc_in_band_group_);
}

int ProcessTopology::pool_root_rank(int pool) const
{
    if (pool < 0 || pool >= kpar_)
    {
        return -1;
    }
    int offset = 0;
    for (int i = 0; i < pool; ++i)
    {
        offset += nproc_in_pool_[i];
    }
    return offset;
}

int ProcessTopology::band_group_root_rank(int band_group) const
{
    if (band_group < 0 || band_group >= bndpar_)
    {
        return -1;
    }
    // In ABACUS divide_pools the band-group layout over world ranks is
    // stripe-contiguous inside each k pool: within pool P the first
    // (nproc_in_pool[P]/bndpar) ranks belong to band-group 0, the next
    // slice to band-group 1, and so on.  The first rank of the
    // concatenated "same band-group across all pools" set (i.e. the
    // root of bsame_kdiff_world for that band-group) is therefore the
    // first occurrence in pool 0, which falls at offset band_group *
    // (nproc_in_pool[0]/bndpar) from pool_root_rank(0).  The invariant
    // bndpar_ * nproc_in_band_group_ == world_nproc_ + the even split
    // enforced by MPICommGroup::divide_group_comm make that offset
    // equal to (band_group * nproc_in_band_group_) directly because
    // each band-group contains exactly nproc_in_band_group_ processes
    // globally and they appear in ascending band-group id order in
    // world rank when scanned pool by pool.
    const int per_bg_in_pool0 = nproc_in_pool_[0] / bndpar_;
    assert(per_bg_in_pool0 * bndpar_ == nproc_in_pool_[0]);
    const int via_pool0 = band_group * per_bg_in_pool0;
    const int via_global = band_group * nproc_in_band_group_;
    assert(via_pool0 == via_global);
    return via_global;
}
