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
    // slice to band-group 1, and so on.  Pool 0 always starts at world
    // rank 0 (divide_mpi_groups assigns pools as contiguous blocks in
    // ascending world-rank order), so the first member of the
    // "same band-group across all pools" set -- i.e. the root of
    // bsame_kdiff_world for that band-group -- is the first occurrence
    // in pool 0, which falls at offset band_group * (nproc_in_pool[0]
    // / bndpar) from world rank 0.
    //
    // NOTE: this is generally NOT equal to band_group *
    // nproc_in_band_group_.  The global band-group size is
    // kpar * (nproc_in_pool[0] / bndpar), so the two formulas coincide
    // only when kpar == 1 (or band_group == 0).
    const int per_bg_in_pool0 = nproc_in_pool_[0] / bndpar_;
    assert(per_bg_in_pool0 * bndpar_ == nproc_in_pool_[0]);
    return band_group * per_bg_in_pool0;
}

#ifdef __MPI
// Immutable builders: each returns a copy of *this with exactly one
// communicator replaced. Scalar topology fields are never touched.
ProcessTopology ProcessTopology::with_pw_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.pw_world_comm_ = comm;
    return t;
}

ProcessTopology ProcessTopology::with_kmesh_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.kmesh_world_comm_ = comm;
    return t;
}

ProcessTopology ProcessTopology::with_bsame_kdiff_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.bsame_kdiff_world_comm_ = comm;
    return t;
}

ProcessTopology ProcessTopology::with_bdiff_ksame_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.bdiff_ksame_world_comm_ = comm;
    return t;
}

ProcessTopology ProcessTopology::with_rgrid_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.rgrid_world_comm_ = comm;
    return t;
}

ProcessTopology ProcessTopology::with_diag_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.diag_world_comm_ = comm;
    return t;
}

ProcessTopology ProcessTopology::with_matrix_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.matrix_world_comm_ = comm;
    return t;
}

ProcessTopology ProcessTopology::with_atom_world_comm(MPI_Comm comm) const
{
    ProcessTopology t = *this;
    t.atom_world_comm_ = comm;
    return t;
}
#endif
