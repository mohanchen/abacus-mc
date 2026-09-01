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
      my_bndgroup_(0),
      rank_in_bgroup_(0),
      nproc_in_bgroup_(1)
#ifdef __MPI
      ,
      pool_comm_(MPI_COMM_SELF),
      kp_world_comm_(MPI_COMM_NULL),
      int_bgroup_comm_(MPI_COMM_SELF),
      bp_world_comm_(MPI_COMM_NULL),
      grid_comm_(MPI_COMM_NULL),
      diag_comm_(MPI_COMM_NULL)
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
                                 )
    : world_nproc_(world_nproc_in),
      my_rank_(my_rank_in),
      kpar_(kpar_in),
      my_pool_(my_pool_in),
      rank_in_pool_(rank_in_pool_in),
      nproc_in_pool_(nproc_in_pool_in),
      bndpar_(bndpar_in),
      my_bndgroup_(my_bndgroup_in),
      rank_in_bgroup_(rank_in_bgroup_in),
      nproc_in_bgroup_(nproc_in_bgroup_in)
#ifdef __MPI
      ,
      pool_comm_(pool_comm_in),
      kp_world_comm_(kp_world_comm_in),
      int_bgroup_comm_(int_bgroup_comm_in),
      bp_world_comm_(bp_world_comm_in),
      grid_comm_(grid_comm_in),
      diag_comm_(diag_comm_in)
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
    assert(my_bndgroup_ >= 0 && my_bndgroup_ < bndpar_);
    assert(nproc_in_bgroup_ >= 1);
    assert(bndpar_ * nproc_in_bgroup_ == world_nproc_);
    assert(rank_in_bgroup_ >= 0 && rank_in_bgroup_ < nproc_in_bgroup_);
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
