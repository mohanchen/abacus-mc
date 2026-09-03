#include "para_rgrid_world.h"

#include <cassert>

namespace Parallel
{

ParaRgridWorld::ParaRgridWorld(int ncx, int ncy, int ncz)
    : ParaWorld("rgrid"), ncx_(ncx), ncy_(ncy), ncz_(ncz)
{
    assert(ncx > 0 && ncy > 0 && ncz > 0);
    distribute_z();
    nczp_ = numz_[0];
}

#ifdef __MPI
ParaRgridWorld::ParaRgridWorld(const MPI_Comm& comm, int ncx, int ncy, int ncz)
    : ParaWorld("rgrid", comm), ncx_(ncx), ncy_(ncy), ncz_(ncz)
{
    assert(ncx > 0 && ncy > 0 && ncz > 0);
    distribute_z();
    nczp_ = numz_[rank()];
}
#endif

void ParaRgridWorld::distribute_z()
{
    const int np = size();
    numz_.resize(np, 0);
    startz_.resize(np, 0);
    whichpro_.resize(ncz_, 0);

    // Evenly distribute z-planes, remainder to front processes
    const int base = ncz_ / np;
    const int rem = ncz_ % np;
    int acc = 0;
    for (int p = 0; p < np; ++p)
    {
        numz_[p] = base + (p < rem ? 1 : 0);
        startz_[p] = acc;
        acc += numz_[p];
    }

    // Build owner table
    for (int p = 0; p < np; ++p)
    {
        for (int iz = 0; iz < numz_[p]; ++iz)
        {
            whichpro_[startz_[p] + iz] = p;
        }
    }
}

int ParaRgridWorld::numz(int p) const
{
    assert(p >= 0 && p < static_cast<int>(numz_.size()));
    return numz_[p];
}

int ParaRgridWorld::startz(int p) const
{
    assert(p >= 0 && p < static_cast<int>(startz_.size()));
    return startz_[p];
}

int ParaRgridWorld::whichpro(int iz) const
{
    assert(iz >= 0 && iz < ncz_);
    return whichpro_[iz];
}

} // namespace Parallel
