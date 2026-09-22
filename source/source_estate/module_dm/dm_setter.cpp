#include "density_matrix.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/libm/libm.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/klist.h"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>

namespace module_dm
{

// set DMK using a pointer
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK_pointer(const int ik, TK* DMK_in)
{
#ifdef __DEBUG
    assert(ik < this->_nk * this->_nspin);
#endif
    this->_DMK[ik].assign(DMK_in, DMK_in + this->pv->nrow * this->pv->ncol);
}

// set _DMK element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK(const int ispin, const int ik, const int i, const int j, const TK value)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
    assert(ik >= 0 && ik < this->_nk);
#endif
    // consider transpose col=>row
    this->_DMK[ik + this->_nk * (ispin - 1)][i * this->pv->nrow + j] = value;
}

// set _DMK element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK_zero()
{
    for (int ik = 0; ik < _nspin * _nk; ik++)
    {
        std::fill(this->_DMK[ik].begin(), this->_DMK[ik].end(), TK{});
    }
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::save_DMR()
{
    ModuleBase::TITLE("DensityMatrix", "save_DMR");
    ModuleBase::timer::start("DensityMatrix", "save_DMR");

    const int nnr = this->_DMR[0]->get_nnr();
    // allocate if _DMR_save is empty
    if (_DMR_save.size() == 0)
    {
        _DMR_save.resize(this->_DMR.size());
    }
    // resize if _DMR_save[is].size is not equal to _DMR.size
    for (int is = 0; is < _DMR_save.size(); is++)
    {
        if (_DMR_save[is].size() != nnr)
        {
            _DMR_save[is].resize(nnr);
        }
    }
    // save _DMR to _DMR_save
    for (int is = 0; is < this->_DMR.size(); is++)
    {
        TR* DMR_pointer = this->_DMR[is]->get_wrapper();
        TR* DMR_save_pointer = _DMR_save[is].data();
        // The resize above value-initializes newly added elements, and the
        // whole [0, nnr) range is overwritten by the copy, so a prior
        // zeroing of the destination would be a dead store.
        std::copy(DMR_pointer, DMR_pointer + nnr, DMR_save_pointer);
    }

    ModuleBase::timer::end("DensityMatrix", "save_DMR");
}

// T of HContainer can be double or std::complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace module_dm
