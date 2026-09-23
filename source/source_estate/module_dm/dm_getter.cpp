#include "density_matrix.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/libm/libm.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/klist.h"

#include <cstddef>
#include <memory>
#include <stdexcept>

namespace module_dm
{

// get _DMR pointer
template <typename TK, typename TR>
hamilt::HContainer<TR>* DensityMatrix<TK, TR>::get_DMR_pointer(const int ispin) const
{
    if (ispin <= 0 || ispin > this->spin_mult)
    {
        throw std::out_of_range("DensityMatrix::get_DMR_pointer: DMR spin index is out of range");
    }
    if (this->_DMR.size() != static_cast<std::size_t>(this->spin_mult))
    {
        throw std::logic_error("DensityMatrix::get_DMR_pointer: DMR has not been initialized");
    }
    return this->_DMR[ispin - 1];
}

// get _DMK[ik] pointer
template <typename TK, typename TR>
TK* DensityMatrix<TK, TR>::get_DMK_pointer(const int ik) const
{
#ifdef __DEBUG
    assert(ik < this->_nk * this->spin_mult);
#endif
    return const_cast<TK*>(this->_DMK[ik].data());
}

// get a matrix element of density matrix dm(k)
template <typename TK, typename TR>
TK DensityMatrix<TK, TR>::get_DMK(const int ispin, const int ik, const int i, const int j) const
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->spin_mult);
#endif
    // consider transpose col=>row
    return this->_DMK[ik + this->_nk * (ispin - 1)][i * this->pv->nrow + j];
}

// get _DMK nks, nrow, ncol
template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_nks() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() == _nk * spin_mult);
#endif
    return _nk * spin_mult;
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_size() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_DMK.size();
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_nrow() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->pv->nrow;
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_ncol() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->pv->ncol;
}

// T of HContainer can be double or std::complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace module_dm
