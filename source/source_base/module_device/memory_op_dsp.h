#ifndef MODULE_DEVICE_MEMORY_DSP_H_
#define MODULE_DEVICE_MEMORY_DSP_H_

#include "types.h"

#include <complex>
#include <cstddef>

namespace base_device
{

namespace memory
{

#ifdef __DSP

template <typename FPTYPE, typename Device>
struct resize_memory_op_mt
{
    /// @brief Allocate memory for a given pointer. Note this op will free the pointer first.
    ///
    /// Input Parameters
    /// \param size : array size
    /// \param record_string : label for memory record
    ///
    /// Output Parameters
    /// \param arr : allocated array
    void operator()(FPTYPE*& arr, const size_t size, const char* record_in = nullptr);
};

template <typename FPTYPE, typename Device>
struct set_memory_op_mt
{
    /// @brief memset for DSP memory allocated by mt allocator.
    ///
    /// Input Parameters
    /// \param var : the specified constant byte value
    /// \param size : array size
    ///
    /// Output Parameters
    /// \param arr : output array initialized by the input value
    void operator()(FPTYPE* arr, const int var, const size_t size);
};

template <typename FPTYPE, typename Device>
struct delete_memory_op_mt
{
    /// @brief free memory for multi-device
    ///
    /// Input Parameters
    /// \param arr : the input array
    void operator()(FPTYPE* arr);
};

#endif // __DSP

} // end of namespace memory
} // end of namespace base_device

#endif // MODULE_DEVICE_MEMORY_DSP_H_