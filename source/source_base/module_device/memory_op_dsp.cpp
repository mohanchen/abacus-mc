#include "memory_op_dsp.h"
#include "dsp_selector.h"

#include "source_base/memory.h"
#include "source_base/tool_threading.h"
#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#endif

#include <complex>
#include <cstring>

namespace base_device
{
namespace memory
{

#ifdef __DSP

template <typename FPTYPE>
struct resize_memory_op_mt<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(FPTYPE*& arr, const size_t size, const char* record_in)
    {
        if (arr != nullptr)
        {
            mtfunc::free_ht(arr);
        }
        int rank = get_dsp_selector()->get_rank();
        arr = (FPTYPE*)mtfunc::malloc_ht(sizeof(FPTYPE) * size, rank);
        std::string record_string;
        if (record_in != nullptr)
        {
            record_string = record_in;
        }
        else
        {
            record_string = "no_record";
        }

        if (record_string != "no_record")
        {
            ModuleBase::Memory::record(record_string, sizeof(FPTYPE) * size);
        }
    }
};

template <typename FPTYPE>
struct set_memory_op_mt<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(FPTYPE* arr, const int var, const size_t size)
    {
        ModuleBase::OMP_PARALLEL([&](int num_thread, int thread_id)
        {
            int beg = 0, len = 0;
            ModuleBase::BLOCK_TASK_DIST_1D(num_thread, thread_id, size, (size_t)4096 / sizeof(FPTYPE), beg, len);
            memset(arr + beg, var, sizeof(FPTYPE) * len);
        });
    }
};

template <typename FPTYPE>
struct delete_memory_op_mt<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(FPTYPE* arr)
    {
        mtfunc::free_ht(arr);
    }
};


template struct resize_memory_op_mt<int, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<float, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<double, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<std::complex<float>, base_device::DEVICE_CPU>;
template struct resize_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>;

template struct set_memory_op_mt<int, base_device::DEVICE_CPU>;
template struct set_memory_op_mt<float, base_device::DEVICE_CPU>;
template struct set_memory_op_mt<double, base_device::DEVICE_CPU>;
template struct set_memory_op_mt<std::complex<float>, base_device::DEVICE_CPU>;
template struct set_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>;

template struct delete_memory_op_mt<int, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<float, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<double, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<std::complex<float>, base_device::DEVICE_CPU>;
template struct delete_memory_op_mt<std::complex<double>, base_device::DEVICE_CPU>;
#endif

} // namespace memory
} // namespace base_device