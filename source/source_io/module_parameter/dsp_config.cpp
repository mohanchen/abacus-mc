#include "source_base/module_device/dsp_selector.h"
#include <memory>
#include <stdexcept>

#ifdef __DSP

namespace ModuleIO {

// Initialize DSP selector
void init_dsp_selector(const int my_rank, const int dsp_count)
{
    // Validate parameters
    if (my_rank < 0)
    {
        throw std::runtime_error(
            "ModuleIO::init_dsp_selector: "
            "my_rank must be non-negative"
        );
    }
    if (dsp_count <= 0)
    {
        throw std::runtime_error(
            "ModuleIO::init_dsp_selector: "
            "dsp_count must be positive"
        );
    }
    
    // Calculate DSP rank
    const int rank = my_rank % dsp_count;
    
    // Create default DSP selector and set it as global
    base_device::memory::create_default_selector(rank);
}

} // namespace ModuleIO

#endif
