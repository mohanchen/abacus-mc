#ifndef MODULE_DEVICE_DSP_SELECTOR_H_
#define MODULE_DEVICE_DSP_SELECTOR_H_

#ifdef __DSP

#include <memory>

namespace base_device {
namespace memory {

// DSP selector interface
class DspSelector {
public:
    virtual ~DspSelector() = default;
    // Get DSP rank
    virtual int get_rank() const = 0;
    // Set DSP rank
    virtual void set_rank(const int rank) = 0;
};

// Global selector instance
extern std::unique_ptr<DspSelector> dsp_selector;

// Get current DSP selector
DspSelector* get_dsp_selector();

// Create default DSP selector and set it as global
void create_default_selector(const int rank);

} // namespace memory
} // namespace base_device

#endif // end __DSP

#endif // MODULE_DEVICE_DSP_SELECTOR_H_
