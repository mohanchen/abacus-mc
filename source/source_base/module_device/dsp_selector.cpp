#include "dsp_selector.h"
#include <string>
#include <stdexcept>

#ifdef __DSP

namespace base_device
{
namespace memory
{

// Global selector instance
std::unique_ptr<DspSelector> dsp_selector = nullptr;

// Get current DSP selector
DspSelector* get_dsp_selector()
{
    if (!dsp_selector)
    {
        throw std::runtime_error(
            "ModuleBase::memory::get_dsp_selector: "
            "DSP selector not initialized. Call init_dsp_selector first."
        );
    }
    return dsp_selector.get();
}

// Default DSP selector implementation
class DefaultDspSelector : public DspSelector
{
private:
    int rank_ = 0;

public:
    int get_rank() const override
    {
        return rank_;
    }

    void set_rank(const int rank) override
    {
        if (rank < 0)
        {
            throw std::runtime_error(
                "ModuleBase::memory::DspSelector: "
                "DSP rank must be non-negative"
            );
        }
        rank_ = rank;
    }
};


// Create default DSP selector and set it as global
void create_default_selector(const int rank)
{
    dsp_selector = std::unique_ptr<DefaultDspSelector>(new DefaultDspSelector());
    dsp_selector->set_rank(rank);
}

} // namespace memory
} // namespace base_device

#endif 
