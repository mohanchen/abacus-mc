#include "charge_mixing.h"
#include "chg_dmr.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"

void Charge_Mixing::allocate_mixing_dmr(const int nnr)
{
    ModuleBase::TITLE("Charge_Mixing", "allocate_mixing_dmr");
    ModuleBase::timer::start("Charge_Mixing", "allocate_mixing_dmr");
    module_charge::init_mixing_dmr(this->mixing, this->dmr_mdata, nnr, this->cfg_);
    ModuleBase::timer::end("Charge_Mixing", "allocate_mixing_dmr");

    return;
}

void Charge_Mixing::mix_dmr(elecstate::DensityMatrix<double, double>* DM)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_dmr");
    ModuleBase::timer::start("Charge_Mixing", "mix_dmr");
    module_charge::mix_dmr(DM, this->mixing, this->dmr_mdata, this->cfg_);
    ModuleBase::timer::end("Charge_Mixing", "mix_dmr");

    return;
}

void Charge_Mixing::mix_dmr(elecstate::DensityMatrix<std::complex<double>, double>* DM)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_dmr");
    ModuleBase::timer::start("Charge_Mixing", "mix_dmr");
    module_charge::mix_dmr(DM, this->mixing, this->dmr_mdata, this->cfg_);
    ModuleBase::timer::end("Charge_Mixing", "mix_dmr");

    return;
}
