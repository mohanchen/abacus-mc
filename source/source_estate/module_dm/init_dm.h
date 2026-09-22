#ifndef INIT_DM_H
#define INIT_DM_H

#include <complex>
#include <functional>
#include <map>
#include <vector>

#include "source_base/vector3.h"
#include "source_cell/unitcell.h" // use unitcell
#include "source_estate/elecstate.h"// use ElecState
#include "source_psi/psi.h" // use electronic wave functions
#include "source_estate/module_charge/charge.h" // use charge
#include "source_lcao/setup_dm.h" // define Setup_DM
#include "source_hamilt/module_hcontainer/hcontainer.h"

namespace module_dm
{

struct Init_DM_Config
{
    std::string esolver_type;
    int td_stype;
    int nspin;
    double nelec;
    // RT-TDDFT (td_stype==2, esolver_type!="tddft") parameters, passed explicitly
    // to avoid a reverse dependency of source_estate on source_lcao/module_rt.
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>* td_phase_hybrid = nullptr;
    ModuleBase::Vector3<double> td_cart_At;
    // dm2rho lives in source_lcao; pass it as a callback so source_estate does
    // not depend on source_lcao.
    std::function<void(std::vector<hamilt::HContainer<double>*>&, int, Charge*, double, double, bool)> dm2rho_func;
};

template <typename TK>
void init_dm(UnitCell& ucell,
        elecstate::ElecState* pelec,
        LCAO_domain::Setup_DM<TK> &dmat,
        psi::Psi<TK>* psi,
        Charge &chr,
        const int iter,
        const int exx_two_level_step,
        const Init_DM_Config& cfg);

}

#endif
