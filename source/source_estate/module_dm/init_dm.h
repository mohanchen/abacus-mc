#ifndef INIT_DM_H
#define INIT_DM_H

#include "source_cell/unitcell.h" // use unitcell
#include "source_estate/elecstate.h"// use ElecState
#include "source_psi/psi.h" // use electronic wave functions
#include "source_estate/module_charge/charge.h" // use charge
#include "source_lcao/setup_dm.h" // define Setup_DM

namespace elecstate
{

struct Init_DM_Config
{
    std::string esolver_type;
    int td_stype;
    int nspin;
    double nelec;
};

template <typename TK>
void init_dm(UnitCell& ucell,
        ElecState* pelec,
        LCAO_domain::Setup_DM<TK> &dmat,
        psi::Psi<TK>* psi,
        Charge &chr,
        const int iter,
        const int exx_two_level_step,
        const Init_DM_Config& cfg);

}

#endif
