#include "chg_parallel.h"

#ifdef __MPI

#include <cassert>

#include "charge.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

namespace module_charge
{

void reduce_diff_pools(double* array_rho, const Charge& chr)
{
    ModuleBase::TITLE("Charge", "reduce_diff_pools");
    ModuleBase::timer::start("Charge", "reduce_diff_pools");
    assert(array_rho != nullptr);
    if (GlobalV::KPAR > 1)
    {
        assert(chr.pgrid != nullptr);
        chr.pgrid->reduce_across_pools(array_rho);
    }
    if (PARAM.globalv.all_ks_run && PARAM.inp.bndpar > 1)
    {
        assert(chr.nrxx > 0);
        MPI_Allreduce(MPI_IN_PLACE, array_rho, chr.nrxx, MPI_DOUBLE, MPI_SUM, BP_WORLD);
    }
    ModuleBase::timer::end("Charge", "reduce_diff_pools");
}

void rho_mpi(Charge& chr)
{
    ModuleBase::TITLE("Charge", "rho_mpi");
    if (GlobalV::KPAR * PARAM.inp.bndpar <= 1)
    {
        return;
    }
    ModuleBase::timer::start("Charge", "rho_mpi");

    assert(chr.rho != nullptr);
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        reduce_diff_pools(chr.rho[is], chr);
        if (XC_Functional::get_ked_flag() || PARAM.inp.out_elf[0] > 0)
        {
            assert(chr.kin_r != nullptr);
            reduce_diff_pools(chr.kin_r[is], chr);
        }
    }

    ModuleBase::timer::end("Charge", "rho_mpi");
    return;
}

void kin_r_mpi(Charge& chr)
{
    ModuleBase::TITLE("Charge", "kin_r_mpi");
    if (GlobalV::KPAR * PARAM.inp.bndpar <= 1)
    {
        return;
    }
    ModuleBase::timer::start("Charge", "kin_r_mpi");

    if (XC_Functional::get_ked_flag() || PARAM.inp.out_elf[0] > 0)
    {
        assert(chr.kin_r != nullptr);
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            reduce_diff_pools(chr.kin_r[is], chr);
        }
    }

    ModuleBase::timer::end("Charge", "kin_r_mpi");
    return;
}

} // namespace module_charge

#endif
