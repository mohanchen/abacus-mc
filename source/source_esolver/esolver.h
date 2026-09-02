#ifndef ESOLVER_H
#define ESOLVER_H

#include "source_base/matrix.h"
#include "source_base/parallel_partition.h"
#include "source_cell/base_cell.h"
#include "source_cell/unitcell.h"

struct Input_para;

namespace ModuleESolver
{
class ESolver
{
  public:
    ESolver()
    {
        classname = "ESolver";
    }

    virtual ~ESolver()
    {
        //****************************************************
        // do not add any codes in this deconstructor funcion
        //****************************************************
    }

    //! initialize the energy solver by using input parameters and cell modules
    virtual void before_all_runners(BaseCell& cell, const Input_para& inp) = 0;

    //! run energy solver
    virtual void runner(BaseCell& cell, const int istep) = 0;

    //! perform post processing calculations
    virtual void after_all_runners(BaseCell& cell) = 0;

    //! deal with exx and other calculation than scf/md/relax/cell-relax:
    //! such as nscf, get_wf and get_pchg
    virtual void others(BaseCell&, const int) {}

    //! calculate total energy of a given system
    virtual double cal_energy() = 0;

    //! calcualte forces for the atoms in the given cell
    virtual void cal_force(BaseCell& cell, ModuleBase::matrix& force) = 0;

    //! calcualte stress of given cell
    virtual void cal_stress(BaseCell& cell, ModuleBase::matrix& stress) = 0;

    virtual bool supports_mdcell() const
    {
        return false;
    }

    virtual double mdcell_cutoff(const Input_para& inp) const
    {
        static_cast<void>(inp);
        return 0.0;
    }

    bool conv_esolver = true; // whether esolver is converged

    std::string classname;

    /// Inject the parallel partition built by the driver (value copy).
    /// Called before before_all_runners; the solver may afterwards derive
    /// solver-specific domains via topo_ = topo_.with_*_world_comm(...).
    void set_topology(const ParallelPartition& topo)
    {
        topo_ = topo;
    }

  protected:
    /// Bound in before_all_runners; members use inp_->xxx instead of PARAM.inp.xxx
    const Input_para* inp_ = nullptr;

    /**
     * @brief Parallel partition for this solver, injected by the driver
     *        before before_all_runners.
     *
     * The world-level domains (pw / kmesh / bsame_kdiff / bdiff_ksame /
     * rgrid / diag) come from Parallel_Global::create_partition. An
     * ESolver variant that needs solver-specific domains (matrix, atom,
     * or a custom re-split of a world domain) derives its own view in
     * before_all_runners via the immutable with_*_world_comm builders,
     * e.g.:
     *     topo_ = topo_.with_matrix_world_comm(blacs_grid_comm);
     * The partition does not own communicators; the owner of a split
     * communicator must free it (e.g. in the ESolver destructor).
     */
    ParallelPartition topo_;
};

} // namespace ModuleESolver

#endif
