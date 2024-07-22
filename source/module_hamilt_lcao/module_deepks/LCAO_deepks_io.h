#ifndef LCAO_DEEPKS_IO 
#define LCAO_DEEPKS_IO 

#ifdef __DEEPKS

/*
#include "module_base/complexmatrix.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_io/winput.h"

#include <torch/script.h>
#include <torch/torch.h>
#include <unordered_map>
*/

namespace LCAO_deepks_io
{

    /// print density matrices
    void print_dm(const std::vector<double>& dm);

    void print_dm_k(const int nks, const std::vector<std::vector<std::complex<double>>>& dm);

    void load_npy_gedm(const int nat);

    ///----------------------------------------------------------------------
    /// The following 4 functions save the `[dm_eig], [e_base], [f_base], [grad_vx]`
    /// of current configuration as `.npy` file, when `deepks_scf = 1`.
    /// After a full group of consfigurations are calculated,
    /// we need a python script to `load` and `torch.cat` these `.npy` files,
    /// and get `l_e_delta,npy` and `l_f_delta.npy` corresponding to the exact E, F data.
    ///
    /// Unit of energy: Ry
    /// Unit of force: Ry/Bohr
    ///----------------------------------------------------------------------
    void save_npy_d(const int nat);
    void save_npy_gvx(const int nat);
    void save_npy_gvepsl(const int nat);

    void save_npy_e(const double& e /**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/, const std::string& e_file);
    void save_npy_f(const ModuleBase::matrix& fbase /**<[in] \f$F_{base}\f$ or \f$F_{tot}\f$, in Ry/Bohr*/,
                    const std::string& f_file,
                    const int nat);

    void save_npy_s(const ModuleBase::matrix& sbase /**<[in] \f$S_{base}\f$ or \f$S_{tot}\f$, in Ry/Bohr^3*/,
                    const std::string& s_file,
                    const double omega);

    // QO added on 2021-12-15
    void save_npy_o(const ModuleBase::matrix& bandgap /**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/,
                    const std::string& o_file,
                    const int nks);
    void save_npy_orbital_precalc(const int nat, const int nks);


    //xinyuan added on 2023-2-20
    void save_npy_h(const ModuleBase::matrix &H,const std::string &h_file,const int nlocal);//just for gamma only
    void save_npy_v_delta_precalc(const int nat, const int nks,const int nlocal);
    void save_npy_psialpha(const int nat, const int nks,const int nlocal);
    void save_npy_gevdm(const int nat);
};

#endif
#endif
