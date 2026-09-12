#ifndef DFTU_LCAO_H
#define DFTU_LCAO_H

#include "source_pw/module_pwdft/dftu_base.h"

#include <complex>
#include <string>
#include <vector>


class UnitCell;
class Parallel_Orbitals;

#ifdef __LCAO
class LCAO_Orbitals;
#endif


class Plus_U : public Plus_U_Base
{

  public:
    Plus_U();
    ~Plus_U();

    // allocate relevant data strcutures
    void init(UnitCell& cell,
                const Parallel_Orbitals* pv,
                const int npol,
                const int nspin,
                const std::vector<int>& l_channel,
                const bool yukawa_potential,
                const double yukawa_lambda,
                const std::string& global_readin_dir,
                const std::string& global_out_dir,
                const std::string& init_chg,
                const int nlocal,
                const std::string& ks_solver,
                const std::string& device,
                const std::vector<double>& hubbard_u,
                const double uramping,
                const int occ_mat_ctrl,
                const int mixing_dftu
#ifdef __LCAO
                , const LCAO_Orbitals* orb = nullptr
#endif
                );

  private:

#ifdef __LCAO
    const LCAO_Orbitals* ptr_orb_ = nullptr;
    std::vector<double> orb_cutoff_;

  public:
    /// read-only accessors for state needed by DFTU_LCAO free functions
    const std::vector<double>& get_orb_cutoff() const { return orb_cutoff_; }
    const LCAO_Orbitals* get_ptr_orb() const { return ptr_orb_; }
#endif
};

#endif
