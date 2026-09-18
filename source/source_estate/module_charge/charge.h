#ifndef CHARGE_H
#define CHARGE_H

#include <vector>

#include "source_base/complexmatrix.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"
// #include "source_estate/fp_energy.h"
#include "source_base/parallel_grid.h"

//a forward declaration of UnitCell
class UnitCell;

// Electron Charge Density
class Charge
{

  public:

    Charge();
    ~Charge();

    // rho/rhog/kin_r views alias the vector-backed _space_* storage, so
    // copying a Charge would duplicate dangling pointers into another
    // object's vector buffer. Forbid copies until a deep copy is needed.
    Charge(const Charge&) = delete;
    Charge& operator=(const Charge&) = delete;

    //==========================================================
    // MEMBER VARIABLES :
    // init_chg : "atomic" or "file"
    // NAME : total number of electrons
    // NAME : rho (nspin,ncxyz), the charge density in real space
    // NAME : rho_save (nspin,ncxyz), for charge mixing
    // NAME : rhog, charge density in G space
    // NAME : rhog_save, chage density in G space
    // NAME : rho_core [nrxx], the core charge in real space
    // NAME : rhog_core [ngm], the core charge in reciprocal space
    //==========================================================

    double **rho = nullptr;
    double **rho_save = nullptr;

    std::complex<double> **rhog = nullptr;
    std::complex<double> **rhog_save = nullptr;

    double **kin_r = nullptr; // kinetic energy density in real space, for meta-GGA
    double **kin_r_save = nullptr; // kinetic energy density in real space, for meta-GGA
    const Parallel_Grid* pgrid = nullptr;

  private:

    // Underlying contiguous storage backing the public rho/rhog/kin_r views.
    // Each buffer holds nspin rows; rho[is] points at _space_rho.data()+is*nrxx.
    // Owned here as std::vector so the storage self-manages (no raw new/delete).
    std::vector<double> _space_rho;
    std::vector<double> _space_rho_save;
    std::vector<std::complex<double>> _space_rhog;
    std::vector<std::complex<double>> _space_rhog_save;
    std::vector<double> _space_kin_r;
    std::vector<double> _space_kin_r_save;

  public:

    double **nhat = nullptr; //compensation charge for PAW
    double **nhat_save = nullptr; //compensation charge for PAW
                                 // wenfei 2023-09-05

    double *rho_core = nullptr;
    std::complex<double> *rhog_core = nullptr;

    int prenspin = 1;

    void set_rhopw(ModulePW::PW_Basis* rhopw_in);

    /**
     * @brief Init charge density from file or atomic pseudo-wave-functions
     *
     * @param eferm_iout [out] fermi energy to be initialized
     * @param ucell [in] unit cell
     * @param strucFac [in] structure factor
     * @param symm [in] symmetry
     * @param klist [in] k points list if needed
     * @param wfcpw [in] PW basis for wave function if needed
     */
    void init_rho(const UnitCell& ucell,
                  const Parallel_Grid& pgrid,
                  const ModuleBase::ComplexMatrix& strucFac,
                  ModuleSymmetry::Symmetry& symm,
                  const void* klist = nullptr,
                  const void* wfcpw = nullptr);

    // mohan add 2025-12-02
    bool kin_density() const;

    void allocate(const int &nspin_in, const bool kin_den);

    void renormalize_rho();

    double sum_rho() const;

    void save_rho_before_sum_band();

    double cal_rho2ne(const double *rho_in) const;

    void check_rho(); // to check whether the charge density is normal

    void init_final_scf(); //LiuXh add 20180619

    void set_omega(double* omega_in){this->omega_ = omega_in;};

    // mohan add 2021-02-20
    int nrxx=0; // number of r vectors in this processor
    int nxyz = 0; // total number of r vectors
    int ngmc=0; // number of g vectors in this processor
    int nspin=0; // number of spins
    ModulePW::PW_Basis* rhopw = nullptr;// When double_grid is used, rhopw = rhodpw (dense grid)
    bool cal_elf = false; // whether to calculate electron localization function (ELF)

  private:

    void destroy();    // free arrays  liuyu 2023-03-12

    double* omega_ = nullptr; // omega for non-linear core correction

    bool allocate_rho;

    bool allocate_rho_final_scf; // LiuXh add 20180606
};

#endif // charge
