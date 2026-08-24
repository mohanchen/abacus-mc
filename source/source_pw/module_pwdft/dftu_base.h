#ifndef DFTU_BASE_H
#define DFTU_BASE_H

#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge_mixing.h"

#include <string>
#include <vector>


class Plus_U_Base
{
  public:
    Plus_U_Base();
    ~Plus_U_Base();

  public:
    /// allocate relevant data structures (base part, no LCAO types)
    void init_base(UnitCell& cell,
                   const int npol,
                   const int nspin,
                   const std::vector<int>& orbital_corr,
                   const bool yukawa_potential,
                   const double yukawa_lambda,
                   const std::string& global_readin_dir,
                   const std::string& global_out_dir,
                   const std::string& init_chg,
                   const int nlocal,
                   const bool gamma_only_local,
                   const std::string& ks_solver,
                   const bool cal_force,
                   const bool cal_stress,
                   const std::string& device,
                   const int kpar);

    void uramping_update();
    bool u_converged();

    // mohan change these parameters to static, 2025-11-07
    static std::vector<double> U;
    static std::vector<double> U0;
    static std::vector<int> orbital_corr;
    static double uramping;
    static int omc;
    static int mixing_dftu;
    static int nspin;

    // --- Accessors for static data ---

    static double get_hubbard_u(int it) { return U[it]; }
    static double get_hubbard_u0(int it) { return U0[it]; }
    static int get_num_u_types() { return static_cast<int>(U.size()); }
    static int get_orbital_corr(int it) { return orbital_corr[it]; }
    static bool has_correlated_orbital(int it) { return orbital_corr[it] != -1; }
    static const int* get_orbital_corr_data() { return orbital_corr.data(); }

    static double get_energy() { return energy_u; }
    static void set_energy(const double &e) { energy_u = e; }
    static void set_double_energy() { energy_u *= 2.0; }

  private:
    static double energy_u;

  protected:
    int cal_type = 3;
    std::string global_readin_dir;
    std::string global_out_dir;
    double yukawa_lambda = 0.0;
    std::string init_chg;
    int npol = 1;
    int nlocal = 0;
    bool gamma_only_local = false;
    std::string ks_solver;
    bool cal_force = false;
    bool cal_stress = false;
    std::string device;
    int kpar = 1;

    // transform between iwt index and it, ia, L, N and m index
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>
        iatlnmipol2iwt;

  public:
    /// interface for PW base
    /// calculate the local occupation number matrix for PW based wave functions
    void cal_occ_pw(const int iter,
                    const void* psi_in,
                    const ModuleBase::matrix& wg_in,
                    const UnitCell& cell,
                    Charge_Mixing* p_chgmix,
                    const int* isk);

    /// get effective potential pointer for the given spin channel (PW basis)
    ///
    /// nspin=1: isk is ignored, returns &eff_pot_pw[0]
    /// nspin=2: isk selects spin-up (0) or spin-down (1) half of the
    ///          split layout [all_up | all_dn]
    /// nspin=4: isk is ignored, returns &eff_pot_pw[0] (all Pauli blocks)
    const std::complex<double>* get_eff_pot_pw_spin(const int isk) const
    {
        if (nspin == 2 && isk == 1)
        {
            return eff_pot_pw.data() + eff_pot_pw.size() / 2;
        }
        return eff_pot_pw.data();
    }

    /// get size of effective potential for a single spin channel (PW basis)
    ///
    /// nspin=1: full array size
    /// nspin=2: half of the total (one spin channel in split layout)
    /// nspin=4: full array size (all Pauli blocks are packed together)
    int get_size_eff_pot_pw_spin() const
    {
        return (nspin == 2) ? static_cast<int>(eff_pot_pw.size() / 2)
                            : static_cast<int>(eff_pot_pw.size());
    }

    /// get effective potential matrix for PW base (per-atom, raw index)
    /// @deprecated Use get_eff_pot_pw_spin() for nspin-aware access.
    [[deprecated("Use get_eff_pot_pw_spin() for nspin-aware access")]]
    const std::complex<double>* get_eff_pot_pw(const int iat) const
    {
        return &(eff_pot_pw[eff_pot_pw_index[iat]]);
    }

    int get_size_eff_pot_pw() const
    {
        return eff_pot_pw.size();
    }

    // dftu can be calculated only after locale has been initialed
    bool initialed_locale = false;

    bool is_locale_initialized() const { return initialed_locale; }
    void mark_locale_initialized() { initialed_locale = true; }
    void mark_locale_dirty() { initialed_locale = false; }

    static bool is_mixing_enabled() { return mixing_dftu != 0; }
    static void enable_mixing() { mixing_dftu = 1; }

  private:
    void copy_locale(const UnitCell& ucell);
    void zero_locale(const UnitCell& ucell);
    void mix_locale(const UnitCell& ucell, const double& mixing_beta);
    void set_locale(const UnitCell& ucell);

    std::vector<std::complex<double>> eff_pot_pw;
    std::vector<int> eff_pot_pw_index;
    std::vector<double> uom_array;
    std::vector<double> uom_save;

    // Yukawa-related members (base part, no LCAO dependency)
    double lambda = 0.0;
    std::vector<std::vector<std::vector<std::vector<double>>>> Fk;
    std::vector<std::vector<std::vector<double>>> U_Yukawa;
    std::vector<std::vector<std::vector<double>>> J_Yukawa;

  public:
    /// get occupation matrix element locale[iat][l][n][spin](m1,m2)
    double get_locale(const int iat, const int l, const int n, const int spin,
                     const int m1, const int m2) const
    {
        return locale[iat][l][n][spin](m1, m2);
    }

    /// set occupation matrix element locale[iat][l][n][spin](m1,m2)
    void set_locale(const int iat, const int l, const int n, const int spin,
                   const int m1, const int m2, const double val)
    {
        locale[iat][l][n][spin](m1, m2) = val;
    }

    /// get flat occupation matrix for an atom's correlated orbital.
    /// nspin=1: fills occ with locale[iat][l][0][0] data
    /// nspin=2: fills occ with interleaved locale[iat][l][0][0] and [1] data
    /// nspin=4: fills occ with locale[iat][l][0][0] data (all 4 Pauli blocks)
    void get_locale_flat(const int iat, const int l, std::vector<double>& occ) const;

    /// set flat occupation matrix for an atom's correlated orbital (write-back)
    void set_locale_flat(const int iat, const int l, const int spin,
                        const std::vector<double>& occ);

    // local occupancy matrix of the correlated subspace
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> locale;
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> locale_save;

    //=============================================================
    // In dftu_io.cpp
    // For reading/writing/broadcasting/copying relevant data structures
    //=============================================================
  public:
    void output(const UnitCell& ucell,
                bool out_chg,
                const std::string& global_out_dir,
                int nspin,
                int npol);

  private:
    void write_occup_m(const UnitCell& ucell,
                       std::ofstream& ofs,
                       bool diag,
                       int nspin,
                       int npol);
    void read_occup_m(const UnitCell& ucell,
                      const std::string& fn,
                      const std::string& init_chg,
                      int nspin,
                      int npol);
    void local_occup_bcast(const UnitCell& ucell,
                           int nspin,
                           int npol);

    //=============================================================
    // In dftu_yukawa.cpp
    // Relevant for calculating U using Yukawa potential
    //=============================================================

  public:
    static bool Yukawa;
};


#endif
