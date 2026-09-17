#ifndef CHARGE_MIXING_H
#define CHARGE_MIXING_H
#include "charge.h"
#include "chg_mix_cfg.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_base/module_mixing/mixing.h"
#include "source_base/module_mixing/plain_mixing.h"
#include <functional>

class Charge_Mixing
{
  /// Charge_Mixing class
  /// This class is used to mix charge density, kinetic energy density and real-space density matrix
  /// This Charge_Mixing class offers the following interfaces:
  /// 1. set_mixing() to set all private mixing parameters
  /// 2. init_mixing() to initialize mixing, including allocating memory for mixing data and reset mixing
  /// 3. mix_rho() to mix charge density
  /// 4. mix_dmr() to mix real-space density matrix
  /// how to use it:
  /// you can (re)start a mixing by calling set_mixing() and init_mixing() before calling mix_rho() or mix_dmr()

  public:
    Charge_Mixing();
    ~Charge_Mixing();

    /**
     * @brief Set all private mixing parameters from an aggregated config
     * @param cfg mixing parameters and runtime globals (nspin, scf_thr_type, double_grid)
     * @param omega_in omega for non-linear core correction
     * @param tpiba_in 2*pi/beta for non-linear core correction
     */
    void set_mixing(const MixingConfig& cfg,
                    double& omega_in,
                    double& tpiba_in);

    void close_kerker_gg0() { mixing_gg0 = 0.0; mixing_gg0_mag = 0.0; }
    void conserve_setting() { mixing_beta = 0.01; mixing_beta_mag = 0.04; }
    /**
     * @brief initialize mixing, including constructing mixing and allocating memory for mixing data
     * @brief this function should be called at eachiterinit()
     */
    void init_mixing();

    /**
     * @brief allocate memory of dmr_mdata
     * @param nnr size of real-space density matrix
     */
    void allocate_mixing_dmr(const int nnr);

    /**
     * @brief charge mixing
     * @param chr pointer of Charge object
     */
    void mix_rho(Charge* chr);

    /**
     * @brief density matrix mixing, only for LCAO
     * @param DM pointer of DensityMatrix object
     */
    void mix_dmr(elecstate::DensityMatrix<double, double>* DM);
    void mix_dmr(elecstate::DensityMatrix<std::complex<double>, double>* DM);

    /**
     * @brief allocate memory of uom_mdata
     * @param uom_size size of DFT+U occupation matrix
     */
    void allocate_mixing_uom(int size_uom);

    /**
     * @brief DFT+U occupation matrix mixing
     * @param uom_in output occupation matrix
     * @param uom_save_in input occupation matrix
     */
    void mix_uom(std::vector<double>& uom_in, std::vector<double>& uom_save_in);

    /**
     * @brief reset mixing, actually we only call init_mixing() to reset mixing instead of this function
     */
    void mix_reset();
    
    /**
     * @brief Set the smooth and dense grids
     * @param rhopw_in smooth grid
     * @param rhodpw_in dense grid when double grid is used, otherwise same as rhopw
     */
    void set_rhopw(ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis* rhodpw_in);

    // extracting parameters normally these parameters will not be used outside charge mixing
    // while Exx is using them as well as some other places
    const std::string& get_mixing_mode() const {return mixing_mode;}
    double get_mixing_beta() const {return mixing_beta;}
    int get_mixing_ndim() const {return mixing_ndim;}
    double get_mixing_gg0() const {return mixing_gg0;}
    Base_Mixing::Mixing* get_mixing() const {return mixing;}

    /**
     * @brief read-only access to the aggregated mixing config set by set_mixing()
     */
    const MixingConfig& get_mixing_config() const {return cfg_;}

    // for mixing restart
    int mixing_restart_step = 0; //which step to restart mixing during SCF, always equal to scf_namx except for the mixing restart
    int mixing_restart_count = 0; // the number of restart mixing during SCF. Do not set mixing_restart_count as bool since I want to keep some flexibility in the future
    int mixing_restart_last = 0; // the label of mixing restart step, store the step number of the last mixing restart

    // to calculate the slope of drho curve during SCF, which is used to determine if SCF oscillate
    bool if_scf_oscillate(const int iteration, const double drho, const int iternum_used, const double threshold);
    
  private:
  
    // mixing_data
    Base_Mixing::Mixing* mixing = nullptr; ///< Mixing object to mix charge density, kinetic energy density and compensation density
    Base_Mixing::Mixing_Data rho_mdata;    ///< Mixing data for charge density
    Base_Mixing::Mixing_Data tau_mdata;    ///< Mixing data for kinetic energy density
    Base_Mixing::Mixing_Data nhat_mdata;   ///< Mixing data for compensation density
    Base_Mixing::Mixing_Data dmr_mdata;    ///< Mixing data for real space density matrix
    Base_Mixing::Mixing_Data uom_mdata;    ///< Mixing data for DFT+U occupation matrix
    Base_Mixing::Plain_Mixing* mixing_highf = nullptr; ///< The high_frequency part is mixed by plain mixing method.

    //======================================
    // private mixing parameters
    //======================================
    MixingConfig cfg_;                 ///< aggregated mixing config, also holds nspin/scf_thr_type/double_grid
    std::string mixing_mode = "broyden"; ///< mixing mode: "plain", "broyden", "pulay"
    double mixing_beta = 0.8;            ///< mixing beta for density
    double mixing_beta_mag = 1.6;        ///< mixing beta for magnetism
    int mixing_ndim = 8;                 ///< mixing ndim for broyden and pulay
    double mixing_gg0 = 0.0;             ///< mixing gg0 for Kerker screen
    bool mixing_tau = false;             ///< whether to use tau mixing
    double mixing_gg0_mag = 0.0;         ///< mixing gg0 for Kerker screen for magnetism
    double mixing_gg0_min = 0.1;         ///< minimum kerker coefficient
    double mixing_angle = 0.0;           ///< mixing angle for nspin=4
    bool mixing_dmr = false;             ///< whether to mixing real space density matrix
    double* omega = nullptr;                  ///< omega for non-linear core correction
    double* tpiba = nullptr;                  ///< 2*pi/beta for non-linear core correction
    double* tpiba2 = nullptr;                 ///< 2*pi/beta^2 for non-linear core correction
    std::vector<double> _drho_history; ///< history of drho used to determine the oscillation, size is scf_nmax
    
    bool new_e_iteration = true;

    ModulePW::PW_Basis* rhopw = nullptr;  ///< smooth grid
    ModulePW::PW_Basis* rhodpw = nullptr; ///< dense grid, same as rhopw for ncpp.

    /**
     * @brief charge mixing for reciprocal space
     * @param chr pointer of Charge object
     */
    void mix_rho_recip(Charge* chr);

    /**
     * @brief charge mixing for real space
     * @param chr pointer of Charge object
     */
    void mix_rho_real(Charge* chr);

    /**
     * @brief two-beta mixing functor: mix the first `nunit` elements with
     * mixing_beta and the rest (nunit..total) with mixing_beta_mag. Used for
     * magnetic cases (nspin==2/4) where the charge channel and the magnetism
     * channels use different betas. Replaces the duplicated local lambdas.
     * @tparam T element type, double (real space) or std::complex<double> (reciprocal)
     */
    template <typename T>
    std::function<void(T*, const T*, const T*)> make_twobeta_mix(const int total, const int nunit)
    {
        return [this, total, nunit](T* out, const T* in, const T* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < nunit; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = nunit; i < total; ++i)
            {
                out[i] = in[i] + this->mixing_beta_mag * sres[i];
            }
        };
    }

    /**
     * @brief divide rho/tau to smooth and high frequency parts
     * @param data_d dense data
     * @param data_s smooth data
     * @param data_hf high frequency data = dense data - smooth data
     *
     */
    void divide_data(std::complex<double>* data_d, std::complex<double>*& data_s, std::complex<double>*& data_hf);
    /**
     * @brief gather smooth and high frequency parts to rho/tau
     * @param data_d dense data
     * @param data_s smooth data
     * @param data_hf high frequency data = dense data - smooth data
     *  
     */
    void combine_data(std::complex<double>* data_d, std::complex<double>*& data_s, std::complex<double>*& data_hf);
    /**
     * @brief clean smooth and high frequency parts
     * @param data_d dense data
     * @param data_s smooth data
     * @param data_hf high frequency data = dense data - smooth data
     *
     */
    void clean_data(std::complex<double>*& data_s, std::complex<double>*& data_hf);
};

#endif
