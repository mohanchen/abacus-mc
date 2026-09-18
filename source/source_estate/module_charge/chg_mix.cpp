#include "chg_mix.h"
#include "chg_drho.h"
#include "chg_precond.h"
#include "chg_rho_detail.h"
#include "chg_tau.h"
#include "chg_uspp.h"

#include <functional>
#include <memory>

#include "source_io/module_parameter/parameter.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_base/module_mixing/pulay_mixing.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_hamilt/module_xc/xc_functional.h"

Charge_Mixing::Charge_Mixing()
{
    // unique_ptr members default-construct to nullptr
}

Charge_Mixing::~Charge_Mixing()
{
    // unique_ptr members (mixing, mixing_highf) are released automatically
}

void Charge_Mixing::set_mixing(const MixingConfig& cfg,
                               double& omega_in,
                               double& tpiba_in)
{
    // store the aggregated config; init_mixing/mix_rho read nspin,
    // scf_thr_type and double_grid from it instead of PARAM/GlobalV.
    this->cfg_ = cfg;
    // get private mixing parameters
    this->mixing_mode = cfg.mixing_mode;
    this->mixing_beta = cfg.mixing_beta;
    this->mixing_beta_mag = cfg.mixing_beta_mag;
    this->mixing_ndim = cfg.mixing_ndim;
    this->mixing_gg0 = cfg.mixing_gg0;

    this->mixing_gg0_mag = cfg.mixing_gg0_mag;
    this->mixing_gg0_min = cfg.mixing_gg0_min;
    this->mixing_angle = cfg.mixing_angle;
    this->mixing_dmr = cfg.mixing_dmr;
    this->omega = &omega_in;
    this->tpiba = &tpiba_in;
    // check the paramters
    if (this->mixing_beta > 1.0 || this->mixing_beta < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta to [0.0, 1.0]!");
    }
    if (cfg.nspin >= 2 && this->mixing_beta_mag < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta_mag >= 0.0!");
    }

    if (!(this->mixing_mode == "plain" || this->mixing_mode == "broyden" || this->mixing_mode == "pulay"))
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    // print into running.log
    //GlobalV::ofs_running << "\n\n";
    GlobalV::ofs_running << "\n";
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                 "
           "   |" << std::endl;
    GlobalV::ofs_running << " | Setup charge mixing parameters                                  "
           "   |" << std::endl;
    GlobalV::ofs_running << " |                                                                 "
           "   |" << std::endl;
    GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<" << std::endl;
    GlobalV::ofs_running << "\n";


    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_type", this->mixing_mode);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_beta", this->mixing_beta);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0", this->mixing_gg0);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0_min", this->mixing_gg0_min);

    if (cfg.nspin==2 || cfg.nspin==4)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_beta_mag", this->mixing_beta_mag);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0_mag", this->mixing_gg0_mag);
    }
    if (this->mixing_angle > 0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_angle", this->mixing_angle);
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_ndim", this->mixing_ndim);

    return;
}

void Charge_Mixing::init_mixing()
{
    // this init should be called at the 1-st iteration of each scf loop

    ModuleBase::TITLE("Charge_Mixing", "init_mixing");
    ModuleBase::timer::start("Charge_Mixing", "init_mixing");

    // (re)construct mixing object
    if (this->mixing_mode == "broyden")
    {
        this->mixing = std::make_unique<Base_Mixing::Broyden_Mixing>(this->mixing_ndim, this->mixing_beta);
    }
    else if (this->mixing_mode == "plain")
    {
        this->mixing = std::make_unique<Base_Mixing::Plain_Mixing>(this->mixing_beta);
    }
    else if (this->mixing_mode == "pulay")
    {
        this->mixing = std::make_unique<Base_Mixing::Pulay_Mixing>(this->mixing_ndim, this->mixing_beta);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    if ( this->cfg_.double_grid)
    {
        // ONLY smooth part of charge density is mixed by specific mixing method
        // The high_frequency part is mixed by plain mixing method.
        this->mixing_highf = std::make_unique<Base_Mixing::Plain_Mixing>(this->mixing_beta);
    }

    // allocate memory for mixing data, if exists, free it first and then allocate new memory
    // initailize rho_mdata
    if (this->cfg_.scf_thr_type == 1)
    {
        if (this->cfg_.nspin == 4 && this->mixing_angle > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * 2,
                                        sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * this->cfg_.nspin,
                                        sizeof(std::complex<double>));
        }
    }
    else
    {
        if (this->cfg_.nspin == 4 && this->mixing_angle > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * 2, sizeof(double));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * this->cfg_.nspin, sizeof(double));
        }
    }

    // initailize tau_mdata
    if ((XC_Functional::get_ked_flag()) && cfg_.mixing_tau)
    {
        if (this->cfg_.scf_thr_type == 1)
        {
            this->mixing->init_mixing_data(this->tau_mdata,
                                           this->rhopw->npw * this->cfg_.nspin,
                                           sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->tau_mdata, this->rhopw->nrxx * this->cfg_.nspin, sizeof(double));
        }
    }

    ModuleBase::timer::end("Charge_Mixing", "init_mixing");

    return;
}

void Charge_Mixing::set_rhopw(ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis* rhodpw_in)
{
    this->rhopw = rhopw_in;
    this->rhodpw = rhodpw_in;
}

void Charge_Mixing::mix_reset()
{
    this->mixing->reset();
    this->rho_mdata.reset();
    // initailize tau_mdata
    if ((XC_Functional::get_ked_flag()) && cfg_.mixing_tau)
    {
        this->tau_mdata.reset();
    }
}

bool Charge_Mixing::if_scf_oscillate(const int iteration, const double drho,
                                     const int iternum_used, const double threshold)
{
    ModuleBase::TITLE("Charge_Mixing", "if_scf_oscillate");

    if(this->_drho_history.size() == 0)
    {
        this->_drho_history.resize(PARAM.inp.scf_nmax);
    }

    // add drho into history
    this->_drho_history[iteration - 1] = drho;

    if(threshold >= 0) // close the function
    {
        return false;
    }

    // check if the history is long enough
    if(iteration < iternum_used + this->mixing_restart_last)
    {
        return false;
    }

    // calculate the slope of the last iternum_used iterations' drho
    double slope = 0.0;

    // Least Squares Method
    // this part is too short, so I do not design it as a free function in principle
    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    for (int i = iteration - iternum_used; i < iteration; i++)
    {
        sumX += i;
        sumY += std::log10(this->_drho_history[i]);
        sumXY += i * std::log10(this->_drho_history[i]);
        sumXX += i * i;
    }
    double numerator = iternum_used * sumXY - sumX * sumY;
    double denominator = iternum_used * sumXX - sumX * sumX;
    if (denominator == 0) {
        return false;
    }
    slope =  numerator / denominator;

    // if the slope is less than the threshold, return true
    if(slope > threshold)
    {
        return true;
    }

    return false;
}

void Charge_Mixing::allocate_mixing_uom(int uom_size)
{
    ModuleBase::TITLE("Charge_Mixing", "allocate_mixing_uom");
    ModuleBase::timer::start("Charge_Mixing", "allocate_mixing_uom");
    // For nspin=2, uom_size already includes both spin channels
    // (uterm_mat.size() = pot_index * 2 for nspin=2)
    // So uom_fold should always be 1
    this->mixing->init_mixing_data(this->uom_mdata, uom_size, sizeof(double));
    this->uom_mdata.reset();
    ModuleBase::timer::end("Charge_Mixing", "allocate_mixing_uom");
    return;
}

void Charge_Mixing::mix_uom(std::vector<double>& uom_in, std::vector<double>& uom_save_in)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_uom");
    ModuleBase::timer::start("Charge_Mixing", "mix_uom");
    double* uom_value_out = uom_in.data();
    double* uom_value_in = uom_save_in.data();
    // For all nspin cases, uom_array layout is already fully sized
    // and mixing operates on the entire array
    this->mixing->push_data(this->uom_mdata, uom_value_in, uom_value_out, nullptr, false);
    this->mixing->mix_data(this->uom_mdata, uom_value_out);
    ModuleBase::timer::end("Charge_Mixing", "mix_uom");
#ifdef __MPI
    // Synchronize mixed uom across all ranks to prevent divergence
    // after multiple Pulay steps (same pattern as mix_dmr)
    Parallel_Common::bcast_double(uom_in.data(), uom_in.size());
#endif
    return;
}

void Charge_Mixing::mix_rho_recip(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_rho_recip");
    ModuleBase::timer::start("Charge_Mixing", "mix_rho_recip");

    const int nspin = cfg_.nspin;
    assert(nspin==1 || nspin==2 || nspin==4);

    std::complex<double>* rhog_in = nullptr;
    std::complex<double>* rhog_out = nullptr;
    // RAII owners for the smooth / high-frequency parts on the double grid.
    // The raw pointers below alias these vectors when double_grid is on,
    // or alias chr->rhog[_save][0] directly when double_grid is off so the
    // mixing still mutates chr in place.
    std::vector<std::complex<double>> rho_sg_in;
    std::vector<std::complex<double>> rho_sg_out;
    std::vector<std::complex<double>> rho_hf_in;
    std::vector<std::complex<double>> rho_hf_out;
    // for smooth part
    std::complex<double>* rhogs_in = nullptr;
    std::complex<double>* rhogs_out = nullptr;
    // for high_frequency part
    std::complex<double>* rhoghf_in = nullptr;
    std::complex<double>* rhoghf_out = nullptr;

    if ( cfg_.double_grid)
    {
        // divide into smooth part and high_frequency part
        const int npw_smooth = this->rhopw->npw;
        const int npw_dense = this->rhodpw->npw;
        rho_sg_in.resize(nspin * npw_smooth);
        rho_hf_in.resize(nspin * (npw_dense - npw_smooth));
        rho_sg_out.resize(nspin * npw_smooth);
        rho_hf_out.resize(nspin * (npw_dense - npw_smooth));
        module_charge::split_dgrid(chr->rhog_save[0], rho_sg_in, rho_hf_in,
                                    nspin, npw_smooth, npw_dense);
        module_charge::split_dgrid(chr->rhog[0], rho_sg_out, rho_hf_out,
                                    nspin, npw_smooth, npw_dense);
        rhogs_in = rho_sg_in.data();
        rhoghf_in = rho_hf_in.data();
        rhogs_out = rho_sg_out.data();
        rhoghf_out = rho_hf_out.data();
    }
    else
    {
        rhogs_in = chr->rhog_save[0];
        rhogs_out = chr->rhog[0];
    }

    //  inner_product_recip_hartree is a hartree-like sum, unit is Ry
    std::function<double(std::complex<double>*, std::complex<double>*)> inner_product
        = [this](std::complex<double>* rhog1, std::complex<double>* rhog2)
    {
        return module_charge::inner_product_recip_hartree(
            rhog1, rhog2, *this->rhopw, this->cfg_, *this->omega, *this->tpiba);
    };

    // Kerker screening functor, shared by all nspin branches
    std::function<void(std::complex<double>*)> screen = [this](std::complex<double>* p) {
        module_charge::kerker_screen_recip(this->cfg_, this->rhopw, *this->tpiba, p);
    };

    // DIIS Mixing Only for smooth part, while high_frequency part is mixed by plain mixing method.
    if (nspin == 1)
    {
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
    }
    else if (nspin == 2)
    {
        // magnetic density
        const int npw = this->rhopw->npw;
        std::vector<std::complex<double>> rhog_mag(npw * nspin);
        std::vector<std::complex<double>> rhog_mag_save(npw * nspin);
        module_charge::detail::pack_rho_mag(rhog_mag.data(), chr->rhog[0], chr->rhog[1], npw);
        module_charge::detail::pack_rho_mag(rhog_mag_save.data(), chr->rhog_save[0], chr->rhog_save[1], npw);
        //
        rhog_in = rhog_mag_save.data();
        rhog_out = rhog_mag.data();
        std::function<void(std::complex<double>*, const std::complex<double>*,
            const std::complex<double>*)> twobeta_mix
            = module_charge::detail::make_twobeta_mix<std::complex<double>>(
                2 * npw, npw, this->mixing_beta, this->mixing_beta_mag);
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
        // get rhog[is][ngmc] from rhog_mag[is*ngmc]
        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(chr->rhog[is], npw);
        }
        module_charge::detail::unpack_rho_mag(chr->rhog[0], chr->rhog[1], rhog_mag.data(), npw);
        // sync rhogs_out so merge_dgrid has the latest smooth part
        if ( cfg_.double_grid)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                rhogs_out[ig] = chr->rhog[0][ig];
                rhogs_out[ig + npw] = chr->rhog[1][ig];
            }
        }
    }
    else if (nspin == 4 && cfg_.mixing_angle <= 0)
    {
        // normal broyden mixing for {rho, mx, my, mz}
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        const int npw = this->rhopw->npw;
        std::function<void(std::complex<double>*, const std::complex<double>*,
            const std::complex<double>*)> twobeta_mix
            = module_charge::detail::make_twobeta_mix<std::complex<double>>(
                4 * npw, npw, this->mixing_beta, this->mixing_beta_mag);
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
    }
    else if (nspin == 4 && cfg_.mixing_angle > 0)
    {
        // special broyden mixing for {rho, |m|} proposed by J. Phys. Soc. Jpn. 82 (2013) 114706
        // here only consider the case of mixing_angle = 1, which mean only change |m| and keep angle fixed
        // old support see mix_rho_recip()
        if ( cfg_.double_grid)
        {
            ModuleBase::WARNING_QUIT("Charge_Mixing", "double_grid is not supported for new mixing method yet.");
        }
        // rho_magabs and rho_magabs_save, zero-initialized
        const int nrxx = this->rhopw->nrxx;
        std::vector<double> rho_magabs(nrxx);
        std::vector<double> rho_magabs_save(nrxx);
        // calculate rho_magabs and rho_magabs_save
        for (int ir = 0; ir < nrxx; ir++)
        {
            // |m| for rho
            rho_magabs[ir] = std::sqrt(chr->rho[1][ir] * chr->rho[1][ir]
            + chr->rho[2][ir] * chr->rho[2][ir]
            + chr->rho[3][ir] * chr->rho[3][ir]);
            // |m| for rho_save
            rho_magabs_save[ir] = std::sqrt(chr->rho_save[1][ir] * chr->rho_save[1][ir]
            + chr->rho_save[2][ir] * chr->rho_save[2][ir]
            + chr->rho_save[3][ir] * chr->rho_save[3][ir]);
        }
        // rhog_magabs and rhog_magabs_save, zero-initialized
        const int npw = this->rhopw->npw;
        std::vector<std::complex<double>> rhog_magabs(npw * 2);
        std::vector<std::complex<double>> rhog_magabs_save(npw * 2);
        // calculate rhog_magabs and rhog_magabs_save
        for (int ig = 0; ig < npw; ig++)
        {
            rhog_magabs[ig] = chr->rhog[0][ig]; // rho
            rhog_magabs_save[ig] = chr->rhog_save[0][ig]; // rho_save
        }
        // FT to get rhog_magabs and rhog_magabs_save
        this->rhopw->real2recip(rho_magabs.data(), rhog_magabs.data() + this->rhopw->npw);
        this->rhopw->real2recip(rho_magabs_save.data(), rhog_magabs_save.data() + this->rhopw->npw);
        //
        rhog_in = rhog_magabs_save.data();
        rhog_out = rhog_magabs.data();
        std::function<void(std::complex<double>*, const std::complex<double>*,
            const std::complex<double>*)> twobeta_mix
            = module_charge::detail::make_twobeta_mix<std::complex<double>>(
                2 * npw, npw, this->mixing_beta, this->mixing_beta_mag);
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
        // get new |m| in real space using FT
        this->rhopw->recip2real(rhog_magabs.data() + this->rhopw->npw, rho_magabs.data());
        // use new |m| and angle to update {mx, my, mz}
        for (int ig = 0; ig < npw; ig++)
        {
            chr->rhog[0][ig] = rhog_magabs[ig]; // rhog
            double norm = std::sqrt(chr->rho[1][ig] * chr->rho[1][ig]
                    + chr->rho[2][ig] * chr->rho[2][ig]
                    + chr->rho[3][ig] * chr->rho[3][ig]);
            if (std::abs(norm) < 1e-10)
            {
                continue;
            }
            double rescale_tmp = rho_magabs[npw + ig] / norm;
            chr->rho[1][ig] *= rescale_tmp;
            chr->rho[2][ig] *= rescale_tmp;
            chr->rho[3][ig] *= rescale_tmp;
        }
    }

    if ( cfg_.double_grid)
    {
        // plain mixing for high_frequencies
        const int ndimhf = (this->rhodpw->npw - this->rhopw->npw) * nspin;
        this->mixing_highf->plain_mix(rhoghf_out, rhoghf_in, rhoghf_out, ndimhf, nullptr);

        // combine smooth part and high_frequency part;
        // rho_sg_* / rho_hf_* vectors are released automatically at scope exit
        module_charge::merge_dgrid(chr->rhog[0], rho_sg_out, rho_hf_out,
                                    nspin, this->rhopw->npw, this->rhodpw->npw);
    }

    // rhog to rho
    if (nspin == 4 && cfg_.mixing_angle > 0)
    {
        // only tranfer rhog[0]
        // do not support double_grid, use rhopw directly
        chr->rhopw->recip2real(chr->rhog[0], chr->rho[0]);
    }
    else
    {
        for (int is = 0; is < nspin; is++)
        {
            // use rhodpw for double_grid
            // rhodpw is the same as rhopw for ! cfg_.double_grid
            this->rhodpw->recip_to_real<std::complex<double>, double,
                base_device::DEVICE_CPU>(chr->rhog[is], chr->rho[is]);
        }
    }
    // For kinetic energy density
    if ((XC_Functional::get_ked_flag()) && cfg_.mixing_tau)
    {
        module_charge::detail::mix_tau_recip(chr, nspin, cfg_.double_grid,
                      this->rhopw, this->rhodpw,
                      this->mixing, this->tau_mdata, this->mixing_highf);
    }

    ModuleBase::timer::end("Charge_Mixing", "mix_rho_recip");
    return;
}

void Charge_Mixing::mix_rho_real(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_rho_real");
    ModuleBase::timer::start("Charge_Mixing", "mix_rho_real");

    const int nspin = cfg_.nspin;
    assert(nspin==1 || nspin==2 || nspin==4);

    double* rhor_in=nullptr;
    double* rhor_out=nullptr;

    std::function<void(double*)> screen = [this](double* p) {
        module_charge::kerker_screen_real(this->cfg_, this->rhopw, *this->tpiba, p);
    };
    std::function<double(double*, double*)> inner_product = [this](double* rho1, double* rho2)
    {
        return module_charge::inner_product_real(rho1, rho2, *this->rhopw, this->cfg_);
    };

    if (nspin == 1)
    {
        rhor_in = chr->rho_save[0];
        rhor_out = chr->rho[0];
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
    }
    else if (nspin == 2)
    {
        // magnetic density
        const int nrxx = this->rhopw->nrxx;
        std::vector<double> rho_mag(nrxx * nspin);
        std::vector<double> rho_mag_save(nrxx * nspin);
        module_charge::detail::pack_rho_mag(rho_mag.data(), chr->rho[0], chr->rho[1], nrxx);
        module_charge::detail::pack_rho_mag(rho_mag_save.data(), chr->rho_save[0], chr->rho_save[1], nrxx);
        //
        rhor_in = rho_mag_save.data();
        rhor_out = rho_mag.data();
        std::function<void(double*, const double*, const double*)> twobeta_mix
            = module_charge::detail::make_twobeta_mix<double>(2 * nrxx, nrxx, this->mixing_beta, this->mixing_beta_mag);
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
        // get new rho[is][nrxx] from rho_mag[is*nrxx]
        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(chr->rho[is], nrxx);
        }
        module_charge::detail::unpack_rho_mag(chr->rho[0], chr->rho[1], rho_mag.data(), nrxx);
    }
    else if (nspin == 4 && cfg_.mixing_angle <= 0)
    {
        rhor_in = chr->rho_save[0];
        rhor_out = chr->rho[0];
        const int nrxx = this->rhopw->nrxx;
        std::function<void(double*, const double*, const double*)> twobeta_mix
            = module_charge::detail::make_twobeta_mix<double>(4 * nrxx, nrxx, this->mixing_beta, this->mixing_beta_mag);
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
    }
    else if (nspin == 4 && cfg_.mixing_angle > 0)
    {
        // real-space version of the {rho, |m|} broyden mixing
        const int nrxx = this->rhopw->nrxx;
        std::vector<double> rho_magabs(nrxx * 2);
        std::vector<double> rho_magabs_save(nrxx * 2);
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho_magabs[ir] = chr->rho[0][ir]; // rho
            rho_magabs_save[ir] = chr->rho_save[0][ir]; // rho_save
            // |m| for rho
            rho_magabs[nrxx + ir] = std::sqrt(chr->rho[1][ir] * chr->rho[1][ir]
                    + chr->rho[2][ir] * chr->rho[2][ir]
                    + chr->rho[3][ir] * chr->rho[3][ir]);
            // |m| for rho_save
            rho_magabs_save[nrxx + ir] = std::sqrt(chr->rho_save[1][ir] * chr->rho_save[1][ir]
                    + chr->rho_save[2][ir] * chr->rho_save[2][ir]
                    + chr->rho_save[3][ir] * chr->rho_save[3][ir]);
        }
        rhor_in = rho_magabs_save.data();
        rhor_out = rho_magabs.data();

        std::function<void(double*, const double*, const double*)> twobeta_mix
            = module_charge::detail::make_twobeta_mix<double>(2 * nrxx, nrxx, this->mixing_beta, this->mixing_beta_mag);
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);

        for (int ir = 0; ir < nrxx; ir++)
        {
            chr->rho[0][ir] = rho_magabs[ir]; // rho
            double norm = std::sqrt(chr->rho[1][ir] * chr->rho[1][ir]
                    + chr->rho[2][ir] * chr->rho[2][ir]
                    + chr->rho[3][ir] * chr->rho[3][ir]);

            if (norm < 1e-10)
            {
                continue;
            }
            double rescale_tmp = rho_magabs[nrxx + ir] / norm;
            chr->rho[1][ir] *= rescale_tmp;
            chr->rho[2][ir] *= rescale_tmp;
            chr->rho[3][ir] *= rescale_tmp;
        }
    }

    double *taur_out=nullptr;
    double *taur_in=nullptr;
    if ((XC_Functional::get_ked_flag()) && cfg_.mixing_tau)
    {
        taur_in = chr->kin_r_save[0];
        taur_out = chr->kin_r[0];
        // Note: there is no kerker modification for tau because I'm not sure
        // if we should have it. If necessary we can try it in the future.
        this->mixing->push_data(this->tau_mdata, taur_in, taur_out, nullptr, false);

        this->mixing->mix_data(this->tau_mdata, taur_out);
    }

    ModuleBase::timer::end("Charge_Mixing", "mix_rho_real");
    return;
}


void Charge_Mixing::mix_rho(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_rho");
    ModuleBase::timer::start("Charge_Mixing", "mix_rho");

    const int nspin = cfg_.nspin;
    assert(nspin==1 || nspin==2 || nspin==4);

    // the charge before mixing.
    const int nrxx = chr->rhopw->nrxx;
    std::vector<double> rho123(nspin * nrxx);
    for (int is = 0; is < nspin; ++is)
    {
        if (is == 0 || is == 3 || !cfg_.domag_z)
        {
            double* rho123_is = rho123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                rho123_is[ir] = chr->rho[is][ir];
            }
        }
    }
    std::vector<double> kin_r123;
    if ((XC_Functional::get_ked_flag()) && cfg_.mixing_tau)
    {
        kin_r123.resize(nspin * nrxx);
        for (int is = 0; is < nspin; ++is)
        {
            double* kin_r123_is = kin_r123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                kin_r123_is[ir] = chr->kin_r[is][ir];
            }
        }
    }
    // --------------------Mixing Body--------------------
    if (cfg_.scf_thr_type == 1)
    {
        mix_rho_recip(chr);
    }
    else if (cfg_.scf_thr_type == 2)
    {
        mix_rho_real(chr);
    }
    // ---------------------------------------------------

    // mohan add 2012-06-05
    // rho_save is the charge before mixing
    for (int is = 0; is < nspin; ++is)
    {
        if (is == 0 || is == 3 || !cfg_.domag_z)
        {
            double* rho123_is = rho123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                chr->rho_save[is][ir] = rho123_is[ir];
            }
        }
    }

    if ((XC_Functional::get_ked_flag()) && cfg_.mixing_tau)
    {
        for (int is = 0; is < nspin; ++is)
        {
            double* kin_r123_is = kin_r123.data() + is * nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                chr->kin_r_save[is][ir] = kin_r123_is[ir];
            }
        }
    }

    ModuleBase::timer::end("Charge_Mixing", "mix_rho");
    return;
}
