#include "charge_mixing.h"
#include "chg_drho.h"
#include "chg_precond.h"
#include "chg_uspp.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"

void Charge_Mixing::mix_rho_recip(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_rho_recip");
    ModuleBase::timer::start("Charge_Mixing", "mix_rho_recip");

    const int nspin = PARAM.inp.nspin;
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

    if ( PARAM.globalv.double_grid)
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
    auto inner_product = [this](std::complex<double>* rhog1, std::complex<double>* rhog2)
    {
        return module_charge::inner_product_recip_hartree(
            rhog1, rhog2, *this->rhopw, this->cfg_, *this->omega, *this->tpiba);
    };

    // DIIS Mixing Only for smooth part, while high_frequency part is mixed by plain mixing method.
    if (nspin == 1)
    {
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        auto screen = [this](std::complex<double>* p) {
            module_charge::kerker_screen_recip(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
    }
    else if (nspin == 2)
    {
        // magnetic density
        const int npw = this->rhopw->npw;
        // rhog_mag[is*ngmc] and rhog_mag_save[is*ngmc], zero-initialized
        std::vector<std::complex<double>> rhog_mag(npw * nspin);
        std::vector<std::complex<double>> rhog_mag_save(npw * nspin);
        // get rhog_mag[is*ngmc] and rhog_mag_save[is*ngmc]
        for (int ig = 0; ig < npw; ig++)
        {
            rhog_mag[ig] = chr->rhog[0][ig] + chr->rhog[1][ig];
            rhog_mag_save[ig] = chr->rhog_save[0][ig] + chr->rhog_save[1][ig];
        }
        for (int ig = 0; ig < npw; ig++)
        {
            rhog_mag[ig + npw] = chr->rhog[0][ig] - chr->rhog[1][ig];
            rhog_mag_save[ig + npw] = chr->rhog_save[0][ig] - chr->rhog_save[1][ig];
        }
        //
        rhog_in = rhog_mag_save.data();
        rhog_out = rhog_mag.data();
        //
        auto screen = [this](std::complex<double>* p) {
            module_charge::kerker_screen_recip(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        auto twobeta_mix = this->make_twobeta_mix<std::complex<double>>(2 * npw, npw);
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
        // get rhog[is][ngmc] from rhog_mag[is*ngmc]
        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(chr->rhog[is], npw);
        }
        for (int ig = 0; ig < npw; ig++)
        {
            chr->rhog[0][ig] = 0.5 * (rhog_mag[ig] + rhog_mag[ig+npw]);
            chr->rhog[1][ig] = 0.5 * (rhog_mag[ig] - rhog_mag[ig+npw]);
        }
        // sync rhogs_out so merge_dgrid has the latest smooth part
        if ( PARAM.globalv.double_grid)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                rhogs_out[ig] = chr->rhog[0][ig];
                rhogs_out[ig + npw] = chr->rhog[1][ig];
            }
        }
    }
    else if (nspin == 4 && PARAM.inp.mixing_angle <= 0)
    {
        // normal broyden mixing for {rho, mx, my, mz}
        rhog_in = rhogs_in;
        rhog_out = rhogs_out;
        const int npw = this->rhopw->npw;
        auto screen = [this](std::complex<double>* p) { // use old one
            module_charge::kerker_screen_recip(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        auto twobeta_mix = this->make_twobeta_mix<std::complex<double>>(4 * npw, npw);
        this->mixing->push_data(this->rho_mdata, rhog_in, rhog_out, screen, twobeta_mix, true);
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhog_out);
    }
    else if (nspin == 4 && PARAM.inp.mixing_angle > 0)
    {
        // special broyden mixing for {rho, |m|} proposed by J. Phys. Soc. Jpn. 82 (2013) 114706
        // here only consider the case of mixing_angle = 1, which mean only change |m| and keep angle fixed
        // old support see mix_rho_recip()
        if ( PARAM.globalv.double_grid)
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
        auto screen = [this](std::complex<double>* p) { // use old one
            module_charge::kerker_screen_recip(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        auto twobeta_mix = this->make_twobeta_mix<std::complex<double>>(2 * npw, npw);
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

    if ( PARAM.globalv.double_grid)
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
    if (nspin == 4 && PARAM.inp.mixing_angle > 0)
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
            // rhodpw is the same as rhopw for ! PARAM.globalv.double_grid
            this->rhodpw->recip_to_real<std::complex<double>,double,base_device::DEVICE_CPU>(chr->rhog[is], chr->rho[is]);
        }
    }
    // For kinetic energy density
    if ((XC_Functional::get_ked_flag()) && mixing_tau)
    {
        std::vector<std::complex<double>> kin_g(nspin * rhodpw->npw);
        std::vector<std::complex<double>> kin_g_save(nspin * rhodpw->npw);
        // FFT to get kin_g and kin_g_save
        for (int is = 0; is < nspin; ++is)
        {
            rhodpw->real2recip(chr->kin_r[is], &kin_g[is * rhodpw->npw]);
            rhodpw->real2recip(chr->kin_r_save[is], &kin_g_save[is * rhodpw->npw]);
        }
        // RAII owners for the smooth / high-frequency parts on the double grid;
        // raw pointers below alias these vectors when double_grid is on, or
        // alias kin_g[_save] directly when double_grid is off so the mixing
        // mutates the dense buffer in place.
        std::vector<std::complex<double>> tau_sg_in;
        std::vector<std::complex<double>> tau_sg_out;
        std::vector<std::complex<double>> tau_hf_in;
        std::vector<std::complex<double>> tau_hf_out;
        // for smooth part, for ! PARAM.globalv.double_grid only have this part
        std::complex<double>* taugs_in = nullptr;
        std::complex<double>* taugs_out = nullptr;
        // for high frequency part
        std::complex<double>* taughf_in = nullptr;
        std::complex<double>* taughf_out = nullptr;
        if ( PARAM.globalv.double_grid)
        {
            // divide into smooth part and high_frequency part
            const int npw_smooth = this->rhopw->npw;
            const int npw_dense = this->rhodpw->npw;
            tau_sg_in.resize(nspin * npw_smooth);
            tau_hf_in.resize(nspin * (npw_dense - npw_smooth));
            tau_sg_out.resize(nspin * npw_smooth);
            tau_hf_out.resize(nspin * (npw_dense - npw_smooth));
            module_charge::split_dgrid(kin_g_save.data(), tau_sg_in, tau_hf_in,
                                        nspin, npw_smooth, npw_dense);
            module_charge::split_dgrid(kin_g.data(), tau_sg_out, tau_hf_out,
                                        nspin, npw_smooth, npw_dense);
            taugs_in = tau_sg_in.data();
            taughf_in = tau_hf_in.data();
            taugs_out = tau_sg_out.data();
            taughf_out = tau_hf_out.data();
        }
        else
        {
            taugs_in = kin_g_save.data();
            taugs_out = kin_g.data();
        }

        // Note: there is no kerker modification for tau because I'm not sure
        // if we should have it. If necessary we can try it in the future.
        this->mixing->push_data(this->tau_mdata, taugs_in, taugs_out, nullptr, false);

        this->mixing->mix_data(this->tau_mdata, taugs_out);

        if ( PARAM.globalv.double_grid)
        {
            // simple mixing for high_frequencies
            const int ndimhf = (this->rhodpw->npw - this->rhopw->npw) * nspin;
            this->mixing_highf->plain_mix(taughf_out, taughf_in, taughf_out, ndimhf, nullptr);

            // combine smooth part and high_frequency part;
            // tau_sg_* / tau_hf_* vectors are released automatically at scope exit
            module_charge::merge_dgrid(kin_g.data(), tau_sg_out, tau_hf_out,
                                        nspin, this->rhopw->npw, this->rhodpw->npw);
        }

        // kin_g to kin_r
        for (int is = 0; is < nspin; is++)
        {
            rhodpw->recip2real(&kin_g[is * rhodpw->npw], chr->kin_r[is]);
        }
    }

    ModuleBase::timer::end("Charge_Mixing", "mix_rho_recip");
    return;
}

void Charge_Mixing::mix_rho_real(Charge* chr)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_rho_real");
    ModuleBase::timer::start("Charge_Mixing", "mix_rho_real");

    const int nspin = PARAM.inp.nspin;
    assert(nspin==1 || nspin==2 || nspin==4);

    double* rhor_in=nullptr;
    double* rhor_out=nullptr;

    if (nspin == 1)
    {
        rhor_in = chr->rho_save[0];
        rhor_out = chr->rho[0];
        auto screen = [this](double* p) {
            module_charge::kerker_screen_real(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, true);
        auto inner_product = [this](double* rho1, double* rho2)
        {
            return module_charge::inner_product_real(rho1, rho2, *this->rhopw, this->cfg_);
        };
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
    }
    else if (nspin == 2)
    {
        // magnetic density
        const int nrxx = this->rhopw->nrxx;
        // rho_mag[is*nnrx] and rho_mag_save[is*nnrx], zero-initialized
        std::vector<double> rho_mag(nrxx * nspin);
        std::vector<double> rho_mag_save(nrxx * nspin);
        // get rho_mag[is*nnrx] and rho_mag_save[is*nnrx]
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho_mag[ir] = chr->rho[0][ir] + chr->rho[1][ir];
            rho_mag_save[ir] = chr->rho_save[0][ir] + chr->rho_save[1][ir];
        }
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho_mag[ir + nrxx] = chr->rho[0][ir] - chr->rho[1][ir];
            rho_mag_save[ir + nrxx] = chr->rho_save[0][ir] - chr->rho_save[1][ir];
        }
        //
        rhor_in = rho_mag_save.data();
        rhor_out = rho_mag.data();
        auto screen = [this](double* p) {
            module_charge::kerker_screen_real(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        auto twobeta_mix = this->make_twobeta_mix<double>(2 * nrxx, nrxx);
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        auto inner_product = [this](double* rho1, double* rho2)
        {
            return module_charge::inner_product_real(rho1, rho2, *this->rhopw, this->cfg_);
        };
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
        // get new rho[is][nrxx] from rho_mag[is*nrxx]
        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(chr->rho[is], nrxx);
            //ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
        }
        for (int ir = 0; ir < nrxx; ir++)
        {
            chr->rho[0][ir] = 0.5 * (rho_mag[ir] + rho_mag[ir+nrxx]);
            chr->rho[1][ir] = 0.5 * (rho_mag[ir] - rho_mag[ir+nrxx]);
        }
    }
    else if (nspin == 4 && PARAM.inp.mixing_angle <= 0)
    {
        // normal broyden mixing for {rho, mx, my, mz}
        rhor_in = chr->rho_save[0];
        rhor_out = chr->rho[0];
        const int nrxx = this->rhopw->nrxx;
        auto screen = [this](double* p) {
            module_charge::kerker_screen_real(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        auto twobeta_mix = this->make_twobeta_mix<double>(4 * nrxx, nrxx);
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        auto inner_product = [this](double* rho1, double* rho2)
        {
            return module_charge::inner_product_real(rho1, rho2, *this->rhopw, this->cfg_);
        };
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);
    }
    else if (nspin == 4 && PARAM.inp.mixing_angle > 0)
    {
        // special broyden mixing for {rho, |m|} proposed by J. Phys. Soc. Jpn. 82 (2013) 114706
        // here only consider the case of mixing_angle = 1, which mean only change |m| and keep angle fixed
        const int nrxx = this->rhopw->nrxx;
        // rho_magabs and rho_magabs_save, zero-initialized
        std::vector<double> rho_magabs(nrxx * 2);
        std::vector<double> rho_magabs_save(nrxx * 2);
        // calculate rho_magabs and rho_magabs_save
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

        auto screen = [this](double* p) {
            module_charge::kerker_screen_real(this->cfg_, this->rhopw, *this->tpiba, p);
        };
        auto twobeta_mix = this->make_twobeta_mix<double>(2 * nrxx, nrxx);
        this->mixing->push_data(this->rho_mdata, rhor_in, rhor_out, screen, twobeta_mix, true);
        auto inner_product = [this](double* rho1, double* rho2)
        {
            return module_charge::inner_product_real(rho1, rho2, *this->rhopw, this->cfg_);
        };
        this->mixing->cal_coef(this->rho_mdata, inner_product);
        this->mixing->mix_data(this->rho_mdata, rhor_out);

        // use new |m| and angle to update {mx, my, mz}
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
    if ((XC_Functional::get_ked_flag()) && mixing_tau)
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

    const int nspin = PARAM.inp.nspin;
    assert(nspin==1 || nspin==2 || nspin==4);

    // the charge before mixing.
    const int nrxx = chr->rhopw->nrxx;
    std::vector<double> rho123(nspin * nrxx);
    for (int is = 0; is < nspin; ++is)
    {
        if (is == 0 || is == 3 || !PARAM.globalv.domag_z)
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
    if ((XC_Functional::get_ked_flag()) && mixing_tau)
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
    if (PARAM.inp.scf_thr_type == 1)
    {
        mix_rho_recip(chr);
    }
    else if (PARAM.inp.scf_thr_type == 2)
    {
        mix_rho_real(chr);
    }
    // ---------------------------------------------------

    // mohan add 2012-06-05
    // rho_save is the charge before mixing
    for (int is = 0; is < nspin; ++is)
    {
        if (is == 0 || is == 3 || !PARAM.globalv.domag_z)
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

    if ((XC_Functional::get_ked_flag()) && mixing_tau)
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

    if (new_e_iteration) 
    {
        new_e_iteration = false;
    }

    ModuleBase::timer::end("Charge_Mixing", "mix_rho");
    return;
}
