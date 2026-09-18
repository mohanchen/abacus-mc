#include "chg_dmr.h"

#include <functional>
#include <vector>

#include "source_base/global_function.h"
#include "source_base/module_mixing/mixing.h"
#include "source_base/tool_quit.h"
#include "source_estate/module_dm/density_matrix.h"

namespace module_charge
{
namespace
{

/**
 * @brief Two-beta plain step for the magnetic density matrix.
 *
 * The first nnr elements are the charge channel (beta), the next nnr
 * elements are the magnetization channel (beta_mag).
 *
 * @param out      mixed output, length 2 * nnr
 * @param in       mixed input, length 2 * nnr
 * @param sres     residual, length 2 * nnr
 * @param nnr      number of DMR elements per spin channel
 * @param beta     mixing beta for the charge channel
 * @param beta_mag mixing beta for the magnetization channel
 */
void twobeta_step(double* out,
                  const double* in,
                  const double* sres,
                  const int nnr,
                  const double beta,
                  const double beta_mag)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
    for (int i = 0; i < nnr; ++i)
    {
        out[i] = in[i] + beta * sres[i];
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
    for (int i = nnr; i < 2 * nnr; ++i)
    {
        out[i] = in[i] + beta_mag * sres[i];
    }
}

} // namespace

void init_mixing_dmr(Base_Mixing::Mixing* mixing,
                     Base_Mixing::Mixing_Data& mdata,
                     const int nnr,
                     const MixingConfig& cfg)
{
    if (mixing == nullptr)
    {
        ModuleBase::WARNING_QUIT("module_charge::init_mixing_dmr", "mixing pointer is null");
    }
    if (nnr <= 0)
    {
        ModuleBase::WARNING_QUIT("module_charge::init_mixing_dmr", "nnr must be > 0");
    }

    const int dmr_nspin = (cfg.nspin == 2) ? 2 : 1;
    // DMR mixing currently supports only the real-space convergence threshold.
    if (cfg.scf_thr_type == 1)
    {
        ModuleBase::WARNING_QUIT("module_charge::init_mixing_dmr",
                                 "This Mixing of Density Matrix is not supported for PW basis yet");
    }
    else if (cfg.scf_thr_type == 2)
    {
        mixing->init_mixing_data(mdata, nnr * dmr_nspin, sizeof(double));
    }

    // Clear the history counters while keeping the allocated storage.
    mdata.reset();
}

template <typename TK>
void mix_dmr(elecstate::DensityMatrix<TK, double>* dm,
             Base_Mixing::Mixing* mixing,
             Base_Mixing::Mixing_Data& mdata,
             const MixingConfig& cfg)
{
    if (dm == nullptr)
    {
        ModuleBase::WARNING_QUIT("module_charge::mix_dmr", "density matrix pointer is null");
    }
    if (mixing == nullptr)
    {
        ModuleBase::WARNING_QUIT("module_charge::mix_dmr", "mixing pointer is null");
    }
    if (cfg.nspin != 1 && cfg.nspin != 2 && cfg.nspin != 4)
    {
        ModuleBase::WARNING_QUIT("module_charge::mix_dmr", "unsupported nspin, require 1, 2 or 4");
    }

    std::vector<hamilt::HContainer<double>*> dmr = dm->get_DMR_vector();
    std::vector<std::vector<double>>& dmr_save = dm->get_DMR_save();

    double* dmr_in = nullptr;
    double* dmr_out = nullptr;
    if (cfg.nspin == 1 || cfg.nspin == 4)
    {
        dmr_in = dmr_save[0].data();
        dmr_out = dmr[0]->get_wrapper();
        mixing->push_data(mdata, dmr_in, dmr_out, nullptr, false);
        mixing->mix_data(mdata, dmr_out);
    }
    else // cfg.nspin == 2
    {
        // Magnetic density matrix: up/down channels are transformed into
        // charge/magnetization channels before mixing and back afterwards.
        const int nnr = dmr[0]->get_nnr();
        std::vector<double> dmr_mag(nnr * cfg.nspin, 0.0);
        std::vector<double> dmr_mag_save(nnr * cfg.nspin, 0.0);

        // Transfer the current DMR into the charge/magnetization layout.
        double* dmr_up = dmr[0]->get_wrapper();
        double* dmr_down = dmr[1]->get_wrapper();
        for (int ir = 0; ir < nnr; ++ir)
        {
            dmr_mag[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }
        // Transfer the saved DMR into the charge/magnetization layout.
        dmr_up = dmr_save[0].data();
        dmr_down = dmr_save[1].data();
        for (int ir = 0; ir < nnr; ++ir)
        {
            dmr_mag_save[ir] = dmr_up[ir] + dmr_down[ir];
            dmr_mag_save[ir + nnr] = dmr_up[ir] - dmr_down[ir];
        }

        dmr_in = dmr_mag_save.data();
        dmr_out = dmr_mag.data();
        const double beta = cfg.mixing_beta;
        const double beta_mag = cfg.mixing_beta_mag;
        std::function<void(double*, const double*, const double*)> twobeta
            = [nnr, beta, beta_mag](double* out, const double* in, const double* sres) {
                  twobeta_step(out, in, sres, nnr, beta, beta_mag);
              };
        // No Kerker screening in DMR mixing.
        mixing->push_data(mdata, dmr_in, dmr_out, nullptr, twobeta, false);
        mixing->mix_data(mdata, dmr_out);

        // Transform the mixed charge/magnetization channels back to up/down.
        dmr_up = dmr[0]->get_wrapper();
        dmr_down = dmr[1]->get_wrapper();
        ModuleBase::GlobalFunc::ZEROS(dmr_up, nnr);
        ModuleBase::GlobalFunc::ZEROS(dmr_down, nnr);
        for (int ir = 0; ir < nnr; ++ir)
        {
            dmr_up[ir] = 0.5 * (dmr_mag[ir] + dmr_mag[ir + nnr]);
            dmr_down[ir] = 0.5 * (dmr_mag[ir] - dmr_mag[ir + nnr]);
        }
    }
}

template void mix_dmr<double>(elecstate::DensityMatrix<double, double>* dm,
                              Base_Mixing::Mixing* mixing,
                              Base_Mixing::Mixing_Data& mdata,
                              const MixingConfig& cfg);
template void mix_dmr<std::complex<double>>(elecstate::DensityMatrix<std::complex<double>, double>* dm,
                                            Base_Mixing::Mixing* mixing,
                                            Base_Mixing::Mixing_Data& mdata,
                                            const MixingConfig& cfg);

} // namespace module_charge
