#include "charge.h"
#include "chg_drho.h"
#include "chg_drho_detail.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cassert>

namespace module_charge
{

// Charge residual between chr->rho and chr->rho_save, normalized per electron.
double cal_drho(Charge* chr,
                const double nelec,
                const ModulePW::PW_Basis& rhopw,
                const MixingConfig& cfg,
                const double omega,
                const double tpiba)
{
    assert(chr != nullptr);
    ModuleBase::TITLE("module_charge", "cal_drho");
    ModuleBase::timer::start("module_charge", "cal_drho");
    const int nspin = cfg.nspin;
    assert(nspin==1 || nspin==2 || nspin==4);
    double drho = 0.0;

    if (cfg.scf_thr_type == 1)
    {
        for (int is = 0; is < nspin; ++is)
        {
            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
            chr->rhopw->real2recip(chr->rho[is], chr->rhog[is]);

            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
            chr->rhopw->real2recip(chr->rho_save[is], chr->rhog_save[is]);
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        std::vector<std::complex<double>> drhog(nspin * rhopw.npw);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
        for (int is = 0; is < nspin; ++is)
        {
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                drhog[is * rhopw.npw + ig] = chr->rhog[is][ig] - chr->rhog_save[is][ig];
            }
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
        drho = module_charge::detail::inner_product_recip_rho(
            drhog.data(), drhog.data(), rhopw, cfg, omega, tpiba);
    }
    else
    {
        // Note: Maybe it is wrong.
        //       The inner_product_real function (L1-norm) is different from that (L2-norm) in mixing.
        for (int is = 0; is < nspin; is++)
        {
            if (is != 0 && is != 3 && cfg.domag_z)
            {
                continue;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : drho)
#endif
            for (int ir = 0; ir < rhopw.nrxx; ir++)
            {
                drho += std::abs(chr->rho[is][ir] - chr->rho_save[is][ir]);
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(drho);
#endif
        assert(nelec != 0);
        assert(omega > 0);
        assert(rhopw.nxyz > 0);
        drho *= omega / static_cast<double>(rhopw.nxyz);
        drho /= nelec;
    }

    ModuleBase::timer::end("module_charge", "cal_drho");
    return drho;
}

// Kinetic-energy-density residual between chr->kin_r and chr->kin_r_save.
double cal_dkin(Charge* chr,
                const double nelec,
                const ModulePW::PW_Basis& rhopw,
                const MixingConfig& cfg,
                const double omega)
{
    assert(chr != nullptr);
    if (!(XC_Functional::get_ked_flag()))
    {
        return 0.0;
    };
    ModuleBase::TITLE("module_charge", "cal_dkin");
    ModuleBase::timer::start("module_charge", "cal_dkin");
    double dkin = 0.0;

    // Get dkin from kin_r and kin_r_save for PW and LCAO both, which is different from drho.
    for (int is = 0; is < cfg.nspin; is++)
    {
        if (is != 0 && is != 3 && cfg.domag_z)
        {
            continue;
        }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : dkin)
#endif
        for (int ir = 0; ir < rhopw.nrxx; ir++)
        {
            dkin += std::abs(chr->kin_r[is][ir] - chr->kin_r_save[is][ir]);
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(dkin);
#endif
    assert(nelec != 0);
    assert(omega > 0);
    assert(rhopw.nxyz > 0);
    dkin *= omega / static_cast<double>(rhopw.nxyz);
    dkin /= nelec;

    ModuleBase::timer::end("module_charge", "cal_dkin");
    return dkin;
}

} // namespace module_charge

namespace module_charge
{
namespace detail
{

double inner_product_recip_rho(const std::complex<double>* rho1,
                               const std::complex<double>* rho2,
                               const ModulePW::PW_Basis& rhopw,
                               const MixingConfig& cfg,
                               const double omega,
                               const double tpiba)
{
    assert(rho1 != nullptr);
    assert(rho2 != nullptr);
    assert(cfg.nspin == 1 || cfg.nspin == 2 || cfg.nspin == 4);
    ModuleBase::TITLE("Charge_Mixing", "recip_rho");
    ModuleBase::timer::start("Charge_Mixing", "recip_rho");

    const int nspin = cfg.nspin;
    std::vector<const std::complex<double>*> rhog1(nspin);
    std::vector<const std::complex<double>*> rhog2(nspin);
    for (int is = 0; is < nspin; is++)
    {
        rhog1[is] = rho1 + is * rhopw.npw;
        rhog2[is] = rho2 + is * rhopw.npw;
    }

    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (tpiba * tpiba);
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;

    auto part_of_noncolin = [&]()
    {
        double sum = 0.0;
        const int ig0 = rhopw.ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            if (ig == ig0) {continue;}
            sum += (conj(rhog1[0][ig]) * rhog2[0][ig]).real() / rhopw.gg[ig];
        }
        sum *= fac;
        return sum;
    };

    switch (nspin)
    {
    case 1:
        sum += part_of_noncolin();
        break;

    case 2: {
        // (1) First part of density error.
        const int ig0 = rhopw.ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            if (ig == ig0) {continue;}
            sum += (conj(rhog1[0][ig] + rhog1[1][ig]) * (rhog2[0][ig] + rhog2[1][ig])).real() / rhopw.gg[ig];
        }
        sum *= fac;

        if (cfg.gamma_only_pw)
        {
            sum *= 2.0;
        }

        // (2) Second part of density error.
        // including |G|=0 term.
        double sum2 = 0.0;

        // The G=0 component is the ig_gge0-th element of the local G-list on the
        // rank that owns it, not necessarily element 0: the local G-list is built
        // by scanning (x,y) sticks in grid order, so element 0 is the first plane
        // wave of the first owned stick. Using a hardcoded index 0 made the inner
        // product partition-dependent for pools with more than one rank.
        if (ig0 >= 0)
        {
            sum2 += fac2 * (conj(rhog1[0][ig0] - rhog1[1][ig0]) * (rhog2[0][ig0] - rhog2[1][ig0])).real();
        }

        double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
        for (int ig = 0; ig < rhopw.npw; ig++)
        {
            if (ig == ig0) { continue; }
            mag += (conj(rhog1[0][ig] - rhog1[1][ig]) * (rhog2[0][ig] - rhog2[1][ig])).real();
        }
        mag *= fac2;

        if (cfg.gamma_only_pw)
        {
            mag *= 2.0;
        }

        sum2 += mag;
        sum += sum2;
        break;
    }
    case 4:
        // non-collinear spin, added by zhengdy
        if (!cfg.domag && !cfg.domag_z) {
            sum += part_of_noncolin();
        } else
        {
            // another part with magnetization
            const int ig0 = rhopw.ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                if (ig == ig0)
                {
                    continue;
                }
                sum += (conj(rhog1[0][ig]) * rhog2[0][ig]).real() / rhopw.gg[ig];
            }
            sum *= fac;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[1][ig0]) * rhog2[1][ig0]).real() + (conj(rhog1[2][ig0]) * rhog2[2][ig0]).real()
                          + (conj(rhog1[3][ig0]) * rhog2[3][ig0]).real());
            }
            double fac3 = fac2;
            if (cfg.gamma_only_pw)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                if (ig == ig0) {
                    continue;
}
                sum += fac3
                       * ((conj(rhog1[1][ig]) * rhog2[1][ig]).real() + (conj(rhog1[2][ig]) * rhog2[2][ig]).real()
                          + (conj(rhog1[3][ig]) * rhog2[3][ig]).real());
            }
        }
        break;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif

    sum *= omega * 0.5;

    ModuleBase::timer::end("Charge_Mixing", "recip_rho");
    return sum;
}

} // namespace detail

double inner_product_real(const double* rho1,
                          const double* rho2,
                          const ModulePW::PW_Basis& rhopw,
                          const MixingConfig& cfg)
{
    assert(rho1 != nullptr);
    assert(rho2 != nullptr);
    double rnorm = 0.0;
    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (cfg.nspin == 4 && cfg.mixing_angle > 0)
    {
        resize_tmp = 2;
    }

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : rnorm)
#endif
    for (int ir = 0; ir < rhopw.nrxx * cfg.nspin / resize_tmp; ++ir)
    {
        rnorm += rho1[ir] * rho2[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(rnorm);
#endif
    return rnorm;
}

// a Hartree-like inner product
double inner_product_recip_hartree(const std::complex<double>* rhog1,
                                   const std::complex<double>* rhog2,
                                   const ModulePW::PW_Basis& rhopw,
                                   const MixingConfig& cfg,
                                   const double omega,
                                   const double tpiba)
{
    assert(rhog1 != nullptr);
    assert(rhog2 != nullptr);
    assert(cfg.nspin == 1 || cfg.nspin == 2 || cfg.nspin == 4);
    ModuleBase::TITLE("Charge_Mixing", "recip_hartree");
    ModuleBase::timer::start("Charge_Mixing", "recip_hartree");

    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (tpiba * tpiba);
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;
    const int npw = rhopw.npw;

    // a lambda function for summing the charge density
    auto part_of_rho = [&]()
    {
        double sum = 0.0;
        const int ig0 = rhopw.ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            if (ig == ig0)
            {
                continue;
            }
            sum += (conj(rhog1[ig]) * rhog2[ig]).real() / rhopw.gg[ig];
        }
        sum *= fac;
        return sum;
    };

    if (cfg.nspin==1)
    {
        sum += part_of_rho();
    }
    else if (cfg.nspin==2)
    {
        // charge density part
        const int ig0 = rhopw.ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int ig = 0; ig < rhopw.npw; ++ig)
        {
            if (ig == ig0)
            {
                continue;
            }
            sum += (conj(rhog1[ig]) * (rhog2[ig])).real() / rhopw.gg[ig];
        }
        sum *= fac;

        if (cfg.gamma_only_pw)
        {
            sum *= 2.0;
        }

        // (2) Second part of density error.
        // including |G|=0 term.
        double sum2 = 0.0;

        // Same G=0 indexing remark as in inner_product_recip_rho: use ig_gge0
        // instead of a hardcoded index 0, otherwise the inner product (and hence
        // the DIIS mixing coefficients) depends on how the pool is divided.
        if (ig0 >= 0)
        {
            sum2 += fac2 * (conj(rhog1[ig0 + rhopw.npw]) * rhog2[ig0 + rhopw.npw]).real();
        }

        double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
        for (int ig = 0; ig < rhopw.npw; ig++)
        {
            if (ig == ig0) { continue; }
            mag += (conj(rhog1[ig + rhopw.npw]) * rhog2[ig + rhopw.npw]).real();
        }
        mag *= fac2;

        if (cfg.gamma_only_pw)
        {
            mag *= 2.0;
        }

        sum2 += mag;
        sum += sum2;
    }
    else if (cfg.nspin==4)
    {
        if (!cfg.domag && !cfg.domag_z)
        {
            sum += part_of_rho();
        }
        else if (cfg.mixing_angle <= 0)
        {
            // sum for tradtional mixing
            const int ig0 = rhopw.ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                if (ig == ig0) {continue;}
                sum += (conj(rhog1[ig]) * rhog2[ig]).real() / rhopw.gg[ig];
            }
            sum *= fac;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[ig0 + npw]) * rhog2[ig0 + npw]).real() + (conj(rhog1[ig0 + 2*npw]) * rhog2[ig0 + 2*npw]).real()
                          + (conj(rhog1[ig0 + 3*npw]) * rhog2[ig0 + 3*npw]).real());
            }
            double fac3 = fac2;
            if (cfg.gamma_only_pw)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                if (ig == ig0) {
                    continue;
}
                sum += fac3
                       * ((conj(rhog1[ig + npw]) * rhog2[ig + npw]).real() + (conj(rhog1[ig + 2*npw]) * rhog2[ig + 2*npw]).real()
                          + (conj(rhog1[ig + 3*npw]) * rhog2[ig + 3*npw]).real());
            }
        }
        else if (cfg.mixing_angle > 0)
        {
            // sum for angle mixing
            const int ig0 = rhopw.ig_gge0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                if (ig == ig0)
                {
                    continue;
                }
                sum += (conj(rhog1[ig]) * rhog2[ig]).real() / rhopw.gg[ig];
            }
            sum *= fac;
            if (ig0 > 0)
            {
                sum += fac2
                       * ((conj(rhog1[ig0 + rhopw.npw]) * rhog2[ig0 + rhopw.npw]).real());
            }
            double fac3 = fac2;
            if (cfg.gamma_only_pw)
            {
                fac3 *= 2.0;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                if (ig == ig0) {
                    continue;
}
                sum += fac3
                       * ((conj(rhog1[ig + rhopw.npw]) * rhog2[ig + rhopw.npw]).real());
            }
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif

    sum *= omega * 0.5;

    ModuleBase::timer::end("Charge_Mixing", "recip_hartree");

    return sum;
}

} // namespace module_charge
