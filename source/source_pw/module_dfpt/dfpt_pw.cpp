#include "dfpt_pw.h"
#include "dfpt_pw_impl.h"

#include "dfpt_hamilt_shift.h"
#include "dfpt_kq_basis.h"
#include "dfpt_metal.h"
#include "dfpt_pert.h"
#include "dfpt_phon.h"
#include "dfpt_pw_data.h"
#include "dfpt_q0.h"
#include "dfpt_rho.h"
#include "dfpt_stern.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/qlist.h"
#include "source_pw/module_pwdft/stru_fac.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace ModuleDFPT
{

DFPT_PW::Impl::Impl()
{
}

DFPT_PW::Impl::~Impl()
{
}

bool DFPT_PW::Impl::wired() const
{
    return pw_rho_ != nullptr && pw_wfc_ != nullptr;
}

DFPT_PW::DFPT_PW() : pimpl_(std::unique_ptr<Impl>(new Impl()))
{
}

DFPT_PW::~DFPT_PW()
{
}

bool DFPT_PW::get_with_u() const
{
    return pimpl_->data_.with_u();
}

bool DFPT_PW::get_u_active() const
{
    return pimpl_->data_.u_active();
}

double DFPT_PW::Impl::solve_displacement(int q_idx, int iat, int idir)
{
    ModuleBase::TITLE("DFPT_PW", "solve_displacement");
    ModuleBase::timer::start("DFPT_PW", "solve_displacement");
    if (!wired() || hamilt_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_PW", "solve_displacement");
        return 0.0;
    }
    const ModuleBase::Vector3<double> q_frac = data_.get_qvec(q_idx);
    const ModuleBase::Vector3<double> q_cart = q_frac * ucell_->G;
    const int nrxx = pw_rho_->nrxx;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();

    pert_.build_dv(q_idx, iat, idir, data_);
    rho_.reset_mixing(q_idx);
    // the previous perturbation's stored response must not leak into the
    // first iteration of this one
    data_.set_drho_g(q_idx, 0, std::vector<std::complex<double>>(pw_rho_->npw, std::complex<double>(0.0, 0.0)));

    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    bool converged = false;
    double residual = 0.0;
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);
    // last screened response potential (hoisted out of the loop: the 2n+1
    // accumulation below needs the converged v_sc of this displacement)
    std::vector<std::complex<double>> v_sc_r_last;
    for (int iter = 0; iter < max_iter_ && !converged; ++iter)
    {
        // the per-displacement SCF state (iter / residual / converged) is
        // local to this solve: the DFPT_PW_Data ledger is the per-(q,irrep)
        // outer-pass record kept by run(), and the final residual is
        // returned to the caller for that aggregation (B4)

        // ---- 1. screened response potential from the mixed input density:
        // q-shifted complex periodic amplitude on the shared grid, i.e. the
        // same convention as dv_rc (v_hartree_q acts on the q-shifted
        // coefficients; the XC kernel responds to Re/Im of the amplitude)
        std::vector<std::complex<double>> v_sc_r(nrxx, std::complex<double>(0.0, 0.0));
        const std::vector<std::complex<double>> drho_in_g = data_.get_drho_g(q_idx, 0);
        if (!drho_in_g.empty() && static_cast<int>(drho_in_g.size()) == pw_rho_->npw)
        {
            std::vector<std::complex<double>> dv_ha_g;
            rho_.v_hartree_q(q_cart, drho_in_g, dv_ha_g);
            std::vector<std::complex<double>> vh_r(nrxx);
            pw_rho_->recip2real(dv_ha_g.data(), vh_r.data());
            for (int ir = 0; ir < nrxx; ++ir)
            {
                v_sc_r[ir] = vh_r[ir];
            }
            if (xc_ != nullptr)
            {
                std::vector<std::complex<double>> a_r(nrxx);
                pw_rho_->recip2real(drho_in_g.data(), a_r.data());
                std::vector<std::complex<double>> b_r;
                xc_->apply(a_r, b_r);
                if (static_cast<int>(b_r.size()) == nrxx)
                {
                    for (int ir = 0; ir < nrxx; ++ir)
                    {
                        v_sc_r[ir] += b_r[ir];
                    }
                }
            }
            if (dbg)
            {
                double dh = 0.0;
                double dv = 0.0;
                for (int ig = 0; ig < pw_rho_->npw; ++ig)
                {
                    dh += std::norm(drho_in_g[ig]);
                }
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    dv += std::norm(v_sc_r[ir]);
                }
                std::cout << "DBG iter=" << iter << " |drho_in_g|=" << std::sqrt(dh) << " |v_sc_r|=" << std::sqrt(dv)
                          << std::endl;
            }
        }
        v_sc_r_last = v_sc_r;

        // ---- 2. Sternheimer solve of every occupied (k, band)
        for (int ik = 0; ik < nk; ++ik)
        {
            if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty())
            {
                if (dbg)
                {
                    std::cout << "DBG skip ik=" << ik << " no occ_kq" << std::endl;
                }
                continue; // no occupied states at k+q: nothing to solve
            }
            // dV_ext |psi_n> for all bands (dVloc convolution + dVnl_dtau)
            pert_.apply_dv(q_idx, ik, gs_psi_, data_);
            // screened response part |v_sc psi_n>
            std::vector<std::vector<std::complex<double>>> dv_sc;
            pert_.apply_vr(q_idx, ik, v_sc_r, gs_psi_, q_cart, dv_sc);
            if (ik != last_ik_ || last_q_ != q_idx)
            {
                hamilt_->set_context(q_cart, ik);
                last_ik_ = ik;
                if (dbg)
                {
                    std::cout << "DBG occ_kq nstates=" << occ_kq_[ik].size() << std::endl;
                    for (size_t m = 0; m < occ_kq_[ik].size(); ++m)
                    {
                        double nrm = 0.0;
                        for (size_t i = 0; i < occ_kq_[ik][m].size(); ++i)
                        {
                            nrm += std::norm(occ_kq_[ik][m][i]);
                        }
                        std::cout << "DBG   occ[" << m << "] |psi|^2=" << nrm << std::endl;
                    }
                    // kernel consistency: <psi_m|H(k+q)|psi_m> must equal
                    // eig(ikq, m); the eigenvalue used by set_shift below is
                    // the k-side one (equal only when H is assembled right)
                    for (size_t m = 0; m < occ_kq_[ik].size(); ++m)
                    {
                        hamilt_->set_shift(0.0);
                        std::vector<std::complex<double>> hp(occ_kq_[ik][m].size());
                        hamilt_->apply(occ_kq_[ik][m].data(), hp.data());
                        std::complex<double> dot(0.0, 0.0);
                        for (size_t i = 0; i < hp.size(); ++i)
                        {
                            dot += std::conj(occ_kq_[ik][m][i]) * hp[i];
                        }
                        std::cout << "DBG   <psi_" << m << "|H|psi_" << m << "> = " << dot.real() << " + i "
                                  << dot.imag() << "  (GS eig " << eig_(ikq_of_k_[ik], static_cast<int>(m)) << ")"
                                  << std::endl;
                        std::cout << "DBG   <psi_" << m << "|T+Vnl|psi_" << m
                                  << "> = " << hamilt_->debug_t_vnl(occ_kq_[ik][m]) << std::endl;
                        std::cout << "DBG   <psi_" << m << "|V(wfc-path)|psi_" << m
                                  << "> = " << hamilt_->debug_v_wfc(occ_kq_[ik][m]) << std::endl;
                    }
                }
            }
            for (int ib = 0; ib < nbands; ++ib)
            {
                if (!dfpt_band_occupied(wg_, ik, ib))
                {
                    continue; // unoccupied: no Sternheimer equation
                }
                std::vector<std::complex<double>> rhs = data_.get_dpsi(q_idx, ik, ib);
                if (rhs.empty() || static_cast<int>(dv_sc.size()) != nbands || rhs.size() != dv_sc[ib].size())
                {
                    if (dbg)
                    {
                        std::cout << "DBG skip solve ik=" << ik << " ib=" << ib << " rhs.size=" << rhs.size()
                                  << " dv_sc.size=" << dv_sc.size() << " dv_sc[ib].size="
                                  << (dv_sc.size() > static_cast<size_t>(ib) ? dv_sc[ib].size() : 999999) << std::endl;
                    }
                    continue;
                }
                // b = -(dV_ext + dV_sc)|psi_n>
                for (size_t i = 0; i < rhs.size(); ++i)
                {
                    rhs[i] = -(rhs[i] + dv_sc[ib][i]);
                }
                hamilt_->set_shift(eig_(ik, ib));
                std::vector<std::complex<double>> dpsi_out;
                double res = 0.0;
                stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr, dpsi_out, res);
                if (dbg)
                {
                    double nr = 0.0, nb2 = 0.0;
                    for (size_t i = 0; i < dpsi_out.size(); ++i)
                    {
                        nr += std::norm(dpsi_out[i]);
                        nb2 += std::norm(rhs[i]);
                    }
                    std::cout << "DBG solve ik=" << ik << " ib=" << ib << " eps=" << eig_(ik, ib) << " res=" << res
                              << " |dpsi|=" << std::sqrt(nr) << " |rhs|=" << std::sqrt(nb2)
                              << " finite=" << (std::isfinite(std::sqrt(nr)) ? 1 : 0) << std::endl;
                }
                data_.set_dpsi(q_idx, ik, ib, dpsi_out);
            }
        }

        // ---- 3. first-order density and mixing
        rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
        rho_.mix_drho(q_idx, data_);
        residual = rho_.get_residual(q_idx, data_);
        if (dbg)
        {
            std::cout << "DBG iter=" << iter << " residual=" << residual << " conv_thr=" << conv_thr_ << std::endl;
        }
        converged = (residual < conv_thr_);
    }
    // stash the converged screened potential and dpsi of this displacement
    // for the two-pass 2n+1 accumulation (term2 cross section needs
    // dV_ext^b + dV_sc^b and dpsi^b of every displacement)
    data_.set_vsc_r(iat, idir, v_sc_r_last);
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> disp(
            nk,
            std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                disp[ik][ib] = data_.get_dpsi(q_idx, ik, ib);
            }
        }
        data_.set_dpsi_disp(iat, idir, disp);
    }

    ModuleBase::timer::end("DFPT_PW", "solve_displacement");
    return residual;
}

void DFPT_PW::Impl::solve_pos_resp(int q_idx)
{
    ModuleBase::TITLE("DFPT_PW", "solve_pos_resp");
    ModuleBase::timer::start("DFPT_PW", "solve_pos_resp");
    // Y^a_{k,v} = P_c x_a|psi_{k,v}> through the Sternheimer equation
    //   (H(k) - eps_v) Y^a_v = P_c [H, x_a]|psi_v>,
    //   [H, x_a]|psi> = -(i/tpiba) dH/dk_a|psi>   (velocity form),
    // exactly the linear solve of QE dvpsi_e (whose rhs negation restores
    // P_c[H,x]psi from commutator_Hx_psi's [x,H] convention). dH/dk_a is the
    // pos_matrix velocity operator: the diagonal kinetic 2 tpiba^2 (k+G)_a
    // plus the separable projector derivative (build_vkb/build_vkb_dk). The
    // solved vector carries the complete conduction-space position response
    // and replaces the empty-eigenvector-truncated r-matrix contraction.
    if (!wired() || hamilt_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_PW", "solve_pos_resp");
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data_.get_qvec(q_idx) * ucell_->G;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    const double tpiba = ucell_->tpiba;
    const double tpiba2 = tpiba * tpiba;
    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();
    const bool dbg = (getenv("DFPT_DEBUG") != nullptr);

    for (int a = 0; a < 3; ++a)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> yvec(
            nk,
            std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik)
        {
            if (occ_kq_[ik].empty())
            {
                continue; // matches the displacement solve guard
            }
            if (last_q_ != q_idx || last_ik_ != ik)
            {
                hamilt_->set_context(q_cart, ik);
                last_q_ = q_idx;
                last_ik_ = ik;
            }
            const int npwk = pw_wfc_->npwk[ik];
            std::vector<ModuleBase::Vector3<double>> gk(npwk);
            for (int ig = 0; ig < npwk; ++ig)
            {
                gk[ig] = pw_wfc_->getgpluskcar(ik, ig);
            }
            // dH/dk_a|psi_b> for every band: diagonal kinetic part
            std::vector<std::vector<std::complex<double>>> vel(
                nbands,
                std::vector<std::complex<double>>(npwk, std::complex<double>(0.0, 0.0)));
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int ig = 0; ig < npwk; ++ig)
                {
                    vel[ib][ig] = 2.0 * tpiba2 * gk[ig][a] * gs_psi_(ik, ib, ig);
                }
            }
            // nonlocal derivative part (pos_matrix velocity form; NCPP
            // separable projectors only)
            for (int it = 0; it < ucell_->ntype; ++it)
            {
                const pseudo& ncpp = ucell_->atoms[it].ncpp;
                const int nh = ncpp.nh;
                if (nh == 0)
                {
                    continue;
                }
                // projector -> (radial beta index, m channel) table
                std::vector<int> mu_ib(nh, 0);
                std::vector<int> mu_m(nh, 0);
                int mu_idx = 0;
                for (int ib = 0; ib < ncpp.nbeta; ++ib)
                {
                    const int l = ncpp.lll[ib];
                    for (int m = 0; m < 2 * l + 1; ++m)
                    {
                        if (mu_idx < nh)
                        {
                            mu_ib[mu_idx] = ib;
                            mu_m[mu_idx] = m;
                        }
                        ++mu_idx;
                    }
                }
                for (int ia = 0; ia < ucell_->atoms[it].na; ++ia)
                {
                    std::vector<std::vector<std::complex<double>>> vkb;
                    pert_.build_vkb(it, ia, gk, vkb);
                    // becp_b[mu] = <vkb_mu|psi_b>
                    std::vector<std::vector<std::complex<double>>> becp(nbands);
                    for (int b = 0; b < nbands; ++b)
                    {
                        becp[b].assign(nh, std::complex<double>(0.0, 0.0));
                        for (int mu = 0; mu < nh; ++mu)
                        {
                            for (int ig = 0; ig < npwk; ++ig)
                            {
                                becp[b][mu] += std::conj(vkb[mu][ig]) * gs_psi_(ik, b, ig);
                            }
                        }
                    }
                    std::vector<std::vector<std::complex<double>>> dvkb;
                    pert_.build_vkb_dk(it, ia, a, gk, vkb, dvkb);
                    // dbecp_b[mu] = <dvkb_mu|psi_b>
                    std::vector<std::vector<std::complex<double>>> dbecp(nbands);
                    for (int b = 0; b < nbands; ++b)
                    {
                        dbecp[b].assign(nh, std::complex<double>(0.0, 0.0));
                        for (int mu = 0; mu < nh; ++mu)
                        {
                            for (int ig = 0; ig < npwk; ++ig)
                            {
                                dbecp[b][mu] += std::conj(dvkb[mu][ig]) * gs_psi_(ik, b, ig);
                            }
                        }
                    }
                    // dV_nl/dk_a|psi_b> = sum_mu |dvkb_mu> (D becp_b)_mu
                    //                              + |vkb_mu> (D dbecp_b)_mu
                    for (int b = 0; b < nbands; ++b)
                    {
                        for (int mu = 0; mu < nh; ++mu)
                        {
                            std::complex<double> out_b(0.0, 0.0);
                            std::complex<double> in_b(0.0, 0.0);
                            for (int nu = 0; nu < nh; ++nu)
                            {
                                if (mu_m[mu] != mu_m[nu])
                                {
                                    continue;
                                }
                                const double dij = ncpp.dion(mu_ib[mu], mu_ib[nu]);
                                out_b += dij * becp[b][nu];
                                in_b += dij * dbecp[b][nu];
                            }
                            for (int ig = 0; ig < npwk; ++ig)
                            {
                                vel[b][ig] += dvkb[mu][ig] * out_b + vkb[mu][ig] * in_b;
                            }
                        }
                    }
                }
            }
            // solve (H - eps_v) Y = -(i/tpiba) vel for every occupied band
            for (int ib = 0; ib < nbands; ++ib)
            {
                if (!dfpt_band_occupied(wg_, ik, ib))
                {
                    continue;
                }
                std::vector<std::complex<double>> rhs(npwk, std::complex<double>(0.0, 0.0));
                const std::complex<double> fac(0.0, -1.0 / tpiba);
                for (int ig = 0; ig < npwk; ++ig)
                {
                    rhs[ig] = fac * vel[ib][ig];
                }
                hamilt_->set_shift(eig_(ik, ib));
                double res = 0.0;
                stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr, yvec[ik][ib], res);
                if (dbg)
                {
                    std::cout << "DBG posresp a=" << a << " ik=" << ik << " ib=" << ib << " eps=" << eig_(ik, ib)
                              << " res=" << res << std::endl;
                }
            }
        }
        data_.set_pos_resp(a, yvec);
    }
    ModuleBase::timer::end("DFPT_PW", "solve_pos_resp");
}

void DFPT_PW::Impl::solve_efield_resp(int q_idx)
{
    ModuleBase::TITLE("DFPT_PW", "solve_efield_resp");
    ModuleBase::timer::start("DFPT_PW", "solve_efield_resp");
    // E-field SCF response (QE solve_e + dfpt_kernel form): the bare legs
    // Y^a stashed by solve_pos_resp are the field rhs base and the fixed
    // point adds the screened response potential of the mixed drho^E
    // exactly like solve_displacement. The converged dpsi^E,a feeds the
    // SCF dielectric tensor (DFPT_Q0::compute_eps) and the zstar_eu
    // cross-check probe (DFPT_ALEG).
    if (!wired() || hamilt_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_PW", "solve_efield_resp");
        return;
    }
    const ModuleBase::Vector3<double> q_cart = data_.get_qvec(q_idx) * ucell_->G;
    const int nrxx = pw_rho_->nrxx;
    const int nk = gs_psi_.get_nk();
    const int nbands = gs_psi_.get_nbands();
    const int lin_max = data_.get_max_iter();
    const double lin_thr = data_.get_conv_thr();

    for (int a = 0; a < 3; ++a)
    {
        const std::vector<std::vector<std::vector<std::complex<double>>>> yr = data_.get_pos_resp(a);
        if (static_cast<int>(yr.size()) != nk)
        {
            continue; // bare legs not solved: no E response either
        }
        rho_.reset_mixing(q_idx);
        data_.set_drho_g(q_idx, 0, std::vector<std::complex<double>>(pw_rho_->npw, std::complex<double>(0.0, 0.0)));
        bool converged = false;
        for (int iter = 0; iter < max_iter_ && !converged; ++iter)
        {
            // screened response potential of the mixed input density
            // (identical assembly to solve_displacement)
            std::vector<std::complex<double>> v_sc_r(nrxx, std::complex<double>(0.0, 0.0));
            const std::vector<std::complex<double>> drho_in_g = data_.get_drho_g(q_idx, 0);
            if (!drho_in_g.empty() && static_cast<int>(drho_in_g.size()) == pw_rho_->npw)
            {
                std::vector<std::complex<double>> dv_ha_g;
                rho_.v_hartree_q(q_cart, drho_in_g, dv_ha_g);
                pw_rho_->recip2real(dv_ha_g.data(), v_sc_r.data());
                if (xc_ != nullptr)
                {
                    std::vector<std::complex<double>> a_r(nrxx);
                    pw_rho_->recip2real(drho_in_g.data(), a_r.data());
                    std::vector<std::complex<double>> b_r;
                    xc_->apply(a_r, b_r);
                    if (static_cast<int>(b_r.size()) == nrxx)
                    {
                        for (int ir = 0; ir < nrxx; ++ir)
                        {
                            v_sc_r[ir] += b_r[ir];
                        }
                    }
                }
            }
            for (int ik = 0; ik < nk; ++ik)
            {
                if (static_cast<int>(occ_kq_.size()) <= ik || occ_kq_[ik].empty())
                {
                    continue;
                }
                std::vector<std::vector<std::complex<double>>> dv_sc;
                pert_.apply_vr(q_idx, ik, v_sc_r, gs_psi_, q_cart, dv_sc);
                if (last_q_ != q_idx || last_ik_ != ik)
                {
                    hamilt_->set_context(q_cart, ik);
                    last_q_ = q_idx;
                    last_ik_ = ik;
                }
                for (int ib = 0; ib < nbands; ++ib)
                {
                    if (!dfpt_band_occupied(wg_, ik, ib))
                    {
                        continue;
                    }
                    if (static_cast<int>(yr[ik][ib].size()) == 0 || static_cast<int>(dv_sc.size()) != nbands
                        || yr[ik][ib].size() != dv_sc[ib].size())
                    {
                        continue;
                    }
                    std::vector<std::complex<double>> rhs(yr[ik][ib].size());
                    for (size_t i = 0; i < rhs.size(); ++i)
                    {
                        rhs[i] = -(yr[ik][ib][i] + dv_sc[ib][i]);
                    }
                    hamilt_->set_shift(eig_(ik, ib));
                    std::vector<std::complex<double>> dpsi_out;
                    double res = 0.0;
                    stern_.solve(*hamilt_, occ_kq_[ik], rhs, lin_max, lin_thr, dpsi_out, res);
                    data_.set_dpsi(q_idx, ik, ib, dpsi_out);
                }
            }
            rho_.compute_drho(gs_psi_, wg_, q_idx, data_);
            rho_.mix_drho(q_idx, data_);
            const double residual = rho_.get_residual(q_idx, data_);
            converged = (residual < conv_thr_);
            if (converged)
            {
                std::cout << "DFPT efield dir=" << a << " converged, residual=" << residual << " (iter=" << iter << ")"
                          << std::endl;
            }
        }
        // stash dpsi^E,a before any later solve reuses the slots
        std::vector<std::vector<std::vector<std::complex<double>>>> de(
            nk,
            std::vector<std::vector<std::complex<double>>>(nbands));
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                de[ik][ib] = data_.get_dpsi(q_idx, ik, ib);
            }
        }
        data_.set_dpsi_efield(a, de);
    }
    ModuleBase::timer::end("DFPT_PW", "solve_efield_resp");
}

int DFPT_PW::get_nq() const
{
    return pimpl_->qlist_.get_nq();
}

ModuleBase::Vector3<double> DFPT_PW::get_qvec(int q_idx) const
{
    return pimpl_->data_.get_qvec(q_idx);
}

std::vector<double> DFPT_PW::get_phonon_freq(int q_idx) const
{
    return pimpl_->data_.get_phon_freq(q_idx);
}

std::vector<double> DFPT_PW::get_phon_freq_loto() const
{
    return pimpl_->data_.get_phon_freq_loto();
}

ModuleBase::Vector3<double> DFPT_PW::get_loto_dir() const
{
    return pimpl_->data_.get_loto_dir();
}

std::string DFPT_PW::format_q_report(int q_idx) const
{
    return pimpl_->phon_.format_q_report(q_idx, pimpl_->data_);
}

std::string DFPT_PW::format_loto_report() const
{
    return pimpl_->phon_.format_loto_report(pimpl_->data_);
}

ModuleBase::matrix DFPT_PW::get_dielectric_tensor() const
{
    return pimpl_->data_.get_dielectric();
}

ModuleBase::matrix DFPT_PW::get_born_charges(int atom_idx) const
{
    return pimpl_->data_.get_born(atom_idx);
}

void DFPT_PW::set_qfile(const std::string& filename)
{
    pimpl_->qfile_ = filename;
}

void DFPT_PW::set_qmesh(int nqx, int nqy, int nqz)
{
    pimpl_->nqx_ = nqx;
    pimpl_->nqy_ = nqy;
    pimpl_->nqz_ = nqz;
}

void DFPT_PW::set_conv_thr(double thr)
{
    pimpl_->conv_thr_ = thr;
    pimpl_->data_.set_conv_thr(thr);
}

void DFPT_PW::set_max_iter(int max_iter)
{
    pimpl_->max_iter_ = max_iter;
    pimpl_->data_.set_max_iter(max_iter);
}

void DFPT_PW::set_mix_beta(double beta)
{
    if (beta > 0.0 && beta <= 1.0)
    {
        pimpl_->mix_beta_ = beta;
    }
}

void DFPT_PW::set_compute_q0(bool flag)
{
    pimpl_->data_.set_compute_q0(flag);
}

void DFPT_PW::set_loto(bool flag)
{
    pimpl_->data_.set_loto(flag);
}

void DFPT_PW::set_loto_dir(const ModuleBase::Vector3<double>& dir)
{
    pimpl_->data_.set_loto_dir(dir);
}

} // namespace ModuleDFPT
