#include "to_wannier90_pw.h"
#include "source_base/parallel_comm.h" // use POOL_WORLD

#include "source_io/module_parameter/parameter.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_out/binstream.h"

toWannier90_PW::toWannier90_PW(
    const bool &out_wannier_mmn, 
    const bool &out_wannier_amn, 
    const bool &out_wannier_unk, 
    const bool &out_wannier_eig,
    const bool &out_wannier_wvfn_formatted, 
    const std::string &nnkpfile,
    const std::string &wannier_spin
):toWannier90(out_wannier_mmn, out_wannier_amn, out_wannier_unk, out_wannier_eig, out_wannier_wvfn_formatted, nnkpfile, wannier_spin)
{

}

toWannier90_PW::~toWannier90_PW()
{
    
}

void toWannier90_PW::calculate(
    const UnitCell& ucell,
    const ModuleBase::matrix& ekb,
    const ModulePW::PW_Basis_K* wfcpw,
    const ModulePW::PW_Basis_Big* bigpw,
    const K_Vectors& kv,
    const psi::Psi<std::complex<double>>* psi
)
{
    read_nnkp(ucell,kv);

    if (PARAM.inp.nspin == 2)
    {
        if (wannier_spin == "up")
        {
            start_k_index = 0;
        }
        else if (wannier_spin == "down")
        {
            start_k_index = num_kpts / 2;
        }
        else
        {
            ModuleBase::WARNING_QUIT("toWannier90::calculate", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    if (out_wannier_eig)
    {
        out_eig(ekb);
    }

    if (out_wannier_mmn)
    {
        cal_Mmn(*psi, wfcpw);
    }

    if (out_wannier_amn)
    {
        cal_Amn(*psi, wfcpw);
    }

    if (out_wannier_unk)
    {
        out_unk(*psi, wfcpw, bigpw);
    }

}
void toWannier90_PW::set_tpiba_omega(const double& tpiba, const double& omega)
{
    this->tpiba = &tpiba;
    this->omega = &omega;
    
}




void toWannier90_PW::unkdotkb(
    const psi::Psi<std::complex<double>>& psi_pw, 
    const ModulePW::PW_Basis_K* wfcpw,
    const int& cal_ik,
    const int& cal_ikb,
    const ModuleBase::Vector3<double> G,
    ModuleBase::ComplexMatrix &Mmn
)
{
    Mmn.create(num_bands, num_bands);

    for (int m = 0; m < num_bands; m++)
    {
        int im = cal_band_index[m];
        std::complex<double>* phase = new std::complex<double>[wfcpw->nmaxgr];
        ModuleBase::GlobalFunc::ZEROS(phase, wfcpw->nmaxgr);

        // get the phase value in realspace
        for (int ig = 0; ig < wfcpw->npwk[cal_ik]; ig++)
        {
            // if wfcpw->getgdirect(cal_ik, ig) == phase_G
            ModuleBase::Vector3<double> temp_G = wfcpw->getgdirect(cal_ik, ig) - G;
            if (temp_G.norm() < 1e-9)
            {
                phase[ig] = std::complex<double>(1.0, 0.0);
                break;
            }
        }

        wfcpw->recip2real(phase, phase, cal_ik);

        if (PARAM.inp.nspin == 4)
        {
            // (1) set value
            std::complex<double>* psir_up = new std::complex<double>[wfcpw->nmaxgr];
            std::complex<double>* psir_dn = new std::complex<double>[wfcpw->nmaxgr];
            ModuleBase::GlobalFunc::ZEROS(psir_up, wfcpw->nmaxgr);
            ModuleBase::GlobalFunc::ZEROS(psir_dn, wfcpw->nmaxgr);

            // (2) fft and get value
            // int npw_ik = wfcpw->npwk[cal_ik];
            int npwx = wfcpw->npwk_max;
            wfcpw->recip2real(&psi_pw(cal_ik, im, 0), psir_up, cal_ik);
            // wfcpw->recip2real(&psi_pw(cal_ik, im, npw_ik), psir_dn, cal_ik);
            wfcpw->recip2real(&psi_pw(cal_ik, im, npwx), psir_dn, cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir_up[ir] *= phase[ir];
                psir_dn[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir_up, psir_up, cal_ikb);
            wfcpw->real2recip(psir_dn, psir_dn, cal_ikb);

            for (int n = 0; n < num_bands; n++)
            {
                int in = cal_band_index[n];

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // int npw_ikb = wfcpw->npwk[cal_ikb];
                    int pwNumberMax = wfcpw->npwk_max;

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < pwNumberMax; ig++)
                    {
                        result_tem = result_tem + conj(psir_up[ig]) * psi_pw(cal_ikb, in, ig) + conj(psir_dn[ig]) * psi_pw(cal_ikb, in, ig+pwNumberMax);
                    }
#ifdef __MPI
                    Parallel_Reduce::reduce_all(result_tem);
#endif
                    Mmn(m, n) = result_tem;
                }
                else
                {
                    // GlobalV::ofs_running << "gamma only test" << std::endl;
                }

            }

            delete[] psir_up;
            delete[] psir_dn;

        }
        else
        {
            // (1) set value
            std::complex<double>* psir = new std::complex<double>[wfcpw->nmaxgr];
            ModuleBase::GlobalFunc::ZEROS(psir, wfcpw->nmaxgr);

            // (2) fft and get value
            wfcpw->recip2real(&psi_pw(cal_ik, im, 0), psir, cal_ik);
            for (int ir = 0; ir < wfcpw->nrxx; ir++)
            {
                psir[ir] *= phase[ir];
            }

            wfcpw->real2recip(psir, psir, cal_ikb);

            for (int n = 0; n < num_bands; n++)
            {
                int in = cal_band_index[n];

                if (!gamma_only_wannier)
                {
                    std::complex<double> result_tem(0.0, 0.0);

                    // Can be accelerated using lapack
                    for (int ig = 0; ig < wfcpw->npwk[cal_ikb]; ig++)
                    {
                        result_tem = result_tem + conj(psir[ig]) * psi_pw(cal_ikb, in, ig);
                    }

#ifdef __MPI
                    Parallel_Reduce::reduce_all(result_tem);
#endif
                    Mmn(m, n) = result_tem;
                }
                else
                {
                    // GlobalV::ofs_running << "gamma only test" << std::endl;
                }

            }

            delete[] psir;
            
        }

        delete[] phase;

    }

}

void toWannier90_PW::gen_radial_function_in_q(std::vector<ModuleBase::matrix> &radial_in_q)
{
    // The radial function is Fourier transformed into the q-space.
    radial_in_q.resize(num_wannier);

    double *r = new double[mesh_r];
    double *dr = new double[mesh_r];
    double *psi = new double[mesh_r];
    double *psir = new double[mesh_r];

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        double x = 0;
        for (int ir = 0; ir < mesh_r; ir++)
        {
            x = x_min + ir * dx;
            r[ir] = exp(x) / alfa[wannier_index];
            dr[ir] = dx * r[ir];
        }

        double alfa32 = pow(alfa[wannier_index], 3.0 / 2.0);
        double alfa_new = alfa[wannier_index];

        if (rvalue[wannier_index] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 2.0 * alfa32 * exp(-alfa_new * r[ir]);
            }
        }

        if (rvalue[wannier_index] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r[ir]) * exp(-alfa_new * r[ir] * 0.5);
            }
        }

        if (rvalue[wannier_index] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi[ir] = sqrt(4.0 / 27.0) * alfa32
                        * (1.0 - 2.0 / 3.0 * alfa_new * r[ir] + 2.0 / 27.0 * pow(alfa_new, 2.0) * r[ir] * r[ir])
                        * exp(-alfa_new * r[ir] * 1.0 / 3.0);
            }
        }

        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir[ir] = psi[ir] * r[ir];
        }

        auto &tmp_radial = radial_in_q[wannier_index];
        if (L[wannier_index] >= 0)
        {
            tmp_radial.create(1, PARAM.globalv.nqx);
            integral(mesh_r, psir, r, dr, L[wannier_index], tmp_radial.c);
        }
        else
        {
            int tmp_size = 0;

            if (L[wannier_index] == -1 || L[wannier_index] == -2 || L[wannier_index] == -3) tmp_size = 2;

            if (L[wannier_index] == -4 || L[wannier_index] == -5) tmp_size = 3;

            tmp_radial.create(tmp_size, PARAM.globalv.nqx);

            for (int tmp_L = 0; tmp_L < tmp_size; tmp_L++)
            {
                integral(mesh_r, psir, r, dr, tmp_L, tmp_radial.c+tmp_L*PARAM.globalv.nqx);
            }
        }

    }

    delete[] r;
    delete[] dr;
    delete[] psi;
    delete[] psir;

}

void toWannier90_PW::produce_trial_in_pw(
    const psi::Psi<std::complex<double>>& psi_pw,
    const int& ik,
    const ModulePW::PW_Basis_K* wfcpw,
    const std::vector<ModuleBase::matrix> &radial_in_q,
    ModuleBase::ComplexMatrix& trial_orbitals_k
)
{
    const int npw = wfcpw->npwk[ik];
    const int npwx = wfcpw->npwk_max;

    trial_orbitals_k.create(num_wannier, npwx);

    const int total_lm = 16;
    ModuleBase::matrix ylm(total_lm, npw);

    ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik, ig);
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk, ylm);

    // Keep it consistent with the definition of spherical harmonic functions in Wannier90 
    std::vector<int> need_inv = {2, 3, 5, 6, 14, 15};
    for (auto index : need_inv)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            ylm(index, ig) = -1.0 * ylm(index, ig);
        }
    }

    double bs2, bs3, bs6, bs12;
    bs2 = 1.0 / sqrt(2.0);
    bs3 = 1.0 / sqrt(3.0);
    bs6 = 1.0 / sqrt(6.0);
    bs12 = 1.0 / sqrt(12.0);

    std::complex<double> *sf = new std::complex<double>[npw];
    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        if (L[wannier_index] >= 0)
        {
            get_trial_orbitals_lm_k(L[wannier_index], m[wannier_index], ylm, gk, npw, radial_in_q[wannier_index].c, trial_orbitals_k.c + wannier_index*npwx);

            for (int ig = 0; ig < npw; ig++)
            {
                const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                trial_orbitals_k(wannier_index, ig) = sf[ig] * trial_orbitals_k(wannier_index, ig);
            }

        }
        else
        {
            if (L[wannier_index] == -1)
            {
                double tmp_bs2 = 0;
                if (m[wannier_index] == 0) tmp_bs2 = bs2;
                if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                std::complex<double> *orb_s = new std::complex<double>[npw];
                std::complex<double> *orb_px = new std::complex<double>[npw];

                get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_px);

                for (int ig = 0; ig < npw; ig++)
                {
                    const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                    sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                    trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs2 * orb_s[ig]  + tmp_bs2 * orb_px[ig]);
                }

                delete[] orb_s;
                delete[] orb_px;
            }

            if (L[wannier_index] == -2)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = bs2;
                    if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    std::complex<double> *orb_py = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_px);
                    get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_py);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] - bs6 * orb_px[ig] + tmp_bs2 * orb_py[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                    delete[] orb_py;
                }
                else if (m[wannier_index] == 2)
                {
                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_px);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] + 2.0 * bs6 * orb_px[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                }
            }

            if (L[wannier_index] == -3)
            {
                double m_px = 1.0;
                double m_py = 1.0;
                double m_pz = 1.0;

                if (m[wannier_index] == 1)
                {
                    m_py = -1.0;
                    m_pz = -1.0;
                }
                else if (m[wannier_index] == 2)
                {
                    m_px = -1.0;
                    m_pz = -1.0;
                }
                else if (m[wannier_index] == 3)
                {
                    m_px = -1.0;
                    m_py = -1.0;
                }

                std::complex<double> *orb_s = new std::complex<double>[npw];
                std::complex<double> *orb_px = new std::complex<double>[npw];
                std::complex<double> *orb_py = new std::complex<double>[npw];
                std::complex<double> *orb_pz = new std::complex<double>[npw];

                get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_px);
                get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_py);
                get_trial_orbitals_lm_k(1, 0, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_pz);

                for (int ig = 0; ig < npw; ig++)
                {
                    const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                    sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                    trial_orbitals_k(wannier_index, ig) = sf[ig] * 0.5 * (orb_s[ig] + m_px * orb_px[ig] + m_py * orb_py[ig] + m_pz * orb_pz[ig]);
                }

                delete[] orb_s;
                delete[] orb_px;
                delete[] orb_py;
                delete[] orb_pz;
            }
            
            if (L[wannier_index] == -4)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = bs2;
                    if (m[wannier_index] == 1) tmp_bs2 = -bs2;

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    std::complex<double> *orb_py = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_px);
                    get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_py);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] - bs6 * orb_px[ig] + tmp_bs2 * orb_py[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                    delete[] orb_py;
                }
                else if (m[wannier_index] == 2)
                {
                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_px);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs3 * orb_s[ig] + 2.0 * bs6 * orb_px[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px; 
                }
                else if (m[wannier_index] == 3 || m[wannier_index] == 4)
                {
                    double m_pz = 1.0;
                    if (m[wannier_index] == 4) m_pz = -1.0;

                    std::complex<double> *orb_pz = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(1, 0, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_pz);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*PARAM.globalv.nqx, orb_dz2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * bs2 * (m_pz * orb_pz[ig] + orb_dz2[ig]);
                    }

                    delete[] orb_pz;
                    delete[] orb_dz2;
                }
            }

            if (L[wannier_index] == -5)
            {
                if (m[wannier_index] == 0 || m[wannier_index] == 1)
                {
                    double tmp_bs2 = -bs2;
                    double tmp_bs12 = -bs12;
                    double tmp_d = 0.5;

                    if (m[wannier_index] == 1)
                    {
                        tmp_bs2 = bs2;
                    }

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_px = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];
                    std::complex<double> *orb_dx2_y2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 1, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_px);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*PARAM.globalv.nqx, orb_dz2);
                    get_trial_orbitals_lm_k(2, 3, ylm, gk, npw, radial_in_q[wannier_index].c+2*PARAM.globalv.nqx, orb_dx2_y2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs6 * orb_s[ig] + tmp_bs2 * orb_px[ig] + tmp_bs12 * orb_dz2[ig] + tmp_d * orb_dx2_y2[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_px;
                    delete[] orb_dz2;
                    delete[] orb_dx2_y2;
                }
                else if (m[wannier_index] == 2 || m[wannier_index] == 3)
                {
                    double tmp_bs2 = -bs2;
                    double tmp_bs12 = -bs12;
                    double tmp_d = -0.5;

                    if (m[wannier_index] == 3)
                    {
                        tmp_bs2 = bs2;
                    }

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_py = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];
                    std::complex<double> *orb_dx2_y2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 2, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_py);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*PARAM.globalv.nqx, orb_dz2);
                    get_trial_orbitals_lm_k(2, 3, ylm, gk, npw, radial_in_q[wannier_index].c+2*PARAM.globalv.nqx, orb_dx2_y2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs6 * orb_s[ig] + tmp_bs2 * orb_py[ig] + tmp_bs12 * orb_dz2[ig] + tmp_d * orb_dx2_y2[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_py;
                    delete[] orb_dz2;
                    delete[] orb_dx2_y2;
                }
                else if (m[wannier_index] == 4 || m[wannier_index] == 5)
                {
                    double tmp_pz = -1.0;

                    if (m[wannier_index] == 5) tmp_pz = 1.0;

                    std::complex<double> *orb_s = new std::complex<double>[npw];
                    std::complex<double> *orb_pz = new std::complex<double>[npw];
                    std::complex<double> *orb_dz2 = new std::complex<double>[npw];

                    get_trial_orbitals_lm_k(0, 0, ylm, gk, npw, radial_in_q[wannier_index].c, orb_s);
                    get_trial_orbitals_lm_k(1, 0, ylm, gk, npw, radial_in_q[wannier_index].c+PARAM.globalv.nqx, orb_pz);
                    get_trial_orbitals_lm_k(2, 0, ylm, gk, npw, radial_in_q[wannier_index].c+2*PARAM.globalv.nqx, orb_dz2);

                    for (int ig = 0; ig < npw; ig++)
                    {
                        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
                        sf[ig] = std::complex<double>(cos(arg), -sin(arg));

                        trial_orbitals_k(wannier_index, ig) = sf[ig] * (bs6 * orb_s[ig] + tmp_pz * bs2 * orb_pz[ig] + bs3 * orb_dz2[ig]);
                    }

                    delete[] orb_s;
                    delete[] orb_pz;
                    delete[] orb_dz2;
                }
            }

        }
    }

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        std::complex<double> anorm(0.0, 0.0);
        for (int ig = 0; ig < npw; ig++)
        {
            anorm += conj(trial_orbitals_k(wannier_index, ig)) * trial_orbitals_k(wannier_index, ig);
        }

#ifdef __MPI
        Parallel_Reduce::reduce_all(anorm);
#endif
        for (int ig = 0; ig < npw; ig++)
        {
            trial_orbitals_k(wannier_index, ig) = trial_orbitals_k(wannier_index, ig) / sqrt(anorm);
        }
    }

    delete[] gk;
    delete[] sf;
}

void toWannier90_PW::get_trial_orbitals_lm_k(
    const int &orbital_L,
    const int &orbital_m,
    const ModuleBase::matrix &ylm,
    const ModuleBase::Vector3<double> *gk,
    const int &npw,
    double *radial_in_q_single,
    std::complex<double> *orbital_in_G_single
)
{
    const double tpiba = *this->tpiba;
    for (int ig = 0; ig < npw; ig++)
    {
        orbital_in_G_single[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(radial_in_q_single, PARAM.globalv.nqx, PARAM.globalv.dq, gk[ig].norm() * tpiba);
    }

    std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, orbital_L);
    int index = orbital_L * orbital_L + orbital_m;
    for (int ig = 0; ig < npw; ig++)
    {
        orbital_in_G_single[ig] = lphase * ylm(index, ig) * orbital_in_G_single[ig];
    }

    return;
}

void toWannier90_PW::integral(
    const int meshr,
    const double *psir,
    const double *r,
    const double *rab,
    const int &l,
    double *table
)
{
    const double pref = ModuleBase::FOUR_PI / sqrt(*this->omega);

    double *inner_part = new double[meshr];
    for (int ir = 0; ir < meshr; ir++)
    {
        inner_part[ir] = psir[ir] * psir[ir];
    }

    double unit = 0.0;
    ModuleBase::Integral::Simpson_Integral(meshr, inner_part, rab, unit);
    delete[] inner_part;

    double *aux = new double[meshr];
    double *vchi = new double[meshr];
    for (int iq = 0; iq < PARAM.globalv.nqx; iq++)
    {
        const double q = PARAM.globalv.dq * iq;
        ModuleBase::Sphbes::Spherical_Bessel(meshr, r, q, l, aux);
        for (int ir = 0; ir < meshr; ir++)
        {
            vchi[ir] = psir[ir] * aux[ir] * r[ir];
        }

        double vqint = 0.0;
        ModuleBase::Integral::Simpson_Integral(meshr, vchi, rab, vqint);

        table[iq] = vqint * pref;
    }
    delete[] aux;
    delete[] vchi;
    return;
}

void toWannier90_PW::unkdotW_A(
    const psi::Psi<std::complex<double>>& psi_pw,
    const ModulePW::PW_Basis_K* wfcpw,
    const int& ik,
    const std::vector<ModuleBase::matrix> &radial_in_q,
    ModuleBase::ComplexMatrix &Amn
)
{
    Amn.create(num_bands, num_wannier);

    int npw = wfcpw->npwk[ik];
    int npwx = wfcpw->npwk_max;
    ModuleBase::ComplexMatrix trial_orbitals;
    produce_trial_in_pw(psi_pw, ik, wfcpw, radial_in_q, trial_orbitals);

    for (int iw = 0; iw < num_wannier; iw++)
    {
        for (int ib_w = 0; ib_w < num_bands; ib_w++)
        {
            int ib = cal_band_index[ib_w];

            if (PARAM.inp.nspin != 4)
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    Amn(ib_w, iw) += conj(psi_pw(ik, ib, ig)) * trial_orbitals(iw, ig);
                }
            }
            else
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    Amn(ib_w, iw) += up_con[iw] * conj(psi_pw(ik, ib, ig)) * trial_orbitals(iw, ig)
                                   + dn_con[iw] * conj(psi_pw(ik, ib, ig+npwx)) * trial_orbitals(iw, ig);
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_all(Amn.c, Amn.size);
#endif

}
