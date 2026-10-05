#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// local helpers for eigenvalue calculation
// migrated from dftu_base.cpp, mohan 2025-11-08
namespace
{

inline void JacobiRotate(std::vector<std::vector<double>>& A, int p, int q, int n)
{
    if (std::abs(A[p][q]) > 1e-10)
    {
        double r = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t = 0.0;
        if (r >= 0)
        {
            t = 1.0 / (r + sqrt(1.0 + r * r));
        }
        else
        {
            t = -1.0 / (-r + sqrt(1.0 + r * r));
        }
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;

        A[p][p] -= t * A[p][q];
        A[q][q] += t * A[p][q];
        A[p][q] = A[q][p] = 0.0;

        for (int k = 0; k < n; k++)
        {
            if (k != p && k != q)
            {
                double Akp = c * A[k][p] - s * A[k][q];
                double Akq = s * A[k][p] + c * A[k][q];
                A[k][p] = A[p][k] = Akp;
                A[k][q] = A[q][k] = Akq;
            }
        }
    }
}

inline std::vector<double> CalculateEigenvalues(std::vector<std::vector<double>>& A, int n)
{
    std::vector<double> eigenvalues(n);
    while (true)
    {
        int p = 0, q = 1;
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (std::abs(A[i][j]) > std::abs(A[p][q]))
                {
                    p = i;
                    q = j;
                }
            }
        }

        if (std::abs(A[p][q]) < 1e-10)
        {
            for (int i = 0; i < n; i++)
            {
                eigenvalues[i] = A[i][i];
            }
            break;
        }

        JacobiRotate(A, p, q, n);
    }
    return eigenvalues;
}

/// @brief Extract the four Pauli-component blocks of one SOC correlated shell.
///
/// The blocks, each (2l+1)x(2l+1) stored row-major with index m0*m+m1, are
/// the charge channel and the sigma_x/y/z spin channels. The PW storage keeps
/// the four blocks contiguous; the LCAO storage is a real symmetric matrix in
/// spin basis, whose spin-off-diagonal imaginary channel is not stored.
///
/// @param dftu DFT+U object holding the occupation matrix
/// @param iat atom index
/// @param l angular-momentum channel
/// @param layout storage layout of the nspin == 4 occupation matrix
/// @param blocks [out] four Pauli-component blocks
inline void extract_soc_pauli_blocks(const Plus_U_Base& dftu,
                                     const int iat,
                                     const int l,
                                     const OccmatSocLayout layout,
                                     std::vector<std::vector<double>>& blocks)
{
    const int m = 2 * l + 1;
    const int m2 = m * m;
    blocks.assign(4, std::vector<double>(m2, 0.0));

    const OccupationMatrix& occmat = dftu.occmat();
    if (layout == SOC_LAYOUT_PAULI)
    {
        const ModuleBase::matrix& occ = occmat.mat(iat, l, 0);
        for (int is = 0; is < 4; ++is)
        {
            for (int k = 0; k < m2; ++k)
            {
                blocks[is][k] = occ.c[is * m2 + k];
            }
        }
    }
    else // SOC_LAYOUT_SPIN_BASIS_REAL
    {
        for (int m0 = 0; m0 < m; ++m0)
        {
            for (int m1 = 0; m1 < m; ++m1)
            {
                const int k = m0 * m + m1;
                const double suu = occmat.get(iat, l, 0, m0, m1);
                const double sud = occmat.get(iat, l, 0, m0, m1 + m);
                const double sdu = occmat.get(iat, l, 0, m0 + m, m1);
                const double sdd = occmat.get(iat, l, 0, m0 + m, m1 + m);
                blocks[0][k] = suu + sdd;
                blocks[1][k] = sud + sdu;
                // The spin-basis storage drops imaginary parts, so the
                // sigma_y channel has no contribution.
                blocks[2][k] = 0.0;
                blocks[3][k] = suu - sdd;
            }
        }
    }
}

/// @brief Diagonalize the SOC occupation matrix reconstructed in spinor space.
///
/// The Hermitian spinor matrix is built from the Pauli blocks as
///   n_upup   = (b0 + b3)/2, n_dndn = (b0 - b3)/2,
///   n_updown = (b1 + i*b2)/2.
/// Its 2*(2l+1) real eigenvalues are computed with LAPACK zheev.
///
/// @param blocks four Pauli-component blocks
/// @param m 2l+1
/// @return eigenvalues sorted ascending by zheev
inline std::vector<double> calculate_spinor_eigenvalues(
    const std::vector<std::vector<double>>& blocks,
    const int m)
{
    const int n = 2 * m;
    ModuleBase::ComplexMatrix spinor_n(n, n);
    for (int m0 = 0; m0 < m; ++m0)
    {
        for (int m1 = 0; m1 < m; ++m1)
        {
            const int k = m0 * m + m1;
            const double b0 = blocks[0][k];
            const double b1 = blocks[1][k];
            const double b2 = blocks[2][k];
            const double b3 = blocks[3][k];

            spinor_n(m0, m1) = std::complex<double>(0.5 * (b0 + b3), 0.0);
            spinor_n(m0, m1 + m) = std::complex<double>(0.5 * b1, 0.5 * b2);
            spinor_n(m0 + m, m1) = std::complex<double>(0.5 * b1, -0.5 * b2);
            spinor_n(m0 + m, m1 + m) = std::complex<double>(0.5 * (b0 - b3), 0.0);
        }
    }

    std::vector<double> eigenvalues(n, 0.0);
    std::vector<std::complex<double>> work(1, std::complex<double>(0.0, 0.0));
    std::vector<double> rwork(std::max(1, 3 * n - 2), 0.0);
    int info = 0;
    int lwork = -1;

    LapackConnector::zheev('N', 'U', n, spinor_n, n, eigenvalues.data(),
                           work.data(), lwork, rwork.data(), &info);
    if (info != 0)
    {
        ModuleBase::WARNING_QUIT("calculate_spinor_eigenvalues",
                                 "LAPACK zheev workspace query failed");
    }
    lwork = std::max(1, static_cast<int>(work[0].real()));
    work.resize(lwork);
    LapackConnector::zheev('N', 'U', n, spinor_n, n, eigenvalues.data(),
                           work.data(), lwork, rwork.data(), &info);
    if (info != 0)
    {
        ModuleBase::WARNING_QUIT("calculate_spinor_eigenvalues",
                                 "LAPACK zheev diagonalization failed");
    }
    return eigenvalues;
}

} // namespace


namespace DFTU_BASE
{

bool is_ion_step_output_step(const int istep, const OccmatOutputCfg& cfg)
{
    if (cfg.out_freq_ion <= 0)
    {
        return false;
    }
    return (istep % cfg.out_freq_ion == 0);
}

bool is_elec_snapshot_trigger(const int iter,
                              const bool conv_esolver,
                              const OccmatOutputCfg& cfg)
{
    const bool periodic = (iter % cfg.out_freq_elec == 0);
    const bool last_step = (iter == cfg.scf_nmax);
    return periodic || last_step || conv_esolver;
}

std::string gen_ion_step_dm_onsite_filename(const std::string& out_dir, const int istep)
{
    std::stringstream ss;
    ss << out_dir << "dm_onsiteg" << (istep + 1) << ".txt";
    return ss.str();
}

void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& l_channel,
                  const int occ_mat_ctrl,
                  const std::string& fn,
                  const std::string& init_chg,
                  int nspin,
                  int npol)
{
    ModuleBase::TITLE("DFTU_BASE", "read_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ifstream ifdftu(fn.c_str(), std::ios::in);

    if (!ifdftu)
    {
        if (occ_mat_ctrl > 0)
        {
            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "Can not find the file dm_onsite_ini.txt. Please check your dm_onsite_ini.txt");
        }
        else
        {
            if (init_chg == "file")
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "Can not find the file dm_onsite.txt. Please do scf calculation first");
            }
        }
        ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "Can not open dm_onsite.txt file");
    }

    ifdftu.clear();
    ifdftu.seekg(0);

    char word[20];

    int T = 0;
    int iat = 0;
    int spin = 0;
    int L = 0;
    int zeta = 0;

    ifdftu.rdstate();

    while (ifdftu.good())
    {
        ifdftu >> word;
        if (ifdftu.eof())
        {
            break;
        }

        if (strcmp("Atom=", word) == 0)
        {
            ifdftu >> iat;
            iat -= 1;
            ifdftu >> word;

            if (strcmp("L=", word) != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> L;
            ifdftu >> word;

            if (strcmp("ORBITAL=", word) != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> zeta;
            ifdftu.ignore(150, '\n');

            if (zeta != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m",
                                         "only the first radial channel (ORBITAL=0) is supported");
            }

            T = ucell.iat2it[iat];
            const int NL = ucell.atoms[T].nwl + 1;

            for (int l = 0; l < NL; l++)
            {
                if (l != l_channel[T])
                {
                    continue;
                }

                if (nspin == 1 || nspin == 2)
                {
                    for (int is = 0; is < 2; is++)
                    {
                        ifdftu >> word;
                        if (strcmp("spin=", word) == 0)
                        {
                            ifdftu >> spin;
                            spin -= 1;
                            ifdftu.ignore(150, '\n');

                            double value = 0.0;
                            for (int m0 = 0; m0 < 2 * L + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * L + 1; m1++)
                                {
                                    ifdftu >> value;
                                    occ.set(iat, L, spin, m0, m1, value);
                                }
                                ifdftu.ignore(150, '\n');
                            }
                        }
                        else
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                    }
                }
                else if (nspin == 4) // SOC
                {
                    double value = 0.0;
                    for (int m0 = 0; m0 < 2 * L + 1; m0++)
                    {
                        for (int ipol0 = 0; ipol0 < npol; ipol0++)
                        {
                            const int m0_all = m0 + (2 * L + 1) * ipol0;

                            for (int m1 = 0; m1 < 2 * L + 1; m1++)
                            {
                                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                {
                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    ifdftu >> value;
                                    occ.set(iat, L, 0, m0_all, m1_all, value);
                                }
                            }
                            ifdftu.ignore(150, '\n');
                        }
                    }
                }
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
        }

        ifdftu.rdstate();

        if (ifdftu.eof() != 0)
        {
            break;
        }
    }

    return;
}

#ifdef __MPI
/// Broadcast the local occupation number matrices from rank 0 to all ranks.
///
/// Each occupation matrix is broadcast as one contiguous block
/// (matrix::c stores nr * nc consecutive doubles) instead of element by
/// element.
void local_occup_bcast(const UnitCell& ucell,
                       OccupationMatrix& occ,
                       const std::vector<int>& l_channel,
                       int nspin,
                       int npol)
{
    ModuleBase::TITLE("DFTU_BASE", "local_occup_bcast");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (l_channel[T] == -1)
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            const int L = l_channel[T];

            for (int l = 0; l <= ucell.atoms[T].nwl; l++)
            {
                if (l != l_channel[T])
                {
                    continue;
                }

                if (nspin == 1 || nspin == 2)
                {
                    for (int spin = 0; spin < 2; spin++)
                    {
                        Parallel_Common::bcast_double(occ.mat(iat, l, spin).c,
                                                      occ.mat(iat, l, spin).nr * occ.mat(iat, l, spin).nc);
                    }
                }
                else if (nspin == 4) // SOC
                {
                    Parallel_Common::bcast_double(occ.mat(iat, l, 0).c,
                                                  occ.mat(iat, l, 0).nr * occ.mat(iat, l, 0).nc);
                }
            }
        }
    }
    return;
}
#endif


void output(const Plus_U_Base& dftu,
            const UnitCell& ucell,
            bool out_chg,
            const std::string& global_out_dir,
            int nspin,
            int npol,
            int istep,
            int iter,
            const OccmatOutputCfg& cfg,
            OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "output");

    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output", "istep must be >= 0");
    }
    if (iter < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output", "iter must be >= 1");
    }
    if (cfg.out_freq_ion < 0 || cfg.out_freq_elec < 1 || cfg.scf_nmax < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::output",
                                "invalid occupation-matrix output frequency configuration");
    }

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | #DFT+U INFORMATION# |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;

        for (int L = 0; L < NL; L++)
        {
            const int N = ucell.atoms[T].l_nchi[L];

            if (L >= dftu.get_l_channel(T) && dftu.has_l_channel(T))
            {
                if (L != dftu.get_l_channel(T))
                {
                    continue;
                }

                if (!dftu.use_yukawa())
                {
                    GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << " ORBITAL=" << 0
                                         << " U=" << dftu.get_u_current(T) * ModuleBase::Ry_to_eV << " eV" << std::endl;
                }
                else
                {
                    double Ueff = (dftu.yukawa().get_U(T, L) - dftu.yukawa().get_J(T, L)) * ModuleBase::Ry_to_eV;
                    GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << "  ORBITAL=" << 0
                                         << " U=" << dftu.yukawa().get_U(T, L) * ModuleBase::Ry_to_eV << " eV"
                                         << " J=" << dftu.yukawa().get_J(T, L) * ModuleBase::Ry_to_eV << " eV"
                                         << std::endl;
                }
            }
        }
    }

    GlobalV::ofs_running << " Local Occupation Matrices for each atom" << std::endl;
    write_occup_m(dftu, ucell, GlobalV::ofs_running, true, nspin, npol,
                  OCMAT_FMT_LEGACY, soc_layout);

    // dm_onsite.txt is always overwritten with the latest occupation matrix;
    // it is the entry file of init_chg=file and NSCF restarts.
    if (out_chg && GlobalV::MY_RANK == 0)
    {
        const std::string latest_fn = global_out_dir + "dm_onsite.txt";
        std::ofstream ofdftu;
        ofdftu.open(latest_fn);
        if (!ofdftu)
        {
            ModuleBase::WARNING_QUIT("DFTU_BASE::output", "Can't create file dm_onsite.txt");
        }
        write_occup_m(dftu, ucell, ofdftu, false, nspin, npol,
                      OCMAT_FMT_LEGACY, soc_layout);
        ofdftu.close();
    }

    // At the first electronic step of an output ionic step, create the
    // per-ionic-step file and write its header. Electronic-step sections are
    // appended later by append_ion_step_snapshot().
    const bool ion_step_output = is_ion_step_output_step(istep, cfg);
    if (out_chg && ion_step_output && iter == 1 && GlobalV::MY_RANK == 0)
    {
        const std::string ion_step_fn = gen_ion_step_dm_onsite_filename(global_out_dir, istep);
        std::ofstream ofs_ion_step;
        ofs_ion_step.open(ion_step_fn, std::ios::out | std::ios::trunc);
        if (!ofs_ion_step)
        {
            ModuleBase::WARNING_QUIT("DFTU_BASE::output",
                                     "Can't create per-ionic-step occupation-matrix file");
        }
        ofs_ion_step << "# ################################################################\n";
        ofs_ion_step << "# # DFT+U occupation-matrix snapshots\n";
        ofs_ion_step << "# # Ionic (geometry) step g" << (istep + 1) << "\n";
        ofs_ion_step << "# # Electronic steps recorded every out_freq_elec = "
                     << cfg.out_freq_elec << " iterations\n";
        ofs_ion_step << "# ################################################################\n";
        ofs_ion_step.close();
    }

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | # END DFT+U INFO    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl << std::endl;

    return;
}

void append_ion_step_snapshot(const Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const std::string& global_out_dir,
                              int nspin,
                              int npol,
                              int istep,
                              int iter,
                              bool conv_esolver,
                              double etot_ry,
                              double tot_mag,
                              const double* tot_mag_nc,
                              const OccmatOutputCfg& cfg,
                              OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "append_ion_step_snapshot");

    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "nspin must be 1, 2 or 4");
    }
    if (istep < 0)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "istep must be >= 0");
    }
    if (iter < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot", "iter must be >= 1");
    }
    if (nspin == 4 && tot_mag_nc == nullptr)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                 "tot_mag_nc must be provided for non-collinear calculations");
    }
    // out_freq_ion == 0 is the valid default (no numbered file) and is
    // handled by the gates below; only negative values are an error.
    if (cfg.out_freq_ion < 0 || cfg.out_freq_elec < 1 || cfg.scf_nmax < 1)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                "invalid occupation-matrix output frequency configuration");
    }

    const bool ion_step_output = is_ion_step_output_step(istep, cfg);
    const bool elec_trigger = is_elec_snapshot_trigger(iter, conv_esolver, cfg);
    if (!ion_step_output || !elec_trigger || GlobalV::MY_RANK != 0)
    {
        return;
    }

    // Human-readable convergence status of this electronic step.
    std::string status = "not_converged";
    if (conv_esolver)
    {
        status = "converged";
    }
    else if (iter == cfg.scf_nmax)
    {
        status = "reached_scf_nmax";
    }

    const std::string ion_step_fn = gen_ion_step_dm_onsite_filename(global_out_dir, istep);
    std::ofstream ofs_ion_step;
    ofs_ion_step.open(ion_step_fn, std::ios::app);
    if (!ofs_ion_step)
    {
        ModuleBase::WARNING_QUIT("DFTU_BASE::append_ion_step_snapshot",
                                 "Can't open per-ionic-step occupation-matrix file");
    }

    // In append mode the put pointer sits at end of file; position 0 means
    // the file was just created and the section must not start with a blank line.
    const bool file_was_empty = (ofs_ion_step.tellp() == std::streampos(0));

    const double etot_ev = etot_ry * ModuleBase::Ry_to_eV;

    // Separate consecutive electronic-step sections by a blank line.
    if (!file_was_empty)
    {
        ofs_ion_step << "\n";
    }
    ofs_ion_step << "# ================================================================\n";
    ofs_ion_step << "# Electronic step " << iter << "\n";
    ofs_ion_step << "# Total energy (Kohn-Sham): " << std::setw(20)
                 << std::setprecision(8) << std::fixed << etot_ev << " eV\n";
    if (nspin == 4)
    {
        ofs_ion_step << "# Total magnetism (Bohr mag/cell) mx, my, mz:"
                     << std::setw(14) << tot_mag_nc[0]
                     << std::setw(14) << tot_mag_nc[1]
                     << std::setw(14) << tot_mag_nc[2] << "\n";
    }
    else
    {
        ofs_ion_step << "# Total magnetism (Bohr mag/cell):"
                     << std::setw(14) << tot_mag << "\n";
    }
    ofs_ion_step << "# Status: " << status << "\n";
    ofs_ion_step << "# ================================================================\n";

    // Full occupation matrices together with eigenvalues and per-atom magnetism.
    write_occup_m(dftu, ucell, ofs_ion_step, true, nspin, npol,
                  OCMAT_FMT_READABLE, soc_layout);

    ofs_ion_step.close();

    GlobalV::ofs_running << " Append local occupation matrices of electronic step " << iter
                         << " to file: " << ion_step_fn << std::endl;
    return;
}


void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ofstream& ofs,
                   bool diag,
                   int nspin,
                   int npol,
                   OccmatTextFormat fmt,
                   OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "write_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!dftu.has_l_channel(T))
        {
            continue;
        }
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = dftu.get_l_channel(T);

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < NL; l++)
            {
                if (l != dftu.get_l_channel(T))
                {
                    continue;
                }

                if (fmt == OCMAT_FMT_READABLE)
                {
                    ofs << "\nAtom=" << iat + 1;
                    ofs << " L=" << l << std::endl;
                }
                else
                {
                    ofs << "\n Atom= " << iat + 1;
                    ofs << " L= " << l;
                    ofs << " ORBITAL= " << 0 << std::endl;
                }

                if (nspin == 1 || nspin == 2)
                {
                    double sum0[2];
                    const std::string eigen_label = (fmt == OCMAT_FMT_READABLE)
                        ? " Eigenvalues for spin " : " Eigenvalues for spin=";
                    const std::string trace_label = (fmt == OCMAT_FMT_READABLE)
                        ? " Trace (electrons) = " : " sum is ";
                    const std::string matrix_label = (fmt == OCMAT_FMT_READABLE)
                        ? " Occupation matrix for spin " : " spin= ";
                    for (int is = 0; is < 2; is++)
                    {
                        if (diag)
                        {
                            std::vector<std::vector<double>> A(2 * l + 1, std::vector<double>(2 * l + 1));
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    A[m0][m1] = dftu.occmat().get(iat, l, is, m0, m1);
                                }
                            }
                            std::vector<double> eigenvalues = CalculateEigenvalues(A, 2 * l + 1);
                            sum0[is] = 0.0;
                            ofs << eigen_label << is + 1 << std::endl;
                            ofs << std::setprecision(8) << std::fixed;
                            for (int i = 0; i < 2 * l + 1; i++)
                            {
                                ofs << std::setw(12) << eigenvalues[i];
                                sum0[is] += eigenvalues[i];
                            }
                            ofs << std::endl;
                            ofs << trace_label << std::setw(12) << sum0[is] << std::endl;
                        }
                        ofs << matrix_label << is + 1 << std::endl;
                        ofs << std::setprecision(8) << std::fixed;
                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int m1 = 0; m1 < 2 * l + 1; m1++)
                            {
                                ofs << std::setw(12)
                                    << dftu.occmat().get(iat, l, is, m0, m1);
                            }
                            ofs << std::endl;
                        }
                    }
                    if (diag)
                    {
                        ofs << std::setw(12) << std::setprecision(8)
                            << std::fixed << " Magnetism for atom " << iat+1 << ": " << sum0[0] - sum0[1]
                            << std::endl;
                    }
                }
                else if (nspin == 4) // SOC
                {
                    const int m = 2 * l + 1;
                    std::vector<std::vector<double>> blocks;
                    extract_soc_pauli_blocks(dftu, iat, l, soc_layout, blocks);

                    if (diag)
                    {
                        // Occupation numbers of the correlated spinor shell.
                        std::vector<double> eigenvalues = calculate_spinor_eigenvalues(blocks, m);
                        ofs << " Eigenvalues in spinor space:" << std::endl;
                        ofs << std::setprecision(8) << std::fixed;
                        for (int i = 0; i < 2 * m; ++i)
                        {
                            ofs << std::setw(12) << eigenvalues[i];
                            if ((i + 1) % m == 0)
                            {
                                ofs << std::endl;
                            }
                        }

                        // Per-atom charge and spin-moment traces of the shell.
                        double n_tot = 0.0;
                        double mag_x = 0.0;
                        double mag_y = 0.0;
                        double mag_z = 0.0;
                        for (int mm = 0; mm < m; ++mm)
                        {
                            const int k = mm * m + mm;
                            n_tot += blocks[0][k];
                            mag_x += blocks[1][k];
                            mag_y += blocks[2][k];
                            mag_z += blocks[3][k];
                        }
                        ofs << " Trace (electrons) = " << std::setw(12) << n_tot << std::endl;
                        ofs << " Magnetism for atom " << iat + 1 << " (mx, my, mz):"
                            << std::setw(12) << mag_x << std::setw(12) << mag_y
                            << std::setw(12) << mag_z << std::endl;

                        // Hermitian occupation matrix shown spin block by spin block.
                        ofs << " Occupation matrix, spin up/up (Re):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * (blocks[0][k] + blocks[3][k]);
                                ofs << std::setw(12) << val;
                            }
                            ofs << std::endl;
                        }
                        ofs << " Occupation matrix, spin up/down (Re):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * blocks[1][k];
                                ofs << std::setw(12) << val;
                            }
                            ofs << std::endl;
                        }
                        ofs << " Occupation matrix, spin up/down (Im):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * blocks[2][k];
                                ofs << std::setw(12) << val;
                            }
                            ofs << std::endl;
                        }
                        ofs << " Occupation matrix, spin down/down (Re):" << std::endl;
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int m1 = 0; m1 < m; ++m1)
                            {
                                const int k = m0 * m + m1;
                                const double val = 0.5 * (blocks[0][k] - blocks[3][k]);
                                ofs << std::setw(12) << val;
                            }
                            ofs << std::endl;
                        }
                    }
                    else
                    {
                        // Real 2m x 2m matrix in spin basis; this is the
                        // layout parsed by read_occup_m() on restart.
                        for (int m0 = 0; m0 < m; ++m0)
                        {
                            for (int ipol0 = 0; ipol0 < npol; ++ipol0)
                            {
                                for (int m1 = 0; m1 < m; ++m1)
                                {
                                    for (int ipol1 = 0; ipol1 < npol; ++ipol1)
                                    {
                                        const int k = m0 * m + m1;
                                        double val = 0.0;
                                        if (ipol0 == 0 && ipol1 == 0)
                                        {
                                            val = 0.5 * (blocks[0][k] + blocks[3][k]);
                                        }
                                        else if (ipol0 == 1 && ipol1 == 1)
                                        {
                                            val = 0.5 * (blocks[0][k] - blocks[3][k]);
                                        }
                                        else
                                        {
                                            val = 0.5 * blocks[1][k];
                                        }
                                        ofs << std::setw(12) << std::setprecision(8)
                                            << std::fixed << val;
                                    }
                                }
                                ofs << std::endl;
                            }
                        }
                    }
                }
            }     // l
        }         // I
    }             // T

    return;
}


} // namespace DFTU_BASE
