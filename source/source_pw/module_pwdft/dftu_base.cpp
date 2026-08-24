#include "source_pw/module_pwdft/dftu_base.h"

#include "source_base/parallel_global.h"
#include "source_base/timer.h"

#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <vector>

// local inline helpers for eigenvalue calculation (mirrors dftu_io.cpp)
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

// static member definitions (mohan add 2025-11-06)
double Plus_U_Base::energy_u = 0.0;

std::vector<double> Plus_U_Base::U = {};

std::vector<double> Plus_U_Base::U0 = {};

std::vector<int> Plus_U_Base::orbital_corr = {};

double Plus_U_Base::uramping = 0.0;

int Plus_U_Base::omc = 0;

int Plus_U_Base::mixing_dftu = 0;
int Plus_U_Base::nspin = 0;

bool Plus_U_Base::Yukawa = false;


Plus_U_Base::Plus_U_Base()
{}


Plus_U_Base::~Plus_U_Base()
{}


// init_base will be fully implemented in Phase 4
void Plus_U_Base::init_base(UnitCell& cell,
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
                             const int kpar)
{
    (void)cell;
    (void)npol;
    (void)nspin;
    (void)orbital_corr;
    (void)yukawa_potential;
    (void)yukawa_lambda;
    (void)global_readin_dir;
    (void)global_out_dir;
    (void)init_chg;
    (void)nlocal;
    (void)gamma_only_local;
    (void)ks_solver;
    (void)cal_force;
    (void)cal_stress;
    (void)device;
    (void)kpar;
}


void Plus_U_Base::uramping_update()
{
    // Yukawa calculates U directly every iteration, no need for ramping
    if (Yukawa) {
        return;
    }
    // if uramping < 0.1, use the original U
    if (this->uramping < 0.01) {
        return;
    }
    // loop to change U
    for (int i = 0; i < static_cast<int>(this->U0.size()); i++)
    {
        if (this->U[i] + this->uramping < this->U0[i])
        {
            this->U[i] += this->uramping;
        }
        else
        {
            this->U[i] = this->U0[i];
        }
    }
}


bool Plus_U_Base::u_converged()
{
    // Yukawa calculates U directly every iteration, always considered converged
    if (Yukawa) {
        return true;
    }
    for (int i = 0; i < static_cast<int>(this->U0.size()); i++)
    {
        if (this->U[i] != this->U0[i])
        {
            return false;
        }
    }
    return true;
}


// copy_locale — save current locale to locale_save and uom_save
void Plus_U_Base::copy_locale(const UnitCell& ucell)
{
    ModuleBase::TITLE("Plus_U_Base", "copy_locale");
    ModuleBase::timer::start("Plus_U_Base", "copy_locale");

    for (int T = 0; T < ucell.ntype; T++)
    {
        int target_l = get_orbital_corr(T);
        if (target_l == -1)
            continue;

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            if (Plus_U_Base::nspin == 4)
            {
                locale_save[iat][target_l][0][0] = locale[iat][target_l][0][0];
                if(this->uom_save.size() != 0)
                {
                    const int size = locale[iat][target_l][0][0].nr * locale[iat][target_l][0][0].nc;
                    for(int mm=0; mm<size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat]+mm] = locale[iat][target_l][0][0].c[mm];
                    }
                }
            }
            else if (Plus_U_Base::nspin == 1 || Plus_U_Base::nspin == 2)
            {
                locale_save[iat][target_l][0][0] = locale[iat][target_l][0][0];
                locale_save[iat][target_l][0][1] = locale[iat][target_l][0][1];
                if(this->uom_save.size() != 0)
                {
                    const int size = locale[iat][target_l][0][0].nr * locale[iat][target_l][0][0].nc;
                    const int half_size = this->uom_save.size() / 2;
                    for(int mm=0; mm<size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat]+mm] = locale[iat][target_l][0][0].c[mm];
                        this->uom_save[half_size + eff_pot_pw_index[iat]+mm] = locale[iat][target_l][0][1].c[mm];
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("Plus_U_Base", "copy_locale");
}


void Plus_U_Base::zero_locale(const UnitCell& ucell)
{
    ModuleBase::TITLE("Plus_U_Base", "zero_locale");
    ModuleBase::timer::start("Plus_U_Base", "zero_locale");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T))
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < ucell.atoms[T].nwl + 1; l++)
            {
                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (Plus_U_Base::nspin == 4)
                    {
                        locale[iat][l][n][0].zero_out();
                    }
                    else if (Plus_U_Base::nspin == 1 || Plus_U_Base::nspin == 2)
                    {
                        locale[iat][l][n][0].zero_out();
                        locale[iat][l][n][1].zero_out();
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("Plus_U_Base", "zero_locale");
}


void Plus_U_Base::mix_locale(const UnitCell& ucell,
                              const double& mixing_beta)
{
    ModuleBase::TITLE("Plus_U_Base", "mix_locale");
    ModuleBase::timer::start("Plus_U_Base", "mix_locale");

    double beta = mixing_beta;

    for (int T = 0; T < ucell.ntype; T++)
    {
        int target_l = get_orbital_corr(T);
        if (target_l == -1)
            continue;

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            if (Plus_U_Base::nspin == 4)
            {
                const int size = locale[iat][target_l][0][0].nr * locale[iat][target_l][0][0].nc;
                for (int mm = 0; mm < size; mm++)
                {
                    locale[iat][target_l][0][0].c[mm] = locale[iat][target_l][0][0].c[mm] * beta + locale_save[iat][target_l][0][0].c[mm] * (1.0 - beta);
                }
                if (this->uom_save.size() != 0)
                {
                    for (int mm = 0; mm < size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat] + mm] = locale[iat][target_l][0][0].c[mm];
                    }
                }
            }
            else if (Plus_U_Base::nspin == 1 || Plus_U_Base::nspin == 2)
            {
                const int size = locale[iat][target_l][0][0].nr * locale[iat][target_l][0][0].nc;
                const int half_size = this->uom_save.size() / 2;
                for (int mm = 0; mm < size; mm++)
                {
                    locale[iat][target_l][0][0].c[mm] = locale[iat][target_l][0][0].c[mm] * beta + locale_save[iat][target_l][0][0].c[mm] * (1.0 - beta);
                    locale[iat][target_l][0][1].c[mm] = locale[iat][target_l][0][1].c[mm] * beta + locale_save[iat][target_l][0][1].c[mm] * (1.0 - beta);
                }
                if (this->uom_save.size() != 0)
                {
                    for (int mm = 0; mm < size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat] + mm] = locale[iat][target_l][0][0].c[mm];
                        this->uom_save[half_size + eff_pot_pw_index[iat] + mm] = locale[iat][target_l][0][1].c[mm];
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("Plus_U_Base", "mix_locale");
}


void Plus_U_Base::set_locale(const UnitCell& ucell)
{
    ModuleBase::TITLE("Plus_U_Base", "set_locale");
    ModuleBase::timer::start("Plus_U_Base", "set_locale");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T)) continue;
        const int l = get_orbital_corr(T);
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            if (Plus_U_Base::nspin == 4)
            {
                for(int mm = 0; mm < locale[iat][l][0][0].nr * locale[iat][l][0][0].nc; mm++)
                    locale[iat][l][0][0].c[mm] = this->uom_array[eff_pot_pw_index[iat] + mm];
            }
            else if (Plus_U_Base::nspin == 1 || Plus_U_Base::nspin == 2)
            {
                const int half_size = this->uom_array.size() / 2;
                for(int mm = 0; mm < locale[iat][l][0][0].nr * locale[iat][l][0][0].nc; mm++)
                {
                    locale[iat][l][0][0].c[mm] = this->uom_array[eff_pot_pw_index[iat] + mm];
                    if (Plus_U_Base::nspin == 2)
                    {
                        locale[iat][l][0][1].c[mm] = this->uom_array[half_size + eff_pot_pw_index[iat] + mm];
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("Plus_U_Base", "set_locale");
}


void Plus_U_Base::get_locale_flat(const int iat, const int l, std::vector<double>& occ) const
{
    const int tlp1 = 2 * l + 1;
    const int size = tlp1 * tlp1;
    if (nspin == 2)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int i = 0; i < size; i++)
            {
                occ[is * size + i] = locale[iat][l][0][is].c[i];
            }
        }
    }
    else
    {
        for (int i = 0; i < static_cast<int>(occ.size()); i++)
        {
            occ[i] = locale[iat][l][0][0].c[i];
        }
    }
}


void Plus_U_Base::set_locale_flat(const int iat, const int l, const int spin,
                                   const std::vector<double>& occ)
{
    for (int i = 0; i < static_cast<int>(occ.size()); i++)
    {
        locale[iat][l][0][spin].c[i] = occ[i];
    }
}


void Plus_U_Base::output(const UnitCell& ucell,
                         bool out_chg,
                         const std::string& global_out_dir,
                         int nspin,
                         int npol)
{
    ModuleBase::TITLE("Plus_U_Base", "output");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | #DFT+U INFORMATION# |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    for (int T = 0; T < ucell.ntype; T++)
    {
        const int NL = ucell.atoms[T].nwl + 1;

        for (int L = 0; L < NL; L++)
        {
            const int N = ucell.atoms[T].l_nchi[L];

            if (L >= get_orbital_corr(T) && has_correlated_orbital(T))
            {
                if (L != get_orbital_corr(T))
                {
                    continue;
                }

                if (!Yukawa)
                {
                    GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << " ORBITAL=" << 0
                                         << " U=" << this->U[T] * ModuleBase::Ry_to_eV << " eV" << std::endl;
                }
                else
                {
                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                        {
                            continue;
                        }
                        double Ueff = (this->U_Yukawa[T][L][n] - this->J_Yukawa[T][L][n]) * ModuleBase::Ry_to_eV;
                        GlobalV::ofs_running << " Type=" << T+1 << " L=" << L << "  ORBITAL=" << n
                                             << " U=" << this->U_Yukawa[T][L][n] * ModuleBase::Ry_to_eV << " eV"
                                             << " J=" << this->J_Yukawa[T][L][n] * ModuleBase::Ry_to_eV << " eV"
                                             << std::endl;
                    }
                }
            }
        }
    }

    GlobalV::ofs_running << " Local Occupation Matrices for each atom" << std::endl;
    this->write_occup_m(ucell, GlobalV::ofs_running, true, nspin, npol);

    // Write dm_onsite.txt
    if (out_chg && GlobalV::MY_RANK == 0)
    {
        std::ofstream ofdftu;
        ofdftu.open(global_out_dir + "dm_onsite.txt");
        if (!ofdftu)
        {
            ModuleBase::WARNING_QUIT("Plus_U_Base::output", "Can't create file dm_onsite.txt");
        }
        this->write_occup_m(ucell, ofdftu, false, nspin, npol);
        ofdftu.close();
    }

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " | # END DFT+U INFO    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>" << std::endl << std::endl;

    return;
}


void Plus_U_Base::write_occup_m(const UnitCell& ucell,
                                std::ofstream& ofs,
                                bool diag,
                                int nspin,
                                int npol)
{
    ModuleBase::TITLE("Plus_U_Base", "write_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T))
        {
            continue;
        }
        const int NL = ucell.atoms[T].nwl + 1;
        const int LC = get_orbital_corr(T);

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < NL; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    ofs << "\n Atom= " << iat+1;
                    ofs << " L= " << l;
                    ofs << " ORBITAL= " << n << std::endl;

                    if (nspin == 1 || nspin == 2)
                    {
                        double sum0[2];
                        for (int is = 0; is < 2; is++)
                        {
                            if (diag)
                            {
                                std::vector<std::vector<double>> A(2 * l + 1, std::vector<double>(2 * l + 1));
                                for (int m0 = 0; m0 < 2 * l + 1; m0++)
                                {
                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        A[m0][m1] = locale[iat][l][n][is](m0, m1);
                                    }
                                }
                                std::vector<double> eigenvalues = CalculateEigenvalues(A, 2 * l + 1);
                                sum0[is] = 0.0;
                                ofs << " Eigenvalues for spin=" << is+1 << std::endl;
                                ofs << std::setprecision(8) << std::fixed;
                                for (int i = 0; i < 2 * l + 1; i++)
                                {
                                    ofs << std::setw(12) << eigenvalues[i];
                                    sum0[is] += eigenvalues[i];
                                }
                                ofs << std::endl;
                                ofs << " sum is " << std::setw(12) << sum0[is] << std::endl;
                            }
                            ofs << " spin= " << is+1 << std::endl;
                            ofs << std::setprecision(8) << std::fixed;
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    ofs << std::setw(12)
                                        << locale[iat][l][n][is](m0, m1);
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
                        if (diag)
                        {
                            double sum0[4];
                            std::vector<std::vector<double>> A(2 * l + 1, std::vector<double>(2 * l + 1));
                            int index = 0;
                            for (int is = 0; is < 4; is++)
                            {
                                for (int m0 = 0; m0 < 2 * l + 1; m0++)
                                {
                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        A[m0][m1] = locale[iat][l][n][0].c[index];
                                        index++;
                                    }
                                }
                                std::vector<double> eigenvalues = CalculateEigenvalues(A, 2 * l + 1);
                                sum0[is] = 0.0;
                                ofs << " Eigenvalues for is=" << is << std::endl;
                                ofs << std::setprecision(8) << std::fixed;
                                for (int i = 0; i < 2 * l + 1; i++)
                                {
                                    ofs << std::setw(12) << eigenvalues[i];
                                    sum0[is] += eigenvalues[i];
                                }
                                ofs << std::endl;
                                ofs << " sum is " << std::setw(12) << sum0[is] << std::endl;
                            }
                            ofs << std::setw(12) << std::setprecision(8)
                                << std::fixed << " Magnetism for atom " << iat + 1 << ": "
                                << sum0[1] << " " << sum0[2] << " " << sum0[3] << std::endl;
                        }
                        else
                        {
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int ipol0 = 0; ipol0 < npol; ipol0++)
                                {
                                    const int m0_all = m0 + (2 * l + 1) * ipol0;

                                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                    {
                                        for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                        {
                                            int m1_all = m1 + (2 * l + 1) * ipol1;
                                            ofs << std::setw(12) << std::setprecision(8) << std::fixed
                                                << locale[iat][l][n][0](m0_all, m1_all);
                                        }
                                    }
                                    ofs << std::endl;
                                }
                            }
                        }
                    }
                } // n
            }     // l
        }         // I
    }             // T

    return;
}


void Plus_U_Base::read_occup_m(const UnitCell& ucell,
                               const std::string& fn,
                               const std::string& init_chg,
                               int nspin,
                               int npol)
{
    ModuleBase::TITLE("Plus_U_Base", "read_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ifstream ifdftu(fn.c_str(), std::ios::in);

    if (!ifdftu)
    {
        if (omc > 0)
        {
            ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "Can not find the file dm_onsite_ini.txt. Please check your dm_onsite_ini.txt");
        }
        else
        {
            if (init_chg == "file")
            {
                ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "Can not find the file dm_onsite.txt. Please do scf calculation first");
            }
        }
        ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "Can not open dm_onsite.txt file");
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
                ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> L;
            ifdftu >> word;

            if (strcmp("ORBITAL=", word) != 0)
            {
                ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> zeta;
            ifdftu.ignore(150, '\n');

            T = ucell.iat2it[iat];
            const int NL = ucell.atoms[T].nwl + 1;
            const int LC = get_orbital_corr(T);

            for (int l = 0; l < NL; l++)
            {
                if (l != get_orbital_corr(T))
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
                                    locale[iat][L][zeta][spin](m0, m1) = value;
                                }
                                ifdftu.ignore(150, '\n');
                            }
                        }
                        else
                        {
                            ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
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
                                    locale[iat][L][zeta][0](m0_all, m1_all) = value;
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
            ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
        }

        ifdftu.rdstate();

        if (ifdftu.eof() != 0)
        {
            break;
        }
    }

    return;
}


void Plus_U_Base::local_occup_bcast(const UnitCell& ucell,
                                    int nspin,
                                    int npol)
{
    ModuleBase::TITLE("Plus_U_Base", "local_occup_bcast");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T))
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            const int L = get_orbital_corr(T);

            for (int l = 0; l <= ucell.atoms[T].nwl; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                for (int n = 0; n < ucell.atoms[T].l_nchi[l]; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    if (nspin == 1 || nspin == 2)
                    {
                        for (int spin = 0; spin < 2; spin++)
                        {
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
#ifdef __MPI
                                    MPI_Bcast(&locale[iat][l][n][spin](m0, m1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
                                }
                            }
                        }
                    }
                    else if (nspin == 4) // SOC
                    {
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
#ifdef __MPI
                                        MPI_Bcast(&locale[iat][l][n][0](m0_all, m1_all),
                                                  1,
                                                  MPI_DOUBLE,
                                                  0,
                                                  MPI_COMM_WORLD);
#endif
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}


// cal_occ_pw() is implemented in source_lcao/module_dftu/dftu_pw.cpp
// as a Plus_U_Base method. It will be relocated to this directory
// in Phase 5 of the class-split refactor.
