#include "source_pw/module_pwdft/dftu_base.h"

#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base_io.h"
#include "source_base/global_function.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_global.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

// local inline helpers for eigenvalue calculation (JacobiRotate, CalculateEigenvalues)
// have been migrated to dftu_base_io.cpp, where they are used by DFTU_BASE::write_occup_m.
// mohan refactored 2025-11-08
// All members are now non-static; default values are in the header.


Plus_U_Base::Plus_U_Base()
{}


Plus_U_Base::~Plus_U_Base()
{}


// init_base: base-only initialization shared by PW and LCAO paths.
// LCAO-specific setup (paraV, orb, ucell pointer) stays in Plus_U::init().
void Plus_U_Base::init_base(UnitCell& cell,
                             const int npol,
                             const int nspin,
                             const std::vector<int>& orbital_corr,
                             const bool yukawa_potential,
                             const double yukawa_lambda,
                             const std::string& global_readin_dir,
                             const std::string& global_out_dir,
                             const std::string& init_chg,
                             const std::string& device,
                             const int kpar,
                             const std::vector<double>& hubbard_u,
                             const double uramping,
                             const int occ_mat_ctrl,
                             const int mixing_dftu)
{
    ModuleBase::TITLE("Plus_U_Base", "init_base");

#ifndef __MPI
    std::cout << "DFT+U module is only accessible in mpi version" << std::endl;
    exit(0);
#endif

    this->nspin = nspin;
    this->orbital_corr = orbital_corr;
    this->uramping = uramping;
    this->occ_mat_ctrl = occ_mat_ctrl;
    this->mixing_dftu = mixing_dftu;
    this->u_target = hubbard_u;
    this->u_current = hubbard_u;
    if (uramping > 0.01)
    {
        std::fill(this->u_current.begin(),
                  this->u_current.end(),
                  0.0);
    }
    this->device = device;
    this->kpar = kpar;

    this->energy_u = 0.0;

    this->occmat_.init(cell, orbital_corr, nspin, npol);

    this->occ_mat.resize(cell.nat);
    this->occ_mat_save.resize(cell.nat);
    this->pot_uterm_pw_index.resize(cell.nat);
    int pot_index = 0;

    this->iatlnmipol2iwt.resize(cell.nat);

    int num_locale = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            const int iat = cell.itia2iat(it, ia);

            occ_mat[iat].resize(cell.atoms[it].nwl + 1);
            occ_mat_save[iat].resize(cell.atoms[it].nwl + 1);

            this->iatlnmipol2iwt[iat].resize(cell.atoms[it].nwl + 1);

            if(!has_correlated_orbital(it))
            {
                continue;
            }

            const int tlp1_npol = (get_orbital_corr(it)*2+1)*npol;
            const int tlp1 = 2 * get_orbital_corr(it) + 1;
            const int elem_size = tlp1 * tlp1;
            if(nspin == 4)
            {
                this->pot_uterm_pw_index[iat] = pot_index;
                pot_index += tlp1_npol * tlp1_npol;
            }
            else
            {
                this->pot_uterm_pw_index[iat] = pot_index;
                pot_index += elem_size;
            }

            for (int l = 0; l <= cell.atoms[it].nwl; l++)
            {
                const int N = cell.atoms[it].l_nchi[l];

                occ_mat[iat][l].resize(N);
                occ_mat_save[iat][l].resize(N);

                for (int n = 0; n < N; n++)
                {
                    if (nspin == 1 || nspin == 2)
                    {
                        occ_mat[iat][l][n].resize(2);
                        occ_mat_save[iat][l][n].resize(2);

                        occ_mat[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        occ_mat[iat][l][n][1].create(2 * l + 1, 2 * l + 1);

                        occ_mat_save[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        occ_mat_save[iat][l][n][1].create(2 * l + 1, 2 * l + 1);
                        num_locale += (2 * l + 1) * (2 * l + 1) * 2;
                    }
                    else if (nspin == 4)
                    {
                        occ_mat[iat][l][n].resize(1);
                        occ_mat_save[iat][l][n].resize(1);

                        occ_mat[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                        occ_mat_save[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                        num_locale += (2 * l + 1) * (2 * l + 1) * npol * npol;
                    }
                }
            }

            this->iatlnmipol2iwt[iat].resize(cell.atoms[it].nwl + 1);
            for (int L = 0; L <= cell.atoms[it].nwl; L++)
            {
                this->iatlnmipol2iwt[iat][L].resize(cell.atoms[it].l_nchi[L]);

                for (int n = 0; n < cell.atoms[it].l_nchi[L]; n++)
                {
                    this->iatlnmipol2iwt[iat][L][n].resize(2 * L + 1);

                    for (int m = 0; m < 2 * L + 1; m++)
                    {
                        this->iatlnmipol2iwt[iat][L][n][m].resize(npol);
                    }
                }
            }

            for (int iw = 0; iw < cell.atoms[it].nw * npol; iw++)
            {
                int iw0 = iw / npol;
                int ipol = iw % npol;
                int iwt = cell.itiaiw2iwt(it, ia, iw);
                int l = cell.atoms[it].iw2l[iw0];
                int n = cell.atoms[it].iw2n[iw0];
                int m = cell.atoms[it].iw2m[iw0];

                this->iatlnmipol2iwt[iat][l][n][m][ipol] = iwt;
            }
        }
    }

    if (nspin == 2) pot_index *= 2;

    this->pot_uterm_pw.resize(pot_index, 0.0);
    this->uom_array.resize(pot_index, 0.0);
    this->uom_save.resize(pot_index, 0.0);

    if (yukawa_potential)
    {
        this->yukawa_.reset(new YukawaScreening());
        this->yukawa_->init(cell, orbital_corr, yukawa_lambda);
    }

    if (occ_mat_ctrl != 0)
    {
        std::stringstream sst;
        sst << global_readin_dir << "dm_onsite_ini.txt";
        DFTU_BASE::read_occup_m(cell, this->occmat_.data(), this->orbital_corr, this->occ_mat_ctrl,
                                sst.str(), init_chg, nspin, npol);
#ifdef __MPI
        DFTU_BASE::local_occup_bcast(cell, this->occmat_.data(), this->orbital_corr, nspin, npol);
#endif

        mark_occ_mat_initialized();
        this->copy_occ_mat(cell);
    }
    else
    {
        if (init_chg == "file")
        {
            std::stringstream sst;
            sst << global_readin_dir << "dm_onsite.txt";
            DFTU_BASE::read_occup_m(cell, this->occmat_.data(), this->orbital_corr, this->occ_mat_ctrl,
                                    sst.str(), init_chg, nspin, npol);
#ifdef __MPI
            DFTU_BASE::local_occup_bcast(cell, this->occmat_.data(), this->orbital_corr, nspin, npol);
#endif
            mark_occ_mat_initialized();
        }
        else
        {
            this->zero_occ_mat(cell);
        }
    }

    ModuleBase::Memory::record("Plus_U_Base::occ_mat", sizeof(double) * num_locale);
    return;
}


void Plus_U_Base::uramping_update()
{
    // Yukawa calculates U directly every iteration, no need for ramping
    if (use_yukawa()) {
        return;
    }
    // if uramping < 0.1, use the original U
    if (this->uramping < 0.01) {
        return;
    }
    // loop to change U
    for (int i = 0; i < static_cast<int>(this->u_target.size()); i++)
    {
        if (this->u_current[i] + this->uramping < this->u_target[i])
        {
            this->u_current[i] += this->uramping;
        }
        else
        {
            this->u_current[i] = this->u_target[i];
        }
    }
}


bool Plus_U_Base::u_converged()
{
    // Yukawa calculates U directly every iteration, always considered converged
    if (use_yukawa()) {
        return true;
    }
    for (int i = 0; i < static_cast<int>(this->u_target.size()); i++)
    {
        if (this->u_current[i] != this->u_target[i])
        {
            return false;
        }
    }
    return true;
}


// copy_occ_mat — save current occ to occ_save and uom_save
void Plus_U_Base::copy_occ_mat(const UnitCell& ucell)
{
    this->occmat_.copy_to_save(ucell, this->orbital_corr);
    this->occmat_.write_save_to_flat(ucell, this->orbital_corr,
                                     this->pot_uterm_pw_index, this->uom_save);
}


void Plus_U_Base::zero_occ_mat(const UnitCell& ucell)
{
    this->occmat_.zero(ucell, this->orbital_corr);
}


void Plus_U_Base::set_occ_mat(const UnitCell& ucell)
{
    this->occmat_.read_from_flat(ucell, this->orbital_corr,
                                 this->pot_uterm_pw_index, this->uom_array);
}


void Plus_U_Base::get_occ_mat_flat(const int iat, const int l, std::vector<double>& occ) const
{
    this->occmat_.get_flat(iat, l, occ);
}


void Plus_U_Base::set_occ_mat_flat(const int iat, const int l, const int spin,
                                   const std::vector<double>& occ)
{
    this->occmat_.set_flat(iat, l, spin, occ);
}


// cal_occ_pw() is implemented in source_pw/module_pwdft/dftu_base_occ.cpp
// as a Plus_U_Base method. Pure per-atom kernels live in dftu_base_tools.{h,cpp}
// as free functions in namespace DFTU_BASE.
