//----------------------------------------------------------
// EXPLAIN : This routine calculates rhoa as the
// superposition of atomic charges.
//
// nspina is the number of spin components to be calculated
//
// nspina = 1 the total atomic charge density is calculated
// nspina = 2 the spin up and spin down atomic charge
// densities are calculated assuming an uniform atomic
// spin-polarization equal to starting_mag(nt)
// nspina = 4 noncollinear case. The total density is
// calculated in the first component and the magnetization
// std::vector in the other three.
//
// NB: nspina may not be equal to nspin because in some cases
// (as in update) the total charge only could be needed,
// even in a LSDA calculation.
//----------------------------------------------------------
#include "charge.h"
#include "chg_tools.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/libm/libm.h"
#include "source_base/math_integral.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_threading.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

#include <vector>

Charge::Charge()
{
    allocate_rho = false;
    allocate_rho_final_scf = false; // LiuXh add 20180619
}

Charge::~Charge()
{
    this->destroy();
}

void Charge::set_rhopw(ModulePW::PW_Basis* rhopw_in)
{
    this->rhopw = rhopw_in;
}

// mohan add 2025-12-02
bool Charge::kin_density() const
{
    if (XC_Functional::get_ked_flag() || PARAM.inp.out_elf[0] > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Charge::destroy()
{
    if (allocate_rho || allocate_rho_final_scf) // LiuXh add 20180619
    {
        // All storage (rho, rhog, rho_core, etc.) is backed by std::vector
        // members that self-manage; just clear the vectors.
        _ptrs_rho.clear();
        _ptrs_rhog.clear();
        _ptrs_rho_save.clear();
        _ptrs_rhog_save.clear();
        _ptrs_kin_r.clear();
        _ptrs_kin_r_save.clear();
        _space_rho_core.clear();
        _space_rhog_core.clear();
        rho = nullptr;
        rhog = nullptr;
        rho_save = nullptr;
        rhog_save = nullptr;
        rho_core = nullptr;
        rhog_core = nullptr;
        kin_r = nullptr;
        kin_r_save = nullptr;
    }
}

void Charge::allocate(const int& nspin_in, const bool kin_den)
{
    ModuleBase::TITLE("Charge", "allocate");

    if (this->rhopw == nullptr)
    {
        ModuleBase::WARNING_QUIT("Charge::allocate","rhopw is nullptr.");
    }

    this->nrxx = this->rhopw->nrxx;
    this->nxyz = this->rhopw->nxyz;
    this->ngmc = this->rhopw->npw;


    if (allocate_rho == true)
    {
        this->destroy();
        allocate_rho = false;
    }

    assert(allocate_rho == false);

    //  mohan add 2021-02-20
    this->nspin = nspin_in;

    if (PARAM.inp.test_charge > 1)
    {
        std::cout << "\n spin_number = " << nspin << " real_point_number = " << nrxx << std::endl;
    }

    // allocate memory (std::vector self-manages the storage)
    _space_rho.resize(nspin * nrxx);
    _space_rho_save.resize(nspin * nrxx);
    _space_rhog.resize(nspin * ngmc);
    _space_rhog_save.resize(nspin * ngmc);
    if(kin_den)
    {
        _space_kin_r.resize(nspin * nrxx);
        _space_kin_r_save.resize(nspin * nrxx);
    }
    _ptrs_rho.resize(nspin);
    _ptrs_rhog.resize(nspin);
    _ptrs_rho_save.resize(nspin);
    _ptrs_rhog_save.resize(nspin);
    rho = _ptrs_rho.data();
    rhog = _ptrs_rhog.data();
    rho_save = _ptrs_rho_save.data();
    rhog_save = _ptrs_rhog_save.data();
    if(kin_den)
    {
        _ptrs_kin_r.resize(nspin);
        _ptrs_kin_r_save.resize(nspin);
        kin_r = _ptrs_kin_r.data();
        kin_r_save = _ptrs_kin_r_save.data();
    }
    for (int is = 0; is < nspin; is++)
    {
        rho[is] = _space_rho.data() + is * nrxx;
        rhog[is] = _space_rhog.data() + is * ngmc;
        rho_save[is] = _space_rho_save.data() + is * nrxx;
        rhog_save[is] = _space_rhog_save.data() + is * ngmc;
        ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog[is], ngmc);
        ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog_save[is], ngmc);
        if(kin_den) 
        {
            kin_r[is] = _space_kin_r.data() + is * nrxx;
            ModuleBase::GlobalFunc::ZEROS(kin_r[is], nrxx);
            kin_r_save[is] = _space_kin_r_save.data() + is * nrxx;
            ModuleBase::GlobalFunc::ZEROS(kin_r_save[is], nrxx);
        }
    }

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * nspin * nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * nspin * nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * nspin * ngmc);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * nspin * ngmc);
    if(kin_den)
    {
        ModuleBase::Memory::record("Chg::kin_r", sizeof(double) * nspin * ngmc);
        ModuleBase::Memory::record("Chg::kin_r_save", sizeof(double) * nspin * ngmc);
    }

    _space_rho_core.resize(nrxx);
    this->rho_core = _space_rho_core.data();
    ModuleBase::GlobalFunc::ZEROS(rho_core, nrxx);

    _space_rhog_core.resize(ngmc);
    this->rhog_core = _space_rhog_core.data();
    ModuleBase::GlobalFunc::ZEROS(rhog_core, ngmc);

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * ngmc);

    this->allocate_rho = true;
    return;
}

double Charge::sum_rho() const
{
    const int nspin0 = (nspin == 2) ? 2 : 1;
    return module_charge::sum_rho(this->rho, nspin0, this->nrxx, *this->omega_, this->rhopw->nxyz);
}

void Charge::renormalize_rho()
{
    ModuleBase::TITLE("Charge", "renormalize_rho");

    const double sr = this->sum_rho();
    GlobalV::ofs_warning << std::setprecision(15);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge before normalized", sr);
    const double normalize_factor = PARAM.inp.nelec / sr;

    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            rho[is][ir] *= normalize_factor;
        }
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning, "charge after normalized", this->sum_rho());

    GlobalV::ofs_running << std::setprecision(6);
    return;
}

void Charge::save_rho_before_sum_band()
{
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::DCOPY(rho[is], rho_save[is], this->rhopw->nrxx);
        if (XC_Functional::get_ked_flag())
        {
            ModuleBase::GlobalFunc::DCOPY(kin_r[is], kin_r_save[is], this->rhopw->nrxx);
        }
    }
    return;
}

double Charge::cal_rho2ne(const double* rho_in) const
{
    return module_charge::cal_rho2ne(rho_in, this->rhopw->nrxx, *this->omega_, this->rhopw->nxyz);
}

void Charge::check_rho()
{
    if (this->nspin==1 || this->nspin==4)
    {
        double ne = 0.0;
        ne = this->cal_rho2ne(rho[0]);
        if (std::abs(ne - PARAM.inp.nelec) > 1.0e-6)
        {
            ModuleBase::WARNING("Charge", "Charge is not equal to the number of electrons!");
        }
    }
    else if (this->nspin == 2)
    {
        // for spin up
        double ne_up = 0.0;
        ne_up = this->cal_rho2ne(rho[0]);
        if (ne_up < 0.0)
        {
            ModuleBase::WARNING_QUIT("Charge",
                "Number of spin-down electrons set in starting magnetization exceeds all available.");
        }
        // for spin down
        double ne_dn = 0.0;
        ne_dn = this->cal_rho2ne(rho[1]);
        if (ne_dn < 0.0)
        {
            ModuleBase::WARNING_QUIT("Charge",
                "Number of spin-up electrons set in starting magnetization exceeds all available.");
        }
        // for total charge
        if (std::abs(ne_up + ne_dn - PARAM.inp.nelec) > 1.0e-6)
        {
            ModuleBase::WARNING("Charge", "Charge is not equal to the number of electrons!");
        }
    }
}

// LiuXh add 20180619
void Charge::init_final_scf()
{
    ModuleBase::TITLE("Charge", "init_after_scf");

    assert(allocate_rho_final_scf == false);
    if (PARAM.inp.test_charge > 1)
    {
        std::cout << "\n spin_number = " << PARAM.inp.nspin
                  << " real_point_number = " << this->rhopw->nrxx << std::endl;
    }

    // allocate memory for final SCF (std::vector self-manages storage)
    const int ns = PARAM.inp.nspin;
    const int nrxx = this->rhopw->nrxx;
    const int ngmc = this->rhopw->npw;
    _space_rho.resize(ns * nrxx);
    _space_rho_save.resize(ns * nrxx);
    _space_rhog.resize(ns * ngmc);
    _space_rhog_save.resize(ns * ngmc);
    _ptrs_rho.resize(ns);
    _ptrs_rhog.resize(ns);
    _ptrs_rho_save.resize(ns);
    _ptrs_rhog_save.resize(ns);
    rho = _ptrs_rho.data();
    rhog = _ptrs_rhog.data();
    rho_save = _ptrs_rho_save.data();
    rhog_save = _ptrs_rhog_save.data();

    for (int is = 0; is < ns; is++)
    {
        rho[is] = _space_rho.data() + is * nrxx;
        rhog[is] = _space_rhog.data() + is * ngmc;
        rho_save[is] = _space_rho_save.data() + is * nrxx;
        rhog_save[is] = _space_rhog_save.data() + is * ngmc;
        ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog[is], ngmc);
        ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
        ModuleBase::GlobalFunc::ZEROS(rhog_save[is], ngmc);
    }

    ModuleBase::Memory::record("Chg::rho", sizeof(double) * ns * nrxx);
    ModuleBase::Memory::record("Chg::rho_save", sizeof(double) * ns * nrxx);
    ModuleBase::Memory::record("Chg::rhog", sizeof(double) * ns * ngmc);
    ModuleBase::Memory::record("Chg::rhog_save", sizeof(double) * ns * ngmc);

    _space_rho_core.resize(nrxx);
    this->rho_core = _space_rho_core.data();
    ModuleBase::GlobalFunc::ZEROS(rho_core, nrxx);

    _space_rhog_core.resize(ngmc);
    this->rhog_core = _space_rhog_core.data();
    ModuleBase::GlobalFunc::ZEROS(rhog_core, ngmc);

    ModuleBase::Memory::record("Chg::rho_core", sizeof(double) * this->rhopw->nrxx);
    ModuleBase::Memory::record("Chg::rhog_core", sizeof(double) * this->rhopw->npw);

    this->allocate_rho_final_scf = true;
    return;
}
