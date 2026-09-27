#include "hs_r_legacy.h"

#include "source_base/timer.h"
#include "source_lcao/lcao_hs_arrays.h"
#include "source_lcao/spar_dh.h"
#include "source_lcao/spar_hsr.h"
#include "source_lcao/spar_st.h"
#include "dhs_sparse_writer.h"
#include "hs_sparse_io.h"

#include <complex>
#include <fstream>
#include <sstream>
#include <vector>

// if 'binary=true', output binary file.
// The 'sparse_thr' is the accuracy of the sparse matrix.
// If the absolute value of the matrix element is less than or equal to the
// 'sparse_thr', it will be ignored.

void ModuleIO::output_dSR(const int& istep,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const bool& binary,
                          const double& sparse_thr,
                          const int precision,
                          const std::string& global_out_dir,
                          const std::string& global_matrix_dir,
                          const std::string& calculation,
                          const bool out_app_flag,
                          const int nspin,
                          const bool gamma_only_local,
                          const int npol,
                          const int nlocal)
{
    ModuleBase::TITLE("ModuleIO", "output_dSR");
    ModuleBase::timer::start("ModuleIO", "output_dSR");

    sparse_format::cal_dS(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr,
                          gamma_only_local, nspin, npol);

    // mohan update 2024-04-01
    const std::string fileflag_s = "s";
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary, fileflag_s, precision,
                             global_out_dir, global_matrix_dir, calculation, out_app_flag,
                             nspin, nlocal);

    sparse_format::destroy_dH_R_sparse(HS_Arrays, nspin);

    ModuleBase::timer::end("ModuleIO", "output_dSR");
    return;
}

void ModuleIO::output_dHR(const int& istep,
                          const ModuleBase::matrix& v_eff,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const bool& binary,
                          const double& sparse_thr,
                          const int precision,
                          const std::string& global_out_dir,
                          const std::string& global_matrix_dir,
                          const std::string& calculation,
                          const bool out_app_flag,
                          const int nspin,
                          const bool gamma_only_local,
                          const int npol,
                          const int nlocal)
{
    ModuleBase::TITLE("ModuleIO", "output_dHR");
    ModuleBase::timer::start("ModuleIO", "output_dHR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                         #Print out dH/dR#                          |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    if (nspin == 1 || nspin == 4)
    {
        // mohan add 2024-04-01
        const int cspin = 0;

        sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff,
                              gamma_only_local, nspin, npol);
    }
    else if (nspin == 2)
    {
        for (int cspin = 0; cspin < 2; cspin++)
        {
            sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff,
                                  gamma_only_local, nspin, npol);
        }
    }
    // mohan update 2024-04-01
    const std::string fileflag_h = "h";
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary, fileflag_h, precision,
                             global_out_dir, global_matrix_dir, calculation, out_app_flag,
                             nspin, nlocal);

    sparse_format::destroy_dH_R_sparse(HS_Arrays, nspin);

    ModuleBase::timer::end("ModuleIO", "output_dHR");
    return;
}

template <typename TK>
void ModuleIO::output_SR(Parallel_Orbitals& pv,
                         const Grid_Driver& grid,
                         hamilt::Hamilt<TK>* p_ham,
                         const std::string& SR_filename,
                         const bool& binary,
                         const double& sparse_thr,
                         const int precision,
                         const std::string& global_out_dir,
                         const std::string& global_matrix_dir,
                         const std::string& calculation,
                         const bool out_app_flag,
                         const int nspin)
{
    ModuleBase::TITLE("ModuleIO", "output_SR");
    ModuleBase::timer::start("ModuleIO", "output_SR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                 #Print out overlap matrix S(R)#                    |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::cout << " Overlap matrix file is in " << SR_filename << std::endl;
    GlobalV::ofs_running << " Overlap matrix file is in " << SR_filename << std::endl;

    LCAO_HS_Arrays HS_Arrays;

    sparse_format::cal_SR(pv,
                          HS_Arrays.all_R_coor,
                          HS_Arrays.SR_sparse,
                          HS_Arrays.SR_soc_sparse,
                          grid,
                          sparse_thr,
                          p_ham);

    const int istep = 0;
    ModuleIO::SparseWriteOptions options;
    options.filename = SR_filename;
    options.label = "S";
    options.threshold = sparse_thr;
    options.binary = binary;
    options.precision = precision;
    options.istep = istep;
    options.reduce = true;
    options.calculation = calculation;
    options.out_app_flag = out_app_flag;

    if (nspin == 4)
    {
        ModuleIO::save_sparse(HS_Arrays.SR_soc_sparse,
                              HS_Arrays.all_R_coor,
                              pv,
                              options);
    }
    else
    {
        ModuleIO::save_sparse(HS_Arrays.SR_sparse,
                              HS_Arrays.all_R_coor,
                              pv,
                              options);
    }

    sparse_format::destroy_HS_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_SR");
    return;
}

void ModuleIO::output_TR(const int istep,
                         const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         LCAO_HS_Arrays& HS_Arrays,
                         const Grid_Driver& grid,
                         const TwoCenterBundle& two_center_bundle,
                         const LCAO_Orbitals& orb,
                         const std::string& TR_filename,
                         const bool& binary,
                         const double& sparse_thr,
                         const int precision,
                         const std::string& global_out_dir,
                         const std::string& global_matrix_dir,
                         const std::string& calculation,
                         const bool out_app_flag,
                         const int nspin)
{
    ModuleBase::TITLE("ModuleIO", "output_TR");
    ModuleBase::timer::start("ModuleIO", "output_TR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |           #Print out kinetic energy term matrix T(R)#              |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::stringstream sst;
    const bool md_no_append = (calculation == "md") && !out_app_flag;
    if (md_no_append)
    {
        sst << global_matrix_dir << TR_filename << "g" << istep;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }
    else
    {
        sst << global_out_dir << TR_filename;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }

    sparse_format::cal_TR(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr);
    ModuleIO::SparseWriteOptions options;
    options.filename = sst.str();
    options.label = "T";
    options.threshold = sparse_thr;
    options.binary = binary;
    options.precision = precision;
    options.istep = istep;
    options.reduce = true;
    options.calculation = calculation;
    options.out_app_flag = out_app_flag;

    ModuleIO::save_sparse(HS_Arrays.TR_sparse,
                          HS_Arrays.all_R_coor,
                          pv,
                          options);

    sparse_format::destroy_T_R_sparse(HS_Arrays);

    ModuleBase::timer::end("ModuleIO", "output_TR");
    return;
}

template void ModuleIO::output_SR<double>(Parallel_Orbitals& pv,
                                          const Grid_Driver& grid,
                                          hamilt::Hamilt<double>* p_ham,
                                          const std::string& SR_filename,
                                          const bool& binary,
                                          const double& sparse_thr,
                                          const int precision,
                                          const std::string& global_out_dir,
                                          const std::string& global_matrix_dir,
                                          const std::string& calculation,
                                          const bool out_app_flag,
                                          const int nspin);
template void ModuleIO::output_SR<std::complex<double>>(Parallel_Orbitals& pv,
                                                        const Grid_Driver& grid,
                                                        hamilt::Hamilt<std::complex<double>>* p_ham,
                                                        const std::string& SR_filename,
                                                        const bool& binary,
                                                        const double& sparse_thr,
                                                        const int precision,
                                                        const std::string& global_out_dir,
                                                        const std::string& global_matrix_dir,
                                                        const std::string& calculation,
                                                        const bool out_app_flag,
                                                        const int nspin);
