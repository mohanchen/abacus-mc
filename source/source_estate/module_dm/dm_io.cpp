#include "dm_io.h"

#include "density_matrix.h"

#include "source_base/tool_title.h"

#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <string>

namespace module_dm
{

// read *.dmk into density matrix dm(k)
template <typename TK, typename TR>
void read_DMK_file(DensityMatrix<TK, TR>& dm,
                   const std::string& directory,
                   const int ispin,
                   const int ik)
{
    ModuleBase::TITLE("DensityMatrix", "read_DMK");
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= dm._nspin);
#endif
    // read
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
    //
    bool quit_abacus = false;

    std::ifstream ifs;

    ifs.open(fn.c_str());
    if (!ifs)
    {
        quit_abacus = true;
    }
    else
    {
        // if the number is not match,
        // quit the program or not.
        bool quit = false;

        ModuleBase::CHECK_DOUBLE(ifs, dm._kvec_d[ik].x, quit);
        ModuleBase::CHECK_DOUBLE(ifs, dm._kvec_d[ik].y, quit);
        ModuleBase::CHECK_DOUBLE(ifs, dm._kvec_d[ik].z, quit);
        ModuleBase::CHECK_INT(ifs, dm._paraV->nrow);
        ModuleBase::CHECK_INT(ifs, dm._paraV->ncol);
    } // If file exist, read in data.
    // Finish reading the first part of density matrix.

    for (int i = 0; i < dm._paraV->nrow; ++i)
    {
        for (int j = 0; j < dm._paraV->ncol; ++j)
        {
            ifs >> dm._DMK[ik + dm._nk * (ispin - 1)][i * dm._paraV->ncol + j];
        }
    }
    ifs.close();
}

// output density matrix dm(k) into *.dmk
template <typename TK, typename TR>
void write_DMK_file(const DensityMatrix<TK, TR>& dm,
                    const std::string& directory,
                    const int ispin,
                    const int ik)
{
    ModuleBase::TITLE("DensityMatrix", "write_DMK");
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= dm._nspin);
#endif
    // write
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
    std::ofstream ofs;
    ofs.open(fn.c_str());
    if (!ofs)
    {
        ModuleBase::WARNING("elecstate::write_dmk", "Can't create DENSITY MATRIX File!");
    }
    ofs << dm._kvec_d[ik].x << " " << dm._kvec_d[ik].y << " " << dm._kvec_d[ik].z << std::endl;
    ofs << "\n  " << dm._paraV->nrow << " " << dm._paraV->ncol << std::endl;

    ofs << std::setprecision(3);
    ofs << std::scientific;

    for (int i = 0; i < dm._paraV->nrow; ++i)
    {
        for (int j = 0; j < dm._paraV->ncol; ++j)
        {
            if (j % 8 == 0)
            {
                ofs << "\n";
            }
            ofs << " " << dm._DMK[ik + dm._nk * (ispin - 1)][i * dm._paraV->ncol + j];
        }
    }

    ofs.close();
}

template <>
void write_DMK_file<std::complex<double>, double>(
    const DensityMatrix<std::complex<double>, double>& dm,
    const std::string& directory,
    const int ispin,
    const int ik)
{
    ModuleBase::TITLE("DensityMatrix", "write_DMK");
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= dm._nspin);
#endif
    // write
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
    std::ofstream ofs;
    ofs.open(fn.c_str());
    if (!ofs)
    {
        ModuleBase::WARNING("elecstate::write_dmk", "Can't create DENSITY MATRIX File!");
    }
    ofs << dm._kvec_d[ik].x << " " << dm._kvec_d[ik].y << " " << dm._kvec_d[ik].z << std::endl;
    ofs << "\n  " << dm._paraV->nrow << " " << dm._paraV->ncol << std::endl;

    ofs << std::setprecision(3);
    ofs << std::scientific;

    for (int i = 0; i < dm._paraV->nrow; ++i)
    {
        for (int j = 0; j < dm._paraV->ncol; ++j)
        {
            if (j % 8 == 0)
            {
                ofs << "\n";
            }
            ofs << " " << dm._DMK[ik + dm._nk * (ispin - 1)][i * dm._paraV->ncol + j].real();
        }
    }

    ofs.close();
}

// explicit instantiation
template void read_DMK_file<double, double>(DensityMatrix<double, double>&,
                                            const std::string&,
                                            const int,
                                            const int);
template void read_DMK_file<std::complex<double>, double>(
    DensityMatrix<std::complex<double>, double>&,
    const std::string&,
    const int,
    const int);
template void read_DMK_file<std::complex<double>, std::complex<double>>(
    DensityMatrix<std::complex<double>, std::complex<double>>&,
    const std::string&,
    const int,
    const int);
template void write_DMK_file<double, double>(const DensityMatrix<double, double>&,
                                             const std::string&,
                                             const int,
                                             const int);
// write_DMK_file<std::complex<double>, double> has an explicit specialization above
template void write_DMK_file<std::complex<double>, std::complex<double>>(
    const DensityMatrix<std::complex<double>, std::complex<double>>&,
    const std::string&,
    const int,
    const int);

} // namespace module_dm
