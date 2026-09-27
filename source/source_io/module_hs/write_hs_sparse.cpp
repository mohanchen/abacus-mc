#include "write_hs_sparse.h"

#include "source_base/global_function.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_lcao/module_rt/td_info.h"
#include "single_r_io.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <vector>

namespace
{
template <typename Tdata>
std::vector<long long> count_nonzeros_by_R(
    const ModuleIO::SparseRMatrix<Tdata>& smat,
    const std::set<ModuleIO::RCoordinate>& all_R_coor,
    const double threshold,
    const bool reduce)
{
    std::vector<long long> nonzero_num(all_R_coor.size(), 0);
    int count = 0;
    for (const auto& R_coor: all_R_coor)
    {
        const auto iter = smat.find(R_coor);
        if (iter != smat.end())
        {
            for (const auto& row_loop: iter->second)
            {
                for (const auto& col_value: row_loop.second)
                {
                    if (std::abs(col_value.second) > threshold)
                    {
                        ++nonzero_num[count];
                    }
                }
            }
        }
        ++count;
    }

    if (reduce)
    {
        Parallel_Reduce::reduce_all(nonzero_num.data(), static_cast<int>(nonzero_num.size()));
    }
    return nonzero_num;
}

int count_output_R(const std::vector<long long>& nonzero_num)
{
    int output_R_number = 0;
    for (const long long count: nonzero_num)
    {
        if (count != 0)
        {
            ++output_R_number;
        }
    }
    return output_R_number;
}

void open_sparse_file(std::ofstream& ofs, const ModuleIO::SparseWriteOptions& options)
{
    std::ios_base::openmode mode = std::ios::out;
    if (options.binary)
    {
        mode |= std::ios::binary;
    }
    const bool append_on_restart = (options.calculation == "md") && options.out_app_flag && options.istep;
    if (append_on_restart)
    {
        mode |= std::ios::app;
    }
    ofs.open(options.filename.c_str(), mode);
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::open_sparse_file",
                                 "Cannot open sparse matrix file: " + options.filename);
    }
}

void write_sparse_header(std::ofstream& ofs,
                         const ModuleIO::SparseWriteOptions& options,
                         const int nlocal,
                         const int output_R_number)
{
    const int step = std::max(options.istep, 0);
    if (options.binary)
    {
        ofs.write(reinterpret_cast<const char*>(&step), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&nlocal), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&output_R_number), sizeof(int));
    }
    else
    {
        ofs << "STEP: " << step << std::endl;
        ofs << "Matrix Dimension of " + options.label + "(R): " << nlocal
            << std::endl;
        ofs << "Matrix number of " + options.label + "(R): "
            << output_R_number << std::endl;
    }
}

void write_R_record(std::ofstream& ofs,
                    const ModuleIO::RCoordinate& R_coor,
                    const long long nonzero_count,
                    const bool binary)
{
    int dRx = R_coor.x;
    int dRy = R_coor.y;
    int dRz = R_coor.z;
    if (binary)
    {
        const int count = static_cast<int>(nonzero_count);
        ofs.write(reinterpret_cast<char*>(&dRx), sizeof(int));
        ofs.write(reinterpret_cast<char*>(&dRy), sizeof(int));
        ofs.write(reinterpret_cast<char*>(&dRz), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&count), sizeof(int));
    }
    else
    {
        ofs << dRx << " " << dRy << " " << dRz << " " << nonzero_count
            << std::endl;
    }
}

void check_output_file_open(const std::ofstream& ofs,
                            const std::string& filename,
                            const std::string& context)
{
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT(context, "Cannot open sparse matrix file: " + filename);
    }
}
} // namespace

void ModuleIO::save_dH_sparse(const int& istep,
                              const Parallel_Orbitals& pv,
                              LCAO_HS_Arrays& HS_Arrays,
                              const double& sparse_thr,
                              const bool& binary,
                              const std::string& fileflag,
                              const int precision,
                              const std::string& global_out_dir,
                              const std::string& global_matrix_dir,
                              const std::string& calculation,
                              const bool out_app_flag,
                              const int nspin,
                              const int nlocal) {
    ModuleBase::TITLE("ModuleIO", "save_dH_sparse");
    ModuleBase::timer::start("ModuleIO", "save_dH_sparse");
    SparseWriteOptions single_R_options;
    single_R_options.threshold = sparse_thr;
    single_R_options.binary = binary;
    single_R_options.precision = precision;
    single_R_options.reduce = true;
    single_R_options.calculation = calculation;
    single_R_options.out_app_flag = out_app_flag;

    auto& all_R_coor_ptr = HS_Arrays.all_R_coor;
    auto& output_R_coor_ptr = HS_Arrays.output_R_coor;

    // The three Cartesian derivative components (x, y, z) share identical
    // file handling; only the label and the sparse matrices differ.
    struct Component
    {
        char axis;                                              // 'x' / 'y' / 'z'
        SparseRMatrix<double>* sparse;                          // nspin != 4 (array of 2)
        SparseRMatrix<std::complex<double>>* soc_sparse;        // nspin == 4
        std::vector<long long> nonzero_num[2];
        std::stringstream fname[2];
        std::ofstream ofs[2];
    };
    Component comps[3] = {{'x', HS_Arrays.dHRx_sparse, &HS_Arrays.dHRx_soc_sparse, {}, {}, {}},
                          {'y', HS_Arrays.dHRy_sparse, &HS_Arrays.dHRy_soc_sparse, {}, {}, {}},
                          {'z', HS_Arrays.dHRz_sparse, &HS_Arrays.dHRz_soc_sparse, {}, {}, {}}};

    const int total_R_num = static_cast<int>(all_R_coor_ptr.size());
    int output_R_number = 0;
    int step = istep;

    int spin_loop = 1;
    if (nspin == 2) {
        spin_loop = 2;
    }

    for (auto& comp: comps)
    {
        if (nspin != 4)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                comp.nonzero_num[ispin] = count_nonzeros_by_R(comp.sparse[ispin], all_R_coor_ptr, sparse_thr, true);
            }
        }
        else
        {
            comp.nonzero_num[0] = count_nonzeros_by_R(*comp.soc_sparse, all_R_coor_ptr, sparse_thr, true);
        }
    }

    const auto has_output_R = [&](const int index) {
        for (const auto& comp: comps)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                if (comp.nonzero_num[ispin][index] != 0)
                {
                    return true;
                }
            }
        }
        return false;
    };

    for (int index = 0; index < total_R_num; ++index)
    {
        if (has_output_R(index))
        {
            output_R_number++;
        }
    }

    const bool md_no_append = (calculation == "md") && !out_app_flag;
    for (auto& comp: comps)
    {
        for (int ispin = 0; ispin < 2; ++ispin)
        {
            if (md_no_append)
            {
                comp.fname[ispin] << global_matrix_dir
                                  << "d" << fileflag << "r" << comp.axis
                                  << "s" << (ispin + 1) << "g" << step << "_nao.csr";
            }
            else
            {
                comp.fname[ispin] << global_out_dir
                                  << "d" << fileflag << "r" << comp.axis
                                  << "s" << (ispin + 1) << "_nao.csr";
            }
        }
    }

    if (GlobalV::DRANK == 0)
    {
        const bool open_in_append = (calculation == "md") && out_app_flag && step;
        for (auto& comp: comps)
        {
            const std::string label = std::string("dH") + comp.axis;
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                std::ios_base::openmode mode = std::ios::out;
                if (binary)
                {
                    mode |= std::ios::binary;
                }
                if (open_in_append)
                {
                    mode |= std::ios::app;
                }
                else if (!binary)
                {
                    GlobalV::ofs_running << " " << label << " data are in file: "
                                         << comp.fname[ispin].str() << std::endl;
                }
                comp.ofs[ispin].open(comp.fname[ispin].str().c_str(), mode);
                check_output_file_open(comp.ofs[ispin], comp.fname[ispin].str(), "ModuleIO::save_dH_sparse");

                if (binary)
                {
                    comp.ofs[ispin].write(reinterpret_cast<char*>(&step), sizeof(int));
                    comp.ofs[ispin].write(reinterpret_cast<char*>(const_cast<int*>(&nlocal)), sizeof(int));
                    comp.ofs[ispin].write(reinterpret_cast<char*>(&output_R_number), sizeof(int));
                }
                else
                {
                    comp.ofs[ispin] << "STEP: " << step << std::endl;
                    comp.ofs[ispin] << "Matrix Dimension of " << label << "(R): " << nlocal << std::endl;
                    comp.ofs[ispin] << "Matrix number of " << label << "(R): " << output_R_number << std::endl;
                }
            }
        }
    }

    output_R_coor_ptr.clear();

    int count = 0;
    for (auto& R_coor: all_R_coor_ptr) {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (!has_output_R(count))
        {
            count++;
            continue;
        }

        output_R_coor_ptr.insert(R_coor);

        if (GlobalV::DRANK == 0) {
            for (auto& comp: comps)
            {
                for (int ispin = 0; ispin < spin_loop; ++ispin) {
                    if (binary) {
                        const int comp_count = static_cast<int>(comp.nonzero_num[ispin][count]);
                        comp.ofs[ispin].write(reinterpret_cast<char*>(&dRx), sizeof(int));
                        comp.ofs[ispin].write(reinterpret_cast<char*>(&dRy), sizeof(int));
                        comp.ofs[ispin].write(reinterpret_cast<char*>(&dRz), sizeof(int));
                        comp.ofs[ispin].write(reinterpret_cast<const char*>(&comp_count), sizeof(int));
                    } else {
                        comp.ofs[ispin] << dRx << " " << dRy << " " << dRz << " "
                                        << comp.nonzero_num[ispin][count] << std::endl;
                    }
                }
            }
        }

        for (int ispin = 0; ispin < spin_loop; ++ispin) {
            for (auto& comp: comps)
            {
                if (comp.nonzero_num[ispin][count] > 0) {
                    if (nspin != 4) {
                        output_single_R(comp.ofs[ispin], comp.sparse[ispin][R_coor], pv, single_R_options);
                    } else {
                        output_single_R(comp.ofs[ispin], (*comp.soc_sparse)[R_coor], pv, single_R_options);
                    }
                }
            }
        }

        count++;
    }

    if (GlobalV::DRANK == 0) {
        for (auto& comp: comps)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin) {
                comp.ofs[ispin].close();
            }
        }
    }

    ModuleBase::timer::end("ModuleIO", "save_dH_sparse");
    return;
}

template <typename Tdata>
void ModuleIO::save_sparse(
    const SparseRMatrix<Tdata>& smat,
    const std::set<RCoordinate>& all_R_coor,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options) {
    ModuleBase::TITLE("ModuleIO", "save_sparse");
    ModuleBase::timer::start("ModuleIO", "save_sparse");
    const int nlocal = pv.get_global_row_size();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_sparse",
                                 "Parallel_Orbitals global row size must be positive.");
    }

    const std::vector<long long> nonzero_num
        = count_nonzeros_by_R(smat, all_R_coor, options.threshold, options.reduce);
    const int output_R_number = count_output_R(nonzero_num);
    std::ofstream ofs;
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        open_sparse_file(ofs, options);
        write_sparse_header(ofs, options, nlocal, output_R_number);
    }

    int count = 0;
    for (const auto& R_coor: all_R_coor)
    {
        if (nonzero_num[count] == 0)
        {
            count++;
            continue;
        }

        if (!options.reduce || GlobalV::DRANK == 0)
        {
            write_R_record(ofs, R_coor, nonzero_num[count], options.binary);
        }

        if (smat.count(R_coor))
        {
            output_single_R(ofs, smat.at(R_coor), pv, options);
        }
        else
        {
            SparseRBlock<Tdata> empty_map;
            output_single_R(ofs, empty_map, pv, options);
        }
        ++count;
    }
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        ofs.close();
    }

    ModuleBase::timer::end("ModuleIO", "save_sparse");
}

template void ModuleIO::save_sparse<double>(
    const SparseRMatrix<double>&,
    const std::set<RCoordinate>&,
    const Parallel_Orbitals&,
    const SparseWriteOptions&);

template void ModuleIO::save_sparse<std::complex<double>>(
    const SparseRMatrix<std::complex<double>>&,
    const std::set<RCoordinate>&,
    const Parallel_Orbitals&,
    const SparseWriteOptions&);
