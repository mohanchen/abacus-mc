#include "dhs_sparse_writer.h"

#include "hs_sparse_io.h"
#include "single_r_io.h"
#include "source_base/global_function.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"

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
