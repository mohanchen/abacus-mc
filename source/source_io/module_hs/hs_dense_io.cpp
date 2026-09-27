#include "hs_dense_io.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/global_function.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <type_traits>

// output a square matrix
template <typename T>
void ModuleIO::save_mat(const int istep,
    const T* mat,
    const int dim,
    const bool bit,
    const int precision,
    const bool tri,
    const bool app,
    const std::string& filename,
    const Parallel_2D& pv,
    const int drank,
    const std::string& ks_solver,
    const bool reduce)
{
    ModuleBase::TITLE("ModuleIO", "save_mat");
    ModuleBase::timer::start("ModuleIO", "save_mat");

    const bool gamma_only = std::is_same<T, double>::value;

    // write .dat file
    if (bit)
    {
// write .dat file with MPI
#ifdef __MPI
        FILE* out_matrix = nullptr;

        if (drank == 0)
        {
            const char* mode = (app && istep > 0) ? "ab" : "wb";
            out_matrix = fopen(filename.c_str(), mode);
            if (out_matrix == nullptr)
            {
                ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
            }
            fwrite(&dim, sizeof(int), 1, out_matrix);
        }

        int ir=0;
        int ic=0;
        for (int i = 0; i < dim; ++i)
        {
            T* line = new T[tri ? dim - i : dim];
            ModuleBase::GlobalFunc::ZEROS(line, tri ? dim - i : dim);

            ir = pv.global2local_row(i);
            if (ir >= 0)
            {
                // data collection
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    ic = pv.global2local_col(j);
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(ks_solver))
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        line[tri ? j - i : j] = mat[iic];
                    }
                }
            }

            if (reduce)
            {
                Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);
            }

            if (drank == 0)
            {
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    fwrite(&line[tri ? j - i : j], sizeof(T), 1, out_matrix);
                }
            }
            delete[] line;

            MPI_Barrier(DIAG_WORLD);
        }

        if (drank == 0)
        {
            fclose(out_matrix);
        }
// write .dat file without MPI
#else
        const char* mode = (app && istep > 0) ? "ab" : "wb";
        FILE* out_matrix = fopen(filename.c_str(), mode);
        if (out_matrix == nullptr)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
        }

        fwrite(&dim, sizeof(int), 1, out_matrix);

        for (int i = 0; i < dim; i++)
        {
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                fwrite(&mat[i * dim + j], sizeof(T), 1, out_matrix);
            }
        }
        fclose(out_matrix);
#endif
    } // end writing .dat file
    else // write .txt file
    {
        std::ofstream out_matrix;
        out_matrix << std::scientific << std::setprecision(precision);
#ifdef __MPI
        if (drank == 0)
        {
            if (app && istep > 0)
            {
                out_matrix.open(filename.c_str(), std::ofstream::app);
            }
            else
            {
                out_matrix.open(filename.c_str());
            }
            if (!out_matrix.is_open())
            {
                ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
            }
            out_matrix << "#------------------------------------------------------------------------" << std::endl;
            out_matrix << "# ionic step " << istep+1 << std::endl; // istep starts from 0
            out_matrix << "# filename " << filename << std::endl;
            out_matrix << "# gamma only " << gamma_only << std::endl;
            out_matrix << "# rows " << dim << std::endl;
            out_matrix << "# columns " << dim << std::endl;
            out_matrix << "#------------------------------------------------------------------------" << std::endl;

        }

        int ir=0;
        int ic=0;
        for (int i = 0; i < dim; i++)
        {
            T* line = new T[tri ? dim - i : dim];
            ModuleBase::GlobalFunc::ZEROS(line, tri ? dim - i : dim);

            ir = pv.global2local_row(i);
            if (ir >= 0)
            {
                // data collection
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    ic = pv.global2local_col(j);
                    if (ic >= 0)
                    {
                        int iic=0;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(ks_solver))
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        line[tri ? j - i : j] = mat[iic];
                    }
                }
            }

            if (reduce)
            {
                Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);
            }

            if (drank == 0)
            {
                out_matrix << "Row " << i+1 << std::endl;
                size_t count = 0;
                for (int j = (tri ? i : 0); j < dim; j++)
                {
                    out_matrix << " " << line[tri ? j - i : j];
                    ++count;
                    if(count%8==0)
                    {
                        if(j!=dim-1)
                        {
                            out_matrix << std::endl;
                        }
                    }
                }
                out_matrix << std::endl;
            }
            delete[] line;
        }

        if (drank == 0)
        {
            out_matrix.close();
        }
#else
        if (app)
        {
            out_matrix.open(filename.c_str(), std::ofstream::app);
        }
        else
        {
            out_matrix.open(filename.c_str());
        }
        if (!out_matrix.is_open())
        {
            ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
        }

        out_matrix << dim;
        out_matrix << std::setprecision(precision);
        for (int i = 0; i < dim; i++)
        {
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                out_matrix << " " << mat[i * dim + j];
            }
            out_matrix << std::endl;
        }
        out_matrix.close();
#endif
    }
    ModuleBase::timer::end("ModuleIO", "save_mat");
    return;
}

// Explicit instantiations
template void ModuleIO::save_mat<double>(const int, const double*, const int, const bool,
    const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
template void ModuleIO::save_mat<std::complex<double>>(const int, const std::complex<double>*, const int,
    const bool, const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
template void ModuleIO::save_mat<float>(const int, const float*, const int, const bool,
    const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
template void ModuleIO::save_mat<std::complex<float>>(const int, const std::complex<float>*, const int,
    const bool, const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
