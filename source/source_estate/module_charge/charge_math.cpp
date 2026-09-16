#include "charge_math.h"

#include "source_base/global_function.h"
#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_threading.h"
#include "source_base/tool_title.h"

#include <cassert>
#include <cmath>
#include <new>

namespace charge_math
{

double sum_rho(double* const* rho,
               const int nspin0,
               const int nrxx,
               const double omega,
               const int nxyz)
{
    ModuleBase::TITLE("charge_math", "sum_rho");

    double sum_rho = 0.0;

    for (int is = 0; is < nspin0; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            sum_rho += rho[is][ir];
        }
    }

    // multiply the sum of charge density by a factor
    sum_rho *= omega / static_cast<double>(nxyz);

#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum_rho);
#endif

    // mohan fixed bug 2010-01-18,
    // sum_rho may be smaller than 1, like Na bcc.
    if (sum_rho <= 0.1)
    {
        ModuleBase::WARNING_QUIT("charge_math::sum_rho", "Can't find even an electron!");
    }

    return sum_rho;
}

double cal_rho2ne(const double* rho_in,
                  const int nrxx,
                  const double omega,
                  const int nxyz)
{
    assert(nxyz > 0); // mohan add 2025-12-02

    double ne = 0.0;
    for (int ir = 0; ir < nrxx; ir++)
    {
        ne += rho_in[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(ne);
#endif
    ne = ne * omega / static_cast<double>(nxyz);

    return ne;
}

void non_linear_core_correction(const bool numeric,
                                const double omega,
                                const double tpiba2,
                                const int mesh,
                                const double* r,
                                const double* rab,
                                const double* rhoc,
                                double* rhocg,
                                const double* gg_uniq,
                                const int ngg)
{
    ModuleBase::TITLE("charge_math", "drhoc");

    // use labmda instead of repeating codes
    const auto kernel = [&](int num_threads, int thread_id)
    {

    double gx = 0.0;
    double rhocg1 = 0.0;
    double *aux = nullptr;

    // here we compute the fourier transform is the charge in numeric form
    if (numeric)
    {
        aux = new double [mesh];
        // G=0 term

        int igl0 = 0;
        if (gg_uniq [0] < 1.0e-8)
        {
            // single thread term
            if (thread_id == 0)
            {
                for (int ir = 0;ir < mesh; ir++)
                {
                    aux [ir] = r [ir] * r [ir] * rhoc [ir];
                }
                ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
                //rhocg [1] = fpi * rhocg1 / omega;
                rhocg [0] = ModuleBase::FOUR_PI * rhocg1 / omega;//mohan modify 2008-01-19
            }
            igl0 = 1;
        }

        int igl_beg, igl_end;
        // exclude igl0
        ModuleBase::TASK_DIST_1D(num_threads, thread_id, ngg - igl0, igl_beg, igl_end);
        igl_beg += igl0;
        igl_end += igl_beg;

        // G <> 0 term
        for (int igl = igl_beg; igl < igl_end;igl++)
        {
            gx = sqrt(gg_uniq[igl] * tpiba2);
            ModuleBase::Sphbes::Spherical_Bessel(mesh, r, gx, 0, aux);
            for (int ir = 0;ir < mesh; ir++)
            {
                aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
            } //  enddo
            ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
            rhocg [igl] = ModuleBase::FOUR_PI * rhocg1 / omega;
        } //  enddo
        delete [] aux;
    }
    else
    {
        // here the case where the charge is in analytic form,
        // check old version before 2008-12-9
    }

    }; // end kernel

    // do not use omp parallel when this function is already in parallel block
    //
    // it is called in parallel block in Forces::cal_force_cc,
    // but not in other funtcion such as Stress_Func::stress_cc.
    ModuleBase::TRY_OMP_PARALLEL(kernel);

    return;
}

} // namespace charge_math
