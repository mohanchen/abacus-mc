#ifndef DM_TOOLS_H
#define DM_TOOLS_H

#include <complex>
#include <map>
#include <vector>

#include "source_base/vector3.h"

namespace hamilt
{
template <typename T>
class HContainer;
}

namespace elecstate
{
template <typename TK, typename TR>
class DensityMatrix;

// DensityMatrix<complex<double>,TR>::cal_DMR() is illegal in C++, so DensityMatrix_Tools is used instead.
namespace DensityMatrix_Tools
{
    template <typename TK, typename TR_in, typename TR_out>
    extern void cal_DMR(
        const DensityMatrix<TK, TR_in> &dm,
        std::vector<hamilt::HContainer<TR_out>*> &dmR_out,
        const int ik_in);

    template <typename TK, typename TR_in, typename TR_out>
    extern void cal_DMR_td(
        const DensityMatrix<TK, TR_in> &dm,
        std::vector<hamilt::HContainer<TR_out>*> &dmR_out,
        const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
        const ModuleBase::Vector3<double> At,
        const int ik_in);

    template <typename TK, typename TR_in, typename TR_out>
    extern void cal_DMR_full(
        const DensityMatrix<TK, TR_in> &dm,
        hamilt::HContainer<TR_out>* dmR_out,
        const int ik_in);

    template <typename TR>
    extern void func_exp_mul_dmk(const std::complex<double> kphase,
                                const std::vector<std::complex<double>>& DMK_mat_trans,
                                TR* target_DMR_mat);

    template <typename TR>
    extern void func_xyz_to_updown(const std::complex<double> tmp[4],
                                  const int icol,
                                  const int step_trace[4],
                                  TR* target_DMR_mat);
}

} // namespace elecstate

#endif
