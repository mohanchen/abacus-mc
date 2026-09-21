#ifndef DM_TOOLS_H
#define DM_TOOLS_H

#include <complex>
#include <map>
#include <string>
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
        DensityMatrix<TK, TR_in> &dm,
        std::vector<hamilt::HContainer<TR_out>*> &dmR_out,
        const int ik_in);

    template <typename TK, typename TR_in, typename TR_out>
    extern void cal_DMR_td(
        DensityMatrix<TK, TR_in> &dm,
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

    /// read a DMK file (SPIN<is>_<ik>.dmk) into dm's DMK block
    template <typename TK, typename TR>
    extern void read_DMK_file(DensityMatrix<TK, TR>& dm,
                              const std::string& directory,
                              const int ispin,
                              const int ik);

    /// write dm's DMK block to a DMK file (SPIN<is>_<ik>.dmk)
    template <typename TK, typename TR>
    extern void write_DMK_file(const DensityMatrix<TK, TR>& dm,
                               const std::string& directory,
                               const int ispin,
                               const int ik);
}

} // namespace elecstate

#endif
