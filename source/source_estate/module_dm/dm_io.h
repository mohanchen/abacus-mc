#ifndef DM_IO_H
#define DM_IO_H

#include <string>

namespace module_dm
{
template <typename TK, typename TR>
class DensityMatrix;

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
} // namespace module_dm

#endif
