#ifndef KLIST_IO_H
#define KLIST_IO_H

#include "source_base/vector3.h"

#include <fstream>
#include <string>
#include <vector>

/// this-free helpers extracted from K_Vectors, kept in a separate TU so they
/// can be unit-tested and reused without dragging in the K_Vectors class.
namespace KListIO
{
/// Render the IBZ reduction table ("IBZ" k-point -> originating k-point).
std::string ibz_kpt_table(int nkstot,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                          const std::vector<int>& ibz_index,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz);

/// Render the IBZ weight table (IBZ k-point, weight, multiplicity, origin index).
std::string ibz_wk_table(int nkstot_ibz,
                         const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                         const std::vector<double>& wk_ibz,
                         const std::vector<int>& ibz2bz);

/// Result of line-mode interpolation between special k-points.
struct LineK
{
    std::vector<ModuleBase::Vector3<double>> kpts; ///< interpolated k points
    std::vector<int> segids;                       ///< segment id per k point (ISSUE#3482)
    int nks_total = 0;                             ///< total interpolated k-point count
};

/// Read the special k points and per-point interpolation counts from `ifk`,
/// then linearly interpolate the line-mode k points. Pure function of the
/// stream and `nks_special`; dies via WARNING_QUIT on malformed input.
LineK interp_line(std::ifstream& ifk, int nks_special);
} // namespace KListIO

#endif // KLIST_IO_H
