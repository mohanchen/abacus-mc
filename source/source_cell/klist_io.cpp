/**
 * @file klist_io.cpp
 * @brief this-free helpers extracted from K_Vectors (IBZ table formatting and
 *        line-mode k-point interpolation). Kept separate from klist.cpp so the
 *        logic is testable in isolation; klist.cpp only keeps thin wrappers.
 */
#include "klist_io.h"

#include "source_base/formatter.h"
#include "source_base/global_function.h"

#include <sstream>

namespace KListIO
{

std::string ibz_kpt_table(const int nkstot,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                          const std::vector<int>& ibz_index,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz)
{
    std::stringstream ss;
    ss << " " << std::setw(40) << "nkstot"
       << " = " << nkstot << std::setw(66) << "ibzkpt" << std::endl;
    std::string table;
    table += "K-POINTS REDUCTION ACCORDING TO SYMMETRY\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s%12s%12s%12s\n",
                             "KPT",
                             "DIRECT_X",
                             "DIRECT_Y",
                             "DIRECT_Z",
                             "IBZ",
                             "DIRECT_X",
                             "DIRECT_Y",
                             "DIRECT_Z");
    for (int i = 0; i < nkstot; ++i)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8d%12.8f%12.8f%12.8f\n",
                                 i + 1,
                                 kvec_d[i].x,
                                 kvec_d[i].y,
                                 kvec_d[i].z,
                                 ibz_index[i] + 1,
                                 kvec_d_ibz[ibz_index[i]].x,
                                 kvec_d_ibz[ibz_index[i]].y,
                                 kvec_d_ibz[ibz_index[i]].z);
    }
    ss << table << std::endl;
    return ss.str();
}

std::string ibz_wk_table(const int nkstot_ibz,
                         const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                         const std::vector<double>& wk_ibz,
                         const std::vector<int>& ibz2bz)
{
    std::string table;
    table += "\n K-POINTS REDUCTION ACCORDING TO SYMMETRY\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s%8s\n", "IBZ", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT", "ibz2bz");
    for (int ik = 0; ik < nkstot_ibz; ik++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f%8d\n",
                                 ik + 1,
                                 kvec_d_ibz[ik].x,
                                 kvec_d_ibz[ik].y,
                                 kvec_d_ibz[ik].z,
                                 wk_ibz[ik],
                                 ibz2bz[ik]);
    }
    return table;
}

LineK interp_line(std::ifstream& ifk, const int nks_special)
{
    // number of points to the next k points
    std::vector<int> nkl(nks_special, 0);

    // coordinates of special points.
    std::vector<ModuleBase::Vector3<double>> ks(nks_special);

    LineK out;
    std::vector<int> kpt_segids;
    int kpt_segid = 0;
    for (int iks = 0; iks < nks_special; iks++)
    {
        ifk >> ks[iks].x;
        ifk >> ks[iks].y;
        ifk >> ks[iks].z;
        ModuleBase::GlobalFunc::READ_VALUE(ifk, nkl[iks]);

        if (nkl[iks] <= 0)
        {
            ModuleBase::WARNING_QUIT("KListIO::interp_line",
                                     "Line-mode interpolation counts must be positive.");
        }
        out.nks_total += nkl[iks];
        /* ISSUE#3482: to distinguish different kline segments */
        if ((nkl[iks] == 1) && (iks != (nks_special - 1))) {
            kpt_segid++;
        }
        kpt_segids.push_back(kpt_segid);
    }
    if (nkl[nks_special - 1] != 1)
    {
        ModuleBase::WARNING_QUIT("KListIO::interp_line",
                                 "The final line-mode k-point must have an interpolation count of 1.");
    }

    out.kpts.resize(out.nks_total);
    out.segids.reserve(out.nks_total);

    int count = 0;
    for (int iks = 1; iks < nks_special; iks++)
    {
        double dxs = (ks[iks].x - ks[iks - 1].x) / nkl[iks - 1];
        double dys = (ks[iks].y - ks[iks - 1].y) / nkl[iks - 1];
        double dzs = (ks[iks].z - ks[iks - 1].z) / nkl[iks - 1];
        for (int is = 0; is < nkl[iks - 1]; is++)
        {
            out.kpts[count].x = ks[iks - 1].x + is * dxs;
            out.kpts[count].y = ks[iks - 1].y + is * dys;
            out.kpts[count].z = ks[iks - 1].z + is * dzs;
            out.segids.push_back(kpt_segids[iks - 1]); /* ISSUE#3482 */
            ++count;
        }
    }

    // deal with the last special k point.
    out.kpts[count].x = ks[nks_special - 1].x;
    out.kpts[count].y = ks[nks_special - 1].y;
    out.kpts[count].z = ks[nks_special - 1].z;
    out.segids.push_back(kpt_segids[nks_special - 1]); /* ISSUE#3482 */
    ++count;

    assert(count == out.nks_total);
    assert(out.segids.size() == static_cast<size_t>(out.nks_total)); /* ISSUE#3482 */
    return out;
}

} // namespace KListIO
