/**
 * @file klist_io.cpp
 * @brief this-free helpers extracted from K_Vectors (IBZ table formatting and
 *        line-mode k-point interpolation). Kept separate from klist.cpp so the
 *        logic is testable in isolation; klist.cpp only keeps thin wrappers.
 */
#include "klist_io.h"

#include "source_base/formatter.h"
#include "source_base/global_function.h"
#include "source_cell/reciprocal_grid.h"

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

bool find_kpoints_header(std::ifstream& ifk)
{
    std::string word;
    while (ifk.good())
    {
        ifk >> word;
        // LiuXh add 20180416, fix bug in k-point file when the first line with comments
        ifk.ignore(150, '\n');
        if (word == "K_POINTS" || word == "KPOINTS" || word == "K")
        {
            return true;
        }
    }
    return false;
}

void read_kpt_list(std::ifstream& ifk,
                   const int nkstot,
                   std::vector<ModuleBase::Vector3<double>>& kvec,
                   std::vector<double>& wk)
{
    for (int i = 0; i < nkstot; i++)
    {
        ifk >> kvec[i].x >> kvec[i].y >> kvec[i].z;
        ModuleBase::GlobalFunc::READ_VALUE(ifk, wk[i]);
    }
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

void build_kstars(const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                  const std::vector<ModuleBase::Matrix3>& kgmatrix,
                  const int nrotkm,
                  const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                  const double epsilon,
                  const std::function<bool(double, double)>& equal,
                  std::vector<std::map<int, ModuleBase::Vector3<double>>>& kstars)
{
    const int nkstot = static_cast<int>(kvec_d.size());
    const int nkstot_ibz = static_cast<int>(kvec_d_ibz.size());
    kstars.resize(nkstot_ibz);

    ModuleBase::Vector3<double> kvec_rot;
    for (int i = 0; i < nkstot; ++i)
    {
        int exist_number = -1;
        int isym = 0;
        for (int j = 0; j < nrotkm; ++j)
        {
            kvec_rot = kvec_d[i] * kgmatrix[j];
            ModuleCell::restrict_kpt(kvec_rot, epsilon);
            for (int k = 0; k < nkstot_ibz; ++k)
            {
                if (equal(kvec_rot.x, kvec_d_ibz[k].x) && equal(kvec_rot.y, kvec_d_ibz[k].y)
                    && equal(kvec_rot.z, kvec_d_ibz[k].z))
                {
                    isym = j;
                    exist_number = k;
                    break;
                }
            }
            if (exist_number != -1)
            {
                break;
            }
        }
        kstars[exist_number].insert(std::make_pair(isym, kvec_d[i]));
    }
}

void pack_kpts(const std::vector<int>& isk,
               const std::vector<double>& wk,
               const std::vector<ModuleBase::Vector3<double>>& kvec_c,
               const std::vector<ModuleBase::Vector3<double>>& kvec_d,
               const std::vector<ModuleBase::Vector3<double>>& kvec_c_full,
               const int nkstot,
               std::vector<int>& isk_aux,
               std::vector<double>& wk_aux,
               std::vector<double>& kvec_c_aux,
               std::vector<double>& kvec_d_aux,
               std::vector<double>& kvec_c_full_aux)
{
    for (int ik = 0; ik < nkstot; ik++)
    {
        isk_aux[ik] = isk[ik];
        wk_aux[ik] = wk[ik];
        kvec_c_aux[3 * ik] = kvec_c[ik].x;
        kvec_c_aux[3 * ik + 1] = kvec_c[ik].y;
        kvec_c_aux[3 * ik + 2] = kvec_c[ik].z;
        kvec_d_aux[3 * ik] = kvec_d[ik].x;
        kvec_d_aux[3 * ik + 1] = kvec_d[ik].y;
        kvec_d_aux[3 * ik + 2] = kvec_d[ik].z;
    }
    const int nkstot_full = static_cast<int>(kvec_c_full.size());
    for (int ik = 0; ik < nkstot_full; ik++)
    {
        kvec_c_full_aux[3 * ik] = kvec_c_full[ik].x;
        kvec_c_full_aux[3 * ik + 1] = kvec_c_full[ik].y;
        kvec_c_full_aux[3 * ik + 2] = kvec_c_full[ik].z;
    }
}

void unpack_kpts(const std::vector<int>& isk_aux,
                 const std::vector<double>& wk_aux,
                 const std::vector<double>& kvec_c_aux,
                 const std::vector<double>& kvec_d_aux,
                 const std::vector<double>& kvec_c_full_aux,
                 const int nks,
                 const int startk,
                 std::vector<int>& isk,
                 std::vector<double>& wk,
                 std::vector<ModuleBase::Vector3<double>>& kvec_c,
                 std::vector<ModuleBase::Vector3<double>>& kvec_d,
                 std::vector<ModuleBase::Vector3<double>>& kvec_c_full)
{
    for (int i = 0; i < nks; i++)
    {
        // 3 is because each k point has three value:kx, ky, kz
        const int k_index = i + startk;
        kvec_c[i].x = kvec_c_aux[k_index * 3];
        kvec_c[i].y = kvec_c_aux[k_index * 3 + 1];
        kvec_c[i].z = kvec_c_aux[k_index * 3 + 2];
        kvec_d[i].x = kvec_d_aux[k_index * 3];
        kvec_d[i].y = kvec_d_aux[k_index * 3 + 1];
        kvec_d[i].z = kvec_d_aux[k_index * 3 + 2];
        kvec_c_full[i].x = kvec_c_full_aux[k_index * 3];
        kvec_c_full[i].y = kvec_c_full_aux[k_index * 3 + 1];
        kvec_c_full[i].z = kvec_c_full_aux[k_index * 3 + 2];
        wk[i] = wk_aux[k_index];
        isk[i] = isk_aux[k_index];
    }
}

} // namespace KListIO
