#include "dftu_nao_adj.h"

#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base.h"

namespace DFTU_LCAO
{

std::vector<AdjacentAtomInfo> build_adjacent_atoms(const UnitCell* ucell,
                                                   Plus_U_Base* dftu,
                                                   const Grid_Driver* gridD,
                                                   const std::vector<double>& orb_cutoff,
                                                   const double onsite_radius)
{
    std::vector<AdjacentAtomInfo> adjs_all;
    adjs_all.reserve(ucell->nat);
    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0 = 0;
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!dftu->has_l_channel(T0))
        {
            continue;
        }

        AdjacentAtomInfo adjs;
        gridD->Find_atom(*ucell, tau0, T0, I0, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell->cal_dtau(iat0, iat1, R_index1).norm() * ucell->lat0
                < orb_cutoff[T1] + onsite_radius)
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        adjs_all.push_back(adjs);
    }
    return adjs_all;
}

} // namespace DFTU_LCAO
