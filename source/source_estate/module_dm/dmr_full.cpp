#include "density_matrix.h"

#include "source_base/libm/libm.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"

namespace module_dm
{

template <typename TK, typename TR_in, typename TR_out>
void DensityMatrix_Tools::cal_DMR_full(
    const DensityMatrix<TK, TR_in> &dm,
    hamilt::HContainer<TR_out>* dmR_out,
    const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR_full");

    ModuleBase::timer::start("DensityMatrix", "cal_DMR_full");
    const int ld_hk = dm._paraV->nrow;
    hamilt::HContainer<TR_out>* target_DMR = dmR_out;
    target_DMR->set_zero();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR_out>& target_ap = target_DMR->get_atom_pair(i);
        const int iat1 = target_ap.get_atom_i();
        const int iat2 = target_ap.get_atom_j();
        const int row_ap = dm._paraV->atom_begin_row[iat1];
        const int col_ap = dm._paraV->atom_begin_col[iat2];
        const int row_size = dm._paraV->get_nrow_atom(iat1);
        const int col_size = dm._paraV->get_ncol_atom(iat2);
        const int mat_size = row_size * col_size;
        const int R_size = target_ap.get_R_size();
        assert(row_ap != -1 && col_ap != -1 && "Atom-pair not belong this process");

        // calculate kphase and target_mat_ptr
        std::vector<std::vector<TK>> kphase_vec(dm._nk, std::vector<TK>(R_size));
        std::vector<TR_out*> target_DMR_mat_vec(R_size);
        for(int iR = 0; iR < R_size; ++iR)
        {
            const ModuleBase::Vector3<int> R_index = target_ap.get_R_index(iR);
            hamilt::BaseMatrix<TR_out>*const target_mat = target_ap.find_matrix(R_index);
#ifdef __DEBUG
            if (target_mat == nullptr)
            {
                std::cout << "target_mat is nullptr" << std::endl;
                continue;
            }
#endif
            target_DMR_mat_vec[iR] = target_mat->get_pointer();
            for(int ik = 0; ik < dm._nk; ++ik)
            {
                if(ik_in >= 0 && ik_in != ik) { continue; }
                // cal k_phase
                const ModuleBase::Vector3<double> dR(R_index[0], R_index[1], R_index[2]);
                const double arg = (dm._kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                kphase_vec[ik][iR] = TK(cosp, sinp);
            }
        }

        std::vector<TK> DMK_mat_trans(mat_size);
        for(int ik = 0; ik < dm._nk; ++ik)
        {
            if(ik_in >= 0 && ik_in != ik) { continue; }
            const TK*const DMK_mat_ptr
                = dm._DMK[ik].data()
                  + col_ap * dm._paraV->nrow + row_ap;
            for(int icol = 0; icol < col_size; ++icol) {
                for(int irow = 0; irow < row_size; ++irow) {
                    DMK_mat_trans[irow * col_size + icol] = DMK_mat_ptr[icol * ld_hk + irow];
            }}

            for(int iR = 0; iR < R_size; ++iR)
            {
                const TK kphase = kphase_vec[ik][iR];
                BlasConnector::axpy(mat_size,
                                    kphase,
                                    DMK_mat_trans.data(),
                                    1,
                                    target_DMR_mat_vec[iR],
                                    1);
            }
        }
    }
    ModuleBase::timer::end("DensityMatrix", "cal_DMR_full");
}

template <>
void DensityMatrix<double, double>::cal_DMR_full(
    hamilt::HContainer<std::complex<double>>* dmR_out,
    const int ik_in) const{}
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR_full(
    hamilt::HContainer<std::complex<double>>* dmR_out,
    const int ik_in) const
{
    DensityMatrix_Tools::cal_DMR_full(*this, dmR_out, ik_in);
}

} // namespace module_dm
