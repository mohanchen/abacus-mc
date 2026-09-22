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
void DensityMatrix_Tools::cal_DMR_td(
    DensityMatrix<TK, TR_in> &dm,
    std::vector<hamilt::HContainer<TR_out>*> &dmR_out,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR_td");
    assert(dmR_out.size()==dm._nspin && "DMR has not been initialized!");

    ModuleBase::timer::start("DensityMatrix", "cal_DMR_td");
    const int ld_hk = dm.pv->nrow;
    for (int is = 1; is <= dm._nspin; ++is)
    {
        const int ik_begin = dm._nk * (is - 1); // jump dm._nk for spin_down if nspin==2
        hamilt::HContainer<TR_out>*const target_DMR = dmR_out[is - 1];
        target_DMR->set_zero();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<TR_out>& target_ap = target_DMR->get_atom_pair(i);
            const int iat1 = target_ap.get_atom_i();
            const int iat2 = target_ap.get_atom_j();
            const int row_ap = dm.pv->atom_begin_row[iat1];
            const int col_ap = dm.pv->atom_begin_col[iat2];
            const int row_size = dm.pv->get_nrow_atom(iat1);
            const int col_size = dm.pv->get_ncol_atom(iat2);
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
                    if(ik_in >= 0 && ik_in != ik)
                    {
                        continue;
                    }
                    // cal k_phase
                    const ModuleBase::Vector3<double> dR(R_index[0], R_index[1], R_index[2]);
                    const double arg = (dm._kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    kphase_vec[ik][iR] = TK(cosp, sinp);
                    if(!phase_hybrid.empty())
                    {
                        //phase for hybrid gauge tddft
                        kphase_vec[ik][iR] *= phase_hybrid.at(R_index);
                    }
                }
            }

            std::vector<TK> DMK_mat_trans(mat_size);
            std::vector<TK> tmp_DMR( (dm._nspin==4) ? mat_size*R_size : 0);
            for(int ik = 0; ik < dm._nk; ++ik)
            {
                if(ik_in >= 0 && ik_in != ik)
                {
                    continue;
                }
                const TK*const DMK_mat_ptr
                    = dm._DMK[ik + ik_begin].data()
                      + col_ap * dm.pv->nrow + row_ap;
                for(int icol = 0; icol < col_size; ++icol)
                {
                    for(int irow = 0; irow < row_size; ++irow)
                    {
                        DMK_mat_trans[irow * col_size + icol] = DMK_mat_ptr[icol * ld_hk + irow];
                    }
                }

                // if nspin != 4, fill DMR
                // if nspin == 4, fill tmp_DMR
                for(int iR = 0; iR < R_size; ++iR)
                {
                    // (kr+i*ki) * (Dr+i*Di) = (kr*Dr-ki*Di) + i*(kr*Di+ki*Dr)
                    const TK kphase = kphase_vec[ik][iR];
                    if(dm._nspin != 4)                // only save real kr*Dr-ki*Di
                {
                    func_exp_mul_dmk(kphase, DMK_mat_trans, target_DMR_mat_vec[iR]);
                }
                else if(dm._nspin == 4)
                {
                    BlasConnector::axpy(mat_size,
                                        kphase,
                                        DMK_mat_trans.data(),
                                        1,
                                        &tmp_DMR[iR * mat_size],
                                        1);
                }
                }
            }

            // if nspin == 4
            // copy tmp_DMR to fill target_DMR
            if(dm._nspin == 4)
            {
                int step_trace[4]{};
                constexpr int npol = 2;
                for (int is = 0; is < npol; is++)
                {
                    for (int is2 = 0; is2 < npol; is2++)
                    {
                        step_trace[is * npol + is2] = target_ap.get_col_size() * is + is2;
                    }
                }

                TK tmp[4]{};
                for(int iR = 0; iR < R_size; ++iR)
                {
                    const TK* tmp_DMR_mat = &tmp_DMR[iR * mat_size];
                    TR_out* target_DMR_mat = target_DMR_mat_vec[iR];
                    for (int irow = 0; irow < row_size; irow += 2)
                    {
                        for (int icol = 0; icol < col_size; icol += 2)
                        {
                            tmp[0] = tmp_DMR_mat[icol + step_trace[0]];
                            tmp[1] = tmp_DMR_mat[icol + step_trace[1]];
                            tmp[2] = tmp_DMR_mat[icol + step_trace[2]];
                            tmp[3] = tmp_DMR_mat[icol + step_trace[3]];

                            func_xyz_to_updown(tmp, icol, step_trace, target_DMR_mat);
                        }
                        tmp_DMR_mat += col_size * 2;
                        target_DMR_mat += col_size * 2;
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("DensityMatrix", "cal_DMR_td");
    dm._dmr_ready = true;
}

template <>
void DensityMatrix<double, double>::cal_DMR_td(
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    return;
}
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR_td(
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    DensityMatrix_Tools::cal_DMR_td(*this, this->_DMR, phase_hybrid, At, ik_in);
}

template <>
void DensityMatrix<std::complex<double>, std::complex<double>>::cal_DMR_td(
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    DensityMatrix_Tools::cal_DMR_td(*this, this->_DMR, phase_hybrid, At, ik_in);
}

} // namespace module_dm
