#include "density_matrix.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/libm/libm.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_cell/klist.h"

namespace elecstate
{

//----------------------------------------------------
// density matrix class
//----------------------------------------------------

// destructor
template <typename TK, typename TR>
DensityMatrix<TK, TR>::~DensityMatrix()
{
    this->clear_DMR();
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::clear_DMR()
{
    for (hamilt::HContainer<TR>*& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    this->_dmr_ready = false;
}

template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const Parallel_Orbitals* paraV_in,
                                     const int nspin,
                                     const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                                     const int nk)
    : _paraV(paraV_in), _nspin(nspin), _kvec_d(kvec_d), _nk((nk > 0 && nk <= _kvec_d.size()) ? nk : _kvec_d.size())
{
    ModuleBase::TITLE("DensityMatrix", "resize_DMK");
    const int nks = _nk * _nspin;
    this->_DMK.resize(nks);
    for (int ik = 0; ik < nks; ik++)
    {
        this->_DMK[ik].resize(this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->_DMK.size() * this->_DMK[0].size() * sizeof(TK));
}

template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const Parallel_Orbitals* paraV_in, const int nspin)
    : _paraV(paraV_in), _nspin(nspin),
      _kvec_d({ModuleBase::Vector3<double>(0, 0, 0)}), _nk(1)
{
    ModuleBase::TITLE("DensityMatrix", "resize_gamma");
    this->_DMK.resize(_nspin);
    for (int ik = 0; ik < this->_nspin; ik++)
    {
        this->_DMK[ik].resize(this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->_DMK.size() * this->_DMK[0].size() * sizeof(TK));
}




// calculate DMR from DMK using blas for gamma-only calculation
template <>
void DensityMatrix<double, double>::cal_DMR(const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR");
    using TK = double;
    using TR = double;

    assert(ik_in == -1 || ik_in == 0);
    assert(this->_nk == 1);

    assert(this->_DMR.size()==this->_nspin && "DMR has not been initialized!");

    ModuleBase::timer::start("DensityMatrix", "cal_DMR");
    const int ld_hk = this->_paraV->nrow;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        const int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
        hamilt::HContainer<TR>*const target_DMR = this->_DMR[is - 1];
        target_DMR->set_zero();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<TR>& target_ap = target_DMR->get_atom_pair(i);
            const int iat1 = target_ap.get_atom_i();
            const int iat2 = target_ap.get_atom_j();
            const int row_ap = this->_paraV->atom_begin_row[iat1];
            const int col_ap = this->_paraV->atom_begin_col[iat2];
            const int row_size = this->_paraV->get_nrow_atom(iat1);
            const int col_size = this->_paraV->get_ncol_atom(iat2);
            const int R_size = target_ap.get_R_size();
            assert(row_ap != -1 && col_ap != -1 && "Atom-pair not belong this process");
            assert(R_size == 1);
            const ModuleBase::Vector3<int> R_index = target_ap.get_R_index(0);
            assert(R_index.x == 0 && R_index.y == 0 && R_index.z == 0);
            hamilt::BaseMatrix<TR>*const target_mat = target_ap.find_matrix(R_index);
#ifdef __DEBUG
            if (target_mat == nullptr)
            {
                std::cout << "target_mat is nullptr" << std::endl;
                continue;
            }
#endif
            // k index
            constexpr TK kphase = 1;
            // transpose DMK col=>row
            const TK* DMK_mat_ptr
                = this->_DMK[0 + ik_begin].data()
                  + col_ap * this->_paraV->nrow + row_ap;
            // set DMR element
            TR* target_DMR_ptr = target_mat->get_pointer();
            for (int mu = 0; mu < row_size; ++mu)
            {
                BlasConnector::axpy(col_size,
                                    kphase,
                                    DMK_mat_ptr,
                                    ld_hk,
                                    target_DMR_ptr,
                                    1);
                DMK_mat_ptr += 1;
                target_DMR_ptr += col_size;
            }
        }
    }
    ModuleBase::timer::end("DensityMatrix", "cal_DMR");
    this->_dmr_ready = true;
}



// switch_dmr
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::switch_dmr(const int mode)
{
    ModuleBase::TITLE("DensityMatrix", "switch_dmr");
    if (this->_nspin != 2)
    {
        return;
    }
    else
    {
        ModuleBase::timer::start("DensityMatrix", "switch_dmr");
        switch(mode)
        {
        case 0:
            // switch to original density matrix
            if (!this->dmr_tmp_.empty() && this->dmr_origin_.size() != 0)
            {
                this->_DMR[0]->allocate(this->dmr_origin_.data(), false);
                this->dmr_tmp_.clear();
            }
            // else: do nothing
            break;
        case 1:
            // switch to total magnetization density matrix, dmr_up + dmr_down
            if(this->dmr_tmp_.empty())
            {
                const size_t size = this->_DMR[0]->get_nnr();
                this->dmr_tmp_.resize(size);
                this->dmr_origin_.resize(size);
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_origin_[i] = this->_DMR[0]->get_wrapper()[i];
                    this->dmr_tmp_[i] = this->dmr_origin_[i] + this->_DMR[1]->get_wrapper()[i];
                }
                this->_DMR[0]->allocate(this->dmr_tmp_.data(), false);
            }
            else
            {
                const size_t size = this->_DMR[0]->get_nnr();
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_tmp_[i] = this->dmr_origin_[i] + this->_DMR[1]->get_wrapper()[i];
                }
            }
            break;
        case 2:
            // switch to magnetization density matrix, dmr_up - dmr_down
            if(this->dmr_tmp_.empty())
            {
                const size_t size = this->_DMR[0]->get_nnr();
                this->dmr_tmp_.resize(size);
                this->dmr_origin_.resize(size);
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_origin_[i] = this->_DMR[0]->get_wrapper()[i];
                    this->dmr_tmp_[i] = this->dmr_origin_[i] - this->_DMR[1]->get_wrapper()[i];
                }
                this->_DMR[0]->allocate(this->dmr_tmp_.data(), false);
            }
            else
            {
                const size_t size = this->_DMR[0]->get_nnr();
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_tmp_[i] = this->dmr_origin_[i] - this->_DMR[1]->get_wrapper()[i];
                }
            }
            break;
        default:
            ModuleBase::WARNING_QUIT("density_matrix.cpp", "Unknown mode in switch_dmr");
        }
        ModuleBase::timer::end("DensityMatrix", "switch_dmr");
    }
}



template <>
void DensityMatrix_Tools::func_exp_mul_dmk<double>(
    const std::complex<double> kphase,
    const std::vector<std::complex<double>>& DMK_mat_trans,
    double* target_DMR_mat)
{
    const std::size_t mat_size = DMK_mat_trans.size();
    for(std::size_t i = 0; i < mat_size; i++)
    {
        target_DMR_mat[i]
            += kphase.real() * DMK_mat_trans[i].real() 
            - kphase.imag() * DMK_mat_trans[i].imag();
    }
}

template <>
void DensityMatrix_Tools::func_exp_mul_dmk<std::complex<double>>(
    const std::complex<double> kphase,
    const std::vector<std::complex<double>>& DMK_mat_trans,
    std::complex<double>* target_DMR_mat)
{
    BlasConnector::axpy(DMK_mat_trans.size(),
                        kphase,
                        DMK_mat_trans.data(),
                        1,
                        target_DMR_mat,
                        1);
}

template <>
void DensityMatrix_Tools::func_xyz_to_updown<double>(
    const std::complex<double> tmp[4],
    const int icol,
    const int step_trace[4],
    double* target_DMR_mat)
{
    target_DMR_mat[icol + step_trace[0]] = tmp[0].real() + tmp[3].real();  // rho_0 = (rho_upup + rho_downdown).real()
    target_DMR_mat[icol + step_trace[1]] = tmp[1].real() + tmp[2].real();  // rho_x = (rho_updown + rho_downup).real()
    // rho_y: the stored DM block is the complex conjugate of the physical 1-RDM P (cal_dm_psi builds
    // DM_{ab}=sum conj(c_a) c_b = conj(P), so tmp[1]=DM_{ud}=conj(P_{ud})). Extracting m_y from the
    // CONJUGATED block therefore carries the opposite sign of the bare-textbook formula; m_x/m_z read
    // Re() and are conjugation-invariant. Using the bare formula (PR #7664) sign-flips m_y and quenches
    // in-plane non-collinear moments (e.g. Mn3Sn 120-deg AFM); see issue #7831.
    target_DMR_mat[icol + step_trace[2]] = tmp[1].imag() - tmp[2].imag();  // rho_y = Im(P_updown) - Im(P_downup)
    target_DMR_mat[icol + step_trace[3]] = tmp[0].real() - tmp[3].real();  // rho_z = (rho_upup - rho_downdown).real()
}

template <>
void DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(
    const std::complex<double> tmp[4],
    const int icol,
    const int step_trace[4],
    std::complex<double>* target_DMR_mat)
{
    target_DMR_mat[icol + step_trace[0]] = tmp[0] + tmp[3];  // rho_0 = (rho_upup + rho_downdown)
    target_DMR_mat[icol + step_trace[1]] = tmp[1] + tmp[2];  // rho_x = (rho_updown + rho_downup)
    // rho_y sign accounts for the conjugated stored DM block (conj(P)); see the <double> specialization above.
    target_DMR_mat[icol + step_trace[2]]
        = -ModuleBase::IMAG_UNIT * (tmp[1] - tmp[2]);  // rho_y = -i*(rho_updown - rho_downup)
    target_DMR_mat[icol + step_trace[3]] = tmp[0] - tmp[3];  // rho_z = (rho_upup - rho_downdown)
}



// T of HContainer can be double or complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace elecstate
