#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include <complex>
#include <map>
#include <string>
#include <vector>

#include "dm_shift.h"
#include "source_base/vector3.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/record_adj.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

namespace module_dm
{
/**
 * @brief DensityMatrix Class
 * <TK,TR> = <double,double> for Gamma-only calculation
 * <TK,TR> = <std::complex<double>,double> for multi-k calculation
 */
    template <typename TK, typename TR> class DensityMatrix;

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

}

} // namespace module_dm

namespace module_dm
{

template <typename TK, typename TR>
class DensityMatrix
{
    using TRShift = typename ShiftRealComplex<TR>::type;

    public:
    /**
     * @brief Destructor of class DensityMatrix
     */
    ~DensityMatrix();

    /**
     * @brief Constructor of class DensityMatrix for multi-k calculation
     * @param pv pointer of Parallel_Orbitals object
     * @param spin_mult spin multiplicity used to size the DM: 1 for input nspin 1 or 4
     *  (non-collinear k points are not doubled), 2 for input nspin 2 (LSDA up/down).
     *  This is NOT the physical nspin (1/2/4); it matches K_Vectors::spin_mult.
     * @param kvec_d direct coordinates of kpoints
     * @param nk number of k-points, not always equal to K_Vectors::get_nks()/spin_mult.
     *               it will be set to kvec_d.size() if the value is invalid
     * @param nspin the global physical nspin from INPUT (1/2/4); defaults to spin_mult for
     *               non-SOC cases where they coincide. Pass 4 explicitly for SOC/noncollinear
     *               calculations so that cal_DMR selects the spin-resolved (Pauli) branch.
     */
    DensityMatrix(const Parallel_Orbitals* pv,
            const int spin_mult,
            const std::vector<ModuleBase::Vector3<double>>& kvec_d,
            const int nk,
            const int nspin = 0);

    /**
     * @brief Constructor of class DensityMatrix for gamma-only calculation, where kvector is not required
     * @param pv pointer of Parallel_Orbitals object
     * @param spin_mult spin multiplicity of the density matrix (1 or 2); NOT the physical nspin
     * @param nspin the global physical nspin from INPUT (1/2/4); defaults to spin_mult.
     */
    DensityMatrix(const Parallel_Orbitals* pv, const int spin_mult, const int nspin = 0);

    /**
     * @brief initialize density matrix DMR from UnitCell
     * @param GridD_in pointer of Grid_Driver object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_DMR(const Grid_Driver* GridD_in, const UnitCell* ucell);

    /**
     * @brief initialize density matrix DMR from UnitCell and RA
     * @param ra pointer of Record_adj object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_DMR(Record_adj& ra, const UnitCell* ucell);

    /**
     * @brief initialize density matrix DMR from another HContainer
     * now only support HContainer<double>
     * @param _DMR_in pointer of another HContainer object
     */
    void init_DMR(const hamilt::HContainer<TR>& _DMR_in);

    /// @brief initialize density matrix DMR from another HContainer
    /// this is a temprory function for NSPIN=4 case 
    /// since copy HContainer from another HContainer with different TR is not supported yet
    /// would be refactor in the future
    /// @param _DMR_in 
    // the old input type ``:HContainer<complex<double>` causes redefination error if TR = complex<double>
    void init_DMR(const hamilt::HContainer<TRShift>& _DMR_in);

    /**
     * @brief set _DMK element directly
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @param value value to be set
     */
    void set_DMK(const int ispin, const int ik, const int i, const int j, const TK value);

    /**
     * @brief set _DMK element to zero
    */
    void set_DMK_zero();
    
    /**
     * @brief get a matrix element of density matrix dm(k)
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @return T a matrix element of density matrix dm(k)
     */
    TK get_DMK(const int ispin, const int ik, const int i, const int j) const;

    /**
     * @brief get total number of k-points of density matrix dm(k)
     */
    int get_DMK_nks() const;
    int get_DMK_size() const;

    /**
     * @brief get number of rows of density matrix dm(k)
     */
    int get_DMK_nrow() const;

    /**
     * @brief get number of columns of density matrix dm(k)
     */
    int get_DMK_ncol() const;

    /**
     * @brief get pointer of DMR
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @return HContainer<TR>* pointer of DMR
     */
    hamilt::HContainer<TR>* get_DMR_pointer(const int ispin) const;

    /**
     * @brief check whether the stored DMR is a valid density matrix calculated from DMK
     * init_DMR() resets the flag and cal_DMR()/cal_DMR_td() set it, so a freshly
     * allocated, zeroed or file-read DMR is reported as not ready until the first
     * wavefunction-derived calculation
     * @return true if DMR is ready for Hamiltonian construction
     */
    bool is_dmr_ready() const
    {
        return this->_dmr_ready;
    }

    /**
     * @brief get pointer vector of DMR
     * @return HContainer<TR>* vector of DMR
     */
    const std::vector<hamilt::HContainer<TR>*>& get_DMR_vector() const
    {
        return this->_DMR;
    }
    std::vector<hamilt::HContainer<TR>*>& get_DMR_vector()
    {
        return this->_DMR;
    }

    const std::vector<std::vector<TR>>& get_DMR_save() const
    {
        return this->_DMR_save;
    }
    std::vector<std::vector<TR>>& get_DMR_save()
    {
        return this->_DMR_save;
    }

    /**
     * @brief get pointer of DMK
     * @param ik k-point index, which is the index of _DMK
     * @return TK* pointer of DMK
     */
    TK* get_DMK_pointer(const int ik) const;

    /**
     * @brief get pointer vector of DMK
    */
    const std::vector<std::vector<TK>>& get_DMK_vector() const
    {
        return this->_DMK;
    }
    std::vector<std::vector<TK>>& get_DMK_vector()
    {
        return this->_DMK;
    }

    /**
     * @brief set _DMK using a input TK* pointer
     * please make sure the size of TK* is correct
    */
    void set_DMK_pointer(const int ik, TK* DMK_in);

    /**
     * @brief get pointer of paraV
     */
    const Parallel_Orbitals* get_paraV_pointer() const
    {
        return this->pv;
    }

    const std::vector<ModuleBase::Vector3<double>>& get_kvec_d() const
    {
        return this->_kvec_d;
    }

    /**
     * @brief calculate density matrix DMR from dm(k) using blas::axpy
     * @param ik_in
     * if ik_in < 0, calculate all k-points
     * if ik_in >= 0, calculate only one k-point without summing over k-points
     */
    void cal_DMR(const int ik_in);

    /**
     * @brief calculate density matrix DMR with additional vector potential phase, used for hybrid gauge tddft
     * @param ik_in
     * if ik_in < 0, calculate all k-points
     * if ik_in >= 0, calculate only one k-point
     */
    void cal_DMR_td(const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
                    const ModuleBase::Vector3<double> At,
                    const int ik_in);

    /**
     * @brief calculate complex density matrix DMR with both real and imaginary part for noncollinear-spin calculation
     * the stored dm(k) has been used to calculate the passin DMR
     * @param dmR_out pointer of HContainer object to store the calculated complex DMR
     * @param ik_in
     * if ik_in < 0, calculate all k-points
     * if ik_in >= 0, calculate only one k-point
     */
    void cal_DMR_full(hamilt::HContainer<std::complex<double>>* dmR_out, const int ik_in) const;

    /**
     * @brief (Only nspin=2) switch DMR to total density matrix or magnetization density matrix
     * @param mode 0 - original density matrix; 1 - total density matrix; 2 - magnetization density matrix
     */
    void switch_dmr(const int mode);

    /**
     * @brief save _DMR into _DMR_save
     */
    void save_DMR();
    
    std::vector<ModuleBase::ComplexMatrix> EDMK; // for TD-DFT

#ifdef __PEXSI
    /**
     * @brief EDM storage for PEXSI
     * used in MD calculation
     */
    std::vector<TK*> pexsi_EDM;
#endif

  private:
    /**
     * @brief delete all HContainer objects in _DMR and clear the vector
     */
    void clear_DMR();

    /**
     * @brief HContainer for density matrix in real space for 2D parallelization
     * vector.size() = 1 for non-polarization and SOC
     * vector.size() = 2 for spin-polarization
     */
    std::vector<hamilt::HContainer<TR>*> _DMR;
    std::vector<std::vector<TR>> _DMR_save;

    /// @brief whether _DMR holds a density matrix calculated from DMK (reset by init_DMR, set by cal_DMR)
    bool _dmr_ready = false;

    /**
     * @brief HContainer for density matrix in real space for grid parallelization
     * same size semantics as _DMR
     */
    std::vector<hamilt::HContainer<TR>*> _DMR_grid;

    /**
     * @brief density matrix in k space, which is a vector[ik]
     * DMK should be a [spin_mult][_nk][i][j] matrix,
     * whose size is spin_mult * _nk * pv->get_nrow() * pv->get_ncol()
     */
    // std::vector<ModuleBase::ComplexMatrix> _DMK;
    std::vector<std::vector<TK>> _DMK;

    /**
     * @brief K_Vectors object, which is used to get k-point information
     */
    const std::vector<ModuleBase::Vector3<double>> _kvec_d;

    /**
     * @brief Parallel_Orbitals object, which contain all information of 2D block cyclic distribution
     */
    const Parallel_Orbitals* pv = nullptr;

    /**
     * @brief spin multiplicity used to size the density matrix (1 - none spin and SOC ;
     * 2 - spin polarization). This is NOT the physical nspin (1/2/4); it matches
     * K_Vectors::spin_mult.
     */
    int spin_mult = 1;

    /**
     * @brief the global physical nspin from INPUT (1/2/4).
     * For SOC/noncollinear (nspin==4) the density matrix is stored with spin_mult==1
     * (a single 2x2 spin-block matrix), but cal_DMR/cal_DMR_td must still take the
     * spin-resolved (Pauli) branch. spin_mult cannot distinguish this, so keep the
     * global value here. Equals spin_mult for non-SOC cases.
     */
    int nspin = 1;

    /**
     * @brief real number of k-points
     * _nk is not equal to _kv->get_nks() when spin-polarization is considered
     * _nk = kv->get_nks() / nspin when nspin=2
     */
    int _nk = 0;

    /// temporary pointers for switch DMR, only used with nspin=2
    std::vector<TR> dmr_origin_;
    std::vector<TR> dmr_tmp_;

    friend void DensityMatrix_Tools::cal_DMR<TK, TR>(
        DensityMatrix<TK, TR>& dm,
        std::vector<hamilt::HContainer<TR>*>& dmR_out,
        const int ik_in);
    friend void DensityMatrix_Tools::cal_DMR_td<TK, TR>(
        DensityMatrix<TK, TR>& dm,
        std::vector<hamilt::HContainer<TR>*>& dmR_out,
        const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
        const ModuleBase::Vector3<double> At,
        const int ik_in);
    friend void DensityMatrix_Tools::cal_DMR_full<TK, TR>(
        const DensityMatrix<TK, TR>& dm,
        hamilt::HContainer<std::complex<double>>* dmR_out,
        const int ik_in);
};

} // namespace module_dm

#endif
