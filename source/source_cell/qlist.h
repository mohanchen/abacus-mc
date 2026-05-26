// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================
// Added by Shengjun Chen on 2026-05-26.
// ============================================================


#ifndef Q_VECTORS_H
#define Q_VECTORS_H

#include "source_base/global_function.h"
#include "source_base/matrix3.h"
#include "source_cell/unitcell.h"
#include "parallel_kpoints.h"
#include "k_vector_utils.h"
#include <vector>

class Q_Vectors
{
public:
    std::vector<ModuleBase::Vector3<double>> qvec_c; /// Cartesian coordinates of q points
    std::vector<ModuleBase::Vector3<double>> qvec_d; /// Direct coordinates of q points
    std::vector<ModuleBase::Vector3<double>> qvec_c_full; // Cartesian coordinates of full q mesh match with nqstot_full

    std::vector<double> wq; /// wq, weight of q points

    std::vector<int> ngq; /// ngq, number of plane waves for each q point

    int nmp[3]={0};                 /// Number of Monkhorst-Pack
    std::vector<int> ql_segids; /// index of qline segment

    /// @brief equal q points to each ibz-qpoint, corresponding to a certain symmetry operations. 
    /// dim: [iqs_ibz][(isym, qvec_d)]
    std::vector<std::map<int, ModuleBase::Vector3<double>>> qstars;

    bool qc_done = false;
    bool qd_done = false;

    Q_Vectors(){};
    ~Q_Vectors(){};
    Q_Vectors& operator=(const Q_Vectors&) = default;
    Q_Vectors& operator=(Q_Vectors&& rhs) = default;

    Parallel_Kpoints para_q; ///< parallel for qpoints

    /**
     * @brief Set up the q-points for phonon calculation.
     *
     * This function sets up the q-points according to the input parameters and symmetry operations.
     *
     * @param symm The symmetry of the system.
     * @param q_file_name The name of the file containing the q-points.
     * @param reciprocal_vec The reciprocal vector of the system.
     * @param latvec The lattice vector of the system.
     *
     * @return void
     *
     * @note This function will quit with a warning if something goes wrong while reading the QPOINTS file.
     * @note If the optimized lattice type of the reciprocal lattice cannot match the optimized real lattice,
     *       it will output a warning and suggest possible solutions.
     */
    void set(const UnitCell& ucell,
        const ModuleSymmetry::Symmetry& symm,
        const std::string& q_file_name,
        const ModuleBase::Matrix3& reciprocal_vec,
        const ModuleBase::Matrix3& latvec,
        std::ofstream& ofs);

    int get_nqs() const
    {
        return this->nqs;
    }

    int get_nqstot() const
    {
        return this->nqstot;
    }

    int get_nqstot_full() const
    {
        return this->nqstot_full;
    }

    double get_qoffset(const int i) const
    {
        return this->qoffset[i];
    }

    int get_q_nqstot() const
    {
        return this->q_nqstot;
    }

    std::string get_q_qword() const
    {
        return this->q_qword;
    }

    void set_nqs(int value)
    {
        this->nqs = value;
    }

    void set_nqstot(int value)
    {
        this->nqstot = value;
    }

    void set_nqstot_full(int value)
    {
        this->nqstot_full = value;
    }

    bool get_is_mp() const
    {
        return is_mp;
    }

    std::vector<int> iq2iqtot; ///<[nqs] map iq to the global index of q points
    std::vector<int> ibz_index; ///< map q points (before symmetry reduction) to irreducible q-points

    /**
     * @brief Updates the q-points to use the irreducible Brillouin zone (IBZ).
     *
     * This function updates the q-points to use the irreducible Brillouin zone (IBZ) instead of the full Brillouin
     * zone.
     *
     * @return void
     *
     * @note This function should only be called by the master process (MY_RANK == 0).
     * @note This function assumes that the number of q-points in the IBZ (nqstot_ibz) is greater than 0.
     * @note This function updates the total number of q-points (nqstot) to be the number of q-points in the IBZ.
     * @note This function resizes the vector of q-points (qvec_d) and updates its values to be the q-points in the IBZ.
     * @note This function also updates the weights of the q-points (wq) to be the weights in the IBZ.
     * @note After this function is called, the flag qd_done is set to true to indicate that the q-points have been
     * updated, and the flag qc_done is set to false to indicate that the Cartesian coordinates of the q-points need to
     * be recalculated.
     */
    void update_use_ibz(const int& nqstot_ibz,
                        const std::vector<ModuleBase::Vector3<double>>& qvec_d_ibz,
                        const std::vector<double>& wq_ibz);

  private:
    int nqs = 0;         ///< number of symmetry-reduced q points in this pool(processor)
    int nqstot = 0;      ///< number of symmetry-reduced q points in full q mesh
    int nqstot_full = 0; ///< number of q points before symmetry reduction in full q mesh

    double qoffset[3] = {0.0}; // used only in automatic q-points.
    std::string q_qword;       // coordinate type keyword in QPOINTS file
    int q_nqstot = 0;          // original number of q points declared in input
    bool is_mp = false;        // Monkhorst-Pack

    /**
     * @brief Resize the q-point related vectors according to the new q-point number.
     *
     * This function resizes the vectors that store the q-point information,
     * including the Cartesian and Direct coordinates of q-points,
     * the weights of q-points, the index of q-points, and the number of plane waves for each q-point.
     *
     * @param qpoint_number The new number of q-points.
     *
     * @return void
     */
    void renew(const int& qpoint_number);

    // step 1 : generate qpoints

    /**
     * @brief Reads the q-points from a file.
     *
     * This function reads the q-points from a file specified by the filename.
     * It supports both Cartesian and Direct coordinates, and can handle different types of q-points,
     * including Gamma, Monkhorst-Pack, and Line mode. It also supports automatic generation of q-points
     * file if the file does not exist.
     *
     * @param fn The name of the file containing the q-points.
     *
     * @return bool Returns true if the q-points are successfully read from the file,
     *              false otherwise.
     *
     * @note It will generate a q-points file automatically according to the global variables GAMMA_ONLY_LOCAL and QSPACING (if applicable).
     * @note If the q-points type is neither Gamma nor Monkhorst-Pack, it will quit with a warning.
     * @note If the q-points type is Line mode and the symmetry flag is 1, it will quit with a warning.
     * @note If the number of q-points is greater than 100000, it will quit with a warning.
     */
    bool read_qpoints(const UnitCell& ucell,
                      const std::string& fn); // return 0: something wrong.

    /**
     * @brief Adds q-points linearly between special points.
     *
     * This function adds q-points linearly between special points in the Brillouin zone.
     * The special points and the number of q-points between them are read from an input file.
     *
     * @param ifq The input file stream from which the special points and the number of q-points between them are read.
     * @param qvec A vector to store the generated q-points.
     *
     * @return void
     *
     * @note The function first reads the number of special points (nqs_special) and the number of q-points between them
     * (nql) from the input file.
     * @note The function then recalculates the total number of q-points (nqstot) based on the number of q-points
     * between the special points.
     * @note The function generates the q-points by linearly interpolating between the special points.
     * @note The function also assigns a segment ID to each q-point to distinguish different q-line segments.
     * @note The function checks that the total number of generated q-points matches the calculated total number of
     * q-points.
     * @note The function checks that the size of the segment ID vector matches the total number of q-points.
     */
    void interpolate_q_between(std::ifstream& ifq, std::vector<ModuleBase::Vector3<double>>& qvec);

    /**
     * @brief Generates q-points using the Monkhorst-Pack scheme.
     *
     * This function generates q-points in the reciprocal space using the Monkhorst-Pack scheme.
     *
     * @param nmp_in the number of q-points in each dimension.
     * @param qoffset_in the offset for the q-points in each dimension.
     * @param q_type The type of q-point.  1 means without Gamma point, 0 means with Gamma.
     *
     * @return void
     *
     * @note The function assumes that the q-points are evenly distributed in the reciprocal space.
     * @note The function sets the weight of each q-point to be equal, so that the total weight of all q-points is 1.
     * @note The function sets the flag qd_done to true to indicate that the q-points have been generated.
     */
    void Monkhorst_Pack(const int* nmp_in, const double* qoffset_in, const int tipo);

    /**
     * @brief Calculates the coordinate of a q-point using the Monkhorst-Pack scheme.
     *
     * This function calculates the coordinate of a q-point in the reciprocal space using the Monkhorst-Pack scheme.
     * The Monkhorst-Pack scheme is a method for generating q-points in the Brillouin zone.
     *
     * @param q_type The type of q-point. 1 means without Gamma point, 0 means with Gamma.
     * @param offset The offset for the q-point.
     * @param n The index of the q-point in the current dimension.
     * @param dim The total number of q-points in the current dimension.
     *
     * @return double Returns the coordinate of the q-point.
     *
     * @note The function assumes that the q-points are evenly distributed in the reciprocal space.
     */
    double Monkhorst_Pack_formula(const int& q_type, const double& offset, const int& n, const int& dim);

    // step 2 : set both qvec and qved; normalize weight

    /**
     * @brief Normalizes the weights of the q-points.
     *
     * This function normalizes the weights of the q-points so that their sum is equal to 1.
     *
     * @return void
     *
     * @note This function should only be called by the master process (MY_RANK == 0).
     * @note If the sum of the weights is zero or very small (< 1e-10), the function will set equal weights for all
     * q-points and issue a warning.
     */
    void normalize_wq();

    /**
     * @brief Gets the global index of a q-point.
     * @return this->iq2iqtot[iq]
     */
    void cal_iq_global();
#ifdef __MPI
    friend void KVectorUtils::kvec_mpi_k(Q_Vectors& qvec);
#endif
};
#endif // QVECT_H