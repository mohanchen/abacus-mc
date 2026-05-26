// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================
// Added by Shengjun Chen on 2026-05-26.
// NOT COMPLETE YET. DO NOT USE THIS CODE.
// ============================================================

#include "qlist.h"

#include "k_vector_utils.h"
#include "source_base/formatter.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_io/module_parameter/parameter.h"

void Q_Vectors::cal_iq_global()
{
    const int my_pool = this->para_q.my_pool;
    this->iq2iqtot.resize(this->nqs);
#ifdef __MPI
    for (int iq = 0; iq < this->nqs; ++iq)
    {
        this->iq2iqtot[iq] = this->para_q.startk_pool[my_pool] + iq;
    }
#else
    for (int iq = 0; iq < this->nqs; ++iq)
    {
        this->iq2iqtot[iq] = iq;
    }
#endif
}

void Q_Vectors::set(const UnitCell& ucell,
                    const ModuleSymmetry::Symmetry& symm,
                    const std::string& q_file_name,
                    const ModuleBase::Matrix3& reciprocal_vec,
                    const ModuleBase::Matrix3& latvec,
                    std::ofstream& ofs)
{
    ModuleBase::TITLE("Q_Vectors", "set");

    ofs << "\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " |                       #Setup Q-points#                             |" << std::endl;
    ofs << " | We setup the q-points according to input parameters.               |" << std::endl;
    ofs << " | The reduced q-points are set according to symmetry operations.     |" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    ofs << "\n";

    ofs << "\n SETUP Q-POINTS" << std::endl;

    // (1) read qpoints.
    bool read_succesfully = this->read_qpoints(ucell, q_file_name);
#ifdef __MPI
    Parallel_Common::bcast_bool(read_succesfully);
#endif
    if (!read_succesfully)
    {
        ModuleBase::WARNING_QUIT("Q_Vectors::set", "Something wrong while reading QPOINTS.");
    }

    // output qpoints file
    std::string sqpt1;
    std::string sqpt2;

    if (!this->qc_done && this->qd_done)
    {
        for (size_t iq = 0; iq != this->nqstot_full; ++iq)
            this->qvec_c_full[iq] = this->qvec_d[iq] * reciprocal_vec;
    }
    else if (this->qc_done && !this->qd_done)
    {
        for (size_t iq = 0; iq != this->nqstot_full; ++iq)
            this->qvec_c_full[iq] = this->qvec_c[iq];
    }

    // (2)
    // q-point symmetry reduction differs from k-points:
    // q-points do not involve berry phase or time-reversal-specific handling.
    // Acoustic sum rules should be enforced separately in phonon calculations.
    if (ModuleSymmetry::Symmetry::symm_flag != -1)
    {
        bool match = true;
        // calculate qpoints in IBZ and reduce qpoints according to symmetry
        KVectorUtils::kvec_ibz_kpoint(*this, symm, ModuleSymmetry::Symmetry::symm_flag, sqpt1, ucell, match);
#ifdef __MPI
        Parallel_Common::bcast_bool(match);
#endif
        if (!match)
        {
            std::cout << "Optimized lattice type of reciprocal lattice cannot match the optimized real lattice. "
                      << std::endl;
            std::cout << "It is often because the inaccuracy of lattice parameters in STRU." << std::endl;
            if (ModuleSymmetry::Symmetry::symm_autoclose)
            {
                ModuleBase::WARNING("Q_Vectors::ibz_qpoint", "Automatically set symmetry to 0 and continue ...");
                std::cout << "Automatically set symmetry to 0 and continue ..." << std::endl;
                ModuleSymmetry::Symmetry::symm_flag = 0;
                match = true;
                KVectorUtils::kvec_ibz_kpoint(*this, symm, ModuleSymmetry::Symmetry::symm_flag, sqpt1, ucell, match);
            } else {
                ModuleBase::WARNING_QUIT("Q_Vectors::ibz_qpoint",
                                         "Possible solutions: \n \
1. Refine the lattice parameters in STRU;\n \
2. Use a different`symmetry_prec`.  \n \
3. Close symemtry: set `symmetry` to 0 in INPUT. \n \
4. Set `symmetry_autoclose` to 1 in INPUT to automatically close symmetry when this error occurs.");
            }
        }
    }

    // (3)
    // Improve q point information

    // Complement the coordinates of q point
    KVectorUtils::set_both_kvec(*this, reciprocal_vec, latvec, sqpt2);

    if (GlobalV::MY_RANK == 0)
    {
        // output qpoints file
        std::stringstream sqpt;
        sqpt << PARAM.globalv.global_out_dir << "QPT.info";
        std::ofstream ofqpt(sqpt.str().c_str()); // clear qpoints
        ofqpt << sqpt2 << sqpt1;
        ofqpt.close();
    }

    // normalize q points weights (phonons have no spin degeneracy)
    this->normalize_wq();

    // It's very important in parallel case,
    // firstly do the mpi_q() and then
    // do other operations
    this->para_q.kinfo(nqstot,
                       GlobalV::KPAR,
                       GlobalV::MY_POOL,
                       GlobalV::RANK_IN_POOL,
                       GlobalV::NPROC,
                       1); // assign q points to several process pools (no spin for phonons)
#ifdef __MPI
    // distribute Q point data to the corresponding process
    KVectorUtils::kvec_mpi_k(*this);
#endif

    // initialize ibz_index
    this->ibz_index.resize(this->nqstot_full);
    for (int iq = 0; iq < this->nqstot_full; iq++)
    {
        this->ibz_index[iq] = iq;
    }
    
    // get iq2iqtot
    this->cal_iq_global();

    KVectorUtils::print_klists(*this, ofs);

    return;
}

// 1.reset the size of the Q-point container according to nqstot
void Q_Vectors::renew(const int& qpoint_number)
{
    qvec_c.resize(qpoint_number);
    qvec_d.resize(qpoint_number);
    qvec_c_full.resize(qpoint_number);
    wq.resize(qpoint_number);
    ngq.resize(qpoint_number);

    return;
}

// Read the QPT file, which contains Q-point coordinates, weights, and grid size information
// Generate Q-point grid according to different parameters of the QPT file
bool Q_Vectors::read_qpoints(const UnitCell& ucell,
                             const std::string& fn)
{
    ModuleBase::TITLE("Q_Vectors", "read_qpoints");
    if (GlobalV::MY_RANK != 0)
    {
        return true;
    }

    // 1. Overwrite the QPT file and default Q-point information if needed
    if (PARAM.globalv.gamma_only_local)
    {
        GlobalV::ofs_warning << " Auto generating q-points file: " << fn << std::endl;
        std::ofstream ofs(fn.c_str());
        ofs << "Q_POINTS" << std::endl;
        ofs << "0" << std::endl;
        ofs << "Gamma" << std::endl;
        ofs << "1 1 1 0 0 0" << std::endl;
        ofs.close();
    }
    else if (PARAM.inp.kspacing[0] > 0.0)
    {
        if (PARAM.inp.kspacing[1] <= 0 || PARAM.inp.kspacing[2] <= 0)
        {
            ModuleBase::WARNING_QUIT("Q_Vectors", "kspacing should > 0");
        };
        // number of Q points = max(1,int(|bi|/KSPACING+1))
        ModuleBase::Matrix3 btmp = ucell.G;
        double b1 = sqrt(btmp.e11 * btmp.e11 + btmp.e12 * btmp.e12 + btmp.e13 * btmp.e13);
        double b2 = sqrt(btmp.e21 * btmp.e21 + btmp.e22 * btmp.e22 + btmp.e23 * btmp.e23);
        double b3 = sqrt(btmp.e31 * btmp.e31 + btmp.e32 * btmp.e32 + btmp.e33 * btmp.e33);
        int nq1
            = std::max(1, static_cast<int>(b1 * ModuleBase::TWO_PI / PARAM.inp.kspacing[0] / ucell.lat0 + 1));
        int nq2
            = std::max(1, static_cast<int>(b2 * ModuleBase::TWO_PI / PARAM.inp.kspacing[1] / ucell.lat0 + 1));
        int nq3
            = std::max(1, static_cast<int>(b3 * ModuleBase::TWO_PI / PARAM.inp.kspacing[2] / ucell.lat0 + 1));

        GlobalV::ofs_warning << " Generate q-points file according to KSPACING: " << fn << std::endl;
        std::ofstream ofs(fn.c_str());
        ofs << "Q_POINTS" << std::endl;
        ofs << "0" << std::endl;
        if (PARAM.inp.kmesh_type == "mp")
        {
            ofs << "Monkhorst-Pack" << std::endl;
        }
        else
        {
            ofs << "Gamma" << std::endl;
        }
        ofs << nq1 << " " << nq2 << " " << nq3 << " " << PARAM.inp.koffset[0] << " " << PARAM.inp.koffset[1] << " "
            << PARAM.inp.koffset[2] << std::endl;
        ofs.close();
    }

    // 2. Generate the Q-point grid automatically according to the QPT file
    // 2.1 read the QPT file
    std::ifstream ifq(fn.c_str());
    if (!ifq)
    {
        GlobalV::ofs_warning << " Can't find File name : " << fn << std::endl;
        return false;
    }

    ifq >> std::setiosflags(std::ios::uppercase);

    ifq.clear();
    ifq.seekg(0);

    std::string word;
    std::string qword;

    int ierr = 0;

    ifq.rdstate();

    while (ifq.good())
    {
        ifq >> word;
        ifq.ignore(150, '\n'); // LiuXh add 20180416, fix bug in q-point file when the first line with comments
        if (word == "Q_POINTS" || word == "QPOINTS" || word == "Q")
        {
            ierr = 1;
            break;
        }

        ifq.rdstate();
    }

    if (ierr == 0)
    {
        GlobalV::ofs_warning << " symbol Q_POINTS not found." << std::endl;
        return false;
    }

    // input q-points are in 2pi/a units
    ModuleBase::GlobalFunc::READ_VALUE(ifq, nqstot);

    this->q_nqstot = nqstot;

    ModuleBase::GlobalFunc::READ_VALUE(ifq, qword);

    this->q_qword = qword;

    const int max_qpoints = 100000;
    if (nqstot > max_qpoints)
    {
        GlobalV::ofs_warning << " nqstot > MAX_QPOINTS" << std::endl;
        return false;
    }

    // 2.2 Select different methods and generate Q-point grid
    int q_type = 0;
    if (nqstot == 0) // nqstot==0, use monkhorst_pack.
    {
        if (qword == "Gamma") // MP(Gamma)
        {
            is_mp = true;
            q_type = 0;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Input type of q points", "Monkhorst-Pack(Gamma)");
        }
        else if (qword == "Monkhorst-Pack" || qword == "MP" || qword == "mp")
        {
            is_mp = true;
            q_type = 1;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Input type of q points", "Monkhorst-Pack");
        }
        else
        {
            GlobalV::ofs_warning << " Error: neither Gamma nor Monkhorst-Pack." << std::endl;
            return false;
        }

        ifq >> nmp[0] >> nmp[1] >> nmp[2];

        qoffset[0] = 0;
        qoffset[1] = 0;
        qoffset[2] = 0;
        if (!(ifq >> qoffset[0] >> qoffset[1] >> qoffset[2]))
        {
            ModuleBase::WARNING("Q_Vectors::read_qpoints", "Missing q-point offsets in the q-points file.");
        }

        this->Monkhorst_Pack(nmp, qoffset, q_type);
    }
    else if (nqstot > 0) // nqstot>0, the Q-point information is clearly set
    {
        if (qword == "Cartesian" || qword == "C") // Cartesian coordinates
        {
            this->renew(nqstot);
            for (int i = 0; i < nqstot; i++)
            {
                ifq >> qvec_c[i].x >> qvec_c[i].y >> qvec_c[i].z;
                ModuleBase::GlobalFunc::READ_VALUE(ifq, wq[i]);
            }

            this->qc_done = true;
        }
        else if (qword == "Direct" || qword == "D") // Direct coordinates
        {
            this->renew(nqstot);
            for (int i = 0; i < nqstot; i++)
            {
                ifq >> qvec_d[i].x >> qvec_d[i].y >> qvec_d[i].z;
                ModuleBase::GlobalFunc::READ_VALUE(ifq, wq[i]);
            }
            this->qd_done = true;
        }
        else if (qword == "Line_Cartesian")
        {
            if (ModuleSymmetry::Symmetry::symm_flag == 1)
            {
                ModuleBase::WARNING("Q_Vectors::read_qpoints",
                                    "Line mode of q-points is open, please set symmetry to 0 or -1.");
                return false;
            }

            interpolate_q_between(ifq, qvec_c);

            std::for_each(wq.begin(), wq.end(), [](double& d) { d = 1.0; });

            this->qc_done = true;
        }

        else if (qword == "Line_Direct" || qword == "L" || qword == "Line")
        {
            if (ModuleSymmetry::Symmetry::symm_flag == 1)
            {
                ModuleBase::WARNING("Q_Vectors::read_qpoints",
                                    "Line mode of q-points is open, please set symmetry to 0 or -1.");
                return false;
            }

            interpolate_q_between(ifq, qvec_d);

            std::for_each(wq.begin(), wq.end(), [](double& d) { d = 1.0; });

            this->qd_done = true;
        }

        else
        {
            GlobalV::ofs_warning << " Error : neither Cartesian nor Direct qpoint." << std::endl;
            return false;
        }
    }

    this->nqstot_full = this->nqs = this->nqstot;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nqstot", nqstot);
    return true;
} // END SUBROUTINE

void Q_Vectors::interpolate_q_between(std::ifstream& ifq, std::vector<ModuleBase::Vector3<double>>& qvec)
{
    // how many special points.
    int nqs_special = this->nqstot;

    // number of points to the next q points
    std::vector<int> nql(nqs_special, 0);

    // coordinates of special points.
    std::vector<ModuleBase::Vector3<double>> qs(nqs_special);

    // recalculate nqstot.
    nqstot = 0;
    /* ISSUE#3482: to distinguish different qline segments */
    std::vector<int> qpt_segids;
    ql_segids.clear();
    ql_segids.shrink_to_fit();
    int qpt_segid = 0;
    for (int iqs = 0; iqs < nqs_special; iqs++)
    {
        ifq >> qs[iqs].x;
        ifq >> qs[iqs].y;
        ifq >> qs[iqs].z;
        ModuleBase::GlobalFunc::READ_VALUE(ifq, nql[iqs]);

        assert(nql[iqs] >= 0);
        nqstot += nql[iqs];
        /* ISSUE#3482: to distinguish different qline segments */
        if ((nql[iqs] == 1) && (iqs != (nqs_special - 1))) {
            qpt_segid++;
        }
        qpt_segids.push_back(qpt_segid);
    }
    assert(nql[nqs_special - 1] == 1);

    this->renew(nqstot);

    int count = 0;
    for (int iqs = 1; iqs < nqs_special; iqs++)
    {
        double dxs = (qs[iqs].x - qs[iqs - 1].x) / nql[iqs - 1];
        double dys = (qs[iqs].y - qs[iqs - 1].y) / nql[iqs - 1];
        double dzs = (qs[iqs].z - qs[iqs - 1].z) / nql[iqs - 1];
        for (int is = 0; is < nql[iqs - 1]; is++)
        {
            qvec[count].x = qs[iqs - 1].x + is * dxs;
            qvec[count].y = qs[iqs - 1].y + is * dys;
            qvec[count].z = qs[iqs - 1].z + is * dzs;
            ql_segids.push_back(qpt_segids[iqs - 1]); /* ISSUE#3482: to distinguish different qline segments */
            ++count;
        }
    }

    // deal with the last special q point.
    qvec[count].x = qs[nqs_special - 1].x;
    qvec[count].y = qs[nqs_special - 1].y;
    qvec[count].z = qs[nqs_special - 1].z;
    ql_segids.push_back(qpt_segids[nqs_special - 1]); /* ISSUE#3482: to distinguish different qline segments */
    ++count;

    assert(count == nqstot);
    assert(ql_segids.size() == nqstot); /* ISSUE#3482: to distinguish different qline segments */
}

double Q_Vectors::Monkhorst_Pack_formula(const int& q_type, const double& offset, const int& n, const int& dim)
{
    double coordinate = 0.0;
    if (q_type == 1)
    {
        coordinate = (offset + 2.0 * (double)n - (double)dim - 1.0) / (2.0 * (double)dim);
    }
    else
    {
        coordinate = (offset + (double)n - 1.0) / (double)dim;
    }
    return coordinate;
}

void Q_Vectors::Monkhorst_Pack(const int* nmp_in, const double* qoffset_in, const int q_type)
{
    const int mpnx = nmp_in[0];
    const int mpny = nmp_in[1];
    const int mpnz = nmp_in[2];

    this->nqstot = mpnx * mpny * mpnz;
    // only can renew after nqstot is estimated.
    this->renew(nqstot);

    for (int x = 1; x <= mpnx; x++)
    {
        double v1 = Monkhorst_Pack_formula(q_type, qoffset_in[0], x, mpnx);
        if (std::abs(v1) < 1.0e-10) {
            v1 = 0.0;
        }
        for (int y = 1; y <= mpny; y++)
        {
            double v2 = Monkhorst_Pack_formula(q_type, qoffset_in[1], y, mpny);
            if (std::abs(v2) < 1.0e-10) {
                v2 = 0.0;
            }
            for (int z = 1; z <= mpnz; z++)
            {
                double v3 = Monkhorst_Pack_formula(q_type, qoffset_in[2], z, mpnz);
                if (std::abs(v3) < 1.0e-10) {
                    v3 = 0.0;
                }
                // index of nqs qpoint
                const int i = mpnx * mpny * (z - 1) + mpnx * (y - 1) + (x - 1);
                qvec_d[i].set(v1, v2, v3);
            }
        }
    }

    const double weight = 1.0 / static_cast<double>(nqstot);
    for (int iq = 0; iq < nqstot; iq++)
    {
        wq[iq] = weight;
    }
    this->qd_done = true;

    return;
}

void Q_Vectors::update_use_ibz(const int& nqstot_ibz,
                               const std::vector<ModuleBase::Vector3<double>>& qvec_d_ibz,
                               const std::vector<double>& wq_ibz)
{
    if (GlobalV::MY_RANK != 0) {
        return;
    }
    ModuleBase::TITLE("Q_Vectors", "update_use_ibz");
    assert(nqstot_ibz > 0);
    assert(nqstot_ibz <= qvec_d_ibz.size());
    // update nqstot
    this->nqs = this->nqstot = nqstot_ibz;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nqstot now", nqstot);

    this->qvec_d.resize(this->nqstot);

    for (int i = 0; i < this->nqstot; ++i)
    {
        this->qvec_d[i] = qvec_d_ibz[i];

        // update weight.
        this->wq[i] = wq_ibz[i];
    }

    this->qd_done = true;
    this->qc_done = false;
    return;
}

void Q_Vectors::normalize_wq()
{
    if (GlobalV::MY_RANK != 0) {
        return;
    }
    double sum = 0.0;

    for (int iq = 0; iq < nqstot; iq++)
    {
        sum += this->wq[iq];
    }

    // If sum of weights is zero or very small, set equal weights
    if (sum < 1e-10)
    {
        ModuleBase::WARNING("Q_Vectors::normalize_wq",
                            "Sum of q-point weights is zero or very small. "
                            "Setting equal weights for all q-points.");
        for (int iq = 0; iq < nqstot; iq++)
        {
            this->wq[iq] = 1.0 / double(nqstot);
        }
        sum = 1.0;
    }

    for (int iq = 0; iq < nqstot; iq++)
    {
        this->wq[iq] /= sum;
    }

    // phonons have no spin degeneracy; weights already sum to 1.0

    return;
}