#include "setup_nonlocal.h"

#include "source_base/parallel_common.h"

#ifdef __LCAO
#include "source_pw/module_pwdft/soc.h"
#include "../source_base/complexmatrix.h"
// mohan add 2013-08-02
// In order to get rid of the read in file .NONLOCAL.

InfoNonlocal::InfoNonlocal()
{
    this->Beta.resize(1);
    this->nprojmax = 0;
    this->rcutmax_Beta = 0.0;
}
InfoNonlocal::~InfoNonlocal() = default;

void InfoNonlocal::build_soc_coefficients(const Atom* atom,
                                          const int& n_projectors,
                                          ModuleBase::ComplexMatrix& coefficient_D_nc_in)
{
    const int nh = atom->ncpp.nh;

    int lmaxkb = -1;
    for (int ibeta = 0; ibeta < atom->ncpp.nbeta; ibeta++)
    {
        lmaxkb = std::max(lmaxkb, atom->ncpp.lll[ibeta]);
    }

    Soc soc;
    if (atom->ncpp.has_so)
    {
        soc.rot_ylm(lmaxkb);
        soc.fcoef.create(1, atom->ncpp.nh, atom->ncpp.nh);
    }

    int ip1 = 0;
    for (int p1 = 0; p1 < n_projectors; p1++)
    {
        const int l1 = atom->ncpp.lll[p1];
        const double j1 = atom->ncpp.jjj[p1];
        for (int m1 = 0; m1 < 2 * l1 + 1; m1++)
        {
            int ip2 = 0;
            for (int p2 = 0; p2 < n_projectors; p2++)
            {
                const int l2 = atom->ncpp.lll[p2];
                const double j2 = atom->ncpp.jjj[p2];
                for (int m2 = 0; m2 < 2 * l2 + 1; m2++)
                {
                    if (l1 == l2 && fabs(j1 - j2) < 1e-7)
                    {
                        for (int is1 = 0; is1 < 2; is1++)
                        {
                            for (int is2 = 0; is2 < 2; is2++)
                            {
                                if (atom->ncpp.has_so)
                                {
                                    soc.set_fcoef(l1, l2, is1, is2, m1, m2, j1, j2, 0, ip1, ip2);

                                    coefficient_D_nc_in(ip1 + nh * is1, ip2 + nh * is2)
                                        = atom->ncpp.dion(p1, p2) * soc.fcoef(0, is1, is2, ip1, ip2);
                                    if (p1 != p2)
                                    {
                                        soc.fcoef(0, is1, is2, ip1, ip2) = std::complex<double>(0.0, 0.0);
                                    }
                                }
                                else
                                {
                                    if (is1 == is2 && m1 == m2)
                                    {
                                        coefficient_D_nc_in(ip1 + nh * is1, ip2 + nh * is2) = atom->ncpp.dion(p1, p2);
                                    }
                                }
                            } // end is2
                        }     // end is1
                    }         // end l1==l2
                    ip2++;
                } // end m2
            }     // end p2
            assert(ip2 == nh);
            ip1++;
        } // end m1
    }
}

void InfoNonlocal::build_beta_r(const Atom* atom,
                                const int& p1,
                                std::vector<double>& beta_r,
                                int& cut_mesh)
{
    cut_mesh = atom->ncpp.mesh;
    for (int ir = atom->ncpp.mesh - 1; ir >= 0; --ir)
    {
        if (std::abs(atom->ncpp.betar(p1, ir)) > 1.0e-10)
        {
            cut_mesh = ir;
            break;
        }
    }
    if (cut_mesh % 2 == 0)
    {
        ++cut_mesh;
    }

    beta_r.resize(cut_mesh, 0.0);
    for (int ir = 0; ir < cut_mesh; ++ir)
    {
        beta_r[ir] = atom->ncpp.betar(p1, ir);
    }
}

void InfoNonlocal::read_header(std::ifstream& ifs,
                               const int& my_rank,
                               std::string& label,
                               std::string& ps_type,
                               int& nlmax)
{
    if (my_rank == 0)
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<HEADER>"))
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, label);
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ps_type);
            if (ps_type != "NC")
            {
                ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal",
                                         "Only available for NC nonlocal pseudopotential");
            }
            ModuleBase::GlobalFunc::READ_VALUE(ifs, nlmax);
            assert(nlmax >= -1);
            ModuleBase::GlobalFunc::SCAN_END(ifs, "</HEADER>");
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_string(label);
    Parallel_Common::bcast_string(ps_type);
    Parallel_Common::bcast_int(nlmax);
#endif
}

void InfoNonlocal::read_dij(std::ifstream& ifs,
                            const int& my_rank,
                            const int& nlmax,
                            int& n_projectors,
                            std::ofstream& log)
{
    if (my_rank == 0)
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<DIJ>"))
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, n_projectors);
            ModuleBase::GlobalFunc::OUT(log, "n_projectors", n_projectors);

            for (int p1 = 0; p1 < n_projectors; p1++)
            {
                for (int p2 = 0; p2 < n_projectors; p2++)
                {
                    int L1_read, L2_read;
                    ifs >> L1_read >> L2_read;
                    assert(L1_read <= nlmax);
                    assert(L2_read <= nlmax);
                    double dion_read;
                    ifs >> dion_read;
                }
            }
            ModuleBase::GlobalFunc::SCAN_END(ifs, "</DIJ>");
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(n_projectors);
#endif
}

void InfoNonlocal::read_projector(std::ifstream& ifs,
                                  const int& my_rank,
                                  const int& p1,
                                  const int& nlmax,
                                  int& meshr_ps,
                                  int& lfrombeta,
                                  std::vector<double>& radial_ps,
                                  std::vector<double>& rab_ps,
                                  std::vector<double>& beta_r)
{
    meshr_ps = 0;
    if (my_rank == 0)
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_BETA>", false))
        {
            int iproj = 0;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, iproj);
            if (iproj != p1)
            {
                std::cout << " iproj=" << iproj << " p1=" << p1 << std::endl;
                ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal", "Check non-local projector index.");
            }

            ModuleBase::GlobalFunc::READ_VALUE(ifs, lfrombeta);
            assert(lfrombeta >= 0);
            assert(lfrombeta <= nlmax);

            ModuleBase::GlobalFunc::READ_VALUE(ifs, meshr_ps);
            if (meshr_ps % 2 == 0)
            {
                std::cout << " meshr_ps = " << meshr_ps << std::endl;
                ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal", "meshr_ps must be odd!");
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal", "<PP_BETA> doesn't match!");
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(meshr_ps);
    Parallel_Common::bcast_int(lfrombeta);
#endif

    radial_ps.assign(meshr_ps, 0.0);
    rab_ps.assign(meshr_ps, 0.0);
    beta_r.assign(meshr_ps, 0.0);

    if (my_rank == 0)
    {
        for (int ir = 0; ir < meshr_ps; ir++)
        {
            ifs >> radial_ps[ir];
            ifs >> beta_r[ir];
            ifs >> rab_ps[ir];
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_double(radial_ps.data(), meshr_ps);
    Parallel_Common::bcast_double(beta_r.data(), meshr_ps);
    Parallel_Common::bcast_double(rab_ps.data(), meshr_ps);
#endif

    if (my_rank == 0)
    {
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_BETA>");
    }
}

void InfoNonlocal::Set_NonLocal(const int& it,
                                Atom* atom,
                                int& n_projectors,
                                const int& kmesh,
                                const double& dk,
                                const double& dr_uniform,
                                std::ofstream& log,
                                const bool& out_element_info,
                                const bool& lspinorb,
                                const int& nspin,
                                const int& my_rank)
{
    ModuleBase::TITLE("InfoNonlocal", "Set_NonLocal");

    // get the number of non-local projectors
    n_projectors = atom->ncpp.nbeta;
    const int nh = atom->ncpp.nh; // zhengdy-soc

    // set the nonlocal projector objects
    std::vector<Numerical_Nonlocal_Lm> tmpBeta_lm(n_projectors);
    ModuleBase::ComplexMatrix coefficient_D_nc_in(nh * 2, nh * 2); // zhengdy-soc

    build_soc_coefficients(atom, n_projectors, coefficient_D_nc_in);

    for (int p1 = 0; p1 < n_projectors; p1++)
    {
        const int lnow = atom->ncpp.lll[p1];

        int cut_mesh = 0;
        std::vector<double> beta_r;
        build_beta_r(atom, p1, beta_r, cut_mesh);

        tmpBeta_lm[p1].set_NL_proj(atom->label,
                                   it,       // type
                                   lnow,     // angular momentum L
                                   cut_mesh, // number of radial mesh
                                   atom->ncpp.rab.data(),
                                   atom->ncpp.r.data(), // radial mesh value (a.u.)
                                   beta_r.data(),
                                   kmesh,
                                   dk,
                                   dr_uniform); // delta k mesh in reciprocal space

        if (out_element_info)
        {
            tmpBeta_lm[p1].plot(my_rank);
        }
    }

    this->Beta[it].set_type_info(it,
                                 atom->label,
                                 atom->ncpp.pp_type,
                                 atom->ncpp.lmax,
                                 n_projectors,
                                 tmpBeta_lm.data()); // zhengdy-soc 2018-09-10

    // mohan add 2021-05-07
    atom->ncpp.set_d_so(coefficient_D_nc_in, n_projectors, nh, atom->ncpp.has_so, lspinorb, nspin);

    log << " SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS FOR ELEMENT " << atom->label << std::endl;
}

void InfoNonlocal::Read_NonLocal(const int& it,
                                 Atom* atom,
                                 int& n_projectors,
                                 const int& my_rank,
                                 const int& kmesh,
                                 const double& dk,
                                 const double& dr_uniform,
                                 const std::string& nonlocalFile,
                                 const bool& out_element_info,
                                 std::ofstream& log)
{
    ModuleBase::TITLE("InfoNonlocal", "Read_NonLocal");

    std::ifstream ifs;

    // check if the non-local pseudopotential file exist
    bool open = false;
    if (my_rank == 0)
    {
        ifs.open(nonlocalFile.c_str());
        if (ifs)
        {
            open = true;
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_bool(open);
#endif
    if (!open)
    {
        std::cout << " Non-local File : " << nonlocalFile << std::endl;
        ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal", "Can not find the NONLOCAL file.");
    }

    std::string label;
    std::string ps_type;
    int nlmax = 0;

    read_header(ifs, my_rank, label, ps_type, nlmax);

    // mohan add 2012-06-09
    if (nlmax != -1)
    {
        bool find_lmax = false;
        for (int ic = 0; ic < atom->ncpp.nbeta; ic++)
        {
            if (nlmax == atom->ncpp.lll[ic])
            {
                find_lmax = true;
                break;
            }
        }

        if (!find_lmax)
        {
            std::cout << " For element " << label << std::endl;
            std::cout << " Max L Read in from NONLOCAL = " << nlmax << std::endl;
            for (int ib = 0; ib < atom->ncpp.nbeta; ++ib)
            {
                std::cout << " Max L Read in from pseudopotential file = " << atom->ncpp.lll[ib] << std::endl;
            }
            ModuleBase::WARNING_QUIT("InfoNonlocal::Read_NonLocal", "nlmax != atom->lll");
        }
    }

    ModuleBase::GlobalFunc::OUT(log, "label", label);
    ModuleBase::GlobalFunc::OUT(log, "nlmax", nlmax);

    read_dij(ifs, my_rank, nlmax, n_projectors, log);

    std::vector<Numerical_Nonlocal_Lm> tmpBeta_lm(n_projectors);
    std::vector<int> LfromBeta(n_projectors, 0);

    for (int p1 = 0; p1 < n_projectors; p1++)
    {
        int meshr_ps = 0;
        int lfrombeta = 0;
        std::vector<double> radial_ps;
        std::vector<double> rab_ps;
        std::vector<double> beta_r;

        read_projector(ifs, my_rank, p1, nlmax, meshr_ps, lfrombeta, radial_ps, rab_ps, beta_r);
        LfromBeta[p1] = lfrombeta;

        tmpBeta_lm[p1].set_NL_proj(label,
                                   it,            // type
                                   LfromBeta[p1], // angular momentum L
                                   meshr_ps,      // number of radial mesh
                                   rab_ps.data(),
                                   radial_ps.data(), // radial mesh value(a.u.)
                                   beta_r.data(),
                                   kmesh,
                                   dk,
                                   dr_uniform); // delta k mesh in reciprocal space

        if (out_element_info)
        {
            tmpBeta_lm[p1].plot(my_rank);
        }
    }

    this->Beta[it].set_type_info(it, label, ps_type, nlmax, n_projectors, tmpBeta_lm.data());

    ifs.close();
}

void InfoNonlocal::setupNonlocal(const int& ntype, Atom* atoms, std::ofstream& log, LCAO_Orbitals& orb,
                                 const std::string& basis_type,
                                 const bool& out_element_info,
                                 const bool& lspinorb,
                                 const int& nspin,
                                 const int& my_rank)
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //~~~~~~~~~~~~~~~~~~~~~~   2    ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Read in non-local projector for each atom type.
    // In fact this should be improved,
    // the non-local projector should be transferred
    // from .UPF file directly.
    // mohan note 2011-03-04
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (basis_type == "lcao" || basis_type == "lcao_in_pw")
    {
        this->Beta.resize(ntype);
        this->nproj.assign(ntype, 0);

        this->nprojmax = 0;

        // if true: read in the nonlocal file from file.
        // if false: get nonlocal information from .upf or .vwr directly
        bool readin_nonlocal = false;

        for (int it = 0; it < ntype; it++)
        {
            Atom* atom = &atoms[it];
            if (readin_nonlocal)
            {
                this->Read_NonLocal(it,
                                    atom,
                                    this->nproj[it],
                                    my_rank,
                                    orb.get_kmesh(),
                                    orb.get_dk(),
                                    orb.get_dr_uniform(),
                                    orb.orbital_file[it],
                                    out_element_info,
                                    log);
            }
            else
            {
                this->Set_NonLocal(it, atom, this->nproj[it], orb.get_kmesh(), orb.get_dk(), orb.get_dr_uniform(), log,
                                   out_element_info, lspinorb, nspin, my_rank);
            }
            this->nprojmax = std::max(this->nprojmax, this->nproj[it]);
            // caoyu add 2021-05-24 to reconstruct atom_arrange::set_sr_NL
            this->rcutmax_Beta = std::max(this->rcutmax_Beta, this->Beta[it].get_rcut_max());
        }

        ModuleBase::GlobalFunc::OUT(log, "Max number of nonlocal projectors (all elements)", this->nprojmax);
        
    }
    return;
}

#endif
