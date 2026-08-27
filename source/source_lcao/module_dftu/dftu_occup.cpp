#include "dftu_lcao.h"
#include "dftu_folding.h"
#include "source_base/timer.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __LCAO
#include "source_lcao/hamilt_lcao.h"
#endif

// copy_occ_mat(), zero_occ_mat(), mix_occ_mat(), set_occ_mat(ucell),
// get_occ_mat_flat(), set_occ_mat_flat()
// are now implemented in dftu_base.cpp as Plus_U_Base methods (inherited by Plus_U).

#ifdef __LCAO

void Plus_U::cal_occ_mat_k(const Parallel_Orbitals* pv,
                         const int iter,
                         const UnitCell& ucell,
                         const std::vector<std::vector<std::complex<double>>>& dm_k,
                         const K_Vectors& kv,
                         const double& mixing_beta,
                         hamilt::Hamilt<std::complex<double>>* p_ham,
                         const bool gamma_only_local,
                         const int nspin,
                         const int npol,
                         const int nlocal,
                         const std::string& ks_solver,
                         const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>& iatlnmipol2iwt,
                         const std::vector<int>& orbital_corr)
{
    ModuleBase::TITLE("Plus_U", "cal_occ_mat_k");
    ModuleBase::timer::start("Plus_U", "cal_occ_mat_k");

    this->copy_occ_mat(ucell);
    this->zero_occ_mat(ucell);

    //=================Part 1======================
    // call SCALAPACK routine to calculate the product of the S and density matrix
    const char transN = 'N';
    const char transT = 'T';
    const int one_int = 1;
    const std::complex<double> beta(0.0,0.0), alpha(1.0,0.0);

    std::vector<std::complex<double>> srho(pv->nloc);

    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        // srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_k(iw,nu)
        DFTU_LCAO::folding_matrix_k_new(ks_solver, gamma_only_local, nspin, ik, p_ham);

        std::complex<double>* s_k_pointer = nullptr;

        if(nspin != 4)
        {
            s_k_pointer = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham)->getSk();
        }
        else
        {
            s_k_pointer = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham)->getSk();
        }

#ifdef __MPI
        ScalapackConnector::gemm(transN,
            transT,
            nlocal,
            nlocal,
            nlocal,
            alpha,
            s_k_pointer,
            one_int,
            one_int,
            &pv->desc[0],
            dm_k[ik].data(),
            one_int,
            one_int,
            &pv->desc[0],
            beta,
            srho.data(),
            one_int,
            one_int,
            &pv->desc[0]);
#endif

        const int spin = kv.isk[ik];
        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

			if (LC == -1)
			{
				continue;
			}

            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
					if (l != orbital_corr[it])
					{
						continue;
					}

                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        // if(!Yukawa && n!=0) continue;
						if (n != 0)
						{
							continue;
						}

                        // Calculate the local occupation number matrix
                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < npol; ipol0++)
                            {
                                const int iwt0 = iatlnmipol2iwt[iat][l][n][m0][ipol0];
                                const int mu = pv->global2local_row(iwt0);
                                const int mu_prime = pv->global2local_col(iwt0);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                    {
                                        const int iwt1 = iatlnmipol2iwt[iat][l][n][m1][ipol1];
                                        const int nu = pv->global2local_col(iwt1);
                                        const int nu_prime = pv->global2local_row(iwt1);

                                        const int irc = nu * pv->nrow + mu;
                                        const int irc_prime = mu_prime * pv->nrow + nu_prime;

                                        const int m0_all = m0 + ipol0 * (2 * l + 1);
                                        const int m1_all = m1 + ipol1 * (2 * l + 1);

										if ((nu >= 0) && (mu >= 0))
										{
											occ_mat[iat][l][n][spin](m0_all, m1_all) += (srho[irc]).real() / 4.0;
										}

										if ((nu_prime >= 0) && (mu_prime >= 0))
										{
											occ_mat[iat][l][n][spin](m0_all, m1_all)
												+= (std::conj(srho[irc_prime])).real() / 4.0;
										}
                                    } // ipol1
                                } // m1
                            } // ipol0
                        } // m0
                    } // end n
                } // end l
            } // end ia
        } // end it
    } // ik

    for (int it = 0; it < ucell.ntype; it++)
    {
        const int NL = ucell.atoms[it].nwl + 1;
        const int LC = orbital_corr[it];

		if (LC == -1)
		{
			continue;
		}

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const int iat = ucell.itia2iat(it, ia);

            for (int l = 0; l < NL; l++)
            {
				if (l != orbital_corr[it])
				{
					continue;
				}

                const int N = ucell.atoms[it].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    // if(!Yukawa && n!=0) continue;
					if (n != 0)
					{
						continue;
					}
					// set the local occupation mumber matrix of spin up and down zeros

#ifdef __MPI
                    if (nspin == 1 || nspin == 4)
                    {
                        ModuleBase::matrix temp(occ_mat[iat][l][n][0]);
                        MPI_Allreduce(&temp(0, 0),
                                      &occ_mat[iat][l][n][0](0, 0),
                                      (2 * l + 1) * npol * (2 * l + 1) * npol,
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
                    }
                    else if (nspin == 2)
                    {
                        ModuleBase::matrix temp0(occ_mat[iat][l][n][0]);
                        MPI_Allreduce(&temp0(0, 0),
                                      &occ_mat[iat][l][n][0](0, 0),
                                      (2 * l + 1) * (2 * l + 1),
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);

                        ModuleBase::matrix temp1(occ_mat[iat][l][n][1]);
                        MPI_Allreduce(&temp1(0, 0),
                                      &occ_mat[iat][l][n][1](0, 0),
                                      (2 * l + 1) * (2 * l + 1),
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
                    }
#endif

                    switch (nspin)
                    {
                    case 1:
                        occ_mat[iat][l][n][0] += transpose(occ_mat[iat][l][n][0]);
                        occ_mat[iat][l][n][0] *= 0.5;
                        occ_mat[iat][l][n][1] += occ_mat[iat][l][n][0];
                        break;

                    case 2:
                        for (int is = 0; is < nspin; is++)
                            occ_mat[iat][l][n][is] += transpose(occ_mat[iat][l][n][is]);
                        break;

                    case 4:
                        occ_mat[iat][l][n][0] += transpose(occ_mat[iat][l][n][0]);
                        break;

                    default:
                        std::cout << "Not supported NSPIN parameter" << std::endl;
                        exit(0);
                    }
                } // end n
            } // end l
        } // end ia
    } // end it

    if(is_mixing_enabled() && is_occ_mat_initialized())
    {
        this->mix_occ_mat(ucell,mixing_beta);
    }

    mark_occ_mat_initialized();
    ModuleBase::timer::end("Plus_U", "cal_occ_mat_k");
    return;
}

void Plus_U::cal_occ_mat_gamma(const Parallel_Orbitals* pv,
                             const int iter,
                             const UnitCell &ucell,
                             const std::vector<std::vector<double>> &dm_gamma,
                             const double& mixing_beta,
                             hamilt::Hamilt<double>* p_ham,
                             const int nspin,
                             const int npol,
                             const int nlocal,
                             const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>& iatlnmipol2iwt,
                             const std::vector<int>& orbital_corr)
{
    ModuleBase::TITLE("Plus_U", "cal_occ_mat_gamma");
    ModuleBase::timer::start("Plus_U", "cal_occ_mat_gamma");
    this->copy_occ_mat(ucell);
    this->zero_occ_mat(ucell);

    //=================Part 1======================
    // call PBLAS routine to calculate the product of the S and density matrix
    char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double alpha = 1.0, beta = 0.0;

    std::vector<double> srho(pv->nloc);
    for (int is = 0; is < nspin; is++)
    {
        double* s_gamma_pointer = dynamic_cast<hamilt::HamiltLCAO<double, double>*>(p_ham)->getSk();

#ifdef __MPI
        ScalapackConnector::gemm(transN,
            transT,
            nlocal,
            nlocal,
            nlocal,
            alpha,
            s_gamma_pointer,
            one_int,
            one_int,
            &pv->desc[0],
            dm_gamma[is].data(),
            //dm_gamma[is].c,
            one_int,
            one_int,
            &pv->desc[0],
            beta,
            srho.data(),
            one_int,
            one_int,
            &pv->desc[0]);
#endif

        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
			if (!orbital_corr[it] != -1)
			{
				continue;
			}
			for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
					if (l != orbital_corr[it])
					{
						continue;
					}

                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
						if (n != 0)
						{
							continue;
						}

                        // Calculate the local occupation number matrix
                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < npol; ipol0++)
                            {
                                const int iwt0 = iatlnmipol2iwt[iat][l][n][m0][ipol0];
                                const int mu = pv->global2local_row(iwt0);
                                const int mu_prime = pv->global2local_col(iwt0);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                    {
                                        const int iwt1 = iatlnmipol2iwt[iat][l][n][m1][ipol1];
                                        const int nu = pv->global2local_col(iwt1);
                                        const int nu_prime = pv->global2local_row(iwt1);

                                        const int irc = nu * pv->nrow + mu;
                                        const int irc_prime = mu_prime * pv->nrow + nu_prime;

                                        if ((nu >= 0) && (mu >= 0))
                                        {
                                            int m0_all = m0 + (2 * l + 1) * ipol0;
                                            int m1_all = m0 + (2 * l + 1) * ipol1;

                                            occ_mat[iat][l][n][is](m0, m1) += srho[irc] / 4.0;
                                        }

                                        if ((nu_prime >= 0) && (mu_prime >= 0))
                                        {
                                            int m0_all = m0 + (2 * l + 1) * ipol0;
                                            int m1_all = m0 + (2 * l + 1) * ipol1;

                                            occ_mat[iat][l][n][is](m0, m1) += srho[irc_prime] / 4.0;
                                        }
                                    }
                                }
                            }
                        }

                        ModuleBase::matrix temp(occ_mat[iat][l][n][is]);

#ifdef __MPI
                        MPI_Allreduce(&temp(0, 0),
                                      &occ_mat[iat][l][n][is](0, 0),
                                      (2 * l + 1) * npol * (2 * l + 1) * npol,
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
#endif

                        // for the case spin independent calculation
                        switch (nspin)
                        {
                        case 1:
                            occ_mat[iat][l][n][0] += transpose(occ_mat[iat][l][n][0]);
                            occ_mat[iat][l][n][0] *= 0.5;
                            occ_mat[iat][l][n][1] += occ_mat[iat][l][n][0];
                            break;

                        case 2:
                            occ_mat[iat][l][n][is] += transpose(occ_mat[iat][l][n][is]);
                            break;

                        default:
                            std::cout << "Not supported NSPIN parameter" << std::endl;
                            exit(0);
                        }

                    } // end for(n)
                } // L
            } // ia
        } // it
    } // is

    if(is_mixing_enabled() && is_occ_mat_initialized())
    {
        this->mix_occ_mat(ucell,mixing_beta);
    }

    mark_occ_mat_initialized();
    ModuleBase::timer::end("Plus_U", "cal_occ_mat_gamma");
    return;
}
#endif
