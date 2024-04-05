#include "LCAO_hamilt.h"

#include "module_base/parallel_reduce.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"	//caoyu add 2021-07-26
#endif
#include "module_base/timer.h"

#ifdef __EXX
#include "LCAO_hamilt.hpp"
#endif

#include "sparse_format.h"

using namespace sparse_format;

LCAO_Hamilt::LCAO_Hamilt()
{
}

LCAO_Hamilt::~LCAO_Hamilt()
{
    if(GlobalV::test_deconstructor)
    {
        std::cout << " ~LCAO_Hamilt()" << std::endl;
    }
}


void LCAO_Hamilt::cal_STN_R_for_T(const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","cal_STN_R_for_T");

    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    double tmp=0.0;
    std::complex<double> tmpc=complex<double>(0.0,0.0);

    for(int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if(adj)
                {
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);

					Abfs::Vector3_Order<int> dR(
							GlobalC::GridD.getBox(ad).x, 
							GlobalC::GridD.getBox(ad).y, 
							GlobalC::GridD.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = this->LM->ParaV->global2local_row(iw1_all);

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->global2local_col(iw2_all);

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
                                tmp = this->LM->Hloc_fixedR[index];
                                if (std::abs(tmp) > sparse_threshold)
                                {
                                    this->LM->TR_sparse[dR][iw1_all][iw2_all] = tmp;
                                }
                            }
                            else
                            {
                                tmpc = this->LM->Hloc_fixedR_soc[index];
                                if(std::abs(tmpc) > sparse_threshold)
                                {
                                    this->LM->TR_soc_sparse[dR][iw1_all][iw2_all] = tmpc;
                                }
                            }

                            ++index;
                        }
                    }
                }
            }
        }
    }

    return;
}


void LCAO_Hamilt::cal_SR(const double &sparse_threshold, hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("LCAO_Hamilt","cal_SR");
    sparse_format::set_R_range(*this->LM);
    //cal_STN_R_sparse(current_spin, sparse_threshold);
    if(GlobalV::NSPIN!=4)
    {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao 
        = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham);
        this->cal_HContainer_sparse_d(0, sparse_threshold, *(p_ham_lcao->getSR()), this->LM->SR_sparse);
    }
    else
    {
        hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_ham_lcao 
        = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham);
        this->cal_HContainer_sparse_cd(0, sparse_threshold, *(p_ham_lcao->getSR()), this->LM->SR_soc_sparse);
    }
}


void LCAO_Hamilt::cal_TR(
		LCAO_gen_fixedH &gen_h,
		const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","cal_TR");
    
    //need to rebuild T(R)
    this->LM->Hloc_fixedR.resize(this->LM->ParaV->nnr);
    this->LM->zeros_HSR('T');

    gen_h.build_ST_new('T', 0, GlobalC::ucell, this->LM->Hloc_fixedR.data());

    sparse_format::set_R_range(*this->LM);
    this->cal_STN_R_sparse_for_T(sparse_threshold);

    return;
}


// in case there are elements smaller than the threshold
void LCAO_Hamilt::clear_zero_elements(
		const int &current_spin, 
		const double &sparse_threshold)
{
    if(GlobalV::NSPIN != 4)
    {
        for (auto &R_loop : this->LM->HR_sparse[current_spin])
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

        for (auto &R_loop : this->LM->SR_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

    }
    else
    {
        for (auto &R_loop : this->LM->HR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

        for (auto &R_loop : this->LM->SR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }
    }
}


void LCAO_Hamilt::destroy_all_HSR_sparse(void)
{
	this->LM->destroy_HS_R_sparse();
}

void LCAO_Hamilt::destroy_TR_sparse(void)
{
	this->LM->destroy_T_R_sparse();
}

void LCAO_Hamilt::destroy_dH_R_sparse(void)
{
	this->LM->destroy_dH_R_sparse();
}
