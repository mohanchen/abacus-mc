#include "spar_hsr.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

void sparse_format::cal_HSR(
		const int &current_spin, 
		const double &sparse_threshold, 
		const int (&nmp)[3], 
		hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("sparse_format","cal_HSR");

    sparse_format::set_R_range(*this->LM);

    const int nspin = GlobalV::NSPIN;

    //cal_STN_R_sparse(current_spin, sparse_threshold);
    if(nspin==1 || nspin==2)
    {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao = 
        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham);

		this->cal_HContainer_sparse_d(current_spin, 
				sparse_threshold, 
				*(p_ham_lcao->getHR()), 
				this->LM->HR_sparse[current_spin]);

		this->cal_HContainer_sparse_d(current_spin, 
				sparse_threshold, 
				*(p_ham_lcao->getSR()), 
				this->LM->SR_sparse);
    }
	else if(nspin==4)
	{
		hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_ham_lcao = 
			dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham);

		this->cal_HContainer_sparse_cd(current_spin, 
				sparse_threshold, 
				*(p_ham_lcao->getHR()), 
				this->LM->HR_soc_sparse);

		this->cal_HContainer_sparse_cd(current_spin, 
				sparse_threshold, 
				*(p_ham_lcao->getSR()), 
				this->LM->SR_soc_sparse);
	}
	else
	{
        ModuleBase::WARNING_QUIT("cal_HSR","check the value of nspin.");       
    }

    // only old DFT+U method need to cal extra contribution to HR
    if (GlobalV::dft_plus_u == 2)
    {
        if(nspin==1 || nspin==2)
        {
            cal_HR_dftu(current_spin, sparse_threshold);
        }
	    else if(nspin==4)
        {
            cal_HR_dftu_soc(current_spin, sparse_threshold);
        }
		else
		{
			ModuleBase::WARNING_QUIT("cal_HSR","check the value of nspin.");       
		}
	}

#ifdef __EXX
#ifdef __MPI
    // if EXX is considered
    if( GlobalC::exx_info.info_global.cal_exx )
    {
		if(GlobalC::exx_info.info_ri.real_number)
		{
			this->cal_HR_exx_sparse(current_spin, sparse_threshold, nmp, *this->LM->Hexxd);
		}
		else
		{
			this->cal_HR_exx_sparse(current_spin, sparse_threshold, nmp, *this->LM->Hexxc);
		}
	}
#endif // __MPI
#endif // __EXX

    sparse_matrix::clear_zero_elements(current_spin, sparse_threshold);

    return;
}


void sparse_format::cal_HContainer_d(
		const int &current_spin, 
		const double &sparse_threshold, 
		const hamilt::HContainer<double>& hR, 
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target)
{
    ModuleBase::TITLE("sparse_format","cal_HContainer_d");

    const Parallel_Orbitals* paraV = this->LM->ParaV;
    auto row_indexes = paraV->get_indexes_row();
    auto col_indexes = paraV->get_indexes_col();
    for(int iap=0;iap<hR.size_atom_pairs();++iap)
    {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = paraV->atom_begin_row[atom_i];
        int start_j = paraV->atom_begin_col[atom_j];
        int row_size = paraV->get_row_size(atom_i);
        int col_size = paraV->get_col_size(atom_j);
        for(int iR=0;iR<hR.get_atom_pair(iap).get_R_size();++iR)
        {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            int* r_index = hR.get_atom_pair(iap).get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index[0], r_index[1], r_index[2]);
            for(int i=0;i<row_size;++i)
            {
                int mu = row_indexes[start_i+i];
                for(int j=0;j<col_size;++j)
                {
                    int nu = col_indexes[start_j+j];
                    const auto& value_tmp = matrix.get_value(i,j);
                    if(std::abs(value_tmp)>sparse_threshold)
                    {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}

void sparse_format::cal_HContainer_cd(
		const int &current_spin, 
		const double &sparse_threshold, 
		const hamilt::HContainer<std::complex<double>>& hR, 
		std::map<Abfs::Vector3_Order<int>, 
		std::map<size_t, std::map<size_t, std::complex<double>>>>& target)
{
    ModuleBase::TITLE("sparse_format","cal_HContainer_cd");

    const Parallel_Orbitals* paraV = this->LM->ParaV;
    auto row_indexes = paraV->get_indexes_row();
    auto col_indexes = paraV->get_indexes_col();
    for(int iap=0;iap<hR.size_atom_pairs();++iap)
    {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = paraV->atom_begin_row[atom_i];
        int start_j = paraV->atom_begin_col[atom_j];
        int row_size = paraV->get_row_size(atom_i);
        int col_size = paraV->get_col_size(atom_j);
        for(int iR=0;iR<hR.get_atom_pair(iap).get_R_size();++iR)
        {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            int* r_index = hR.get_atom_pair(iap).get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index[0], r_index[1], r_index[2]);
            for(int i=0;i<row_size;++i)
            {
                int mu = row_indexes[start_i+i];
                for(int j=0;j<col_size;++j)
                {
                    int nu = col_indexes[start_j+j];
                    const auto& value_tmp = matrix.get_value(i,j);
                    if(std::abs(value_tmp)>sparse_threshold)
                    {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}


// in case there are elements smaller than the threshold
void sparse_format::clear_zero_elements(
		const int &current_spin, 
		const double &sparse_threshold)
{
    ModuleBase::TITLE("sparse_format","clear_zero_elements");

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
                }// end while iter
            }// end row loop
        }// end R loop

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
                }// end while iter
            }// end row_loop
        }// end R_loop
    }

    return;
}
