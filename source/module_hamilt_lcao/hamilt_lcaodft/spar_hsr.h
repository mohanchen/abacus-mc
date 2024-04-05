

void LCAO_Hamilt::cal_HSR(
		const int &current_spin, 
		const double &sparse_threshold, 
		const int (&nmp)[3], 
		hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("sparse_matrix","cal_HSR");

    sparse_format::set_R_range(*this->LM);

    //cal_STN_R_sparse(current_spin, sparse_threshold);
    if(GlobalV::NSPIN!=4)
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
    else
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

    // only old DFT+U method need to cal extra contribution to HR
    if (GlobalV::dft_plus_u == 2)
    {
        if (GlobalV::NSPIN != 4)
        {
            cal_HR_dftu(current_spin, sparse_threshold);
        }
        else
        {
            cal_HR_dftu_soc(current_spin, sparse_threshold);
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

    clear_zero_elements(current_spin, sparse_threshold);

    return;
}

