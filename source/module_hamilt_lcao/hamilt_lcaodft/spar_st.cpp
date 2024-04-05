
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
