#include "spar_st.h"

void sparse_format::cal_SR(
		std::set<Abfs::Vector3_Order<int>> &all_R_coor,
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> SR_sparse;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> SR_soc_sparse;
		const double &sparse_thr, 
		hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("sparse_format","cal_SR");

    sparse_format::set_R_range(all_R_coor, grid);

    //cal_STN_R_sparse(current_spin, sparse_thr);
    if(nspin==1 || nspin==2)
    {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao 
        = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham);
        this->cal_HContainer_d(0, sparse_thr, *(p_ham_lcao->getSR()), SR_sparse);
    }
    else if(nspin==4)
    {
        hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_ham_lcao 
        = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham);
        this->cal_HContainer_cd(0, sparse_thr, *(p_ham_lcao->getSR()), SR_soc_sparse);
    }

    return;
}


void LCAO_Hamilt::cal_TR(
		LCAO_gen_fixedH &gen_h,
		const double &sparse_thr)
{
    ModuleBase::TITLE("sparse_format","cal_TR");
    
    //need to rebuild T(R)
    this->LM->Hloc_fixedR.resize(this->LM->ParaV->nnr);
    this->LM->zeros_HSR('T');

    gen_h.build_ST_new('T', 0, GlobalC::ucell, this->LM->Hloc_fixedR.data());

    sparse_format::set_R_range(all_R_coor, grid);

    this->cal_STN_R_for_T(sparse_thr);

    return;
}
