
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

void sparse_matrix_h::cal_HContainer_d(
		const int &current_spin, 
		const double &sparse_threshold, 
		const hamilt::HContainer<double>& hR, 
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target)
{
    ModuleBase::TITLE("sparse_matrix_h","cal_HContainer_d");

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

void sparse_matrix_h::cal_HContainer_cd(
		const int &current_spin, 
		const double &sparse_threshold, 
		const hamilt::HContainer<std::complex<double>>& hR, 
		std::map<Abfs::Vector3_Order<int>, 
		std::map<size_t, std::map<size_t, std::complex<double>>>>& target)
{
    ModuleBase::TITLE("sparse_matrix_h","cal_HContainer_cd");

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
