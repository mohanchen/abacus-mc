#ifndef SPARSE_FORMAT_ST_H 
#define SPARSE_FORMAT_ST_H

#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

namespace sparse_format
{
	void cal_SR(
			std::set<Abfs::Vector3_Order<int>> &all_R_coor,
			std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> SR_sparse,
			std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> SR_soc_sparse,
			Grid_Driver &grid,
			const double &sparse_thr, 
			hamilt::Hamilt<std::complex<double>>* p_ham);

	void cal_TR(
			const UnitCell &ucell,
			Parallel_Orbitals &pv,
			LCAO_Matrix &lm,
			Grid_Driver &grid,
			LCAO_gen_fixedH &gen_h,
			const double &sparse_thr);

	void cal_STN_R_for_T(
			const UnitCell &ucell,
			Parallel_Orbitals &pv,
			LCAO_Matrix &lm,
			Grid_Driver &grid,
			const double &sparse_thr);
}

#endif
