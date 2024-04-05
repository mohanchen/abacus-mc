#ifndef SPARSE_FORMAT_ST_H 
#define SPARSE_FORMAT_ST_H

#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

namespace sparse_format
{
	void cal_SR(
			const double &sparse_thr, 
			hamilt::Hamilt<std::complex<double>>* p_ham);

	void cal_TR(
			LCAO_gen_fixedH &gen_h,
			const double &sparse_thr);

    void cal_STN_R_for_T(const double &sparse_thr);
}

#endif
