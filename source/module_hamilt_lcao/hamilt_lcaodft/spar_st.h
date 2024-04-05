#ifndef LCAO_HAMILT_H
#define LCAO_HAMILT_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"


#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

class LCAO_Hamilt
{
    public:

    LCAO_Hamilt();

    ~LCAO_Hamilt();


	void cal_SR(
			const double &sparse_thr, 
			hamilt::Hamilt<std::complex<double>>* p_ham);

	void cal_TR(
			LCAO_gen_fixedH &gen_h,
			const double &sparse_thr);

    LCAO_Matrix* LM;
};

#endif
