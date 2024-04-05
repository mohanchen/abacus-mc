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


    void cal_STN_R_for_T(const double &sparse_thr);


#ifdef __EXX
    template<typename Tdata> void cal_HR_exx_sparse(
            const int &current_spin,
            const double &sparse_thr,
            const int (&nmp)[3],
            const std::vector< std::map <int, std::map < std::pair<int, std::array<int,3>>, RI::Tensor<Tdata> > >>& Hexxs);
#endif


    void cal_SR(const double &sparse_thr, hamilt::Hamilt<std::complex<double>>* p_ham);

    void cal_TR(
			LCAO_gen_fixedH &gen_h,
			const double &sparse_thr);


    void clear_zero_elements(const int &current_spin, const double &sparse_thr);


    void destroy_all_HSR(void);

    void destroy_TR(void);

    void destroy_dH_R(void);

    LCAO_Matrix* LM;
};

#endif
