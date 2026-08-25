#include "source_base/matrix.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_io/module_parameter/parameter.h"

namespace spinconstrain {

/**
 * @brief Calculate atomic magnetic moments using projector overlap (PW basis).
 *
 * @details For each k-point:
 *   1. Tabulate atomic projectors: set up |alpha_{l,m}> for each atom
 *   2. Compute becp = <alpha_{l,m}|psi_{k,i}> via overlap_proj_psi
 *   3. Decompose becp into magnetic moments via accumulate_Mi_from_becp
 *
 * The magnetic moment is computed as:
 *   Mi = sum_{k,i} w_{k,i} * <psi_{k,i}|P_at|sigma|psi_{k,i}>
 * where P_at is the atomic projector and sigma are the Pauli matrices.
 *
 * Finally, Mi is summed across all MPI k-pool ranks since each pool only
 * has a subset of k-points.
 */
template <>
void SpinConstrain<std::complex<double>>::cal_mi_pw()
{
    ModuleBase::TITLE("module_deltaspin", "cal_mi_pw");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "cal_mi_pw");

    this->zero_Mi();
    if(PARAM.inp.device == "cpu")
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
        // Loop over k-points to calculate Mi of sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_CPU>*>(this->psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->get_npol();
        for(int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik); // Set up atomic projectors for this k-point
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer); // Compute becp = <alpha|psi>
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_tot_nproj();
            this->accumulate_Mi_from_becp(becp, nkb, nbands, npol, ik,
                &this->pelec->wg(ik, 0), &onsite_p->get_nh(0));
        }
    }
#if ((defined __CUDA) || (defined __ROCM))
    else
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->get_npol();
        for(int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik);
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer);
            const std::complex<double>* becp = onsite_p->get_h_becp();
            int nkb = onsite_p->get_size_becp() / nbands / npol;
            this->accumulate_Mi_from_becp(becp, nkb, nbands, npol, ik,
                &this->pelec->wg(ik, 0), &onsite_p->get_nh(0));
        }
    }
#endif
    // MPI reduction: sum Mi across all k-pool ranks
    Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, GlobalV::NPROC_IN_POOL, &(this->Mi_[0][0]), 3 * this->Mi_.size());

    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mi_pw");
}

} // namespace spinconstrain
