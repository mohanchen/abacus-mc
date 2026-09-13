#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_hamilt.h"
#include "dftu_nao_pots.h"
#include "source_base/global_function.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_ao/parallel_orbitals.h"

namespace DFTU_LCAO {

template <typename T>
void cal_pot_uterm(Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   const Parallel_Orbitals* pv,
                   const int spin,
                   T* pot_uterm,
                   const T* sk)
{
    ModuleBase::TITLE("DFTU_LCAO", "cal_pot_uterm");
    ModuleBase::timer::start("DFTU_LCAO", "cal_pot_uterm");

    const int nlocal = pv->get_global_row_size();
    ModuleBase::GlobalFunc::ZEROS(pot_uterm, pv->nloc);

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const T half = static_cast<T>(0.5);
    const T one = static_cast<T>(1.0);
    const T zero = static_cast<T>(0.0);

    std::vector<T> pot_onsite(pv->nloc);
    DFTU_LCAO::cal_pot_onsite(dftu, ucell, pv, spin, true, &pot_onsite[0]);

#ifdef __MPI
    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), one_int, one_int, pv->desc,
            sk, one_int, one_int, pv->desc,
            zero,
            pot_uterm, one_int, one_int, pv->desc);
#endif

    for (int irc = 0; irc < pv->nloc; irc++)
    {
        pot_onsite[irc] = pot_uterm[irc];
    }

#ifdef __MPI
    ScalapackConnector::tranu(nlocal, nlocal,
            one,
            &pot_onsite[0], one_int, one_int, pv->desc,
            one,
            pot_uterm, one_int, one_int, pv->desc);
#endif

    ModuleBase::timer::end("DFTU_LCAO", "cal_pot_uterm");
    return;
}

// Explicit instantiation
template void cal_pot_uterm<double>(Plus_U_Base& dftu,
                                    const UnitCell& ucell,
                                    const Parallel_Orbitals* pv,
                                    const int spin,
                                    double* pot_uterm,
                                    const double* sk);
template void cal_pot_uterm<std::complex<double>>(Plus_U_Base& dftu,
                                                  const UnitCell& ucell,
                                                  const Parallel_Orbitals* pv,
                                                  const int spin,
                                                  std::complex<double>* pot_uterm,
                                                  const std::complex<double>* sk);

} // namespace DFTU_LCAO
