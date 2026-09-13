#include "source_pw/module_pwdft/dftu_base.h"
#include "dftu_hamilt.h"
#include "dftu_nao_pots.h"
#include "source_base/global_function.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_ao/parallel_orbitals.h"

namespace DFTU_LCAO {

void pot_uterm_complex(Plus_U_Base& dftu,
                       const UnitCell& ucell,
                       const Parallel_Orbitals* pv,
                       const int ik,
                       std::complex<double>* pot_uterm,
                       const std::vector<int>& isk,
                       const std::complex<double>* sk)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_uterm_complex");
    ModuleBase::timer::start("DFTU_LCAO", "pot_uterm_complex");

    int spin = isk[ik];

    const int nlocal = pv->get_global_row_size();
    ModuleBase::GlobalFunc::ZEROS(pot_uterm, pv->nloc);

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const std::complex<double> one(1.0, 0.0);
    const std::complex<double> half = 0.5;
    const std::complex<double> zero = 0.0;

    std::vector<std::complex<double>> pot_onsite(pv->nloc);
    DFTU_LCAO::pot_onsite_complex(dftu, ucell, pv, spin, true, &pot_onsite[0]);

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

    ModuleBase::timer::end("DFTU_LCAO", "pot_uterm_complex");
    return;
}

void pot_uterm_real(Plus_U_Base& dftu,
                    const UnitCell& ucell,
                    const Parallel_Orbitals* pv,
                    const int ik,
                    double* pot_uterm,
                    const std::vector<int>& isk,
                    const double* sk)
{
    ModuleBase::TITLE("DFTU_LCAO", "pot_uterm_real");
    ModuleBase::timer::start("DFTU_LCAO", "pot_uterm_real");

    int spin = isk[ik];

    const int nlocal = pv->get_global_row_size();
    ModuleBase::GlobalFunc::ZEROS(pot_uterm, pv->nloc);

    //=============================================================
    //   PART2: call pblas to calculate effective potential matrix
    //=============================================================
    const char transN = 'N', transT = 'T';
    int one_int = 1;
    double alpha = 1.0, beta = 0.0, half = 0.5, one = 1.0;

    std::vector<double> pot_onsite(pv->nloc);
    DFTU_LCAO::pot_onsite_real(dftu, ucell, pv, spin, true, &pot_onsite[0]);

#ifdef __MPI
    ScalapackConnector::gemm(transN, transN,
            nlocal, nlocal, nlocal,
            half,
            ModuleBase::GlobalFunc::VECTOR_TO_PTR(pot_onsite), 1, 1, pv->desc,
            sk, 1, 1, pv->desc,
            beta,
            pot_uterm, 1, 1, pv->desc);
#endif

    for (int irc = 0; irc < pv->nloc; irc++)
        pot_onsite[irc] = pot_uterm[irc];

#ifdef __MPI
    pdtran_(&nlocal, &nlocal,
            &one,
            &pot_onsite[0], &one_int, &one_int, const_cast<int*>(pv->desc),
            &one,
            pot_uterm, &one_int, &one_int, const_cast<int*>(pv->desc));
#endif

    ModuleBase::timer::end("DFTU_LCAO", "pot_uterm_real");
    return;
}

} // namespace DFTU_LCAO
