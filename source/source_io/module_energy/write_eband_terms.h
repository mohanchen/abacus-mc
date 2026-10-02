#ifndef WRITE_EBAND_TERMS_H
#define WRITE_EBAND_TERMS_H

#include "source_psi/psi.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_basis/module_ao/parallel_orbitals.h"

#include <complex>
#include <map>
#include <vector>

#ifdef __EXX
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#endif

class Grid_Driver;

namespace ModuleIO
{

/// @brief Write the band-decomposed energy of each Hamiltonian term
/// (kinetic, local/nonlocal pp, Hartree, XC) in KS orbital representation.
/// @tparam TK K-point data type (double or std::complex<double>)
/// @tparam TR Real-space data type (double or std::complex<double>)
template <typename TK, typename TR>
void write_eband_terms(const int nspin,
                       const int nbasis,
                       const int drank,
                       const Parallel_Orbitals* pv,
                       const psi::Psi<TK>& psi,
                       const UnitCell& ucell,
                       Structure_Factor& sf,
                       surchem& solvent,
                       const ModulePW::PW_Basis& rho_basis,
                       const ModulePW::PW_Basis& rhod_basis,
                       const ModuleBase::matrix& vloc,
                       const Charge& chg,
                       const K_Vectors& kv,
                       const ModuleBase::matrix& wg,
                       Grid_Driver& gd,
                       const std::vector<double>& orb_cutoff,
                       const TwoCenterBundle& two_center_bundle,
                       const Exx_Info& exx_info
#ifdef __EXX
                       ,
                       std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>* Hexxd,
                       std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
);

} // namespace ModuleIO

#endif // WRITE_EBAND_TERMS_H
