#ifndef LCAO_SETUP_DM_H
#define LCAO_SETUP_DM_H

#include "source_estate/module_dm/setup_dm.h"

class K_Vectors;
class Parallel_Orbitals;

namespace LCAO_domain
{
// Backward-compatible alias: Setup_DM now lives in module_dm (source_estate).
using module_dm::Setup_DM;

// Free function replacing the former member Setup_DM::allocate_dm.
template <typename TK>
void allocate_dm(module_dm::Setup_DM<TK>& dmat,
                 const K_Vectors* kv,
                 const Parallel_Orbitals* pv,
                 const int nspin);
} // namespace LCAO_domain

#endif
