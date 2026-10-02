
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_basis/module_ao/parallel_orbitals.h"

#include <tuple>
#include <vector>
#include <set>

// constructor of Atom
Atom::Atom()
{
}
Atom::~Atom()
{
}

Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}

Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}

pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}

SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}

Sep_Cell::Sep_Cell() noexcept : ntype(0), omega(0.0), tpiba2(0.0)
{
}
Sep_Cell::~Sep_Cell() noexcept = default;

// constructor of UnitCell
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}

Parallel_Orbitals::Parallel_Orbitals()
{
}
Parallel_Orbitals::~Parallel_Orbitals()
{
}

namespace RI_2D_Comm
{
int get_is_block(const int is_k, const int is_row_b, const int is_col_b)
{
    return 0;
}
std::tuple<int, int, int> get_iat_iw_is_block(const UnitCell& ucell, const int& iwt)
{
    return std::make_tuple(0, 0, 0);
}
std::tuple<int, int> split_is_block(const int is_b)
{
    return std::make_tuple(0, 0);
}
int get_iwt(const UnitCell& ucell, const int iat, const int iw_b, const int is_b)
{
    return 0;
}
std::vector<int> get_ik_list(const K_Vectors& kv, const int is_k)
{
    return {};
}
std::vector<std::tuple<std::set<int>, std::set<int>>> get_2D_judge(const UnitCell& ucell, const Parallel_2D& pv)
{
    return {};
}
} // namespace RI_2D_Comm
