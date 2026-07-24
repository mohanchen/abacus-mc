#ifndef CAL_UX_H
#define CAL_UX_H

#include "source_cell/unitcell.h"

namespace unitcell {

void cal_ux(UnitCell& ucell, const int nspin);

bool judge_parallel(double a[3], ModuleBase::Vector3<double> b);

}

#endif