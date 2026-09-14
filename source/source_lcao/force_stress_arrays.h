#ifndef FORCESTRESS_ARRAYS_H
#define FORCESTRESS_ARRAYS_H

#include <vector>

class ForceStressArrays
{
    public:

    ForceStressArrays(){};
    ~ForceStressArrays(){};

    //-----------------------------------------
    // force in LCAO
    // used in gamma only algorithm.
    //-----------------------------------------
    double* DSloc_x = nullptr;
    double* DSloc_y = nullptr;
    double* DSloc_z = nullptr;

    //-----------------------------------------
    // force in LCAO
    // used in k-points algorithm.
    //-----------------------------------------
    double* DSloc_Rx = nullptr;
    double* DSloc_Ry = nullptr;
    double* DSloc_Rz = nullptr;

    //-----------------------------------------
    // dT + part of dVNL
    // used in gamma only algorithm.
    //-----------------------------------------
    double* DHloc_fixed_x = nullptr;
    double* DHloc_fixed_y = nullptr;
    double* DHloc_fixed_z = nullptr;

    //-----------------------------------------
    // dT + part of dVNL
    // used in kpoint algorithm.
    //-----------------------------------------
    double* DHloc_fixedR_x = nullptr;
    double* DHloc_fixedR_y = nullptr;
    double* DHloc_fixedR_z = nullptr;

    //----------------------------------------
    // r_mu - r_nu
    //----------------------------------------

    std::vector<double> DH_r;//zhengdy added 2017-07

    std::vector<double> stvnl11;
    std::vector<double> stvnl12;
    std::vector<double> stvnl13;
    std::vector<double> stvnl22;
    std::vector<double> stvnl23;
    std::vector<double> stvnl33;

    std::vector<double> DSloc_11;
    std::vector<double> DSloc_12;
    std::vector<double> DSloc_13;
    std::vector<double> DSloc_22;
    std::vector<double> DSloc_23;
    std::vector<double> DSloc_33;

    std::vector<double> DHloc_fixed_11;
    std::vector<double> DHloc_fixed_12;
    std::vector<double> DHloc_fixed_13;
    std::vector<double> DHloc_fixed_22;
    std::vector<double> DHloc_fixed_23;
    std::vector<double> DHloc_fixed_33;

};

#endif
