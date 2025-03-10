#ifndef MAGNETISM_H
#define MAGNETISM_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "module_elecstate/module_charge/charge.h"

class Magnetism
{
public:
    // constructor and deconstructor
    Magnetism();
    ~Magnetism();

    // notice : bcast (MPI operation) is done in unitcell
    double *start_magnetization=nullptr;

    // tot_magnetization : majority spin - minority spin (nelup - neldw).
    double tot_magnetization;
    double tot_magnetization_nc[3]={0.0};
    double abs_magnetization;

    void compute_magnetization(const double& omega,
                               const int& nrxx, 
                               const int& nxyz, 
                               const double* const * rho, 
                               
                               double* nelec_spin = nullptr);

    ModuleBase::Vector3<double> *m_loc_=nullptr; //magnetization for each element along c-axis
	double *angle1_=nullptr;                     //angle between c-axis and real spin std::vector
	double *angle2_=nullptr;                     //angle between a-axis and real spin std::vector projection in ab-plane
    double ux_[3]={0.0};
	bool lsign_=false;

private:
    bool judge_parallel(double a[3],ModuleBase::Vector3<double> b);

};

/*
 A comment about variables nelup, neldw, multiplicity and tot_magnetization:
 All these variables contain the same information and must be kept harmonized.
 Variables nelup and neldw will be removed in future versions of the code.
 Variables multiplicity and tot_magnetization, though redundent will probably
 coexist since multiplicity is the more natural way (?)for defining the spin
 configuratio in the quantum-chemistry community while tot_magnetization is
 more natural (?) when dealing with extended systems.
*/

#endif
