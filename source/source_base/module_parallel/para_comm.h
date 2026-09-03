#ifndef PARA_COMM_H
#define PARA_COMM_H

#include <complex>
#include <string>

#include "para_world.h"

namespace Parallel
{

// ========== Broadcast ==========

void bcast_bool(bool& v, const ParaWorld& world, int root = 0);
void bcast_int(int& v, const ParaWorld& world, int root = 0);
void bcast_double(double& v, const ParaWorld& world, int root = 0);
void bcast_complex(std::complex<double>& v, const ParaWorld& world, int root = 0);
void bcast_string(std::string& s, const ParaWorld& world, int root = 0);

void bcast_int(int* v, int n, const ParaWorld& world, int root = 0);
void bcast_double(double* v, int n, const ParaWorld& world, int root = 0);
void bcast_complex(std::complex<double>* v, int n, const ParaWorld& world, int root = 0);
void bcast_char(char* v, int n, const ParaWorld& world, int root = 0);
void bcast_string(std::string* v, int n, const ParaWorld& world, int root = 0);

// ========== Reduce (Allreduce, result on all ranks) ==========

void reduce_all(double& v, const ParaWorld& world);
void reduce_all(int& v, const ParaWorld& world);
void reduce_all(double* v, int n, const ParaWorld& world);

// ========== Reduce min/max ==========

void reduce_min(double& v, const ParaWorld& world);
void reduce_max(double& v, const ParaWorld& world);
void reduce_min(int& v, const ParaWorld& world);
void reduce_max(int& v, const ParaWorld& world);

// ========== Gather ==========

void gather_int(int& v, int* all, const ParaWorld& world);

} // namespace Parallel

#endif // PARA_COMM_H
