#ifndef DRIVER_H
#define DRIVER_H

#include "source_base/parallel_partition.h"

class Driver
{
  public:
    Driver();
    ~Driver();

    /**
     * @brief An interface function.
     * This function calls "this->reading()", "this->atomic_world()" in order.
     */
    void init();

  private:
    /**
     * @brief Print the start information.
     *
     */
    void print_start_info();
    /**
     * @brief Reading the parameters and split the MPI world.
     * This function read the parameter in "INPUT", "STRU" etc,
     * and split the MPI world into different groups.
     */
    void reading();

    /**
     * @brief An interface function.
     * This function calls "this->driver_run()" to do calculation,
     * and log the time and  memory consumed during calculation.
     */
    void atomic_world();

    // the actual calculations
    void driver_run();

    // Init harewares according to Input parameters
    void init_hardware();
    void finalize_hardware();

    /// Parallel partition built once in reading() and injected into
    /// the ESolver in driver_run().
    ParallelPartition topo_;
};

#endif
