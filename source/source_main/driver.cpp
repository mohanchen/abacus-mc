#include "source_main/driver.h"

#include "source_base/global_file.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_esolver/esolver.h"
#include "source_io/module_output/cal_test.h"
#include "source_io/module_parameter/input_conv.h"
#include "source_io/module_json/para_json.h"
#include "source_io/module_output/print_info.h"
#include "source_io/module_parameter/read_input.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/version.h"
#include "source_base/parallel_global.h"
#ifdef __DSP
#include "source_base/module_device/memory_op.h"
#include "source_base/module_external/blas_connector.h"
#endif

Driver::Driver()
{
}

Driver::~Driver()
{
}

void Driver::init()
{
    // 1) Let's start by printing a title.
    ModuleBase::TITLE("Driver", "ABACUS_begins");

    // 2) Print the current time, since it may run a long time.
    time_t time_start = std::time(nullptr);
    ModuleBase::timer::start();

    // 3) Welcome to the atomic world! Let's do some fancy stuff here.
    this->atomic_world();

    // 4) All timers recorders are printed.
    ModuleBase::timer::finish(GlobalV::ofs_running);

    // 5) All memory recorders are printed.
    ModuleBase::Memory::print_all(GlobalV::ofs_running);

    // 6) Print the final time, hopefully it will not cost too long. 
    time_t time_finish = std::time(nullptr);
    ModuleIO::print_time(time_start, time_finish);

    // 7) Clean up: close all of the running logs
    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK, PARAM.inp.out_alllog,PARAM.inp.calculation);

}

void Driver::print_start_info()
{
    ModuleBase::TITLE("Driver", "print_start_info");
#ifdef VERSION
    const char* version = VERSION;
#else
    const char* version = "unknown";
#endif
#ifdef COMMIT_INFO
#include "commit.h"
    const char* commit = COMMIT;
#else
    const char* commit = "unknown";
#endif
    time_t time_now = time(nullptr);

    PARAM.set_start_time(time_now);
    GlobalV::ofs_running << "                                                  "
                            "                                   "
                         << std::endl;
    GlobalV::ofs_running << "                              ABACUS " << version << std::endl << std::endl;
    GlobalV::ofs_running << "               Atomic-orbital Based Ab-initio "
                            "Computation at UStc                    "
                         << std::endl
                         << std::endl;
    GlobalV::ofs_running << "                     Website: http://abacus.ustc.edu.cn/           "
                            "                  "
                         << std::endl;
    GlobalV::ofs_running << "               Documentation: https://abacus.deepmodeling.com/     "
                            "                  "
                         << std::endl;
    GlobalV::ofs_running << "                  Repository: "
                            "https://github.com/abacusmodeling/abacus-develop       "
                         << std::endl;
    GlobalV::ofs_running << "                              "
                            "https://github.com/deepmodeling/abacus-develop         "
                         << std::endl;
    GlobalV::ofs_running << "                      Commit: " << commit << std::endl << std::endl;
    GlobalV::ofs_running << std::setiosflags(std::ios::right);

#ifndef __MPI
    GlobalV::ofs_running << "    This is SERIES version." << std::endl;
#endif
    GlobalV::ofs_running << "    Start Time is " << ctime(&time_now);
    GlobalV::ofs_running << "                                                  "
                            "                                   "
                         << std::endl;
    GlobalV::ofs_running << " -------------------------------------------------"
                            "-----------------------------------"
                         << std::endl;

    GlobalV::ofs_running << std::setiosflags(std::ios::left);
    std::cout << std::setiosflags(std::ios::left);

    GlobalV::ofs_running << "\n READING GENERAL INFORMATION" << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "global_out_dir", PARAM.globalv.global_out_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "global_in_card",  PARAM.globalv.global_in_card);
}

void Driver::reading()
{
    ModuleBase::TITLE("Driver", "reading");
    ModuleBase::timer::start("Driver", "reading");
    // temperarily
    GlobalV::MY_RANK = PARAM.globalv.myrank;
    GlobalV::NPROC = PARAM.globalv.nproc;

    // (1) read the input file
    ModuleIO::ReadInput input(PARAM.globalv.myrank);
    input.read_parameters(PARAM, PARAM.globalv.global_in_card);

    ModuleBase::set_quit_out_dir(PARAM.globalv.global_out_dir);
    ModuleBase::set_quit_calculation(PARAM.inp.calculation);

#if defined(__CUDA) && defined(__USE_NVTX)
    ModuleBase::timer::set_nvtx_enabled(PARAM.inp.timer_enable_nvtx);
#endif

#ifdef __DSP
    if (PARAM.inp.dsp_count <= 0)
    {
        ModuleBase::WARNING_QUIT("driver", "dsp_count must be > 0");
    }
    base_device::memory::set_dsp_cluster_id(GlobalV::MY_RANK % PARAM.inp.dsp_count);
    BlasConnector::set_dsp_cluster_id(GlobalV::MY_RANK % PARAM.inp.dsp_count);
#endif

    // (2) create the output directory, running_*.log and print info
    input.create_directory(PARAM);
    this->print_start_info();

    // (3) write the input file
    std::stringstream ss1;
    ss1 << PARAM.globalv.global_out_dir <<  PARAM.globalv.global_in_card << ".info";
    input.write_parameters(PARAM, ss1.str());

    // (*temp*) copy the variables from INPUT to each class
    Input_Conv::Convert();

    // (4)+(5) Build the full parallel partition in one factory call.
    // The factory internally reproduces the exact legacy call order:
    //   1. split_diag_world  -> DIAG_WORLD
    //   2. split_grid_world -> GRID_WORLD
    //   3. divide_pools     -> POOL_WORLD/KP_WORLD/INT_BGROUP/BP_WORLD
    // and writes the same GlobalV scalars via the legacy helpers.
    // split_diag_world / split_grid_world write DRANK/DSIZE/DCOLOR/
    // GRANK/GSIZE through their reference parameters; divide_pools
    // writes NPROC_IN_POOL/RANK_IN_POOL/MY_POOL/etc through
    // init_pools' reference parameters into GlobalV. We then read
    // those scalars back for the OUT prints, keeping behavior identical.
    this->topo_ = Parallel_Global::create_partition(GlobalV::NPROC,
                                                     GlobalV::MY_RANK,
                                                     GlobalV::KPAR,
                                                     PARAM.inp.bndpar,
                                                     PARAM.inp.diago_proc,
                                                     PARAM.inp.diago_proc);

    // The factory's split_diag_world / split_grid_world write to local
    // variables, not to GlobalV. Read the values back from the topology
    // handles so the GlobalV scalars and printed output match the old
    // driver exactly.
    GlobalV::DRANK = 0;
    GlobalV::DSIZE = 1;
    GlobalV::DCOLOR = 0;
    GlobalV::GRANK = 0;
    GlobalV::GSIZE = 1;
#ifdef __MPI
    if (this->topo_.diag_world_comm() != MPI_COMM_NULL)
    {
        MPI_Comm_rank(this->topo_.diag_world_comm(), &GlobalV::DRANK);
        MPI_Comm_size(this->topo_.diag_world_comm(), &GlobalV::DSIZE);
    }
    if (this->topo_.rgrid_world_comm() != MPI_COMM_NULL)
    {
        MPI_Comm_rank(this->topo_.rgrid_world_comm(), &GlobalV::GRANK);
        MPI_Comm_size(this->topo_.rgrid_world_comm(), &GlobalV::GSIZE);
        GlobalV::DCOLOR = GlobalV::GRANK;
    }
#endif
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DRANK", GlobalV::DRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DSIZE", GlobalV::DSIZE);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "DCOLOR", GlobalV::DCOLOR + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GRANK", GlobalV::GRANK + 1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GSIZE", GlobalV::GSIZE);
    ModuleBase::timer::end("Driver", "reading");
}

void Driver::atomic_world()
{
    ModuleBase::TITLE("Driver", "atomic_world");
    ModuleBase::timer::start("Driver", "atomic_world");

    // reading information 
    this->reading();

    // where the actual stuff is done
    this->driver_run();

    ModuleBase::timer::end("Driver", "atomic_world");
}
