#include "write_istate_info.h"

#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"

void ModuleIO::write_istate_info(const ModuleBase::matrix &ekb,const ModuleBase::matrix &wg, const K_Vectors& kv)
{
	ModuleBase::TITLE("ModuleIO","write_istate_info");
	ModuleBase::timer::tick("ModuleIO", "write_istate_info");

	std::stringstream ss;
    ss << PARAM.globalv.global_out_dir << "istate.info";
    if (GlobalV::MY_RANK == 0)
    {
        std::ofstream ofsi(ss.str().c_str()); // clear istate.info
        ofsi.close();
    }

    for (int ip = 0; ip < GlobalV::KPAR; ip++)
    {
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        if (GlobalV::MY_POOL == ip)
        {
            if (GlobalV::RANK_IN_POOL != 0 || GlobalV::MY_BNDGROUP != 0 ) { continue;
}
#endif
            std::ofstream ofsi2(ss.str().c_str(), std::ios::app);
            if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 4)
            {
                for (int ik = 0; ik < kv.get_nks(); ik++)
                {
#ifdef __MPI
                    int ik_global = kv.para_k.startk_pool[ip] + ik + 1;
#else
                    int ik_global = ik + 1;
#endif
                    ofsi2 << "BAND" << std::setw(25) << "Energy(ev)" << std::setw(25) << "Occupation"
                          << std::setw(25) << "Kpoint = " << ik_global
                          << std::setw(25) << "(" << kv.kvec_d[ik].x << " " << kv.kvec_d[ik].y
                          << " " << kv.kvec_d[ik].z << ")" << std::endl;
                    for (int ib = 0; ib < PARAM.globalv.nbands_l; ib++)
                    {
                        ofsi2.precision(16);
                        ofsi2 << std::setw(6) << ib + 1 << std::setw(25)
                              << ekb(ik, ib) * ModuleBase::Ry_to_eV << std::setw(25) << wg(ik, ib)
                              << std::endl;
                    }
                    ofsi2 << std::endl;
                    ofsi2 << std::endl;
                }
            }
            else
            {
                for (int ik = 0; ik < kv.get_nks() / 2; ik++)
                {
#ifdef __MPI
                    int ik_global = kv.para_k.startk_pool[ip] + ik + 1;
#else
                    int ik_global = ik + 1;
#endif
                    ofsi2 << "BAND" << std::setw(25) << "Spin up Energy(ev)" << std::setw(25) << "Occupation"
                          << std::setw(25) << "Spin down Energy(ev)" << std::setw(25) << "Occupation"
                          << std::setw(25) << "Kpoint = " << ik_global
                          << std::setw(25) << "(" << kv.kvec_d[ik].x << " " << kv.kvec_d[ik].y
                          << " " << kv.kvec_d[ik].z << ")" << std::endl;

                    for (int ib = 0; ib < PARAM.inp.nbands; ib++)
                    {
                        ofsi2 << std::setw(6) << ib + 1 << std::setw(25)
                              << ekb(ik, ib) * ModuleBase::Ry_to_eV << std::setw(25) << wg(ik, ib)
                              << std::setw(25) << ekb((ik + kv.get_nks() / 2), ib) * ModuleBase::Ry_to_eV
                              << std::setw(25) << wg(ik + kv.get_nks() / 2, ib) << std::endl;
                    }
                    ofsi2 << std::endl;
                    ofsi2 << std::endl;
                }
            }
            ofsi2.close();
#ifdef __MPI
        }
#endif
    }
	ModuleBase::timer::tick("ModuleIO", "write_istate_info");
	return;
}
