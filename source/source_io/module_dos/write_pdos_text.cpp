#include "write_pdos_text.h"

#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

void ModuleIO::write_pdos_text(
        const UnitCell& ucell,
        const ModuleBase::matrix* pdos,
        const int nspin,
        const int nlocal,
        const int npoints,
        const double& emin,
        const double& dos_edelta_ev,
        const std::string& out_dir,
        const std::string& basis,
        const int istep)
{
    ModuleBase::TITLE("ModuleIO", "write_pdos_text");

    // number of spin channels to write
    const int nspin0 = (nspin == 2) ? 2 : 1;

    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;
        ss << out_dir << "pdoss" << is + 1;
        if (istep >= 0)
        {
            ss << "g" << istep + 1;
        }
        ss << "_" << basis << ".txt";

        std::ofstream ofs(ss.str().c_str());

        ofs << "# energy(eV)  atom  species  pdos(s,py,pz,px,dxy,dyz,dz2,dxz,dx2,f..., 1/eV)" << std::endl;

        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int ia = ucell.iat2ia[iat];
            const int it = ucell.iat2it[iat];
            const Atom* atom = &ucell.atoms[it];
            const int s0 = ucell.itiaiw2iwt(it, ia, 0);
            const int max_l = atom->nwl;

            // collect all (l, m) for this atom, ordered s,p,d,f
            struct OrbInfo
            {
                int l;
                int m;
            };
            std::vector<OrbInfo> orb_list;

            for (int L = 0; L <= max_l; ++L)
            {
                for (int m = -L; m <= L; ++m)
                {
                    for (int j = 0; j < atom->nw; ++j)
                    {
                        if (atom->iw2l[j] == L && atom->iw2m[j] == m)
                        {
                            orb_list.push_back({L, m});
                            break;
                        }
                    }
                }
            }

            for (int n = 0; n < npoints; ++n)
            {
                const double en = emin + n * dos_edelta_ev;

                ofs << std::setw(12) << std::fixed << std::setprecision(6) << en
                    << "  " << std::setw(3) << iat + 1
                    << "  " << std::setw(4) << ucell.atoms[it].label;

                for (size_t io = 0; io < orb_list.size(); ++io)
                {
                    const int L = orb_list[io].l;
                    const int m = orb_list[io].m;
                    double pdos_val = 0.0;

                    // sum over all zeta for this (l, m)
                    for (int j = 0; j < atom->nw; ++j)
                    {
                        if (atom->iw2l[j] != L || atom->iw2m[j] != m)
                        {
                            continue;
                        }
                        const int w = ucell.itiaiw2iwt(it, ia, j);

                        if (nspin == 4)
                        {
                            const int w0 = w - s0;
                            pdos_val += pdos[0](s0 + 2 * w0, n) + pdos[0](s0 + 2 * w0 + 1, n);
                        }
                        else
                        {
                            pdos_val += pdos[is](w, n);
                        }
                    }

                    // zero out negligible values
                    if (std::abs(pdos_val) < 1e-6)
                    {
                        pdos_val = 0.0;
                    }

                    ofs << "  " << std::setw(8) << std::fixed << std::setprecision(6) << pdos_val;
                }
                ofs << std::endl;
            }
        }

        ofs.close();
    }
}
