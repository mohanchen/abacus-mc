#include "write_pdos_text.h"

#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"

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

        // Build the orbital list once: (w, atom_index, species, l, z, m)
        struct OrbitalInfo
        {
            int w;       // global orbital index (into pdos matrix row)
            int iat;     // 1-based atom index
            std::string species;
            int l;
            int z;       // zeta index (1-based)
            int m;
        };
        std::vector<OrbitalInfo> orb_list;

        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int ia = ucell.iat2ia[iat];
            const int it = ucell.iat2it[iat];
            const Atom* atom = &ucell.atoms[it];
            const int s0 = ucell.itiaiw2iwt(it, ia, 0);

            for (int j = 0; j < atom->nw; ++j)
            {
                OrbitalInfo oi;
                oi.iat = iat + 1;
                oi.species = ucell.atoms[it].label;
                oi.l = atom->iw2l[j];
                oi.z = atom->iw2n[j] + 1;
                oi.m = atom->iw2m[j];
                oi.w = ucell.itiaiw2iwt(it, ia, j);
                orb_list.push_back(oi);
            }
        }

        // Header: column-to-orbital mapping
        ofs << "# PDOS: energy(eV) in column 1, pdos(1/eV) in columns 2+"
            << std::endl;
        ofs << "# col  atom  species  l  z  m" << std::endl;
        for (size_t icol = 0; icol < orb_list.size(); ++icol)
        {
            ofs << "#" << std::setw(4) << icol + 2
                << std::setw(6) << orb_list[icol].iat
                << std::setw(8) << orb_list[icol].species
                << std::setw(3) << orb_list[icol].l
                << std::setw(3) << orb_list[icol].z
                << std::setw(3) << orb_list[icol].m
                << std::endl;
        }

        // Data: one row per energy point
        for (int n = 0; n < npoints; ++n)
        {
            const double en = emin + n * dos_edelta_ev;
            ofs << std::setw(14) << std::fixed << std::setprecision(6) << en;

            for (size_t icol = 0; icol < orb_list.size(); ++icol)
            {
                const int w = orb_list[icol].w;
                double pdos_val = 0.0;

                if (nspin == 4)
                {
                    // sum the two spinor components
                    const int iat = orb_list[icol].iat - 1;
                    const int ia = ucell.iat2ia[iat];
                    const int it = ucell.iat2it[iat];
                    const int s0 = ucell.itiaiw2iwt(it, ia, 0);
                    const int w0 = w - s0;
                    pdos_val = pdos[0](s0 + 2 * w0, n) + pdos[0](s0 + 2 * w0 + 1, n);
                }
                else
                {
                    pdos_val = pdos[is](w, n);
                }

                ofs << std::setw(16) << std::scientific << std::setprecision(6) << pdos_val;
            }
            ofs << std::endl;
        }

        ofs.close();
    }
}
