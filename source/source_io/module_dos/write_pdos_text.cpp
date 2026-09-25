#include "write_pdos_text.h"

#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"

#include <fstream>
#include <iomanip>
#include <sstream>

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

        ofs << "# energy(eV)  atom  species  l  z  m  pdos(1/eV)" << std::endl;

        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int ia = ucell.iat2ia[iat];
            const int it = ucell.iat2it[iat];
            const Atom* atom = &ucell.atoms[it];
            const int s0 = ucell.itiaiw2iwt(it, ia, 0);

            for (int j = 0; j < atom->nw; ++j)
            {
                const int L = atom->iw2l[j];
                const int N = atom->iw2n[j];
                const int m = atom->iw2m[j];
                const int w = ucell.itiaiw2iwt(it, ia, j);

                for (int n = 0; n < npoints; ++n)
                {
                    const double en = emin + n * dos_edelta_ev;
                    double pdos_val = 0.0;

                    if (nspin == 4)
                    {
                        // sum the two spinor components
                        const int w0 = w - s0;
                        pdos_val = pdos[0](s0 + 2 * w0, n) + pdos[0](s0 + 2 * w0 + 1, n);
                    }
                    else
                    {
                        pdos_val = pdos[is](w, n);
                    }

                    ofs << std::setw(12) << std::fixed << std::setprecision(6) << en
                        << std::setw(6) << iat + 1
                        << std::setw(8) << ucell.atoms[it].label
                        << std::setw(3) << L
                        << std::setw(3) << N + 1
                        << std::setw(3) << m
                        << std::setw(16) << std::scientific << std::setprecision(6) << pdos_val
                        << std::endl;
                }
            }
        }

        ofs.close();
    }
}
