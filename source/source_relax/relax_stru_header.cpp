#include "relax_driver.h"

#include "source_base/formatter.h"
#include "source_main/version.h"

#include <ctime>

std::string Relax_Driver::build_stru_header(const int istep,
                                            const double etot,
                                            const ModuleBase::matrix& stress,
                                            const Input_para& inp,
                                            const bool is_final,
                                            const bool geometry_evaluated)
{
    std::time_t now = std::time(nullptr);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    const char* final_tag = is_final ? " (FINAL)" : "";
    std::string header = FmtCore::format("# ABACUS version: %s\n# Written at %s\n# RELAX STEP %d%s, Energy: %.8f eV\n",
                                          VERSION,
                                          time_buf,
                                          istep + 1,
                                          final_tag,
                                          etot * ModuleBase::Ry_to_eV);

    // stress in kbar: Ry/Bohr^3 -> kbar, always 3 lines for downstream parsers.
    // N/A marks uncomputed stress; the reason is appended on the first line.
    const double stress_transform = ModuleBase::RYDBERG_SI
                                    / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI)
                                    * 1.0e-8;
    const bool stress_valid = inp.cal_stress && geometry_evaluated;
    if (stress_valid)
    {
        for (int i = 0; i < 3; i++)
        {
            header += FmtCore::format("# Stress (kbar): %.6f %.6f %.6f\n",
                                      stress(i, 0) * stress_transform,
                                      stress(i, 1) * stress_transform,
                                      stress(i, 2) * stress_transform);
        }
    }
    else
    {
        const char* reason = !inp.cal_stress ? " (N/A = not computed, cal_stress=0)"
                                             : " (N/A = proposed geometry not evaluated)";
        header += std::string("# Stress (kbar): N/A N/A N/A") + reason + "\n"
                  "# Stress (kbar): N/A N/A N/A\n"
                  "# Stress (kbar): N/A N/A N/A\n";
    }

    if (!geometry_evaluated)
    {
        header += "# NOTE: geometry proposed by optimizer but not evaluated; forces omitted, energy above belongs to the last evaluated geometry\n";
    }
    else if (!inp.cal_force)
    {
        header += "# Forces not computed (cal_force=0); per-atom f fields omitted intentionally\n";
    }

    return header;
}
