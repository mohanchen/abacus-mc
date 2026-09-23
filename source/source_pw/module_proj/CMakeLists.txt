# module_proj: projection of plane-wave wavefunctions onto local (atomic) orbitals.
#
# This module hosts the onsite-projector code extracted from module_pwdft:
# radial projectors (SBT), onsite projector init/tabulate/overlap, and the
# Onsite_Proj_tools computation backend.
#
# Known temporary dependency (recorded during extraction):
# - Onsite_Proj_tools includes source_pw/module_pwdft/vnl_pw.h and
#   source_pw/module_pwdft/kernels/{force_op,stress_op}.h, i.e. module_proj
#   currently depends on module_pwdft headers. This is tolerated short-term;
#   the planned cleanup is to split the pseudopotential-specific constructor
#   of Onsite_Proj_tools so this reverse dependency can be removed.

list(APPEND objects
    radial_proj.cpp
)

add_library(
    module_proj
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(module_proj)
endif()

# Note: no test/ subdirectory yet. It will be added together with the first
# migrated test (radial_proj_test) when the corresponding sources move in.
