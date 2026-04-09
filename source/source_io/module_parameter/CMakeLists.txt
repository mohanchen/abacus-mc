
add_library(
    parameter
    OBJECT
    parameter.cpp
    dsp_config.cpp
)

# Link with base library for WARNING_QUIT
target_link_libraries(parameter PUBLIC base)

# Ensure device library is built before parameter
add_dependencies(parameter device)

# if(ENABLE_COVERAGE)
#   add_coverage(parameter)
# endif()

# if(BUILD_TESTING)
#   add_subdirectory(test)
# endif()