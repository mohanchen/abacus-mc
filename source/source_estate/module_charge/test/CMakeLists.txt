abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(_OPENMP)

if (ENABLE_MPI)

# Copy at configure time so a plain `make` + `ctest` run finds the data.
# install() only runs during `cmake --install`, which local test runs skip.
file(COPY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_CHARGE_rhog_io
  LIBS parameter base device planewave
  SOURCES test_chg_rhog_io.cpp ../chg_rhog_io.cpp ../../../source_basis/module_pw/test/test_tool.cpp
  # This test drives PW_Basis::initmpi and read/write_rhog's MPI collectives,
  # so it must keep __MPI even though this directory disables it. Its main()
  # calls MPI_Init via test_tool.cpp's setupmpi().
  KEEP_FEATURE_DEFINITIONS __MPI
)

endif()
