abacus_disable_feature_definitions(__MPI)
AddTest(
  TARGET MODULE_BASE_para_world
  SOURCES para_world_test.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_collection
  SOURCES para_collection_test.cpp ../para_collection.cpp ../para_world.cpp ../para_kmesh_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_kmesh_world
  SOURCES para_kmesh_world_test.cpp ../para_kmesh_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_pw_world
  SOURCES para_pw_world_test.cpp ../para_pw_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_diag_world
  SOURCES para_diag_world_test.cpp ../para_diag_world.cpp ../para_world.cpp
)

AddTest(
  TARGET MODULE_BASE_para_world_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_world_mpi_test.cpp ../para_world.cpp
)
target_compile_definitions(MODULE_BASE_para_world_mpi PRIVATE __MPI)

AddTest(
  TARGET MODULE_BASE_para_collection_mpi
  LIBS MPI::MPI_CXX
  SOURCES para_collection_mpi_test.cpp ../para_collection.cpp ../para_world.cpp ../para_kmesh_world.cpp
)
target_compile_definitions(MODULE_BASE_para_collection_mpi PRIVATE __MPI)

file(COPY para_world_mpi_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY para_collection_mpi_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME MODULE_BASE_para_world_mpi_test
      COMMAND ${BASH} para_world_mpi_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME MODULE_BASE_para_collection_mpi_test
      COMMAND ${BASH} para_collection_mpi_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
