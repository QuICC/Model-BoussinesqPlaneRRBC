#include(ConfigureBenchmarks)
#quicc_add_benchmark(Explicit
#  MODEL "BoussinesqPlaneRRBC"
#  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Benchmarks"
#  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
#  )

include(ConfigureStability)
quicc_add_stability(Explicit
  MODEL "BoussinesqPlaneRRBC"
  WORKDIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/TestSuite/Stability"
  ARCHIVEDIR "${CMAKE_BINARY_DIR}/Models/_refdata"
  TIMEOUT 500
  VARIANTS "framework/parallel/algorithm:tubular"
  MPIRANKS 1
  GITTAG "new/rrbc_stability"
  )
