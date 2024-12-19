add_library(rtfmm INTERFACE)

file(GLOB RT_SOURCES
    rtfmm_solver.cpp
    body.cpp
    surface.cpp
    fmm.cpp
    kernel.cpp
    mathfunc.cpp
    tree.cpp
    traverser.cpp
    ewald.cpp
    argument.cpp
    epkifmm.cpp
    ewald_solver.cpp
)

find_package(BLAS REQUIRED)
message(STATUS BLAS found: ${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)
message(STATUS LAPACK found: ${LAPACK_LIBRARIES})

find_package(OpenMP REQUIRED)

target_include_directories(rtfmm INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(rtfmm INTERFACE
    legacy_api
    gpu_utils
    math
    fftw3
    OpenMP::OpenMP_CXX
    LAPACK::LAPACK
    BLAS::BLAS
)



target_link_libraries(mdtypes INTERFACE utility)

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${RT_SOURCES} PARENT_SCOPE)
