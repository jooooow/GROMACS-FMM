add_library(zeta_fmm INTERFACE)

file(GLOB ZETA_SOURCES
    zeta_fmm_solver.cpp
    fmm.cpp
    tree.cpp
    traverser.cpp
    solver.cpp  
)

target_include_directories(zeta_fmm INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(zeta_fmm INTERFACE
                      legacy_api
                      gpu_utils
                      math
                      )



target_link_libraries(mdtypes INTERFACE utility)

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${ZETA_SOURCES} PARENT_SCOPE)


if(GMX_GPU_CUDA)
gmx_add_libgromacs_sources(
    gpu_solver.cu
    body.cu
    kernel.cu
)
endif()
