add_library(rtfmm INTERFACE)

file(GLOB RT_SOURCES
    rtfmm_solver.cpp
)

target_include_directories(rtfmm INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(rtfmm INTERFACE
                      legacy_api
                      gpu_utils
                      math
                      )



target_link_libraries(mdtypes INTERFACE utility)

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${RT_SOURCES} PARENT_SCOPE)
