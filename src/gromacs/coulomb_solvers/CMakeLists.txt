add_library(coulomb_solvers INTERFACE)

#file(GLOB SOLVER_SOURCES solver.cpp)

target_include_directories(coulomb_solvers INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

                           

target_link_libraries(coulomb_solvers INTERFACE
                      legacy_api
                      gpu_utils
                      math
                      )

target_link_libraries(mdtypes INTERFACE utility)

#set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${SOLVER_SOURCES} PARENT_SCOPE)
