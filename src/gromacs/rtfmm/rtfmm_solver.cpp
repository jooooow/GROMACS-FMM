#include "rtfmm_solver.h"
#include "type.h"

namespace rtfmm
{
rtFmmSolver::rtFmmSolver(FILE* f, const t_mdatoms& md) : coulomb_solver::BaseSolver(f, md)
{
    std::cout<<"rtFmmSolver() with f and md"<<std::endl;
}
void rtFmmSolver::execute()
{
    std::cout<<"rtFmmSolver::execute()"<<std::endl;
    std::cout<<"natoms="<<natoms<<std::endl;
    test();
}
}

