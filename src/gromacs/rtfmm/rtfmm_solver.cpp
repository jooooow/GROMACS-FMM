#include "rtfmm_solver.h"
#include "type.h"

namespace rtfmm
{
rtFmmSolver::rtFmmSolver(FILE* f, const t_mdatoms& md) : coulomb_solver::BaseSolver(f, md)
{
    std::cout<<"rtFmmSolver() with f and md"<<std::endl;
}

void rtFmmSolver::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)
{
    std::cout<<"rtFmmSolver::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)"<<std::endl;
    std::cout<<"natoms="<<natoms<<std::endl;
    test();
}

coulomb_solver::Energy rtFmmSolver::get_energy()
{
    coulomb_solver::Energy energy;
    std::cout<<"rtFmmSolver::get_energy"<<std::endl;
    return energy;
}

}

