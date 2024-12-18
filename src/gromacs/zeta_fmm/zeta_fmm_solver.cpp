#include "zeta_fmm_solver.h"
#include "type.h"

namespace zeta_fmm
{
ZetaFmmSolver::ZetaFmmSolver(FILE* f, const t_mdatoms& md) : coulomb_solver::BaseSolver(f, md)
{
    std::cout<<"ZetaFmmSolver() with f and md"<<std::endl;
}

void ZetaFmmSolver::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)
{
    std::cout<<"ZetaFmmSolver::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)"<<std::endl;
    std::cout<<"ZetaFmmSolver::natoms="<<natoms<<std::endl;
}

coulomb_solver::Energy ZetaFmmSolver::get_energy()
{
    coulomb_solver::Energy energy;
    std::cout<<"ZetaFmmSolver::get_energy"<<std::endl;
    return energy;
}

}

