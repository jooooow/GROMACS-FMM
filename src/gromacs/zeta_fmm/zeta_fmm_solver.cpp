#include "zeta_fmm_solver.h"
#include "type.h"

namespace zeta_fmm
{
ZetaFmmSolver::ZetaFmmSolver(FILE* f, const t_mdatoms& md) : coulomb_solver::BaseSolver(f, md)
{
    std::cout<<"ZetaFmmSolver() with f and md"<<std::endl;
}
void ZetaFmmSolver::execute()
{
    std::cout<<"ZetaFmmSolver::execute()"<<std::endl;
    std::cout<<"ZetaFmmSolver::natoms="<<natoms<<std::endl;
}
}

