#pragma once
#include <iostream>

#include "gromacs/coulomb_solvers/basesolver.h"

namespace zeta_fmm
{
class ZetaFmmSolver : public coulomb_solver::BaseSolver
{
public:
    ZetaFmmSolver(FILE* f, const t_mdatoms& md);
    void execute();
};
}

