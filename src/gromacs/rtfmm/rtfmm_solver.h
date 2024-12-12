#pragma once
#include <iostream>

#include "gromacs/coulomb_solvers/basesolver.h"

namespace rtfmm
{
class rtFmmSolver : public coulomb_solver::BaseSolver
{
public:
    rtFmmSolver();
    void execute();
};
}

