#pragma once
#include <iostream>

#include "gromacs/coulomb_solvers/basesolver.h"

namespace zeta_fmm
{
class ZetaFmmSolver : public coulomb_solver::BaseSolver
{
public:
    ZetaFmmSolver(FILE* f, const t_mdatoms& md);

    void execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy) override;

    coulomb_solver::Energy get_energy() override;
};
}

