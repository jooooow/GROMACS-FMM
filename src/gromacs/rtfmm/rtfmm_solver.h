#pragma once
#include <iostream>

#include "gromacs/coulomb_solvers/basesolver.h"
#include "body.h"

namespace rtfmm
{
class rtFmmSolver : public coulomb_solver::BaseSolver
{
public:
    rtFmmSolver(FILE* f, const t_mdatoms& md);

    void execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy) override;

    void add_force(gmx::ForceWithVirial* forcewithvirial) override;
    
    void compare_forces(const gmx::ArrayRef<gmx::RVec>& fs) override;

    coulomb_solver::Energy get_energy() override;
    
private:
    Bodies3 bs;
};
}

