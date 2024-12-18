#pragma once
#include <iostream>

#include "gromacs/coulomb_solvers/basesolver.h"
#include "solver.h"
#include "type.h"

namespace zeta_fmm
{
class ZetaFmmSolver : public coulomb_solver::BaseSolver
{
public:
    ZetaFmmSolver(FILE* f, const t_mdatoms& md);

    void execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy) override;

    void add_force(gmx::ForceWithVirial* forcewithvirial) override;

    coulomb_solver::Energy get_energy() override;

    void compare_forces(const gmx::ArrayRef<gmx::RVec>& fs) override;

private:
    int P;
    int max_depth;
    int images;
    real box_r;
    vec3r center;
    real cycle;
    real reg_size;
    std::vector<Body3> bodies;
};
}

