#include "zeta_fmm_solver.h"
#include "type.h"

namespace zeta_fmm
{
ZetaFmmSolver::ZetaFmmSolver(FILE* f, const t_mdatoms& md) : coulomb_solver::BaseSolver(f, md)
{
    std::cout<<"ZetaFmmSolver() with f and md"<<std::endl;

    P = 10;
    max_depth = 3;
    images = 4;
    reg_size = 0.00f;
}

void ZetaFmmSolver::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)
{
    std::cout<<"ZetaFmmSolver::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)"<<std::endl;
    std::cout<<"ZetaFmmSolver::natoms="<<natoms<<std::endl;

    if(bodies.empty())
    {
        std::cout<<"resize bodies"<<std::endl;
        bodies.resize(natoms);
    }
    for(int i = 0; i < natoms;i++)
    {
        bodies[i].idx = i;
        bodies[i].p = 0;
        bodies[i].f = zeta_fmm::vec3r(0,0,0);
        for(int m = 0; m < 3; m++)
        {
            bodies[i].loc[m] = coord.at(i)[m];
        }
        bodies[i].q = chargeA_[i];
    }

    cycle = box[XX][XX];
    box_r = cycle / 2;
    center = zeta_fmm::vec3r(box[XX][XX] / 2, box[YY][YY] / 2, box[ZZ][ZZ] / 2);
    if(verbose)
    {
        printf("[verbose] cycle = %.8f, box_r = %.8f, center = [%.8f, %.8f, %.8f]\n", cycle, box_r, center[0], center[1], center[2]);
    }
    std::cout<<"init GPUSolver"<<std::endl;
    GPUSolver gpu_solver(
        fplog,
        bodies,
        P,
        box_r,
        center,
        max_depth,
        images,
        cycle,
        reg_size,
        verbose,
        dummy
    );
    std::cout<<"gpu_solver.solve()"<<std::endl;
    gpu_solver.solve();
    gpu_solver.dipole_correction(cycle);
    bodies = gpu_solver.get_result(true);
    std::cout<<"finished"<<std::endl;
}

void ZetaFmmSolver::add_force(gmx::ForceWithVirial* forcewithvirial)
{
    std::cout<<"ZetaFmmSolver::add_force"<<std::endl;
    for(int i = 0; i < natoms; i++)
    {
        zeta_fmm::Body3& b = bodies[i];
        forcewithvirial->force_[i][0] += b.f[0];
        forcewithvirial->force_[i][1] += b.f[1];
        forcewithvirial->force_[i][2] += b.f[2];
    }
}

void ZetaFmmSolver::compare_forces(const gmx::ArrayRef<gmx::RVec>& fs)
{
    zeta_fmm::real FDif = 0, FNrm = 0;
    for(int i = 0; i < natoms; i++)
    {
        zeta_fmm::Body3& b = bodies[i];
        if(isnan(b.f[0]) || isnan(b.f[0]) || isnan(b.f[0]) || isnan(b.p))
        {
            std::cerr<<"have nan!"<<std::endl;
            break;
        }
        zeta_fmm::Body3 pb;
        pb.f = zeta_fmm::vec3r(fs[i][0],fs[i][1],fs[i][2]);
        zeta_fmm::vec3r diff = (b.f - pb.f);
        FDif += (b.f - pb.f).norm();
        FNrm += pb.f.norm();
        if(i < 3)
        {
            printf("%d [%.6f,%.6f,%.6f], (%.6f,%.6f,%.6f), (%.6f,%.6f,%.6f), (%.6f,%.6f,%.6f)\n", 
                b.idx, 
                b.loc[0], b.loc[1], b.loc[2], 
                b.f[0], b.f[1], b.f[2], 
                pb.f[0], pb.f[1], pb.f[2], 
                diff[0], diff[1], diff[2]
            );
        }
    }
    zeta_fmm::real L2f = std::sqrt(FDif / FNrm);
    printf("L2f = %.8f\n", L2f);
}

coulomb_solver::Energy ZetaFmmSolver::get_energy()
{
    coulomb_solver::Energy energy;
    std::cout<<"ZetaFmmSolver::get_energy"<<std::endl;
    return energy;
}

}

