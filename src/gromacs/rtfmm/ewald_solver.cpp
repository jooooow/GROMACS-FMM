#include "ewald_solver.h"
#include "type.h"
#include "fmm.h"
#include "argument.h"
#include "gromacs/math/units.h"
#include "ewald.h"

namespace rtfmm
{
EwaldSolverWrapper::EwaldSolverWrapper(FILE* f, const t_mdatoms& md) : coulomb_solver::BaseSolver(f, md)
{
    std::cout<<"EwaldSolverWrapper() with f and md"<<std::endl;
}

void EwaldSolverWrapper::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)
{
    std::cout<<"EwaldSolverWrapper::execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)"<<std::endl;
    std::cout<<"natoms="<<natoms<<std::endl;

    char* fakeArgv[] = {""};
    Argument args(1,fakeArgv);
    args.P = 10;
    args.n = natoms;
    args.images = 4;
    
    args.cycle = box[XX][XX];
    args.r = args.cycle / 2;
    args.x = vec3r(box[XX][XX] / 2, box[YY][YY] / 2, box[ZZ][ZZ] / 2);

    args.show();

    if(bs.empty())
    {
        bs.resize(args.n);
    }
    for(int i = 0; i < natoms;i++)
    {
        bs[i].idx = i;
        bs[i].p = 0;
        bs[i].f = vec3r(0,0,0);
        for(int m = 0; m < 3; m++)
        {
            bs[i].x[m] = coord.at(i)[m];
        }
        bs[i].q = chargeA_[i];
    }
    
    
    std::cout<<"Ewald begin"<<std::endl;

    rtfmm::EwaldSolver ewald(bs, args);
    bs = ewald.solve();

    std::cout<<"Ewald end"<<std::endl;

    for(size_t i = 0; i < bs.size(); i++)
    {
        Body3& b = bs[i];
        b.p *= gmx::c_one4PiEps0;
        b.f *= -b.q * gmx::c_one4PiEps0;
    }
}

void EwaldSolverWrapper::add_force(gmx::ForceWithVirial* forcewithvirial)
{
    std::cout<<"EwaldSolverWrapper::add_force"<<std::endl;
    for(int i = 0; i < natoms; i++)
    {
        Body3& b = bs[i];
        forcewithvirial->force_[i][0] += b.f[0];
        forcewithvirial->force_[i][1] += b.f[1];
        forcewithvirial->force_[i][2] += b.f[2];
    }
}

void EwaldSolverWrapper::compare_forces(const gmx::ArrayRef<gmx::RVec>& fs)
{
    real FDif = 0, FNrm = 0;
    for(int i = 0; i < natoms; i++)
    {
        Body3& b = bs[i];
        if(isnan(b.f[0]) || isnan(b.f[0]) || isnan(b.f[0]) || isnan(b.p))
        {
            std::cerr<<"have nan!"<<std::endl;
            break;
        }
        Body3 pb;
        pb.f = vec3r(fs[i][0],fs[i][1],fs[i][2]);
        vec3r diff = (b.f - pb.f);
        FDif += (b.f - pb.f).norm();
        FNrm += pb.f.norm();
        if(i < 3)
        {
            printf("%d [%.6f,%.6f,%.6f], (%.6f,%.6f,%.6f), (%.6f,%.6f,%.6f), (%.6f,%.6f,%.6f)\n", 
                b.idx, 
                b.x[0], b.x[1], b.x[2], 
                b.f[0], b.f[1], b.f[2], 
                pb.f[0], pb.f[1], pb.f[2], 
                diff[0], diff[1], diff[2]
            );
        }
    }
    real L2f = std::sqrt(FDif / FNrm);
    printf("L2f = %.8f\n", L2f);
}

coulomb_solver::Energy EwaldSolverWrapper::get_energy()
{
    coulomb_solver::Energy energy;
    std::cout<<"EwaldSolverWrapper::get_energy"<<std::endl;
    return energy;
}

}

