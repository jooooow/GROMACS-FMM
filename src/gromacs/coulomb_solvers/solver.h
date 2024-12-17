#pragma once

#include <memory>
#include <iostream>
#include "gromacs/zeta_fmm/zeta_fmm_solver.h"
#include "gromacs/rtfmm/rtfmm_solver.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "basesolver.h"

namespace zeta_fmm {
    class ZetaFmmSolver;
}

namespace rtfmm {
    class rtFmmSolver;
}


namespace coulomb_solver
{

class SolverExecutor
{
private:
    std::shared_ptr<BaseSolver> solver;
public:
    enum class SolverType
    {
        ZETA,
        RT,
        EWALD    
    };

    void set_solver(SolverType solver_type, FILE* f, const t_mdatoms& md) 
    {
        std::cout<<"set solver : ";
        if(solver_type == SolverType::ZETA)
        {
            std::cout<<"ZETA"<<std::endl;
            solver = std::make_shared<zeta_fmm::ZetaFmmSolver>(f, md);
        }
        else if(solver_type == SolverType::RT)
        {
            std::cout<<"RT"<<std::endl;
            solver = std::make_shared<rtfmm::rtFmmSolver>(f, md);
        }
    }

    void execute() 
    {
        if (solver) 
        {
            solver->execute();
        } 
        else 
        {
            std::cerr << "No solver set!" << std::endl;
        }
    }
};

}