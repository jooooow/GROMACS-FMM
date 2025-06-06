#pragma once

#include <memory>
#include <iostream>
#include "gromacs/zeta_fmm/zeta_fmm_solver.h"
#include "gromacs/rtfmm/rtfmm_solver.h"
#include "gromacs/rtfmm/ewald_solver.h"
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

class SolverExecuter
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
        else if(solver_type == SolverType::EWALD)
        {
            std::cout<<"EWALD"<<std::endl;
            solver = std::make_shared<rtfmm::EwaldSolverWrapper>(f, md);
        }
        else
        {
            std::cerr<<"invalid solver type !"<<static_cast<int>(solver_type)<<std::endl;
            exit(0);
        }
    }

    void execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy) 
    {
        if (solver) 
        {
            solver->execute(coord, box, verbose, dummy);
        } 
        else 
        {
            std::cerr << "No solver set!" << std::endl;
        }
    }

    void set_parameter(const t_mdatoms& md)
    {
        if (solver) 
        {
            solver->set_parameter(md);
        } 
        else 
        {
            std::cerr << "No solver set!" << std::endl;
        }
    }

    void add_force(gmx::ForceWithVirial* forcewithvirial)
    {
        if (solver)
        {
            solver->add_force(forcewithvirial);
        }
        else
        {
            std::cerr << "No solver set!" << std::endl;
        }
    }

    void compare_forces(const gmx::ArrayRef<gmx::RVec>& fs)
    {
        if (solver)
        {
            solver->compare_forces(fs);
        }
        else
        {
            std::cerr << "No solver set!" << std::endl;
        }
    }

    coulomb_solver::Energy get_energy()
    {
        if (solver)
        {
            return solver->get_energy();
        }
        else
        {
            std::cerr << "No solver set!" << std::endl;
        }

        return coulomb_solver::Energy();
    }
};

}