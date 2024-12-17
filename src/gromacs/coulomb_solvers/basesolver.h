#pragma once
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/forceoutput.h"

namespace coulomb_solver
{

class BaseSolver
{
public:
    BaseSolver(FILE* f, const t_mdatoms& md)
    {
        std::cout<<"BaseSolver(FILE* f, const t_mdatoms& md)"<<std::endl;
        fplog = f;
        natoms = md.homenr;
        chargeA_ = md.chargeA;
        chargeB_ = md.chargeB;
        std::cout<<"BaseSolver set log and charge ok"<<std::endl;
    }
    virtual ~BaseSolver() = default;
    virtual void execute() = 0;

protected:
    FILE* fplog;
    int natoms;
    gmx::ArrayRef<const float> chargeA_; //! State A charge
    gmx::ArrayRef<const float> chargeB_; //! State B charge
};

}