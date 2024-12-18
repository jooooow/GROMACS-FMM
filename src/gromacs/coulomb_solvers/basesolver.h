#pragma once
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/forceoutput.h"
#include <iostream>

namespace coulomb_solver
{

class Energy
{
public: 
    Energy()
    {
        coulomb_energy = 0;
        lj_energy = 0;
        potential_energy = 0;
        kinetic_energy = 0;
        total_energy = 0;
        conserved_energy = 0;
    }

    double coulomb_energy;
    double lj_energy;
    double potential_energy;
    double kinetic_energy;
    double total_energy;
    double conserved_energy;

    friend std::ostream& operator<<(std::ostream& os, const Energy& energy) {
        os << "Energy Details = [" << std::endl;
        os << "\tCoulomb Energy       = " << energy.coulomb_energy << std::endl;
        os << "\tLennard-Jones Energy = " << energy.lj_energy << std::endl;
        os << "\tPotential Energy     = " << energy.potential_energy << std::endl;
        os << "\tKinetic Energy       = " << energy.kinetic_energy << std::endl;
        os << "\tTotal Energy         = " << energy.total_energy << std::endl;
        os << "\tConserved Energy     = " << energy.conserved_energy << std::endl;
        os << "]" << std::endl;
        return os;
    }
};

class BaseSolver
{
public:
    BaseSolver(FILE* f, const t_mdatoms& md)
    {
        std::cout<<"BaseSolver(FILE* f, const t_mdatoms& md)"<<std::endl;
        fplog = f;
        set_parameter(md);
        std::cout<<"BaseSolver set log and charge ok"<<std::endl;
    }

    virtual ~BaseSolver() = default;

    virtual void set_parameter(const t_mdatoms& md)
    {
        std::cout<<"BaseSolver::set_parameter(const t_mdatoms&)"<<std::endl;
        natoms = md.homenr;
        chargeA_ = md.chargeA;
        chargeB_ = md.chargeB;
    }

    virtual void execute(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose = 0, int dummy = 0) = 0;

    virtual Energy get_energy() = 0;

protected:
    FILE* fplog;
    int natoms;
    gmx::ArrayRef<const float> chargeA_; //! State A charge
    gmx::ArrayRef<const float> chargeB_; //! State B charge
};

}