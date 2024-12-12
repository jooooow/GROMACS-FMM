#pragma once

namespace coulomb_solver
{

class BaseSolver
{
public:
    virtual ~BaseSolver() = default;
    virtual void execute() = 0;
};

}