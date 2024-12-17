#include "solver.h"
#include <algorithm>
#include "gromacs/math/units.h"

zeta_fmm::Solver::Solver(FILE* f)
{
    fplog = f;
    //fprintf(fplog, "init Solver hjm\n");
}

zeta_fmm::Solver::Solver(FILE* f, std::vector<Body3>& bs, std::string n)
{
    fplog = f;
    //fprintf(fplog, "init Solver hjm\n");
    bodies = bs;
    for(size_t i = 0; i < bodies.size(); i++)
    {
        zeta_fmm::Body3& bi = bodies[i];
        bi.p = 0;
        bi.f = zeta_fmm::vec3r(0,0,0);
    }
    name = n;
}

zeta_fmm::Solver::~Solver()
{

}

void zeta_fmm::Solver::dipole_correction(zeta_fmm::real cycle)
{
    int num_body = bodies.size();
    zeta_fmm::real coef = 4 * M_PI / (3 * cycle * cycle * cycle);
    zeta_fmm::vec3r dipole(0,0,0);
    zeta_fmm::real Q = 0;
    for (int i = 0; i < num_body; i++) 
    {
        Q += bodies[i].q;
        dipole += bodies[i].loc * bodies[i].q;
    }
    //printf("Q = %.20f\n",Q);
    zeta_fmm::real dnorm = dipole.norm();
    for (int i = 0; i < num_body; i++) 
    { 
        //bodies[i].p -= - bodies[i].q * Q * gmx::c_one4PiEps0; 
        bodies[i].p -= coef * dnorm / num_body / bodies[i].q * gmx::c_one4PiEps0; 
        bodies[i].f += coef * dipole * bodies[i].q * gmx::c_one4PiEps0;
    }
}

std::vector<zeta_fmm::Body3> zeta_fmm::Solver::get_result(bool is_sort)
{
    if(is_sort)
        sort_bodies();
    return bodies;
}

void zeta_fmm::Solver::sort_bodies()
{
    std::sort(
        bodies.begin(), 
        bodies.end(),
        [](const Body3& a, const Body3& b)
        {
            return a.idx < b.idx; 
        }
    );
}