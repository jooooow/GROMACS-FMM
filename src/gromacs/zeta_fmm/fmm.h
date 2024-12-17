#ifndef ZETA_FMM_H
#define ZETA_FMM_H

#include <cstdio>
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "solver.h"
#include "type.h"

namespace zeta_fmm
{
    /*
    Wrapper of GPUSolver
    */
class FMM
{
public:

    enum class Solver
    {
        FMM
    };

    FMM(FILE* f, const t_mdatoms& md);
    ~FMM();
    void set_data(const t_mdatoms& md);
    void calculate(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose = 0, int dummy = 0, Solver solver_type = Solver::FMM);
    void calculate_fmm(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose = 0, int dummy = 0);
    void calculate_ewald(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose = 0, int dummy = 0);
    zeta_fmm::real get_energy();
    
    void add_force(gmx::ForceWithVirial* forcewithvirial);
    int check_force(gmx::ForceWithShiftForces* pme_result_fs, int print = 0);
    int check_force(gmx::ForceWithVirial* pme_result_fv, int print = 0);
    void check_force(gmx::ForceWithVirial* forcewithvirial);
    void clear_fs(gmx::ForceWithShiftForces* pme_result_fs);
    void clear_mtsforce(gmx::ArrayRef<gmx::RVec> f);
    void print_one();

    zeta_fmm::real coulomb_energy;
    zeta_fmm::real lj_energy;
    zeta_fmm::real potential_energy;
    zeta_fmm::real kinetic_energy;
    zeta_fmm::real total_energy;
    zeta_fmm::real conserved_energy;
private:
    FILE* fplog;
    int natoms;

    //! State A charge
    gmx::ArrayRef<const float> chargeA_;
    //! State B charge
    gmx::ArrayRef<const float> chargeB_;

    int P;
    int max_depth;
    int images;
    zeta_fmm::real box_r;
    zeta_fmm::vec3r center;
    zeta_fmm::real cycle;
    zeta_fmm::real reg_size;

    std::vector<Body3> bodies;
    std::vector<Body3> pme_bodies;

    void set_bodies(gmx::ArrayRef<const gmx::RVec> coord);
};
}

#endif