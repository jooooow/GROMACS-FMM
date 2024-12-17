#include "fmm.h"
#include "gromacs/math/units.h"

zeta_fmm::FMM::FMM(FILE* f, const t_mdatoms& md)
{
    fplog = f;
    P = 10;
    max_depth = 3;
    images = 4;
    reg_size = 0.00f;

    set_data(md);

    printf("init FMM(P=%d, max_depth=%d, images=%d, reg_size=%.8f, USE_DOUBLE=%d, FMM_REG_USE_DW=%d)\n", P, max_depth, images, reg_size, (typeid(zeta_fmm::real) == typeid(double)) ,FMM_REG_USE_DW);
}

zeta_fmm::FMM::~FMM()
{

}

void zeta_fmm::FMM::set_data(const t_mdatoms& md)
{
    natoms = md.homenr;
    chargeA_ = md.chargeA;
    chargeB_ = md.chargeB;
    
    bodies.resize(natoms);
    pme_bodies.resize(natoms);
}

void zeta_fmm::FMM::set_bodies(gmx::ArrayRef<const gmx::RVec> coord)
{
    zeta_fmm::real Q = 0;
    zeta_fmm::vec3r dipole(0,0,0);
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
        Q += chargeA_[i];
        dipole += bodies[i].loc * bodies[i].q;
    }
}

void zeta_fmm::FMM::calculate(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy, Solver solver_type)
{
    calculate_fmm(coord, box, verbose, dummy);
}

void zeta_fmm::FMM::calculate_fmm(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)
{
    set_bodies(coord);
    cycle = box[XX][XX];
    box_r = cycle / 2;
    center = zeta_fmm::vec3r(box[XX][XX] / 2, box[YY][YY] / 2, box[ZZ][ZZ] / 2);
    if(verbose)
    {
        printf("[verbose] cycle = %.8f, box_r = %.8f, center = [%.8f, %.8f, %.8f]\n", cycle, box_r, center[0], center[1], center[2]);
    }
    zeta_fmm::GPUSolver gpu_solver(
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
    gpu_solver.solve();
    gpu_solver.dipole_correction(cycle);
    bodies = gpu_solver.get_result(true);
}

void zeta_fmm::FMM::calculate_ewald(gmx::ArrayRef<const gmx::RVec> coord, const matrix box, int verbose, int dummy)
{

}

zeta_fmm::real zeta_fmm::FMM::get_energy()
{
    zeta_fmm::real energy = 0;
    for(int i = 0; i < natoms; i++)
    {
        zeta_fmm::Body3& b = bodies[i];
        energy += b.p * b.q;
    }
    coulomb_energy = energy / 2;
    return energy / 2;
}

void zeta_fmm::FMM::add_force(gmx::ForceWithVirial* forcewithvirial)
{
    for(int i = 0; i < natoms; i++)
    {
        zeta_fmm::Body3& b = bodies[i];
        forcewithvirial->force_[i][0] += b.f[0];
        forcewithvirial->force_[i][1] += b.f[1];
        forcewithvirial->force_[i][2] += b.f[2];
    }
}

void zeta_fmm::FMM::check_force(gmx::ForceWithVirial* forcewithvirial)
{
    zeta_fmm::Body3& b = bodies[0];
    printf("%.8f,%.8f,%.8f\n", forcewithvirial->force_[0][0],forcewithvirial->force_[0][1],forcewithvirial->force_[0][2]);
}

void zeta_fmm::FMM::print_one()
{
    zeta_fmm::Body3& b = bodies[0];
    printf("%d (%.8f,%.8f,%.8f) %.8f\n", b.idx, b.f[0], b.f[1], b.f[2], b.p);
}

int zeta_fmm::FMM::check_force(gmx::ForceWithShiftForces* pme_result_fs, int print)
{
    zeta_fmm::real FDif = 0, FNrm = 0;
    const gmx::ArrayRef<gmx::RVec> fs = pme_result_fs->force();
    for(int i = 0; i < natoms; i++)
    {
        zeta_fmm::Body3& pb = pme_bodies[i];
        pb.f = zeta_fmm::vec3r(fs[i][0],fs[i][1],fs[i][2]);
    }

    int have_nan = 0;
    for(int i = 0; i < natoms; i++)
    {
        zeta_fmm::Body3& b = bodies[i];
        if(isnan(b.f[0]) || isnan(b.f[0]) || isnan(b.f[0]) || isnan(b.p))
        {
            have_nan = 1;
            break;
        }
        zeta_fmm::Body3& pb = pme_bodies[i];
        zeta_fmm::vec3r diff = (b.f - pb.f);
        FDif += (b.f - pb.f).norm();
        FNrm += pb.f.norm();
        if(i < print)
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
    coulomb_energy = get_energy();
    printf("L2f = %.8f\n", L2f);
    if(isnan(L2f) || have_nan)
    {
        for(int i = 0; i < natoms; i++)
        {
            zeta_fmm::Body3& b = bodies[i];
            zeta_fmm::Body3& pb = pme_bodies[i];
            zeta_fmm::vec3r diff = (b.f - pb.f);
            FDif += (b.f - pb.f).norm();
            FNrm += pb.f.norm();
            printf("%d (%.4f,%.4f,%.4f)%.4f, (%.4f,%.4f,%.4f), (%.4f,%.4f,%.4f)\n", b.idx, b.f[0], b.f[1], b.f[2], b.p, pb.f[0], pb.f[1], pb.f[2], diff[0], diff[1], diff[2]);
        }
        fprintf(fplog, "body info : \n");
        for(int i = 0; i < natoms;i++)
        {
            fprintf(fplog, "%d %.32f %.32f %.32f %.32f\n", i, bodies[i].loc[0], bodies[i].loc[1], bodies[i].loc[2], bodies[i].q);
        }
        return 0;
    }
    return 1;
}

void zeta_fmm::FMM::clear_fs(gmx::ForceWithShiftForces* pme_result_fs)
{
    const gmx::ArrayRef<gmx::RVec> fs = pme_result_fs->force();
    for(int i = 0; i < natoms; i++)
    {
        fs[i][0] = fs[i][1] = fs[i][2] = 0;
    }
}

void zeta_fmm::FMM::clear_mtsforce(gmx::ArrayRef<gmx::RVec> f)
{
    for(int i = 0; i < natoms; i++)
    {
        f[i][0] = f[i][1] = f[i][2] = 0;
    }
}

int zeta_fmm::FMM::check_force(gmx::ForceWithVirial* pme_result_fv, int print)
{
    zeta_fmm::real FDif = 0, FNrm = 0;
    const gmx::ArrayRef<gmx::RVec> fv = pme_result_fv->force_;
    for(int i = 0; i < natoms; i++)
    {
        zeta_fmm::Body3& b = bodies[i];
        zeta_fmm::Body3 pb;
        pb.f = zeta_fmm::vec3r(fv[i][0],fv[i][1],fv[i][2]);
        zeta_fmm::vec3r diff = b.f - pb.f;
        if(i < print)
        {
            printf("%d (%.6f,%.6f,%.6f), (%.6f,%.6f,%.6f), (%.6f,%.6f,%.6f)\n", b.idx, b.f[0], b.f[1], b.f[2], pb.f[0], pb.f[1], pb.f[2], diff[0], diff[1], diff[2]);
        }
    }

    return 1;
}
