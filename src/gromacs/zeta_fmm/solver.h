#pragma once
#include <vector>
#include <cstdio>
#include "type.h"
#include "body.h"
#include "traverser.h"
#include "kernel.h"

namespace zeta_fmm
{
class Solver
{
public:
    Solver(FILE* f);
    Solver(FILE* f, std::vector<Body3>& bs, std::string n);
    virtual ~Solver();
    virtual void solve() = 0;
    void dipole_correction(zeta_fmm::real cycle);
    std::vector<zeta_fmm::Body3> get_result(bool is_sort = 0);
    void sort_bodies();
protected:
    FILE* fplog;
    std::vector<Body3> bodies;
    std::string name;
};

class GPUSolver : public Solver
{
public:
    GPUSolver(FILE* f);
    GPUSolver(
        FILE* f,
        std::vector<Body3>& bs,
        int P_,
        zeta_fmm::real box_r_,
        zeta_fmm::vec3r center_,
        int max_depth_,
        int images_,
        zeta_fmm::real cycle_,
        zeta_fmm::real rega_ = 0.0,
        int verbose_ = 0,
        int dummy_ = 0
    );
    ~GPUSolver();
    virtual void solve();
private:
    std::vector<int> get_leaf_list();
    int get_real_max_depth();
    void prepare_maps();
private:
    GPUKernel gpu_kernel;
    int P;
    zeta_fmm::real box_r;
    zeta_fmm::vec3r center;
    int images;
    zeta_fmm::real cycle;
    zeta_fmm::real rega;
    int max_depth;
    std::vector<Cell3> cells;
    Traverser::CellMap depth_map;
    Traverser::CellMap p2p_map;
    Traverser::CellMap m2l_map;
    int merged_num_body;
    int verbose;
    int dummy;

    std::vector<Body3> merged_bodies;
};
}