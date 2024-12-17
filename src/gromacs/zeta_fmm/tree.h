#pragma once
#include "type.h"
#include "body.h"
#include "cell.h"
#include <map>
#include <vector>

namespace zeta_fmm
{

class CompleteBalancedOctree
{
public:
    CompleteBalancedOctree();
    void build(std::vector<Body3>& bodies, zeta_fmm::real box_r, zeta_fmm::vec3r center, int max_depth);
    void build_reg(std::vector<Body3>& bodies, zeta_fmm::real rega, int max_depth, zeta_fmm::real cycle);
    std::vector<Cell3> get_cells();
    int get_max_depth();
private:
    std::vector<Cell3> cells;
    int num_body;
    int root_cell_idx;
    int reg0_cell_idx;

    void add_regcells(int max_depth, zeta_fmm::real cycle);
    void add_regbodies(std::vector<Body3>& bodies, zeta_fmm::real rega);
    void insert_regbodies_into_bodies(std::vector<Body3>& bodies, std::map<int, std::vector<zeta_fmm::Body3>>& reg_map);
};

}