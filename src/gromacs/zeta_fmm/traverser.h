#pragma once
#include <map>
#include "tree.h"

namespace zeta_fmm
{

class NearFieldSelector
{
public:
    enum class SelectorType
    {
        MAC, NEARALL, NEARNONE
    };
public:
    NearFieldSelector(SelectorType t, real s = 0);
    bool is_near(const Cell3& a, const Cell3& b, vec3r offset = vec3r(0,0,0));
private:
    SelectorType type;
    real scale;
};

class Traverser
{
public:
    using CellMap = std::map<int, std::vector<IndexAndOffset3r>>;

    Traverser();
    ~Traverser();

    void traverse(CompleteBalancedOctree& tree, NearFieldSelector& nearfield_selector, int images = 0, zeta_fmm::real cycle = 2 * M_PI);
    CellMap get_depth_map();
    CellMap get_p2p_map();
    CellMap get_m2l_map();

    std::vector<Cell3> get_cells();

private:
    void upward(
        int this_cell_idx = 0
    );
    void horizontal(
        NearFieldSelector& nearfield_selector,
        int this_cell_idx, 
        int that_cell_idx, 
        vec3r offset = vec3r(0,0,0)
    );
    void horizontal_periodic_near(
        NearFieldSelector& nearfield_selector,
        int this_cell_idx, 
        int that_cell_idx,
        zeta_fmm::real cycle
    );
    void horizontal_periodic_far(
        NearFieldSelector& nearfield_selector,
        int this_cell_idx, 
        int that_cell_idx,
        zeta_fmm::real cycle,
        int images
    );
private:
    CellMap depth_map;
    CellMap p2p_map;
    CellMap m2l_map;
    std::vector<zeta_fmm::Cell3> cells;
};


}
