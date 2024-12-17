#pragma once 
#include "type.h"
#include "body.h"
#include <vector>

namespace zeta_fmm
{
struct Cell3
{
    int depth;
    real r;
    vec3r center;
    OffsetAndNumber child_info;
    OffsetAndNumber body_info;
    
    //std::vector<Body3*> reg_body;
    //std::vector<complexr> M;
    //std::vector<complexr> L;
    friend std::ostream &operator<<(std::ostream & os, const Cell3 & cell) 
    {
        os << "---[cell]--- "
        <<"r=" << cell.r
        << ",depth=" << cell.depth
        <<",center=" << cell.center
        << ",child_info=" << cell.child_info
        << ",body_info=" << cell.body_info;
        return os;
    }
};

}