#pragma once 
#include "type.h"
#include <random>
#include <algorithm>

namespace zeta_fmm
{
struct Body3
{
    int idx;
    real q;
    real p;
    vec<3, real> loc;
    vec<3, real> f;

    Body3():q(0),p(0),loc(0),f(0){}

    friend std::ostream &operator<<(std::ostream & os, const Body3 & body) 
    {
        os << "---[body]--- "
        <<"idx=" << body.idx
        <<",loc=" << body.loc 
        << ",q=" << body.q 
        << ",p=" << body.p 
        << ",f=" << body.f;
        return os;
    }
};

typedef std::vector<Body3> Bodies3;

void print_one_body(Bodies3& bodies, int idx);

}

