#include "body.h"

void zeta_fmm::print_one_body(Bodies3& bodies, int idx)
{
    zeta_fmm::Body3 b = bodies[idx];
    printf("---[body]--- idx=%d,loc=[%.8f,%.8f,%.8f],q=%.8f,p=%.10f,f=[%.10f,%.10f,%.10f]\n",b.idx,b.loc[0],b.loc[1],b.loc[2],b.q,b.p,b.f[0],b.f[1],b.f[2]);
}
