#include "ewald.h"
#include <omp.h>

//#define EWALD_USE_HALF

rtfmm::EwaldSolver::EwaldSolver(const Bodies3& tar_bs_, const Argument& args_) : EwaldSolver(tar_bs_, tar_bs_, args_)
{
    std::cout<<"EwaldSolver(const Bodies3& tar_bs_, const Argument& args_)"<<std::endl;
}

rtfmm::EwaldSolver::EwaldSolver(const Bodies3& tar_bs_, const Bodies3& src_bs_, const Argument& args_) 
: tar_bs(tar_bs_), src_bs(src_bs_), args(args_)
{
    ksize = args.ewald_ksize;
    if(verbose) printf("ksize = %d\n", ksize);
    scale = 2 * M_PI / args.cycle;
    if(verbose) printf("scale = %.4f\n", scale);
    alpha = ksize / args.cycle;
    cutoff = args.cycle / 2;
    #ifndef EWALD_USE_HALF
    for (int k = -ksize; k <= ksize; k++)
    {                               
        for (int m = -ksize; m <= ksize; m++)
        {                          
            for (int p = -ksize; p <= ksize; p++)
            {                        
                if(k == 0 && m == 0 && p == 0) continue;
                Wave wave;
                wave.K = rtfmm::vec3r(k,m,p);
                wave.val = std::complex<double>(0,0);
                //if(wave.K.norm() <= ksize * ksize)
                waves.push_back(wave);
            }
        }
    }
    #else
    int kmaxsq = ksize * ksize;
    int kmax = ksize;
    for (int l=0; l<=kmax; l++)
    {
        int mmin = -kmax;
        if (l==0) mmin = 0;
        for (int m=mmin; m<=kmax; m++)
        {
            int nmin = -kmax;
            if (l==0 && m==0) nmin=1;
            for (int n=nmin; n<=kmax; n++) 
            {
                real ksq = l * l + m * m + n * n;
                if (ksq <= kmaxsq) 
                {
                    Wave wave;
                    wave.K = vec3r(l,m,n);
                    wave.val = std::complex<double>(0,0);
                    waves.push_back(wave);
                }
            }
        }
    }
    #endif
}

rtfmm::Bodies3 rtfmm::EwaldSolver::solve()
{
    rtfmm::Tree tar_tree;
    tar_tree.build(tar_bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    tar_cells = tar_tree.get_cells();

    rtfmm::Tree src_tree;
    src_tree.build(src_bs, args.x, args.r, args.ncrit, Tree::TreeType::nonuniform);
    src_cells = src_tree.get_cells();

    TIME_BEGIN(real_part);
    RealMap real_map;
    real_part(0, 0, real_map);
    real_p2p(real_map);
    if(args.timing) {TIME_END(real_part);}
    TIME_BEGIN(fourier_part);
    fourier_part();
    if(args.timing) {TIME_END(fourier_part);}
    TIME_BEGIN(self_correction);
    self_correction();
    if(args.timing) {TIME_END(self_correction);}

    return sort_bodies_by_idx(tar_bs);
}

void rtfmm::EwaldSolver::real_part(int this_cell_idx, int that_cell_idx, RealMap& real_map)
{
    rtfmm::Cell3& this_cell = tar_cells[this_cell_idx];
    if (this_cell.crange.number == 0) 
    {
        child(this_cell_idx, that_cell_idx, real_map);
    }
    else
    {
        for(int i = 0; i < this_cell.crange.number; i++)
        {
            real_part(this_cell.crange.offset + i, that_cell_idx, real_map);
        }
    }
}


void rtfmm::EwaldSolver::fourier_part()
{
    DFT(waves,src_bs);
    #ifndef EWALD_USE_HALF
    rtfmm::real coef = 4 * M_PI / std::pow(args.cycle,3);
    #else
    rtfmm::real coef = 8 * M_PI / std::pow(args.cycle,3);
    #endif
    #pragma omp parallel for
    for (int w = 0; w < waves.size(); w++) 
    {
        rtfmm::real k2 = waves[w].K.norm() * scale * scale; //according to the definition of reciprocal vector, ai*bj=2pi*delta_ij. when T is peroid, ai=T, therefore, bj=2pi/T
        waves[w].val *= coef * std::exp(-k2 / (4 * alpha * alpha)) / k2;
    }
    IDFT(waves,tar_bs);
}


void rtfmm::EwaldSolver::self_correction()
{
    for (int i = 0; i < tar_bs.size(); i++) 
    {
        tar_bs[i].p -= M_2_SQRTPI * tar_bs[i].q * alpha;
    }
}


void rtfmm::EwaldSolver::child(int this_cell_idx, int that_cell_idx, RealMap& real_map)
{
    rtfmm::Cell3& this_cell = tar_cells[this_cell_idx];
    rtfmm::Cell3& that_cell = src_cells[that_cell_idx];
    rtfmm::vec3r offset(0,0,0);
    rtfmm::vec3r dx = this_cell.x - that_cell.x;
    for (int d = 0; d < 3; d++)
    {
        if(dx[d] < -args.cycle / 2)
        {
            offset[d]-=args.cycle;
        }
        else if(dx[d] >  args.cycle / 2)
        {
            offset[d]+=args.cycle;
        }
    }
    if ((dx - offset).r() - this_cell.r - that_cell.r < sqrtf(3) * cutoff)
    {
        if(that_cell.crange.number == 0)
        {
            real_map.append(this_cell_idx, that_cell_idx, offset);
        }
        else
        {
            for (int j = 0; j < that_cell.crange.number; j++)
            {
                child(this_cell_idx, that_cell.crange.offset + j, real_map);
            }
        }
    }
}

void rtfmm::EwaldSolver::real_p2p(RealMap& real_map)
{
    int ss = real_map.real_map.size();
    #pragma omp parallel for
    for(int i = 0; i < ss; i++)
    {
        auto real_list = real_map.real_map.begin();
        std::advance(real_list, i);
        int tar = real_list->first;
        std::vector<IndexAndOffset> list = real_list->second;
        for(int j = 0; j < list.size(); j++)
        {
            int src = list[j].index;
            vec3r offset = list[j].offset;
            real_p2p_kernel(tar, src, offset);
        }
    }
}

void rtfmm::EwaldSolver::real_p2p_kernel(int this_cell_idx, int that_cell_idx, rtfmm::vec3r offset)
{
    rtfmm::Cell3& this_cell = tar_cells[this_cell_idx];
    rtfmm::Cell3& that_cell = src_cells[that_cell_idx];
    for (int i = 0; i < this_cell.brange.number; i++)
    {
        rtfmm::Body3& bi = tar_bs[this_cell.brange.offset + i];
        real psum = 0;
        vec3r fsum(0,0,0);
        for (int j = 0; j < that_cell.brange.number; j++)
        {  
            rtfmm::Body3& bj = src_bs[that_cell.brange.offset + j];
            rtfmm::vec3r dx = bi.x - bj.x - offset;
            rtfmm::real r = dx.r();
            //printf("r = %.8f\n",r);
            if (r > 0 && r < cutoff)
            {
                psum += bj.q * std::erfc(alpha * r) / r;
                fsum += -bj.q * dx / std::pow(r,3) * (std::erfc(alpha * r) + (2 * alpha * r * std::pow(M_E, -alpha * alpha * r * r)) / std::sqrt(M_PI));
            }
        }
        //printf("real = %.8f\n", psum);
        bi.p += psum;
        bi.f += fsum;
    }
}


void rtfmm::EwaldSolver::DFT(std::vector<Wave>& ws, std::vector<rtfmm::Body3>& bs)
{
    #pragma omp parallel for
    for (int w = 0; w < ws.size(); w++) 
    {                     
        for (int i = 0; i < bs.size(); i++)  
        {                 
            rtfmm::real ph = (ws[w].K * bs[i].x).sum() * scale;                                
            ws[w].val += bs[i].q * std::pow(rtfmm::real(M_E), complexr(0, -ph));
        }                                                       
    } 
}


void rtfmm::EwaldSolver::IDFT(std::vector<Wave>& ws, std::vector<rtfmm::Body3>& bs)
{
    #pragma omp parallel for
    for (int i = 0; i < bs.size(); i++)
    {
        for (int w = 0; w < ws.size(); w++)
        {
            rtfmm::real ph = (ws[w].K * bs[i].x).sum() * scale;
            bs[i].p += std::real(ws[w].val * std::pow(rtfmm::real(M_E), complexr(0, ph)));
            bs[i].f += std::real(ws[w].val * std::pow(rtfmm::real(M_E), complexr(0, ph)) * complexr(0, 1) * scale) * ws[w].K; // actually this is dp/dx, while real f is -dp/dx*q
        }              
    }   
}
