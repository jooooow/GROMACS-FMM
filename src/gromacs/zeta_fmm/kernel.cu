#include "kernel.h"

#define P2P_KERNEL_TH_NUM 64

__device__ inline zeta_fmm::real reg_w(zeta_fmm::real x)
{
    return 0.25 * (2 + 3 * x - x * x * x);
}

__device__ inline zeta_fmm::real reg_dw(zeta_fmm::real x)
{
    return 0.25 * (3 - 3 * x * x);
}

__device__ inline zeta_fmm::real get_w_single(zeta_fmm::real dx, zeta_fmm::real R, zeta_fmm::real rega)
{
    zeta_fmm::real r = std::abs(dx) - R + rega;
    zeta_fmm::real t;
    if(r <= 0) t = 1;
    else if(r > 2 * rega) t = -1;
    else t = 1 - r / rega;
    return reg_w(t);
}

__device__ inline zeta_fmm::real get_dw_single(zeta_fmm::real dx, zeta_fmm::real R, zeta_fmm::real rega)
{
    zeta_fmm::real r = std::abs(dx) - R + rega;
    zeta_fmm::real dw;
    if(r > 0 && r <= 2 * rega)
    {
        dw = reg_dw(1 - r / rega) * (-1 / rega) * (dx > 0 ? 1 : -1);
    }
    else
    {
        dw = 0;
    }

    return dw;
}

__device__ inline zeta_fmm::real get_w(zeta_fmm::vec3r dx, zeta_fmm::real R, zeta_fmm::real rega)
{
    zeta_fmm::real w = 1;
    for(int d = 0; d < 3; d++)
    {
        w *= get_w_single(dx[d], R, rega);
    }
    return w;
}

__device__ inline zeta_fmm::vec3r get_dw(zeta_fmm::vec3r dx, zeta_fmm::real R, zeta_fmm::real rega)
{
    zeta_fmm::vec3r dw;
    dw[0] = get_dw_single(dx[0], R, rega) * get_w_single(dx[1], R, rega) * get_w_single(dx[2], R, rega);
    dw[1] = get_w_single(dx[0], R, rega) * get_dw_single(dx[1], R, rega) * get_w_single(dx[2], R, rega);
    dw[2] = get_w_single(dx[0], R, rega) * get_w_single(dx[1], R, rega) * get_dw_single(dx[2], R, rega);
    return dw;
}

__device__ inline zeta_fmm::vec3r add(zeta_fmm::vec3r a, zeta_fmm::Offset3rPadding b)
{
    zeta_fmm::vec3r c;
    c[0] = a[0] + b.x;
    c[1] = a[1] + b.y;
    c[2] = a[2] + b.z;
    return c;
}

__device__ int oddOrEven(int n)
{
    return (((n) & 1) == 1) ? -1 : 1;
}

__device__ int ipow2n(int n) 
{
    return (n >= 0) ? 1 : oddOrEven(n);
}

//d/d(r) -> d/d(x)
__device__ zeta_fmm::vec3r derivate_sph2cart(zeta_fmm::Sph dev_sph_cart, zeta_fmm::vec3r dev_p_sph) 
{
    zeta_fmm::real sin_theta = std::sin(dev_sph_cart.theta);
    zeta_fmm::real inv_sin_theta = sin_theta == 0 ? 0 : 1 / sin_theta;
    zeta_fmm::real inv_rho = dev_sph_cart.rho == 0 ? 0 : 1 / dev_sph_cart.rho;

    zeta_fmm::vec3r cart;
    cart[0] = sin_theta * std::cos(dev_sph_cart.phi) * dev_p_sph[0]
            + std::cos(dev_sph_cart.theta) * std::cos(dev_sph_cart.phi) * inv_rho * dev_p_sph[1]
            - std::sin(dev_sph_cart.phi) * inv_rho * inv_sin_theta * dev_p_sph[2];
    cart[1] = sin_theta * std::sin(dev_sph_cart.phi) * dev_p_sph[0]
            + std::cos(dev_sph_cart.theta) * std::sin(dev_sph_cart.phi) * inv_rho * dev_p_sph[1]
            + std::cos(dev_sph_cart.phi) * inv_rho * inv_sin_theta * dev_p_sph[2];
    cart[2] = std::cos(dev_sph_cart.theta) * dev_p_sph[0]
            - sin_theta * inv_rho * dev_p_sph[1];
    return cart;
}

__device__ void calc_Ynm(zeta_fmm::complex* Ynm, int P, zeta_fmm::real rho, zeta_fmm::real theta, zeta_fmm::real phi)
{
    zeta_fmm::real x = std::cos(theta);
    zeta_fmm::real y = std::sin(theta);
    zeta_fmm::real fact = 1;
    zeta_fmm::real pn = 1;
    zeta_fmm::real rhom = 1;
    zeta_fmm::complex ei;
    ei.rel = 0; ei.img = phi;
    ei = ei.exp();
    zeta_fmm::complex eim;
    eim.rel = 1; eim.img = 0;
    for (int m=0; m<=P; m++)
    {
        zeta_fmm::real p = pn;
        int npn = m * m + 2 * m;
        int nmn = m * m;
        Ynm[npn] = rhom * p * eim;
        Ynm[nmn] = Ynm[npn].conj();
        zeta_fmm::real p1 = p;
        p = x * (2 * m + 1) * p1;
        rhom *= rho;
        zeta_fmm::real rhon = rhom;
        for (int n=m+1; n<=P; n++) 
        {
            int npm = n * n + n + m;
            int nmm = n * n + n - m;
            rhon /= -(n + m);
            Ynm[npm] = rhon * p * eim;
            Ynm[nmm] = Ynm[npm].conj();
            zeta_fmm::real p2 = p1;
            p1 = p;
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
            rhon *= rho;
        }
        rhom /= -(2 * m + 2) * (2 * m + 1);
        pn = -pn * fact * y;
        fact += 2;
        eim *= ei;
    } 
}

__device__ void calc_YnmD(zeta_fmm::complex* YnmD, int P, zeta_fmm::real rho, zeta_fmm::real theta, zeta_fmm::real phi)
{
    zeta_fmm::real x = std::cos(theta);
    zeta_fmm::real y = std::sin(theta);
    zeta_fmm::real invY = y == 0 ? 0 : 1 / y;
    zeta_fmm::real fact = 1;
    zeta_fmm::real pn = 1;
    zeta_fmm::real rhom = 1;
    zeta_fmm::complex ei;
    ei.rel = 0; ei.img = phi;
    ei = ei.exp();
    zeta_fmm::complex eim;
    eim.rel = 1; eim.img = 0;
    for (int m=0; m<=P; m++)
    {
        zeta_fmm::real p = pn;
        int npn = m * m + 2 * m;
        zeta_fmm::real p1 = p;
        p = x * (2 * m + 1) * p1;
        YnmD[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim;
        rhom *= rho;
        zeta_fmm::real rhon = rhom;
        for (int n=m+1; n<=P; n++) 
        {
            int npm = n * n + n + m;
            rhon /= -(n + m);
            zeta_fmm::real p2 = p1;
            p1 = p;
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
            YnmD[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;
            rhon *= rho;
        }
        rhom /= -(2 * m + 2) * (2 * m + 1);
        pn = -pn * fact * y;
        fact += 2;
        eim *= ei;
    } 
}

__device__ void calc_Ynm2(zeta_fmm::complex* Ynm, int P, zeta_fmm::real rho, zeta_fmm::real theta, zeta_fmm::real phi)
{
    zeta_fmm::real x = std::cos(theta);                                
    zeta_fmm::real y = std::sin(theta);                                
    zeta_fmm::real fact = 1;                                           
    zeta_fmm::real pn = 1;                                             
    zeta_fmm::real invR = -1.0 / rho;                                  
    zeta_fmm::real rhom = -invR;                                       
    zeta_fmm::complex ei;
    ei.rel = 0; ei.img = phi;
    ei = ei.exp();
    zeta_fmm::complex eim;
    eim.rel = 1; eim.img = 0;                                
    for (int m=0; m<=P; m++) 
    {                                  
        zeta_fmm::real p = pn;                                           
        int npn = m * m + 2 * m;                                 
        int nmn = m * m;                                         
        Ynm[npn] = rhom * p * eim;                               
        Ynm[nmn] = Ynm[npn].conj();                          
        zeta_fmm::real p1 = p;                                           
        p = x * (2 * m + 1) * p1;                                
        rhom *= invR;                                            
        zeta_fmm::real rhon = rhom;                                      
        for (int n=m+1; n<=P; n++) {                              
            int npm = n * n + n + m;                               
            int nmm = n * n + n - m;                               
            Ynm[npm] = rhon * p * eim;                             
            Ynm[nmm] = Ynm[npm].conj();                        
            zeta_fmm::real p2 = p1;                                        
            p1 = p;                                                
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
            rhon *= invR * (n - m + 1);                            
        }                                                        
        pn = -pn * fact * y;                                     
        fact += 2;                                               
        eim *= ei;                                               
    }                                                          
}

__global__ void p2m(zeta_fmm::Body3* bodies, zeta_fmm::Cell3* cells, int* leaf_cells, int P, zeta_fmm::complex* Ms)
{
    int cell_idx = leaf_cells[blockIdx.x];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* M = Ms + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];

    for(int i = 0; i < cell->body_info.number; i++)
    {
        zeta_fmm::Body3* body = bodies + cell->body_info.offset + i;
        zeta_fmm::vec3r dx = body->loc - cell->center;
        zeta_fmm::Sph sph = dx.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, -sph.phi);
        for (int n = 0; n <= P; n++) 
        {
            for (int m = 0; m <= n; m++) 
            {
                M[n * (n + 1) / 2 + m] += body->q * Ynm[n * n + n + m];
            }
        }
    }
}

__global__ void p2m_dummy(zeta_fmm::Body3* bodies, zeta_fmm::Cell3* cells, int* leaf_cells, int P, zeta_fmm::complex* Ms)
{
    int cell_idx = leaf_cells[blockIdx.x];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* M = Ms + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];

    for(int i = 0; i < cell->body_info.number; i++)
    {
        zeta_fmm::Body3* body = bodies + cell->body_info.offset + i;
        zeta_fmm::vec3r dx = body->loc - cell->center;
        zeta_fmm::Sph sph = dx.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, -sph.phi);
        int have_nan = 0;
        for (int n = 0; n <= P; n++) 
        {
            for (int m = 0; m <= n; m++) 
            {
                M[n * (n + 1) / 2 + m] = Ynm[n * n + n + m];
                if(isnan(Ynm[n * n + n + m].rel) || isnan(Ynm[n * n + n + m].rel))
                {
                    have_nan = 1;
                }
            }
        }
        if(have_nan)
        {
            M[0].rel = dx[0];
            M[0].img = dx[1];
            M[1].rel = dx[2];
            M[2].rel = sph.rho;
            M[2].img = sph.theta;
            M[3].rel = sph.phi;
            return;
        }
    }
}

__global__ void p2m_reg(zeta_fmm::Body3* bodies, zeta_fmm::Cell3* cells, int* leaf_cells, int P, zeta_fmm::complex* Ms, zeta_fmm::real rega)
{
    int cell_idx = leaf_cells[blockIdx.x];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* M = Ms + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];

    for(int i = 0; i < cell->body_info.number; i++)
    {
        zeta_fmm::Body3* body = bodies + cell->body_info.offset + i;
        zeta_fmm::vec3r dx = body->loc - cell->center;
        zeta_fmm::Sph sph = dx.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, -sph.phi);
        zeta_fmm::real w = get_w(dx, cell->r, rega);
        for (int n = 0; n <= P; n++) 
        {
            for (int m = 0; m <= n; m++) 
            {
                M[n * (n + 1) / 2 + m] += w * body->q * Ynm[n * n + n + m];
            }
        }
    }
}

__global__ void m2m(zeta_fmm::Cell3* cells, int* branch_cells, int offset, int P, zeta_fmm::complex* Ms)
{
    int cell_idx = branch_cells[blockIdx.x + offset];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* M = Ms + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];

    for(int i = 0; i < cell->child_info.number; i++)
    {
        int child_cell_idx = cell->child_info.offset + i;
        zeta_fmm::Cell3* child_cell = cells + child_cell_idx;
        zeta_fmm::complex* child_M = Ms + child_cell_idx * pole_size_eachcell;
        
        zeta_fmm::vec3r dx = cell->center - child_cell->center;
        zeta_fmm::Sph sph = dx.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, sph.phi);
        for(int j = 0; j <= P; j++) 
        {
            for (int k = 0; k <= j; k++) 
            {
                int jks = j * (j + 1) / 2 + k;
                zeta_fmm::complex temp_M;
                temp_M.rel = 0; temp_M.img = 0;
                for (int n=0; n<=j; n++) 
                {
                    for (int m=max(-n,-j+k+n); m<=min(k-1,n); m++) 
                    {
                        int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                        int nm    = n * n + n - m;
                        temp_M += child_M[jnkms] * Ynm[nm] * zeta_fmm::real(ipow2n(m) * oddOrEven(n));
                    }
                    for (int m=k; m<=min(n,j+k-n); m++) 
                    {
                        int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                        int nm    = n * n + n - m;
                        temp_M += child_M[jnkms].conj() * Ynm[nm] * zeta_fmm::real(oddOrEven(k+n+m));
                    }
                }
                M[jks] += temp_M;
            }
        }
    }
}

__global__ void m2m_image(zeta_fmm::Cell3* cells, int* image_cells, int offset, int P, zeta_fmm::complex* Ms, zeta_fmm::real cycle)
{
    int cell_idx = image_cells[blockIdx.x + offset];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* M = Ms + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];

    int child_cell_idx = cell->child_info.offset;
    zeta_fmm::Cell3* child_cell = cells + child_cell_idx;
    zeta_fmm::complex* child_M = Ms + child_cell_idx * pole_size_eachcell;
    
    for(int pz = -1; pz <= 1; pz++)
    {
        for(int py = -1; py <= 1; py++)
        {
            for(int px = -1; px <= 1; px++)
            {
                zeta_fmm::vec3r dx = cell->center - child_cell->center - zeta_fmm::vec3r(px,py,pz) * cycle;
                zeta_fmm::Sph sph = dx.toSph2();
                calc_Ynm(Ynm, P, sph.rho, sph.theta, sph.phi);
                for(int j = 0; j <= P; j++) 
                {
                    for (int k = 0; k <= j; k++) 
                    {
                        int jks = j * (j + 1) / 2 + k;
                        zeta_fmm::complex temp_M;
                        temp_M.rel = 0; temp_M.img = 0;
                        for (int n=0; n<=j; n++) 
                        {
                            for (int m=max(-n,-j+k+n); m<=min(k-1,n); m++) 
                            {
                                int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                                int nm    = n * n + n - m;
                                temp_M += child_M[jnkms] * Ynm[nm] * zeta_fmm::real(ipow2n(m) * oddOrEven(n));
                            }
                            for (int m=k; m<=min(n,j+k-n); m++) 
                            {
                                int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                                int nm    = n * n + n - m;
                                temp_M += child_M[jnkms].conj() * Ynm[nm] * zeta_fmm::real(oddOrEven(k+n+m));
                            }
                        }
                        M[jks] += temp_M;
                    }
                }
            }
        }
    }
}

__global__ void l2l(zeta_fmm::Cell3* cells, int* branch_cells, int offset, int P, zeta_fmm::complex* Ls)
{
    int cell_idx = branch_cells[blockIdx.x + offset];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* L = Ls + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];

    for(int i = 0; i < cell->child_info.number; i++)
    {
        int child_cell_idx = cell->child_info.offset + i;
        zeta_fmm::Cell3* child_cell = cells + child_cell_idx;
        zeta_fmm::complex* child_L = Ls + child_cell_idx * pole_size_eachcell;        
        zeta_fmm::vec3r dX = child_cell->center - cell->center;
        zeta_fmm::Sph sph = dX.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, sph.phi);
        for (int j=0; j<=P; j++) 
        {
            for (int k=0; k<=j; k++) 
            {
                int jks = j * (j + 1) / 2 + k;
                zeta_fmm::complex temp_L;
                temp_L.rel = 0; temp_L.img = 0;
                for (int n=j; n<=P; n++) 
                {
                    for (int m=j+k-n; m<0; m++) 
                    {
                        int jnkm = (n - j) * (n - j) + n - j + m - k;
                        int nms  = n * (n + 1) / 2 - m;
                        temp_L += L[nms].conj() * Ynm[jnkm] * zeta_fmm::real(oddOrEven(k));
                    }
                    for (int m=0; m<=n; m++) 
                    {
                        if (n-j >= abs(m-k)) 
                        {
                            int jnkm = (n - j) * (n - j) + n - j + m - k;
                            int nms  = n * (n + 1) / 2 + m;
                            temp_L += L[nms] * Ynm[jnkm] * zeta_fmm::real(oddOrEven((m-k)*(m<k)));
                        }
                    }
                }
                child_L[jks] += temp_L;
            }
        }
    }
}

__global__ void p2p(
    zeta_fmm::Body3* bodies, 
    zeta_fmm::Cell3* cells, 
    int* p2p_matrix, zeta_fmm::Offset3rPadding* p2p_offset_matrix, int p2p_matrix_col
)
{
    int tar_cell_idx = p2p_matrix[blockIdx.x * p2p_matrix_col];
    int src_cell_num = p2p_matrix[blockIdx.x * p2p_matrix_col + 1];
    zeta_fmm::Cell3* tar_cell = cells + tar_cell_idx;

    for(int offset = 0; offset < src_cell_num; offset++)
    {
        int src_cell_idx = p2p_matrix[blockIdx.x * p2p_matrix_col + 2 + offset];
        zeta_fmm::Offset3rPadding src_cell_offset = p2p_offset_matrix[blockIdx.x * (p2p_matrix_col - 2) + offset];
        zeta_fmm::Cell3* src_cell = cells + src_cell_idx;

        for(int tar_offset = 0; tar_offset < tar_cell->body_info.number; tar_offset += P2P_KERNEL_TH_NUM)
        {
            if(tar_offset + threadIdx.x < tar_cell->body_info.number)
            {
                int tar_body_idx = tar_cell->body_info.offset + tar_offset + threadIdx.x;
                zeta_fmm::Body3* b_tar = bodies + tar_body_idx;
                zeta_fmm::real p = 0.0;
                zeta_fmm::vec3r f(0,0,0);
                for(int j = src_cell->body_info.offset; j < src_cell->body_info.offset + src_cell->body_info.number; j++)
                {
                    zeta_fmm::Body3* b_src = bodies + j;
                    zeta_fmm::vec3r dx = b_tar->loc - add(b_src->loc, src_cell_offset);
                    zeta_fmm::real r = dx.r();
                    if(r > 0)
                    {
                        p += b_src->q / r;
                        f += b_src->q / r / r * (-dx) / r;
                    }
                }
                b_tar->p += p;
                b_tar->f += f;
            }
        }
    }
}

__global__ void p2p_reg(
    zeta_fmm::Body3* bodies, 
    zeta_fmm::Cell3* cells, 
    int* p2p_matrix, zeta_fmm::Offset3rPadding* p2p_offset_matrix, int p2p_matrix_col,
    zeta_fmm::real rega
)
{
    int tar_cell_idx = p2p_matrix[blockIdx.x * p2p_matrix_col];
    int src_cell_num = p2p_matrix[blockIdx.x * p2p_matrix_col + 1];
    zeta_fmm::Cell3* tar_cell = cells + tar_cell_idx;

    for(int offset = 0; offset < src_cell_num; offset++)
    {
        int src_cell_idx = p2p_matrix[blockIdx.x * p2p_matrix_col + 2 + offset];
        zeta_fmm::Offset3rPadding src_cell_offset = p2p_offset_matrix[blockIdx.x * (p2p_matrix_col - 2) + offset];
        zeta_fmm::Cell3* src_cell = cells + src_cell_idx;

        for(int tar_offset = 0; tar_offset < tar_cell->body_info.number; tar_offset += P2P_KERNEL_TH_NUM)
        {
            if(tar_offset + threadIdx.x < tar_cell->body_info.number)
            {
                int tar_body_idx = tar_cell->body_info.offset + tar_offset + threadIdx.x;
                zeta_fmm::Body3* b_tar = bodies + tar_body_idx;
                zeta_fmm::real wi = get_w(b_tar->loc - tar_cell->center, tar_cell->r, rega);
                #if FMM_REG_USE_DW
                zeta_fmm::vec3r dwi = get_dw(b_tar->loc - tar_cell->center, tar_cell->r, rega);
                #endif
                zeta_fmm::real p = 0.0;
                zeta_fmm::vec3r f(0,0,0);
                for(int j = src_cell->body_info.offset; j < src_cell->body_info.offset + src_cell->body_info.number; j++)
                {
                    zeta_fmm::Body3* b_src = bodies + j;
                    zeta_fmm::real wj = get_w(b_src->loc - src_cell->center, src_cell->r, rega);
                    zeta_fmm::vec3r dx = b_tar->loc - add(b_src->loc, src_cell_offset);
                    zeta_fmm::real r = dx.r();
                    if(r > 0)
                    {
                        p += wj * b_src->q / r;
                        f += wj * b_src->q / r / r * (-dx) / r;
                    }
                }
                b_tar->p += wi * p;
                #if FMM_REG_USE_DW
                b_tar->f += wi * f + dwi * p;
                #else
                b_tar->f += wi * f;
                #endif
            }
        }
    }
}

__global__ void m2l(
    zeta_fmm::Cell3* cells, 
    int* m2l_matrix, zeta_fmm::Offset3rPadding* m2l_offset_matrix, int m2l_matrix_col, 
    int P, zeta_fmm::complex* Ms, zeta_fmm::complex* Ls
)
{
    int tar_cell_idx = m2l_matrix[blockIdx.x * m2l_matrix_col];
    int src_cell_num = m2l_matrix[blockIdx.x * m2l_matrix_col + 1];
    int pole_size_eachcell = get_pole_size_eachcell(P);

    __shared__ zeta_fmm::complex Ynm[512];
    
    zeta_fmm::Cell3* tar_cell = cells + tar_cell_idx;
    zeta_fmm::complex* tar_L = Ls + tar_cell_idx * pole_size_eachcell;
    
    for(int i = 0; i < src_cell_num; i++)
    {
        int src_cell_idx = m2l_matrix[blockIdx.x * m2l_matrix_col + 2 + i];
        zeta_fmm::Offset3rPadding src_cell_offset = m2l_offset_matrix[blockIdx.x * (m2l_matrix_col - 2) + i];
        zeta_fmm::Cell3* src_cell = cells + src_cell_idx;
        zeta_fmm::complex* src_M = Ms + src_cell_idx * pole_size_eachcell;
        
        zeta_fmm::vec3r dX = tar_cell->center - add(src_cell->center, src_cell_offset);
        zeta_fmm::Sph sph = dX.toSph2();
        calc_Ynm2(Ynm, P, sph.rho, sph.theta, sph.phi);
        for (int j=0; j<=P; j++) 
        {
            zeta_fmm::real Cnm = oddOrEven(j);
            for (int k=0; k<=j; k++) 
            {
                int jks = j * (j + 1) / 2 + k;
                zeta_fmm::complex temp_L;
                temp_L.rel = 0; temp_L.img = 0;
                for (int n=0; n<=P; n++) 
                {
                    if(j + n <= P)   // WARNING : this line is forgot in exafmm/minial(if lack, jnkm for Ynm will be out of bounds)
                    {
                        for (int m=-n; m<0; m++) 
                        {
                            int nms  = n * (n + 1) / 2 - m;
                            int jnkm = (j + n) * (j + n) + j + n + m - k;
                            temp_L += src_M[nms].conj() * Cnm * Ynm[jnkm];
                        }
                        for (int m=0; m<=n; m++) 
                        {
                            int nms  = n * (n + 1) / 2 + m;
                            int jnkm = (j + n) * (j + n) + j + n + m - k;
                            zeta_fmm::real Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
                            temp_L += src_M[nms] * Cnm2 * Ynm[jnkm];
                        }
                    }
                }
                tar_L[jks] += temp_L;
            }
        }
    }
}

__global__ void l2p(
    zeta_fmm::Body3* bodies, 
    zeta_fmm::Cell3* cells, 
    int* leaf_cells, 
    int P, 
    zeta_fmm::complex* Ls
)
{
    int cell_idx = leaf_cells[blockIdx.x];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* L = Ls + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];
    __shared__ zeta_fmm::complex YnmD[512];

    for(int i = 0; i < cell->body_info.number; i++)
    {
        zeta_fmm::Body3* body = bodies + cell->body_info.offset + i;
        zeta_fmm::vec3r dx = body->loc - cell->center;
        zeta_fmm::Sph sph = dx.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, sph.phi);
        calc_YnmD(YnmD, P, sph.rho, sph.theta, sph.phi);
        zeta_fmm::vec3r temp = {0,0,0};
        zeta_fmm::complex I;
        I.rel = 0; I.img = 1;
        zeta_fmm::real p = 0;
        for (int n=0; n<=P; n++)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;
            p += (L[nms] * Ynm[nm]).rel;
            temp[0] += (L[nms] * Ynm[nm]).rel / sph.rho * n;
            temp[1] += (L[nms] * YnmD[nm]).rel;
            for (int m=1; m<=n; m++) 
            {
                nm  = n * n + n + m;
                nms = n * (n + 1) / 2 + m;
                p += 2 * (L[nms] * Ynm[nm]).rel;
                temp[0] += 2 * (L[nms] * Ynm[nm]).rel / sph.rho * n;
                temp[1] += 2 * (L[nms] * YnmD[nm]).rel;
                temp[2] += 2 * (L[nms] * Ynm[nm] * I).rel * m;
            }
        }
        body->p += p;
        zeta_fmm::vec3r cart = derivate_sph2cart(sph, temp);
        body->f += cart;
    }
}

__global__ void l2p_dummy(
    zeta_fmm::Body3* bodies, 
    zeta_fmm::Cell3* cells, 
    int* leaf_cells, 
    int P, 
    zeta_fmm::complex* Ls
)
{
    int cell_idx = leaf_cells[blockIdx.x];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3* cell = cells + cell_idx;
    zeta_fmm::complex* L = Ls + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];
    __shared__ zeta_fmm::complex YnmD[512];

    for(int i = 0; i < cell->body_info.number; i++)
    {
        zeta_fmm::Body3* body = bodies + cell->body_info.offset + i;
        zeta_fmm::vec3r dx = body->loc - cell->center;
        zeta_fmm::Sph sph = dx.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, sph.phi);
        calc_YnmD(YnmD, P, sph.rho, sph.theta, sph.phi);
        zeta_fmm::vec3r temp = {0,0,0};
        zeta_fmm::complex I;
        I.rel = 0; I.img = 1;
        zeta_fmm::real p = 0;
        for (int n=0; n<=P; n++)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;
            p += (L[nms] * Ynm[nm]).rel;
            temp[0] += (L[nms] * Ynm[nm]).rel / sph.rho * n;
            temp[1] += (L[nms] * YnmD[nm]).rel;
            for (int m=1; m<=n; m++) 
            {
                nm  = n * n + n + m;
                nms = n * (n + 1) / 2 + m;
                p += 2 * (L[nms] * Ynm[nm]).rel;
                temp[0] += 2 * (L[nms] * Ynm[nm]).rel / sph.rho * n;
                temp[1] += 2 * (L[nms] * YnmD[nm]).rel;
                temp[2] += 2 * (L[nms] * Ynm[nm] * I).rel * m;
            }
        }
        body->p += p;
        zeta_fmm::vec3r cart = derivate_sph2cart(sph, temp);
        body->f += cart;
    }
}

__global__ void l2p_reg(
    zeta_fmm::Body3* bodies,
    zeta_fmm::Cell3* cells, 
    int* leaf_cells, 
    int P, 
    zeta_fmm::complex* Ls, 
    zeta_fmm::real rega
)
{
    int cell_idx = leaf_cells[blockIdx.x];
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::Cell3 cell = cells[cell_idx];
    zeta_fmm::complex* L = Ls + cell_idx * pole_size_eachcell;

    __shared__ zeta_fmm::complex Ynm[512];
    __shared__ zeta_fmm::complex YnmD[512];

    for(int i = 0; i < cell.body_info.number; i++)
    {
        int bidx = cell.body_info.offset + i;
        zeta_fmm::Body3* body = bodies + bidx;
        zeta_fmm::vec3r dx = body->loc - cell.center;
        zeta_fmm::real w = get_w(dx, cell.r, rega);
        #ifdef FMM_REG_USE_DW
        zeta_fmm::vec3r dw = get_dw(dx, cell.r, rega);
        #endif
        zeta_fmm::Sph sph = dx.toSph2();
        calc_Ynm(Ynm, P, sph.rho, sph.theta, sph.phi);
        calc_YnmD(YnmD, P, sph.rho, sph.theta, sph.phi);
        zeta_fmm::vec3r temp = {0,0,0};
        zeta_fmm::complex I;
        I.rel = 0; I.img = 1;
        zeta_fmm::real p = 0;
        for (int n=0; n<=P; n++)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;
            p += (L[nms] * Ynm[nm]).rel;
            temp[0] += (L[nms] * Ynm[nm]).rel / sph.rho * n;
            temp[1] += (L[nms] * YnmD[nm]).rel;
            for (int m=1; m<=n; m++) 
            {
                nm  = n * n + n + m;
                nms = n * (n + 1) / 2 + m;
                p += 2 * (L[nms] * Ynm[nm]).rel;
                temp[0] += 2 * (L[nms] * Ynm[nm]).rel / sph.rho * n;
                temp[1] += 2 * (L[nms] * YnmD[nm]).rel;
                temp[2] += 2 * (L[nms] * Ynm[nm] * I).rel * m;
            }
        }
        zeta_fmm::vec3r cart = derivate_sph2cart(sph, temp);
        #if FMM_REG_USE_DW
        body->f += w * cart + p * dw;
        #else
        body->f += w * cart;
        #endif
        body->p += w * p;
    }
}

__global__ void arrange_bodies(
    zeta_fmm::Body3* bodies, 
    zeta_fmm::Cell3* cells, 
    int* leaf_cells, 
    int merged_num_body, zeta_fmm::real* premerged_body_p, zeta_fmm::vec3r* premerged_body_f
)
{
    int cell_idx = leaf_cells[blockIdx.x];
    zeta_fmm::Cell3 cell = cells[cell_idx];

    for(int i = 0; i < cell.body_info.number; i++)
    {
        int bidx = cell.body_info.offset + i;
        zeta_fmm::Body3 body = bodies[bidx];
        int partion_idx = ((body.loc[0] > cell.center[0]) << 2) 
                        + ((body.loc[1] > cell.center[1]) << 1) 
                        + ((body.loc[2] > cell.center[2]) << 0);
        premerged_body_p[merged_num_body * partion_idx + body.idx] = body.p;
        premerged_body_f[merged_num_body * partion_idx + body.idx] = body.f;
    }
}

__global__ void merge_bodies(
    zeta_fmm::Body3* merged_bodies, 
    int merged_num_body, zeta_fmm::real* premerged_body_p, zeta_fmm::vec3r* premerged_body_f
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < merged_num_body)
    {
        zeta_fmm::real p = 0;
        zeta_fmm::vec3r f(0,0,0);
        for(int i = 0; i < 8; i++)
        {
            p += premerged_body_p[merged_num_body * i + idx];
            f += premerged_body_f[merged_num_body * i + idx];
        }
        merged_bodies[idx].p = p;
        merged_bodies[idx].f = f;
    }
}

zeta_fmm::GPUKernel::GPUKernel()
{
    //std::cout<<"init GPUKernel"<<std::endl;
}

void zeta_fmm::GPUKernel::solve(
    int P, zeta_fmm::real rega,
    Body3* bodies, int real_num_body, 
    Body3* merged_bodies, int merged_num_body,
    Cell3* cells, int num_cell, 
    zeta_fmm::complex* Ms, zeta_fmm::complex* Ls,
    int* leaf_cells, int leaf_cell_num,
    int* branch_cells, std::vector<zeta_fmm::OffsetAndNumber> level_infos,
    int* image_cells, int image_cell_num, zeta_fmm::real cycle,
    int* p2p_matrix, Offset3rPadding* p2p_offset_matrix, int p2p_matrix_rows, int p2p_matrix_col,
    int* m2l_matrix, Offset3rPadding* m2l_offset_matrix, int m2l_matrix_rows, int m2l_matrix_col,
    int dummy
)
{
    // P2M
    if(!dummy)
    {
        if(rega == 0)
            p2m<<<leaf_cell_num, 1>>>(bodies, cells, leaf_cells, P, Ms);
        else
            p2m_reg<<<leaf_cell_num, 1>>>(bodies, cells, leaf_cells, P, Ms, rega);
    }
    else
    {
        p2m_dummy<<<leaf_cell_num, 1>>>(bodies, cells, leaf_cells, P, Ms);
        int pole_size_eachcell = get_pole_size_eachcell(P);
        zeta_fmm::complex* c_Ls = new zeta_fmm::complex[num_cell * pole_size_eachcell];
        cudaMemcpy(c_Ls, Ms, num_cell * pole_size_eachcell * sizeof(zeta_fmm::complex), cudaMemcpyDeviceToHost);
        for(int i = 0; i < num_cell; i++)
        {
            printf("M[%d]= ", i);
            for(int p = 0; p < pole_size_eachcell; p++)
            {
                zeta_fmm::complex L = c_Ls[i * pole_size_eachcell + p];
                printf("(%.8f,%.8f)  ", L.rel, L.img);
            }
            printf("\n");
        }
        delete[] c_Ls;
    }

    // M2M
    for(int i = 0; i < level_infos.size(); i++)
    {
        zeta_fmm::OffsetAndNumber offset_number = level_infos[i];
        m2m<<<offset_number.number, 1>>>(cells, branch_cells, offset_number.offset, P, Ms);
    }

    // M2M for image cells
    zeta_fmm::real c = cycle;
    for(int i = 0; i < image_cell_num - 1; i++) // the last one is not necessary
    {
        m2m_image<<<1,1>>>(cells, image_cells, i, P, Ms, c);
        c *= 3;
    }

    // M2L
    m2l<<<m2l_matrix_rows,1>>>(cells, m2l_matrix, m2l_offset_matrix, m2l_matrix_col, P, Ms, Ls);

    // P2P  
    if(rega == 0)
    {
        p2p<<<p2p_matrix_rows,P2P_KERNEL_TH_NUM>>>(bodies, cells, p2p_matrix, p2p_offset_matrix, p2p_matrix_col);
    }
    else
    {
        p2p_reg<<<p2p_matrix_rows,P2P_KERNEL_TH_NUM>>>(bodies, cells, p2p_matrix, p2p_offset_matrix, p2p_matrix_col, rega);
    }

    // L2L
    for(int i = level_infos.size() - 1; i >= 0; i--)
    {
        zeta_fmm::OffsetAndNumber offset_number = level_infos[i];
        l2l<<<offset_number.number, 1>>>(cells, branch_cells, offset_number.offset, P, Ls);
    }

    // L2P
    
    if(rega == 0)
        l2p<<<leaf_cell_num, 1>>>(bodies, cells, leaf_cells, P, Ls);
    else
        l2p_reg<<<leaf_cell_num, 1>>>(bodies, cells, leaf_cells, P, Ls, rega);

    // arrange and merge reg bodies
    if(rega > 0)
    {
        zeta_fmm::real* g_premerged_body_p;
        cudaMalloc(&g_premerged_body_p, merged_num_body * 8 * sizeof(zeta_fmm::real));
        cudaMemset(g_premerged_body_p, 0.0, merged_num_body * 8 * sizeof(zeta_fmm::real));
        zeta_fmm::vec3r* g_premerged_body_f;
        cudaMalloc(&g_premerged_body_f, merged_num_body * 8 * sizeof(zeta_fmm::vec3r));
        cudaMemset(g_premerged_body_f, 0.0, merged_num_body * 8 * sizeof(zeta_fmm::vec3r));

        arrange_bodies<<<leaf_cell_num, 1>>>(bodies, cells, leaf_cells, merged_num_body, g_premerged_body_p, g_premerged_body_f);

        int merge_block_num = (merged_num_body + 1023) / 1024;
        merge_bodies<<<merge_block_num, 1024>>>(merged_bodies, merged_num_body, g_premerged_body_p, g_premerged_body_f);

        cudaFree(g_premerged_body_p);
        cudaFree(g_premerged_body_f);
    }
}