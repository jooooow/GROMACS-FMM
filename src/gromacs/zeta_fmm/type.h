#pragma once
#include <iostream>
#include <math.h>
#include <complex>
#include <chrono>
#include <cuda.h>
#include <cuda_runtime.h>

#define TIME_VERBOSE


#ifdef TIME_VERBOSE
#define TIME_BEGIN(a) auto time_begin_##a = std::chrono::high_resolution_clock::now()

#define TIME_END(a)   auto time_end_##a = std::chrono::high_resolution_clock::now();\
					  auto elapse_##a = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_##a - time_begin_##a);\
                      printf("[%s time measured : %.5f seconds.]\n", #a, elapse_##a.count() * 1e-9)
#else
#define TIME_BEGIN(a)
#define TIME_END(a)
#endif

#define CHECK_CUDA_ERR(a) cudaDeviceSynchronize();\
                          cudaError_t cuda_err_##a = cudaGetLastError();\
						  if (cuda_err_##a != cudaSuccess){\
                              printf("CUDA Error at %s : %s\n", #a, cudaGetErrorString(cuda_err_##a));\
							  exit(-1);\
						  }

#define get_pole_size_eachcell(P) ((P + 1) * (P + 2) / 2)


namespace zeta_fmm
{

using real = double;
using complexr = std::complex<real>;

struct Sph
{
    real rho;
    real theta;
    real phi;
};

class complex
{
public:
    __host__ __device__
    complex exp()
    {
        complex c;
        c.rel = std::exp(this->rel) * std::cos(this->img);
        c.img = std::exp(this->rel) * std::sin(this->img);
        return c;
    }
    __host__ __device__
    complex conj()
    {
        complex c;
        c.rel = this->rel;
        c.img = -this->img;
        return c;
    }
    __host__ __device__
    complex operator*(const double v) const 
    {
        return complex(*this) *= v;
    }
    __host__ __device__
    const complex &operator+=(const complex v) 
    {
        this->rel += v.rel;
        this->img += v.img;
        return *this;
    }
    __host__ __device__
    const complex &operator*=(const double v) 
    {
        this->rel *= v;
        this->img *= v;
        return *this;
    }
    __host__ __device__
    const complex &operator*=(const complex v) 
    {
        real temp_rel = v.rel * this->rel - v.img * this->img;
        real temp_img = v.rel * this->img + v.img * this->rel;
        this->rel = temp_rel;
        this->img = temp_img;
        return *this;
    }
    __host__ __device__
    friend complex operator*(const double v, const complex& _c)  // value * this
    {
        complex res;
        res.rel = v * _c.rel;
        res.img = v * _c.img;
        return res;
    }
    __host__ __device__
    friend complex operator*(const complex v, const complex& _c)  // complex * this
    {
        complex res;
        res.rel = v.rel * _c.rel - v.img * _c.img;
        res.img = v.rel * _c.img + v.img * _c.rel;
        return res;
    }
public:
    real rel;
    real img;
};


using complexd = std::complex<double>;

template<int N, typename T>
class vec
{
private:
    T data[N];
public:
    __host__  __device__
    vec(){}
    __host__  __device__
    vec(T v)
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] = v;
        }
    }
    __host__  __device__
    vec(T v0, T v1)
    {  
        static_assert(N == 2);
        data[0] = v0;
        data[1] = v1;
    }
    __host__  __device__
    vec(T v0, T v1, T v2)
    {  
        static_assert(N == 3);
        data[0] = v0;
        data[1] = v1;
        data[2] = v2;
    }
    __host__  __device__
    vec(const vec &v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] = v[i];
        }
    }
    __host__  __device__
    T &operator[](int i) 
    {
        return data[i];
    }
    __host__  __device__
    const T &operator[](int i) const
    {
        return data[i];
    }
    __host__  __device__
    vec operator+(const T v) const 
    {
        return vec(*this) += v;
    }
    __host__  __device__
    vec operator-(const T v) const 
    {
        return vec(*this) -= v;
    }
    __host__  __device__
    vec operator-() const
    {
        vec res;
        for(int i = 0; i < N; i++)
        {
            res.data[i] = -this->data[i];
        }
        return res;
    }
    __host__  __device__
    vec operator*(const T v) const 
    {
        return vec(*this) *= v;
    }
    __host__  __device__
    vec operator/(const T v) const 
    {
        return vec(*this) /= v;
    }
    __host__  __device__
    vec operator+(const vec v) const 
    {
        return vec(*this) += v;
    }
    __host__  __device__
    vec operator-(const vec v) const 
    {
        return vec(*this) -= v;
    }
    __host__  __device__
    vec operator*(const vec v) const 
    {
        return vec(*this) *= v;
    }
    __host__  __device__
    vec operator/(const vec v) const 
    {
        return vec(*this) /= v;
    }
    __host__  __device__
    const vec &operator+=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] += v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator-=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] -= v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator*=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] *= v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator/=(const T v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] /= v;
        }
        return *this;
    }
    __host__  __device__
    const vec &operator+=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] += v.data[i];
        }
        return *this;
    }
    __host__  __device__
    const vec &operator-=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] -= v.data[i];
        }
        return *this;
    }
    __host__  __device__
    const vec &operator*=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] *= v.data[i];
        }
        return *this;
    }
    __host__  __device__
    const vec &operator/=(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] /= v.data[i];
        }
        return *this;
    }
    __host__  __device__
    constexpr vec& operator=(const vec& v) 
    {
        for(int i = 0; i < N; i++)
        {
            this->data[i] = v.data[i];
        }
        return *this;
    }
    __host__  __device__
    bool operator==(const vec v) 
    {
        for(int i = 0; i < N; i++)
        {
            if(this->data[i] != v.data[i])
                return false;
        }
        return true;
    }
    __host__  __device__
    bool operator!=(const vec v) 
    {
        return !(*this == v);
    }
    __host__  __device__
    friend vec operator*(const T v, const vec& _vec)  // value * this
    {
        vec res;
        for(int i = 0; i < N; i++)
        {
            res[i] = v * _vec[i];
        }
        return res;
    }

    friend std::ostream &operator<<(std::ostream & os, const vec & v) 
    {
        os << "[";
        for(int i = 0; i < N; i++)
        {
            os << v.data[i];
            if(i < N - 1) os << ",";
        }
        os << "]";
        return os;
    }
    __host__  __device__
    T norm()
    {
        T s = 0;
        for(int i = 0; i < N; i++)
        {
            s += data[i] * data[i];
        }
        return s;
    }
    vec abs()
    {
        vec res;
        for(int i = 0; i < N; i++)
        {
            res.data[i] = fabs(this->data[i]);
        }
        return res;
    }
    __host__  __device__
    T r()
    {
        return sqrt(norm());
    }
    T sum()
    {
        T res = 0;
        for(int i = 0; i < N; i++)
        {
            res += this->data[i];
        }
        return res;
    }
    __host__  __device__
    Sph toSph()
    {
        real rho = this->r();
        real theta = rho == 0 ? 0 : std::acos(this->data[2] / rho);
        real phi = std::atan2(this->data[1], this->data[0]);
        return {rho, theta, phi};
    }
    __host__  __device__
    Sph toSph2()
    {
        real rho = this->r();
        real d2r = this->data[2] / rho;
        if(d2r < -1) d2r = -1;
        else if(d2r > 1) d2r = 1;
        real theta = rho == 0 ? 0 : std::acos(d2r);
        real phi = std::atan2(this->data[1], this->data[0]);
        return {rho, theta, phi};
    }
};

typedef vec<2, real> vec2r;
typedef vec<3, real> vec3r;

struct Offset3rPadding
{
    zeta_fmm::real x;
    zeta_fmm::real y;
    zeta_fmm::real z;
    zeta_fmm::real padding;
};

struct Offset3r
{
    zeta_fmm::real x;
    zeta_fmm::real y;
    zeta_fmm::real z;

    Offset3r() : x(0), y(0), z(0) {}
    Offset3r(zeta_fmm::real x_, zeta_fmm::real y_, zeta_fmm::real z_) : x(x_), y(y_), z(z_) {}
};

struct IndexAndOffset3r
{
    int idx;
    Offset3r offset;
    IndexAndOffset3r() : idx(0) {}
    IndexAndOffset3r(int i, zeta_fmm::real x, zeta_fmm::real y, zeta_fmm::real z) : idx(i), offset(x,y,z) {}
    IndexAndOffset3r(int i, zeta_fmm::vec3r v) : idx(i), offset(v[0],v[1],v[2]) {}
};

struct OffsetAndNumber
{
    int offset;
    int number;

    OffsetAndNumber(){}
    OffsetAndNumber(int offset_, int number_) : offset(offset_), number(number_){}

    friend std::ostream &operator<<(std::ostream & os, const OffsetAndNumber & on) 
    {
        os << "[" << on.offset << "," << on.number << "]";
        return os;
    }
};

}