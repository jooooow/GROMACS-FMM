#pragma once
#include "type.h"
#include <vector>

namespace rtfmm
{

struct Body3
{
    int idx;
    real q;
    real p;
    vec3r x;
    vec3r f;
    Body3(){}
    Body3(int idx_, real q_, vec3r x_) : idx(idx_), q(q_), p(0), x(x_), f(0) {}
};

struct ManyBody
{
    int num;
    std::vector<int> idxs;
    std::vector<real> qs;
    std::vector<real> ps;
    std::vector<real> xs;
    std::vector<real> ys;
    std::vector<real> zs;
    std::vector<real> fxs;
    std::vector<real> fys;
    std::vector<real> fzs;
};

struct BodyCompareResult
{
    int num_compared;
    real rmsp;
    real rmsf;
    real l2p;
    real l2f;
    real epot1;
    real epot2;
    real l2e;
    std::string name1;
    std::string name2;
    void show(std::string output_path);
};

using Bodies3 = std::vector<Body3>;

/**
 * @brief Generate random bodies in 3d box.
 * 
 * @param num number of bodies
 * @param r half size of cubic box
 * @param offset box center
 * 
 * @return bodies with zero net charge 
 */
Bodies3 generate_random_bodies(int num, real r, vec3r offset = vec3r(0,0,0), int seed = 0, int zero_netcharge = 1);

/**
 * @brief extract x from bodies
 * 
 * @param bs bodies
 * 
 * @return vector of x
 */
std::vector<vec3r> get_bodies_x(Bodies3& bs, Range range, vec3r offset = vec3r(0,0,0));

/**
 * @brief extract q from bodies
 * 
 * @param bs bodies
 * 
 * @return vector of q
 */
Matrix get_bodies_q(Bodies3& bs, Range range);


void set_boides_p(Bodies3& bs, Matrix& ps, Range range);

void set_boides_f(Bodies3& bs, Matriv& fs, Range range);

void add_boides_p(Bodies3& bs, Matrix& ps, Range range);

void add_boides_f(Bodies3& bs, Matriv& fs, Range range);

void scale_bodies(Bodies3& bs, real scale = 1.0f / (4 * M_PI));

/**
 * @brief Print one body.
 * 
 * @param b body
 */
void print_body(const Body3& b);


/**
 * @brief Print many bodies.
 * 
 * @param bs bodies
 * @param num if -1, print all; else, print min(num,bs.size()) bodies
 * @param offset offset
 */
void print_bodies(const Bodies3& bs, int num = -1, int offset = 0, std::string name = "bodies");


/**
 * @brief compare two Bodies's potential and force to get information(rms-error, l2-error, potential-energy)
 * 
 * @param bs1 bodies1
 * @param bs2 bodies2
 * @param name1 name of bs1
 * @param name2 name of bs2
 * @param num_compare number of comparison(-1 means all)
 * 
 * @return result of body data comparison
 * 
 * @warning bs2 is used to calculate norm
 */
BodyCompareResult compare(const Bodies3& bs1, const Bodies3& bs2, std::string name1, std::string name2, int num_compare = -1);

/**
 * @brief RT
 */
Bodies3 sort_bodies_by_idx(const Bodies3& bs);

ManyBody Bodies2Manybody(const Bodies3& bs);

Bodies3 Manybody2Bodies(const ManyBody& bs);

void dipole_correction(Bodies3& bs, real cycle, int dir = -1);

}