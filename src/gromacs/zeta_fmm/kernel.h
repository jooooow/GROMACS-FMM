#pragma once
#include "type.h"
#include "body.h"
#include "traverser.h"

#define FMM_REG_USE_DW 0

namespace zeta_fmm
{
class GPUKernel
{
public:
    GPUKernel();

    void solve(
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
    );
    
};
}