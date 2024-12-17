#include "solver.h"
#include "gromacs/math/units.h"

zeta_fmm::GPUSolver::GPUSolver(FILE* f) : Solver(f)
{
    //fprintf(fplog, "init non-para GPUSolver hjm\n");
}

zeta_fmm::GPUSolver::GPUSolver(
    FILE* f,
    std::vector<Body3>& bs,
    int P_,
    zeta_fmm::real box_r_,
    zeta_fmm::vec3r center_,
    int max_depth_,
    int images_,
    zeta_fmm::real cycle_,
    zeta_fmm::real rega_,
    int verbose_,
    int dummy_
) : Solver(f, bs, "gpuFMM"), P(P_), box_r(box_r_), center(center_), max_depth(max_depth_), images(images_), cycle(cycle_), rega(rega_), verbose(verbose_), dummy(dummy_)
{
    merged_num_body = bs.size();
    merged_bodies = bs;
    if(verbose)
    {
        printf("[verbose] bs.size = %d\n", merged_num_body);
    }
}

zeta_fmm::GPUSolver::~GPUSolver()
{
    
}

void zeta_fmm::GPUSolver::prepare_maps()
{
    cudaFree(0);
    
    // build tree from bodies
    zeta_fmm::CompleteBalancedOctree tree;
    tree.build(bodies, box_r, center, max_depth);
    tree.build_reg(bodies, rega, max_depth, cycle); // caution : this step will extend the bodies(thus bodies's real size changes)

    // set near field policy
    zeta_fmm::NearFieldSelector nearfield_selector(zeta_fmm::NearFieldSelector::SelectorType::MAC, 0.4);

    // traverse tree to get interaction map
    zeta_fmm::Traverser traverser;
    traverser.traverse(tree, nearfield_selector, images, cycle); // caution : when images > 0, this function will add extra cells
    cells = traverser.get_cells();
    depth_map = traverser.get_depth_map();
    p2p_map = traverser.get_p2p_map();
    m2l_map = traverser.get_m2l_map();

    if(verbose)
    {
        printf("p2p_map.size = %ld\n", p2p_map.size());
        printf("m2l_map.size = %ld\n", m2l_map.size());
    }
}

void zeta_fmm::GPUSolver::solve()
{
    //fprintf(fplog, "GPUSolver solve hjm\n");

    // prepare interaction maps
    prepare_maps();

    // prepare body & cells
    int real_num_body = bodies.size();
    int real_num_cell = cells.size();
    
    if(verbose)
    {
        printf("merged_num_body = %d, real_num_body = %d, real_num_cell = %d\n", merged_num_body, real_num_body, real_num_cell);
    }
    
    Body3* g_bodies;
    cudaMalloc(&g_bodies, real_num_body * sizeof(Body3));
    cudaMemcpy(g_bodies, bodies.data(), real_num_body * sizeof(Body3), cudaMemcpyHostToDevice);

    Body3* g_merged_bodies;
    cudaMalloc(&g_merged_bodies, merged_num_body * sizeof(Body3));
    cudaMemcpy(g_merged_bodies, merged_bodies.data(), merged_num_body * sizeof(Body3), cudaMemcpyHostToDevice);

    Cell3* g_cells;
    cudaMalloc(&g_cells, real_num_cell * sizeof(Cell3));
    cudaMemcpy(g_cells, cells.data(), real_num_cell * sizeof(Cell3), cudaMemcpyHostToDevice);

    // prepare pole
    int pole_size_eachcell = get_pole_size_eachcell(P);
    zeta_fmm::complex* g_Ms;
    cudaMalloc(&g_Ms, real_num_cell * pole_size_eachcell * sizeof(complex));
    cudaMemset(g_Ms, 0, real_num_cell * pole_size_eachcell * sizeof(complex));
    zeta_fmm::complex* g_Ls;
    cudaMalloc(&g_Ls, real_num_cell * pole_size_eachcell * sizeof(complex));
    cudaMemset(g_Ls, 0, real_num_cell * pole_size_eachcell * sizeof(complex));

    // find max_depth
    int real_max_depth = get_real_max_depth();
    if(verbose)
    {
        std::cout<<"real_max_depth = "<<real_max_depth<<std::endl;
    }

    // get leaf cell (P2M & L2P)
    std::vector<int> leaf_list = get_leaf_list(); // max-depth-cells
    int leaf_cell_num = leaf_list.size();
    if(verbose)
    {
        std::cout<<"leaf_cell_num = "<<leaf_cell_num<<std::endl;
    }
    int* g_leaf_cells;
    cudaMalloc(&g_leaf_cells, leaf_cell_num * sizeof(int));
    cudaMemcpy(g_leaf_cells, leaf_list.data(), leaf_cell_num * sizeof(int), cudaMemcpyHostToDevice);

    // store sorted branch cell (M2M & L2L)
    std::vector<int> branch_cells; // cells except max-depth-cells and image-cells (sorted from max-depth to 0)
    std::vector<OffsetAndNumber> level_infos;
    int level_offset = 0;
    for(auto map = depth_map.rbegin(); map != depth_map.rend(); map++)
    {
        if(map->first == real_max_depth || map->first < 0) continue;
        std::vector<IndexAndOffset3r>& cs = map->second;
        level_infos.push_back(OffsetAndNumber(level_offset, cs.size()));
        for(size_t i = 0; i < cs.size(); i++)
        {
            branch_cells.push_back(cs[i].idx);
        }
        level_offset += cs.size();
    }
    int* g_branch_cells;
    cudaMalloc(&g_branch_cells, branch_cells.size() * sizeof(int));
    cudaMemcpy(g_branch_cells, branch_cells.data(), branch_cells.size() * sizeof(int), cudaMemcpyHostToDevice);

    // get image cells
    std::vector<int> img_cells; // image-cells
    for(auto map = depth_map.rbegin(); map != depth_map.rend(); map++)
    {
        if(map->first < 0)
        {
            img_cells.push_back(map->second[0].idx);
        }
    }
    int img_cells_num = img_cells.size();
    if(verbose)
    {
       printf("img_cells_num = %d\n", img_cells_num); 
    }
    int* g_img_cells;
    cudaMalloc(&g_img_cells, img_cells_num * sizeof(int));
    cudaMemcpy(g_img_cells, img_cells.data(), img_cells_num * sizeof(int), cudaMemcpyHostToDevice);

    // generate p2p matrix (P2P)
    int max_p2p_source_num = 0;
    for(auto map : p2p_map)
    {
        max_p2p_source_num = std::max(max_p2p_source_num, (int)map.second.size());
    }
    int p2p_matrix_row = p2p_map.size();
    int p2p_matrix_col = max_p2p_source_num + 2;
    int* c_p2p_matrix = new int[p2p_matrix_row * p2p_matrix_col];
    Offset3rPadding* c_p2p_offset_matrix = new Offset3rPadding[p2p_matrix_row * (p2p_matrix_col - 2)];
    int i = 0;
    for(auto map : p2p_map)
    {
        c_p2p_matrix[i * p2p_matrix_col + 0] = map.first;
        c_p2p_matrix[i * p2p_matrix_col + 1] = map.second.size();
        for(size_t j = 0; j < map.second.size(); j++)
        {
            c_p2p_matrix[i * p2p_matrix_col + 2 + j] = map.second[j].idx;
            Offset3rPadding offset;
            offset.x = map.second[j].offset.x;
            offset.y = map.second[j].offset.y;
            offset.z = map.second[j].offset.z;
            c_p2p_offset_matrix[i * (p2p_matrix_col - 2) + j] = offset;
        }
        i++;
    }
    int* g_p2p_matrix;
    cudaMalloc(&g_p2p_matrix, p2p_matrix_row * p2p_matrix_col * sizeof(int));
    cudaMemcpy(g_p2p_matrix, c_p2p_matrix, p2p_matrix_row * p2p_matrix_col * sizeof(int),cudaMemcpyHostToDevice);
    Offset3rPadding* g_p2p_offset_matrix;
    cudaMalloc(&g_p2p_offset_matrix, p2p_matrix_row * (p2p_matrix_col - 2) * sizeof(Offset3rPadding));
    cudaMemcpy(g_p2p_offset_matrix, c_p2p_offset_matrix, p2p_matrix_row * (p2p_matrix_col - 2) * sizeof(Offset3rPadding),cudaMemcpyHostToDevice);

    // prepare m2l matrix (M2L)
    int max_m2l_source_num = 0;
    for(auto map : m2l_map)
    {
        max_m2l_source_num = std::max(max_m2l_source_num, (int)map.second.size());
    }
    int m2l_matrix_row = m2l_map.size();
    int m2l_matrix_col = max_m2l_source_num + 2;
    int* c_m2l_matrix = new int[m2l_matrix_row * m2l_matrix_col];
    Offset3rPadding* c_m2l_offset_matrix = new Offset3rPadding[m2l_matrix_row * (m2l_matrix_col - 2)];
    i = 0;
    for(auto map : m2l_map)
    {
        c_m2l_matrix[i * m2l_matrix_col + 0] = map.first;
        c_m2l_matrix[i * m2l_matrix_col + 1] = map.second.size();
        for(size_t j = 0; j < map.second.size(); j++)
        {
            c_m2l_matrix[i * m2l_matrix_col + 2 + j] = map.second[j].idx;
            Offset3rPadding offset;
            offset.x = map.second[j].offset.x;
            offset.y = map.second[j].offset.y;
            offset.z = map.second[j].offset.z;
            c_m2l_offset_matrix[i * (m2l_matrix_col - 2) + j] = offset;
        }
        i++;
    }
    int* g_m2l_matrix;
    cudaMalloc(&g_m2l_matrix, m2l_matrix_row * m2l_matrix_col * sizeof(int));
    cudaMemcpy(g_m2l_matrix, c_m2l_matrix, m2l_matrix_row * m2l_matrix_col * sizeof(int),cudaMemcpyHostToDevice);
    Offset3rPadding* g_m2l_offset_matrix;
    cudaMalloc(&g_m2l_offset_matrix, m2l_matrix_row * (m2l_matrix_col - 2) * sizeof(Offset3rPadding));
    cudaMemcpy(g_m2l_offset_matrix, c_m2l_offset_matrix, m2l_matrix_row * (m2l_matrix_col - 2) * sizeof(Offset3rPadding),cudaMemcpyHostToDevice);

    // solve FMM
    //TIME_BEGIN(gpuFMM_kernels);
    gpu_kernel.solve(
        P, rega,
        g_bodies, real_num_body, 
        g_merged_bodies, merged_num_body,
        g_cells, real_num_cell,
        g_Ms, g_Ls,
        g_leaf_cells, leaf_cell_num, 
        g_branch_cells, level_infos,
        g_img_cells, img_cells_num, cycle,
        g_p2p_matrix, g_p2p_offset_matrix, p2p_matrix_row, p2p_matrix_col,
        g_m2l_matrix, g_m2l_offset_matrix, m2l_matrix_row, m2l_matrix_col,
        dummy
    );
    cudaDeviceSynchronize();
    //TIME_END(gpuFMM_kernels);

    // store result
    if(rega == 0)
    {
        cudaMemcpy(bodies.data(), g_bodies, real_num_body * sizeof(Body3), cudaMemcpyDeviceToHost);
    }
    else
    {
        cudaMemcpy(merged_bodies.data(), g_merged_bodies, merged_num_body * sizeof(Body3), cudaMemcpyDeviceToHost);
        bodies = merged_bodies;
    }

    // convert -E to F
    for(size_t i = 0; i < bodies.size(); i++)
    {
        zeta_fmm::Body3& b = bodies[i];
        b.p *= gmx::c_one4PiEps0;
        b.f *= -b.q * gmx::c_one4PiEps0;
    }

    // memory dellocate
    cudaFree(g_bodies);
    cudaFree(g_merged_bodies);
    cudaFree(g_cells);
    cudaFree(g_leaf_cells);
    cudaFree(g_branch_cells);
    cudaFree(g_img_cells);
    cudaFree(g_p2p_matrix);
    cudaFree(g_p2p_offset_matrix);
    cudaFree(g_m2l_matrix);
    cudaFree(g_m2l_offset_matrix);
    cudaFree(g_Ms);
    cudaFree(g_Ls);

    // check error
    CHECK_CUDA_ERR(overall);

    delete[] c_p2p_matrix;
    delete[] c_p2p_offset_matrix;
    delete[] c_m2l_matrix;
    delete[] c_m2l_offset_matrix;
}

std::vector<int> zeta_fmm::GPUSolver::get_leaf_list()
{
    std::vector<IndexAndOffset3r>& leaf_cell_list = depth_map[get_real_max_depth()];
    std::vector<int> leaf_list;
    for(size_t i = 0; i < leaf_cell_list.size(); i++)
    {
        IndexAndOffset3r cell = leaf_cell_list[i];
        leaf_list.push_back(cell.idx);
    }
    return leaf_list;
}

int zeta_fmm::GPUSolver::get_real_max_depth()
{
    int m = 0;
    for(auto map = depth_map.rbegin(); map != depth_map.rend(); map++)
    {
        m = std::max(m, map->first);
    }
    return m;
}
