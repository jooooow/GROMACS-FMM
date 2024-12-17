#include "tree.h"
#include <queue>
#include <vector>

zeta_fmm::CompleteBalancedOctree::CompleteBalancedOctree()
{
    //std::cout<<"init complete balance octree"<<std::endl;
    num_body = 0;
    root_cell_idx = -1;
    reg0_cell_idx = -1;
}

void zeta_fmm::CompleteBalancedOctree::build(std::vector<Body3>& bodies, zeta_fmm::real box_r, zeta_fmm::vec3r center, int max_depth)
{
    //std::cout<<"build complete balanced octree"<<std::endl;
    std::queue<int> big_cells;
    num_body = bodies.size();
    Cell3 root;
    root.depth = 0;
    root.r = box_r;
    root.center = center;
    root.child_info = OffsetAndNumber(0,0);
    root.body_info = OffsetAndNumber(0, bodies.size());
    this->cells.push_back(root);
    big_cells.push(0);
    root_cell_idx = 0;

    while(!big_cells.empty())
    {
        int branch_cell_idx = big_cells.front();
        Cell3 branch_cell = this->cells[branch_cell_idx];
        //std::cout<<branch_cell_idx<<","<<branch_cell<<std::endl;
        big_cells.pop();

        if(branch_cell.depth < max_depth)
        {
            int begin = branch_cell.body_info.offset;
            int num = branch_cell.body_info.number;
            int end = branch_cell.body_info.offset + branch_cell.body_info.number - 1;
            //printf("begin=%d, num=%d, end=%d\n", begin, num, end);

            //divide
            int quad_num[8] = {0,0,0,0,0,0,0,0};
            int offset[8] = {0,0,0,0,0,0,0,0};
            int offset_sum = 0;
            for(int i = begin; i <= end; i++)
            {
                //std::cout<<bodies[i]<<std::endl;
                int idx = ((bodies[i].loc[0] > branch_cell.center[0]) << 2) 
                        + ((bodies[i].loc[1] > branch_cell.center[1]) << 1) 
                        + ((bodies[i].loc[2] > branch_cell.center[2]) << 0);
                quad_num[idx]++;
            }
            for(int i = 0; i < 8; i++)
            {
                offset[i] = offset_sum;
                offset_sum += quad_num[i];
            }
            //sort
            std::vector<Body3> bodies_buffer;
            bodies_buffer.resize(num);
            for(int i = begin; i <= end; i++)
            {
                int idx = ((bodies[i].loc[0] > branch_cell.center[0]) << 2) 
                        + ((bodies[i].loc[1] > branch_cell.center[1]) << 1) 
                        + ((bodies[i].loc[2] > branch_cell.center[2]) << 0);
                bodies_buffer[offset[idx]] = bodies[i];
                offset[idx]++;
            }
            //store
            for(int i = 0; i < num; i++)
            {
                bodies[begin + i] = bodies_buffer[i];
            }
            //create leaves
            int insert_offset = this->cells.size();
            this->cells[branch_cell_idx].child_info = OffsetAndNumber(insert_offset,8);
            for(int i = 0; i < 8; i++)
            {
                Cell3 leaf_cell;
                leaf_cell.depth = branch_cell.depth + 1;
                leaf_cell.r = branch_cell.r / 2;
                leaf_cell.child_info = OffsetAndNumber(0,0);
                leaf_cell.body_info = OffsetAndNumber(begin + offset[i] - quad_num[i], quad_num[i]);
                if(i == 0)      leaf_cell.center = branch_cell.center + vec3r(-branch_cell.r/2, -branch_cell.r/2, -branch_cell.r/2);
                else if(i == 1) leaf_cell.center = branch_cell.center + vec3r(-branch_cell.r/2, -branch_cell.r/2,  branch_cell.r/2);
                else if(i == 2) leaf_cell.center = branch_cell.center + vec3r(-branch_cell.r/2,  branch_cell.r/2, -branch_cell.r/2);
                else if(i == 3) leaf_cell.center = branch_cell.center + vec3r(-branch_cell.r/2,  branch_cell.r/2,  branch_cell.r/2);
                else if(i == 4) leaf_cell.center = branch_cell.center + vec3r( branch_cell.r/2, -branch_cell.r/2, -branch_cell.r/2);
                else if(i == 5) leaf_cell.center = branch_cell.center + vec3r( branch_cell.r/2, -branch_cell.r/2,  branch_cell.r/2);
                else if(i == 6) leaf_cell.center = branch_cell.center + vec3r( branch_cell.r/2,  branch_cell.r/2, -branch_cell.r/2);
                else if(i == 7) leaf_cell.center = branch_cell.center + vec3r( branch_cell.r/2,  branch_cell.r/2,  branch_cell.r/2);
                //std::cout<<leaf_cell<<","<<begin + offset[i] - quad_num[i]<<std::endl;
                this->cells.push_back(leaf_cell);
                big_cells.push(insert_offset + i);
            }
        }
    }
}

void zeta_fmm::CompleteBalancedOctree::build_reg(std::vector<Body3>& bodies, zeta_fmm::real rega, int max_depth, zeta_fmm::real cycle)
{
    if(rega > 0)
    {
        add_regcells(max_depth, cycle);
        add_regbodies(bodies, rega);
    }
}

void zeta_fmm::CompleteBalancedOctree::add_regcells(int max_depth, zeta_fmm::real cycle)
{
    int reg_cell_num = std::pow(std::pow(2, max_depth) + 2, 3) - std::pow(2, 3 * max_depth);
    zeta_fmm::real min_cell_r = cells[root_cell_idx].r / std::pow(2, max_depth);
    //printf("reg_cell_num = %d, min_cell_r = %.4f\n", reg_cell_num, min_cell_r);

    std::vector<zeta_fmm::Cell3> reg_cells;

    // add reg0 cell
    zeta_fmm::Cell3 reg0_cell;
    reg0_cell.depth = 0;
    reg0_cell.r = cells[root_cell_idx].r + 2 * min_cell_r;
    reg0_cell.center = cells[root_cell_idx].center;
    reg0_cell.child_info = zeta_fmm::OffsetAndNumber(1, reg_cell_num + 1);
    reg0_cell.body_info = zeta_fmm::OffsetAndNumber(0,num_body);
    reg_cells.push_back(reg0_cell);

    // add reg cells
    int m = get_max_depth();
    //printf("max_depth = %d\n", m);
    for(int cidx = 0; cidx < cells.size(); cidx++)
    {
        zeta_fmm::Cell3* this_cell = &cells[cidx];
        if(this_cell->child_info.number == 0)
        {
            for(int pz = -1; pz <= 1; pz++)
            {
                for(int py = -1; py <= 1; py++)
                {
                    for(int px = -1; px <= 1; px++)
                    {
                        if(px == 0 && py == 0 && pz == 0)
                            continue;
                        zeta_fmm::vec3r reg_cell_center = this_cell->center + zeta_fmm::vec3r(px,py,pz) * cycle;
                        zeta_fmm::vec3r reg_cell_center_dx = (reg_cell_center - cells[root_cell_idx].center).abs() - cells[root_cell_idx].r;
                        if( reg_cell_center_dx[0] <= 1.5 * min_cell_r && 
                            reg_cell_center_dx[1] <= 1.5 * min_cell_r && 
                            reg_cell_center_dx[2] <= 1.5 * min_cell_r)  // need multiple a constant bigger in (1,3), otherwise float error will make this sentence false
                        {
                            zeta_fmm::Cell3 reg_cell;
                            reg_cell.depth = m + 1; // reg cells should be put in the leaf level
                            reg_cell.r = min_cell_r;
                            reg_cell.center = reg_cell_center;
                            reg_cell.child_info = zeta_fmm::OffsetAndNumber(0,0);
                            reg_cell.body_info = zeta_fmm::OffsetAndNumber(num_body,0);
                            reg_cells.push_back(reg_cell);
                        }
                    }
                }
            }
        }
    }

    // shift cells child offset
    for(int i = 0; i < cells.size(); i++)
    {
        zeta_fmm::Cell3& c = cells[i];
        c.depth += 1;
        c.child_info.offset += reg_cell_num + 1;
    }

    // concat reg_cells and cells
    std::vector<zeta_fmm::Cell3> temp_cells;
    temp_cells.insert(temp_cells.end(), reg_cells.begin(), reg_cells.end());
    temp_cells.insert(temp_cells.end(), cells.begin(), cells.end());

    //printf("reg_cells.size = %ld, cells = %ld, temp_cells = %ld\n", reg_cells.size(), cells.size(), temp_cells.size());

    cells = temp_cells;
    root_cell_idx += reg_cell_num + 1;
    reg0_cell_idx = 0;
}

void zeta_fmm::CompleteBalancedOctree::add_regbodies(std::vector<Body3>& bodies, zeta_fmm::real rega)
{
    std::map<int, std::vector<zeta_fmm::Body3>> reg_body_map;
    for(int i = 0; i < cells.size(); i++)
    {
        zeta_fmm::Cell3& ci = cells[i];
        if(ci.child_info.number == 0)
        {
            std::vector<zeta_fmm::Body3> cell_reg_body;
            for(int j = 0; j < cells.size(); j++)
            {
                zeta_fmm::Cell3& cj = cells[j];
                if(i != j && cj.child_info.number == 0)
                {
                    for(int k = cj.body_info.offset; k < cj.body_info.offset + cj.body_info.number; k++)
                    {
                        zeta_fmm::Body3 b = bodies[k];
                        zeta_fmm::vec3r dx_abs = (b.loc - ci.center).abs() - ci.r;
                        if(dx_abs[0] < rega && dx_abs[1] < rega && dx_abs[2] < rega)
                        {
                            cell_reg_body.push_back(b);
                        }
                    }
                }
            }
            if(cell_reg_body.size() > 0)
            {
                reg_body_map[i] = cell_reg_body;
            }
        }
    }
    insert_regbodies_into_bodies(bodies, reg_body_map);
}

void zeta_fmm::CompleteBalancedOctree::insert_regbodies_into_bodies(std::vector<Body3>& bodies, std::map<int, std::vector<zeta_fmm::Body3>>& reg_map)
{
    //std::cout<<"bodies.size() = "<<bodies.size()<<std::endl;
    for(auto m : reg_map)
    {
        int cidx = m.first;
        zeta_fmm::Cell3& c = cells[cidx];
        std::vector<zeta_fmm::Body3>& reg_bodies = m.second;
        int crb_num = reg_bodies.size();
        //printf("cidx = %d, depth = %d, insert %d bodies, offset = %d, number = %d\n", cidx, c.depth, crb_num, c.body_info.offset, c.body_info.number);
        for(int i = 0; i < reg_bodies.size(); i++)
        {
            zeta_fmm::Body3 rb = reg_bodies[i];
            //printf("rb.idx = %d,\n", rb.idx);
        }
        bodies.insert(bodies.begin() + c.body_info.offset + c.body_info.number, reg_bodies.begin(), reg_bodies.end());
        c.body_info.number += crb_num;
        for(int k = 0; k < cells.size(); k++)
        {
            zeta_fmm::Cell3& ck = cells[k];
            if(cidx != k && ck.child_info.number == 0 && ck.body_info.offset >= c.body_info.offset)
            {
                ck.body_info.offset += crb_num;
            }
        }
    }
    //std::cout<<"bodies.size() = "<<bodies.size()<<std::endl;

    // refresh all cells' body, or unnecessary ?
}

int zeta_fmm::CompleteBalancedOctree::get_max_depth()
{
    int m = 0;
    for(int i = 0; i < cells.size(); i++)
    {
        zeta_fmm::Cell3& c = cells[i];
        m = std::max(m, c.depth);  
    }
    return m;
}

std::vector<zeta_fmm::Cell3> zeta_fmm::CompleteBalancedOctree::get_cells()
{
    return cells;
}
