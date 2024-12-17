#include "traverser.h"
#include "type.h"

zeta_fmm::NearFieldSelector::NearFieldSelector(SelectorType t, real s) : type(t), scale(s)
{
    //std::cout<<"init NearFieldSelector"<<std::endl;
}

bool zeta_fmm::NearFieldSelector::is_near(const Cell3& a, const Cell3& b, vec3r offset)
{
    if(type == SelectorType::MAC)
    {
        double r = (a.center - b.center - offset).r() * scale;
        if(r > a.r + b.r)  // well sepatated
        {
            return 0;
        }
        return 1;
    }
    else if(type == SelectorType::NEARALL)
    {
        return 1;
    }
    else if(type == SelectorType::NEARNONE)
    {
        return 0;
    }
    else
    {
        return 0;
    }
}

zeta_fmm::Traverser::Traverser()
{
    //std::cout<<"init Traversaler"<<std::endl;
    p2p_map.clear();
    depth_map.clear();
    m2l_map.clear();
}

zeta_fmm::Traverser::~Traverser()
{

}

void zeta_fmm::Traverser::traverse(zeta_fmm::CompleteBalancedOctree& tree, NearFieldSelector& nearfield_selector, int images, zeta_fmm::real cycle)
{
    cells = tree.get_cells();
    //std::cout<<"traversal start for tree with "<<cells.size()<<" cells"<<std::endl;

    upward(0);
    if(images == 0)
    {
        //std::cout<<"horizontal"<<std::endl;
        horizontal(nearfield_selector, 0, 0);
    }
    else
    {
        //std::cout<<"horizontal_periodic_near"<<std::endl;
        horizontal_periodic_near(nearfield_selector, 0, 0, cycle);
        //std::cout<<"horizontal_periodic_far"<<std::endl;
        horizontal_periodic_far(nearfield_selector, 0, 0, cycle, images);
    }

    //std::cout<<"traversal finish"<<std::endl;
}

void zeta_fmm::Traverser::upward(
        int this_cell_idx
    )
{
    const zeta_fmm::Cell3& this_cell = cells[this_cell_idx];
    for(int i = 0; i < this_cell.child_info.number; i++)
    {
        upward(this_cell.child_info.offset + i);
    }
    if(this_cell.child_info.number == 0)
    {
        depth_map[this_cell.depth].push_back(IndexAndOffset3r(this_cell_idx, 0, 0, 0));
    }
    else
    {
        depth_map[this_cell.depth].push_back(IndexAndOffset3r(this_cell_idx, 0, 0, 0));
    }
}

void zeta_fmm::Traverser::horizontal(
        NearFieldSelector& nearfield_selector,
        int this_cell_idx, 
        int that_cell_idx, 
        vec3r offset
    )
{
    const zeta_fmm::Cell3& this_cell = cells[this_cell_idx];
    const zeta_fmm::Cell3& that_cell = cells[that_cell_idx];
    if(!nearfield_selector.is_near(this_cell, that_cell, offset))  // well sepatated
    {
        //if(that_cell.body_info.number > 0)
        m2l_map[this_cell_idx].push_back(IndexAndOffset3r(that_cell_idx, offset));
    }
    else
    {
        if(this_cell.child_info.number == 0 && that_cell.child_info.number == 0)
        {
            //if(that_cell.body_info.number > 0)
            p2p_map[this_cell_idx].push_back(IndexAndOffset3r(that_cell_idx, offset));
        }
        else if(this_cell.child_info.number == 0)
        {
            for(int i = 0; i <that_cell.child_info.number; i++)
            {
                horizontal(nearfield_selector, this_cell_idx, that_cell.child_info.offset + i, offset);
            }
        }
        else if(that_cell.child_info.number == 0)
        {
            for(int i = 0; i <this_cell.child_info.number; i++)
            {
                horizontal(nearfield_selector, this_cell.child_info.offset + i, that_cell_idx, offset);
            }
        }
        else
        {
            if(this_cell.r >= that_cell.r)
            {
                for(int i = 0; i <this_cell.child_info.number; i++)
                {
                    horizontal(nearfield_selector, this_cell.child_info.offset + i, that_cell_idx, offset);
                }
            }
            else
            {
                for(int i = 0; i <that_cell.child_info.number; i++)
                {
                    horizontal(nearfield_selector, this_cell_idx, that_cell.child_info.offset + i, offset);
                }
            }
        }
    }
}

void zeta_fmm::Traverser::horizontal_periodic_near(
        NearFieldSelector& nearfield_selector,
        int this_cell_idx, 
        int that_cell_idx,
        zeta_fmm::real cycle
)
{
    for(int pz = -1; pz <= 1; pz++)
    {
        for(int py = -1; py <= 1; py++)
        {
            for(int px = -1; px <= 1; px++)
            {
                horizontal(nearfield_selector, this_cell_idx, that_cell_idx, zeta_fmm::vec3r(px,py,pz) * cycle);
            }   
        }
    }
}

void zeta_fmm::Traverser::horizontal_periodic_far(
    NearFieldSelector& nearfield_selector,
    int this_cell_idx, 
    int that_cell_idx,
    zeta_fmm::real cycle,
    int images
)
{
    // add image cells
    int child_idx = that_cell_idx;
    int idx = cells.size();
    for(int m = 0; m < images - 1; m++)
    {
        zeta_fmm::Cell3 c;
        c.depth = -1 - m;
        c.center = cells[child_idx].center;
        c.child_info = zeta_fmm::OffsetAndNumber(child_idx, 1);
        cells.push_back(c);
        child_idx = cells.size() - 1;
        depth_map[c.depth].push_back(IndexAndOffset3r(child_idx, 0, 0, 0));
    }
    for(int m = 0; m < images - 1; m++)
    {
        for(int pz = -1; pz <= 1; pz++)
        {
            for(int py = -1; py <= 1; py++)
            {
                for(int px = -1; px <= 1; px++)
                {
                    if(px == 0 && py == 0 && pz == 0) 
                        continue;
                    for(int cz = -1; cz <= 1; cz++)
                    {
                        for(int cy = -1; cy <= 1; cy++)
                        {
                            for(int cx = -1; cx <= 1; cx++)
                            {
                                int ox = px * 3 + cx;
                                int oy = py * 3 + cy;
                                int oz = pz * 3 + cz;
                                int icell_idx = cells[idx].child_info.offset;
                                //if(cells[icell_idx].body_info.number > 0)
                                m2l_map[this_cell_idx].push_back(IndexAndOffset3r(icell_idx, zeta_fmm::vec3r(ox,oy,oz) * cycle));
                            }
                        }
                    }
                }
            }
        }
        idx++;
        cycle *= 3;
    }
}

zeta_fmm::Traverser::CellMap zeta_fmm::Traverser::get_depth_map()
{
    return depth_map;
}

zeta_fmm::Traverser::CellMap zeta_fmm::Traverser::get_p2p_map()
{
    return p2p_map;
}

zeta_fmm::Traverser::CellMap zeta_fmm::Traverser::get_m2l_map()
{
    return m2l_map;
}

std::vector<zeta_fmm::Cell3> zeta_fmm::Traverser::get_cells()
{
    return cells;
}