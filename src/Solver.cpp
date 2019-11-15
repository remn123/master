#include "../include/Solver.h"


SD::SD(int order)
{
  
}

SD::~SD()
{

}

// 0)
void SD::setup(void)
{
  this->calc_sp_nodes();
  this->calc_fp_nodes();
}

// 1)
void SD::initial_condition (Element& e)
{
  //pass
}


// 2)
void SD::calc_high_order_nodes (Element& e)
{
  //pass
}

// 3)
void SD::interpolate_sp2fp (Element& e)
{
  //pass
}


// 4)
void SD::riemann_solver (Element& e)
{
  //pass
}


// 5)
void SD::interpolate_fp2sp (Element& e)
{
  //pass
}


// 6)
void SD::residue (Element& e)
{
  //pass
}

