#pragma once

#include <DVector.h>
#include <Weights.h>

Weights::Weights(int sp_size, int fp_size)
{
  this->fp = fp_size;
  this->sp = sp_size;
  this->Mxfp = MatrixDouble(fp_size, sp_size);
  this->Myfp = MatrixDouble(fp_size, sp_size);
  this->Mxsp = MatrixDouble(sp_size, fp_size);
  this->Mysp = MatrixDouble(sp_size, fp_size);
}

Weights::~Weights()
{
}

void Weights::prod(
  int type, 
  std::vector<DVector>& Q,
  std::vector<DVector>& results
)
{
  std::size_t num_nodes = Q.size();
  std::size_t num_properties = Q[0].size();

  for (std::size_t p_index=0; p_index<num_properties; p_index++)
  {
    uvector u = uvector(num_nodes);
    // read
    for (std::size_t n_index=0; n_index<num_nodes; n_index++)
      u(n_index) = Q[n_index][p_index];

    //calculate 
    uvector result = this->_prod(type, u);

    // write
    for (std::size_t n_index=0; n_index<results.size(); n_index++)
      results[n_index][p_index] = result(n_index);
  }
  
  //return results;
}

void Weights::prod(
  int type, std::vector<DVector>& F, 
  std::vector<double>& J,
  std::vector<DVector>& results
) {
  std::size_t num_nodes = F.size();
  std::size_t num_properties = F[0].size();

  // std::vector<DVector> results;
  // std::size_t nrows = this->_get_nrows(type);
  // results.clear();
  // results.resize(nrows); // number of nodes
  // for (auto &res : results)
  //   res = DVector(std::vector<double>(num_properties, 0.0));

  for (std::size_t p_index=0; p_index<num_properties; p_index++)
  {
    uvector u = uvector(num_nodes);
    // read
    for (std::size_t n_index=0; n_index<num_nodes; n_index++)
      u(n_index) = F[n_index][p_index]*J[n_index];

    //calculate 
    uvector result = this->_prod(type, u);

    // write
    for (std::size_t n_index=0; n_index<results.size(); n_index++)
      results[n_index][p_index] = result(n_index);
  }
  
  // return results;
}

uvector Weights::_prod(int type, uvector u)
{
  switch (type)
  {
  case 0: // MxFP
    return boost::numeric::ublas::prod(this->Mxfp, u);
  case 1: // MyFP
    return boost::numeric::ublas::prod(this->Myfp, u);
  case 2: // MxSP
    return boost::numeric::ublas::prod(this->Mxsp, u);
  case 3: // MySP
    return boost::numeric::ublas::prod(this->Mysp, u);
  }
  return uvector{};
}

size_t Weights::_get_nrows(int type)
{
  switch (type)
  {
  case 0: // MxFP
    return this->Mxfp.size1();
  case 1: // MyFP
    return this->Myfp.size1();
  case 2: // MxSP
    return this->Mxsp.size1();
  case 3: // MySP
    return this->Mysp.size1();
  }
  return 0;
}

void Weights::set(int type, uint i, uint j, double value)
{
  switch (type)
  {
  case 0: // MxFP
    this->Mxfp(i, j) = value;
    break;
  case 1: // MyFP
    this->Myfp(i, j) = value;
    break;
  case 2: // MxSP
    this->Mxsp(i, j) = value;
    break;
  case 3: // MySP
    this->Mysp(i, j) = value;
    break;
  }
}