#pragma once

#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>


/* Matrix inversion routine.
      Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool invert_matrix (const boost::numeric::ublas::matrix<T>& input, boost::numeric::ublas::matrix<T>& inverse) {
  typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  boost::numeric::ublas::matrix<T> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  // perform LU-factorization
  int res = boost::numeric::ublas::lu_factorize(A,pm);
        if( res != 0 ) return false;
  // create identity matrix of "inverse"
  inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));
  // backsubstitute to get the inverse
  boost::numeric::ublas::lu_substitute(A, pm, inverse);
  return true;
}

/*
  Newton-Raphson root finder
  Params
  ------
    u0   := first guess boost::numeric::ublas::matrix<T>
    f(u) := function with format 
            boost::numeric::ublas::matrix<T> (*foo)(const boost::numeric::ublas::matrix<T>& u)

  u_t+1 = u_t - f(u_t)/df(u_t) 

  f(u) = 0
 */
template<class T, typename U, typename M>
bool newton_raphson (boost::numeric::ublas::vector<T>& guess,
                     U&& f,
                     M&& df,
                     boost::numeric::ublas::vector<T>& root,
                     double MAX_ERROR = 1E-15,
                     double MAX_ITER = 1E+5) {
  //std::cout <<"Initial guess" << guess << "\n";
  boost::numeric::ublas::vector<T> f_mat(guess.size()); // copy constructor
  boost::numeric::ublas::vector<T> guess_0 = guess; // copy constructor

  auto nrow = f_mat.size();
  auto ncol = nrow;
  boost::numeric::ublas::matrix<T> df_mat (nrow, ncol);
  boost::numeric::ublas::matrix<T> df_inv (nrow, ncol);
  
  if (!f(guess, f_mat))
    return false;
  root = guess;
  long counter=0;
  while (boost::numeric::ublas::norm_2(f_mat) > MAX_ERROR && counter < MAX_ITER)
  {
    if (!df(root, df_mat))
      return false;
    if (!invert_matrix(df_mat, df_inv))
      return false;
    boost::numeric::ublas::axpy_prod(df_inv, -f_mat, root, false);
    
    //if (counter % 100 == 0) std::cout << "Iter " << counter << ": " << f_mat << "; " << df_inv << "; " << boost::numeric::ublas::prod(df_inv, f_mat) << "; " << root << "\n";
    counter++;
    if (!f(root, f_mat))
      return false;
  }
  //root = guess;
  // if (counter > MAX_ITER) std::cout << "WARNING: Max number of iterations exceeded! Initial guess: " << guess << "; Current value: " << root << "\n";
  // else std::cout << "INFO: Initial guess: " << guess << "; Converged value: " << root << "\n";
  return true;
}



template<class T, typename U, typename M>
bool newton_raphson (
  boost::numeric::ublas::vector<T>& u1,
  boost::numeric::ublas::vector<T>& u2,
  U&& f,
  M&& df,
  boost::numeric::ublas::vector<T>& root,
  double MAX_ERROR = 1E-15,
  double MAX_ITER = 1E+5
)
{
  boost::numeric::ublas::vector<T> f_mat(u1.size());
  boost::numeric::ublas::vector<T> guess(u1.size());
  boost::numeric::ublas::vector<T> du_old(u1.size());
  boost::numeric::ublas::vector<T> du(u1.size());
  guess = 0.5*(u1+u2);
  du_old = u2-u1;
  du = du_old;

  std::cout <<"Initial guess" << guess << "\n";

  auto nrow = f_mat.size();
  auto ncol = nrow;
  boost::numeric::ublas::matrix<T> df_mat (nrow, ncol);
  boost::numeric::ublas::matrix<T> df_inv (nrow, ncol);
  
  if (!f(guess, f_mat))
    return false;
  if (!df(guess, df_mat))
      return false;

  root = guess;
  long counter=0;
  while (boost::numeric::ublas::norm_2(f_mat) > MAX_ERROR && counter < MAX_ITER)
  {
    if ((
        boost::numeric::ublas::inner_prod(
          (boost::numeric::ublas::prod(
            df_mat, 
            (root-u2)
            )-f_mat
          ),
          (boost::numeric::ublas::prod(
            df_mat, 
            (root-u1)
            )-f_mat
          )
        ) > 0.0) || (boost::numeric::ublas::norm_2(2.0*f_mat) > boost::numeric::ublas::norm_2(boost::numeric::ublas::prod(df_mat, du_old))))
    {
      std::cout << "Bissection!\n";
      du_old=du;
      du=0.5*(u2-u1);
      root=u1+du;
      if (boost::numeric::ublas::norm_2(u1-root) < MAX_ERROR) return true;
    }
    else { 
      std::cout << "Newton-Raphson!\n";
      du_old=du;
      if (!invert_matrix(df_mat, df_inv))
        return false;
      
      du=boost::numeric::ublas::prod(df_inv, f_mat);
      auto temp=root;
      root -= du;
      if (boost::numeric::ublas::norm_2(temp-root) < MAX_ERROR) return false;
    }
    

    if (boost::numeric::ublas::inner_prod((u1-root),(root-u2)) < 0.0)
      throw("Jumped out of brackets in newton_raphson");

    if (!f(root, f_mat))
      return false;
    if (!df(root, df_mat))
      return false;
    if (!invert_matrix(df_mat, df_inv))
      return false;

    //if (counter % 100 == 0) 
    std::cout << "Iter " << counter << ": f = " << f_mat << "; 1/df = " << df_inv << "; f/df = " << boost::numeric::ublas::prod(df_inv, f_mat) << "; u = " << root << "\n";
    counter++;

    if (boost::numeric::ublas::norm_2(f_mat) < 0.0)
     u1=root;
    else
     u2=root;
  }
  if (counter > MAX_ITER) std::cout << "WARNING: Max number of iterations exceeded! Initial guess: " << guess << "; Current value: " << root << "\n";
  else std::cout << "INFO: Initial guess: " << guess << "; Converged value: " << root << "\n";
  return true;
}