#pragma once

#include <algorithm>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas::vector<double> uvector;

// Solution Class
class DVector
{
private:
  std::vector<double> q;

public:
  DVector() { this->q = std::vector<double>(4, 0.0); };
  DVector(const std::vector<double> &sol) : q(sol){};
  DVector(const DVector &rhs) { 
    // this->q = std::vector<double>(4, 0.0);
    // this->q = {0.0, 0.0, 0.0, 0.0};
    //this->q.clear();
    // this->q.resize(rhs.q.size());
    this->q = rhs.q;
    // std::copy(rhs.q.begin(), rhs.q.end(), this->q.begin());
  };
  DVector(const uvector &rhs) { 
    auto size = rhs.size();
    this->q = std::vector<double>(size, 0.0);
    for (auto i=0; i<size; i++)
      this->q[i] = rhs(i);
  };
  ~DVector() { this->q.clear(); };

  uvector as_uvector(void)
  {
    auto size = this->q.size();
    uvector uQ = uvector(size);
    for (auto i=0; i<size; i++)
      uQ(i) = this->q[i];
    return uQ;
  }

  std::size_t size(void)
  {
    return this->q.size();
  }

    /* Iterators (delegation) */
  typedef std::vector<double>::iterator iterator;
  typedef std::vector<double>::const_iterator const_iterator;
  iterator begin() { return q.begin(); }
  iterator end() { return q.end(); }

  /* Access ith element operator */
  double operator[](const unsigned int &i) const
  {
    return this->q[i];
  }

  double &operator[](const unsigned int &i)
  {
    return this->q[i];
  }

  /* Assignment = operator*/
  // Assigning a scalar
  template <typename T>
  DVector &operator=(const T &scalar)
  {
    auto q_size = this->q.size();

    if (q_size > 0)
    {
      // for (unsigned int i = 0; i < this->q.size(); i++)
      // {
      //   this->q[i] = (double)scalar;
      // }
      std::transform(
        this->q.begin(), 
        this->q.end(), 
        this->q.begin(), 
        [&scalar](auto& c){return scalar;}
      );
    }
    else
    {
      this->q.reserve(1);
      this->q[0] = (double)scalar;
    }
    return *this;
  }
  // Assigning a DVector
  DVector &operator=(const DVector &rhs)
  {
    // auto rhs_size = rhs.q.size();

    this->q.clear();
    this->q.resize(rhs.q.size());

    this->q = rhs.q; // copy
    // std::copy(rhs.q.begin(), rhs.q.end(),this->q.begin());
    return *this;
  }

  // Assigning a std::vector
  DVector &operator=(const std::vector<double> &rhs)
  {
    auto rhs_size = rhs.size();

    this->q.clear();
    this->q.resize(rhs_size);

    this->q = rhs; // copy

    return *this;
  }

  /* ADDITION AND SUBTRACTION */
  // +=
  template <typename T>
  DVector &operator+=(const T &scalar)
  {
    // for (unsigned int index = 0; index < this->q.size(); index++)
    // {
    //   this->q[index] += scalar;
    // }
    std::transform(
        this->q.begin(), 
        this->q.end(), 
        this->q.begin(), 
        [&scalar](auto& c){return c+scalar;}
      );
    return *this;
  }

  // -=
  template <typename T>
  DVector &operator-=(const T &scalar)
  {
    // for (unsigned int index = 0; index < this->q.size(); index++)
    // {
    //   this->q[index] -= (double)scalar;
    // }
    std::transform(
      this->q.begin(), this->q.end(), 
      this->q.begin(), [&scalar](auto& c){return c-scalar;}
    );
    return *this;
  }

  DVector &operator+=(const DVector &rhs)
  {
    if (this->size() < rhs.q.size())
    {
      this->q.resize(rhs.q.size());
    }
    // unsigned int index = 0;
    // for (auto &val : rhs.q)
    // {
    //   index++;
    //   this->q[index - 1] += val;
    // }
    std::transform(
      this->q.begin(), this->q.end(), rhs.q.begin(), 
      this->q.begin(), std::plus<double>()
    );
    return *this;
  }

  // -=
  DVector &operator-=(const DVector &rhs)
  {
    unsigned int index = 0;
    for (auto &val : rhs.q)
    {
      index++;
      this->q[index - 1] -= val;
    }

    return *this;
  }

  // DVector addition
  DVector operator+(const DVector &rhs)
  {
    DVector result{*this};
    return result += rhs;
  }

  // DVector subtraction
  DVector operator-(const DVector &rhs)
  {
    DVector result(*this);
    return result -= rhs;
  }

  /* MULTIPLICATIONS */
  // *=
  DVector &operator*=(const DVector &rhs)
  {
    unsigned int index = 0;
    for (auto &val : rhs.q)
    {
      index++;
      this->q[index - 1] *= val;
    }

    return *this;
  }

  template <typename T>
  DVector &operator*=(const T &scalar)
  {
    // for (unsigned int index = 0; index < this->q.size(); index++)
    // {
    //   this->q[index] *= (double)scalar;
    // }
    std::transform(
      this->q.begin(), this->q.end(), 
      this->q.begin(), [&scalar](auto& c){return c*scalar;}
    );
    return *this;
  }

  // DVector multiplication
  DVector operator*(const DVector &rhs)
  {
    DVector result(*this);
    return result *= rhs;
  }

  DVector operator-() const // in order to be able to apply it on const objects
  {
    DVector result{*this};
    return result *= -1;
  }
};

// Scalar multiplication
template <typename T>
DVector operator*(const DVector &rhs, const T &scalar)
{
  DVector result{rhs};
  return result *= scalar;
}

// DVector addition
template <typename T>
DVector operator+(const DVector &rhs, const T &scalar)
{
  DVector result{rhs};
  return result += scalar;
}

// DVector subtraction
template <typename T>
DVector operator-(const DVector &rhs, const T &scalar)
{
  DVector result{rhs};
  return result -= scalar;
}

template <typename T>
DVector operator*(const T &scalar, const DVector &rhs)
{
  DVector lhs{rhs * scalar};
  return lhs;
}

template <typename T>
DVector operator+(const T &scalar, const DVector &rhs)
{
  DVector lhs{rhs + scalar};
  return lhs;
}

template <typename T>
DVector operator-(const T &scalar, const DVector &rhs)
{
  DVector lhs{-rhs + scalar};
  return lhs;
}