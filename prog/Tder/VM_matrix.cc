//---------------------------------------------------------------------
/**
 * @file
 * @brief   matrix class for vandermonde
 * @ingroup Group_Vandermonde
 * @author  Noriyoshi ISHII
 * @since   Thu May 30 16:18:38 JST 2013
 *
 * Copyright (C) 2013 by Noriyoshi Ishii, All Rights Reserved.
 *
 */
//---------------------------------------------------------------------

#ifndef NDEBUG
#define NDEBUG 0
#endif

#include "VM.h"

//---------------------------------------------------------------------
/**
 * @brief a constructor
 */
//---------------------------------------------------------------------

VM::matrix::matrix(int n)
  : dim(n), a(new double[n*n])
{ for(int i = 0; i < dim*dim; i++) a[i] = 0.0; }

//---------------------------------------------------------------------
/**
 * @brief the copy constructor
 */
//---------------------------------------------------------------------

VM::matrix::matrix(const matrix& other)
  : dim(other.dim), a(new double[other.dim*other.dim])
{ for(int i = 0; i < dim*dim; i++) a[i] = other.a[i]; }

//---------------------------------------------------------------------
/**
 * @brief the assingment operator
 */
//---------------------------------------------------------------------

VM::matrix& VM::matrix::operator=(const matrix& other)
{
  if(a) delete[] a;
  dim = other.dim;
  a   = new double[dim*dim];

  for(int i = 0; i < dim*dim; i++) a[i] = other.a[i];
  
  return *this;
}

//---------------------------------------------------------------------
/**
 * @breif += operation
 */
//---------------------------------------------------------------------

VM::matrix& VM::matrix::operator+=(const matrix& other)
{
  check_size(other, "operator+=");
  for(int i = 0; i < dim*dim; i++) a[i] += other.a[i];
  return *this;
}

//---------------------------------------------------------------------
/**
 * @brief -= operation
 */
//---------------------------------------------------------------------

VM::matrix& VM::matrix::operator-=(const matrix& other)
{
  check_size(other, "operator-=");
  for(int i = 0; i < dim*dim; i++) a[i] -= other.a[i];
  return *this;
}

//---------------------------------------------------------------------
/**
 * @brief matrix product *= operation
 */
//---------------------------------------------------------------------

VM::matrix& VM::matrix::operator*=(const matrix& other)
{
  check_size(other, "operator*=");

  matrix self(*this);

  for(  int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      double sum = 0.0;
      for(int k = 0; k < dim; k++){
	sum += self(i,k) * other(k,j);
      }
      (*this)(i,j) = sum;
    }
  }
  return *this;
}

//---------------------------------------------------------------------
/**
 * @brief *= operation (double)
 */
//---------------------------------------------------------------------

VM::matrix& VM::matrix::operator*=(double x)
{
  for(int i = 0; i < dim*dim; i++) a[i] *= x;
  return *this;
}
