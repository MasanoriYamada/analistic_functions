//---------------------------------------------------------------------
/**
 * @file
 * @brief
 * @ingroup Group_Vandermonde
 * @author  Noriyoshi Ishii
 * @since   Thu May 30 16:49:53 JST 2013
 *
 * Copyright (C) 2013 by Noriyoshi Ishii, All Rights Reserved.
 *
 */
//---------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

#include "VM.h"

using namespace VM;

//---------------------------------------------------------------------
// shitauke
//---------------------------------------------------------------------

static double power(double x,int n)
{
  double ret = 1.0;
  for(int i = 0; i < n; i++){
    ret *= x;
  }
  return ret;
}

static double keisu(int n,int length)
{
  double ret=1.0;
  for(int i = 0; i < length; i++){
    ret *= double(n - i);
  }
  return ret;
}

//---------------------------------------------------------------------
// main part
//---------------------------------------------------------------------


//---------------------------------------------------------------------
/**
 * @brief the destructor
 */
//---------------------------------------------------------------------

vandermonde::~vandermonde()
{ if(Mder) delete[] Mder; }

//---------------------------------------------------------------------
/**
 * @brief the default constructor
 */
//---------------------------------------------------------------------

vandermonde::vandermonde()
  : Ndeg(0), Vinv(0), Mder(NULL)
{}

//---------------------------------------------------------------------
/**
 * @brief a constructor
 */
//---------------------------------------------------------------------

vandermonde::vandermonde(int n)
  : Ndeg(n), Vinv(n), Mder(new matrix[n])
{
  for(int i = 0; i < Ndeg; i++){
    Mder[i] = matrix(n);
  }
  double x[Ndeg];
  for(int i = 0; i < Ndeg; i++) x[i] = double(i);
  construct_vinv(x);
  construct_mder(x);
}

//---------------------------------------------------------------------
/**
 * @brief the copy constructor
 */
//---------------------------------------------------------------------

vandermonde::vandermonde(const vandermonde& other)
  : Ndeg(other.Ndeg), Vinv(other.Vinv), Mder(new matrix[other.Ndeg])
{
  for(int i = 0; i < Ndeg; i++){
    Mder[i] = other.Mder[i];
  }
}

//---------------------------------------------------------------------
/**
 * @brief the assignment operator
 */
//---------------------------------------------------------------------

vandermonde& vandermonde::operator=(const vandermonde& other)
{
  Ndeg = other.Ndeg;
  Vinv = other.Vinv;
  if(Mder) delete[] Mder;
  Mder = new matrix[Ndeg];
  for(int i = 0; i < Ndeg; i++){
    Mder[i] = other.Mder[i];
  }
  return *this;
}

//---------------------------------------------------------------------
/**
 * @brief construction of Vinv
 */
//---------------------------------------------------------------------

void vandermonde::construct_vinv(double *x)
{
  matrix Linv(Ndeg);
  matrix Uinv(Ndeg);

  construct_Linv(Linv, x);
  construct_Uinv(Uinv, x);
  Vinv = Uinv * Linv;
}

//---------------------------------------------------------------------
/**
 * @brief construction of Mder
 */
//---------------------------------------------------------------------

void vandermonde::construct_mder(double *x)
{
  for(int nth = 0; nth < Ndeg; nth++){
    matrix lhs(Ndeg);
    construct_d_matrix(lhs, nth, x);
    Mder[nth] = lhs * Vinv;
  }
}

//---------------------------------------------------------------------
/**
 * @brief construction of derivative matrices
 */
//---------------------------------------------------------------------

void vandermonde::construct_d_matrix(matrix& mat, int nth, double *x)
{
  for(  int i = 0; i < Ndeg; i++){
    for(int j = 0; j < Ndeg; j++){
      mat(i,j) = 0.0;
    }
  }

  for(  int i = 0; i < Ndeg; i++){
    for(int j = 0; j < Ndeg; j++){
      if(j < nth) continue;
      else{
	mat(i,j) = power(x[i], j - nth) * keisu(j,nth);
      }
    }
  }
}

//---------------------------------------------------------------------
/**
 * @brief construction of Linv
 */
//---------------------------------------------------------------------

void vandermonde::construct_Linv(matrix& ret, double *x)
{
  for(  int i = 0; i < Ndeg; i++){
    for(int j = 0; j < Ndeg; j++){
      ret(i,j) = 0.0;
    }
  }

  ret(0,0) = 1.0;

  for(  int i = 1; i <  Ndeg; i++){
    for(int j = 0; j <= i   ; j++){
      double tmp = 1.0;
      for(int k = 0; k <=i; k++){
	if (k==j) continue;
	tmp *= (x[j] - x[k]);
      }
      ret(i,j) = 1.0/tmp;
    }
  }
}

//---------------------------------------------------------------------
/**
 *
 */
//---------------------------------------------------------------------

void vandermonde::construct_Uinv(matrix& ret, double *x)
{
  for(  int i = 0; i < Ndeg; i++){
    for(int j = 0; j < Ndeg; j++){
      ret(i,j) = 0.0;
    }
  }

  for(int i = 0; i < Ndeg; i++){
    ret(i,i) = 1.0;
  }

  for(int j = 1; j < Ndeg; j++){
    ret(0,j) = -ret(0,j-1)*x[j-1];
  }
  for(  int i = 1;   i < Ndeg; i++){
    for(int j = i+1; j < Ndeg; j++){
      ret(i,j) = ret(i-1,j-1) - ret(i,j-1)*x[j-1];
    }
  }
}
