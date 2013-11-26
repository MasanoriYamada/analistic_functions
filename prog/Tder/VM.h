//---------------------------------------------------------------------
/**
 * @file
 * @brief   Vandermonde
 * @ingroup Group_Vandermonde
 * @author  Noriyoshi ISHII
 * @since   Thu May 30 16:20:12 JST 2013
 *
 * Copyright (C) 2013 by Noriyoshi Ishii, All Rights Reserved.
 *
 */
//---------------------------------------------------------------------

//---------------------------------------------------------------------
/**
 * @defgroup Group_Vandermonde
 */
//---------------------------------------------------------------------

#ifndef IS_INCLUDED_VM_H_130530
#define IS_INCLUDED_VM_H_130530

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

namespace VM {
  struct matrix;
  struct vandermonde;
};

//---------------------------------------------------------------------
/**
 * @brief matrix lib for van der mode matrix
 */
//---------------------------------------------------------------------

struct VM::matrix {
  int     dim;
  double *a;
  
  ~matrix() { if(a) delete[] a; }
  matrix() : dim(0), a(NULL) {}
  matrix(int n);
  matrix(const matrix& other);
  matrix& operator=(const matrix& other);

  const
  double& operator()(int i,int j) const { return a[range(i,j)]; }
  double& operator()(int i,int j)       { return a[range(i,j)]; }

  const
  double& operator[](int i) const { return a[range(i)]; }
  double& operator[](int i)       { return a[range(i)]; }

  matrix& operator+=(const matrix& other);
  matrix& operator-=(const matrix& other);
  matrix& operator*=(const matrix& other);

  matrix& operator*=(double x);
  matrix& operator/=(double x) { return operator*=(1.0/x); }

private:
  int range(int i) const;
  int range(int i,int j) const;
  void check_size(const matrix& other,const char *msg=NULL,...) const;
};

//---------------------------------------------------------------------
/**
 * @brief range checker (1)
 */
//---------------------------------------------------------------------

inline
int VM::matrix::range(int i) const
{
  if(!NDEBUG){
    if(i < 0 || dim <= i){
      fprintf(stderr, "ERROR(VM::matrix): out of range i=%d\n", i);
      exit(1);
    }
  }
  return i;
}

//---------------------------------------------------------------------
/**
 * @brief range checker (2)
 */
//---------------------------------------------------------------------

inline
int VM::matrix::range(int i,int j) const
{
  if(!NDEBUG){
    if(i < 0 || dim <= i || j < 0 || dim <= j){
      fprintf(stderr, "ERROR(VM:matrix): out of range i=%d, j=%d\n", i,j);
      exit(1);
    }
  }
  return i + dim*j;
}

//---------------------------------------------------------------------
/**
 * @brief size checker
 */
//---------------------------------------------------------------------

inline
void VM::matrix::check_size(const matrix& other,const char *msg,...) const
{
  if(!NDEBUG){
    if(other.dim != dim){
      fprintf(stderr,
	      "ERROR(VM::matrix): inconsistent dim %d != %d\n", dim, other.dim);
      if(msg){
	va_list args;
	va_start(args, msg);
	vfprintf(stderr, msg, args);
	va_end(args);
	fprintf(stderr, "\n");
      }
      exit(1);
    }
  }
}

//---------------------------------------------------------------------
// helper functions of VM::matrix
//---------------------------------------------------------------------

inline VM::matrix operator+(VM::matrix lhs,const VM::matrix& rhs){ return lhs += rhs; }
inline VM::matrix operator-(VM::matrix lhs,const VM::matrix& rhs){ return lhs -= rhs; }
inline VM::matrix operator*(VM::matrix lhs,const VM::matrix& rhs){ return lhs *= rhs; }

inline VM::matrix operator*(VM::matrix mat, double x)       { return mat *= x; }
inline VM::matrix operator*(double x,       VM::matrix mat) { return mat *= x; }
inline VM::matrix operator/(VM::matrix mat, double x)       { return mat /= x; }

namespace VM {
  template <class T>
  inline void fmadd(T ret[], matrix& mat, T vec[], bool is_reset=true)
  {
    int Ndeg = mat.dim;
    for(  int i = 0; i < Ndeg; i++){
      for(int j = 0; j < Ndeg; j++){
	if(true == is_reset && j == 0){
	  ret[i] = mat(i,j)*vec[j];
	}
	else{
	  ret[i] += mat(i,j)*vec[j];
	}
      }
    }
  }
}

//---------------------------------------------------------------------
/**
 * @brief vandermonde matrix
 */
//---------------------------------------------------------------------

struct VM::vandermonde {
  int     Ndeg;
  matrix  Vinv;
  matrix *Mder;

  ~vandermonde();
  vandermonde();
  vandermonde(int n);
  vandermonde(const vandermonde& other);

  vandermonde& operator=(const vandermonde& other);
  
private:
  void construct_vinv(double *x);
  void construct_mder(double *x);

  void construct_d_matrix(matrix& mat,int nth,double *x);

  void construct_Linv(matrix& Linv, double *x);
  void construct_Uinv(matrix& Uinv, double *x);
};

#endif
