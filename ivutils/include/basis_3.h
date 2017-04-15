/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 1995-2009        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   $Revision: 1.2 $
 *   $Date: 2013/01/18 23:52:18 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/basis_3.h,v 1.2 2013/01/18 23:52:18 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/basis_3.h,v $
$Revision: 1.2 $
$Author: lesha $
$Date: 2013/01/18 23:52:18 $
*/
/*s****************************************************************************
 * $Log: basis_3.h,v $
 * Revision 1.2  2013/01/18 23:52:18  lesha
 * eigen_solver is included
 *
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.12  2012/10/16 15:30:00  lesha
 * realization moved from cpp here
 *
 * Revision 1.11  2012/09/25 00:11:11  lesha
 * documentation
 *
 * Revision 1.10  2012/03/29 11:35:09  valuev
 * universal scatter_data() for emSourceWave
 *
*******************************************************************************/
/// \en @file vector_3.h \brief Basis and vector transformations
/// \ru @file vector_3.h \brief базис и преобразования векторов

#ifndef BASIS_3_H

#define BASIS_3_H

#include "vector_3.h"
#include "linsolv.h"

///\en basis
template<int N>
class Basis_N{
public:
  typedef Vector_Nt<vec_type,N> vector_t;
  vector_t vec[N];

  Basis_N(const vector_t &x1, const vector_t &x2, const vector_t &x3){
    if(N>0)vec[0]=x1;
    if(N>1)vec[1]=x2;
    if(N>2)vec[2]=x3;
  }
  Basis_N(vec_type stretch=1.){
    for(int i=0;i<N;i++)
      vec[i][i]=stretch;
  }
  Basis_N(const Basis_N &other){
    for(int i=0;i<N;i++)
      vec[i]=other.vec[i];
  }
  
  vec_type &operator()(int i, int j){
    return vec[j][i]; 
  }

  vec_type operator()(int i, int j) const{
    return vec[j][i]; 
  }

 ///\en transforms a vector from the orthogonal basis to this basis
 /// (contravariant form y=B*x)
  vector_t operator()(const vector_t &x) const{
    vector_t y;
    for(int i=0;i<N;i++){
      y[i]=0;
      for(int j=0;j<N;j++){
        y[i]+=vec[j][i]*x[j];
      }
    }
    return y;
  }

 ///\en transforms a vector from the orthogonal basis to this basis
 /// (covariant form y=transpose(B)*x )
  vector_t cov(const vector_t &x) const{
    vector_t y;
    for(int i=0;i<N;i++){
      y[i]=0;
      for(int j=0;j<N;j++){
        y[i]+=vec[i][j]*x[j];
      }
    }
    return y;
  }
  
 ///\en transforms a vector from this basis into the orthogonal basis
 /// (contravariant, x=B^(-1)*y )
  vector_t inv(const vector_t &x) const{
    vec_type matr[N*N];
    vector_t y, xx;
    for(int i=0;i<N;i++){
      xx[i]=x[i];
      for(int j=0;j<N;j++)
        matr[i*N+j]=vec[i][j];
    }
    linsys(matr,xx.v,N,y.v);
    //Vector_3 chk=(*this)(y);
    return y;
  }

  ///\en returns inverse basis
  Basis_N inv() const{
    Basis_N b;
    for(int i=0;i<N;i++){
      Vector_3 x;
      x[i]=1;
      b.vec[i]=inv(x);
    }
    return b;
  }

  inline vector_t& operator[](int i) {
    return vec[i];
  }

  inline vector_t& operator[](size_t i) {
    return vec[i];
  }

  inline vector_t operator[](int i) const {
    return vec[i];
  }

  inline vector_t operator[](size_t i) const {
    return vec[i];
  }
};

typedef Basis_N<3> Basis_3;

///\en basic class for vector transformation (default is no change)
struct VecTransform{
  virtual Vector_3 operator()(const Vector_3& arg) const{
    return arg;
  }
  virtual VecTransform *copy() const {
    return new VecTransform();
  }
  virtual ~VecTransform(){}
};

///\en vector constant shift
struct VecShift: public VecTransform{
  Vector_3 shift;
  VecShift(const Vector_3 &shift_=Vector_3()):shift(shift_){}
  virtual Vector_3 operator()(const Vector_3& arg) const{
    return arg+shift;
  }
  virtual VecTransform *copy() const {
    return new VecShift(shift);
  }
};

///\en vector homothety
struct VecHomothety: public VecTransform{
  vec_type a;
  VecHomothety(const vec_type a_):a(a_){}
  virtual Vector_3 operator()(const Vector_3& arg) const{
    return a*arg;
  }
  virtual VecTransform *copy() const {
    return new VecHomothety(a);
  }
};

#endif
