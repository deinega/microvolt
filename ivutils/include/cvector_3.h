# ifndef __CVECTOR_3_H
# define __CVECTOR_3_H

/// \en @file vector_3.h \brief complex N-dimensional vectors and some vector operations.

# include "vector_3.h"
# include <complex>

using namespace std;

typedef complex<vec_type> cvec_type;

///\en complex N-dimensional vector of type T with some useful operations
template <class T, size_t N> 
struct cVector_Nt {
  typedef Vector_Nt<T,N> vector_t;
  typedef complex<T> ctype;

  vector_t V[2];

  cVector_Nt(const vector_t &re=vector_t(), const vector_t &im=vector_t()) {
    V[0]=re;
    V[1]=im;
  }

  inline vector_t& re() { 
    return this->V[0];
  }

  inline vector_t& im() { 
    return this->V[1];
  }

  inline vector_t& real() { 
    return this->V[0];
  }

  inline vector_t& imag() { 
    return this->V[1];
  }

  inline const vector_t& real() const { 
    return this->V[0];
  }

  inline const vector_t& imag() const { 
    return this->V[1];
  }

  inline vector_t mod() const {
    vector_t v;
    for(int i=0;i<N;i++)
      v[i]=sqrt(V[0][i]*V[0][i]+V[1][i]*V[1][i]);
    return v;
  }

  inline cVector_Nt conj() const {
    return cVector_Nt(V[0],-V[1]);
  }

  inline ctype operator[](int i) const {
    return ctype(V[0][i],V[1][i]);
  }

  inline void set(int i, const ctype &val){
    V[0][i]=val.real();
    V[1][i]=val.imag();
  }

  inline int operator==(const cVector_Nt &cvect) const{
    return (V[0]==cvect.V[0] && V[1]==cvect.V[1]);
  }

  inline int operator!=(const cVector_Nt &cvect) const{
    return (!(*this==cvect));
  }

  inline cVector_Nt operator+(const cVector_Nt& cvect) const{
    return cVector_Nt(V[0]+cvect.V[0], V[1]+cvect.V[1]);
  }

  inline cVector_Nt operator-(const cVector_Nt &cvect) const {
    return cVector_Nt(V[0]-cvect.V[0], V[1]-cvect.V[1]);
  }

  inline cVector_Nt& operator+=(const cVector_Nt &cvect){
    V[0]+=cvect.V[0];
    V[1]+=cvect.V[1];
    return *this;
  }

  inline ctype operator*(const cVector_Nt &cvect) const {
    return ctype(V[0]*cvect.V[0]-V[1]*cvect.V[1], V[0]*cvect.V[1]+V[1]*cvect.V[0]);
  }

  inline cVector_Nt operator %(const cVector_Nt &cvect) const{
    return cVector_Nt(V[0]%cvect.V[0]-V[1]%cvect.V[1], V[0]%cvect.V[1]+V[1]%cvect.V[0]);
  }

  inline cVector_Nt operator %(const vector_t &vect) const{
    return cVector_Nt(V[0]%vect, V[1]%vect);
  }

  inline cVector_Nt operator*(const T &coeff) const {
    return cVector_Nt(V[0]*coeff, V[1]*coeff);
  }

  inline cVector_Nt operator*(const ctype &coeff) const {
    return cVector_Nt(V[0]*coeff.real() - V[1]*coeff.imag(), V[1]*coeff.real() + V[0]*coeff.imag());
  }

  inline cVector_Nt operator/(const ctype &coeff) const {
    return operator* (1./coeff);
  }

  inline cVector_Nt operator/(const T &coeff){
    return cVector_Nt(V[0]/coeff, V[1]/coeff);
  }

  inline cVector_Nt operator-(){
    return cVector_Nt(-V[0], -V[1]);
  }

  inline cVector_Nt& operator*=(const T &coeff){
    V[0]*=coeff;
    V[1]*=coeff;
    return *this;
  }

  inline cVector_Nt& operator*=(const ctype &coeff){
    *this=operator *(coeff);
    return *this;
  }

  T norm2() const{
    return V[0].norm2()+V[1].norm2();
  }

  T norm() const {
    return sqrt(norm2());
  }

  T normalize() {
    T norm=this->norm();
    if(norm>=VEC_ZERO){
      V[0]/=norm;
      V[1]/=norm;
    }
    return norm;
  }

  ctype cnorm2() const{
    ctype result=0;
    for (int i=0; i<3; i++){
      ctype val=ctype(V[0][i],V[1][i]);
      result+=val*val;
    }
    return result; 

  }

  ctype cnorm() const {
    return sqrt(cnorm2());
  }

  ctype cnormalize() {
    ctype cnorm=this->cnorm();
    T norm=this->norm();
    if(cnorm.real()>=VEC_ZERO || cnorm.imag()>=VEC_ZERO){
      for (int i=0; i<3; i++){
        ctype val=cvec_type(V[0][i],V[1][i]);
        val/=cnorm;
        V[0][i]=val.real();
        V[1][i]=val.imag();
      }
    }
    ctype ctest=this->cnorm();
    T test=this->norm();
    return cnorm;
  }
};

template<class T, size_t N>
cVector_Nt<T, N> operator*(const T &coeff,const cVector_Nt<T, N> &vec){
  return vec*coeff;
}

template<class T, class T2, size_t N>
cVector_Nt<T, N> operator*(const T2 &coeff,const cVector_Nt<T, N> &vec){
  return vec*coeff;
}


template<class T, size_t N>
cVector_Nt<T, N> conj(const cVector_Nt<T, N> &vec){
  return vec.conj();
}


typedef cVector_Nt<vec_type, 3> cVector_3;


template<class T, int N>
Vector_Nt<T, N> real_value(const cVector_Nt<T, N> &a){
  return a.real();
}


template<class T, int N>
struct real_t<const cVector_Nt<T, N> >{
  typedef Vector_Nt<T, N> data;
};


# endif
