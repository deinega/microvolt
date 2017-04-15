/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   $Revision: 1.7 $
 *   $Date: 2013/04/26 13:59:02 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/vector_3.h,v 1.7 2013/04/26 13:59:02 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/vector_3.h,v $
$Revision: 1.7 $
$Author: valuev $
$Date: 2013/04/26 13:59:02 $
*/
/*s****************************************************************************
 * $Log: vector_3.h,v $
 * Revision 1.7  2013/04/26 13:59:02  valuev
 * fabs -> abs
 *
 * Revision 1.6  2013/02/22 11:01:26  valuev
 * thin surface interface
 *
 * Revision 1.5  2013/02/02 19:56:29  valuev
 * fixed new vector product
 *
 * Revision 1.4  2013/02/01 18:01:12  valuev
 * added vector volume (needs testing)
 *
 * Revision 1.3  2013/01/29 16:41:32  valuev
 * added faceted surface and rotation quaternions
 *
 * Revision 1.2  2012/12/06 23:23:25  lesha
 * randdir is moved to region
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.64  2012/11/21 10:00:12  lesha
 * rpcell is modidifed (for transport)
 *
 * Revision 1.63  2012/09/27 23:44:17  lesha
 * bug is fixed
 *
 * Revision 1.62  2012/09/24 19:25:30  lesha
 * documentation
 *
 * Revision 1.61  2012/07/26 21:33:37  valuev
 * added rt_angular experiment
 *
 * Revision 1.60  2012/07/25 16:55:35  belousov
 * modified randdir(int flag2d=0), flag2d specifies randomization type (0 - 3d, 1 - 2d in xy plane)
 *
 * Revision 1.59  2012/06/20 21:41:43  lesha
 * *** empty log message ***
 *
 * Revision 1.58  2012/06/16 01:06:05  valuev
 * sync with GridMD project
 *
 * Revision 1.57  2012/05/30 22:22:58  lesha
 * *** empty log message ***
 *
 * Revision 1.56  2012/04/24 03:12:43  lesha
 * step_fun is added
 *
 * Revision 1.55  2012/04/13 12:41:24  valuev
 * moved stdlib.h back (NULL is defined there under llinux)
 *
 * Revision 1.54  2012/04/13 08:13:21  valuev
 * reverted to M_PI
 *
 * Revision 1.53  2012/04/12 20:11:24  lesha
 * documentation
 *
 * Revision 1.52  2012/02/29 23:41:06  lesha
 * compilation with vs10
 *
*******************************************************************************/

#ifndef VECTOR_3_H
#define VECTOR_3_H

/// \en @file vector_3.h \brief N-dimensional vectors and some useful functions for work with them.
/// \ru @file vector_3.h \brief N-мерные вектора и работа с ними.

#include <stdlib.h> // нужна для NULL,  RAND_MAX в randdir
#include <cmath>
#include <limits>

using namespace std;

#include "math_utils.h"

#ifndef SINGLE_PRECISION
typedef double vec_type;
#else
typedef float vec_type;
#endif


///\en "infinitely large" number
///\ru "бесконечно большое" число
#define VEC_INFTY numeric_limits<vec_type>::max()

///\en "infinitely large" integer number
#define INT_INFTY numeric_limits<int>::max()

///\en how bigger is "ininitely small" number then numeric_limits<...>::epsilon().
///    We use this parameter since numeric_limits<...>::epsilon() is too small.
///\ru во сколько раз "бесконечно малое" число больше numeric_limits<...>::epsilon().
///    Этот множитель необходим потому, что numeric_limits<...>::epsilon() слишком мало.
#define MULT_EPSILON 1024

///\en "infinitely small" number
///\ru "бесконечно малое" число
#define VEC_ZERO MULT_EPSILON*numeric_limits<vec_type>::epsilon()

template <class T, int N> 
struct Vector_Nt;


///\en Traits class for vector product
template<class T, int N>
struct vec_prod_traits {
  typedef void type;
};

template<class T>
struct vec_prod_traits<T, 2> {
  typedef T type;
};

template<class T>
struct vec_prod_traits<T, 3> {
  typedef Vector_Nt<T, 3> type;
};

///\en N-dimensional vector of the type T and some useful operations
///\ru N-мерный вектор типа T с некоторыми полезными операциями
template <class T, int N> 
struct Vector_Nt{
  typedef T value_type;
  static const int dimension=N;
  typedef typename vec_prod_traits<T,N>::type vec_product_t; ///<\en vector product type
  //typedef Vector_Nt<T, N> vec_product_t; ///<\en vector product type

  T v[N];

  ///\en makes all components equal to a
  ///\ru присваивает всем компонентам значение a
  Vector_Nt(const T &a=0){
    for(int i=0;i<N;i++) v[i]=a;
  }

  inline Vector_Nt& operator=(const T &a){
    for(int i=0;i<N;i++) v[i]=a;
    return *this;
  }

  explicit Vector_Nt(const T &a1, const T &a2){
    if(N>0)v[0]=a1;
    if(N>1)v[1]=a2;
    for(int i=2;i<N;i++)v[i]=0;
  }

  explicit Vector_Nt(const T &a1, const T &a2, const T& a3){
    if(N>0)v[0]=a1;
    if(N>1)v[1]=a2;
    if(N>2)v[2]=a3;
    for(int i=3;i<N;i++)v[i]=0;
  }

  ///\en construct from input iterator (or array)
  template <class A>
  Vector_Nt(const A *beg){
    for(int i=0;i<N;i++,++beg)
      v[i]=*beg;
  }

  ///\en construct from another vector
  template <class A>
  Vector_Nt(const Vector_Nt<A,N> &vect){
    for(int i=0;i<N;i++) v[i]=vect[i];
  }

  template <class A>
  inline Vector_Nt& operator=(const Vector_Nt<A,N> &vect){
    for(int i=0;i<N;i++) v[i]=vect[i];
    return *this;
  }

  ///\en Copies vector to iterator
  ///\ru Копирует содержимое вектора в итератор
  template <class A>
  void copy_to(A *beg) const{  
    for(int i=0;i<N;i++,++beg)
      *beg=v[i];
  }

  ///\en obtains element value
  ///\ru получение элемента 
  inline T& operator[](int i) const {return (T&)v[i];}

  inline T& operator[](size_t i) const {return (T&)v[i];}

  ///\en comparison. If the difference is less then VEC_ZERO (MULT_EPSILON*numeric_limits<T>::epsilon()) then components are assumed to be equal
  ///\ru сравнение. При отличии меньше чем на VEC_ZERO (MULT_EPSILON*numeric_limits<T>::epsilon()) компоненты считаются одинаковыми
  inline bool operator==(const Vector_Nt &vect) const{
    for(int i=0;i<N;i++)
      if(fabs(v[i]-vect[i])>MULT_EPSILON*numeric_limits<T>::epsilon())return false;
    return true;
  }

  inline bool operator!=(const Vector_Nt &vect) const{
    return (!(this->operator==(vect)));
  }

  inline Vector_Nt operator+(const Vector_Nt& vect) const{
    Vector_Nt result;
    for (int i=0; i<N ;i++)
      result.v[i]=v[i]+vect.v[i];
    return result;
  }

  inline Vector_Nt operator-(const Vector_Nt &vect) const {
    Vector_Nt result;
    for (int i=0; i<N ;i++)
      result.v[i]=v[i]-vect.v[i];
    return result;
  }
 
  ///\en Scalar product
  ///\ru Скалярное произведение векторов
  inline T operator*(const Vector_Nt& vect) const{
    T result=0;
    for (int i=0; i<N; i++)
      result+=v[i]*vect.v[i];
    return result;
  }

  ///\en Multiplies on coefficient
  ///\ru Покомпонентное умножение на коэффициент
  inline Vector_Nt operator*(const T &coeff) const{
    Vector_Nt result;
    for (int i=0; i<N; i++)
      result[i]=coeff*v[i];
    return result;
  }

  ///\en Multiplies on vector by components
  ///\ru Покомпонентное умножение на вектор
  template <class Vec>
  inline Vector_Nt scale(const Vec &s) const{
    Vector_Nt result;
    for(int i=0;i<N;i++)
      result[i]=s[i]*v[i];
    return result;
  }


  ///\en Divides on coefficient
  ///\ru Покомпонентное деление на коэффициент
  inline Vector_Nt operator/(const T &coeff) const {
    Vector_Nt result;
    for (int i=0; i<N; i++)
      result[i]=v[i]/coeff;
    return result;
  }

  ///\en Vector product (N=3 only)
  ///\ru Векторное произведение
  //inline vec_product_t operator%(const Vector_Nt &vect) const; //reserved for N specializations


  
  ///\en Multiplies on -1
  ///\ru Умножение вектора на -1
  inline Vector_Nt operator-() const {
    Vector_Nt result;
    for(int i=0;i<N;i++)
      result.v[i]=-v[i];
    return result;
  }

  ///\ru Сложение с присваиванием
  inline Vector_Nt& operator+=(const Vector_Nt &vect){
    for (int i=0; i<N; i++)
      v[i]+=vect.v[i];
    return *this;
  }

  ///\ru Вычитание с присваиванием
  inline Vector_Nt& operator-=(const Vector_Nt &vect){
    for (int i=0; i<N; i++)
      v[i]-=vect.v[i];
    return *this;
  }

  ///\ru Умножение на коэффициент с присваиванием
  inline Vector_Nt& operator*=(const T &coeff){
    for (int i=0; i<N; i++)
      v[i]*=coeff;
    return *this;
  }

  ///\ru Деление на скаляр с присваиванием
  inline Vector_Nt& operator/=(const T &coeff){
    for (int i=0; i<N; i++)
      v[i]/=coeff;
    return *this;
  }

  ///\en Norm squared
  ///\ru Квадрат нормы вектора
  T norm2() const {
    T result=0;
    for (int i=0; i<N; i++)
      result+=v[i]*v[i];
    return result; 
  }

  ///\en Norm
  ///\ru Норма вектора
  T norm() const {
    return sqrt(norm2());
  }

  ///\en Returns norm and normalizes vector on newnorm
  ///\ru Возвращает норму и нормализует вектор на newnorm
  T normalize(const T &newnorm=1.){
    T norm=this->norm();
    if(norm>=MULT_EPSILON*numeric_limits<T>::epsilon()){
      T c=newnorm/norm;
      for (int i=0; i<N; i++)
        v[i]*=c;
    }
    return norm;
  }

  ///\en \return a vector which is normal to current, or curren if the norm iz less than \ref VEC_ZERO
  Vector_Nt normal() const {
    T n = norm();
    return n<= MULT_EPSILON*numeric_limits<T>::epsilon() ? (*this): (*this)/n;    
  }

  ///\en Normalizes a vector that may have infinite components
  T inormalize(T newnorm=1., T infty=VEC_INFTY){
    T result=0;
    for(int i=0;i<N;i++){
      if(fabs(v[i])>=infty){
        if(result>=0){
          for(int j=0;j<i;j++)
            v[j]=0.;
          result=0.;
        }
        result-=1;
        v[i]=v[i]>0 ? 1 : -1;
      }
      else if(result>=0) //still summing the normal components
        result+=v[i]*v[i];
      else
        v[i]=0.;
    }
    if(fabs(result)<MULT_EPSILON*numeric_limits<T>::epsilon())
      return 0.;
    newnorm/=sqrt(fabs(result));
    for (int i=0; i<N; i++)
      v[i]*=newnorm;
    return result<0 ? infty : result;
  }

  ///\en return a vector projection on a given axis
  Vector_Nt prj(int i) const {
    Vector_Nt res;
    res[i]=v[i];
    return res;
  }

  ///\en return a vector projection on a given vector (k*v)*k
  Vector_Nt prj(const Vector_Nt &k) const {
    return k*(k*(*this));
  }

  ///\en return a vector of length l parallel to this
  Vector_Nt unit(T l=1) const {
    Vector_Nt res(*this);
    res.normalize(l);
    return res;
  }

  ///\en nearest image distance within rectangular cell (FOR DISTANCE MEASUREMENTS)
  ///    assumes that each coordinate absolute value is in the range [0,cell[i]) 
  ///    returned vector is in the range [-cell[i]/2,cell[i]/2)
  ///    flags indicate the periodicity in specific directions: 0x1 for X, 0x2 for Y, 0x4 for Z
  ///\ru Ближайший образ в прямоугольной ячейке
  ///    Считаем, что все пространство разделено на ячейки - параллелепипеды с ребрами, параллельными
  ///    осям координат и диагональю, заданной вектором rcell, причем начало координат является 
  ///    центром одной из ячеек. Если *this находится центральной ячейке, возвращается копия *this.\n
  ///    Иначе, если *this находится в кубе 3*3 ячейки с центром в начале координат, то создает образ 
  ///    *this в центральной ячейке.\n
  ///    Иначе, возвращает неопределенное значение.  
  Vector_Nt rcell1(const Vector_Nt &cell,int flags=0xffff) const{
    Vector_Nt ret(*this);
    for(int i=0;i<N;i++){
      if(flags&(0x1<<i)){
        if(v[i]>cell[i]/2){
          ret[i]-=cell[i];
        }
        else if(v[i]<-cell[i]/2){
          ret[i]+=cell[i];
        }
      }
    }
    return ret;
  }

  ///\en reduction to elementary cell [0, cell[i]) (FOR REDUCTION TO ELEMENTARY CELL)
  ///    flags indicate the periodicity in specific directions: 0x1 for X, 0x2 for Y, 0x4 for Z
  ///\ru Почти то же, что и rcell1, но без ограничения на положение *this и с другой системой ячеек.
  ///    В начале координат находится не центр ячейки, а ее угол. Может работать медленнее из-за наличия 
  ///    операции деления по модулю с плавающей точкой
  Vector_Nt rcell(const Vector_Nt &cell, int flags=0xffff) const {
    Vector_Nt ret(*this);
    for (int i=0, flag=1; i<N; i++, flag<<=1) {
      if(flags&flag){
        ret.v[i]=fmod(v[i],cell.v[i]);
        if(ret[i]<0)ret.v[i]+=cell.v[i];
      }
    }
    return ret;
  }

  ///\en the same as rcell, but start point of zero cell is p1
  ///\ru то же самое, что и rcell, только начало нулевой ячейки имеет координаты p1
  Vector_Nt rpcell(const Vector_Nt &p1, const Vector_Nt &cell, int flags=0xfff) const {
    Vector_Nt ret(*this);
    for (int i=0, flag=1; i<N; i++, flag<<=1) {
      if(flags&flag){
//        if (ret[i]<p1[i] || ret[i]>p1[i]+cell[i]) {
        if (!acless(p1[i],ret[i]) || !acless(ret[i],p1[i]+cell[i])) {
          ret[i]=fmod(v[i]-p1[i],cell[i])+p1[i];
          if (ret[i]<p1[i]) ret[i]+=cell[i];
        }
      }
    }
    return ret;
  }
  
  ///\en returns maximal vector component and its index
  ///\ru Возвращает максимальную компоненту вектора и ее индекс в ind 
  T maxcoord(int *ind=NULL) const {
    int im=0;
    T vv=v[0];
    for (int i=1; i<N; i++) {
      if(v[i]>vv){
        im=i;
        vv=v[i];
      }
    }
    if(ind)*ind=im;
    return vv;
  }

  ///\en returns minimal vector component and its index
  ///\ru Возвращает минимальную компоненту вектора и ее индекс в ind 
  T mincoord(int *ind=NULL) const {
    int im=0;
    T vv=v[0];
    for (int i=1; i<N; i++) {
      if(v[i]<vv){
        im=i;
        vv=v[i];
      }
    }
    if(ind)*ind=im;
    return vv;
  }

  ///\en returns the coord having maximal absolute value
  T maxabscoord(int *ind=NULL) const {
    int im=0;
    T vv=abs(v[0]);
    for (int i=1; i<N; i++) {
      T vi = abs(v[i]);
      if(vi>vv){
        im=i;
        vv=vi;
      }
    }
    if(ind)*ind=im;
    return v[im];
  }

  ///\en returns the coord having minimal absolute value
  T minabscoord(int *ind=NULL) const {
    int im=0;
    T vv=abs(v[0]);
    for (int i=1; i<N; i++) {
      if(abs(v[i])<vv){
        im=i;
        vv=abs(v[i]);
      }
    }
    if(ind)*ind=im;
    return v[im];
  }

  ///\en returns true if the vector has infinite components
  bool infinite(T infty=VEC_INFTY) const {
    for(int i=0;i<N;i++){
      if(fabs(v[i])>=infty)
        return true;
    }
    return false;
  }
};

template<class T, int N>
Vector_Nt<T, N> operator*(const T &coeff,const Vector_Nt<T, N> &vec){
  return vec*coeff;
}

template<class T, class T2, int N>
Vector_Nt<T, N> operator*(const T2 &coeff,const Vector_Nt<T, N> &vec){
  return vec*coeff;
}

///\en Vector product for N=3
///\ru Векторное произведение
template<class T>
Vector_Nt<T, 3> operator%(const Vector_Nt<T, 3> &x, const Vector_Nt<T, 3> &y){
  return Vector_Nt<T,3>(x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]);
}

///\en Vector product for N=2
///\ru Векторное произведение
template<class T>
vec_type operator%(const Vector_Nt<T, 2> &x, const Vector_Nt<T, 2> &y){
  return x[0]*y[1]-x[1]*y[0];
}



// old Vector_3 compatibility typedefs and functions
typedef Vector_Nt<int,2> iVector_2;
typedef Vector_Nt<int,3> iVector_3;
typedef Vector_Nt<vec_type, 2> Vector_2;
typedef Vector_Nt<vec_type, 3> Vector_3;
typedef Vector_3 *Vector_3P; 
typedef Vector_2 *Vector_2P;


template <int N> 
class  Vector_N: public Vector_Nt<vec_type, N>{
};

///\en returns polygon area based on vectors vect1 and vect2
///\ru возвращает площадь параллелограмма, построенного на векторах vect1 и vect2
inline vec_type vec_area(const Vector_2 &vect1, const Vector_2 &vect2) {
  return fabs(vect1[0]*vect2[1]-vect1[1]*vect2[0]); // ??? MUST BE WITHOUT fabs!
};

inline vec_type vec_area(const Vector_3 &vect1, const Vector_3 &vect2) {
  return (vect1%vect2).norm();
};



///\en finds the maximum distance between vector pairs
///\ru Находит максимальное расстояние между векторами va1[i], va2[i], i=1..n
///    \param va1 - массив Vector_3[n]
///    \param n - длина массивов va1 и va2
vec_type dist_max(Vector_3 *va1,Vector_3 *va2,int n);

///\en finds average distance between vector pairs
///\ru Находит среднее расстояние между векторами va1[i], va2[i], i=1..n
vec_type dist_av(Vector_3 *va1,Vector_3 *va2,int n);

///\en finds the average difference norm between two vector sets of the same length
///    optionally gives the indices for maximal and minimal difference
///    va2 can be NULL, then the norm of va1 is used
///\ru Находит среднее расстояние между va1[i] и va2[i], а также, по желанию, индексы, на которых достигается min и max расстояние
vec_type diff_av(Vector_3 *va1,Vector_3 *va2,int n, int *minind=0, int *maxind=0);

///\en finds suitable perpendicular to a vector
///\ru Находит перпендикуляр к вектору vAB
Vector_3 FindPerp(const Vector_3 &vAB);

///\en Returns the average (center) vector of the vector array
///    and cooordinates of a minimal box in which it is contained
Vector_3 GetScope(const Vector_3 *varr,long n,Vector_3* box_min,Vector_3* box_max);

///\en the same with long index array
Vector_3 GetIScope(const Vector_3 *varr,long *indarr,long n,Vector_3* box_min,Vector_3* box_max);

///\en the same with int index array
Vector_3 GetIScopei(const Vector_3 *varr,int *indarr,int n,Vector_3* box_min,Vector_3* box_max);

// neue Funktionen

///\en clears vector array with optional integer index
///\ru Очистка массива векторов, с поддержкой индексирования 
///    В данном Vector_3 vec[] обнуляет n координат. Если ind==NULL, то 
///    очищает первые n элементов. Если ind!=NULL, то для i=0..n-1
///    очищает vec[ind[i]]
///    См. \ref indexed_calculations.
void clear_vecarri(int n,Vector_3 *vec, int *ind=0);

///\en reflects the vector ini+dir*t+0.5*force*t^2 to be inside a box limited by 0 and box sizes
///    changes dir according to the final state
///    fills crossed dir with bit flags corresponding directions along which the walls were crossed
Vector_3 Reflect(Vector_3& ini, double t,Vector_3 &dir, double *box, int flag=0x7, const Vector_3 &force=Vector_3()); 

///\en Calculates extent of the vector container.
///    \return the center of the vector set, optionally
///    (if arguments are not NULL) fills the bounding box in \a box_min, \a box_max.
template<class vec_inp_it>
Vector_3 get_extent(vec_inp_it beg,vec_inp_it end, Vector_3* box_min=NULL,Vector_3* box_max=NULL){
  if(beg==end)
    return Vector_3();
  Vector_3 center(*beg++);
  Vector_3 cube1(center), cube2(center);
  size_t n=1;
  for(;beg!=end;++beg){
    Vector_3 vec=*beg;
    center+=vec;
    for(size_t j=0;j<3;j++){
      if(cube1[j]>vec[j])
        cube1[j]=vec[j];
     if(cube2[j]<vec[j])
        cube2[j]=vec[j];
    }
    n++;
  }
  if(box_min)
    *box_min=cube1;
  if(box_max)
    *box_max=cube2;
  return center/n;
}


///\en Performs a step of the Stabilized Gramm-Schidt orthonormalization algorithm. 
///    Given a set of orthogonal unit vectors defined by [orth_beg, orth_end)
///    orthonormalizes inp_vec with respect to this set (removes all projections to the set).
///    The result is recorded to out_vec, which is normalized.
///    \return the norm of out_vec before normalization. If return value is <=\ref VEC_ZERO,
///    then inp_vec may be considered as linearly dependent of the supplied set.
template<class inp_it, class vector_tt>
typename vector_tt::value_type gramm_schmidt_project(inp_it orth_set_beg, inp_it orth_set_end, const vector_tt &inp_vec, vector_tt &out_vec){
  out_vec = inp_vec;
  for(;orth_set_beg!=orth_set_end; ++orth_set_beg)
    out_vec -= (*orth_set_beg)*out_vec; 
  return out_vec.normalize();
}

template<class value_t, int N>
int step_fun(const Vector_Nt<value_t, N> &t1, const Vector_Nt<value_t, N> &t2){
  for(int i=0;i<N;i++){
    int res=step_fun(t1[i],t2[i]);
    if(res)return res;
  }
  return 0;
}


inline Vector_3 FindPerp(const Vector_3 &vAB){
  Vector_3 res;
  // finding max. comp
  vec_type vmax=fabs(vAB[0]);
  vec_type vmin=vmax;
  int maxind=0, minind=0,i;
  for(i=1;i<3;i++){
    vec_type v=fabs(vAB[i]);
    if(vmax<v){
      vmax=v;
      maxind=i;
    }
    if(vmin>v){
      vmin=v;
      minind=i;
    }
  }
  vec_type s=0;
  if(minind==maxind){// zero vector
    res=Vector_3(1,1,1);
  }
  else{
    res[minind]=1.;
    res[maxind]=-vAB[minind]/vAB[maxind];
  }
  res.normalize();
  return res;
}


# endif // VECTOR_3_H

