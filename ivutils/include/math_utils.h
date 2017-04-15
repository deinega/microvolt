#ifndef _MATH_UTILS_H
#define _MATH_UTILS_H

/// \en @file math_utils.h \brief Some useful math functions

#include <functional>
#include <complex>
#include <limits>
#include "refobj.h"

/// tests if x is NaN
template<class T>
int test_NaN(T x){
  return (!(x>0 || x<0 || x==0));
}

///\en round given value up to the near number of the type 2^N
///\ru округляет число вверх до ближайшего числа вида 2^N
template<class T>
T ceil2(const T val){
  int iexp;
  frexp(val,&iexp);
  T res=ldexp(.5, iexp+1);
  return res;
}

///\en compare x and y within some accuracy which depends on parameter val.
///\ru сравнивает x и y с точностью до некоторого малого числа (зависящего от параметра val)
///    по умолчанию y=0
template<class T>
bool accomp(const T x, const T y, const T val=T(1), int mult_epsilon=1024) {
  T mult=val==T(1) ? T(2) : ceil2(val);
  T err=mult*T(mult_epsilon)*numeric_limits<T>::epsilon();
  return fabs(x-y)<=err;
}



///\en return floor(x) or ceil(x) [if x is equal to ceil(x) within some accuracy which depends on parameter val]
template<class T>
T acfloor(const T x, T val=T(1), int mult_epsilon=1024) {
  if(accomp(ceil(x),x,val,mult_epsilon))
    return ceil(x);
  else
   return floor(x);
}

template<class T>
bool acless(const T x, const T y, const T val=T(1), int mult_epsilon=1024) {
  if(x<=y || accomp(x,y,val,mult_epsilon))return true;
  return false;
}

///\en Compare predicate based on \ref accomp to use in STL
template<class T> 
class accomp_pred{
  int me; 
public:
  accomp_pred(int mult_epsilon=1024):me(mult_epsilon){}
  ///\en \return true if a<b
  bool operator()(const T& a, const T& b) const {
    if(accomp(a,b,T(1),me)) // equal with given accuracy
      return false;
    return a<b;
  }
};


///\en calculates ratio between x and y.
/// if this ratio is close to some integer number with some accuracy,
/// then returns this integer number.
///\ru делит x на y. Если их частное меньше некоторого малого числа (зависящего от величины частного), то оно округляется.
template<class T>
T acdiv(T x, T y, int mult_epsilon=1024) {
  T fr=x/y;
  T ifr=floor(fr+T(.5));
  int mult = fabs(ifr)<=T(mult_epsilon) ? mult_epsilon : int(ceil2(fabs(ifr)));
  T err=mult*numeric_limits<T>::epsilon();
  if (fabs(fr-ifr)<=err)
    return ifr;
  else
    return fr;
}

///\en Logarithmically reduces the argument:
/// log_reduce(x) = log(x), x>e
/// log_reduce(x) -> -log(-x), x<e
/// linear between [-e,e], first derivative continuous
template <class T>
T log_reduce(const T& x){
  return x<-M_E ? -log(-x) : x<M_E ? x/M_E : log(x) ;
  //return x<0 ? -x*x : x*x ;
  //return x<0 ? -sqrt(-x) : x==0 ? 0 : sqrt(x) ;
  //return x;
}

template<class T>
T real_value(const T &a){
  return a;
}

template<class T>
T real_value(const complex<T> &a){
  return a.real();
}

template<class T>
struct real_t{
  typedef T data;
};

template<class T>
struct real_t<complex<T> >{
  typedef T data;
};



///\en returns -1 if t1<t2; 1 if t2>t1; 0 if t1==t2
/// can be used in implementation of operators < or > for some classes
template<class T>
int step_fun(const T &t1, const T &t2){
  if(t1<t2)return -1;
  else if(t2<t1)return 1;
  return 0;
}

template<class T>
int step_fun(const complex<T> &t1, const complex<T> &t2){
  if(t1.real()<t2.real())return -1;
  else if(t1.real()>t2.real())return 1;
  else if(t1.imag()<t2.imag())return -1;
  else if(t1.imag()>t2.imag())return 1;
  return 0;
}

/*
template<template<class T> class containter_r>
int step_fun(const containter_r<class T> &t1, const containter_r<class T> &t2){
  if(t1.size()<t2.size())return -1;
  if(t1.size()>t2.size())return 1;
  
  for(typename containter_r::const_iterator it1=t1.begin(),it2=t2.begin(),e=t1.end();it1!=e;++it1,++it2){
    int res=step_fun(*it1,*it2);
    if(res)return res;
  }
  return 0;
}
*/
/*
template<class T>
int step_fun(const vector<T> &t1, const vector<T> &t2){
  if(t1.size()<t2.size())return -1;
  if(t1.size()>t2.size())return 1;
  
  for(typename vector<T>::const_iterator it1=t1.begin(),it2=t2.begin(),e=t1.end();it1!=e;++it1,++it2){
    int res=step_fun(*it1,*it2);
    if(res)return res;
  }
  return 0;
}
*/

template<class arg_t, class result_t>
class virt_unary_function: public unary_function<arg_t,result_t>{
public:
  virtual result_t operator()(arg_t x)=0;
  virtual ~virt_unary_function(){}
};

template<class arg1_t, class arg2_t, class result_t>
class virt_binary_function: public binary_function<arg1_t,arg2_t,result_t>{
public:
  virtual result_t operator()(arg1_t x1,arg2_t x2)=0;
  virtual ~virt_binary_function(){}
};

template<class arg_t, class result_t>
class delta_function: public virt_unary_function<arg_t,result_t>{
  arg_t x0;
public:
  delta_function(arg_t x0_):x0(x0_){}
  virtual result_t operator()(arg_t x){
    return accomp(x,x0) ? 1 : 0;
  }
};

///\en used to construct max_function and min_function
template<class arg_t, class result_t, class comp_t=greater<result_t> >
class comp_function: public virt_unary_function<arg_t,result_t>{
protected:
  typedef virt_unary_function<arg_t,result_t> fun_t;
  refvector<fun_t> funs;
public:
  void add_function(fun_t *fun){
    funs.push_back(fun);
  }
  virtual result_t operator()(arg_t x){
    int sign = comp_t()(1,-1) ? 1 : -1;
    result_t res=-sign*numeric_limits<result_t>::max();
    for(size_t i=0;i<funs.size();i++){
      result_t y=funs[i]->operator()(x);
      comp_t()(y,res);
      if(sign>0 ? y>res : y<res)
        res=y;
    }
    return res;
  }
};

///\en maximal value of given functions
template<class arg_t, class result_t>
class max_function: public comp_function<arg_t,result_t,greater<result_t> >{};

///\en minimal value of given functions
template<class arg_t, class result_t>
class min_function: public comp_function<arg_t,result_t,less<result_t> >{};

template<class arg_t, class result_t>
class sum_function: public virt_unary_function<arg_t,result_t>{
protected:
  typedef virt_unary_function<arg_t,result_t> fun_t;
  refvector<fun_t> funs;
public:
  void add_function(fun_t *fun){
    funs.push_back(fun);
  }
  virtual result_t operator()(arg_t x){
    result_t res=0;
    for(size_t i=0;i<funs.size();i++)
      res+=funs[i]->operator()(x);
    return res;
  }
};


/*
///\en Integrates a function according to Simpson's formula with overlapping segments:\n
///    Press, William H., Brian P. Flannery, William T. Vetterling, and Saul A. Teukolsky (1989). Numerical Recipes in Pascal: 
///    The Art of Scientific Computing. Cambridge University Press. ISBN 0-521-37516-9. \n
///    \a func should have operator(const arg_t) const \n
///    
template<class func_t, class value_t, class arg_t>
void simpson_integrate(value_t & result, const func_t &func, arg_t x0, arg_t x1, arg_t dx){ 
  arg_t L= x1-x0;
  int sign = 1;
  if(L<0){
    L=-L;
    std::swap(x0,x1);
    sign = -1;
  }
  if(dx<0)
    dx = -dx;
  if(L<=dx)
    return L*(func(x0)+func(x1))/2.;
  int n = L/dx + 1;
  if(n<8)
    n = 8;
  dx = L/(n-1);
  result = func(x0)*(17./48.);  //0
  x0 += dx;
  result+ = func(x0)*(59./48.);  // 1
  x0 += dx;
  result+ = func(x0)*(43./48.);  // 2
  x0 += dx;
  result+ = func(x0)*(49./48);  // 3
  for(int i =4; i< n-4; i++){
    result += func(x0); // 4 -> (n-4)-1
  }
  result = func(x0)*(49./48.);  // n -4
  x0 += dx;
  result+ = func(x0)*(43./48);  // n-3
  x0 += dx;
  result+ = func(x0)*(59./48);  // n-2
  x0 += dx;
  result+ = func(x0)*(17./48);  // n-1
  return  result*dx*sign;
}
*/
#endif
