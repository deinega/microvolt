#ifndef __LINSOLV_H
#define __LINSOLV_H

/// @file linsolv.h \brief Interface for solving systems of linear equations using LAPACK and PARDISO

#include <vector>
#include <utility>
#include <stdio.h>

using namespace std;

#include "math_utils.h"

#define MATR_INFTY 1e15

#ifdef USE_MKL // MKL includes realization of LAPACK and PARDISO

#include "mkl.h"

#ifndef USE_LAPACK
#define USE_LAPACK
#endif

#ifndef USE_PARDISO
#define USE_PARDISO
#endif

#else // LAPACK or/and PARDISO is linked from some lib file (no MKL)
#ifdef USE_LAPACK
extern "C" void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
extern "C" void zgetrf_(const int *M, const int *N, complex<double> *A, const int *lda, int *ipiv, int *info);
extern "C" void zgetri_(const int *N,complex<double> *A,const int *lda,const int *ipiv,const complex<double> *work,const int *lwork,int *info);
extern "C" void zheev_(const char *jobz,const char *uplo,const int *N,complex<double> *matr,
                       const int *lda,double *x,const complex<double> *work,const int *lwork, double *rwork,int *info);
#endif
#ifdef USE_PARDISO
extern "C" void PARDISO(void **, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *);
#endif
#endif

#ifdef USE_LAPACK
// calls dgesv_ (T is double) or sgesv_ (T is float) LAPACK routine
// parameters are the same as corresponding LAPACK routines
// http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve.html#ga7815e25af1fb6f8ee8fd5fd0fd1dc245
// returns 1 if everything ok, or -1 otherwise
template<class T> void tgesv_(int *N, int *nrhs, T *A, int *lda,int *ipiv, T *b, int *ldb, int *info);
// calls zgesv_ (T is complex) or dsyev_ (T is double) LAPACK routine
// parameters are the same as corresponding LAPACK routines
// http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen.html#ga442c43fca5493590f8f26cf42fed4044
// returns 1 if everything ok, or -1 otherwise
template<class T> void theev_(char *,char *,int *N,T *matr,int *lda,typename real_t<T>::data *x,T *work,int *lwork,
  typename real_t<T>::data *rwork,int *info);
#endif

#ifdef USE_PARDISO
/* calls PARDISO with recommended settings.
n is number of linear equations,
Matrix is stored in arrays a, ia, ja (see http://www.pardiso-project.org/manual/manual.pdf),
RHS column is stored in vect, solution is recorded to x.
returns 1 if everything ok, or -1 otherwise
*/
int pardiso(int n, double *a, int *ia, int *ja, double *vect, double *x);
#endif

using namespace std;

/* This function solves system of linear equations using Gauss method.
T is the variables type (integer, double, complex), matr is matrix data, vect is RHS vector,
num is number of linear equations.

matrix coefficients are stores in the following order:
matr_{ij}, where i is number of equation (matrix raw), j is number of variable (matrix column),
is stored in matr[j*num+i] or matr[j*lda+i] if lda!=0

Function changes values in matr and vect and records the solution to vect.
Returns 1 if unique solution is found, 0 if matrix is degenerate (one possible solution will be recorded to vec),
and -1 if system cannot be soved.
*/
template<class T>
int linsys(T *matr,T *vect,int num,T *x,int lda=0);

#ifdef USE_LAPACK
// inverse NxN matrix A 
int inverse(complex<double>* A, int N);
#endif

/*
Check the matrix for some special case of degeneracy and make it nondegenerate:
If there is i zero column, i zero raw and zero RHS i value, makes matr_{i,i} element as MATR_INFTY.
This methos can be used to exclude this special case of matrix degeneracy. */
template<class matr_t>
int make_nondegenerate(matr_t &matr){
  int deg=0;
  for(int i=0;i<matr.sz();i++){
    if(matr.v(i))
      continue;
    int j=0;
    for(;j<matr.sz();j++){
      if(matr.get_m(i,j) || (matr.get_m(j,i)))
        break;
    }
    if(j==matr.sz()){
      matr.set_m(i,i,MATR_INFTY);
      deg=1;
    }
  }
  return deg;
}


// normalize equations from system of linear equations on RHS values
// it can help to find accurate solution if matrix is well defined
template<class matr_t>
void normalize_val(matr_t &matr){
  for(int j=0;j<matr.sz();j++){
    if(matr.v(j)){
      for(int i=0;i<matr.sz();i++){
        matr.set_m(i,j,matr.get_m(i,j)/matr.v(j));
      }
      matr.v(j)=1;
    }
  }
}

// normalize equation from system of linear equations on maximal coefficient value
// it can help to find accurate solution if matrix is well defined
template<class matr_t>
void normalize_center(matr_t &matr){
  for(int j=0;j<matr.sz();j++){
    if(matr.get_m(j,j)){
      for(int i=0;i<matr.sz();i++){
        if(i==j)
          continue;
        // normalizing coeffs responsible for different variable numbers at the same equation
        matr.set_m(i,j,matr.get_m(i,j)/matr.get_m(j,j));
      }
      matr.v(j)/=matr.get_m(j,j);
      matr.set_m(j,j,1);
    }
  }
}

// normalize equation from system of linear equations on maximal coefficient value in the whole system
// it can help to find accurate solution if matrix is well defined
template<class matr_t>
void normalize_max(matr_t &matr){
  for(int j=0;j<matr.sz();j++){
    double max=0;
    for(int i=0;i<matr.sz();i++){
      if(fabs(matr.get_m(i,j))>max)
        max=fabs(matr.get_m(i,j));
    }
    if(max){
      for(int i=0;i<matr.sz();i++){
        // normalizing coeffs responsible for different variable numbers at the same equation
        matr.set_m(i,j,matr.get_m(i,j)/max);
      }
      matr.v(j)/=max;
    }
  }
}

// print system of linear equations (matrix of coefficients and RHS column)
template<class matr_t>
int print(matr_t &matr){
  FILE *f=fopen("matr.txt","wt");
  for(int i=0;i<matr.sz();i++){
    for(int j=0;j<matr.sz();j++){
//      fprintf(f,"%g\t",matr.get_m(i,j));
      fprintf(f,"%.100g",matr.get_m(i,j));
//      fprintf(f," \t");
      fprintf(f,"\n");
      fflush(f);
    }
    fprintf(f,"\n");
    fflush(f);
  }
  fflush(f);
  fclose(f);
  f=fopen("vect.txt","wt");
  for(int i=0;i<matr.sz();i++){
    fprintf(f,"%.100g",matr.v(i));
    fprintf(f,"\n");
    fflush(f);
  }
  fflush(f);
  fclose(f);
  return 1;
}

/* 
 Below we define some classes (linear solvers) that store coefficients in system of linear equations and solve it.
 This classes are templates from type T which can be integer, double, complex.

 To work with this class you should:
  - start_record
  - define system of linear equations using v (for RHS) and set_m (for matrix);
  - end_record;
  - solve the system using operator ().

  Class has bit flag mng: first bit is responsible for managing memory for RHS, 
  second bit is respondible for managing memory for matrix.

 There are two classes with described interface: linear_solver (uses LAPACK) and pardiso_solver (uses PARDISO).
*/


// This class is responsible only for storing and assigning RHS of system of linear equations.
// Linear solvers could be inherited from this class.
template<class T>
class vect_solver{
protected:
  T *vect; // RHS column
  int SZ; // number of linear equations
  int mng; // first bit is responsible for vect

public:

  vect_solver():vect(NULL),SZ(0),mng(0){}

  ~vect_solver(){
    if(vect && (mng&0x1))delete[]vect;
  }

  // define number of linear equations and allocate memory for RHS if first bit of mng is 1
  int init(int SZ_,int mng_){
    SZ=SZ_;
    mng=mng_;
    if(mng&0x1){
      vect = new T[SZ];
      for(int i=0;i<SZ;i++)vect[i]=0;
    }
    return 1;
  }

  // returns used memory size in bytes
  size_t data_size()const{
    return mng&1 ? SZ*sizeof(T) : 0;
  }

  void start_record(){}

  void end_record(){}

  // returns number of linear equations
  inline int sz()const{
    return SZ;
  }

  // i element of RHS column
  inline T &v(int i){
    return vect[i];
  }

  // get pointer on RHS column
  inline T* get_v(){
    return vect;
  }

  // use a subcolumn (i0, i0+n) of other solver (substituting memory pointer)
  void subvector(vect_solver &other,int i0,int n){
    vect=other.vect+i0;
    SZ=n;
  }

  // use a subcolumn (i0, i0+n) of other solver (direct copying)
  void copy_v(vect_solver &other,int i0,int n){
    for(int i=0;i<n;i++)
      v(i)=other.v(i+i0);
    SZ=n;
  }

  // fill RHS vector with zeros
  inline void zero_vector(){
    for(int i=0;i<SZ;i++)
      vect[i]=0;
  }

  // use the same memory which is already used in other solver (for memory optimization)
  inline void share(vect_solver &other){
    vect=other.vect;
  }

  // test if some of RSH column values is NaN
  int test_NaN(){
    for(int i=0;i<SZ;i++){
      if(::test_NaN(vect[i]))
        return 1;
    }
    return 0;
  }
};

// This solver solves system of linear equations using corresponding LAPACK routine (if USE_LAPACK is defined) 
// or by Gauss method using function linsys (otherwise).
template<class T>
class linear_solver: public vect_solver<T>{
  T *matr; // array for matrix elements
  // lda is responsible for order used to store matrix elements in matr (see realization of function set_m).
  // typically equal to SZ, but can be more than SZ
  int lda;
#ifdef USE_LAPACK
  int *ipiv; // memory reserved for LAPACK routines
#endif

  using vect_solver<T>::SZ;
  using vect_solver<T>::mng;
  using vect_solver<T>::vect;
  
public:

  linear_solver():matr(NULL),lda(0)
#ifdef USE_LAPACK
    ,ipiv(NULL)
#endif
  {}

  // calls vect_solver::init and allocates memory for matrix if second bit of mng_ is used
  int init(int SZ_,int mng_);

  ~linear_solver();

  size_t data_size(int flag)const{
    size_t sz=vect_solver<T>::data_size();

    if(flag&6 && mng&2){
      sz+=SZ*SZ*sizeof(T);
#ifdef USE_LAPACK
      sz+=sizeof(int)*SZ;
#endif
    }
    return sz;
  }

  // set matrix element with some value
  // i - variable number, j - equation number
  int set_m(int i,int j,T value){
    matr[i*lda+j]=value;
    return 1;
  }

  // get matrix element value
  inline T get_m(int i,int j) const{
    return matr[i*lda+j];
  }

  // use the same memory which is already used in other solver (for memory optimization)
  inline void share(linear_solver &other){
    vect_solver<T>::share(other);
    matr=other.matr;
#ifdef USE_LAPACK
    ipiv=other.ipiv;
#endif
  }

  // use submatrix (i0, j0, i0+n, j0+n) from other solver (substituting memory pointer)
  void submatrix(linear_solver &other,int i0,int j0,int n){
    matr=other.matr+i0*other.SZ+j0;
    lda=other.SZ;
#ifdef USE_LAPACK
    ipiv=other.ipiv+j0;
#endif
    SZ=n;
  }

  // use submatrix (i0, j0, i0+n, j0+n) from other solver (direct copying)
  void copy_m(linear_solver &other,int i0,int j0,int n){
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
        set_m(i,j,other.get_m(i0+i,j0+j));
    SZ=n;
  }

  // multiply part of matrix raw (i0, j0, i0+n, j0) on column x
  T multiply(int i0, int j0, int n, T *x){
    T res=0;
    for(int i=0;i<n;i++)
      res+=get_m(i0+i,j0)*x[i];
    return res;
  }

  // fill matrix with zeros
  inline void zero_matrix(){
    for(int i=0;i<SZ;i++)
      for(int j=0;j<SZ;j++)
        set_m(i,j,0);
  }

  // solve system of linear equations and record the solution to x
  // returns 1 if everything ok or -1 otherwise
  int operator()(T *x);

  // see description to ::make_nondegenerate
  int make_nondegenerate(){
    return ::make_nondegenerate(*this);
  }

  // see description to ::normalize_center
  void normalize_center(){
    ::normalize_center(*this);
  }

  // fill a, ia, ja that store matrix in PARDISO format
  // return number of nonzero elements
  int make_pardiso_matrix(T **a, int **ia, int **ja);

};

#ifdef USE_PARDISO

/*
This solver solves system of linear equations using corresponding PARDISO routine.
PARDISO works faster and requires less memory than LAPACK for sparce matrices.
Only nonzero matrix elements could be stored.
Matrix is stored in arrays a, ia, ja (see http://www.pardiso-project.org/manual/manual.pdf)
We omit documentation to some functions since it is the same as for linear_solver
*/
template<class T>
class pardiso_solver: public vect_solver<T>{

  T *a; // matrix values
  // index numeration starts with 1 in ia and ja
  int *ia; // start index in a and ja arrays for each raw
  int *ja; // columns corresponded to values a

  int *ind; // used for next step matrix filling
  int ci; // counter, is zero when matrix starts to be filled

  int sparce; // assumed maximal elements number at each row
  int jasz; // number of stored elements
  typedef pair<int,pair<T,int> > inf_t;

  int vfl; // 0 - array, 1 - vector
  int reuse; // 1 if reuse arrays a,ia,ja at next iteration
  int dfa; // delete_first_array after prepare if reuse

  int *stn; // elements number at each row
  inf_t *st; // there indeces, values and orders

  vector<inf_t> *vst;

  using vect_solver<T>::SZ;
  using vect_solver<T>::mng;
  using vect_solver<T>::vect;

public:

  pardiso_solver():a(NULL),ia(NULL),ja(NULL),sparce(0),jasz(0),
  vfl(0),reuse(1),dfa(0),stn(NULL),st(NULL),vst(NULL),ind(NULL),ci(0){}

  int init(int SZ_,int mem);

  ~pardiso_solver();

  void delete_first_record_arays();

  // set assumed maximal elements number at each row
  void set_sparce(int sp){
    sparce=sp;
  }

  int get_sparce(){
    return sparce;
  }

  // reuse arrays a,ia,ja at next iteration (which is start after start_record)
  // can be used if position of nonzero elements assumed to be the same to optimize filling their values
  void set_reuse(int r){
    reuse=r;
  }

  // 1 - vector, 2 - matrix, 4 - first record arrays
  size_t data_size(int flag)const{
    size_t sz=vect_solver<T>::data_size();

    if(flag&2)
      sz+=sizeof(int)*(SZ+1)+sizeof(int)*jasz+sizeof(T)*(jasz+SZ);

    if(flag&4 && mng&2){
      if(vfl){
        for(int i=0;i<SZ;i++)
          sz+=sizeof(vector<inf_t>)+sizeof(inf_t)*vst[i].size();
      }
      else
        sz+=sizeof(inf_t)*SZ*sparce+sizeof(int)*SZ;
    }

    return sz;
  }

  void start_record(){
    ci=0;
    if(!reuse)
      jasz=0;
  }

  int set_m(int i,int j,T val){
    if(!jasz){
      if(vfl){
        vst[j].push_back(make_pair(i,make_pair(0,ci++)));
        vst[j][vst[j].size()-1].second.first=val;
      }
      else{
        if(stn[j]>=sparce)
          return -1;
        st[j*sparce+stn[j]].first=i;
        st[j*sparce+stn[j]].second.second=ci++;
        st[j*sparce+stn[j]].second.first=val;
        stn[j]++;
      }
    }
    else{ // second cycle of matrix filling, now we know element position in array a
      ci++;
      a[ind[ci-1]]=val;
    }
    return 1;
  }

  T get_m(int i,int j) const;

  inline void share(pardiso_solver &other){
    vect_solver<T>::share(other);
    sparce=other.sparce;
    if(vfl){
      vst=other.vst;
    }
    else{
      stn=other.stn;
      st=other.st;
    }
//    init(SZ,2);
  }

  void submatrix(pardiso_solver &other,int i0,int j0,int n){
    copy_m(other,i0,j0,n);
  }

  void copy_m(pardiso_solver &other,int i0,int j0,int n);

  T multiply(int i0, int j0, int n, T *x);

  inline void zero_matrix(){
    if(jasz)
      return;
    if(vfl){
      for(int i=0;i<SZ;i++)
        vst[i].clear();
    }
    else{
      for(int i=0;i<SZ;i++)
        stn[i]=0;
      for(int i=0;i<SZ*sparce;i++){
        st[i].first=-1;
        st[i].second.first=0;
        st[i].second.second=-1;
      }
    }
  }

  void end_record();

  int operator()(T *x);

  int make_nondegenerate();

  void normalize_center();

  int test_NaN(){
    if(vect_solver<T>::test_NaN())
      return 1;
    for(int i=0;i<jasz;i++){
      if(::test_NaN(a[i]))
        return 1;
    }
    return 0;
  }
};

#endif

/* This function solves system of linear equations using Gauss-Seidel method 
(see http://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method).
 solv is some linear solver,
 snum is number of auxilliary blocks, sz is their sizes (if sz==NULL then sz is filled automatically),
 itnum is number of iterations.
 for common form of Gauss-Seidel method snum = num
 Gauss-Seidel method works as a Gauss methos if snum = 1 */
template<class T,class solver_t>
int gauss_seidel(solver_t &solv,int num,T *x,int snum=1,int *sz=NULL,int itnum=1);

// multidimensional newton-raphson method for solving system of nonlinear equations 
// (see http://en.wikipedia.org/wiki/Newton_Raphson)
template<class T,class solver_t>
class newton_raphson{
  T *dx,*y2; // auxiliary arrays
  int SZ; // size of system of nonlinear equations

  solver_t ls; // solver used to solve system of linear equations defined by Jacobian

public:
  newton_raphson(int SZ_);

  ~newton_raphson();

  // makes one iteration and changes x,
  // fun and function specify for nonlinear equation,
  // rsd is current residue, rsd2 is new residue,
  template<class fun_t>
  int operator()(T *x,fun_t *fun,int (fun_t::*function)(T *y),T &rsd,T &rsd2);

  // make num newton iterations.
  // can stop iterations if ratio between current and new residue is less than 1-stop_ratio
  template<class fun_t>
  int operator()(T *x,fun_t *fun,int (fun_t::*function)(T *y),int num,T stop_ratio=0){
    T rsd=-1,rsd2=-1;
    for(int i=0;i<num;i++){
      if(operator ()(x,fun,function,rsd,rsd2)<0)
        return -1;
      if(rsd2/rsd>=1-stop_ratio)
        break;
    }
    return 1;
  }

  // some sophisticated case
  template<class fun_t>
  int operator()(T **x,fun_t *fun,int (fun_t::**function)(T *y),int snum,int *sz,T &rsd,T &rsd2,int snum_gs=0,int *sz_gs=NULL,int itnum=1);
};


#ifdef USE_LAPACK

// This solver find eigenvalues of symmetric (hermitian) matrix using corresponding LAPACK routine
template<class T>
class eigen_solver{

  int SZ; // the order of matrix
  int mng; // if memory is maneged

  T *matr; // array for matrix elements
  // lda is responsible for order used to store matrix elements in matr (see realization of function set_m).
  // typically equal to SZ, but can be more than SZ
  int lda;

  // memory reserved for LAPACK routines
  T *work;
  typename real_t<T>::data *rwork;

public:

  eigen_solver():SZ(0),mng(0),matr(NULL),work(NULL),rwork(NULL),lda(0){}

  // allocates memory for matrix if mng_
  int init(int SZ_,int mng_=1);

  ~eigen_solver(){
    if(mng){
      delete[]matr;
      delete[]work;
      delete[]rwork;
    }
  }

  // set matrix element with some value
  // i - variable number, j - equation number
  int set_m(int i,int j,T value){
    matr[i*lda+j]=value;
    return 1;
  }

  // get matrix element value
  inline T get_m(int i,int j){
    return matr[i*lda+j];
  }

  // fill matrix with zeros
  inline void zero_matrix(){
    for(int i=0;i<SZ;i++)
      for(int j=0;j<SZ;j++)
        set_m(i,j,0);
  }

  // finds eigenvalues and record them to x
  // if eigenvect==true, stores eigenvectors as columns of matr
  // returns 1 if everything ok or -1 otherwise
  int operator()(typename real_t<T>::data *x, bool eigvect=false);

  // gets pointer on num-th eigenvector
  T *get_eigenvector(int num){
    return matr+SZ*num;
  }

};
#endif

#endif
