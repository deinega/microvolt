#include"linsolv.h"
#include<algorithm>

using namespace std;

///\en Definition of fabs for complex numbers
template<class T>
T fabs(const complex<T> &x){ return abs(x); }

template<class T>
int linsys(T *matr,T *vect,int num,T *x,int lda){
  int res=1;
  int *ind=new int[num];
  if(!ind)return -1; /* mememory allocation failure */
  for(int i=0;i<num;i++)ind[i]=i;

  if(!lda)lda=num;

  // проход по переменным вперед, делаем треугольную матрицу
//  for(int i=0;i<num;i++){
  for(int i=num-1;i>=0;i--){
    x[i]=0;
    int nm=i;
    double max=fabs(matr[i*lda+ind[i]]);
    // ищем уравнение, в котором имеется максимальный коэффициент для данной переменной
//    for(int k=i+1;k<num;k++){
    for(int k=i-1;k>=0;k--){
      double f=fabs(matr[i*lda+ind[k]]);
      if(f>max){
         max=f;
         nm=k;
      }
    }
    /*if(max<1e-20 && )return 0;*/
    /*else{*/
    // меняем местами уравнения так, чтобы на i-ом месте было уравнение с максимальным коэффициентом
    // для данной переменной
    if(nm!=i){
      int ti=ind[i];
      ind[i]=ind[nm];
      ind[nm]=ti;
      T t=vect[i];
      vect[i]=vect[nm];
      vect[nm]=t;
    }
    if(max>1e-20){
      T mul=T(1)/matr[i*lda+ind[i]];
      // вычитаю из всех последующих уравнений текущее уравнение
//      for(int k=i+1;k<num;k++){
      for(int k=i-1;k>=0;k--){
        if(fabs(matr[i*lda+ind[k]])!=0){
          vect[k]-=vect[i]*matr[i*lda+ind[k]]*mul;
//          for(int j=i+1;j<num;j++)
          for(int j=i-1;j>=0;j--)
            matr[j*lda+ind[k]]-=matr[j*lda+ind[i]]*matr[i*lda+ind[k]]*mul;
        }
        matr[i*lda+ind[k]]=0;
      }
    }
    else matr[i*lda+ind[i]]=0;
  }
  // проход по переменным назад
//  for(int i=num-1;i>=0;i--){
  for(int i=0;i<num;i++){
    if(fabs(matr[i*lda+ind[i]])<1e-20){
      if(fabs(vect[i])<1e-20){ // детерминант matr равен нулю
//        x[i]=vect[i]=1;
        x[i]=vect[i]=0;
        res=0;
      }
      else{
        delete[]ind;
        return -1; // система несовместна
      }
    }
    else x[i]=vect[i]=vect[i]/matr[i*lda+ind[i]];
    matr[i*lda+ind[i]]=1;
    // вычитаем решение из предыдущих уравнений
//    for(int j=0;j<i;j++){
    for(int j=num-1;j>i;j--){
      if(fabs(matr[i*lda+ind[j]])!=0){
        vect[j]-=matr[i*lda+ind[j]]*vect[i];
        matr[i*lda+ind[j]]=0;
      }
    }
  }
  delete[]ind;
  return res;
}

template<class T>
int linear_solver<T>::init(int SZ_,int mem){

  vect_solver<T>::init(SZ_,mem);
  lda=SZ;

  if(mng&0x2){
    matr = new T [SZ*SZ];
    for(int i=0;i<SZ*SZ;i++)
      matr[i]=0;
  }

#ifdef USE_LAPACK
  if(mng&0x2){
    ipiv = new int[SZ];
    for(int i=0;i<SZ;i++)ipiv[i]=0;
  }
#endif

  return 1;
}

template<class T>
linear_solver<T>::~linear_solver(){
  if(matr && (mng&0x2))delete[]matr;
#ifdef USE_LAPACK
  if(ipiv && (mng&0x2))delete[]ipiv;
#endif
}

template<class T>
int linear_solver<T>::operator()(T *x){
/*  double *a;
  int *ia,*ja;
  int jasz;
  make_pardiso_matrix(&a,&ia,&ja,jasz);
  if(pardiso(SZ,a,ia,ja,vect,x,jasz)<0)
    return -1;
  delete[]a;
  delete[]ia;
  delete[]ja;
  return 1;*/

#ifndef USE_LAPACK
  return linsys(matr,vect,SZ,x,lda);
#else
  int info;
  int nrhs = 1;
  int ldb = SZ;

  tgesv_(&SZ, &nrhs, matr, &lda, ipiv, vect, &ldb, &info);
  for(int i=0;i<SZ;i++)
    x[i]=vect[i];

  return info ? -1 : 1;
#endif
}

template<class T>
int linear_solver<T>::make_pardiso_matrix(T **a, int **ia, int **ja){
  int jasz=0;
  for(int j=0;j<SZ;j++){
    for(int i=0;i<SZ;i++){
      if(matr[i*lda+j])
        jasz++;
    }
  }

  *ia = new int[SZ+1];
  *ja = new int[jasz];
  *a = new T[jasz+SZ];

  jasz=0;
  (*ia)[0]=1;
  for(int j=0;j<SZ;j++){
    int inum=0;
    for(int i=0;i<SZ;i++){
      if(matr[i*lda+j]){
        (*a)[jasz]=matr[i*lda+j];
        (*ja)[jasz]=i+1;
        jasz++;
        inum++;
      }
    }
    (*ia)[j+1]=(*ia)[j]+inum;
  }
  return jasz;

}

#ifdef USE_PARDISO

template<class T>
int pardiso_solver<T>::init(int SZ_,int mem){

//  reuse=0; // test

  vect_solver<T>::init(SZ_,mem ? 3 : 0);
  if(mng&2){
    if(vfl){
      vst = new vector<inf_t>[SZ];
      for(int i=0;i<SZ;i++)
        vst[i].reserve(9);
    }
    else{
      if(!sparce)
        return -1;
      st = new pair<int,pair<T,int> >[SZ*sparce];
      stn = new int[SZ];
      zero_matrix();
    }
  }
  return 1;
}

template<class T>
pardiso_solver<T>::~pardiso_solver(){
  if(a)delete[]a;
  if(ia)delete[]ia;
  if(ja)delete[]ja;
  delete_first_record_arays();
}

template<class T>
void pardiso_solver<T>::delete_first_record_arays(){
  if(mng&2){
    if(stn)delete[]stn;
    stn=NULL;
    if(st)delete[]st;
    st=NULL;
    if(vst)delete[]vst;
    vst=NULL;
  }
}

template<class T>
T pardiso_solver<T>::get_m(int i,int j) const{
  if(ja[ia[j+1]-1-1] < i+1)
    return 0;
  for(int ind=ia[j]-1;ind<ia[j+1]-1;ind++){
    if(ja[ind] > i+1)
      break;
    if(ja[ind] == i+1)
      return a[ind];
  }
  return 0;
}

template<class T>
void pardiso_solver<T>::copy_m(pardiso_solver &other,int i0,int j0,int n){
  zero_matrix();
/*  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      T val=other.get_m(i0+i,j0+j);
      if(val)
        set_m(i,j,val);
    }
  }*/
  for(int j=0;j<n;j++){
    for(int ind=other.ia[j0+j]-1;ind<other.ia[j0+j+1]-1;ind++){
      int i=other.ja[ind]-1;
      if(i<i0)
        continue;
      else if(i>=i0+n)
        break;
      T val=other.a[ind];
      if(val){
        set_m(i-i0,j,val);
      }
    }
  }
  SZ=n;
}

template<class T>
T pardiso_solver<T>::multiply(int i0, int j0, int n, T *x){
  T res=0;
  for(int ind=ia[j0]-1;ind<ia[j0+1]-1;ind++){
    int i=ja[ind]-1;
    if(i<i0)
      continue;
    else if(i>=i0+n)
      break;
    T val=a[ind];
    res+=val*x[i-i0];
  }
  return res;
}

template<class T>
void pardiso_solver<T>::end_record(){
  
  if(!jasz){ // first matrix filling, define elements position in array a
    vector<int> sizes;
    size_t indsz=0;
    for(int i=0;i<SZ;i++){
      int rstn=0;
      if(vfl){
        for(size_t j=0;j<vst[i].size();j++){
          size_t k=j+1;
          for(;k<vst[i].size();k++){
            if(vst[i][j].first==vst[i][k].first){
              vst[i][j].first=-1;
              break;
            }
          }
          if(k==vst[i].size())
            rstn++;
        }
        indsz+=vst[i].size();
        jasz+=rstn;
        sizes.push_back(rstn);
      }
      else{
        for(int j=0;j<stn[i];j++){
//          if(st[sparce*i+j].second.first!=0){
          if(1){
            int k=j+1;
            for(;k<stn[i];k++){
              if(st[sparce*i+j].first==st[sparce*i+k].first){
                st[sparce*i+j].first=-1;
                break;
              }
            }
            if(k==stn[i])
              rstn++;
          }
          else
            st[sparce*i+j].first=-1;
        }
        indsz+=stn[i];
        stn[i]=rstn;
        jasz+=stn[i];
      }
    }
    if(a)delete[]a;
    if(ia)delete[]ia;
    if(ja)delete[]ja;
    ia = new int[SZ+1];
    ja = new int[jasz];
    a = new T[jasz+SZ];
    ind = new int[indsz];
    for(size_t i=0;i<indsz;i++)
      ind[i]=jasz; // for getting unused matrix 'a' element on default at calling (i,j) operator

    ia[0]=1;
    int ai=0;
    for(int i=0;i<SZ;i++){
      if(vfl){
        ia[i+1]=sizes[i]+ia[i];
        sort(vst[i].begin(),vst[i].end());
        for(size_t j=0;j<vst[i].size();j++){
          if(vst[i][j].first<0){
            continue;
          }
          a[ai]=vst[i][j].second.first;
          ja[ai]=vst[i][j].first+1;
          ind[vst[i][j].second.second]=ai;
          ai++;
        }
      }
      else{
        ia[i+1]=stn[i]+ia[i];
        sort(st+sparce*i,st+sparce*(i+1));
        for(int j=0;j<sparce;j++){
          if(st[sparce*i+j].first<0){
            continue;
          }
          a[ai]=st[sparce*i+j].second.first;
          ja[ai]=st[sparce*i+j].first+1;
          ind[st[sparce*i+j].second.second]=ai;
          ai++;
        }
      }
    }
    if(reuse && dfa)
      delete_first_record_arays();
  }

  for(int i=0;i<SZ;i++){
    a[jasz+i]=0;
  }
}

template<class T>
int pardiso_solver<T>::operator()(T *x){

  if(pardiso(SZ,a,ia,ja,vect,x)<0)
    return -1;

  for(int i=0;i<jasz;i++)
    a[i]=MATR_INFTY; // make nondegenrative

  start_record();

  return 1;
}

template<class T>
int pardiso_solver<T>::make_nondegenerate(){
  if(jasz)
    return 0;

  int deg=0;
  for(int i=0;i<SZ;i++){
    if(vect[i])
      continue;
    if((vfl && vst[i].size()) || (!vfl && stn[i]))
      continue;
    set_m(i,i,MATR_INFTY);
    deg=1;
  }
  return deg;
}

template<class T>
void pardiso_solver<T>::normalize_center(){
//  return;
  if(!jasz){
    for(int j=0;j<SZ;j++){
      T val=0;
      if(vfl){
      }
      else{
        for(int i=stn[j]-1;i>=0;i--){
          if(st[sparce*j+i].first==j){
            val=st[sparce*j+i].second.first;
            break;
          }
        }
        if(!val || stn[j]==1)
          continue;
        for(int i=stn[j]-1;i>=0;i--){
          if(st[sparce*j+i].first==j)
            st[sparce*j+i].second.first=1;
          else
            st[sparce*j+i].second.first/=val;
        }
      }
      vect_solver<T>::v(j)/=val;
    }
  }
  else{
    int ai=0;
    for(int i=0;i<SZ;i++){
      int rsz=ia[i+1]-ia[i];
      T val=0;
      for(int j=0;j<rsz;j++){
        if((ja[ai+j]-1)==i){
          val=a[ai+j];
        }
      }
      if(val && rsz>1){
        for(int j=0;j<rsz;j++){
          if((ja[ai+j]-1)==i)
            a[ai+j]=1;
          else
            a[ai+j]/=val;
        }
        vect[i]/=val;
      }
      ai+=rsz;
    }
  }
}

#endif

template<class T,class solver_t>
int gauss_seidel(solver_t &solv,int num,T *x,int snum,int *sz,int itnum){
  if(snum>num || snum<=0)
    return -1;
  int res=1;
  for(int i=0;i<num;i++)
    x[i]=0;
  int del=0;
  if(!sz){ // submatrices size is not specified
    sz = new int[snum];
    int part=num/snum;
    for(int i=0;i<snum-1;i++)
      sz[i]=part;
    sz[snum-1]=num%snum ? num-(snum-1)*part : part; // the rest
    del=1;
  }
  int szm=0; // maximal submatrix size
  for(int i0=0;i0<snum;i0++){
    if(sz[i0]>szm)
      szm=sz[i0];
  }

  solver_t gsolv;
#ifdef USE_PARDISO
  gsolv.set_reuse(0);
  gsolv.set_sparce(solv.get_sparce());
#endif
  gsolv.init(szm,1+2*(itnum>1));

  for(int it=0;it<itnum;it++){
    int szi=0;
    for(int i0=0;i0<snum;i0++){

      // copy sz[i0] elements starting from szi
      if(itnum==1)
        gsolv.submatrix(solv,szi,szi,sz[i0]);
      else
        gsolv.copy_m(solv,szi,szi,sz[i0]);
      gsolv.copy_v(solv,szi,sz[i0]);

      for(int i=0;i<sz[i0];i++){
        int szj=0;
        for(int j0=0;j0<snum;j0++){
          if(j0!=i0){
            // gauss-seidel changing of vector column
//            for(int j=0;j<sz[j0];j++)
//              gsolv.v(i)-=solv.get_m(szj+j,szi+i)*x[szj+j];
            gsolv.v(i)-=solv.multiply(szj,szi+i,sz[j0],x+szj);
          }
          szj+=sz[j0];
        }
      }
      gsolv.make_nondegenerate();
      gsolv.end_record();
//      int res1=gsolv(x+szi,sz[i0]);
      int res1=gsolv(x+szi);
      if(res1<0){
        if(del)delete[]sz;
        return res1;
      }
      else if(res1==0)
        res=0;
      szi+=sz[i0];
    }
  }

  if(del)delete[]sz;
  return res;
}

template<class T,class solver_t>
newton_raphson<T,solver_t>::newton_raphson(int SZ_):SZ(SZ_){
  ls.init(SZ,3);
  dx = new T[SZ];
  y2 = new T[SZ];
}

template<class T,class solver_t>
newton_raphson<T,solver_t>::~newton_raphson(){
  if(dx)delete[]dx;
  if(y2)delete[]y2;
}

template<class T,class solver_t>
template<class fun_t>
int newton_raphson<T,solver_t>::operator()(T *x,fun_t *fun,int (fun_t::*function)(T *y),T &rsd,T &rsd2){

//  T delta=.001;
  T delta=.0000001;
  T delta_1=1./delta;

  (fun->*function)(ls.get_v());

  T residue=0;
  for(int i=0;i<SZ;i++)
    residue+=ls.v(i)*ls.v(i);

  // calculating Jacobian
  for(int i=0;i<SZ;i++){
    T prev=x[i];
    x[i]+=delta;
    (fun->*function)(y2);
    for(int j=0;j<SZ;j++){
      T dy=y2[j]-ls.v(j);
      ls.set_m(i,j,dy*delta_1);
    }
    x[i]=prev;
  }

  ls.normalize_center();
  ls.make_nondegenerate();
  if(ls(dx)<0)
    return -1;

  // add residue calculation
  T residue2=0;

  T dx_mult=1;

  while(1){

    for(int i=0;i<SZ;i++)
      x[i]-=dx_mult*dx[i];

    (fun->*function)(y2);
    residue2=0;
    for(int i=0;i<SZ;i++)
      residue2+=y2[i]*y2[i];

    if(residue2<=residue)
      break;

    for(int i=0;i<SZ;i++)
      x[i]+=dx_mult*dx[i];

    dx_mult/=2; // make dx step less

  }

  rsd=residue;
  rsd2=residue2;

  return 1;
}

template<class T,class solver_t>
template<class fun_t>
int newton_raphson<T,solver_t>::operator()(T **x,fun_t *fun,int (fun_t::**function)(T *fy),
int snum,int *sz,T &rsd,T &rsd2,int snum_gs,int *sz_gs,int itnum){

  int szi=0;
  for(int i=0;i<snum;i++)
    szi+=sz[i];
  if(SZ!=szi)
    return -1;

  T delta=.00001;
  T delta_1=1./delta;

  szi=0;
  for(int j0=0;j0<snum;j0++){
    (fun->*function[j0])(ls.get_v()+szi);
    szi+=sz[j0];
  }

  szi=0;
  T residue=0;
  for(int i0=0;i0<snum;i0++){
    for(int i=0;i<sz[i0];i++){
      residue+=ls.v(szi+i)*ls.v(szi+i);
    }
    szi+=sz[i0];
  }

  szi=0;
  for(int i0=0;i0<snum;i0++){
    for(int i=0;i<sz[i0];i++){
      T prev=x[i0][i];
      x[i0][i]+=delta;
      int szj=0;
      for(int j0=0;j0<snum;j0++){
        (fun->*function[j0])(y2+szj);
        for(int j=0;j<sz[j0];j++){
          T dy=y2[szj+j]-ls.v(szj+j);
          ls.set_m(szi+i,szj+j,dy*delta_1);
        }
        szj+=sz[j0];
      }
      x[i0][i]=prev;
    }
    szi+=sz[i0];
  }

//  if(snum_gs){
  if(snum_gs>1){
    if(gauss_seidel(ls,SZ,dx,snum_gs,sz_gs,itnum)<0)
      return -1;
  }
  else{
//    ls.prepare();
    ls.prepare_center();
    ls.make_nondegenerate();
    if(ls(dx)<0)
      return -1;
  }

  T dx_mult=1;

  while(1){
    
    szi=0;
    for(int i0=0;i0<snum;i0++){
      for(int i=0;i<sz[i0];i++){
        x[i0][i]-=dx_mult*dx[szi+i];
      }
      szi+=sz[i0];
    }
//    break;

    szi=0;
    for(int j0=0;j0<snum;j0++){
      (fun->*function[j0])(y2+szi);
      szi+=sz[j0];
    }
    szi=0;
    T residue2=0;
    for(int i0=0;i0<snum;i0++){
      for(int i=0;i<sz[i0];i++){
        residue2+=y2[szi+i]*y2[szi+i];
      }
      szi+=sz[i0];
    }
    if(residue2<=residue)
      break;

    szi=0; // move back
    for(int i0=0;i0<snum;i0++){
      for(int i=0;i<sz[i0];i++){
        x[i0][i]+=dx_mult*dx[szi+i];
      }
      szi+=sz[i0];
    }

    dx_mult/=2; // make dx step less
  }

  return 1;
}

#ifdef USE_LAPACK

template<class T>
int eigen_solver<T>::init(int SZ_,int mng_){

  SZ=SZ_;
  mng=mng_;
  lda=SZ;

  if(mng){
    matr = new T [SZ*SZ];
    for(int i=0;i<SZ*SZ;i++)
      matr[i]=0;

    work = new T[2*SZ-1];
    for(int i=0;i<2*SZ-1;i++)work[i]=0;

    rwork = new typename real_t<T>::data[3*SZ-2];
    for(int i=0;i<3*SZ-2;i++)rwork[i]=0;
  }
  return 1;
}

template<class T>
int eigen_solver<T>::operator()(typename real_t<T>::data *x, bool eigvect){

  int lwork=2*SZ-1;
  int info;

  theev_((char *)(eigvect ? "V" : "N"),(char *) "L",&SZ,matr,&lda,x,work,&lwork,rwork,&info);

  return info ? -1 : 1;
}

#endif
