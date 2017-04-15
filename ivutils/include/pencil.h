/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *   $Revision: 1.1 $
 *   $Date: 2012/12/06 02:24:04 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/pencil.h,v 1.1 2012/12/06 02:24:04 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/pencil.h,v $
$Revision: 1.1 $
$Author: lesha $
$Date: 2012/12/06 02:24:04 $
*/
/*s****************************************************************************
 * $Log: pencil.h,v $
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.20  2012/10/16 15:33:11  lesha
 * documentation
 *
 * Revision 1.19  2012/06/16 01:06:04  valuev
 * sync with GridMD project
 *
 * Revision 1.14  2009/07/24 05:08:46  valuev
 * Sync with FDTD, added molecule setup
 *
 * Revision 1.18  2009/05/28 07:49:00  valuev
 * updated chemosensor
 *
 * Revision 1.17  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.11  2008/09/02 22:26:34  valuev
 * chain loader
 *
 * Revision 1.10  2008/02/21 16:35:03  valuev
 * made compilable under Linux
 *
 * Revision 1.9  2007/11/20 22:15:54  valuev
 * Added parametric project
 *
 * Revision 1.8  2007/07/09 21:29:07  valuev
 * plasma with wave packets
 *
 * Revision 1.12  2007/06/25 15:28:02  valuev
 * corrected pencil and buffer transfers
 *
 * Revision 1.11  2007/06/22 08:39:53  valuev
 * New pencils in interpolation
 *
 * Revision 1.9  2007/06/01 11:40:32  valuev
 * corrected oblique incidence
 *
 * Revision 1.8  2007/05/30 21:19:37  lesha
 * 2 errors are fixed
 *
 * Revision 1.7  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.6  2007/02/16 10:16:51  valuev
 * allowed array mngptr
 *
 * Revision 1.7  2007/02/16 09:40:32  valuev
 * Added Nudged Elastic Band saddle point search
 *
 * Revision 1.6  2006/12/20 14:29:33  valuev
 * Updated workflow, sync with FDTD
 *
 * Revision 1.5  2006/11/07 13:50:12  valuev
 * Added geometry dumps for gnuplot
 *
 * Revision 1.4  2006/10/27 20:41:01  valuev
 * Added detectors sceleton. Updated some of ivutils from MD project.
 *
 * Revision 1.4  2006/09/26 10:59:42  valuev
 * Added nonorthogonal TB (Menon-Subbaswamy)
 *
 * Revision 1.3  2006/07/21 16:22:03  valuev
 * Added Tight Binding for graphite+O
 *
 * Revision 1.2  2006/04/26 12:12:01  valuev
 * Fixed Neighbour Lists (double-single counting), added twostep NL scheme, added Step2 to mdtutorial (use of template potentials), added DelAtom to mdStructure
 *
 * Revision 1.1  2006/04/17 07:36:11  valuev
 * Added new files (pairpot, pencil, graspapp/vs7)
 *
*******************************************************************************/

#ifndef PENCIL_H
#define PENCIL_H

# include <algorithm> 

using namespace std;


///\en making container of a pointer (to look like a vector)
template <class T>
class ptr_container {
  T* ptr;
  int master;
public:
  ptr_container(size_t size=0):master(1){
    if(size)ptr= new T[size];
    else ptr=NULL;
  }
  ptr_container(const T* first, const T* last):master(1){
    size_t n=last-first;
    if(n!=0){
      ptr= new T[n];
      memcpy(ptr,first,n*sizeof(T));
    }
    /*if(n){
      ptr= new T[n];
      while(first!=last)
        *ptr++=*first++;
    }*/
  }
  
  
  ptr_container(T* data, int managed=0): ptr(data),master(managed) {}

  void assign(const T* first, const T* last){
    size_t n=last-first;
    if(n!=0)
      memcpy(ptr,first,n*sizeof(T));
  }

  ~ptr_container(){
    if(master && ptr) delete [] ptr;
  }
  typedef T* iterator;
  iterator begin() const { return ptr;}

  T& operator[](int i){ return ptr[i]; }
};

//\en this class is used to convert argument of pencil constructor
/// which can be either T* or container_t &
template<class container_t>
class container_arg{
  container_t *pcont;
public:
  container_arg(container_t *pcont_):pcont(pcont_){}
  container_t &operator*(){return *pcont;} 
  container_t *operator->(){return pcont;} 
  container_t *persistent_ptr() const {return &pcont;}
};

template<class T>
class container_arg< ptr_container<T> >{
  ptr_container<T> cont;
public:  
  typedef ptr_container<T> container_t;
  container_arg(T *data):cont(data){}
  container_t &operator*(){return cont;} 
  container_t *operator->(){return &cont;} 
  container_t *persistent_ptr() const {return new ptr_container<T>(cont.begin(),1);}
};


///\en class for controlling one-dimensional data arrays of type T
/// the array can be managed and copied by reference (ref counting used)
/// compatible with vector<T> used as container
template <class T, class container_t=ptr_container<T> >
class pencil {
protected:  
  
  void clear_rd(){
    if(rd){
      (*rd)--;
      if(*rd==0){
        delete rd;
        if(cont)delete cont;
      }
      rd=NULL;
    }
    else{
      if(cont) delete cont;
    }
    cont=NULL;
  }

  container_t *cont; // container
  int *rd; // how many other pencils manage the same container
  typename container_t::iterator beg;

public:
  typedef typename container_t::iterator container_it;
  int nmax; /// capacity of container
  int n;  /// size of container (actual elements number)
  
  pencil(int snmax=0, int sn=0): nmax(snmax), n(sn){
    if(n>nmax)n=nmax;
    if(nmax){
      cont = new container_t(nmax);
      beg = cont->begin();
      rd = new int(1);
    }
    else{
      cont = NULL;
      rd = NULL;
    }
  }

  /// Constructor form existing container.
  /// Container pointer is to be used as the first argument, f.e.
  /// T*  in case of pointer container
  /// or vector<T> * in case of vector container
  /// Container_arg makes the proper argument conversion
  /// if managed=0 will not delete ptr, 1-deletes, 2-makes a copy of sn elements and does not remember ptr
  pencil(container_arg<container_t> arg,int snmax, int sn=0, int managed=0): nmax(snmax), n(sn){
    if(n>nmax)n=nmax;
    if(managed==2){
      cont = new container_t(nmax);
      beg= cont->begin();
      rd= new int(1);
      // copying
      cont->assign(arg->begin(),arg->begin()+n);
      //memcpy(ptr,ptr_,nmax*sizeof(T));
    }
    else if(managed==1){
      cont = arg.persistent_ptr();
      rd= new int(1);
      beg = cont->begin();
    }
    else{
      cont = NULL; // will not delete cont
      rd = NULL;
      beg = (*arg).begin();
    }
  } 

  /// creates unmanaged reference to the sequence
  pencil(container_it first, container_it last, int snmax):nmax(snmax){
    n=(int)(last-first);
    if(n>nmax)n=nmax;
    cont=NULL;
    rd=NULL;
    beg=first;
  }

  /// resizes the data: if newnmax==0, dereferences the old array
  ///                   if newnmax<=nmax, keeps the old one
  ///                   if newnmax>nmax  dereferences the old array, allocates new, changes to master mode
  int resize(int newnmax, int newn=-1){
    if(newn<0)newn=n;
    if(newn>newnmax)return -1;
    if(newnmax==nmax)return 0;
    if(newnmax<nmax && newnmax!=0){
      nmax=newnmax;
      return 0;
    }

    int *rd1=NULL;
    container_t *cont1=NULL;
    if(newnmax){
      cont1 = new container_t(newnmax);
      rd1 = new int(1);
    }
    int nn=(n>newnmax? newnmax : n); //, i ;
    if(nn)
      cont1->assign(beg,beg+nn);
    /*for(i=0;i<nn;i++){
      (*cont1)[i]=(*cont)[i];
    }*/
    clear_rd();
    cont=cont1;
    if(cont)
      beg=cont->begin();
    nmax=newnmax;
    rd=rd1;
    n=newn;
    return 1;
  }

 
  ///\en copy constructor
  pencil(const pencil &other){
    cont=other.cont;
    rd=other.rd;
    beg=other.beg;
    if(rd)(*rd)++;
    n=other.n;
    nmax=other.nmax;
  }

  ///\en reference copy
  pencil &operator=(const pencil &other){
    if(&other!=this){
      clear_rd();
      cont=other.cont;
      rd=other.rd;
      beg=other.beg;
      if(rd)(*rd)++;
      n=other.n;
      nmax=other.nmax;
    }
    return *this;
  }

  ///\en copies data from other pencil
  /// mode=0:  if this size is too small, creates new reference data!
  ///          in which case other referring pencils will not be updated
  /// mode=1:  if this size is too small, copies the first portion of the data only
  ///          
  pencil &copy_data(const pencil &other, int mode=0){
    n=other.n;
    if(nmax<other.n && mode==0){ // must reallocate
      clear_rd();
      nmax=other.nmax;
      cont = new container_t(other.beg,other.beg+nmax);
      rd = new int(1);
      beg=cont->begin();
    }
    else{// just copying
      if(nmax)
        cont->assign(other.beg,other.beg+nmax);
    }
    return *this;
  }

  ///\en same as begin(), retained for compatibility
  container_it get_ptr() const{
    return beg;
  }

  T& operator[](int i) const {
    return *(beg+i);
  }

  pencil &swap(pencil &other){
    std::swap(n,other.n);
    std::swap(nmax,other.nmax);
    std::swap(rd,other.rd);
    std::swap(cont,other.cont);
    std::swap(beg,other.beg);
    return *this;
  }

  ///\en shifts array by num elements starting from ind
  int shift_arr(int ind, int num){
    int i;
    if(ind<0 || ind>=n)return -1;
    n+=num;
    if(num<0){
      if(n<0)n=0;
      for(i=ind;i<n;i++){
        cont[i]=cont[i-num];
      }
      return 1;
    }
    else if(num>0){
      ind+=num;
      if(n>=nmax)n=nmax;
      for(i=n-1;i>=ind;i--){
        cont[i]=cont[i-num];
      }
      return 2;
    }
    return 0;
  }

  ~pencil(){
    clear_rd();
  }

  int size() const{
    return n;
  }

  ///\en gets the number of copies
  int num_copies(){
    if(rd)return *rd;
    else return 0;
  }

  typedef container_it iterator;

  iterator begin() const { return beg; }
  iterator end() const { return beg+n; }

};


/* old pencil class

//e class for controlling one-dimensional data arrays of type T
template <class T>
class pencil{
protected:  
  int master;

  void clear_rd(){
    if(rd){
      rd->unref(master);
      if(rd->ref==0)delete rd;
      rd=NULL;
    }
  }

public:
  int nmax;
  int n;
  
  struct ref_data{
    int ref;
    T *ptr;
    ref_data(T *sptr, int sref=1):ptr(sptr),ref(sref){}
     
    void unref(int master=1){
      ref--;
      if(!ref && ptr && master){
        delete [] ptr;
        ptr=NULL;
      }
    }
  } *rd;
  
  pencil(int snmax=0, int sn=0): nmax(snmax), n(sn), master(1){
    if(n>nmax)n=nmax;
    if(nmax){
      T *ptr = new T[nmax];
      rd = new ref_data(ptr, 1);
    }
    else{
      rd = NULL;
    }
  }
  /// if managed=0 will not delete ptr, 1-deletes, 2-makes a copy and does not remember ptr
  pencil(T *ptr,int snmax, int sn=0, int managed=0): nmax(snmax), n(sn), master(managed){
    if(n>nmax)n=nmax;
    if(master==2){
      rd= new ref_data(new T[nmax],1);
      // copying
      memcpy(rd->ptr,ptr,nmax*sizeof(T));
      master=1;
    }
    else rd = new ref_data(ptr, 1); // will not delete ptr if master=0
  } 

  /// resizes the data: if newnmax==0, dereferences the old array
  ///                   if newnmax<=nmax, keeps the old one
  ///                   if newnmax>nmax  dereferences the old array, allocates new, changes to master mode
  int resize(int newnmax, int newn=-1){
    if(newn<0)newn=n;
    if(newn>newnmax)return -1;
    if(newnmax==nmax)return 0;
    if(newnmax<nmax && newnmax!=0){
      nmax=newnmax;
      return 0;
    }

    ref_data *rd1=NULL;
    if(newnmax){
      T *ptr = new T[newnmax];
      rd1 = new ref_data(ptr, 1);
    }
    int i, nn=(n>newnmax? newnmax : n);
    for(i=0;i<nn;i++){
      rd1->ptr[i]=rd->ptr[i];
    }
    clear_rd();
    nmax=newnmax;
    rd=rd1;
    master=1;
    n=newn;
    return 1;
  }

 
  //e copy constructor
  pencil(const pencil &other){
    rd=other.rd;
    if(rd)rd->ref++;
    n=other.n;
    nmax=other.nmax;
    master=other.master;
  }

  //e reference copy
  pencil &operator=(const pencil &other){
    if(&other!=this){
      clear_rd();
      rd=other.rd;
      if(rd)rd->ref++;
      n=other.n;
      nmax=other.nmax;
      master=other.master;
    }
    return *this;
  }

  //e copies data from other pencil
  //e warning: if this size is too small, creates new reference data!
  //e          in which case other referring pencils will not be updated
  pencil &copy_data(const pencil &other){
    n=other.n;
    if(nmax<other.n){ // must reallocate
      clear_rd();
      nmax=other.nmax;
      T *ptr = new T[nmax];
      rd = new ref_data(ptr, 1);
      master = 1;
    }
    // copying
    if(n!=0)memcpy(rd->ptr,other.rd->ptr,n*sizeof(T));
    return *this;
  }

  T* get_ptr() const{
    if(!rd)return NULL;
    return rd->ptr;
  }

  T& operator[](int i) const {
    return *(rd->ptr+i);
  }

  pencil &swap(pencil &other){
    int tn=other.n;
    other.n=n;
    n=tn;
    int tnmax=other.nmax;
    other.nmax=nmax;
    nmax=tnmax;
    int tmaster=other.master;
    other.master=master;
    master=tmaster;
    ref_data* trd=other.rd;
    other.rd=rd;
    rd=trd;
    return *this;
  }

  //e shifts array by num elements starting from ind
  int shift_arr(int ind, int num){
    int i;
    if(ind<0 || ind>=n)return -1;
    n+=num;
    if(num<0){
      if(n<0)n=0;
      for(i=ind;i<n;i++){
        rd->ptr[i]=rd->ptr[i-num];
      }
      return 1;
    }
    else if(num>0){
      ind+=num;
      if(n>=nmax)n=nmax;
      for(i=n-1;i>=ind;i--){
        rd->ptr[i]=rd->ptr[i-num];
      }
      return 2;
    }
    return 0;
  }

  ~pencil(){
    clear_rd();
  }

  int size() const{
    return (size_t)n;
  }

};
# endif
*/

# endif
