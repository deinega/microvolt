/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.4 $
 *   $Date: 2014/03/09 05:34:17 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/refobj.h,v 1.4 2014/03/09 05:34:17 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/refobj.h,v $
$Revision: 1.4 $
$Author: lesha $
$Date: 2014/03/09 05:34:17 $
*/
/*s****************************************************************************
 * $Log: refobj.h,v $
 * Revision 1.4  2014/03/09 05:34:17  lesha
 * *** empty log message ***
 *
 * Revision 1.3  2013/07/03 15:41:50  valuev
 * added array_mngptr
 *
 * Revision 1.2  2013/06/19 16:58:23  valuev
 * added vector delete
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.39  2012/07/26 21:33:37  valuev
 * added rt_angular experiment
 *
 * Revision 1.37  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.36  2012/03/23 07:32:34  lesha
 * comments
 *
 * Revision 1.35  2012/03/04 01:50:08  lesha
 * default constructor is added to mngarg
 *
 * Revision 1.34  2009/12/23 10:22:14  lesha
 * comments are added
 *
 * Revision 1.33  2009/05/19 21:50:17  valuev
 * Added TestRay for plane
 *
*******************************************************************************/
#ifndef _REFOBJ_H
#define _REFOBJ_H

///     @file refobj.h \brief Smart pointers.
# include <utility>
#include <stddef.h>
#include <map>
#include <vector>

using namespace std;

///\ru служебный класс, нужный для удобства конструирования умных указателей (mng_ptr и т. д.)
/// содержит указатель и флаг integer
///\en Auxiliary class used as an argument for condtructor and functions set / reset 
/// of smart pointers (mng_ptr, sh_ptr).
/// Contains pointer at dynamic object and integer flag.
template<class T>
class mngarg: public pair<T *, int>{
public:
  typedef pair<T *, int> base_t;
  using base_t::first;
  using base_t::second;

  mngarg(T *ptr=NULL, int managed=0): pair<T*,int>(ptr,managed){}
  template<class A>
  mngarg(const mngarg<A> &arg): pair<T*,int>(arg.first,arg.second){}

  T* ptr() const{
    return first;
  }
  T* operator->() const{
    return first;
  }
  T& operator*() const{
    return *first;
  }
  T& operator[] (int i) const{
    return *(first+i);
  }
};

template<class T>
mngarg<T> make_mngarg(T *ptr, int managed=1){
  return mngarg<T>(ptr,managed);
}

template<class T>
struct delete_ptr{
  void operator()(T *ptr){
    delete ptr;
  }
};

template<class T>
struct delete_arr{
  void operator()(T *ptr){
    delete [] ptr;
  }
};


///\en
/// Managed pointer. Used to control memory allocated for dynamic objects.
/// Contains pointer on dynamic object and integer flag, which define the action
/// with dynamic object while pointer adress is changed:
/// 0 do not delete, 1 delete, 0x8 delete as array.
/// Managed pointers can not be used as elements of STL containers!\n
///\ru
/// Managed pointer. Нужен для автоматического управления динамической памятью.
/// Содержит указатель на динамический объект и флаг integer, говорящий о том, 
/// что делать с динамическим объектом при смене адреса указателя:
/// 0 do not delete, 1 delete, 
/// 2 copy and delete (NOT IMPLEMENTED, requires copy constructor), 0x8 -- delete as array.
template<class T>
class mngptr: public pair<T *, int>{
public:
  typedef pair<T *, int> base_t;
  typedef T *pointer;
  using base_t::first;
  using base_t::second;

  mngptr(T* ptr=NULL, int managed=0): pair<T*,int>(ptr,managed){}
  mngptr(const mngarg<T> &arg): pair<T*,int>(arg.first,arg.second){}
  const mngptr &operator=(const mngarg<T> &arg){
    reset(arg.first,arg.second);
    return *this;
  }
  void reset(T* ptr=NULL, int managed=0){
    if(second && first && first!=ptr){
      if(second&0x8)delete [] first;
      else delete first;
    }
    first=ptr;
    second=managed;
  }
  void reset(const mngarg<T> &arg){
    reset(arg.first,arg.second);
  }
  T* ptr() const {
    return first;
  }
  T* operator->() const {
    return first;
  }
  T& operator*() const{
    return *first;
  }
  T& operator[] (int i) const{
    return *(first+i);
  }
  T& operator[] (size_t i) const{
    return *(first+i);
  }
  int managed() const {
    return second;
  }
  ~mngptr(){
    reset();
  }
};




///\en
/// Array of managed pointers. Used to group managed pointers into variable size array
/// based on std::vector.
template<class T>
class array_mngptr{
  std::vector<mngarg<T> > pointers;
public:
  ///\en get current size
  size_t size() const {
    return pointers.size();
  }

  ///\en add an element to the end of array
  ///    \return index of the added element
  size_t push_back(const mngarg<T> & arg){
    pointers.push_back(arg);
    return size()-1;
  }

  ///\en Get pointer at current index
  const T* operator[](size_t i) const {
    return pointers[i].ptr();
  }

  ///\en Get const pointer at current index
  T* operator[](size_t i){
    return pointers[i].ptr();
  }

  ///\en Change value at current index, deletes existing pointer if necessary
  void reset(size_t i, const mngarg<T> & arg){
    mngptr<T> mptr = pointers[i];
    mptr.reset(arg);
    pointers[i] = mngarg<T>(mptr.ptr(),mptr.managed());
  }

  ///\en Clears array, deleting all pointers if necessary
  void clear(){
    for(size_t i=0;i<size();i++){
      reset(i,NULL);
    }
    pointers.clear();
  }

  ///\en Destructor, deletes all pointers if necessary
  ~array_mngptr(){
    clear();
  }
 
};

/* example for mngarr:
  array_mngptr<int> mngarr;
  mngarr.push_back(make_mngarg(new int(1)));
  int b = 2;
  mngarr.push_back(&b);
  *mngarr[1]=3;
  printf("%d",*mngarr[0]);
  mngarr.clear();
*/

# if 0
template <template<class _type> class cont_tt, class T>
class refcontainer: public cont_tt<T *>{
protected:
  int man;
public:
  typedef cont_tt< T * > base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::const_iterator const_iterator;

  refcontainer(int smanaged=0):man(smanaged){}
  refcontainer(size_t n, int smanaged=0):base_t(n),man(smanaged){}

  void set_managed(int sman){
    man=sman;
  }

  ~refcontainer(){
    if(man){
      size_t i, n=base_t::size();
      for(i=0;i<n;i++)
        if((*this)[i]) delete (*this)[i];
    }
  }
};
# endif

///\en
/// Vector of dynamic objects and flag, which define the action
/// with dynamic objects if function clear is called:
/// 0 - do not delete, nonzero - delete, 8 (bit flag set) - delete as array
///\ru
/// Vector указателей на динамические объекты, который сам следит за их удалением (при man==1),
/// установленный битовый флаг 0x8 означает удаление как массива.
template <class T>
class refvector: public std::vector<T *>{
protected:
  int man;
public:
  typedef vector<T*> base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::const_iterator const_iterator;

  refvector(int man_=0):man(man_){}
  refvector(size_t n, int man_):base_t(n),man(man_){}

  void set_managed(int man_){
    man=man_;
  }

  ~refvector(){
    clear();
  }

  void clear(){
    if(man){
      for(typename base_t::iterator it=base_t::begin(), e=base_t::end();it!=e;++it){
        if(*it){
          if(man&0x8)
            delete [] (*it);
          else 
            delete (*it);
        }
      }
    }
    base_t::clear();
  }

  iterator erase(iterator it){
    if(man && *it){
      if(man&0x8)
        delete [] (*it);
      else 
        delete (*it);
    }
    return base_t::erase(it);
  }

};

///\en
/// Map of dynamic objects and flag, which define the action
/// with dynamic objects if function clear is called:
/// 0 - not delete, 1 - delete
///\ru
/// Map указателей на динамические объекты, который сам следит за их удалением (при man==1)
template <class key_tt, class T>
class refmap: public std::map<key_tt, T *>{
protected:
  int man;
public:
  typedef std::map<key_tt, T * > base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::const_iterator const_iterator;

  refmap(int man_=0):man(man_){}
  refmap(size_t n, int man_):base_t(n),man(man_){}

  void set_managed(int man_){
    man=man_;
  }

  ~refmap(){
    clear();
  }

  void clear() {
    if(man){
      for(typename base_t::iterator it=base_t::begin(), e=base_t::end();it!=e;++it)
        if(it->second) delete it->second;
    }
    base_t::clear();
  }

  iterator erase(iterator it){
    if(man && it->second)
      delete it->second;
    return base_t::erase(it);
  }
};

///\en
/// Smart pointer. Used to control memory allocated for dynamic objects.
/// Contains pointer p on dynamic object and pointer on integer number num, 
/// which specify the number of other smart pointers with refers at the same dynamic object p.
/// If *num becomes zero (no one refers to dynamic object), dynamic object should be deleted.
///\ru
/// Умный указатель, содержащий в себе указатель на динамический объект (p) и указатель на счетчик ссылок на него (num)
/// При нулевом количестве ссылок, указывающих на динамический объект, динамический объект удаляется
template<class T, class delete_t=delete_ptr<T> >
class shptr{
  template<class Y, class Z> friend class shptr;
  T *p;
  int *num; //if num==NULL than p is not managed

  void set(T *p_, int managed){
    p=p_;
    if(p&&managed){
      num=new int;
      *num=1;
    }
    else num=NULL;
  }
  template<class Y,class Z>
  void set(const shptr<Y,Z> &other){
    p=other.p;
    if(p){
      num=other.num;
      if(num)(*num)++;
    }
    else num=NULL;
  }

public:
  shptr(T* p=NULL, int managed=1){
    set(p,managed);
  }
  shptr(const mngarg<T> &arg){
    set(arg.first,arg.second);
  }
  shptr(const shptr &other){
    set(other);
  }
  template<class Y,class Z>
  shptr(const shptr<Y,Z> &other){
    set(other);
  }

  void reset(T *p_, int managed=1) {
    if(p!=p_){
      free();
      set(p_,managed);
    }
  }
  void reset(const shptr &other) {
    if(this!=&other){
      free();
      set(other);
    }
  }

  const shptr &operator=(T *p){
    reset(p,0);
    return *this;
  }
  const shptr &operator=(const mngarg<T> &arg){
    reset(arg.first,arg.second);
    return *this;
  }
  const shptr &operator=(const shptr &other){
    reset(other);
    return *this;
  }
  template<class Y,class Z>
  const shptr &operator=(const shptr<Y,Z> &other){
    reset(other);
    return *this;
  }

  virtual ~shptr(){
    free();
  }
  void free(){
    if(p){
      if(num){
        (*num)--;
        if((*num)==0){
          delete_t()(p);
          delete num;
        }
        num=NULL;
      }
      p=NULL;
    }
  }

  bool valid() const {
    return p!=NULL;
  }

  T* ptr() const {
    return p;
  }
  T* operator->() const {
    return p;
  }
  T& operator*() const{
    return *p;
  }
};

#endif
