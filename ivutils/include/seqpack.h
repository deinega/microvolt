/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005-2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   This source code is Free Software and distributed under the terms of wxWidgets license (www.wxwidgets.org) 
 *
 *   $Revision: 1.7 $
 *   $Date: 2013/11/02 17:14:53 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/seqpack.h,v 1.7 2013/11/02 17:14:53 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/seqpack.h,v $
$Revision: 1.7 $
$Author: lesha $
$Date: 2013/11/02 17:14:53 $
*/
/*s****************************************************************************
 * $Log: seqpack.h,v $
 * Revision 1.7  2013/11/02 17:14:53  lesha
 * nothing important
 *
 * Revision 1.6  2013/09/07 23:56:34  lesha
 * *** empty log message ***
 *
 * Revision 1.5  2013/09/04 16:05:35  valuev
 * continued restructure for 64-bits (added index_t)
 *
 * Revision 1.4  2013/09/04 15:22:56  valuev
 * restructure for 64-bits (added index_t)
 *
 * Revision 1.3  2013/05/14 13:57:10  belousov
 * fixed compatibility with gcc
 *
 * Revision 1.2  2013/04/13 08:44:26  valuev
 * added unfinished sequence message
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.43  2012/11/10 23:59:18  lesha
 * GetTruncatedCone is added
 *
 * Revision 1.42  2012/10/20 00:32:46  lesha
 * operator [] is added to data_unpacked and group_unpacked
 *
 * Revision 1.41  2012/10/17 00:59:38  lesha
 * switching between packing and non-packing by template pack_t
 *
 * Revision 1.40  2012/10/16 19:42:58  lesha
 * documentation
 *
 * Revision 1.39  2012/06/16 01:06:04  valuev
 * sync with GridMD project
 *
 * Revision 1.38  2010/10/07 10:30:22  lesha
 * restorer is commented
 *
 * Revision 1.37  2010/09/29 10:40:49  lesha
 * comments
 *
 * Revision 1.36  2010/04/17 19:32:39  valuev
 * fixed group_pack
 *
*******************************************************************************/


#ifndef SEQPACK_H
#define SEQPACK_H

/*  \en @file seqpack.h
    Classes for packing different types of data (packers).
    All these classes share the same interface. This allows to switch easily from one packer to another.
    Before packing you should call function start_record.
    While packing data you call function record with data to be packed as an argument.
    After you finished packing you call function end_record. After this moment you cannot continue packing.
    However, you can start packing from the beginning, calling function clear and then start_record.

    Access to packed data is provided by STL-like iterators.
    Functions begin and end return iterator at the beginning and end of packed data sequence.
    Iterators have operators ++ (increment), * (take value), 
    != (check for unequality, this operator works for checking unequality with end iterator).

    This is illustrative example of data packing and consecutive access to packed data 
    by some packer of the packer_t type:

    // recording sequence
    packer.start_record(); 
    packer.record(some_data);
    /// ... apply it some times
    packer.end_record();

    // extracting sequence
    for(packer_t::iterator it=packer.begin(), e=packer.end(); it!=e; ++it)
      packer_t::value_t val=*it;

    Packers can also pack groups of data (subsequences of elements) 
    by function next_group with iterator on a first element of subsequence and sequence size as arguments.
    Some packers have a functionality to store size of such subsequences.
    Their iterators have function plus to increment on the size of the whole sequence.
    Besides, one can apply to these iterators functions 
    get_group_count, which return size of the rest of the sequence, and
    get_group_begin, which returns iterator on the first element of the rest of the sequence.

    This is illustrative example:

    // recording sequence of subsequences
    packer.start_record(); 
    packer.next_group(some_iterator, sequence_size);
    /// ... apply it some times
    packer.end_record();

    // extracting sequence of subsequences
    for(packer_t::iterator it=packer.begin(), e=packer.end(); it!=e; it.plus()){
      int sequence_size = get_group_count(it);
      packer_t::iterator::base_t its=get_group_begin(it);
      // extracting subsequence
      for(int i=0;i<sequence_size;i++,++its){
        packer_t::value_t val=*its;
      }
    }

    Following packers are available:

    data_unpacked - just store data to some container (for example vector)
    group_unpacked - derived from data_unpacked for works with groups of data, stores size of each group

    int_pack - optimized packer for integers. 
    data_pack - optimized packer (for arbitrary data) which stores only unique data values and
      remember in which order dublicated data follows in order to extract them in the right way
      this packer has a container which stores only unique data 
      and additional packer of indices to indicate position of dublicated data in this container

    pair_pack - two data_packs
    
    group_pack - derived packer from data_pack optimized for packing of groups of data
    stores size of each group

    Packing data in a optimized way can save lots of memory but at the same time 
    it consumes extra time resources for recording / extracting data.
    Template pack_t allows to switch easily between packers which just collect data 
    (data_unpacked, group_unpacked) and packers which actually pack data (int_pack, data_pack, group_pack).
    Switching between these packers is regulated by macrodefinition USE_PACKERS.


    \ru @file seqpack.h
    Классы для упаковки данных. 
    Они отличаются друг от друга типом данных, для упаковки которых они предназначены, 
    способом их хранения, однако имеют общий интерфейс для работы.
    Начало записи осуществляется с помощью функции start_record.
    Далее последовательно вызывается функция record, аргументом которой является пакуемое значение. 
    По окончании записи вызывается end_record. Производить дозапись после вызова end_record нельзя.
    Однако можно начать запись сначала, предварительно вызвав функцию clear или start_record.

    К данным предоставляется последовательный доступ с помощью итераторов, 
    интерфейс работы с которыми похож на интерфейс работы с итераторами STL.
    Функции begin и end возвращают первый и следующий за последним итератор.
    У итераторов имеются операторы ++ (переход к следующему итератору), * (взятие значения),
    != (сравнение со следующим за последним итератором). Иллюстрацией служит следующий код, 
    в котором присходит перебор всех запакованных значений объекта packer класса packer_t:
    for(packer_t::iterator it=packer.begin(), e=packer.end(); it!=e; ++it)
      packer_t::value_t val=*it;

    В данный момент, похоже, у функций clear и start_record одинаковая реализация, и они дублируют друг друга...
**/

#include <vector>
#include <list>
#include <functional>
#include <map>
#include <algorithm>
#include <string>

#include "refobj.h"
#include "utiltl.h"
#include "logexc.h"  // temporary for messging -- remove this?


using namespace std;

/// packer iterator group categories
/// no groups defined
template <class packer_it>
struct nogroup_it{
  typedef packer_it group_it;
  static group_it get_group_begin(const packer_it &it){
    return it;
  }
  static int get_group_count(const packer_it &it){
    return 1;
  }
};

template< class packer_it>
struct packit_traits{
  typedef typename packer_it::gcategory::group_it  group_it;
  typedef typename packer_it::gcategory gcategory;
  typedef typename packer_it::index_type index_type;
};

/// default for all pointers
template< class T >
struct packit_traits< T* >{
  typedef T*         group_it;
  typedef nogroup_it<T *> gcategory;
  typedef ptrdiff_t index_type;
};


/// returns the number elements left in the sequence grouped with a given iterator
/// the element pointed by given iterator is also counted, so the default is 1
template <class packer_it>
inline typename packit_traits<packer_it>::index_type get_group_count(const packer_it &it){
  typedef typename packit_traits<packer_it>::gcategory cat_t;
  return cat_t::get_group_count(it);
}

/// returns the begining of the group associated with current iteratortemplate <class packer_it>
template <class packer_it>
inline typename packit_traits<packer_it>::group_it get_group_begin(const packer_it &it){
  typedef typename packit_traits<packer_it>::gcategory cat_t;
  return cat_t::get_group_begin(it);
}


/// extracts the next value packed in vector ind
/// at a given position defined by pos and incrcount
/// pos must be less than or equal ind.size()-2
/// @returns the value or -1 if end of sequence or the postion is incorrect
/// modifies pos (the next reading position) and incrcount
template<class index_t>
int get_next_record(const vector<index_t> &ind, index_t &pos, index_t &incrcount, index_t single=0);

/// records next value to vector ind
/// only nonnegative values are possible
/// @retuns    0 if the value was not packed (just stored)
///            1 packed with single increment
///            2 packed with other increment
///            3 packed as continued sequence
///           -1 tried to store after end-of sequence 
template<class index_t>
int put_next_record(vector<index_t> &ind, index_t cur, index_t single=0);

/// puts the end of record indicator 
template<class index_t>
int put_record_end(vector<index_t> &ind, index_t single=0);

/// this packer is optimized for packing integer sequences with the same increment.
/// for example, if this increment is zero, then sequence 2 2 2 2 3 3 3 3
/// will be effectively packed (there are more than one 2 and 3 which follows each other).
template <class index_tt = ptrdiff_t>
class index_pack{
  vector<index_tt> ind;
  /// sequences with this increment will be packed more efficiently
  /// for example, if increment 0 corresponds to the sequence 5 5 5 5 5,
  /// increment 2 corresponds to the sequence 5 7 9 11 13
  index_tt single;
  size_t count; /// number of packed integers
public:
  typedef index_tt index_t;
  typedef index_tt value_t;
  typedef index_tt record_t;

  index_pack(index_t sel_incr=0):single(sel_incr),count(0){
  }
  class iterator{
    friend class index_pack;
    index_t pos, incrcount; /// current position and incrcount in parent->ind (arguments of function get_next_record)
    index_t cur; /// currently extracted integer value
    const index_pack *parent;
    iterator(const index_pack *sparent, index_t end=0):parent(sparent), pos(0),incrcount(0){
      if(!end){
        cur=get_next_record(parent->ind,pos,incrcount,parent->single);
      }
      else{
        pos=(index_t)parent->ind.size()-1;
        cur=-1;
      }
    }
  public:
    typedef index_t index_type;
    typedef nogroup_it<iterator> gcategory;

    iterator(const iterator &other):pos(other.pos),incrcount(other.incrcount),cur(other.cur), parent(other.parent){
    }
    iterator() {} // for iterators array
    index_t operator*() const {
      return cur;
    }
    index_t operator++(){ // prefix
      return (cur=get_next_record(parent->ind,pos,incrcount,parent->single));
    }
    index_t operator++(int){ // postfix
      index_t tmp=cur;
      ++*this;
      return tmp;
    }
    /// checks for the end only
    bool operator!=(const iterator &other) const {
      if(pos!=other.pos)return true;
      else return false;
    }
  };
  int start_record(){
    ind.clear();
    count=0;
    return 1;
  }
  
 
  int next_record(index_t i){
    count++;
    return put_next_record(ind,i,single);
  }

   /// put a group of n elements starting with beg
  template <class int_it>
  index_t next_group(int_it beg, index_t n=1){
    index_t i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }

  int end_record(){
    return put_record_end(ind,single);
  }

  iterator begin() const {
    if(ind.size()<3){
      if(ind.size()>1)
        LOGERR(-1,"int_pack: requested iterator for unfinished sequence!",LINFO);
      return end();
    }
    return iterator(this,0);
  }
  iterator end() const {
    return iterator(this,1);
  }

  /// returns size in bytes ignoring vector overheads
  size_t packed_size() const {
    return sizeof(index_t)*ind.size();
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(index_t)*size();
  }

  size_t size() const {
    return count;
  }

  void clear(){
    count=0;
    ind.clear();
  }
};


typedef index_pack<int> int_pack;

template <class comp_pr>
struct subgroup_test{
  /// tests whether the begining of sequence pointed by beg1 coincides with the sequence 
  /// of n elements pointed by beg2
  /// @return the number of non-coinciding elements (returns 0 if beg2 is a subsequence) 
  template <class g1_it, class g2_it >
  ptrdiff_t operator()(g1_it beg1, g1_it end1, g2_it beg2, ptrdiff_t n, const comp_pr &pr) const {
    while(n>0 && beg1!=end1){
      if(!pr(*beg1,*beg2))break;
      ++beg1;
      ++beg2;
      n--;
    }
    return n;
  }
};

/// container which stores only unique values (or group of values) of some type T
/// by function inset_ind (insert_group).
/// if value is unique, it will be stored in container, otherwise it will not.
/// function insert_ind (insert_group) returns index stored data in container.
/// cointainer type could be vector or list.
/// equality of items is checked by comparison prediacate of type comp_pr
/// which must have bool opertor()(const T&,const T&) returning true for duplicated items.
/// cointainer_set is used as a container for stored data in data_pack
template<class T, class cont_t=vector<T>, class comp_pr=equal_to<T> >
class container_set: public cont_t{
protected:
  comp_pr pr;
public:
  typedef typename cont_t::iterator iterator;
  typedef comp_pr key_compare;

  /// returns the number of elements between begin and the inserted element
  size_t insert_ind(iterator it,const T& elm){
    iterator e=cont_t::end(); 
    size_t i=0;
    for(;it!=e;++it){
      if(pr(*it,elm))return i;
      i++;
    }
    this->push_back(elm);
    return i; 
  }


  /// the same with default starting
  size_t insert_ind(const T& elm){
    return insert_ind(cont_t::begin(),elm);
  }

  template<class inp_it, class group_pr>
  size_t insert_group(iterator it, inp_it beg, ptrdiff_t n, const group_pr & gpr){
    iterator e=cont_t::end(); 
    size_t i=0;
    for(;it!=e;++it){
      if(gpr(it,e,beg,n,pr)==0){ // found subsequence which is exactly n elements long
        return i;
      }
      i++;
    }
    // not found: inserting
    while(n>0){
      push_back(*beg);
      ++beg;
      n--;
    }
    return i;
  }

  template<class inp_it, class group_pr>
  size_t insert_group(inp_it beg, ptrdiff_t n, const group_pr & gpr){
    return insert_group(cont_t::begin(), beg, n, gpr);
  }

  void clear_index(){}
};


template<class T, class comp_pr=equal_to<T> > 
class vector_set: public container_set<T,vector<T>, comp_pr >{
public:
  typedef typename container_set<T,vector<T>, comp_pr >::iterator iterator;
  typedef typename container_set<T,vector<T>, comp_pr >::key_compare key_compare;
};



template<class T, class comp_pr=equal_to<T> > 
class list_set: public container_set<T,list<T>, comp_pr >{
public:
  typedef typename container_set<T,list<T>, comp_pr >::iterator iterator;
  typedef typename container_set<T,list<T>, comp_pr >::key_compare key_compare;
};

/// the functionality is the same as for container_set.
/// however, checking if element is dublicated is done by map which store
/// indices in container_set::container as keys and unique values (groups of values) as values
template<class T, class less_pr=less<T>, class index_t = ptrdiff_t > 
class vector_indexed: public container_set<T,vector<T>, equal_from_less_pr<T,less_pr> >{
public:
  typedef container_set<T,vector<T>, equal_from_less_pr<T,less_pr> > base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::key_compare key_compare;
protected:
  typedef multimap<T, size_t, less_pr> map_t;
  map_t imap;
  typedef typename map_t::iterator map_it;
  typedef typename map_t::value_type map_vt;

public:

  template<class inp_it, class group_pr>
  size_t insert_group(typename base_t::iterator it, inp_it beg, index_t n, const group_pr & gpr){
    size_t i=0;
    map_it mit=imap.find(*beg);
    // check the group in any case
    iterator e=this->end(); 
    while(it!=e){
      if(mit!=imap.end() && mit->first==*beg){
        it=base_t::begin()+mit->second;
        i=mit->second;
        ++mit;
      }
      if(gpr(it,e,beg,n,this->pr)==0){ // found subsequence which is exactly n elements long
        return i;
      }
      i++;
      ++it;
    }
    index_t i0=(index_t)i;
    // not found: inserting
    while(n>0){
      push_back(*beg);
      imap.insert(map_vt(*beg,i0));
      ++beg;
      n--;
      i0++;
    }
    return i;
  }

  template<class inp_it, class group_pr>
  size_t insert_group(inp_it beg, index_t n, const group_pr & gpr){
    return insert_group(this->begin(), beg, n, gpr);
  }


  /// insert with it starting
  size_t insert_ind(typename base_t::iterator it, const T& elm){
    map_it mit=imap.find(elm);
    if(mit!=imap.end())return mit->second;
    // inserting
    imap.insert(map_vt(elm,base_t::size()));
    this->push_back(elm);
    return base_t::size()-1;
  }

  /// the same with default starting
  size_t insert_ind(const T& elm){
    return insert_ind(base_t::begin(),elm);
  }

  void clear_index(){
    imap.clear();
  }

  void clear(){
    clear_index();
    base_t::clear();
  }

};


/// Stores data in some container of the the type cont_t 
template<class T, class cont_t=vector<T>, class index_t = ptrdiff_t >
class data_unpacked{
public:
  typedef typename cont_t::iterator iterator;
  typedef typename cont_t::iterator vset_it;
  typedef T value_t;
  typedef T record_t;

protected: 
  cont_t vset;
public:
  
  data_unpacked(){}

  int start_record(){
    vset.clear();
    return 1;
  }
  int next_record(const T& value){
    vset.push_back(value);
    return 1;
  }

   /// put a group of n elements starting with beg
  template <class rec_it>
  index_t next_group(rec_it beg, index_t n=1){
    index_t i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }

  int end_record(){
    return 1;
  }

  T &operator[] (size_t i) {
    return vset[i];
  }
  
  iterator begin(){
    return vset.begin();
  }
  
  iterator end(){
    return vset.end();
  }

  size_t size() const {
    return vset.size();
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return data_size();
  }
  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size();
  }

  void clear(){
    vset.clear();
  }
};

/// Class for packing data of the type T.
/// Data is placed to some container set_t which stores only unique values.
/// Iterative access to packed data is regulated by packer of integers int_pack
/// which stores indices of packed elements in container set_t.
/// set_t can be container_set<T, container_t> (vector_set or list_set).
/// The comparison predicate is of type comp_pr
/// which must have bool operator()(const T&,const T&) returning true for duplicated items
template<class T, class comp_pr=equal_to<T>, class set_t=vector_set<T, comp_pr>, class index_t = ptrdiff_t >
class data_pack{
public:
  typedef typename set_t::iterator vset_it;
  typedef T value_t;
  typedef T record_t;

protected:
  set_t vset; /// container for unique values of packed data
  index_pack<index_t> ipack; /// packer for indices of stored elements in vset

public:
  class iterator{
  protected:
    friend class data_pack;
    typename index_pack<index_t>::iterator it;
    data_pack *parent;
    iterator(typename index_pack<index_t>::iterator sit, data_pack *sparent): it(sit), parent(sparent){}
  public:  
    typedef nogroup_it<iterator> gcategory;
    typedef index_t index_type;
    /// copy constructor
    iterator(const iterator &other):it(other.it), parent(other.parent){}

    iterator() {} /// empty constructor to keep the possibility specify iterators array
    
    /// WARNING: recording by the given reference will  effectively
    ///          change ALL data entries that are equal to this one 
    T& operator*() const {
      return parent->vset[*it];
    }
    // prefix increment
    iterator& operator++(){
      ++it;
      return *this;
    } 
    // postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++it;
      return tmp;
    } 

    bool operator!=(const iterator &other){
      return (it!=other.it);
    } 

  };
  data_pack(): vset(), ipack(0) {}

  int start_record(){
    vset.clear();
    return ipack.start_record();
  }
  int next_record(const T& value){
    size_t res=vset.insert_ind(value);
    return ipack.next_record((index_t)res);
  }

  /// put a group of n elements starting with beg
  template <class rec_it>
  index_t next_group(rec_it beg, index_t n=1){
    index_t i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }

  int end_record(){
    vset.clear_index();
    return ipack.end_record();
  }
  
  iterator begin() {
    return iterator(ipack.begin(),this);
  }
  
  iterator end()  {
    return iterator(ipack.end(),this);
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return sizeof(value_t)*vset.size()+ipack.packed_size();
  }

  size_t size() const {
    return ipack.size();
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size();
  }
  void clear(){
    vset.clear();
    ipack.clear();
  }
};

/// Class for packing pairs of data of two different types.
/// The types are packed by different packers but recovered simultaneously by an iterator
/// pair_tt is the returned value type, it must have at least one of the following constructors
///          pair_tt(value_t1 &a,value_t2 &b) 
///          pair_tt(const value_t1 &a,const value_t2 &b)
///          pair_tt(value_t1 a,value_t2 b) 
/// which may construct a type passing arguments either by value or by reference. 
/// WARNING: recording by the returned references is possible for some packers, but may lead to changes 
///          in other packed elements,
///          the policy of recording is determined by packers
/// By default the returned value type of pair_pack is pair<value_t1,value_t2> returned by value (no recording possible)
/// see refpair_pack()
template<class T1, class T2, class packer_t1=data_pack<T1>, class packer_t2=data_pack<T2>, 
         class pair_tt=pair<typename packer_t1::value_t,typename packer_t2::value_t >, class index_t = ptrdiff_t >
class pair_pack{
public: 
  packer_t1 *pack1;
  packer_t2 *pack2;
 

  typedef typename packer_t1::record_t record_t1;
  typedef typename packer_t2::record_t record_t2;
  typedef pair_tt value_t;
  typedef pair<record_t1,record_t2> record_t;

  /// constructor with possible packer construction specification
  pair_pack(packer_t1 *sp1=new packer_t1(), packer_t2 *sp2 = new packer_t2()):pack1(sp1),pack2(sp2){}

  class iterator{
    friend class pair_pack;
    typename packer_t1::iterator it1;
    typename packer_t2::iterator it2;
    iterator(typename packer_t1::iterator sit1, typename packer_t2::iterator sit2): it1(sit1),it2(sit2){}
  public:  
    /// the class gcategory must have a function
    /// static int get_group_count(const iterator &) which returns group count
    typedef iterator gcategory;
    typedef index_t index_type;

    static index_t get_group_count(const iterator &it){ 
      index_t n1=get_group_count(it.it1);
      index_t n2=get_group_count(it.it2);
      return (n1>n2 ? n1 : n2);
    }

    /// copy constructor
    iterator(const iterator &other):it1(other.it1),it2(other.it2){}
    iterator() {} // for iterators array
    
    value_t operator*() const {
      return value_t(*it1,*it2);
    }
    /// prefix increment
    iterator& operator++(){
      ++it1;
      ++it2;
      return *this;
    } 
    /// postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    } 
    bool operator!=(const iterator &other){
      return (it1!=other.it1);
    }
  };
  
  int start_record(){
    int res=pack1->start_record();
    if(res<0)return res;
    res=pack2->start_record();
    return res;
  }
  
  int next_record(const record_t1 &v1, const record_t2 &v2){
    int res=pack1->next_record(v1);
    if(res<0)return res;
    res=pack2->next_record(v2);
    return res;
  } 
  
  int next_record(const record_t &p){
    return next_record(p.first,p.second);
  }

  /// put a group of n elements starting with beg
  template <class rec_it>
  index_t next_group(rec_it beg, index_t n=1){
    int res=pack1->next_group(first_it<rec_it>(beg), n);
    if(res<0)return res;
    res=pack2->next_group(second_it<rec_it>(beg), n);
    return res;
  }
  
  int end_record(){
    int res=pack1->end_record();
    if(res<0)return res;
    res=pack2->end_record();
    return res;
  }

  iterator begin(){
    return iterator(pack1->begin(),pack2->begin());
  }
  
  iterator end(){
    return iterator(pack1->end(), pack2->end());
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return pack1->packed_size()+pack2->packed_size();
  }

  size_t size() const {
    return pack1->size(); // must be the same for both
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return pack1->data_size()+pack2->data_size();
  }

  void clear(){
    pack1->clear();
    pack2->clear();
  }

  ~pair_pack(){
    delete pack1;
    delete pack2;
  }
};

/// both types stored by reference
template <class T1, class T2>
struct refpair{
  T1& first;
  T2& second;
  typedef T1 first_type;
  typedef T2 second_type;
  refpair(T1 &a, T2 &b):first(a),second(b){}

  operator pair<T1,T2>(){
    return make_pair(first,second);
  }

  T1 &r1() const { return first; }
  T2 &r2() const { return second; }
  refpair(const refpair<T1,T2> &other): first(other.first), second(other.second) {}
};

template<class T1, class T2, class packer_t1=data_pack<T1>, class packer_t2=data_pack<T2> > 
class refpair_pack: public pair_pack<T1,T2,packer_t1,packer_t2, refpair<typename packer_t1::value_t,typename packer_t2::value_t> >{
public:
  typedef pair_pack<T1,T2,packer_t1,packer_t2, refpair<typename packer_t1::value_t,typename packer_t2::value_t> > base_t;
  typedef typename base_t::record_t1 record_t1;
  typedef typename base_t::record_t2 record_t2;
  typedef typename base_t::value_t value_t;
  typedef typename base_t::record_t record_t;
  typedef typename base_t::iterator iterator;

  /// constructor with possible packer construction specification
  refpair_pack(packer_t1 *sp1=new packer_t1(), packer_t2 *sp2 = new packer_t2()):base_t(sp1,sp2){}

};


/// Packer optimized to pack groups of data.
/// set_t can be container_set (vector_set or list_set).
/// The comparison predicate is of type comp_pr
/// which must have bool opertor()(const T&,const T&) returning true for duplicated sequences
template<class T, class comp_pr=equal_to<T>, class set_t=vector_set<T, comp_pr>, class index_t = ptrdiff_t , class group_pr=subgroup_test<typename set_t::key_compare> >
class group_pack: public data_pack<T,comp_pr,set_t>{
public:
  typedef data_pack<T,comp_pr,set_t, index_t> base_t;
  typedef typename base_t::vset_it vset_it;
  typedef T value_t;
  typedef T record_t;
  typedef group_pr key_compare;

protected: 
  using  base_t::vset;
  using  base_t::ipack;

  /// packer of number of elements in each packed group.
  /// this reference can be shared with some other group_pack
  mngptr<index_pack<index_t> > npack;
  group_pr pr;
  index_t sz; /// total number of packed elements
public:

  class iterator: public base_t::iterator {
    typedef typename base_t::iterator base_it;
    friend class group_pack;
    typename index_pack<index_t>::iterator itn; /// iterator of npack
    index_t n; /// number of already iterated elements inside current group
    using base_it::it; /// actual iterator on a packed data
    iterator(typename index_pack<index_t>::iterator siti, group_pack *sparent, typename index_pack<index_t>::iterator sitn): base_it(siti,sparent),itn(sitn),n(0){}
  public:  
    /// the class gcategory must have a function
    /// static index_t get_group_count(const iterator &) which returns group count
    typedef iterator gcategory;
    typedef vset_it group_it;
    typedef index_t index_type;

    static index_t get_group_count(const iterator &it){  /// number of the rest elements in current group
      return *(it.itn)-it.n;
    }

    static group_it get_group_begin(const iterator &it){ 
/* old version that worked OK: 6 hours of MPI debugging to catch this bug 
      if (it.it!=it.parent->ipack.end()) {
        index_t sh=*(it.it);
        group_pack *ig=(group_pack *)it.parent;
        return (ig)->vset.begin()+sh+it.n;
      }
      return ((group_pack *)it.parent)->vset.end(); */
       
      index_t sh=*(it.it); /// index iterator on data_pack
      if(sh<0) /// if groups with 0 elements only are packed after
        return ((group_pack *)it.parent)->vset.end();
      group_pack *ig=(group_pack *)it.parent;
      return (ig)->vset.begin()+sh+it.n;
    }

    /// copy constructor
    iterator(const iterator &other):base_it(other),itn(other.itn),n(other.n){}
    iterator() {} // for iterators array

    /// WARNING: recording by the given reference will  effectively
    ///          change ALL data entries that are equal to this one 
    T& operator*() const {
      return ((group_pack *)this->parent)->vset[*it+n];
    }
    // prefix increment
    iterator& operator++(){
      n++;
      if(n>=*itn){ // safe also after end of sequence, returns -1
        if(*itn!=0) // not zero group
          ++it;
        ++itn;
        n=0;
      }
      return *this;
    } 
    // postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    } 

    iterator& plus(){
      if(*itn)
        ++it;
      ++itn;
      n=0;
      return *this;
    }

    bool operator!=(const iterator &other) const {
//      return ((it!=other.it)&&(itn!=other.itn)&&(n!=other.n));
      return (it!=other.it);
    } 

  };
  
  group_pack():sz(0),npack(new index_pack<index_t>(1),1){}

  int start_record(){
    sz=0;
    if(npack.managed()){ // if this is not dependent service packer
      npack->start_record();
    }
    return base_t::start_record();
  }

  int next_record(const T& value){
    return next_group(&value,1);
  }

  /// this allows to share auxiliary data with other packer
  /// THE PACKERS MUST BE SYNCRONIZED EXTERNALLY: next_group must be called in sync with parent packer
  /// USE WITH CARE
  template <class packer_t>
  int aux_depends_on(const packer_t &other){
    npack.reset(other.get_grouper(),0);
    return 1;
  }

  index_pack<index_t> *get_grouper() const {
    return npack.ptr();
  }

  /// put a group of n elements starting from beg
  /// packing the group with 0 elements will lead to the following behaviour when unpacked:
  /// get_group_count(it) will return 0, *it will return the next recorded element (if it exists)
  /// then after ++it, *it will return the same element, but get_group_count(it) will return the next group count
  /// THE ELEMENTS POINTED BY BEG WILL NOT BE RECORDED WHEN n=0
  /// DO NOT USE *it if get_group_count(it) returns 0!
  template <class rec_it>
  index_t next_group(rec_it beg, index_t n=1){
    if(npack.managed()){ // if this is not dependent service packer
      if(npack->next_record(n)<0)return -1;
    }
    if(n==0)
      return 0;
    size_t res=vset.insert_group(beg,n,pr);
    ipack.next_record((index_t)res);
    sz+=n;
    return n;
  }

  int end_record(){
    if(npack.managed()){ // if this is not dependent service packer
      npack->end_record();
    }
    return base_t::end_record();
  }
  
  iterator begin() {
    return iterator(ipack.begin(),this, npack->begin());
  }
  
  iterator end()  {
    return iterator(ipack.end(),this, npack->end());
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return sizeof base_t::packed_size()+(npack.managed() ? npack->packed_size() :0);
  }

  size_t size() const {
    return sz;
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size()+(npack.managed() ? npack->data_size() :0);
//    return (npack.managed() ? npack->data_size() :0);
  }

  void clear(){
    npack->clear();
    sz=0;
    return base_t::clear();
  }

};

/// Stores groups of data in some container of the the type cont_t 
/// Number of elements in a group is packed in npack.
template<class T, class cont_t=vector<T>, class index_t = ptrdiff_t >
class group_unpacked: public data_unpacked<T,cont_t>{
public:
  typedef data_unpacked<T> base_t;
  typedef typename base_t::vset_it vset_it;
  typedef T value_t;
  typedef T record_t;

protected: 
  using  base_t::vset;

  /// packer of indices of last elements in cont_t in each packed group (not like in group_pack).
  /// this reference can be shared with some other group_unpacked
  mngptr<data_unpacked<index_t> > npack;
  index_t sz;
public:

  class iterator: public base_t::iterator {
    typedef typename base_t::iterator base_it;
    friend class group_unpacked;
    index_t cur_ind; /// shift between base_t::iterator at the beginning of the group and base_t::vset.begin()
    typename data_unpacked<index_t>::iterator itn; /// iterator of npack
    index_t n; /// number of already iterated elements inside current group
    iterator(base_it siti, typename data_unpacked<index_t>::iterator sitn, index_t cur_ind_): base_it(siti),itn(sitn),cur_ind(cur_ind_),n(0){}
  public:  
    /// the class gcategory must have a function
    /// static index_t get_group_count(const iterator &) which returns group count
    typedef iterator gcategory;
    typedef base_it group_it;
    typedef index_t index_type;

    static index_t get_group_count(const iterator &it){ 
      return *(it.itn)-it.cur_ind-it.n;
    }

    static group_it get_group_begin(const iterator &it){ 
      return base_it(it);
    }

    iterator() {} // for iterators array

    iterator& operator++(){
      n++;
      if(n>=*itn-cur_ind){ // safe also after end of sequence, returns -1
        if((*itn-cur_ind)!=0)
          base_it::operator++();
        cur_ind=*itn;
        ++itn;
        n=0;
      }
      else
        base_it::operator++();
      return *this;
    } 
    // postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    }
    iterator& plus(){
      base_it::operator+=(*itn-cur_ind-n);
      cur_ind=*itn;
      ++itn;
      return *this;
    }

    bool operator!=(const iterator &other) const {
      return itn!=other.itn;
    } 

  };
  
  group_unpacked():sz(0),npack(new data_unpacked<index_t>(),1){}

  int start_record(){
    sz=0;
    if(npack.managed()){ // if this is not dependent service packer
      npack->start_record();
    }
    return base_t::start_record();
  }

  int next_record(const T& value){
    return next_group(&value,1);
  }

  /// this allows to share auxiliary data with other packer
  /// THE PACKERS MUST BE SYNCRONIZED EXTERNALLY: next_group must be called in sync with parent packer
  /// USE WITH CARE
  template <class packer_t>
  index_t aux_depends_on(const packer_t &other){
    npack.reset(other.get_grouper(),0);
    return 1;
  }

  data_unpacked<index_t> *get_grouper() const {
    return npack.ptr();
  }

  /// put a group of n elements starting from beg
  /// packing the group with 0 elements will lead to the following behaviour when unpacked:
  /// get_group_count(it) will return 0, *it will return the next recorded element (if it exists)
  /// then after ++it, *it will return the same element, but get_group_count(it) will return the next group count
  /// THE ELEMENTS POINTED BY BEG WILL NOT BE RECORDED WHEN n=0
  /// DO NOT USE *it if get_group_count(it) returns 0!
  template <class rec_it>
  index_t next_group(rec_it beg, index_t n=1){
    if(npack.managed()){ // if this is not dependent service packer
      index_t cur_n=0;
      size_t size=npack->size();
      if(size)
        cur_n+=(*npack)[size-1];
      if(npack->next_record(cur_n+n)<0)return -1;
    }
    if(n==0)return 0;
    base_t::next_group(beg,n);
    sz+=n;
    return n;
  }

  int end_record(){
    if(npack.managed()){ // if this is not dependent service packer
      npack->end_record();
    }
    return base_t::end_record();
  }
  
  iterator operator[] (size_t i) {
    if(i>=npack->size())
      return end();
    typename data_unpacked<index_t>::iterator itn=npack->begin()+i;
    index_t cur_ind = i ? (*npack)[i-1] : 0;
    return iterator(base_t::begin()+cur_ind,itn,cur_ind);
  }

  iterator begin() {
    return iterator(base_t::begin(), npack->begin(),0);
  }
  
  iterator end() {
    return iterator(base_t::end(), npack->end(),-1);
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return data_size();
  }

  size_t size() const {
    return sz;
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return base_t::data_size()+(npack.managed() ? npack->data_size():0);
  }

  void clear(){
    npack->clear();
    sz=0;
    return base_t::clear();
  }
};


template<class pair_t>
struct pair_sort_t {
  bool operator()(const pair_t& a, const pair_t& b) const {
    return a.first<b.first;
  }
};

/// reorders two packs of the same size in memory with < operator, the first is pack the leading, 
/// the second follows the order of the first
template<class pack_t1, class pack_t2>
int mem_reorder(pack_t1 &p1, pack_t2 &p2){
  size_t n=p1.size(), i;
  if(n!=p2.size())return -1;
  typename pack_t1::iterator it1=p1.begin(), e1=p1.end();
  typename pack_t2::iterator it2=p2.begin(), e2=p2.end();
  typedef typename pack_t1::value_t value_t1;
  typedef typename pack_t2::value_t value_t2;
  typedef pair<value_t1, value_t2> sort_t;

  pair_sort_t<sort_t> pred;
  
  // copying
  sort_t *v= new sort_t[n];
  for(i=0;i<n;i++){
    v[i]=sort_t(*it1,*it2);
    ++it1;
    ++it2;
  }
  // sorting
  sort(v,v+n,pred);
  // copying the sorted records
  p1.clear();
  p1.start_record();
  p2.clear();
  p2.start_record();
  for(i=0;i<n;i++){
    p1.next_record(v[i].first);
    p2.next_record(v[i].second);
  }
  p1.end_record();
  p2.end_record();
  delete [] v;
  return 1;
}

//#define USE_PACKERS

/// this template helps to easily switch between packing and not packing
template<class T, class index_t = ptrdiff_t>
struct pack_t{

#ifndef USE_PACKERS

  typedef data_unpacked<T, std::vector<T>, index_t> data;
  typedef group_unpacked<T, std::vector<T>, index_t> group;

#else

  template<class T>
  struct data_pack_t{
    typedef data_pack<T, index_t> type;
  };

  template<>
  struct data_pack_t<int>{
    typedef index_pack<int> type;
  };

  template<>
  struct data_pack_t<ptrdiff_t>{
    typedef index_pack<ptrdiff_t> type;
  };

  typedef typename data_pack_t<T>::type data;
//  typedef group_pack<T> group;
  typedef group_pack<T, equal_to<T> , vector_indexed<T>, index_t > group;

#endif
};

#endif
