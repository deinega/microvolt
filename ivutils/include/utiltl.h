# ifndef _UTILTL
# define _UTILTL

/*e @file utiltl.h
    @brief Useful templates (iterators, modifiers, etc)
  
**/

#include "stdarg.h"
#include <vector>
#include <iterator>
#include <functional>

using namespace std;

template<class iter>
struct value_type_cast{
  typedef typename iter::value_type value_type;
};

template<class T>
struct value_type_cast<T *>{
  typedef T value_type;
};

///\en make vector with optional amount of elements
template<class T>
vector<T> makevec(int num, ...) {
  vector<T> res(num);
  va_list(args);
  va_start(args, num);
  for(int i=0;i<num;i++)
    res[i]=va_arg(args,T);
  va_end(args);
  return res;
}

///\en make vector with one element
template<class T>
vector<T> makevec_1(T val) {
  vector<T> res;
  res.push_back(val);
  return res;
}

///\en make vector with two elements
template<class T>
vector<T> makevec_2(T val1, T val2) {
  vector<T> res;
  res.push_back(val1);
  res.push_back(val2);
  return res;
}

/// two values are assumed to be equal, if none of them is less than other
template<class T, class less_pr>
struct equal_from_less_pr: public binary_function<T,T,bool>{
  less_pr pr;
  bool operator()(const T& r1, const T& r2) const {
    return (!pr(r1,r2))&&(!pr(r2,r1));
  }
};

struct void_check_t{
  bool operator()(int i) const { return true; }
};

///\en aggregate iterator (two iterators are aggregated in one)
template<class it1, class it2=it1, class value_t= typename value_type_cast<it1>::value_type>
class aggregate_it{
  pair<it1,it1> p1;
  it2 i2;
public:
  aggregate_it(it1 beg1, it1 end1, it2 beg2):p1(beg1,end1),i2(beg2){}

  aggregate_it(const aggregate_it &other):p1(other.p1),i2(other.i2){}

  typedef value_t value_type;

  const value_type operator*() const {
    if(p1.first!=p1.second)return *p1.first;
    else return *i2;
  }
  
  aggregate_it& operator++(){ // prefix
    if(p1.first!=p1.second)p1.first++;
    else i2++;
    return *this;
  }

  aggregate_it operator++(int){ // postfix
    aggregate_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const aggregate_it &other) const{
    return (p1.first!=other.p1.first) || (i2!=other.i2);
  }
};

///\en iterator removing some subsequence (beg2,end2) from the basic sequence (beg1,end1) if comparison predicate returns false
template<class iter1, class iter2=iter1, class pred_t=not_equal_to<iter1> >
class remove_if_it{
  iter1 beg1, end1;
  iter2 beg2, end2;
  pred_t Pr;
  
  /// move beg1 through sequence (beg2,end2), if beg1 doesn't belong there
  void check_beg1(){
    iter2 it2=beg2;
    for(;beg1!=end1 && it2!=end2;++it2){
      if(!Pr(beg1,it2))++beg1;
    }
  }

public:
  remove_if_it(iter1 sbeg1, iter1 send1, iter2 sbeg2, iter2 send2):beg1(sbeg1),beg2(sbeg2),end1(send1),end2(send2){
    check_beg1();
  }
  remove_if_it(const remove_if_it &other):beg1(other.beg1),beg2(other.beg2),end1(other.end1),end2(other.end2){
    check_beg1();
  }

  typedef typename value_type_cast<iter1>::value_type value_type;

  const value_type operator*() const {
    return *beg1;
  }
  
  remove_if_it& operator++(){ // prefix
    ++beg1;
    check_beg1();    
    return *this;
  }

  remove_if_it operator++(int){ // postfix
    remove_if_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const remove_if_it &other) const{
    return (beg1!=other.beg1);
  }
};

///\en  Skips some of the numbers, returned as false by checker
///     Checker MUST return true for out_of range indicies.
template<class iterator, class checker_t=void_check_t>
struct skiper_it{
  checker_t checker;
  iterator it;
  int i;
  
  typedef typename iterator_traits<iterator>::iterator_category iterator_category;
	typedef typename iterator_traits<iterator>::value_type value_type;
	typedef typename iterator_traits<iterator>::difference_type difference_type;
	typedef difference_type distance_type;	// retained
	typedef typename iterator_traits<iterator>::pointer pointer;
	typedef typename iterator_traits<iterator>::reference reference;

  skiper_it(){}

  ///\en Use \a valid=flase to indicate invalid iterators, 
  ///    where no checks are performed (end of sequence, etc.)
  skiper_it(iterator it_, int i_=0, bool valid=true):it(it_),i(i_){
    if(valid){
      while(!checker(i)){
        ++i;
        ++it;
      }
    }
  }

  ///\en Use \a valid=flase to indicate invalid iterators, 
  ///    where no checks are performed (end of sequence, etc.)
  skiper_it(iterator it_, const checker_t &ch, int i_=0, bool valid=true):it(it_),i(i_),checker(ch){
    if(valid){
      while(!checker(i)){
        ++i;
        ++it;
      }
    }
  }
  
  
  value_type &operator*(){
    return *it;
  }
  skiper_it &operator++(){
    do{
      ++i;
      ++it;
    }while(!checker(i));
    return *this;
  }
  
  skiper_it operator++(int){ // postfix
    skiper_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator==(const skiper_it &other) const {
    return (it==other.it);
  }
  bool operator!=(const skiper_it &other) const {
    return !(*this==other);
  }

};



///\en looper: loops the sequence after nd iterations, counts the loop number (iter)
template<class iterator, class checker_t=void_check_t>
struct loop_it{
  size_t nd, i, iter;
  iterator it, beg;
  

  typedef typename iterator_traits<iterator>::iterator_category iterator_category;
	typedef typename iterator_traits<iterator>::value_type value_type;
	typedef typename iterator_traits<iterator>::difference_type difference_type;
	typedef difference_type distance_type;	// retained
	typedef typename iterator_traits<iterator>::pointer pointer;
	typedef typename iterator_traits<iterator>::reference reference;


  loop_it(){}
  loop_it(iterator it_, size_t nd_, size_t iter_=0):it(it_), beg(it_), i(0), iter(iter_), nd(nd_){}
  
  /*dup_it(const dup_it &operator other):it(other.it),i(other.i), nd(other.nd){}
  dup_it& operator=(const dup_it &operator other){
    if(this!=&other)return *this=
  }*/

  value_type &operator*(){
    return *it;
  }
  loop_it &operator++(){
    i++;
    if(i<nd)++it;
    else{
      i=0;
      it=beg;
      iter++;
    }
    return *this;
  }
  loop_it operator++(int){ // postfix
    loop_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator==(const loop_it &other) const {
    return (i==other.i && iter==other.iter);
  }
  bool operator!=(const loop_it &other) const {
    return !(*this==other);
  }

};


///\en duplicator: duplicates each entry in the sequence nd times
template<class iterator, class checker_t=void_check_t>
struct dup_it{
  size_t nd, i;
  iterator it;

  typedef typename iterator_traits<iterator>::iterator_category iterator_category;
  typedef typename iterator_traits<iterator>::value_type value_type;
  typedef typename iterator_traits<iterator>::difference_type difference_type;
  typedef difference_type distance_type;	// retained
  typedef typename iterator_traits<iterator>::pointer pointer;
  typedef typename iterator_traits<iterator>::reference reference;

  dup_it(){}
  dup_it(iterator it_, size_t nd_, size_t i_=0):it(it_), i(i_), nd(nd_){}
  
  
  value_type &operator*(){
    return *it;
  }
  dup_it &operator++(){
    i++;
    if(i>=nd){
      ++it;
      i=0;
    }
    return *this;
  }
  dup_it operator++(int){ // postfix
    dup_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator==(const dup_it &other) const {
    return (i==other.i && it==other.it);
  }
  bool operator!=(const dup_it &other) const {
    return !(*this==other);
  }

};


///\en returns the first pair value of iterate pairs
template<class pair_it>
class first_it:public pair_it{
public:
  typedef typename iterator_traits<pair_it>::iterator_category iterator_category;
  typedef typename iterator_traits<pair_it>::value_type::first_type value_type;
  typedef typename iterator_traits<pair_it>::difference_type difference_type;
  typedef difference_type distance_type;	// retained
  typedef typename iterator_traits<pair_it>::pointer pointer;
  typedef typename iterator_traits<pair_it>::reference reference;

  first_it(){}
  
  first_it(const pair_it &base):pair_it(base){}
  
  value_type &operator*(){
    return ((pair_it *)this)->first;
  }  
};

///\en returns the first pair value of iterate pairs
template<class pair_it>
class second_it:public pair_it{
public:
  typedef typename iterator_traits<pair_it>::iterator_category iterator_category;
  typedef typename iterator_traits<pair_it>::value_type::second_type value_type;
  typedef typename iterator_traits<pair_it>::difference_type difference_type;
  typedef difference_type distance_type;	// retained
  typedef typename iterator_traits<pair_it>::pointer pointer;
  typedef typename iterator_traits<pair_it>::reference reference;

  second_it(){}
  
  second_it(const pair_it &base):pair_it(base){}
  
  value_type &operator*(){
    return ((pair_it *)this)->second;
  }  
};


template <class T>
struct plusplus_t{
  typedef T result_type;
  typedef T argument_type;
  T &operator()(T& v) const { v++; return v; }
};

template <class T>
struct minusminus_t{
  typedef T result_type;
  typedef T argument_type;
  T &operator()(T& v) const { v--; return v; }
};

template <class T>
struct noop_t{
  typedef T result_type;
  typedef T argument_type;
  T &operator()(T& v) const { return v; }
};

template <class T>
struct plus_t{
  T value;
  typedef T result_type;
  typedef T argument_type;
  plus_t(const T& value_=T()):value(value_){}
  T &operator()(T& v) const { v+=value; return v; }
};

///\en incremental iterator (to use together with index sets)
/// increment/decrement operators are applied to the value
/// returns the value itself when calling *it
template<class T, class incr_op=plusplus_t<T>, class decr_op=minusminus_t<T> >
class value_it{
protected:
  T v;
  incr_op incr;
  decr_op decr;
public:
  typedef random_access_iterator_tag iterator_category;
  typedef T value_type;
  typedef T difference_type;
  typedef difference_type distance_type;	// retained
  typedef T* pointer;
  typedef T& reference;

  value_it(){}
  value_it(T v_):v(v_){}
  
  T &operator*(){
    return v;
  }
  
  value_it &operator++(){
    //v++;
    incr(v);
    return *this;
  }

  value_it operator++(int){ // postfix
    value_it tmp=*this;
    ++*this;
    return tmp;
  }

  value_it &operator--(){
    //v--;
    decr(v);
    return *this;
  }

  // difference
  T operator-(const value_it & other) const {
    return v-other.v;
  }

  value_it operator--(int){ // postfix
    value_it tmp=*this;
    --*this;
    return tmp;
  }

  bool operator==(const value_it &other) const {
    return (v==other.v);
  }
  bool operator!=(const value_it &other) const {
    return !(*this==other);
  }

};

///\en forward-iterating a random acces iterator via index set
/// increment, comparison operators  are applied to indexer
/// *it returns the value from rand_it by incrementing begin with the current index
template<class rand_it, class indexer_it>
class indexed_it{
  rand_it beg;
  indexer_it ind;
public:
  typedef forward_iterator_tag iterator_category;
  typedef typename iterator_traits<rand_it>::value_type value_type;
  typedef typename iterator_traits<rand_it>::difference_type difference_type;
//  typedef typename iterator_traits<rand_it>::distance_type distance_type;	// retained
  typedef typename iterator_traits<rand_it>::pointer pointer;
  typedef typename iterator_traits<rand_it>::reference reference;

  indexed_it(){}
  indexed_it(rand_it beg_, indexer_it ind_):beg(beg_),ind(ind_){}

  
  value_type operator*() const {
    return *(beg+*ind);
  }
  
  indexed_it &operator++(){
    ind++;
    return *this;
  }

  indexed_it operator++(int){ // postfix
    indexed_it tmp=*this;
    ++*this;
    return tmp;
  }

 
  bool operator==(const indexed_it &other) const {
    return (ind==other.ind);
  }
  bool operator!=(const indexed_it &other) const {
    return !(*this==other);
  }
 
};


///\en Iterator modifying the returned value by an operation
template<class base_it, class operation_t = noop_t<typename iterator_traits<base_it>::value_type > >
class modified_value_it: public base_it {
  operation_t op;
public:
  typedef typename iterator_traits<base_it>::iterator_category iterator_category;
  typedef typename operation_t::result_type value_type;
  typedef typename iterator_traits<base_it>::difference_type difference_type;
//  typedef typename iterator_traits<base_it>::distance_type distance_type;	// retained
  typedef typename iterator_traits<base_it>::pointer pointer;
  typedef typename iterator_traits<base_it>::reference reference;

  modified_value_it(){}
  modified_value_it(const base_it &it):base_it(it){}
  modified_value_it(const base_it &it, const operation_t op_):base_it(it), op(op_){}

  
  value_type operator*() const {
    typename iterator_traits<base_it>::value_type v= *(*((base_it *)this));
    return op(v);
  }
  
};


///\en aligned contiguous memory allocation with allocation border of mult*sizeof(T) bytes
/// returns original pointer in orig_ptr which must be deleted by delete [] (T*)
template <class T>
inline T *aligned_alloc(size_t mult, size_t size, T **orig_ptr) {
  size_t nmult=size/mult;
  if(size%mult)
    nmult++;
  nmult++;
  *orig_ptr = new T[mult*nmult];
  
  unsigned long long pv=(unsigned long long)*orig_ptr, rest, div=mult*sizeof(T);
  rest=pv%div;
  if(rest)
    pv+=(div-rest);
  return reinterpret_cast<T*>(pv);

  //int mask = mult*sizeof(T)-1;
  
  //return (T *)(((size_t)((char *)*orig_ptr+mask))&(~mask));
}

# if 0
//e incremental iterator reading values from the list
template<class T>
class va_it{
protected:
  T v;
  va_list argptr;
public:
  typedef forward_iterator_tag iterator_category;
  typedef T value_type;
  typedef int difference_type;
  typedef difference_type distance_type;	// retained
  typedef T* pointer;
  typedef T& reference;

  va_it(const T &first, ...){
    va_start(argptr, first);
    v = va_arg( argptr, T);
  }
  
  T operator*(){
    return v;
  }
  
  va_it &operator++(){
    v = va_arg( argptr, T);
    *this;
  }

  va_it operator++(int){ // postfix
    value_it tmp=*this;
    ++*this;
    return tmp;
  }

  //value_it &operator--(){
  //  //v--;
  //  decr(v);
  //  return *this;
 // }

  // difference
  //T operator-(const value_it & other) const {
  //  return v-other.v;
  //}

  va_it operator--(int){ // postfix
    value_it tmp=*this;
    --*this;
    return tmp;
  }

  bool operator==(const va_it &other) const {
    return (argptr==other.argptr);
  }
  bool operator!=(const va_it &other) const {
    return !(*this==other);
  }

};

# endif

# endif
