/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.1 $
 *   $Date: 2012/12/06 02:24:04 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/logger.h,v 1.1 2012/12/06 02:24:04 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/logger.h,v $
$Revision: 1.1 $
$Author: lesha $
$Date: 2012/12/06 02:24:04 $
*/
/*s****************************************************************************
 * $Log: logger.h,v $
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.16  2012/09/25 19:48:30  lesha
 * documentation
 *
 * Revision 1.15  2012/03/13 17:14:50  valuev
 * added some tests, fixed memory leaks
 *
 * Revision 1.14  2011/09/15 00:10:55  lesha
 * emDataProjector is simplified
 *
 * Revision 1.13  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.12  2007/06/14 14:10:08  lesha
 * *** empty log message ***
 *
 * Revision 1.11  2007/06/14 14:09:23  lesha
 * AddField is corrected
 *
 * Revision 1.10  2007/03/27 22:37:58  valuev
 * Added timers
 *
 * Revision 1.9  2007/02/23 14:23:35  valuev
 * Added internal FFT to emSourceWave,
 * corrected transfers for PBC
 *
 * Revision 1.8  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.7  2007/02/16 22:54:18  lesha
 * disconnect is added. one little defect is corrected
 *
 * Revision 1.6  2006/12/26 10:46:55  valuev
 * Added range filters for detector t/f arguments
 *
 * Revision 1.5  2006/12/14 08:42:36  valuev
 * reformulated detector
 *  projectors, corrected open file limit control, tested Fourier sceleton
 *
 * Revision 1.4  2006/11/27 09:39:12  valuev
 * Added separators and plane mode to detectors
 *
 * Revision 1.3  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/

# ifndef _LOGGER_H
# define _LOGGER_H

/// @file loger.h \brief Interface for collecting process variables identified by their symbolic names

# include <map>
# include <string>
# include "refobj.h"

using namespace std;

/// Template base class for collecting process variables
/// identified by their symbolic names into output vector.
/// Process variables are of class out_t;
/// the argument of process function (input variable) is of inp_t.
/// This class is used as base logger class for file loggers (ex. mdlog.h)
template<class inp_t, class out_t>
class process_detail{
protected:
  struct field_t {
    string name;
    int on;
    out_t *ptr;
    field_t(const string &str="",int son=0,out_t *sptr=NULL):name(str),on(son),ptr(sptr){}
  };
  typedef map<string,int> map_t;
  map_t  fmap;
  vector<field_t> fields;
  
  
 
  /// transfers the values which are switched on to the vector
  /// to be used from process() of derivates
  virtual size_t TransferValues(vector<out_t> &descr) const { // my change
    descr.clear();
    int i, n=(int)fields.size();
    for(i=0;i<n;i++){
      if(fields[i].on){
        descr.push_back(*(fields[i].ptr));
      }
    }
    return descr.size();
  }


public:
  /// adds output descriptor, returns the data reference 
  /// to be used in constructors of derivates
  size_t AddField(out_t *ptr, const string& name,int on=1){
    int ind;
    if (fmap.find(name)!=fmap.end()) {
      ind=fmap[name];
      fields[ind].ptr=ptr;
      fields[ind].on=on;
      return fields.size();
    }
    fields.push_back(field_t(name,on,ptr));
    ind=(int)fields.size()-1;
    fmap[name]=ind;
    return ind+1;
  }

  /// swithces the output of specified field on or off
  virtual int SetOutFlag(const string &name, int flag=1){
    if(fmap.find(name)!=fmap.end()){
      fields[fmap[name]].on=flag;
      return 1;
    }
    return -1;
  }
  virtual int GetOutFlag(const string &name){
    if(fmap.find(name)!=fmap.end()){
      return fields[fmap[name]].on;
    }
    return 0;
  }
   
  virtual int SetOutBitFlags(long long bflags){
    int i, n=(int)fields.size();
    for(i=0;i<n;i++){
      long long fl=(long long)(1)<<i;
      if(bflags&fl)fields[i].on=1;
      else fields[i].on=0;
    }
    return 1;
  }

  /// collects the symbolic descriptors of output fields which are switched on
  virtual size_t GetFields(vector<string> &descr) const {
    descr.clear();
    int i, n=(int)fields.size();
    for(i=0;i<n;i++){
      if(fields[i].on){
        descr.push_back(fields[i].name);
      }
    }
    return descr.size();
  }
  /// transforms input value into a set of output values
  /// some of which (thouse switched on) are collected in vout
  /// @return number of fields in vout or <0 if the entry is invalid
  virtual int process(vector<out_t> &vout, const inp_t &v, int refnum=0)=0;


  virtual ~process_detail(){}
};



/// class for building filter chains
template <class inp_t> 
class entry_filter {
protected:
  mngptr<entry_filter> next;
public:
  
  entry_filter *append(entry_filter *subchain, int managed=0){
    entry_filter *cur=this;
    while(cur->next.ptr())cur=next.ptr();
    cur->next.reset(subchain,managed);
    return cur;
  }

  void disconnect() {
    next.reset();
  }

  /// constant test
  /// return values: 0 -- test failed, 1 -- test OK, -1 -- skip filter
  virtual int ctest(const inp_t &entry) const {
    if(next.ptr()){
      return next->ctest(entry);
    }
    return -1;
  }
};

/// links the chain with logical predicate
template <class inp_t, class pred_t= logical_and<bool> > 
class logic_filter: public entry_filter<inp_t> {
protected:
  pred_t comp;
  /// to be overriden by derivates
  virtual int mytest(const inp_t &entry) const {
    return -1;
  }
public:
 
  /// constant test
  /// return values: 0 -- test failed, 1 -- test OK, -1 -- skip filter
  virtual int ctest(const inp_t &entry) const {
    int rnext=-1, rthis=mytest(entry);
    if(this->next.ptr()){
      rnext=this->next->ctest(entry);
    }
    if(rnext==-1)return rthis;
    if(rthis==-1)return rnext;
    return (int)comp(rthis? true: false,rnext ? true : false);
  }

};

/// links the chain with logical predicate
template <class inp_t, class pred_t = logical_and<bool>, class lesseq_t = less_equal<inp_t> > 
class range_filter: public logic_filter<inp_t, pred_t>{
protected:
  lesseq_t leq;
  inp_t rmin, rmax;
  virtual int mytest(const inp_t &entry) const {
    if(leq(rmin,entry) && leq(entry,rmax)) return 1;
    return 0;
  }
public:
  range_filter(const inp_t &rmin_, const inp_t &rmax_):rmin(rmin_),rmax(rmax_){} 
};


# if 0

/// conditional set: makes subset of the initial set
/// based on some condition
template<class set_t, class cond_f>
class cond_set{
public:
  typedef typename set_t::iterator base_it;
  typedef typename set_t::value_type value_type;
  class iterator{
    cond_set *parent;
    base_it it;
    iterator(cond_set *parent_, int send=0):parent(parent_){
      if(send)it=parent->end;
      else it=parent->beg;
    }
  public:
    iterator(){}
   
    value_type &operator*(){
      return *it;
    }
    iterator &operator++(){
      do{
        ++it;
      }while(it!=parent->end && !parent->cond(*it));
      return *this;
    }

    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }
    
    bool operator!=(const iterator &other){
      return (it!=other.it);
    }
  };
protected:
  base_it beg, end;
  cond_f cond;
  /// size
  size_t n; 
public:
  cond_set(base_it beg_, base_it end_, const cond_f &f):beg(beg_),end(end_),cond(f),n(0){
    while(!cond(*beg)){
      if(!(beg!=end))return;
      ++beg;
    }
    n=1;
    base_it pend=beg, it=beg++;
    for(;it!=end;++it){
      if(cond(*it)){
        pend=it;
        n++;
      }
    }
    end=++pend;
  }

  size_t size() const{
    return n;
  }
  iterator begin(){
    return iterator(this);
  }
  iterator end(){
    return iterator(this,1);
  }
};


/// dimension checker, retuns true if the arguments 
/// with frozen dimensions coincide with tolerance
template<class arg_t>
struct dimcheck_f{
  int dim;
  arg_t sampl;
  long frozen;
  arg_t toler;
  dimchec_f(int dim_, const arg_t& sampl_,long frozen_, const arg_t &toler_):dim(dim_),sampl(sampl_),frozen(frozen_),toler(toler_){}
  bool operator()(const arg_t &val){
    int i;
    for(i=0;i<dim;i++){
      long k=1;
      if(frozen&(k<<i)){
        if(fabs(sampl[i]-val[i])>toler[i])return false;
      }
    }
  }
};

template<class inp_it>
class generic_itset{
  size_t n;
  inp_it beg, end;
public:
  typedef inp_it iterator;
  generic_itset(inp_it b, inp_it e):beg(b),end(e),n(0){}
  iterator begin(){ return beg;}
  iterator end(){ return end;}
  size_t size(){
    if(!n){
      inp_it it=b;
      while(it!=end){
        ++it;
        n++;
      }
    }
    return n;
  }
};

template<class set_t>
class dim_restrictor: public cond_set<set_t, dimcheck_f<typename set_t::value_type> >{
public:
  typedef cond_set<set_t, dimcheck_f<typename set_t::value_type> > base_t;
  typedef base_t::base_it base_it;
  typedef base_t::iterator iterator;
  typedef base_t::value_type value_typel
  dim_restrictor(int dim, long frozen,const value_type &toler,base_it beg_, base_it end_):base_t(beg_,end_,dimcheck_f(dim,*beg,frozen,toler)){}
};


// buffered convolution
template<class arg_tt, class val_tt>
class buff_conv{

};

template<class arg_tt, class val_tt>
class func_accessor{
  int dim;
public: 
  
  typedef arg_tt arg_t;
  typedef val_tt val_t;
  

  func_accesor(int dim_):dim(dim_){}

  /// getting internal iterator: any iterator will do
  template<class external_set>
  struct internal_set_traits{
    typedef dim_restrictor<external_set> set_t;
    typedef set_t::iterator iterator;
  };

  /// returns internal  set on the base of external one, freezes some dimensions
  template <class external_set>
  int get_internal_set(external_set &set,func_accessor<arg_tt,val_tt>::internal_set_traits<external_set>::set_t &iset, long freeze, arg_t toler){
    iset.init(dimension(),freeze,toler,set.begin(),set.end());
    return 1;
  }


  /// gets argument dimension
  int dimension() const {
    return dim;
  }

  /// gets the preferred order of iterating dimensions, -1 if not important
  template<class inp_it>
  int get_order(inp_it beg){
    return -1;
  }

  /// prepares a buffer to access data with arg_it
  template<class arg_it>
  int prepare(size_t dim, arg_it beg){
    return 1; // any iterator is OK
  }

  /// gets the value
  int get(const arg_t &arg, val_t &val){
    val=func(arg);
    return 1;
  }

  
};

template<class func_accessor_tt>
class convolution_accessor: public func_accessor_tt{
  size_t max_buf;
  mngptr<func_accessor_tt> func;
public:
  typedef func_accessor_tt base_t;
  typedef typename base_t::arg_t arg_t;
  typedef typename base_t::val_t val_t;
 

  convolution_accessor(const func_accesor_tt *func_=NULL, int managed=0){
    init(func_,managed);
  }
  void init(const func_accesor_tt *func_=NULL, int managed=0){
    func.reset(func,managed);
  }


  
  long get_signature(){}
  long set_signature(){}
  int get_conv_type(){}
  int set_conv_type(){}

   /// prepares a buffer to access convolution data with arg_it
  template<class arg_it>
  int prepare(size_t n, arg_it beg, long signature=0, int type=0){
    //1. getting convoluted dimensions
    int dim=func->dimension(), nconv=0;
    long isign=isignature();
    vector<int> convd(dim);
    int i;
    for(i=0;i<dim;i++){
      long k=1;
      long ec=signature&(k<<i);
      long ic=isign&(k<<i);
      if(ec!=ic){
        convd[i].c=1;
        nconv++;
      }
      else convd[i].c=0;
    }
    //2. determine the strategy of convolution: argument span
    size_t k;
    for(k=0;k<n;k++){
      for(i=0;i<dim;i++){
        arg_t x=(*beg)[i];
        convd[i].;
      
      }
    
    }



    return 1; // any iterator is OK
  }

  /// gets a value of current signature
  int get(const arg_t &arg, val_t &val){
    int recalc=1;
    if(prev_set){ // comparing with previous argument
      // is there a difference in non-convoluted dimensions?
      int i;
      recalc=0;
      for(i=0;i<dim;i++){
        if(!convd[i].c && arg[i]!=prev[i]){
          recalc=1;
        }
      }
    }
    if(recalc){
      if(conv_acc) delete conv_acc;
      // getting iterator form function accessor
      best_it beg=func->int_begin(rbeg,rend,&nrange);
      dim_restr_set<arg_t> rset(func->dimension(),func->int_begin(rbeg,rend),arg,sign);
      conv_acc= convolution(nconv,convit,sign,type);
    }
    if(!conv_acc)return -1;
    return conv_acc->get(arg,val);
  }
  

};
# endif

# endif
