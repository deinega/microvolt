/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
*****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/transfer.h,v $
$Revision: 1.17 $
$Author: valuev $
$Date: 2013/10/28 08:21:56 $
*/
/******************************************************************************
 * $Log: transfer.h,v $
 * Revision 1.17  2013/10/28 08:21:56  valuev
 * removed emtype
 *
 * Revision 1.16  2013/10/27 22:30:16  valuev
 * fixed support for COMPLEX
 *
 * Revision 1.15  2013/10/14 07:53:02  lesha
 * InterpTransfer
 *
 * Revision 1.14  2013/10/08 04:07:48  lesha
 * bufptr is returned back
 *
 * Revision 1.13  2013/10/03 12:45:15  valuev
 * edited comments
 *
 * Revision 1.12  2013/09/30 00:28:35  lesha
 * nothing important
 *
 * Revision 1.11  2013/09/29 06:22:57  lesha
 * slight change
 *
 * Revision 1.10  2013/09/29 06:11:57  lesha
 * transfer documentation
 *
 * Revision 1.9  2013/09/29 01:33:53  lesha
 * documentation
 *
*******************************************************************************/

/** @file transfer.h
    Interpolation forms, packers and interpolation transfers (local and MPI), which can be used for 
    parallel simulation (transfer between CPU domains) or multigrid-approach (transfer between meshes).
**/


#ifndef _TRANSFER_H
#define _TRANSFER_H

#include <map>
#include "pencil.h"
#include "seqpack.h"
#include "component.h"

using namespace std;

///\en Class for interpolation (set of indices in some array and coefficients).
/// Indices and coefficients are recorded in a vector
/// This vector could be stored in this class or in some external class,
/// which records other interpolations to this vector as well (see InterpPacker)
/// coef_tt is type for interpolation coefficients (for example, double or complex)
/// index_t is type for array indexation (int for 32-bit architecture)
template<class coef_tt, class index_tt = ptrdiff_t>
class Interpolation{
public:
  typedef coef_tt coef_t;
  typedef index_tt index_t;
  typedef typename vector<index_t>::iterator shift_it;
  typedef typename vector<coef_t>::iterator coeff_it;

  index_t m_ind; ///\en all indices are shifted from this index. m_ind is always non-negative
  pencil<index_t, vector<index_t> > shifts; ///\en interpolation indices (shifted from m_ind)
  pencil<coef_t, vector<coef_t> > coeffs; ///\en interpolation coefficients

  ///\en default constructor
  Interpolation():m_ind(0){}
  
  ///\en copy constructor
  Interpolation(const Interpolation &other):m_ind(other.m_ind),shifts(other.shifts),coeffs(other.coeffs){}
  
  ///\en reference copy
  Interpolation& operator=(const Interpolation &other){
    if(this!=&other){
      m_ind=other.m_ind;
      shifts=other.shifts;
      coeffs=other.coeffs;
    }
    return *this;
  }

  ///\en constructor from a map {ind, coeff} 
  Interpolation(const map<index_t,coef_t> &cmap){
    typename map<index_t,coef_t>::const_iterator it=cmap.begin(), e=cmap.end();
    index_t n=(index_t)cmap.size();
    shifts.resize(n,n);
    coeffs.resize(n,n);
    m_ind = n ? it->first : 0;
    for(n=0;it!=e;++it,++n){
      shifts[n]=it->first-m_ind;
      coeffs[n]=it->second;
    }
  }

  ///\en constructor from shifts and coeffs input iterators
  /// using data copy: managed is a dummy and assumed to be positive
  template<class sh_it, class c_it>
  Interpolation(index_t m_ind_, index_t n_, sh_it shifts_, c_it coeffs_, int managed):
  m_ind(m_ind_),shifts(n_, n_),coeffs(n_, n_){
    for(int i=0;i<shifts.size();i++){
      shifts[i]=*shifts_++;
      coeffs[i]=*coeffs_++;
    }
  }

  ///\en specialization for vector input iterators
  /// uses pencil referencing (no data copy)
  Interpolation(index_t m_ind_, index_t n_, typename Interpolation::shift_it shifts_, typename Interpolation::coeff_it coeffs_):
  m_ind(m_ind_),shifts(shifts_,(shift_it)(shifts_+n_),n_),coeffs(coeffs_,(coeff_it)(coeffs_+n_),n_){}

  bool valid() const {
    return shifts.size()!=0;
  }

  index_t set_m_ind(index_t m_ind_){
    index_t tmp=m_ind;
    m_ind=m_ind_;
    return tmp;
  }

  index_t get_m_ind() const {
    return m_ind;
  }

  int set_size(const int size) {
    shifts.n=size;
    coeffs.n=size;
    return 1;
  }

  int size() const {
    return shifts.size();
  }

  shift_it get_shifts() const{
    return shifts.get_ptr();
  }

  coeff_it get_coeffs() const{
    return coeffs.get_ptr();
  }

  template<class T>
  Interpolation &operator*=(const T &c){
    for(int i=0;i<shifts.size();i++)
      coeffs[i]*=c;
    return *this;
  }

  /*Interpolation &operator+=(const Interpolation &other){
    *this=linear_mix(1.,*this,1.,other);
    return *this;
  }*/

  ///\en record indices and coefficients to a map (key - index, value - coefficient)
  /// all coefficients are multiplied by c1
  template<class T>
  int add_to(const T& c1, map<index_t,coef_t> &cmap) const{
    for(int i=0;i<shifts.size();i++){
      index_t ind=m_ind+shifts[i];
      if(cmap.find(ind)==cmap.end())cmap[ind]=coeffs[i]*c1;
      else cmap[ind]+=coeffs[i]*c1;
    }
    return shifts.size();
  }

  ///\en record indices and coefficiens to corresponding vectors
  int copy_to(vector<index_t> &vsh, vector<coef_t> &vcf) const {
    for(int i=0;i<shifts.size();i++){
      vsh.push_back(m_ind+shifts[i]);
      vcf.push_back(coeffs[i]);
    }
    return shifts.size();
  }

  ///\en gets the value associated with the interpolation entry
  template<class T>
  T get_value(const T *ptr) const {
    T sum=0;
    for(int k=0;k<shifts.size();k++){
      sum+=ptr[m_ind+shifts[k]]*coeffs[k];
    }
    return sum;
  }
};


///\en prepares the new interpolation which is a linear combination of the arguments
template<class interp_t, class T> 
interp_t linear_mix(const T& c1, const interp_t &ip1, const T& c2, const interp_t &ip2){
  map<typename interp_t::index_t,typename interp_t::coef_t> cmap;
  ip1.add_to(c1,cmap);
  ip2.add_to(c2,cmap);
  return interp_t(cmap);
}
/*template<class coef_t, class T> 
Interpolation<coef_t> linear_mix(const T& c1, const Interpolation<coef_t> &ip1, const T& c2, const Interpolation<coef_t> &ip2){
  map<int,coef_t> cmap;
  ip1.add_to(c1,cmap);
  ip2.add_to(c2,cmap);
  return Interpolation<coef_t>(cmap);
}*/

///\en class to pack Interpolations with STL-like interface (iterators, which are returned by begin and end)
template<class coef_tt, class val_t, class index_tt = ptrdiff_t>
class InterpPacker{
public:
  typedef coef_tt coef_t;
  typedef index_tt index_t;
  typedef Interpolation<coef_t, index_t> value_t;
  typedef Interpolation<coef_t, index_t> record_t;
protected:
  typedef typename pack_t<index_t, index_t>::data mind_t;
  typedef typename pack_t<index_t, index_t>::group sh_t; 
  typedef typename pack_t<coef_t, index_t>::group cf_t; 
  typedef typename mind_t::iterator mind_it;
  typedef typename sh_t::iterator sh_it;
  typedef typename cf_t::iterator cf_it;

  int rec_state; ///<\en record state: 0=initial, 1=recording, 2=recorded
  val_t *ptr; ///<\en data array corresponding to interpolations indices

  mind_t mindpack; ///<\en packer of interpolation main indices
  sh_t shiftpack; ///<\en packer of interpolation indices 
  cf_t coeffpack; ///<\en packer of interpolation coefficients 
public:
  class iterator{
  protected:
    friend class InterpPacker;
    mind_it iti;
    sh_it its;
    cf_it itc;
    iterator(mind_it iti_, sh_it its_, cf_it itc_): iti(iti_), its(its_), itc(itc_){}
  public:  
    typedef nogroup_it<iterator> gcategory;
    iterator() {}
    ///\en copy constructor
    iterator(const iterator &other):iti(other.iti),itc(other.itc),its(other.its){}
    
    Interpolation<coef_t, index_t> operator*() const {
      return Interpolation<coef_t, index_t>(*iti,get_group_count(its),get_group_begin(its),get_group_begin(itc));
    }

    ///\en prefix increment
    iterator& operator++(){
      ++iti;
      its.plus();
      itc.plus();
      return *this;
    } 
    ///\en postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    } 

    bool operator!=(const iterator &other){
      return (iti!=other.iti);
    } 

  };
#ifndef USE_PACKERS
  InterpPacker(val_t *sptr=NULL): ptr(sptr), rec_state(0){
#else
  InterpPacker(val_t *sptr=NULL): ptr(sptr), rec_state(0), mindpack(1){
#endif
    // synchronization of shifpack and coeffpack
    coeffpack.aux_depends_on(shiftpack);
  }

  void SetDataPtr(val_t *sptr){
    ptr=sptr;
  }

  typename pack_t<index_t, index_t>::data *get_grouper() const {
    return shiftpack.get_grouper();
  }

  int start_record(){
    if(rec_state!=1) {
      mindpack.start_record();
      shiftpack.start_record();
      coeffpack.start_record();
      rec_state=1;
    }
    return 1;
  }

  ///\en records new interpolation
  /// interpolation is not necessary not empty
  int next_record(const Interpolation<coef_t, index_t> &value){
    start_record();
    mindpack.next_record(value.get_m_ind());
    shiftpack.next_group(value.get_shifts(),value.size());
    return coeffpack.next_group(value.get_coeffs(),value.size());
  }

  typedef map<int,coef_t> map_t;
  typedef typename map<int,coef_t>::const_iterator map_it;

  template<class pair_it>
  struct first_it: public pair_it{
    typedef typename pair_it::value_type::first_type value_type;
    first_it(const pair_it &it):pair_it(it){}
    value_type operator*(){
      return (this->pair_it::operator*()).first;
    }
  };

  template<class pair_it>
  struct second_it: public pair_it{
    typedef typename pair_it::value_type::second_type value_type;
    second_it(const pair_it &it):pair_it(it){}
    value_type operator*(){
      return (this->pair_it::operator*()).second;
    }
  };

  ///<\en records interpolation from a map (key - index, value - coefficient)
  int next_record(const map_t &coeffmap) {
    start_record();
    mindpack.next_record(0);
    coeffpack.next_group(second_it<map_it>(coeffmap.begin()), (int)coeffmap.size());
    shiftpack.next_group(first_it<map_it>(coeffmap.begin()), (int)coeffmap.size());
    return 1;
  }

  ///<\en records interpolation from indices and coefficients vectors
  int next_record(const vector<index_t> &vind, const vector<coef_t> &vcoeff){
    start_record();
    mindpack.next_record(0);
    coeffpack.next_group(vcoeff.begin(), (int)vcoeff.size());
    shiftpack.next_group(vind.begin(), (int)vind.size());
    return 1;
  }

  int end_record(){
    if(rec_state!=2) {
      mindpack.end_record();
      shiftpack.end_record();
      coeffpack.end_record();
      rec_state=2;
    }
    return 1;
  }
  
  ///<\en get interpolated value
  val_t gather(iterator data) const {
    val_t res=0;
    index_t mind=*data.iti;
    index_t nc=get_group_count(data.itc);
    typename pack_t<coef_t, index_t>::group::iterator::group_it ci=get_group_begin(data.itc);
    for (index_t i=0; i<nc; i++){
      index_t ind=*data.its;
      coef_t coef=*ci;
      res+=coef*ptr[mind+ind];
      ci++;
      data.its++;
    }
    return res;
  }
  
  ///<\en distributes a value: v_new=c1*v_old+val*packed_coeff
  index_t scatter(iterator data, val_t val, val_t c1=1) const{
    index_t mind=*data.iti;
    index_t nc=get_group_count(data.itc);
    typename pack_t<coef_t, index_t>::group::iterator::group_it ci=get_group_begin(data.itc);
    for (index_t i=0; i<nc; i++){
      index_t ind=*data.its;
      coef_t coef=*ci;
      ptr[mind+ind]=c1*ptr[mind+ind]+coef*val;
      ci++;
      data.its++;
    }
    return nc;
  }
#ifndef USE_PACKERS
  ///<\en random access to packed interpolation
  Interpolation<coef_t, index_t> operator[] (size_t i) {
    end_record();
    iterator it(mindpack.begin()+i,shiftpack[i],coeffpack[i]);
    return *it;
  }
#endif
  ///\en returns iterator on a first element of packed interpolations
  iterator begin() {
    end_record();
    return iterator(mindpack.begin(),shiftpack.begin(),coeffpack.begin());
  }
  
  iterator end()  {
    end_record();
    return iterator(mindpack.end(),shiftpack.end(),coeffpack.end());
  }

  ///<\en returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return mindpack.packed_size()+shiftpack.packed_size()+coeffpack.packed_size();
  }

  ///<\en returns the unpacked  size in bytes ignoring vector overheads
  size_t data_size() const {
    return mindpack.data_size()+shiftpack.data_size()+coeffpack.data_size();
  }

  size_t size() const {
    return mindpack.size();
  }

  void clear(){
    mindpack.clear();
    shiftpack.clear();
    coeffpack.clear();
  }
};

///\en Simple implementation of interpolation form (but not efficient in terms of used memory)
template<class coef_t, int N, class index_t = ptrdiff_t >
struct StaticInterpolation{
  index_t shifts[N]; ///<\en interpolation indices
  coef_t coeffs[N]; ///<\en interpolation coefficients
  int n; ///<\en number of terms in interpolation

  StaticInterpolation():n(0){}

  int size() const {
    return n;
  }

  bool valid() const {
    return n!=0;
  }

  StaticInterpolation operator +(const StaticInterpolation &other){
    StaticInterpolation r;
    for(int i=0;i<n;i++){
      r.shifts[i]=shifts[i];
      r.coeffs[i]=coeffs[i];
    }
    for(int i=0;i<other.n;i++){
      r.shifts[n+i]=other.shifts[i];
      r.coeffs[n+i]=other.coeffs[i];
    }
    r.n=n+other.n;
    if(n>=N)
      throw 1;
    return r;
  }

  StaticInterpolation operator*(const coef_t &coeff) const{
    StaticInterpolation r=*this;
    for(int i=0;i<n;i++)
      r.coeffs[i]*=coeff;
    return r;
  }

  template<class T>
  T get_value(const T *ptr) const {
    T sum=0;
    for(int k=0;k<size();k++){
      sum+=ptr[shifts[k]]*coeffs[k];
    }
    return sum;
  }
};

template<class val_t, int N, class index_t, class T>
StaticInterpolation<val_t,N, index_t> linear_mix(const T& c1, const StaticInterpolation<val_t,N, index_t> &ip1, 
const T& c2, const StaticInterpolation<val_t,N, index_t> &ip2){
  return ip1*c1+ip2*c2;
}

/**\en Class for interpolation transfers (for example, between CPU domains in parallel simulation or between meshes in multigrid-approach).
It sends values (for example, electromagnetic field in FDTD method) at chosen space locations to some target memory arrays.
These values are in turn obtained as interpolated by mesh. 
If the target array and the mesh array belong to different CPU domains, 
then the MPI procedures are used to perform the transfer, otherwise the data is simply copied. 
This (transfer) class can be applied for different purposes. 
For example, if the output and input points of some discretized update equation on a mesh belong to different CPU domains, 
the transfer class sends the value at the input point between these two domains to perform the discretized update equation on a mesh. 
Transfer class may be used to group all similar transfers together to benefit from MPI non blocking 
synchronization mechanisms triggered at a right time. 
The transfer class can also be used if some domain collects data values at chosen points of the calculation volume 
in order to record them in a file or for further analysis.

Work with object of this class assumes 2 stages.

During the first stage, all the arguments where data values are required, 
should be registered in the transfer objects along with the target arrays and the indices in these arrays, 
where the data values should be sent (put_request, put_vrequest).
Target arrays are assigned by register_buffer, substitute_buffer.
After all arguments are registered, MPI-procedures should be intialized (prepare_transfer).

During the second stage, the transfer object 
(a) asks meshes which contain the registered arguments to interpolate the data values there, 
(b) sends these values to the corresponding indices of the target arrays.
Local transfers are called by compute_local.
MPI-transfers are called by start_tranfer and complete_transfer.

arg_tt is argument type for which value is collected
(for example, in EMTL this is structure which contains position, field type and field direction;
value is field projection on this direction).
container_t is container (mesh or association of meshes) which manages values array.
interp_packer_t is packer of interpolations (interpolations are structures for extraction of interpolated values from container).
index_t is type for array indexation (int for 32-bit architecture)
*/
// TODO: index_t
template<class arg_tt, class container_tt, class reg_tester_t=container_tt, 
class interp_packer_t=typename container_tt::interp_packer_t, class index_t = typename interp_packer_t::index_t>
class InterpTransfer: public apComponent{

public:

  typedef arg_tt arg_t;
  typedef container_tt container_t;
  typedef typename container_tt::interp_form_t interp_form_t;
  typedef typename container_tt::interp_form_t::coef_t coef_t;
  typedef typename container_tt::val_t val_t; ///<\en type of interpolated value (for example, double or complex)

  typedef pair_pack<typename container_tt::interp_form_t, pair<int,int>, 
    interp_packer_t, pair_pack<int,int,typename pack_t<int, index_t>::data,typename pack_t<int, index_t>::data> > base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::value_t value_t;

protected:

  container_t *cont; ///\en container for collecting interpolated values

  int myrank; ///\en processor rank
  int np; ///\en number of processors
  int rec_state; ///\en record state: 0=initial, 1=recording, 2=recorded

  base_t ip; ///\en packer of mpi-transfers: interpolation - pair(buffer number, index)
  base_t iplocal; ///\en packer for local non-mpi transfers: interpolation - pair(buffer number, index)
  vector<val_t *> buff; ///\en destination buffers. first np buffers are reserved as output buffers to other processors

  int mpidump; ///\en if 1, then logfile will be recorded
  FILE *logfile; ///\en logfile for mpi-transfers
  static int idcount; ///\en number of created object, used in logfile

  typedef void request_t;
  typedef void status_t;

  ///\en first element - bufer for MPI_Request in MPI_Send, second element - buffer for MPI_Request for MPI_Recv
  struct reqpair_t: public  pair<request_t *, request_t *> {
    int dest; ///\en destination processor rank
    reqpair_t(request_t *a, request_t *b, int dest_=0):pair<request_t *, request_t *>(a,b),dest(dest_){}
  };

  ///\en interpolation request
  struct irequest_t{
    int dest; ///<\en destination processor rank
    int ibuf; ///<\en buffer (memory array) number at this processor
    int ind; ///<\en index in memory array
    ///\en auxiliary flag to regulate depth of interpolations:
    /// sometimes, mesh interpolates value from other mesh, which interpolate value from other mesh, etc.
    /// this flag is used to stop possible bad-recursions
    int level;
    arg_t arg; ///\en for which argument value should be collected
    
    irequest_t(int dest_=0,int ibuf_=0, int ind_=0, int level_=0, arg_t arg_=arg_t()):
      dest(dest_),ibuf(ibuf_),ind(ind_),level(level_),arg(arg_){}
  };

  vector<vector<irequest_t> > rankvec; ///\en vectors of interpolation requests for each processor

  typedef pair<int,int> intpair_t;

# ifdef USE_MPI
  void *mycomm; ///\en communicator for this instance
  void *ireq_mpit, *intpair_mpit;   
# endif

  ///\en manages two buffers: to send (receive) value to (from) another processor
  struct comm_t{
# ifdef USE_MPI
    void *mycomm;
# endif
    val_t *bufptr; ///<\en points at outbuf or inbuf (depending on situation)
    val_t *outbuf; ///<\en buffer to send data to another processor
    size_t outsz; ///<\en size of outbuf
    val_t *inbuf; ///<\en buffer to receive data from another processor
    size_t insz; ///<\en size of inbuf
    typedef pair_pack<int,int, typename pack_t<int, index_t>::data,typename pack_t<int, index_t>::data> packer_t;
    ///\en packer of pairs 'buffer number - index in buffer array',
    /// used to distribute values recieved from other processors:
    /// these values are initially received to input buffer,
    /// and then recorded to destination buffers buff
    packer_t destpack; 

    comm_t():outsz(0),insz(0),outbuf(NULL),inbuf(NULL){
      destpack.start_record();
    }

    ///\en allocate memory for out- and in- buffers
    void init(){
      destpack.end_record();
      clear_buffers();
      if(outsz)outbuf = new val_t [outsz];
      if(insz)inbuf = new val_t [insz];
    }

    ///\en switch pointer bufptr
    void reset_buffer(int out){
      bufptr=(out ? outbuf : inbuf );
    }

    ///\en adds inbuf values to destination buffers
    void transfer_dest(vector<val_t *> &buff){
      typename packer_t::iterator it=destpack.begin(), e=destpack.end();
      size_t nb=buff.size();
      for(size_t i=0; it!=e;++it){
        pair<int,int> v=*it;
        if (buff[v.first/2]==NULL)
          continue;
        buff[v.first/2][v.second]+=inbuf[i++];
      }
    }

//    void transfer_dest(vector<val_t *> &buff, int *shifts){
//      packer_t::iterator it=destpack.begin(), e=destpack.end();
//      size_t i=0, nb=buff.size();
//      for(;it!=e;++it){
//        pair<int,int> v=*it;
//        buff[v.first/2][shifts[v.second]]+=inbuf[i++];
//      }
//    }

    ///\en initialize MPI_Send and MPI_Recv procedures
    void make_requests(reqpair_t arg, apMPIComm *comm){
# ifdef USE_MPI
      if(arg.first)
        //MPI_Send_init(outbuf,(int)outsz, EM_MPI_VALTYPE, arg.dest, 0, mycomm, arg.first);
        comm->vec_type_send_init(outbuf, (int)outsz*(sizeof(val_t)/sizeof(vec_type)), arg.dest, 0, mycomm, arg.first);

      if(arg.second)
        //MPI_Recv_init(inbuf,(int)insz, EM_MPI_VALTYPE, arg.dest, 0, mycomm, arg.second);
        comm->vec_type_recv_init(inbuf, (int)insz*(sizeof(val_t)/sizeof(vec_type)), arg.dest, 0, mycomm, arg.second);
# endif
    }

    ///\en remove initialized mpi-procedures
    void free_requests(reqpair_t arg, apMPIComm *comm){
# ifdef USE_MPI
      if(arg.first)
        //MPI_Request_free(arg.first);
        comm->free_request(arg.first);
      if(arg.second)
        //MPI_Request_free(arg.second);
        comm->free_request(arg.second);
# endif
    }

    size_t packed_size() const{
      return destpack.packed_size();
    }

    size_t data_size() const{
      return destpack.data_size();
    }

    void clear_buffers(){
      if(outbuf) delete [] outbuf;
      outbuf=NULL;
      if(inbuf) delete [] inbuf;
      inbuf=NULL;
    }

    ~comm_t(){
      clear_buffers();
    }
  };

  refvector<comm_t> commtbl; ///<\en comm_t structures for each processor

  int nsend, nrecv; ///<\en number of those processors where nonzero size data should be sent (received from)
  void *srequests, *rrequests; ///<\en arrays for MPI_Request of sizes nsend, nrecv
  void *sstatuses, *rstatuses; ///<\en arrays for MPI_Status of sizes nsend, nrecv

  int *rvecsz,*svecsz,*rdispl,*sdispl; ///<\en auxiliary arrays used in transfer_requrest and prepare_transfers

  struct{
    int tr_prep, tr_lcl, tr_gthr, tr_sctr, tr_send, tr_recv;
  } tid;

  void complete_record(){
    ip.end_record();
    iplocal.end_record();
    rec_state=2;
  }

  ///\en allocate output and input arrays in commtbl and initizlize corresponding MPI_Send and MPI_Recv procedures
  void init_mpi_table();

  ///\en deinitizlize MPI_Send and MPI_Recv procedures (initialized in init_mpi_table)
  void free_mpi_requests();

  ///\en fills (destination or mpi-output) buffers with interpolated values
  int fill_buffers(iterator it, iterator e);

//  int fill_buffers(iterator it, iterator e, int *shifts);

  ///\en processes requested interpolations (they are recorded at rankvec).
  /// update iplocal for local interpolation, ip and commtbl for mpi-interpolaton.
  /// secondary interpolation requests are rerecorded to rankvec
  int transfer_requests();

/*  template<class T>
  void reval(T o, T &i){
    i=o;
  }
  
  template<class T>
  void reval(complex<T> o, T &i){
    i=o.real();
  }
*/
public:
  InterpTransfer(container_t *cont_=NULL,interp_packer_t *packer_= new interp_packer_t(),interp_packer_t *packer_loc_= new interp_packer_t() );

  void AddDefaultTimers(const vector<int> &p_id);

  void set_container(container_t *cont_) {
    if (rec_state<0)
      rec_state=0;
    cont=cont_;
  }

  ///\en identify node which mantains given position
  int test_rank(Vector_3 &pos) const{
    return cont->test_rank(pos);
  }

  ///\en returns nonlocal packer iterator
  /// ensures that recording is correctly finished
  iterator begin(){
    if(rec_state!=2)
      complete_record();
    return ip.begin();
  }

  ///\en returns local packer iterator
  /// ensures that recording is correctly finished
  iterator begin_local(){
    if(rec_state!=2)
      complete_record();
    return iplocal.begin();
  }

  ///\en returns nonlocal end
  iterator end(){
    return ip.end();
  }

  ///\en returns local end
  iterator end_local(){
    return iplocal.end();
  }

  ///\en returns packed size in bytes
  size_t packed_size() const {
    size_t sz=ip.packed_size()+iplocal.packed_size();
    for(size_t i=0;i<commtbl.size();i++)
      sz+=commtbl[i]->packed_size();
    return sz;
  }

  ///\en returns full unpacked sequence size in bytes 
  size_t data_size() const {
    size_t sz=ip.data_size()+iplocal.data_size();
    for(size_t i=0;i<commtbl.size();i++)
      sz+=commtbl[i]->data_size();
    return sz;
  }

  size_t size() const {
    return ip.size()+iplocal.size();
  }

  ///\en gets total buffer sizes (for info)
  int mpi_size(int out=0) const {
    int sum=0;
    if(out){
      for(size_t i=0;i<commtbl.size();i++)
        sum+=commtbl[i]->outsz;
    }
    else{
      for(size_t i=0;i<commtbl.size();i++)
        sum+=commtbl[i]->insz;
    }
    return sum*sizeof(val_t);
  }

  void clear(){
    ip.clear();
    iplocal.clear();
  }

  ~InterpTransfer(){
    if(logfile)fclose(logfile);
    free_mpi_requests();
#ifdef USE_MPI
    //if(mycomm!=MPI_COMM_NULL)MPI_Comm_free(&mycomm);
    comm->remove_comm(mycomm);

    //MPI_Type_free(&ireq_mpit); 
    comm->free_type(ireq_mpit);
    //MPI_Type_free(&intpair_mpit); 
    comm->free_type(intpair_mpit);
#endif
  }

#if 0
// THIS IS OLD DUMPER, SHOULD BE MODIFIED 
  class dumper_t{
    emDumpInfo *from, *to;
    InterpTransfer *parent;
  public:
    dumper_t(InterpTransfer *sparent, emDumpInfo *sfrom,emDumpInfo *sto):parent(sparent),from(sfrom),to(sto){}
  
    /// dumps the intialized data
    void dump_on_init(int gind,const interp_form_t &coeff){  
      vector<int> vind;
      vector<coef_t> vcoeff;
      vector<Vector_3> vinp;
      parent->cont->decode_interpolation(&coeff,vind,vcoeff);
      for(int i=0;i<(int)vind.size();i++){
        int sftype;
        Vector_3 sfdir;
        vinp.push_back(parent->cont->get_position(vind[i],sftype,sfdir));
      }
//      theDump.DumpTransfer(gind,vind,vinp,vcoeff,*from,*to);
    }
    /// dumps the data in action (uses operator())
    void dump_on_action(const value_t &data){
      (*this)(data);
    }

    /// operator to be subsituted in Apply(...) 
    /// instead of main operator() to check the action
    int operator()(const value_t &data) const {
      vector<int> vind;
      vector<coef_t> vcoeff;
      vector<Vector_3> vinp;
      parent->cont->decode_interpolation(&data.first,vind,vcoeff);
      for(int i=0;i<(int)vind.size();i++){
        int sftype;
        Vector_3 sfdir;
        from->pos=parent->cont->get_position(vind[i],sftype,sfdir);
        vinp.push_back(from->pos);
        from->ftype=to->ftype=sftype;
        from->fdirv=to->fdirv=sfdir;
      }
//      theDump.DumpTransfer(data.second.second,vind,vinp,vcoeff,*from,*to); // ?
      
      return 1;
    }
  };
#endif

  ///\en register destination buffer array where transfered data should be collected
  /// the pointer bufp is not managed (must be deleted externally)
  /// @return the buffer index
  /// if memory for buffer is not allocated yet, then bufp can be NULL; 
  /// after memory is allocated, pointer at the buffer array should be substituted using substitute_buffer
  int register_buffer(val_t *bufp){
    buff.push_back(bufp);
    return (int)(buff.size()-1);
  }

  ///\en substitute pointer at the buffer array, see also register_buffer
  val_t *substitute_buffer(int bufnum,val_t *ptr){
    val_t *tmp=buff[bufnum];
    buff[bufnum]=ptr;
    return tmp;
  }

  ///\en record argument (arg), index ind at the buffer ibuf, where value (corresponding to arg) should be transfered
  /// action=1 means overwrite existing data at the buffrer, 0 means add
  /// (imlementation: local part is recorded to iplocal and nonlocal part is recorded to rankvec)
  int put_request(int ibuf,int ind, const arg_t &arg, int action/*, emDumpInfo *from=NULL, emDumpInfo *to=NULL*/){
    vector<arg_t> argsum(1,arg);
    return put_vrequest(ibuf,ind,argsum,1.,action/*,from,to*/);
  }

  ///\en the same as put_request, but summed up value for all arguments multiplied on c_sum will be recorded at the buffer
  int put_vrequest(int ibuf, int ind, vector<arg_t> &argsum, coef_t c_sum, int action/*, emDumpInfo *from=NULL, emDumpInfo *to=NULL*/);

  ///\en put requests for all possible arguments at position pos
  int put_full_request(int ibuf,int &ind, const Vector_3& pos, int action=1){
    vector<arg_t> args = arg_t::full_args(pos);
    for(size_t i=0;i<args.size();i++){
      if(put_request(ibuf, ind++, args[i], action)<0)
        return -1;
    }
    return 1;
  }

  ///\en return size of values for all possible arguments
  size_t get_full_data_size(){
    return arg_t::get_full_data_size();
  }

  ///\en use less than 3 components of E and H fields in put_full_request (not implemented)
  void set_full_data(int outtype_){}

  ///\en make transfer negotiation
  /// (implementation: commtbl and MPI-procedures are inintialized)
  int prepare_transfers();

  ///\en makes non-mpi part of transfers
  int compute_local(){
    start(tid.tr_lcl);
    fill_buffers(begin_local(),end_local());
    stop(tid.tr_lcl);
    return 1;
  }

//  int compute_local(int *shifts){
//    return fill_buffers(begin_local(),end_local(), shifts);
//  }

  ///\en prepare output buffers and start non-blocking mpi-transfers
  void start_transfers();

  ///\en complete non-blocking mpi-transfers
//  void complete_transfers(int *shifts=NULL);
  void complete_transfers();

};

#endif
