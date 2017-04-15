/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.12 $
 *   $Date: 2013/11/02 17:14:53 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/grid.h,v 1.12 2013/11/02 17:14:53 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/grid.h,v $
$Revision: 1.12 $
$Author: lesha $
$Date: 2013/11/02 17:14:53 $
*/
/*s****************************************************************************
 * $Log: grid.h,v $
 * Revision 1.12  2013/11/02 17:14:53  lesha
 * nothing important
 *
 * Revision 1.11  2013/10/28 20:46:44  lesha
 * volume_integral_const_t is moved back
 *
 * Revision 1.10  2013/10/26 23:05:22  lesha
 * ptr -> mngptr
 *
 * Revision 1.9  2013/10/26 20:25:19  lesha
 * Integrate is added to grid
 *
 * Revision 1.8  2013/10/15 16:06:01  lesha
 * NonuniformGrid::index_t bug is fixed
 *
 * Revision 1.7  2013/10/15 16:03:41  lesha
 * Nonuniform grid, index_t bugs are fixed
 *
 * Revision 1.6  2013/09/05 15:16:26  zakirov
 * added parametrization on index_t
 *
 * Revision 1.5  2013/06/28 02:37:28  lesha
 * spatial integration is added to Project (emValue is modified)
 *
 * Revision 1.4  2013/06/24 11:42:03  valuev
 * xfel reader
 *
 * Revision 1.3  2013/06/18 16:40:15  valuev
 * added constructor from another grid
 *
 * Revision 1.2  2013/04/19 10:28:27  valuev
 * added [] operator
 *
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.57  2012/10/26 05:43:07  lesha
 * nothing important
 *
 * Revision 1.56  2012/10/17 21:44:26  lesha
 * *** empty log message ***
 *
 * Revision 1.55  2012/09/30 19:32:18  lesha
 * grid is modified
 *
 * Revision 1.54  2012/09/29 20:48:39  lesha
 * NonUniformGrid is simplified
 *
 * Revision 1.53  2012/09/28 21:50:37  lesha
 * documentation
 *
 * Revision 1.52  2012/09/21 20:24:58  lesha
 * renaming some function
 *
*******************************************************************************/
# ifndef _GRID_H
# define _GRID_H

#include "region_3.h"

const int global_il=1;

/**\en uniform grid in space

 grid can interpolate value at its internal points (which are inside box b) by function GetCoeff,
 return iterator on its nodes,
 organize memory access to values at its nodes by functions pack_ind and unpack_ind

 grid size (gr), start node (st), end node (en) along each direction for these these functionalities
 are defined by arrays with names starting with n (interpolation), i (interpolation), m (memory)

 usually they are the same for each functionality, 
 but for some special cases (like parallel implementation) could be different

 memory access is managed by interleave parameter.
 here we explain it using example of 1d grid of 10 elements.
 if interleave = 1, then data is stored in the natural order:
 0 1 2 3 4 5 6 7 8 9
 if interleave = 2, then data is stored in a different way:
 0 2 4 6 8 1 3 5 7 9
 if interleave = 5, then we have:
 0 5 1 6 2 7 3 8 4 9
 total grid size GSZ should be multiple if interleave: GSZ = interleave * nintz, where nintz is some integer number

 in 3d case, index change order is z(2), y(1), x(0), for example for grid 4x4x4 (inerleave=1):
 111 112 113 114   121 122 123 124   131 132 133 134   141 142 143 144
 211 212 213 214   221 222 223 224   231 232 233 234   241 242 243 244
 311 312 313 314   321 322 323 324   331 332 333 334   341 342 343 344
 411 412 413 414   421 422 423 424   431 432 433 434   441 442 443 444

 for interleave!=1 there are two possible cases: global_il!=0 or global_il=0
 which we consider on example interleave=2:

   if global_il!=0, then interleave works "globaly":
   111 113 121 123 131 133 141 143 211 213 221 223, etc
   in this case grid size could be extended by auxiliary cells in order to be multiple of interleave

   if global_il=0, then interleave works only for fastest z-direction:
   111 113 112 114 121 123 122 124, etc
   if this case grid size aling z-direction must be multiple of interleave
*/
template <class value_tt, size_t interleave=1, class index_t=int>
class UniformGrid {
protected:
  mngptr<value_tt> ptr; /// data array
  Vector_3 dx; /// space step
  Vector_3 pref; /// from where grid starts. wierd construction, delete it
  Vector_3 ds; /// surface area of one elementary cell
  Box b; /// grid is confined within this box
  /// interpolation settings: default
  index_t ngr[3];
  index_t nst[3];
  index_t nen[3];
  /// iterator settings. iterator starts, ends and moves correspondingly to this settings
  index_t igr[3];
  index_t ist[3];
  index_t ien[3];
  /// memory settings
  index_t mgr[3];
  index_t mst[3];
  index_t men[3];
  index_t SZ2; /// mgr[1]*mgr[2]
  index_t GSZ; /// total used memory size
  /// if interleave is global then GSZ  = interleave*nintz (GSZ can be extended to satisfy this this)
  /// if interleave is local, then mgr[2] = interleave*nintz *(mgr[2] must be multiple of interleave)
  index_t nintz;

  int drc; // bit flag specifying axis directions
  Vector_3 ch_drc(const Vector_3 &pos)const{
    Vector_3 p1=b.get_p1(),p2=b.get_p2(),p;
    for(int i=0;i<3;i++){
      if(drc&(1<<i))
        p[i]=p2[i]-(pos[i]-p1[i]);
      else
        p[i]=pos[i];
    }
    return p;
  }
  
  /// sets range to st, en and gr
  /// uses start and end for st and en, or nst and nen, if start = end = NULL
  void set_range(index_t *st, index_t *en, index_t *gr, const index_t *start, const index_t *end);

  /// gets 1d interpolation coefficiencts for projection on dir axis
  virtual void get_1d_coeff(const Vector_3 &p, index_t *ind, vec_type *c, bool nonlocal=false) const;

public:
  typedef value_tt value_t;
  /// Iterator going through all points in the grid range corresponding to nst, nen.
  /// Index change order is Z(2), Y(1), X(0), index goes from small to large
  class iterator{
  protected:
    const UniformGrid *parent;
    index_t ind[3]; /// current grid indices
    int ie; /// if this is end iterator 
    /// grid constructor, constructs sequence start 
    iterator(const UniformGrid *sparent, int end_=0):parent(sparent),ie(end_){
      for(int i=0;i<3;i++)ind[i]=parent->ist[i];
    }
  public:
    typedef index_t difference_t;
    friend class UniformGrid;
    /// default constructor constructs sequence end
    iterator():parent(NULL),ie(1){
      for(int i=0;i<3;i++)ind[i]=0;
    }
    /// copy constructor 
    iterator(const iterator& other):parent(other.parent),ie(other.ie){
      for(int i=0;i<3;i++){
        ind[i]=other.ind[i];
      }
    }
    /// iterator difference
    index_t operator-(const iterator& other) const;

    Vector_3 operator*() const {
      return parent->Position(ind[0],ind[1],ind[2]);
    }

    Vector_3 ds() const{
      return parent->ds;
    }

    Vector_3 dx() const{
      return parent->dx;
    }

    /// index change order:  iz, iy, ix
    iterator& operator++();
     
    /// positive increment to an iterator
    iterator& operator+=(index_t incr);

    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }

    /// complete test
    bool operator!=(const iterator &other) const{
      if(ie!=other.ie)
        return true;
      if(ie==1)
        return false; // ends always equal
      for(int i=0;i<3;i++){
        if(ind[i]!=other.ind[i])
          return true;
      }
      return false;
    }

    /// returns iterator with shifted grid indices
    iterator shift(index_t *sh){
      iterator o=*this;
      for(int i=0;i<3;i++){
        o.ind[i]+=sh[i];
      }
      return o;
    }

    index_t get_end(int i) const {
      return parent->ien[i];
    }

    index_t get_start(int i) const {
      return parent->ist[i];
    }

    index_t get_ind(int i) const {
      return ind[i];
    }
  };
  typedef iterator const_iterator;

  iterator begin() const {
    return iterator(this);
  }

  iterator end() const {
    return iterator(this,1);
  }
  
  UniformGrid():ptr(NULL),drc(0){}

  UniformGrid(const Vector_3 &v1, const Vector_3&v2, const int dir, const index_t *sz, const index_t *start=NULL, const index_t *end=NULL):ptr(NULL){
    init(v1,v2,dir,sz,start,end);
  }

  UniformGrid(const Vector_3 &v1, const Vector_3&v2, const index_t *sz, const index_t *start=NULL, const index_t *end=NULL):ptr(NULL){
    init(v1,v2,0,sz,start,end);
  }

  UniformGrid(const Vector_3 &v1, const Vector_3&v2, const Vector_Nt<index_t,3> &sz, const index_t *start=NULL, const index_t *end=NULL):ptr(NULL){
    init(v1,v2,0,sz.v,start,end);
  }

  /// constructs grid with the same geometrical layout
  template <class T, size_t il>
  UniformGrid(const UniformGrid<T, il, index_t> &other):ptr(NULL){
    init(other.v1,other.v2,0,other.ngr);
  }

  /// v1, v2 are opposite points of the box which confines the grid
  /// sz is array form 3 numbers specifying grid nodes numbers along 3 directions
  /// dir is wierd construction. if dir!=0, then for abs(dir)-1 ds is calculated, and sign of dir corresponds
  /// to ds direction
  /// if start and end are not NULL, then interpolation, iterator and memory ranges will be specified by start and end
  void init(const Vector_3 &v1, const Vector_3&v2, const int dir, const index_t *sz, const index_t *start=NULL, const index_t *end=NULL);

  size_t get_interleave() const {
    return interleave;
  }

  int SetPtr(value_tt *sptr, int managed=0){
    ptr.reset(sptr,managed);
    return 1;
  }

  value_tt *GetPtr() const {
    return ptr.ptr();
  }

  /// sets limits for iterator
  void SetIteratorRange(const index_t *start=NULL, const index_t *end=NULL){
    set_range(ist,ien,igr,start,end);
  }

  /// sets limits for interpolation
  void SetInterpolationRange(const index_t *start=NULL, const index_t *end=NULL){
    set_range(nst,nen,ngr,start,end);
  }

  void GetInterpolationRange(index_t *gr){
    for(int i=0;i<3;i++)
      gr[i]=ngr[i];
  }

  /// set memory block ranges to be addressed
  /// if interleave is global and nonunity, memory size is extended to be a multiple of the interleave 
  /// if interleave is local and nonunity, memory size in direction 2 must be a multiple of the interleave 
  void SetMemoryRange(const index_t *start=NULL, const index_t *end=NULL){
    set_range(mst,men,mgr,start,end);
    SZ2=mgr[1]*mgr[2];
    GSZ=mgr[0]*mgr[1]*mgr[2];
    if(global_il){ // interleave the whole grid
      //int aux_cells=2;
      // maximal plane size
      //int szm=max(mgr[0]*mgr[1],SZ2);
      //szm=max(szm,mgr[0]*mgr[2]);
      //GSZ+=aux_cells*szm;
      index_t rest=GSZ%interleave;
      if(rest)
        GSZ+=interleave-rest;
      nintz=GSZ/interleave;
    }
    else{  // interleave the z-direction only 
      nintz=mgr[2]/interleave;
      if(mgr[2]%interleave)  // this will cause divide by zero error in interpolation to indicate the wrong setting
        nintz=0; 
    }
  }

  /// gets active array size in memory
  size_t msize(int dir=-1) const {
    return dir>=0 ? (size_t)mgr[dir] : (size_t)mgr[0]*mgr[1]*mgr[2];
  }

  /// gets full size in memory (including auxilliary cells)
  size_t size() const {
    return GSZ;
  }

  Vector_3 get_dx() const {
    return dx;
  }

  Vector_3 get_ds() const {
    return ds;
  }

  /// get total grid surface, if grid is 2D in 3D space
  Vector_3 get_total_surface() const {
/*    Vector_3 total_ds=ds;
    for(int i=0;i<3;i++)
      if(!ds[i])total_ds*=ngr[i];
    return total_ds;*/
    return ds*ngr[0]*ngr[1]*ngr[2];
  }

  Box get_region() const {
    return b;
  }

  /// @returns the position of a given grid point
  Vector_3 Position(index_t ix, index_t iy, index_t iz) const {
    return pref+Vector_3(ix*dx[0],iy*dx[1],iz*dx[2]);
  }

  inline index_t pack_base_ind(index_t ix, index_t iy, index_t iz) const {
    return ((ix-mst[0])*mgr[1]+iy-mst[1])*mgr[2]+iz-mst[2];
  }

protected:
  inline index_t pack_base_ind(index_t ix, index_t iy) const {
    return ((ix-mst[0])*mgr[1]+iy-mst[1])*mgr[2];
  }

  inline index_t pack_z_ind(index_t iz) const {
    if(interleave==1 || global_il)
      return iz-mst[2];
    else{      
      iz-=mst[2];
      return iz/nintz+interleave*(iz%nintz);
    }
  }

public:

  /// returns index in memory array for given grid node
  inline index_t pack_ind(index_t ix, index_t iy, index_t iz) const {
    if(global_il && interleave!=1){
      index_t ind=pack_base_ind(ix,iy,iz);
      return ind/nintz+interleave*(ind%nintz);
    }
    else
      return pack_base_ind(ix,iy)+pack_z_ind(iz);
  }

  /// returns grid node for given index in memory array
  inline void unpack_ind(index_t ind, index_t &ix, index_t &iy, index_t &iz) const {
    if(global_il && interleave!=1)
      ind=(ind%interleave)*nintz+ind/interleave;
    ix=ind/SZ2+mst[0];
    ind%=SZ2;
    iy=ind/mgr[2]+mst[1];
    if(interleave==1 || global_il)
      iz=ind%mgr[2]+mst[2];
    else{
      iz=ind%mgr[2];
      iz=(iz%interleave)*nintz+iz/interleave+mst[2];
    }
  }

/*  inline int unpack_z_ind(int iz) const {
    if(interleave!=1 && !global_il)
      return (iz%interleave)*nintz+iz/interleave+mst[2];
    else
      return iz+mst[2];
  }*/

/*  inline int incr_z0_ind(int iz, int sh) const {
    if(interleave!=1){
      iz= (iz%interleave)*nintz+iz/interleave; // unpack
      iz+=sh; // shift
      if(global_il){
        //if(iz>=GSZ)
          //return -1;
        iz=iz%GSZ;
      }
      else
        iz=iz%mgr[2];
      return iz/nintz+interleave*(iz%nintz); // pack
    }
    else
      if(global_il)
        return (iz+sh)%GSZ;
      else 
        return (iz+sh)%mgr[2];
  }*/

  /// gets the loop ranges for two loops: [0, range1) and [range1, nintz)
  /// which must be performed for each interleave bank for
  /// convolution with given positive shift
  /// returns: bank section point: range1 
  ///          bank shift for the first loop (the other one is dbank0+1),
  ///          shifts for the first and the second loop
/*  void get_contiguous_ranges(int shift, int &range1, int &dbank0, int &shift0, int &shift1) const {
    dbank0=shift/nintz;
    shift0=shift%nintz;
    range1=nintz-shift0;
    shift1=-range1;
  }*/

  // what is the meaning of this function?
  index_t get_range1(index_t shift) const {
    return nintz-shift%nintz;
  }

  /// validity check for interpolation index
  bool check_interpolation_ind(index_t ix, index_t iy, index_t iz) const {
    if(ix<nst[0] || ix>nen[0])return false;
    if(iy<nst[1] || iy>nen[1])return false;
    if(iz<nst[2] || iz>nen[2])return false;
    return true;
  }

  /// validity check for memory index
  bool check_memory_ind(index_t ix, index_t iy, index_t iz) const {
    if(ix<mst[0] || ix>men[0])return false;
    if(iy<mst[1] || iy>men[1])return false;
    if(iz<mst[2] || iz>men[2])return false;
    return true;
  }

  value_tt &operator()(index_t ix, index_t iy, index_t iz) const {
    return ptr[pack_ind(ix,iy,iz)];
  }

  value_tt &operator()(iterator it) const {
    return ptr[pack_ind(it.ind[0],it.ind[1],it.ind[2])];
  }

  value_tt operator()(const Vector_3 &place) const{
    return Interpolate(place);
  }

  /** Gets interpolation coefficients (int_ind - index ararys, values - coefficients in interpolation).
   the arrays must be at least 8 elements long.
   if force_external, then point which is outside box will be interpolated as well.
   if nonlocal is not NULL,  
   the negative indicies in int_ind correspond to absent grid points,
   the space location of kth absent point is pushed_back to nonlocal[-int_ind[k]-1].
   returns number of elements in interpolation
   */
  int GetCoeff(const Vector_3 &place, index_t *int_ind, vec_type *values, int force_external=0, vector<Vector_3> *nonlocal=NULL, bool all_nonlocal=false) const;

  /// if point indices are inside interpolation region
  int test_local(const Vector_3 &place) const;

  /// interpolate value using array ptr
  value_tt Interpolate(const Vector_3 &place, int force_external=0) const;

  /// distributes the value being added between nearest grid points 
  /// @returns the number of grid points affected
  int AddValue(const Vector_3 &place, value_tt value, int force_external=0);
};

template<size_t interleave>
Vector_3 get_iterator_dx(const typename UniformGrid<Vector_3, interleave, ptrdiff_t>::iterator &it){
  return it.dx();
}

Vector_3 get_iterator_dx(const UniformGrid<Vector_3, 1, ptrdiff_t>::iterator &it);

///\en returns an integral of function recorded on a grid data over a volume enclosed between p1 and p2
class volume_integral_grid_t{
  UniformGrid<vec_type> grid;
public:
  volume_integral_grid_t(const Vector_3 &p1, const Vector_3 &p2, const int *N, vec_type *ptr, int managed=0){
    grid.init(p1,p2,2,N);
    grid.SetPtr(ptr,managed);
  }
  ///\en returns an integral of some function over a volume enclosed between p1 and p2
  vec_type operator()(const Vector_3 &p1, const Vector_3 &p2) const;

};


template <class value_tt, size_t interleave=1, class index_t=ptrdiff_t>
class NonUniformGrid: public UniformGrid<value_tt,interleave,index_t>{

//  int argtype; /// bit flag for working dimensions

  using UniformGrid<value_tt,interleave,index_t>::ptr;
  using UniformGrid<value_tt,interleave,index_t>::dx;
  using UniformGrid<value_tt,interleave,index_t>::pref;
  using UniformGrid<value_tt,interleave,index_t>::ds;
  using UniformGrid<value_tt,interleave,index_t>::b;
  
  using UniformGrid<value_tt,interleave,index_t>::ngr;
  using UniformGrid<value_tt,interleave,index_t>::nst;
  using UniformGrid<value_tt,interleave,index_t>::nen;
  
  using UniformGrid<value_tt,interleave,index_t>::igr;
  using UniformGrid<value_tt,interleave,index_t>::ist;
  using UniformGrid<value_tt,interleave,index_t>::ien;

  using UniformGrid<value_tt,interleave,index_t>::mgr;
  using UniformGrid<value_tt,interleave,index_t>::mst;
  using UniformGrid<value_tt,interleave,index_t>::men;
  using UniformGrid<value_tt,interleave,index_t>::SZ2;
  using UniformGrid<value_tt,interleave,index_t>::GSZ;
  using UniformGrid<value_tt,interleave,index_t>::nintz;

  vector<value_tt> x[3]; /// node sequences along each direction

  /// initialize x as grid is uniform
  int make_uniform(){
    for(int i=0;i<3;i++){
      vec_type *g1d = new value_tt [ngr[i]];
      for(index_t j=0;j<ngr[i];j++)
        g1d[j]=pref[i]+j*dx[i];
      SetMeshSteps(i,g1d,g1d+ngr[i]);
      delete[]g1d;
    }
    return 1;
  }

  virtual void get_1d_coeff(const Vector_3 &p, index_t *ind, vec_type *c, bool nonlocal=false) const;

public:

  void init(const Vector_3 &v1, const Vector_3&v2, const int dir, const index_t *sz){
    UniformGrid<value_tt,interleave,index_t>::init(v1,v2,dir,sz);
    make_uniform();
  }

  /// initialize nonuniform grid direction
  /// beg, end - iterator on nodes coordinates in this direction
  template<class inp_it>
  int SetMeshSteps(size_t dir, inp_it beg, inp_it end){
    size_t nv=0; // nodes amount
    x[dir].clear();
    for(;beg!=end;++beg,++nv)
      x[dir].push_back(*beg);
    if(nv<=1)
      return -1;

    index_t sz[3]={ngr[0],ngr[1],ngr[2]};
    Vector_3 p1=b.get_p1(),p2=b.get_p2();
    sort(x[dir].begin(),x[dir].end());
    pref[dir]=x[dir][0];
    p1[dir]=x[dir][0],p2[dir]=x[dir][nv-1];
    dx[dir]=(p2[dir]-p1[dir])/(nv-1); // average dx
    sz[dir]=nv;
    UniformGrid<value_tt,interleave,index_t>::init(p1,p2,0,sz);
    return 1;
  }

  Vector_3 Position(index_t ix, index_t iy, index_t iz) const{
    return Vector_3(x[0][ix],x[1][iy],x[2][iz]);
  }

  vec_type GetControlVolume(index_t ix, index_t iy, index_t iz) const{
    vec_type v=1;
    index_t ic[3]={ix,iy,iz};
    for(int i=0;i<3;i++){
      vec_type p1,p2;
      p1=p2=x[i][ic[i]];
      if(ic[i]>0)
        p1=x[i][ic[i]-1];
      if(ic[i]<index_t(x[i].size())-1)
        p2=x[i][ic[i]+1];
      if(p2!=p1)
        v*=(p2-p1)/2;
    }
    return v;
  }

};

template <class value_tt, class index_t=ptrdiff_t>
class InterpBox:  public Box{
public:
  value_tt cube[2][2][2];
  
  /// gets the value by linear index 
  value_tt &operator[](index_t i) const {
    return *((value_tt *)cube+i);
  }

  index_t pack_ind(index_t i, index_t j, index_t k) const {
    return 4*i+2*j+k;
  }

  void unpack_ind(index_t ind, index_t &ix, index_t &iy, index_t &iz) const {
    ix=ind/4;
    iy=(ind-ix*4)/2;
    iz=ind%2;
  }

  InterpBox(const Vector_3 &sp1, const Vector_3 &sp2): Box(sp1,sp2){}

  /// Gets interpolation coefficients and indicies
  /// @return the number of nonzero indicies
  int GetCoeff(const Vector_3 &place, index_t *int_ind, vec_type *values, int force_external=0) const;
};

# if 0

struct GrdData{
 vector<Vector_3> points;
 vector<vec_type> coeffs;
 vector<int> shifts;
 int mind;
};

/// given a set of points, creates a box for finding interpolation coefficients
/// to use in 3-linear interpolation
template <class grid_it> 
int InterGridInterpolation(grid_it beg, grid_it end, const Vector_3 &place){
  vector<int> shifts;
  grid_it it=beg;
  for(;it!=end;++it){
    arr[i].mind=it->GetShifts(place,arr[i].shifts,arr[i].coeffs,1);
    // detecting the points
    int ni=arr[i].coeffs.size();
    for(j=0;j<ni;j++){
      points.push_back(it->get_position(arr[i].mind+arr[i].shifts[j]);
    }
  }
  // building the box
  Vector_3 p1, p2;
  get_min_box(points.begin(), points.end(), place,p1,p2);
  InterpBox ibox(p1,p2);
  int bind[8];
  vec_type bcoeff[8];
  int nb=ibox.GetCoeff(place,bind,bcoeff,1);
  for(i=0;i<nb;i++){
    Vector_3 grdpoint=ibox.get_position(bind[i]);
    /// detecting which grid is better to represent this particular point

  }
}

# endif

/*
Generate 1D nonumiform mesh from given uniform mesh section.
Each uniform mesh can be described by three numbers: x1 (left point), x2 (right point), dx (mesh step),
which can be stored in Vector_3.
Vector_3 rough describes initial uniform mesh with small resolution,
dence correspond to higher resolution uniform meshes which are placed somewhere incide the rough mesh.
mesh_generate tries to generate nonuniform which has uniform high resolution dence.dx between dence.x1 and dence.x2,
uniform small resolution rough.dx far outside dence, 
and nonuniform intermediate resolution rough.dx < dx < dence.dx close to dence 
within intermediate region with approximate size (rough.dx/dence.dx)*rough.dx.
If generetation is successful, than nodes position of resulted nonuniform mesh will be recorded to mesh.
Otherwise function will return -1.
*/
int generate_1d_mesh(Vector_3 rough, vector<Vector_3> dence, vector<vec_type> &mesh);

# endif

