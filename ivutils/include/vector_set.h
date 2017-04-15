/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 1995-2012        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   $Revision: 1.6 $
 *   $Date: 2013/12/19 05:45:06 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/vector_set.h,v 1.6 2013/12/19 05:45:06 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/vector_set.h,v $
$Revision: 1.6 $
$Author: lesha $
$Date: 2013/12/19 05:45:06 $
*/
/*s****************************************************************************
 * $Log: vector_set.h,v $
 * Revision 1.6  2013/12/19 05:45:06  lesha
 * load for detectors
 *
 * Revision 1.5  2013/07/23 08:55:32  valuev
 * comments
 *
 * Revision 1.4  2013/06/29 18:49:21  lesha
 * bug is fixed
 *
 * Revision 1.3  2013/06/29 01:55:07  lesha
 * SpaceVectorSet and confining detectors are added
 *
 * Revision 1.2  2013/06/28 02:37:28  lesha
 * spatial integration is added to Project (emValue is modified)
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.37  2012/10/16 15:27:40  lesha
 * filter_vset is moved here
 *
 * Revision 1.36  2012/09/25 01:12:46  lesha
 * documentation
 *
 * Revision 1.35  2012/09/09 16:28:05  lesha
 * get_total_surface is added to surface_set
 *
 * Revision 1.34  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.33  2012/04/13 08:13:21  valuev
 * reverted to M_PI
 *
 * Revision 1.32  2012/04/12 20:11:24  lesha
 * documentation
 *
 * Revision 1.31  2012/03/29 20:17:50  valuev
 * fixed memory leak in base_vset
 *
 * Revision 1.30  2012/03/29 11:35:09  valuev
 * universal scatter_data() for emSourceWave
 *
 * Revision 1.29  2012/03/22 07:31:20  valuev
 * updated increments for iterators, prepared parallel/universal near-to-far (not working)
 *
*******************************************************************************/
#ifndef _VSET_H
#define _VSET_H

/// @file vector_set.h \brief Sequences of vectors. 
/// These sequences can be used in numerical applications to specify 
/// at which points in space fields should be recorded.

/*
 Here there are definition of some classes that can be used
 to store and iterate sequence of vectors and surface elements.
 Surface element is vector normal to some surface, length of this vector is equal to surface area.
 Sequences of vectors can be used in numerical applications to specify 
 at which points in space fields should be recorded.
 Surface elements can be used for integration of some field over chosen surface.
 Basic purpose of such class is to return iterator on start and end element of the sequence.
 Iterators can return corresponding position and surface element.
 SpaceVectorSet is a class which stores all elements of the sequence.
 There are other classes which are intended for some some special cases 
 (like sequence at plane or cylindrical surface)
 where this sequence of vectors can be stored more efficiently.
*/

#include "grid.h"
#include "data_flags.h"

template<class it_t>
Vector_3 get_iterator_dx(const it_t &it){
  return 0;
}

///en This class stores all elements of the sequence.
struct SpaceVectorSet{

  vector<Vector_3> pos; // vector sequence
  vector<Vector_3> dpos; // volume (surface) elements sequence

  class iterator{ // iterator on current element of the sequence
    friend struct SpaceVectorSet;
    const SpaceVectorSet *parent;
    int  i; // current number

    // if e=1 returns the end iterator, otherwise returns iterator on the first element
    iterator(const SpaceVectorSet *sparent,int e=0): parent(sparent){
      i = e ? (int)parent->pos.size() : 0;
    }
  public:

    typedef Vector_3 value_t;

    Vector_3 operator*() const { // returns current position
      return parent->pos[i];
    }

    Vector_3 ds() const{ // returns current surface element
      return parent->dpos.size() ? parent->dpos[i] : 0;
    }

    Vector_3 dx() const{ // returns current volume element
      return parent->dpos.size() ? parent->dpos[i] : 0;
    }

    iterator& operator++(){
      i++;
      return *this;
    }
     
    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }

    bool operator!=(const iterator &other) const{
      return i!=other.i;
    }
  };

  typedef iterator const_iterator;

  SpaceVectorSet(){}

  SpaceVectorSet(const vector<Vector_3> &pos_, const vector<Vector_3> dpos_): pos(pos_), dpos(dpos_){}

  iterator begin() const {
    return iterator(this);
  }

  iterator end() const {
    return iterator(this,1);
  }

  size_t size(){ // size of sequence
    return pos.size();
  }

  // returns sum of oriented surfaces (integral dpos)
  // can be used for calculation of a plane wave flux over total surface
  Vector_3 get_total_surface(){
    Vector_3 surface;
    for(size_t i=0;i<dpos.size();i++)
      surface+=dpos[i];
    return surface;
  }
};

Vector_3 get_iterator_dx(const SpaceVectorSet::iterator &it);

///\en Sequence corresponding to the grid located at some sides of the cube (maximal number of sides is 6).
///\ru Последовательность векторов, пробегающая сетку, натянутую на избранные грани куба.
class BoxSurfaceSet {

  typedef UniformGrid<Vector_3>::iterator grid_it;

  UniformGrid<Vector_3> grd[6]; // grids for 6 sides of the box
  vec_type ds[3]; // surface element area for pairs of opposite grids perpendicular to x, y and z directions
  int used; // which sides are used

public:

  class iterator{
    friend class BoxSurfaceSet;
    const BoxSurfaceSet *parent;
    int gi; // number of current grid
    grid_it it; // curent grid iterator

    iterator(const BoxSurfaceSet *sparent, int e=0);
  public:
    typedef Vector_3 value_t;

    Vector_3 operator*() const {
      return *it;
    }


    Vector_3 ds() const;

    iterator& operator++();
     
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    }

    /// positive increment to an iterator
    iterator& operator+=(size_t incr){
      for(size_t i=0; gi<6 && i<incr;i++)
        ++*this;
      return *this;
    }

    /// complete test
    bool operator!=(const iterator &other) const{
      return gi!=other.gi ? true : gi==6 ? false : it!=other.it;
    }
  };

  typedef iterator const_iterator;

  BoxSurfaceSet() {}

  ///\en
  /// v1, v2 are opposite cube vertices
  /// sz is array of 3 numbers which define grid points number along directions x, y and z
  /// used is bit flag (see BOX_SIDES) which defines used sides (all sides are used by default)
  ///\ru
  /// v1, v2 - противоположные вершины куба, ограничивающего сетку, 
  /// sz - массив из трех чисел, указующих на количество шагов сетки по трем направлениям
  /// used - флаг, отвечающий за грани, на которых натянута сетка
  BoxSurfaceSet(const Vector_3 &v1, const Vector_3&v2, const int *sz, const int used=0xff) {
    init(v1, v2, sz, used);
  }

  void init(const Vector_3 &v1, const Vector_3&v2, const int *sz, const int sused=0xff);

  iterator begin() const {
    return iterator(this);
  }

  iterator end() const {
    return iterator(this,1);
  }

  size_t size();

  // returns sum of oriented surfaces for all used sides (integral ds)
  // f. e. if all surfaces are used, returns 0.
  // can be used for calculation of a plane wave flux over total surface
  Vector_3 get_total_surface();
};

///\en Sequence corresponding to the grid at the cylinder surface.
///\ru Последовательность векторов, пробегающая сетку, натянутую на цилиндрическую поверхность.
class CylinderSurfaceSet {

  Vector_3 origin; // cylinder origin
  Vector_3 n; // direction of cylinder axis
  Vector_3 x, y; // two directions perpendicular to cylinder axis. radial angle is measured in system (x, y)
  vec_type R; // radius
  int znum,finum; // grid points number along axial and radial directions
  vec_type ds; // surface element area (positive if normal to surface is outside in negative otherwise)

public:

  class iterator{
    friend class CylinderSurfaceSet;
    const CylinderSurfaceSet *parent;
    int zi, fi; // axial and radial number of current grid point

    iterator(const CylinderSurfaceSet *sparent,int e=0): parent(sparent), fi(0), zi (e ? parent->znum : 0){}

  public:
    typedef Vector_3 value_t;

    Vector_3 operator*() const;

    Vector_3 ds() const;

    iterator& operator++(){
      fi++;
      if(fi>=parent->finum){
        fi=0;
        zi++;
      }
      return *this;
    }
     
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    }

    /// tests for the end only
    bool operator!=(const iterator &other) const{
      return zi!=other.zi;
    }
  };

  typedef iterator const_iterator;

  CylinderSurfaceSet(){}

  CylinderSurfaceSet(const Vector_3 &origin, const Vector_3 &n, const Vector_3 &x, const Vector_3 &y, 
  vec_type R, const int *sz, int surf_dir=1) {
    init(origin, n, x, y, R, sz, surf_dir);
  }

  // origin is cylinder origin, n is direction of cylinder axis, 
  // x, y are two directions perpendicular to cylinder axis. radial angle is measured in system (x, y)
  // R is radius, znum, finum are grid points number along axial and radial directions.
  // surf_dir defines direction of the normal to cylinder (outside if positive, inside if negative)
  void init(const Vector_3 &origin_, const Vector_3 &n_, const Vector_3 &x_, const Vector_3 &y_, 
  vec_type R_, const int *sz, int surf_dir=1);

  iterator begin() const {
    return iterator(this);
  }

  iterator end() const {
    return iterator(this,1);
  }

  size_t size(){
    return znum*finum;
  }
};

///\en sequence corresponding to the grid at the sphere surface.
/// grid has equal step along azimutal and radial angles
class SphereSurfaceSet {
  Vector_3 origin; // center fo the sphere
  vec_type R; // radius
  vec_type theta0, fi0; // start azimutal and radial angles
  vec_type dtheta, dfi; // azimutal and radial angle grid steps
  int theta_num, fi_num; // azimutal and radial grid points numbers

public:

  class iterator{

    friend class SphereSurfaceSet;
    const SphereSurfaceSet *parent;
    int th, f;
    iterator(const SphereSurfaceSet *sparent,int e=0): parent(sparent), f(0), th(e ? parent->theta_num : 0) {}

  public:
    typedef Vector_3 value_t;

    Vector_3 operator*() const;

    // returns position coordinates in sperical system (R, theta, phi)
    Vector_3 internal_coords() const;

    iterator& operator++();
     
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    }

    /// tests for the end only
    bool operator!=(const iterator &other) const{
      return th!=other.th;
    }
  };

  typedef iterator const_iterator;

  SphereSurfaceSet(): theta_num(0), fi_num(0) {}

  // origin is sphere center, R is radius,
  // theta (fi)_num are azimutal and radial grid points numbers
  // theta (fi) 0 (1) are start and end azimutal and radial angles.
  SphereSurfaceSet(const Vector_3 &origin, vec_type R, int theta_num, int fi_num, 
    vec_type theta0=0, vec_type fi0=0, vec_type theta1=M_PI, vec_type fi1=2*M_PI) {
    init(origin, R, theta_num, fi_num, theta0, fi0, theta1, fi1);
  }

  int init(const Vector_3 &origin_, vec_type R_, int theta_num_, int fi_num_, 
    vec_type theta0_=0, vec_type fi0_=0, vec_type theta1=M_PI, vec_type fi1=2*M_PI);

  iterator begin() const {
    return iterator(this);
  }

  iterator end() const {
    return iterator(this,1);
  }

  size_t size() {
    return theta_num*fi_num;
  }
};

#if 0
template<class vset_t>
class CuttedVectorSet{
  mngptr<vset_t> vset;
  int_pack slices;
  CuttedVectorSet(mngarg<vset_t> vset_): vset(vset_){}
public:
  class iterator{
    friend class SliceVectorSet;
    typename vset_t::iterator vset_it;
    int_pack::iterator it, e;
    int cur;
    
    iterator(typename vset_t::iterator vset_it_, int_pack::iterator it_, int_pack::iterator e_): 
    vset_it(vset_it_), it(it_), e(e_), cur(0) {
      iterate();
    }
    void iterate(){
      while(it!=e){
        int n=*it;
        int num=n/2;
        if(n%2){
          for(int i=0;i<num;i++)
            ++vset_it;
          ++it;
        }
        else{
          cur=num;
          break;
        }
      }
    }

  public:
    typedef Vector_3 value_t;

    Vector_3 operator*() const {
      return *vset_it;
    }

    Vector_3 ds() const{
      return vset_it.ds();
    }

    iterator& operator++(){
      if(cur){
        ++vset;
        cur--;
      }
      else{
        ++it;
        iterate();
      }
    }
     
    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }

    /// tests for the end only
    bool operator!=(const iterator &other) const{
      return it!=other.it;
    }
  };
  void start_record(){
    slices.start_record();
  }
  void next_record(int it){
    slices.next_record(it>0?it*2:-2*it-1);
  }
  void end_record(){
    slices.end_record();
  }
  iterator begin(){
    return iterator(vset->begin(),slices.begin(),slices.end());
  }
  iterator end(){
    return iterator(vset->end(),slices.end(),slices.end());
  }
};
#endif

// returns iterator position coordinates in specific coordinate system
template <class iter>
Vector_3 internal_coords(iter it) {
  return *it;
}

// returns vector coordinates in specific coordinates system  
template <class iter>
Vector_3 internal_coords(const Vector_3 &F, iter it) {
  return F;
}

// returns iterator position coordinates in spherical system (r, phi, theta)
Vector_3 internal_coords(SphereSurfaceSet::iterator it);

// returns vector coordinates in sperical system (r, phi, theta) associated with iterator position,
// see http://fdtd.kintechlab.com/en/fdtd#the_amplitude_scattering_matrix
Vector_3 internal_coords(const Vector_3 &F, SphereSurfaceSet::iterator it);

// if some functions or classes are using vector sequence they are templates from vector sequence.
// it can lead to to many instantiations during compilation.
// if vector sequence iterator operations are not time consuming compare to other operations of the function,
// it is more efficient to use virtual operations.

// virtualizator of iterator at vector sequence
class base_vset_it{
public:
  virtual Vector_3 operator*() const=0;
  virtual void operator++()=0;
  virtual base_vset_it *copy()=0;
  virtual ~base_vset_it(){}
};

// inheritted class for iterator at some chosen sequence vset_it_t
template<class vset_it_t>
class virt_vset_it: public base_vset_it{
  vset_it_t it;
public:
  virt_vset_it(vset_it_t it_): it(it_){}
  virtual Vector_3 operator*() const{return *it;}
  virtual void operator++(){++it;}
  virtual base_vset_it *copy(){
    return new virt_vset_it(it);
  }
};

class base_surf_vset_it{
public:
  virtual Vector_3 operator*() const=0;
  virtual void operator++()=0;
  virtual Vector_3 ds() const=0;
  virtual base_surf_vset_it *copy()=0;
  virtual ~base_surf_vset_it(){}
};

template<class vset_it_t>
class virt_surf_vset_it: public base_surf_vset_it{
  vset_it_t it;
public:
  virt_surf_vset_it(vset_it_t it_): it(it_){}
  virtual Vector_3 operator*() const{return *it;}
  virtual Vector_3 ds() const{return it.ds();}
  virtual void operator++(){++it;}
  virtual base_surf_vset_it *copy(){
    return new virt_surf_vset_it(it);
  }
};

// same as \ref  virt_vset_it, additionally transforming the set to a new shape
template<class vset_it_t>
class transformed_vset_it: public base_surf_vset_it{
  vset_it_t it;
  mngptr<VecTransform> tr;
  mngptr<VecTransform> tr_ds;
public:
  ///\en Constructs the iterator
  transformed_vset_it(vset_it_t it_, mngarg<VecTransform> transform, mngarg<VecTransform> transform_ds ): it(it_),tr(transform), tr_ds(transform_ds){}
  virtual Vector_3 operator*() const{return (*tr)(*it);}
  virtual Vector_3 ds() const{return (*tr_ds)(it.ds());}
  virtual void operator++(){++it;}
  virtual base_surf_vset_it *copy(){
    return new transformed_vset_it(it,make_mngarg(tr->copy()),make_mngarg(tr_ds->copy()));
  }
};

template<class vset_it_t>
base_surf_vset_it *make_virt_vset_it(vset_it_t it, mngarg<VecTransform> transform =NULL, mngarg<VecTransform> transform_ds=NULL ){
  if(transform.ptr() && transform_ds.ptr())
    return new transformed_vset_it<vset_it_t>(it,transform,transform_ds);
  else
    return new virt_surf_vset_it<vset_it_t>(it);
}

// select points of given vector set which belong to some container and confined region
// and record them (and their surface elements) to vectors vect and ds
template<class cont_t,class vset_t>
void filter_vset(cont_t *cont, Region_3 *conf, const vset_t &vset, vector<Vector_3> &vect, vector<Vector_3> *ds=NULL){
  for(typename vset_t::iterator it=vset.begin(),e=vset.end();it!=e;++it){
    Vector_3 pos=*it;
    if(conf && !conf->TestPoint(pos))
      continue;
    if(cont && !cont->TestPoint(pos))
      continue;
    vect.push_back(pos);
    if(ds)
      ds->push_back(it.ds());
  }
}

#endif
