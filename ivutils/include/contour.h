/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.12 $
 *   $Date: 2013/03/05 14:14:44 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/contour.h,v 1.12 2013/03/05 14:14:44 belousov Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/contour.h,v $
$Revision: 1.12 $
$Author: belousov $
$Date: 2013/03/05 14:14:44 $
*/
/*s****************************************************************************
 * $Log: contour.h,v $
 * Revision 1.12  2013/03/05 14:14:44  belousov
 * *** empty log message ***
 *
 * Revision 1.11  2013/02/22 11:01:26  valuev
 * thin surface interface
 *
 * Revision 1.10  2013/02/15 08:08:13  valuev
 * cross edges
 *
 * Revision 1.9  2013/02/14 12:12:16  valuev
 * started adding CrossEdges
 *
 * Revision 1.8  2013/02/02 19:56:29  valuev
 * fixed new vector product
 *
 * Revision 1.7  2013/02/01 18:01:12  valuev
 * added vector volume (needs testing)
 *
 * Revision 1.6  2013/01/29 16:41:32  valuev
 * added faceted surface and rotation quaternions
 *
 * Revision 1.5  2013/01/28 12:26:08  belousov
 * fixed contour compilation under icc
 *
 * Revision 1.4  2013/01/25 15:27:38  valuev
 * restructured contour: got rid of virtual functions and unnecessary storage
 *
 * Revision 1.3  2013/01/25 10:35:35  valuev
 * contour restructure (not working)
 *
 * Revision 1.2  2013/01/24 16:58:46  valuev
 * removed virtual functions from base contour
 *
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.50  2012/09/24 23:10:24  lesha
 * *** empty log message ***
 *
 * Revision 1.49  2012/09/24 19:25:30  lesha
 * documentation
 *
 * Revision 1.48  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.47  2012/03/21 17:01:53  lesha
 * documentation
 *
 * Revision 1.46  2012/02/17 00:10:53  lesha
 * PtrContour is added
 *
*******************************************************************************/
#ifndef CONTOUR_H
#define CONTOUR_H

/// \en @file contour.h \brief Some useful iterators (iterator at points, edges, planes)
/// and classes for 2D and 3D contours.
/// We have different contour implementations that can be used for different purposes

#include "vector_3.h"
#include "plane_3.h"
#include "refobj.h"
#include <vector>

using namespace std;

/// \en records opposite coordinates of maximal rectangular parallelepiped which contains some points
/// defined by point iterators it and end
/// returns average of selected points
template<class point_it,int N>
Vector_Nt<vec_type,N> GetBoundingBox(point_it it, point_it end,Vector_Nt<vec_type,N>* cube1,Vector_Nt<vec_type,N>* cube2);

///\en stores two vertices of the edge
///\ru возращает полезные функции ребра по двум образующим его вершинам
template <class point_it>
class edge_t{
  const point_it p1, p2; // two vertices of the edge
public:

  typedef typename iterator_traits<point_it>::value_type vec_t;

  // sp1, sp2 are vertices of the edge
  edge_t(const point_it &sp1, const point_it &sp2):p1(sp1),p2(sp2){}
  //const vec_t &get_p1() const{return *p1;}
  //const vec_t &get_p2() const{return *p2;}
  
  vec_t get_p1() const {return *p1;}
  vec_t get_p2() const {return *p2;}
  
  vec_t rel() const{
    return *p2-*p1;
  }

  vec_t dir(vec_type *len=NULL) const{
    vec_t v=rel();
    vec_type l=v.normalize();
    if(len)*len=l;
    return v;
  }

  vec_type len() const{
    return (p2-p1).norm();
  }
};

///\en contour edges iterator, based on contour start and end vertices iterators
///\ru итератор по ребрам контура, построенный на основе его начального и конечного итератора по его вершинам
template <class  point_it>
class generic_edge_it{
  point_it beg, end; // start and end iterator at the contour vertex
  point_it cur, cur1; // iterator at two adjoint vertices forming an edge

  void advance_cur1(){
    cur1++;
    if(!(cur1!=end))cur1=beg; // looping the second pointer
  }

public:

  typedef typename iterator_traits<point_it>::value_type vec_t;
  typedef edge_t<point_it> edge;

  // sbeg, send are start and end iterator at the contour vertex,
  // endpoint is 1 if edge iterator is end iterator
  generic_edge_it(point_it sbeg, point_it send, int endpoint=0): beg(sbeg), end(send), cur(sbeg),cur1(sbeg){
    if(!endpoint)advance_cur1();
  }

  generic_edge_it(const generic_edge_it &other): end(other.end), beg(other.beg), cur(other.cur),cur1(other.cur1){}
 
  edge operator*() const {
    return edge(cur,cur1);
  }
  
  generic_edge_it& operator++(){ // prefix
    cur++;
    advance_cur1();
    return *this;
  }

  generic_edge_it operator++(int){ // postfix
    generic_edge_it tmp = *this;
    ++*this;
    return tmp;
  }

  bool operator!=(const generic_edge_it &other) const{
    return cur!=other.cur;
  }
};

///\en iterator at planes normal to the contour plane and passing through edges of the contour
/// Used in ProjectSimplex only.
///\ru итератор по плоскостям, проходящим через ребра контура и нормальным к плоскости контура. 
/// Используется только в ProjectSimplex.
template<class edge_it>
class normplanes_it{
  edge_it cur; // current edge
  Vector_3 n; // normal to the contour
public:

  normplanes_it(const Vector_3 &sn, edge_it beg):cur(beg),n(sn){
    n.normalize();
  }
  normplanes_it(const normplanes_it &other):cur(other.cur),n(other.n){}

  Plane_3 operator*() const {
    typename edge_it::edge cedge=*cur;
    return Plane_3(n%cedge.dir(),cedge.get_p1());
  }
  
  normplanes_it& operator++(){ // prefix
    cur++;
    return *this;
  }

  normplanes_it operator++(int){ // postfix
    normplanes_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const normplanes_it &other) const{
    return cur!=other.cur;
  }
};


///\en Storage class used as template parameter in \ref Contour.
template <class point_it>
struct range_store_t{
  point_it m_beg, m_end;
  ///\en Should have default constructor.
  range_store_t(){}
  ///\en Should have constructor from two \a point_its. 
  range_store_t(point_it beg_, point_it end_): m_beg(beg_), m_end(end_){}
  ///\en \return begin of the range
  point_it begin() const { return m_beg; }
  ///\en \return end of the range
  point_it end() const { return m_end; }
  ///\en \return the size of the store (number of elements)
  //size_t size() const {
  //  return (size_t)(m_end - m_beg); 
  //}
};


/// \en Basic contour interface. Using generic edge iterators helps to construct a contour from any sequence of points.
///     Parameter store_t is used for abstraction from the type of storage for the sequence of points.
template<class point_it, class store_t = range_store_t<point_it>, class edge_it=generic_edge_it<point_it>, 
         class edge_t=typename generic_edge_it<point_it>::edge  >
class Contour{
protected:
  ///\en Knows how to get iterator range.
  store_t store; 
  ///\en Default constructor may be used in derived classes only
  Contour(){}
public:

  typedef typename iterator_traits<point_it>::value_type vector_t;
  typedef point_it point_iterator;
  static const int dimension = vector_t::dimension;
  typedef edge_it edge_iterator;
  typedef edge_t  edge_type;

  Contour(point_it beg_, point_it end_): store(beg_,end_){}

  ///\en Gets begin of point iterator range by calling store_t function
  point_it points_begin() const { return store.begin(); }
  ///\en Gets end of point iterator range by calling store_t function
  point_it points_end() const { return store.end(); }

  edge_iterator edges_begin() const{
    return edge_it(points_begin(),points_end());
  }
  edge_iterator edges_end() const{
    return edge_it(points_end(),points_end(),1); // may be (points_end(),points_begin()) ?
  }

  /*int GetNPoints() const {
    return (int)store.size();
  }

  int GetNEdges()const {
    return GetNPoints()-1;
  }*/

  /// \en gets vertex with given index or infinite vector if index is not found
  vector_t GetPoint(int ind) const {
    point_it it=points_begin(), e=points_end();
    for(int i=0; it!=e; ++it, i++)
      if(i==ind)return *it;
    return vector_t(VEC_INFTY);
  }

  /// \en sets the value for given vertex
  /// returns <0 if this is not possible
  int SetPoint(int ind, const vector_t &vec){
    return -1;
  }

  void Clear(){}

  /// \en calculates the area of the contour (assumed flat)
  vec_type Area() const;

  /// \en Calculates directed area of the contour based from first contour point (for volume calculations).
  ///     Positive direction of the contor is assumed to be anti-clockwise (right-hand rule). 
  ///     Fills the base point if pointer is not NULL.
  typename vector_t::vec_product_t VecArea(vector_t *base=NULL) const{
    typedef typename vector_t::vec_product_t result_t;
    result_t sum = result_t(); // should initialize with 0 for simple types as well
    point_iterator it=points_begin(), e=points_end();
    if(!(it!=e))return sum; // 0 points in the contour?
    vector_t p_base=*it;
    if(base)
      *base=p_base;
    ++it;
    if(!(it!=e))return sum; // 1 point  in the contour?
    vector_t v1=*it-p_base;
    ++it;
    for(;it!=e;++it){
      vector_t v2=*it-p_base;
      sum+=v1%v2;
      v1=v2;
    }
    return sum;
  }

  /// \en gets contour center
  vector_t GetCenter() const;

  vector_t GetBoundingBox(vector_t *v1, vector_t *v2) const{
    return ::GetBoundingBox(points_begin(),points_end(),v1,v2);
  }
};

template<class point_it, class vector_tt=typename iterator_traits<point_it>::value_type, class store_t = range_store_t<point_it>, 
  class edge_it=generic_edge_it<point_it>, 
  class edge_t=typename generic_edge_it<point_it>::edge>
class Contour_N: public Contour<point_it, store_t, edge_it, edge_t>{
protected:
  Contour_N(){}
public:
  typedef Contour<point_it, store_t, edge_it, edge_t> base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  static const int dimension = base_t::dimension;
  typedef typename base_t::edge_iterator edge_iterator;
  typedef typename base_t::edge_type edge_type;

  Contour_N(point_it beg_, point_it end_): Contour<point_it,store_t,edge_it,edge_t>(beg_, end_) {}

};

/// \en Contour_N specialization for 2D
template<class point_it,  class store_t, class edge_it, class edge_t>
class Contour_N<point_it,Vector_2,store_t,edge_it,edge_t>: public Contour<point_it, store_t, edge_it,edge_t>{
protected:
  Contour_N(){}
public:
  typedef Contour<point_it, store_t, edge_it,edge_t> base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  static const int dimension = base_t::dimension;
  typedef typename base_t::edge_iterator edge_iterator;
  typedef typename base_t::edge_type edge_type;

  Contour_N(point_it beg_, point_it end_): base_t(beg_, end_) {}

  bool TestPoint(const Vector_2 p) const;
};

/// \en  Contour_N specialization for 3D
template<class point_it,class store_t, class edge_it, class edge_t>
class Contour_N<point_it,Vector_3,store_t, edge_it,edge_t>: public Contour<point_it,store_t,edge_it,edge_t>{
protected:
  Contour_N(){}
public:


  typedef Contour<point_it,store_t,edge_it,edge_t> base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  static const int dimension = base_t::dimension;
  typedef typename base_t::edge_iterator edge_iterator;
  typedef typename base_t::edge_type edge_type;
  ///\en Iterator of planes containing edges and normal to contour plane
  typedef normplanes_it<edge_it> plane_it;

  Contour_N(point_it beg_, point_it end_): base_t(beg_, end_) {}


  using base_t::edges_begin;
  using base_t::edges_end;

  /// \en tests whether a segment crosses the contour plane inside the contour 
  /// and returns the crossing point (if cross!=NULL)
  bool TestEdge(const Vector_3 &v1, const Vector_3 &v2, Vector_3 *cross=NULL) const;

  /// \en gets begin iterator at normal planes
  plane_it planes_begin() const {
    Plane_3 cpl=GetPlane();
    Vector_3 n=cpl.get_normal();
    return plane_it(n,edges_begin());
  }
  /// \en gets end interator at normal planes
  plane_it planes_end() const {
    return plane_it(Vector_3(0),edges_end());
  }

  /// \en gets the plane of the contour (uses 3 first points only)
  Plane_3 GetPlane() const {
    point_it it=base_t::points_begin();
    return Plane_3(GetNormVect(),*it);
  }

  /// gets normal vector (uses 3 first points only)
  Vector_3 GetNormVect() const;
};


/// \en Contour based on dynamic Vector_Nt array
template <int N=3>
class PtrContour: public Contour_N< Vector_Nt<vec_type,N>* >{
public:

  typedef Contour_N< Vector_Nt<vec_type,N>* > base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  static const int dimension = base_t::dimension;
  typedef typename base_t::edge_iterator edge_iterator;
  typedef typename base_t::edge_type edge_type;


  PtrContour(vector_t * points_=NULL, int num_=0):
    base_t(points_, points_+num_){}


  /// no checks
  vector_t GetPoint(int ind) const {
    //return points[ind];
    return (base_t::store.m_beg)[ind];
  }

  /// checks
  int SetPoint(int i,const vector_t &p){
    if(i<0 || i>=(int)(base_t::store.m_beg-base_t::store.m_end))
      return -1;
    (base_t::store.m_beg)[i]=p;
    return 1;
  }

  int GetNPoints() const {
    return (int)(base_t::store.m_beg-base_t::store.m_end);
  }
};

///\en Storage class for array storage 
template <class T, int num>
struct array_store_t{
  typedef T value_type;
  T data[num];
  ///\en arguments (needed for base contour constructor) are ignored
  array_store_t(T* beg_=NULL, T* end_=NULL){}
  T* begin() const { return (T*)data; }
  T* end() const { return (T*)data+num; }
};

/// \en Contour based on Vector_Nt array.
template <int N, int num>
class ArrayContour: public Contour_N< Vector_Nt<vec_type,N>* ,
                                      Vector_Nt<vec_type,N>,
                                      array_store_t<Vector_Nt<vec_type,N>, num > >{


public:

  typedef Contour_N< Vector_Nt<vec_type,N>* ,
                                      Vector_Nt<vec_type,N>,
                                      array_store_t<Vector_Nt<vec_type,N>, num > > base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  static const int dimension = base_t::dimension;
  typedef typename base_t::edge_iterator edge_iterator;
  typedef typename base_t::edge_type edge_type;
 
  ArrayContour(): base_t(NULL,NULL) {}

  /// no checks
  vector_t GetPoint(int ind) const {
    return base_t::store.data[ind];
  }

  /// checks
  int SetPoint(int i,const vector_t &p){
    if(i<0 || i>=num)
      return -1;
    base_t::store.data[i]=p;
    return 1;
  }

  int GetNPoints() const {
    return num;
  }

  const vector_t* GetPoints() const {
    return base_t::store.data;
  }

  vector_t *GetPoints(){
    return base_t::store.data;
  }

};

///\en Storage class for vector storage
template <class T>
struct vec_store_t{
  typedef T value_type;
  std::vector<T> data;
  ///\en arguments (needed for base contour constructor) are ignored
  vec_store_t(typename std::vector<T>::const_iterator beg_, typename  std::vector<T>::const_iterator end_){}
  vec_store_t(){}
  typename std::vector<T>::const_iterator begin() const { return data.begin(); }
  typename std::vector<T>::const_iterator end() const { return data.end(); }
};

/// \en Contour based on Vector_Nt vector.
template<int N=3>
class VecContour: public Contour_N<typename std::vector<Vector_Nt<vec_type,N> >::const_iterator,
                                   Vector_Nt<vec_type,N>,
                                   vec_store_t<Vector_Nt<vec_type,N> > >{

public:

  typedef Contour_N<typename std::vector<Vector_Nt<vec_type,N> >::const_iterator,
                                   Vector_Nt<vec_type,N>,
                                   vec_store_t<Vector_Nt<vec_type,N> > >  base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  static const int dimension = base_t::dimension;
  typedef typename base_t::edge_iterator edge_iterator;
  typedef typename base_t::edge_type edge_type;

  VecContour(int n=0){
    base_t::store.data.resize(n);
  }

  VecContour(const vector<vector_t> &pvec){
    base_t::store.data=pvec;
  }  

  /// constructs a contour from input point iterator, applies additional transform p=>p*a+shift
  template<class InpIt>
  VecContour(InpIt pbeg, InpIt pend, vec_type a=1.,const vector_t &shift=vector_t()){
    for(;pbeg!=pend;++pbeg){
      base_t::store.data.push_back((*pbeg)*a+shift);
    }   
  }

  void add(const vector_t &vect) {
    base_t::store.data.push_back(vect);
  }


  /// \en sets the ith point, if i is greated than current size
  /// appends the vector with the required amount of default Vector_3
  int SetPoint(int i,const vector_t &p){
    if(i<0)return -1;
    int n=(int)base_t::store.data.size();
    while(i>n){
      base_t::store.data.push_back(vector_t());
      n++;
    }
    if(i==n){
      base_t::store.data.push_back(p);
    }
    else{
      base_t::store.data[i]=p;
    }
    return 1;
  }

  /// no checks
  vector_t GetPoint(int ind) const {
    return base_t::store.data[ind];
  }

  int GetNPoints() const {
    return (int)base_t::store.data.size();
  }

  const vector<vector_t> &GetPoints() const {
    return base_t::store.data;
  }

  vector<vector_t> &GetPoints(){
    return base_t::store.data;
  }

  void Clear(){
    base_t::store.data.clear();
  }

  operator PtrContour<N>()const{
    return PtrContour<N>((vector_t *)&(base_t::store.data[0]),GetNPoints());
  }
};

/*template<class contour_t>
PtrContour<contour_t::dimension> make_PtrContour(const contour_t &cnt){
  int n=cnt.GetNPoints();
  contour_t::vector_t *points = n ? NULL : new contour_t::vector_t[n];
  n=0;
  for(contour_t::point_iterator it=cnt.points_begin(),e=cnt.points_end();it!=e;++it)
    points[n++]=*it;
  return PtrContour<contour_t::dimension>(make_mngarg(points),n);
}

template<int N,int num>
PtrContour<N> make_PtrContour(const ArrayContour<N,num> &cnt){
  return PtrContour<N>(cnt.points,num);
}*/

/// \en projects Vector_3 on 2D coordinate system defined by its origin and axes x, y
inline Vector_2 Vector_3to2(const Vector_3 &v3, const Vector_3 &origin, const Vector_3 &x, const Vector_3 &y){
  Vector_3 tmp=v3-origin;
  return Vector_2(tmp*x, tmp*y);
}

/// \en opposite to Vector_3to2
inline Vector_3 Vector_2to3(const Vector_2 &v2, const Vector_3 &origin, const Vector_3 &x, const Vector_3 &y){
  return origin+v2[0]*x+v2[1]*y;
}

/// \en 2D points iterator as a projection of given 3D points iterator
/// on 2D coordinate system defined by its origin, and axes x, y
template<class point_it>
class contour_projection_it{
  point_it it;
  Vector_3 o,x,y;
public:
  typedef forward_iterator_tag iterator_category;
	typedef void difference_type;
	typedef difference_type distance_type;	// retained
	typedef typename iterator_traits<point_it>::pointer pointer;
	typedef typename iterator_traits<point_it>::reference reference;


  contour_projection_it(){}
  contour_projection_it(point_it it_, const Vector_3 &o_, const Vector_3 &x_, const Vector_3 &y_):it(it_),o(o_),x(x_),y(y_){}

  typedef Vector_2 value_type;
   
  Vector_2 operator*() const {
    return Vector_3to2(*it,o,x,y);
  }
  
  contour_projection_it &operator++(){ // prefix
    ++it;
    return *this;
  }

  contour_projection_it operator++(int){
    contour_projection_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const contour_projection_it &other) const{
    return it!=other.it;
  }
};

/// \en 2D contour as a projection of given 3D contour 
/// on 2D coordinate system defined by its origin, and axes x, y
template<class contour_t>
class Contour_3to2: public Contour_N<contour_projection_it<typename contour_t::point_iterator> >{

public:

  typedef Contour_N<contour_projection_it<typename contour_t::point_iterator> > base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  static const int dimension = base_t::dimension;
  typedef typename base_t::edge_iterator edge_iterator;
  typedef typename base_t::edge_type edge_type;

  Contour_3to2(){}
  Contour_3to2(const contour_t &cnt,const Vector_3 &o,const Vector_3 &x,const Vector_3 &y):
    base_t(point_iterator(cnt.points_begin(),o,x,y),point_iterator(cnt.points_end(),o,x,y)){}

};

///\en record to the contour cnt section of the plane plane0 by polyhedron, 
/// formed by planes iterated from iterator beg till end
///\ru записывает в контур cnt сечение плоскости plane0 многогранником, 
/// образованным плоскостями, по которым проходят итераторы beg, end
template<class plane_it, class contour_t>
int ProjectSimplex(const Plane_3& plane0, plane_it beg, plane_it end, contour_t &cnt, int nplanes=-1);

/// \en gets all faces of the polyhedral region defined by planes from [beg, end) of type plane_it
/// puts the faces as planar contours into container cont using push_back(...)
/// @return the number of faces
template<class out_cont_t, class plane_it>
int GetFaces(out_cont_t &cont,plane_it beg,plane_it end);

/*
template<class contour_t1, class contour_t2>
size_t GetCrossEdges(const contour_t1 &cont1, const contour_t2 &cont2, vector<pair<Vector_3, Vector_3> > *edges = NULL);


template<class contour_t1, class contour_t2>
size_t GetCrossEdges(const contour_t1 &cont1, const contour_t2 &cont2, vector<pair<Vector_3, Vector_3> > *edges = NULL){
  typename contour_t1::point_iterator pi1 = cont1.points_begin(), pi1e = cont1.points_end();
  // finding first contour normal plane 
  Vector_3 p[3], normal1;
  size_t i=0, limit = 3;
  while(pi1!=pi1e && i<limit){
    p[i++]=*pi1++;
    if(i==limit-1){
      normal1 = (p[i]-p[i-1])%(p[i-2]-p[i-1]);
      if(normal1.minabscoord()<VEC_ZERO)
        limit++;
    }
  }
  if(i<limit)
    return 0.; // first contour has less than 3 points or all points are on a single line
  pi1 = cont1.points_begin();
  Plane_3 nplane1(normal1,*pi1); // normal plane for the 1st contour
  typename contour_t2::point_iterator pi2 = cont2.points_begin(), pi2e = cont2.points_end();
  if(pi2 == pie2) // second contour has 0 points
    return 0.;
  Vector_3 v1= *pi2;
  size_t ncross =0;
  vector< pair<Vector_3, Vector_3> > m_edges;
  for(++pi2;pi2!=pi2e;++pi2){
    Vector_3 v2= *pi2, vcross;
    if( TestEdge(v1,v2,&vcross) ){
      ncross++;
      if(ncross==2){ // record edge
        m_edges.push_back(make_pair(cross_prev,vcross));
        ncross = 0;
      } 
      else
        cross_prev = vcross;
      v1= v2;
    }
  }
  // testing edges against contour1
  Vector_3 len;
  if(m_edges.size()){
    normplanes_it<contour1_t::edge_iterator> plb(normal1,cont1.edges_begin()), ple(normal1,cont1.edges_end());
    Polyhedron<plane_it> poly(plb,ple);
    for(size_t i=0;i<m_edges.size();i++){
      vec_type frac[2];
      Vector_3 dir = m_edges[i].second-edges[i].first;
      int res=TestLine(m_edges[i].first,dir,frac,NULL,NULL);
      if(res==2){
        frac[0]=max(0.,frac[0]);
        frac[1]=min(1.,frac[1]);
        Vector_3 e1 = m_edges[i].first+frac[0]*dir;
        Vector_3 e2 = m_edges[i].first+frac[1]*dir;
        if(edges)
          edges->push_bak(make_pair(e1,e2));
        len+=e2-e1;
      }
    }
  }
  return len.norm();
}
*/




#endif
