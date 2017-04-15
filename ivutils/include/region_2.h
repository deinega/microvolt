#ifndef _REGION_2_H
#define _REGION_2_H

/// \ru @file region_2.h \brief Classes for various 2D or 3D bodies.
/// \en @file region_2.h \brief Классы описывающие различные 2D или 3D геометрические тела.

#include "contour.h"
#include "basis_3.h"
#include "refobj.h"
#include "math_utils.h"

///\en Bit flags to specify different visualization formats used in \ref CreateDumper 
enum VISUAL_FORMATS {
  GNUPLOT = 0x1,
  VTK = 0x2,
};

///\en Specifies default visual dump format.
///    Class RegDumper will be implemented correspondingly to this format.
///\ru В зависимости от различных значений параметра DUMP будут подключаться 
/// различные реализации класса RegDumper и функций CreateDumper
#define DUMP GNUPLOT

///\en Class responsible for drawing geometrical objects
/// It should be defined in other file
///\ru Класс, отвечающий за рисование объектов.
/// Переопределяется в отдельно подкючаемых файлах
template<int N>
class RegDumper;

template<class reg_t>
RegDumper<reg_t::dimension> *GetRegionDumper(const reg_t *reg);

///\en Body categories \{
///\en Generic category.
struct unknown_cat{};
///\en Convex linear figure (3D polyhedron).
///    Must have subtype plane_it and functions planes_begin(),
///    planes_end to get all composing planes.
struct polyhedron_cat{};
///\en Sphere.
struct sphere_cat{};
///\}

///\en Bit flags specifying the existing implementations
///    of some \ref Region<> functions, returned by Region::GetFlags(). 
enum RegFlags{
 HAS_TESTPOINT=0x1, ///<\en TestPoint function is implemented for this region
 HAS_TESTLINE=0x2, ///<\en TestLine function (and also TestRay and TestEdge) is implemented for this region
 HAS_TESTCONTOUR=0x4, ///<\en TestContour function is implemented for this region
 HAS_VOLUME=0x8, ///<\en Volume function is implemented for this region
 HAS_CLONE=0x10, ///<\en Clone function is implemented for this region
 NOT_CONVEX=0x20 ///<\en Arbitrary shape
};

///\en Convex n-dimensional body. May have infinite extension in some directions.
///\ru Выпуклое n-мерное тело. Может быть неограниченным в некоторых направлениях.
template<int N>
class Region{
public:
  ///\en Category of the region to determine the correct template behaviour for 
  ///    various algorithms (see \ref Body categories). 
  typedef unknown_cat category;

  typedef Vector_Nt<vec_type,N> vector_t;
  static const int dimension=N;

  ///\en Returns a set of bit flags specifying which functions
  ///    are implemented for this region (see \ref RegFlags).
  ///    \return HAS_TESTPOINT|HAS_CLONE
  virtual int GetFlags() const {
    return HAS_TESTPOINT|HAS_CLONE;
  }

  ///\en Returns true if the point is in the the region.
  virtual bool TestPoint(const vector_t &p) const{
    return true;
  }

  ///\en Records the fractions frac1 and frac2 of the line (p,dir) corresponding to the intersections
  /// surfpi=p+fraci*dir, surfni - coords of intersection and corresponding outer normal.
  /// Intersections with frac<=epsilon are ignored.
  /// Returns the number of intersections: 0, 1 (touching or intersection of infinite region) 
  /// or 2 because region is convex.
  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const{
    return 0;
  }

  /// \en finds the point on the region surface closest to the given one p
  /// returns distance to this point (positive is P is outside and negative is P is inside region),
  /// surfp -- projection, surfn -- norm vector at surfp 
  virtual vec_type SurfProject(const vector_t &p, vector_t *surfp=NULL, vector_t *surfn=NULL) const{
    return 0;
  }

  ///\en Returns the area fraction of the contour part that is inside the region,
  /// subcontour which is inside the region (if pointer is not NULL)
  /// and the center of the subcontour (if pointer is not NULL)
  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<N> *subcont=NULL, vector_t *subcenter=NULL) const{
    return 1;
  }

  /// \en virtual version of TestContour for particulat PtrContour
  virtual vec_type TestPtrContour(const PtrContour<N> &cnt, VecContour<N> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

/*  virtual vector_t GetSize() const{
    vector_t v1, v2;
    GetBoundingBox(&v1,&v2);
    return v2-v1;
  }*/

  /// \en Returns opposite corners of the box which consist the region inside itself
  virtual vector_t GetBoundingBox(vector_t *v1, vector_t *v2) const{
    *v1=vector_t(-VEC_INFTY);
    *v2=vector_t(VEC_INFTY);
    return vector_t();
  }

  ///\en Returns volume of the region. '-1' means that function is not implemented.
  virtual vec_type Volume() const{
    return -1;
  }

  // @return a copy of itself optionally shifted in space with the shift vector
  virtual Region *Clone(const vector_t& shift=vector_t()) const{
    return new Region(*this);
  }

  ///\en return pointer on object RegDumper, which has function Dump drawing the region
  ///\ru Возвращает указатель на объект RegDumper, у которого есть отвечающая за рисование функция Dump
  virtual RegDumper<N>* CreateDumper() const{
    return GetRegionDumper(this);
  }

  virtual ~Region(){}


  ///\en Returns a fraction of the edge [0-1] that belongs to the region
  virtual vec_type GetInsideEdgePart(const vector_t &p1, const vector_t &p2) const{
    vec_type frac[2];
    int res=TestLine(p1,p2-p1,frac,NULL,NULL);
    if(res<=0)
      return 0.;

    if(res==2){
      frac[0]=max(vec_type(0),frac[0]);
      frac[1]=min(vec_type(1),frac[1]);
      return frac[1]-frac[0];
    }
    // res==1
    if(frac[0]<=0.){
      if(TestPoint(p1))
        return 1.;
      else 
        return 0.;
    }
    else{
      if(TestPoint(p1))
        return min(vec_type(1.),frac[0]);
      else 
        return max(vec_type(1.)-frac[0],vec_type(0.));
    }
  }
  ///\en Returns a fraction of the edge from the first point 
  /// till the intersection with the surface of the region.
  /// The intersection point may be calculated as pi=p1+frac*(p2-p1)
  /// frac ==-2 means that the function is not implemented
  /// frac ==-1 means no intersection
  /// fills pointers (if nonzero): surfp -- intersection point on the surface, surfn -- norm vector at surfp 
  virtual vec_type TestEdge(const vector_t &p1, const vector_t &p2, vector_t *surfp=NULL, vector_t *surfn=NULL) const{
    vec_type frac[2];
    vector_t surfp2[2],surfn2[2];
    int res=TestLine(p1,p2-p1,frac,surfp ? surfp2 : NULL,surfn ? surfn2 : NULL);
   
    for(int i=0;i<res;i++){
      if(frac[i]>1 && accomp(frac[i],vec_type(1)))
        frac[i]=1;
      if(frac[i]>=0 && frac[i]<=1){
        if(surfp)
          *surfp=surfp2[i];
        if(surfn)
          *surfn=surfn2[i];
        return frac[i];
      }
    }
    return -1.;
  }

  ///\en Returns the fraction of the ray corresponding to the nearest intersection in positive direction
  ///    frac<0 means no intersection in positive direction
  ///    intersections with frac<=epsilon are ignored.
  ///    Fills pointers (if nonzero): surfp -- intersection point on the surface, surfn -- norm vector at surfp .
  ///    \return fraction>=0. or  -1. if there is no intersection in positive direction, return values less than -1 \n
  ///            may be treated as implementation dependent error codes.            
  virtual vec_type TestRay(const vector_t &p1, const vector_t &dir, vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const{
    vec_type frac[2];
    vector_t surfp2[2],surfn2[2];
    int res=TestLine(p1,dir,frac,surfp ? surfp2 : NULL,surfn ? surfn2 : NULL);
    for(int i=0;i<res;i++){
      if(frac[i]>=0){
        if(surfp)
          *surfp=surfp2[i];
        if(surfn)
          *surfn=surfn2[i];
        return frac[i];
      }
    }
    return -1.;
  }
};

typedef Region<2> Region_2;
typedef Region<3> Region_3;

///\en Retained for compatibility with earlier versions.
typedef Region<3> SpaceRegion;

///\en Returns the area fraction of the contour part that is inside the region reg,
/// subcontour which is inside the region (if pointer is not NULL)
/// and the center of the subcontour (if pointer is not NULL) by Monte Carlo method.
/// nt is number of random points used in Monte Carlo method.
template<int N,class contour_t>
vec_type MonteCarloTestContour(const Region<N> *reg, const contour_t &cnt, Vector_Nt<vec_type,N> *subcenter=NULL, int nt=1000);

///\en Records the fractions frac1 and frac2 of the interval (p,p+dir*length/dir.norm()) corresponding to the intersections
/// surfpi=p+fraci*dir, surfni - coords of intersection and corresponding outer normal.
/// returns the number of intersections: 0, 1 or 2 (if there are more than 2 intersections, skip them).
/// in calculations nt points at the interval (p, p+length*dir/dir.norm()) are checked 
/// if they are inside or outside the region
template<int N>
int MonteCarloTestLine(const Region<N> *reg, const Vector_Nt<vec_type,N> &p, const Vector_Nt<vec_type,N> &dir, vec_type *frac, 
Vector_Nt<vec_type,N> *surfp, Vector_Nt<vec_type,N> *surfn, vec_type length, int nt=1000);

///\en Regular n-dimensional parallelepiped with edges parallel to coordinate axis
///\ru Прямоугольный n-мерный параллелепипед, у которого ребра направлены вдоль координатных осей
template<int N>
class Box_N: public Region<N>{
public:
  typedef typename Region<N>::vector_t vector_t;
protected:
  vector_t p1, p2, sz; // opposite vertices of the box and diagonal

public:

  Box_N(){}

  Box_N(const vector_t &sp1, const vector_t &sp2){
    init(sp1,sp2);
  }

  ///\en Returns a set of bit flags specifying which functions
  ///    are implemented for this region (see \ref RegFlags).
  ///    \return HAS_TESTPOINT|HAS_CLONE|HAS_VOLUME
  virtual int GetFlags() const {
    return HAS_TESTPOINT|HAS_CLONE|HAS_VOLUME;
  }

  virtual vec_type Volume() const{
    if(!valid())
      return 0.;
    vec_type vol=1;
    for(int i=0;i<N;i++){
      if(fabs(sz[i])<VEC_INFTY)
        vol*=sz[i];
      else 
        return VEC_INFTY;
    }
    return vol;
  }

  ///\en sp1, sp2 are opposite vertices of the box 
  void init(const vector_t &sp1, const vector_t &sp2);

  ///\en iterator at 6 vertices of the box
  class point_it {
    friend class Box_N;
    const Box_N *parent;
    int i;
    point_it(const Box_N *sparent, int si=0):parent(sparent), i(si){}
  public:
    typedef vector_t value_type;
    /// default constructor pointing to the end
    point_it():parent(NULL),i(0x1<<N){}

    point_it &operator++(){
      i++;
      return *this;
    }
    point_it operator++(int){
      point_it tmp=*this;
      ++*this;
      return tmp;
    }
    /// point construction
    vector_t operator*() const{
      vector_t p=parent->p1;
      for(size_t j=0;j<N;j++)
        if((i&(1<<j)))
          p[j]=parent->p2[j];
      return p;
    }
    bool operator!=(const point_it &other) const{
      return i!=other.i;
    }
  };

  ///\en returns begin iterator point_it at vertex of the box with minimal x, y and z coordinates
  point_it points_begin() const{
    return point_it(this,0);
  }

  ///\en returns end iterator point_it at vertex of the box
  point_it points_end() const{
    return point_it(this,0x1<<N);
  }

  vector_t get_p1() const{
    return p1;
  }

  vector_t get_p2() const{
    return p2;
  }

  vector_t GetSize() const{
    return sz;
  }

  ///\en this box becomes an intersection of itself with other box
  Box_N &operator&=(const Box_N &other){
    for(size_t i=0;i<N;i++){
      p1[i]=max(p1[i],other.p1[i]);
      p2[i]=min(p2[i],other.p2[i]);
      sz[i]=p2[i]-p1[i];
    }
    return *this;
  }

  ///\en @return true if the box defines nonvoid region
  bool valid() const{
    for(int i=0;i<N;i++){
      if(p1[i]>p2[i])return false;
    }
    return true;
  }

  ///\en @return true of the other box is inside this box 
  bool TestBox(const Box_N &other) const {
    for(int i=0;i<N;i++){
      if(other.p2[i]<p1[i] || p2[i]<other.p1[i])
        return false;
    } 
    return true;
  }

  bool TestPoint(const vector_t &p) const{
    for(int i=0;i<N;i++){
//      if(p[i]<p1[i] || p[i]>=p2[i])
      if(p[i]<p1[i] || p[i]>p2[i])
        return false;
    }
    return true;
  }

  virtual vector_t GetBoundingBox(vector_t *v1, vector_t *v2) const{
    *v1=p1;
    *v2=p2;
    return (p1+p2)/2.;
  }

  ///\en @return a copy of itself optionally shifted in space with the shift vector
  virtual Region<N> *Clone(const vector_t& shift=vector_t()) const{
    return new Box_N(p1+shift,p2+shift);
  }
};

///\en polygon
///\ru Многоугольник
class Polygon_2: public Region_2, public VecContour<2>{
public:
  bool TestPoint(const Vector_2 &p) const{
    return VecContour<2>::TestPoint(p);
  }
};

///\en n-dimensional sphere
template<int N>
class Sphere_N: public Region<N>{
public:
  typedef typename Region<N>::vector_t vector_t;
protected:
 // using typename Region<N>::vector_t; //I.V.: reformulated as typedef to compile with VS8 - this change is not compilable in UNIX
  
  vector_t center;
  vec_type R;

public:

  Sphere_N(vec_type R, const vector_t &center){
    init(R, center);
  }

  ///\en R is radius of the sphere, center is senter of the sphere
  void init(const vec_type &R_, const vector_t &center_) {
    R=R_;
    center=center_;
  }

  vector_t get_center() const{
    return center;
  }

  vec_type get_radius() const{
    return R;
  }

  bool TestPoint(const vector_t &p) const{
    return ((p-center).norm2()<=R*R);
//    return acless((p-center).norm2(),R*R,R*R);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  virtual typename Sphere_N<N>::vector_t GetBoundingBox(vector_t *v1, vector_t *v2) const;

  virtual Region<N> *Clone(const vector_t& shift=vector_t()) const{
    return new Sphere_N(R,center+shift);
  }
};


///\en circle, 2D sphere
///\ru 2D круг
class Circle: public Sphere_N<2>{

public:

  ///\en R is radius, center is center of the circle
  ///\ru R - радиус, center - координаты центра круга
  Circle(const vec_type &R, const Vector_2 &center): Sphere_N<2>(R,center){}

  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<2> *subcont=NULL, Vector_2 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<2> &cnt, VecContour<2> *subcont=NULL, Vector_2 *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual RegDumper<2>* CreateDumper() const{
    return GetRegionDumper(this);
  }
};

///\en Region obtained from regoin reg by shift and linear transformation specified by matrix basis
///\ru Регион, полученный из исходного региона reg путем линейного преобразования, задаваемого матрицей basis,
/// и сдвига на расстояние shift
template<class reg_tt>
class StretchedRegion: public Region<reg_tt::dimension>{

protected:

  typedef reg_tt reg_t;

public:

  typedef typename Region<reg_t::dimension>::vector_t vector_t;

  mngptr<reg_t> reg; // original region (before linear transformation and shift)

  vector_t shift; // distance between this region and region reg
  Basis_N<reg_t::dimension> basis; // linear transformation
  Basis_N<reg_t::dimension> basis_inv; // inverse transformation

  // reg_ is initial region, shift is shift and basis is linear transformation
  StretchedRegion(mngarg<reg_t> reg_,const vector_t &shift_,const Basis_N<reg_t::dimension> &basis_=Basis_N<reg_t::dimension>()): 
    reg(reg_),shift(shift_),basis(basis_),basis_inv(basis_.inv()){
  }

  bool TestPoint(const vector_t &p) const{
    return reg->TestPoint(basis_inv(p-shift));
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  template<class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<reg_tt::dimension> *subcont=NULL, typename reg_tt::vector_t *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<reg_tt::dimension> &cnt, VecContour<reg_tt::dimension> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual RegDumper<reg_t::dimension>* CreateDumper() const{
    return GetRegionDumper(this);
  }
};

///\en union of arbitrary number of regions
template<int N>
class RegionUnion: public Region<N>{
  vector<Region<N> *> regs;
  vec_type mc_length;
public:

  typedef typename Region<N>::vector_t vector_t;

  // mc_length_ is test length used in MonteCarloTestLine function
  RegionUnion(vec_type mc_length_ =1e-7):mc_length(mc_length_){}

  int GetN() const {
    return (int)regs.size();
  }

  ///\en adds new region to union
  int AddRegion(Region<N> *reg){
    regs.push_back(reg);
    return regs.size();
  }

  bool TestPoint(const vector_t &p) const{
    for(size_t i=0;i<regs.size();i++){
      if(regs[i]->TestPoint(p))
        return true;
    }
    return false;
  }

  int TestLine(const Vector_3 &p, const Vector_3 &dir, vec_type *frac,
  Vector_3 *surfp=NULL, Vector_3 *surfn=NULL, vec_type epsilon=0) const{
    return MonteCarloTestLine(this,p,dir,frac,surfp,surfn,mc_length);
  }

};

///\en intersection of arbitrary number of regions
template<int N>
class RegionIntersection: public Region<N>{
  vector<Region<N> *> regs;
  vec_type mc_length;
public:

  typedef typename Region<N>::vector_t vector_t;

  // mc_length_ is test length used in MonteCarloTestLine function
  RegionIntersection(vec_type mc_length_=1e-7):mc_length(mc_length_){}

  ///\en adds new region to intersection
  int AddRegion(Region<N> *reg){
    regs.push_back(reg);
    return regs.size();
  }

  bool TestPoint(const vector_t &p) const{
    for(size_t i=0;i<regs.size();i++){
      if(!regs[i]->TestPoint(p))
        return false;
    }
    return true;
  }

  int TestLine(const Vector_3 &p, const Vector_3 &dir, vec_type *frac,
  Vector_3 *surfp=NULL, Vector_3 *surfn=NULL, vec_type epsilon=0) const{
    return MonteCarloTestLine(this,p,dir,frac,surfp,surfn,mc_length);
  }

};

///\en @returns regular polygon,
/// center is center, d is distance from center to vertex, n is number of vertices
///\ru возвращает правильный многоугольник
/// center - центр, d - расстояние от центра до вершины, n - количество вершин
Polygon_2 *GetRegularPolygon(const Vector_2 &center, vec_type d, int n);



template<int N>
class Wire_N: public Region<N> {
public:
   ///\en Category of the region to determine the correct template behaviour for 
  ///    various algorithms (see \ref Body categories). 
  typedef unknown_cat category;

  typedef typename Region<N>::vector_t vector_t;
  static const int dimension=N;
public:
  vec_type radius;
  std::vector<vector_t> vertices;


  Wire_N(vec_type radius_=0.):radius(radius_){}

  size_t AddVertex(const vector_t &vertex){
    vertices.push_back(vertex);
    return vertices.size();
  }

  ///\en Returns a set of bit flags specifying which functions
  ///    are implemented for this region (see \ref RegFlags).
  ///    \return HAS_TESTPOINT|HAS_CLONE
  virtual int GetFlags() const {
    return HAS_TESTPOINT|HAS_CLONE|NOT_CONVEX;
  }

  ///\en Returns true if the point is in the the region.
  virtual bool TestPoint(const vector_t &p) const{
    size_t sz = vertices.size();
    if(!sz)
      return false;
    if(sz==1){ // sphere ?
      Sphere_N<N> sphere(radius,vertices[0]);
      return sphere.TestPoint(p);
    }
    vec_type dmin = VEC_INFTY;
    for(size_t i=0;i<sz-1;i++){ // finding segment with minimal distance
      vector_t vp = p-vertices[i];
      vector_t vv = vertices[i+1]-vertices[i];
      vec_type nvv = vv.normalize();
      vec_type proj = vv*vp;
      vec_type dist;
      if(proj<0) // distance to first point
        dist = vp.norm();
      else if(proj>nvv) // distance to second point
        dist = (p-vertices[i+1]).norm();
      else // distance to segment
        dist = (vp - proj*vv).norm();
      if(dist<dmin)
        dmin = dist;
    }
    if(dmin<=radius)
      return true;
    else   
      return false;
  }

 
  /// \en finds the point on the region surface closest to the given one p
  /// returns distance to this point (positive is P is outside and negative is P is inside region),
  /// surfp -- projection, surfn -- norm vector at surfp 
  virtual vec_type SurfProject(const vector_t &p, vector_t *surfp=NULL, vector_t *surfn=NULL) const{
    return 0;
  }

  ///\en Returns the area fraction of the contour part that is inside the region,
  /// subcontour which is inside the region (if pointer is not NULL)
  /// and the center of the subcontour (if pointer is not NULL)
  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<N> *subcont=NULL, vector_t *subcenter=NULL) const{
    return 0;
  }

 
  /// \en Returns opposite corners of the box which consist the region inside itself
  virtual vector_t GetBoundingBox(vector_t *v1, vector_t *v2) const{
    vector_t vmin(VEC_INFTY), vmax(-VEC_INFTY);
    for(size_t i=0;i<vertices.size();i++){
      for(size_t j=0;j<N;j++){
        vmin[j] = min(vmin[j], vertices[i][j]);
        vmax[j] = max(vmax[j], vertices[i][j]);
      }
    }
    *v1=vmin;
    *v2=vmax;
    return 0.5*(vmin+vmax);
  }

  ///\en Returns volume of the region. '-1' means that function is not implemented.
  virtual vec_type Volume() const{
    return -1;
  }

  // @return a copy of itself optionally shifted in space with the shift vector
  virtual Region<N> *Clone(const vector_t& shift=vector_t()) const{
    return new Wire_N(*this);
  }

  ///\en return pointer on object RegDumper, which has function Dump drawing the region
  ///\ru Возвращает указатель на объект RegDumper, у которого есть отвечающая за рисование функция Dump
  virtual RegDumper<N>* CreateDumper() const{
    return GetRegionDumper(this);
  }

  virtual ~Wire_N(){}

};

typedef Wire_N<2> Wire_2;
typedef Wire_N<3> Wire_3;


#endif
