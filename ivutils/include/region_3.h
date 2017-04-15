/*s****************************************************************************
 * $Log: region_3.h,v $
 * Revision 1.16  2013/11/23 22:38:44  lesha
 * float is supported
 *
 * Revision 1.15  2013/10/28 20:46:44  lesha
 * volume_integral_const_t is moved back
 *
 * Revision 1.14  2013/10/27 23:13:36  lesha
 * dx is introduced to MakeBoxDomain
 *
 * Revision 1.13  2013/10/26 20:25:19  lesha
 * Integrate is added to grid
 *
 * Revision 1.12  2013/10/26 16:57:40  lesha
 * bugs are fixed
 *
 * Revision 1.11  2013/10/26 14:13:06  valuev
 * modified MAkeBoxDomains
 *
 * Revision 1.10  2013/10/25 06:49:26  lesha
 * balanced domain decomposition interface
 *
 * Revision 1.9  2013/10/25 06:13:47  lesha
 * documentation
 *
 * Revision 1.8  2013/10/15 19:11:02  lesha
 * *** empty log message ***
 *
 * Revision 1.7  2013/10/14 21:46:07  lesha
 * transition is moved here
 *
 * Revision 1.6  2013/10/03 13:21:18  valuev
 * russian comments
 *
 * Revision 1.5  2013/10/01 01:11:28  lesha
 * GetAreaFraction is simplified
 *
 * Revision 1.4  2013/09/30 22:18:13  lesha
 * GetAreaFraction is moved here
 *
*******************************************************************************/
#ifndef _REGION_H
#define _REGION_H

/// \en @file region.h \brief A collection of objects for various 3D bodies.
/// \ru @file region.h \brief Набор объектов, описывающих различные 3D геометрические тела.


#include "region_2.h"
#include "utiltl.h"

// polyhedron traits
template<class poly_t>
struct poly_traits{
  typedef void plane_it; 
};

///\en orthogonal box, see documentation to Box_N
class Box: public Box_N<3>{
  
public:
  typedef polyhedron_cat category;

  Box(){}

  Box(const Vector_3 &sp1, const Vector_3 &sp2){
    init(sp1,sp2);
  }

  Box(vec_type x1, vec_type x2,vec_type y1, vec_type y2,vec_type z1, vec_type z2){
    init(Vector_3(x1,y1,z1),Vector_3(x2,y2,z2));
  }

  ///\en Returns a set of bit flags specifying which functions
  ///    are implemented for this region (see \ref RegFlags).
  ///    \return HAS_TESTPOINT|HAS_TESTLINE|HAS_TESTCONTOUR|HAS_CLONE|HAS_VOLUME
  virtual int GetFlags() const {
    return HAS_TESTPOINT|HAS_TESTLINE|HAS_TESTCONTOUR|HAS_CLONE|HAS_VOLUME;
  }

  class plane_it {
    friend class Box;
    const Box *parent;
    int i;
    plane_it(const Box *sparent, int si=0):parent(sparent), i(si){}
  public:
    typedef Plane_3 value_type;
    /// default constructor pointing to the end
    plane_it():parent(NULL),i(6){}
    
    plane_it(const plane_it &other):parent(other.parent),i(other.i){}
    plane_it &operator++(){
      i++;
      return *this;
    }
    plane_it operator++(int){
      plane_it tmp=*this;
      ++*this;
      return tmp;
    }
    /// plane construction
    Plane_3 operator*() const{
      Vector_3 n;
      n[i%3]= vec_type(i>2? -1: 1); // other coords =0
      return Plane_3(n,(i>2? parent->p2 : parent->p1));
    }
    bool operator!=(const plane_it &other) const{
      return i!=other.i;
    }
  };

  plane_it planes_begin() const{
    return plane_it(this,0);
  }
  plane_it planes_end() const{
    return plane_it(this,6);
  }

  virtual vec_type Volume() const{
    vec_type vol=1;
    for(int i=0;i<3;i++)vol*=sz[i];
    return vol;
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  ///\en Returns the area fraction of the contour part that is inside the region 
  /// and subcontour which is inside the region (if pointer is not NULL)
  /// and the center of the subcontour (if pointer is not NULL)
  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

// polyhedron traits
template<>
struct poly_traits<Box>{
  typedef Box::plane_it plane_it; 
};

///\en polyhedron defined as a set of planes
template <class plane_itt>
class Polyhedron: public Region_3{
protected:
  plane_itt b, e; // start and end iterator at plane sequence which confine the polyhedron
  Box box;
public:
  typedef polyhedron_cat category;
  typedef plane_itt plane_it;

  Polyhedron() {};

  Polyhedron(plane_it beg, plane_it end):box(Vector_3(-VEC_INFTY),Vector_3(VEC_INFTY)){
    init(beg,end);
  }

  plane_it planes_begin() const{
    return b;
  }

  plane_it planes_end() const{
    return e;
  }

  ///\en start and end iterator at plane sequence which confine the polyhedron
  void init(plane_it beg, plane_it end);

  ///\en returns true if the point is in the the region
  /// plane normals are pointing inside the region!!!!
  virtual bool TestPoint(const Vector_3 &p) const;

  ///\en gets the distance between a point and the closest plane
  /// (which is the distance to the polyhedron surface for internal points)
  vec_type MinPlaneDist(const Vector_3 &pos, plane_it *mit=NULL);
  
  ///\en returns the area fraction of the contour part that is inside the region 
  /// and subcontour which is inside the region (if pointer is not NULL)
  /// and the center of the subcontour (if pointer is not NULL)
  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  ///\en find polyhedron surface point surfp which is closest to the given point p and normal to polyhedron surfn at surfp
  /// return distance between p and surfp
  virtual vec_type SurfProject(const vector_t &p, vector_t *surfp=NULL, vector_t *surfn=NULL) const;

  virtual Vector_3 GetBoundingBox(Vector_3 *v1, Vector_3 *v2) const {
    return box.GetBoundingBox(v1, v2);
  }

  ///\en cloning is not possible for arbitrary plane iterator
  virtual Region_3 *Clone(const Vector_3& shift=Vector_3()) const {
    return NULL;
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

// polyhedron traits
template <class plane_itt>
struct poly_traits<Polyhedron<plane_itt> >{
  typedef typename Polyhedron<plane_itt>::plane_it plane_it; 
};


/*template<>
struct value_type_cast<Plane_3 *>{
  typedef Plane_3 value_type;
};*/

///\en polyhedron as array of planes
class Polyhedron_3: public Polyhedron<Plane_3 *>{
protected:
  mngptr<Plane_3> pptr;
public:
  typedef polyhedron_cat category;
  typedef Polyhedron<Plane_3 *> base_t;
  typedef base_t::plane_it plane_it;

  Polyhedron_3(){}

  ///\en copy constructor
  Polyhedron_3(const Polyhedron_3 &other){
    size_t sz=other.e-other.b;
    if(!sz)
      return;
    pptr.reset(new Plane_3[sz],0x8);
    for(size_t i=0;i<sz;i++)
      pptr.ptr()[i]=other.b[i];

    Polyhedron<Plane_3 *>::init(pptr.ptr(),pptr.ptr()+sz);
    box=other.box;
  }

  Polyhedron_3(Plane_3 *beg, Plane_3 *end, bool managed=true):base_t(beg,end),pptr(beg,managed ? 0x8 : 0){}

  ///\en the pointers must be valid for the whole lifetime of the class
  /// beg and end start and end pointers at array of planes
  /// managed: false if this array is not managed, true if this array will be delete if not used anymore
  void init(Plane_3 *beg, Plane_3 *end, bool managed=true){
    base_t::init(beg,end);
    pptr.reset(beg,managed ? 0x8 : 0);
  }
  
  ///\en this constructor copies the planes to the current class
  /// the managed flag is ignored (as if were always true)
  template<class plane_itt>
  void init(plane_itt beg, plane_itt end, bool managed){
    plane_itt it;
    size_t sz=0;
    for(it=beg;it!=end;++it)sz++;
    Plane_3 *data= sz? new Plane_3[sz] : NULL;
      
    pptr.reset(data, sz? 0x8: 0);
    sz=0;
    // copying
    for(it=beg;it!=end;++it){
      data[sz]=*it;
      sz++;
    }
    Polyhedron<Plane_3 *>::init(data,data+sz);
  }

  ///\en this constructor copies the planes to the current class
  /// the managed flag is ignored (as if were always true)
  template<class plane_itt>
  Polyhedron_3(plane_itt beg, plane_itt end, bool managed){
    init(beg,end,managed);
  }

  Polyhedron_3(const vector<Plane_3> &planes){
    init(planes.begin(),planes.end(),true);
  }

  Polyhedron_3(const Box &b){
    init(b.planes_begin(),b.planes_end(),true);
  }

  Polyhedron_3(const Plane_3 &pl){
    init(&pl,&pl+1,true);
  }

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  ///\en @return a copy of itself optionally shifted in space with the shift vector
  virtual Region_3 *Clone(const Vector_3& shift=Vector_3()) const {
    size_t sz=e-b;
    if(!sz)
      return NULL;
    Plane_3 *nplanes= new Plane_3[sz];
    for(size_t i=0;i<sz;i++){
      nplanes[i]=b[i];
      nplanes[i].shift(shift);
    }
    return new Polyhedron_3(nplanes,nplanes+sz);
  }
};

// polyhedron traits
template<>
struct poly_traits<Polyhedron_3>{
  typedef Polyhedron_3::plane_it plane_it; 
};

///\en 3D sphere
class Sphere: public Sphere_N<3>{

public:
  typedef sphere_cat category;

  Sphere(vec_type R, const Vector_3 &center): Sphere_N<3>(R,center){}

  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual Region<3> *Clone(const vector_t& shift=vector_t()) const{
    return new Sphere(R,center+shift);
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en
/// Cylinder with arbitrary 2D base of the type base_t
/// Cylinder axis passes the point origin and directed along vector n
/// Axes of 2D plane which contain the base are directed along vector x and y
/// Base is specified in 2D coordinate system defined by origin, x and y
/// Vector n is perpendicular to x and y (cylinder axis is perpendicular to the base)
/// "Tilted" cylinder with axis not perpendicular to the base
/// can be specified using template StretchedRegion
///\ru
/// Цилиндр с произвольным двумерным основанием типа base_t
/// Ось цилиндра прохожит через точку origin и направлена вдоль оси n
/// Оси x и y двумерной плоскости, на которой лежит основание, направлены вдоль векторов x и y
/// Подразумевается, что вектор n перпендикулярен x и y (основание перпендикулярно оси циллиндра)
/// "Косой" циллиндр, у которого заданное основание не перпендикулярно оси циллиндра, 
/// нужно задавать с помощью шаблона StretchedRegion
template<class base_tt>
class Cylinder: public Region_3{
public:
  typedef base_tt base_t;

  Vector_3 origin, n, x, y;
  ///\ru  основание
  mngptr<base_t> base; ///<\en base

  ///\en project point v at base 2D coordinate system and record to height distance to this base plane
  ///\ru проецирует точку на плоскость основания и записывает в height высоту над этой плоскостью
  Vector_2 Vector_3to2(const Vector_3 &v3, vec_type *height=NULL) const{
    Vector_3 tmp=v3-origin;
    if(height)*height=tmp*n;
    return Vector_2(tmp*x, tmp*y);
  }

  
  ///\en move point from 2D base coordinate system to 3D space and shift from the base plane at distance height
  ///\ru переводит точку обратно в трехмерное пространство и поднимает на высоту height
  Vector_3 Vector_2to3(const Vector_2 &v2, const vec_type height=0) const{
    return origin+v2[0]*x+v2[1]*y+height*n;
  }

  Cylinder(const Vector_3 &origin, const Vector_3 &n, const Vector_3 &x, const Vector_3 &y, mngarg<base_t> base){
    init(origin, n, x, y, base);
  }

//  Cylinder(Cylinder &other): origin(other.origin), n(other.n), x(other.x), y(other.y) {
//    base.reset(make_mngarg(new base_t(*(other.base))));
//  }

  // n, x, y must be orthonormal
  void init(const Vector_3 &origin_, const Vector_3 &n_, const Vector_3 &x_, const Vector_3 &y_, mngarg<base_t> base_){
    origin=origin_, n=n_, x=x_, y=y_, base.reset(base_);
  }

  base_t *get_base() const{
    return base.ptr();
  }

  bool TestPoint(const Vector_3 &p) const{
    return base->TestPoint(Vector_3to2(p));
  }

  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  virtual Region_3 *Clone(const Vector_3& shift=Vector_3()) const {
    return new Cylinder(origin+shift,n,x,y, make_mngarg(new base_t(*base)));
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en Cone with arbitrary 2D base of the type base_tt
/// Cone axis passes the point origin and directed along vector n
/// Axes of 2D plane which contain the base are directed along vector x and y
/// Base is specified in 2D coordinate system defined by origin, x and y
/// Vector n is perpendicular to x and y (cone axis is perpendicular to the base)
/// Top of the cone is origin+L*n
/// "Tilted" cone with axis not perpendicular to the base
/// can be specified using template StretchedRegion
template<class base_tt>
class Cone: public Cylinder<base_tt>{
public:  
  typedef base_tt base_t;
  typedef typename Region<Cylinder<base_tt>::dimension>::vector_t vector_t;
public:
  vec_type L; ///\en distance between origin and top of the cone
  using Cylinder<base_tt>::origin;
  using Cylinder<base_tt>::n;
  using Cylinder<base_tt>::x;
  using Cylinder<base_tt>::y;
  using Cylinder<base_tt>::base;
public:

  Cone(const Vector_3 &origin, const Vector_3 &n, const Vector_3 &x, const Vector_3 &y, mngarg<base_t> base, vec_type L_):
    Cylinder<base_t>(origin,n,x,y,base),L(L_){}

  bool TestPoint(const Vector_3 &p) const{
    vec_type h;
    Vector_2 p2=Cylinder<base_t>::Vector_3to2(p,&h);
    if(h>=L)return false; // behind top
    Vector_2 o2=Cylinder<base_t>::Vector_3to2(origin);
    Vector_2 dist=p2-o2;
    dist*=L/(L-h);
    return base->TestPoint(o2+dist);
  }

  // not realized yet
  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  virtual Region_3 *Clone(const Vector_3& shift=Vector_3()) const {
    return new Cone(origin+shift,n,x,y, make_mngarg(new base_t(*base)),L);
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en Intersection of polyhedron with arbitrary other region.
/// Can be used to specify hemisphere, finite cylinders etc.
///\ru тело, получающееся в результате пересечения многогранника с каким-нибудь другим произвольным телом
/// может использоваться для получения полусфер, ограниченных циллиндров и т.д.
template<class reg_tt>
class ConfinedRegion: public Region_3{
public:
  typedef reg_tt reg_t;

  mngptr<reg_t> reg; ///\en arbitrary region
  mngptr<Polyhedron_3> poly; ///\en polyhedron which confines reg

  ConfinedRegion(mngarg<reg_t> reg_, mngarg<Polyhedron_3> poly_): reg(reg_), poly(poly_){}

  bool TestPoint(const Vector_3 &p) const{
    return reg->TestPoint(p) && poly->TestPoint(p);
  }

  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;


  virtual Region_3 *Clone(const Vector_3& shift=Vector_3()) const {
    return new ConfinedRegion(make_mngarg((reg_t *)reg->Clone(shift)),make_mngarg((Polyhedron_3 *)poly->Clone(shift)));
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en Complement to the region in 3D space
///\ru Дополнение к телу в трехмерном пространстве
template<class reg_tt>
class Inverse: public Region_3{
public:
  typedef reg_tt reg_t;
  mngptr<reg_t> reg;

  Inverse(mngarg<reg_t> reg_=NULL): reg(reg_){}

  bool TestPoint(const Vector_3 &p) const{
    return !(reg->TestPoint(p));
  }

  virtual int TestLine(const Vector_3 &p, const Vector_3 &dir, vec_type *frac,
  Vector_3 *surfp=NULL, Vector_3 *surfn=NULL, vec_type epsilon=0) const{
    int res=reg->TestLine(p,dir,frac,surfp,surfn,epsilon);
    if(surfn){
      for(int i=0;i<res;i++)
        surfn[i]*=-1;
    }
    return res;
  }

  virtual Region_3 *Clone(const Vector_3& shift=Vector_3()) const{
    if(!reg.ptr())return NULL;
    return new Inverse(make_mngarg((reg_t *)reg->Clone(shift)));
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en returns "polyhedron" formed by one plane (half-space)
///\ru возвращает "многогранник", образованный одной плоскостью (полупространство)
Polyhedron_3 *GetHalfSpace(const Vector_3 &n, const Vector_3 &pos);
///\en returns "polyhedron", formed by 2 parallel planes (infinite plane)
///\ru возвращает "многогранник", образованный двумя параллельным плоскостями (бесконечная пластинка)
Polyhedron_3 *GetPlate(const Vector_3 &n, const Vector_3 &pos, vec_type width);
///\en returns polyhedron - regular parallelipiped
///\ru возвращает многогранник - прямоугольный параллелепипед
Polyhedron_3 *GetBox(const Vector_3 &p1, const Vector_3 &p2);
///\en returns polyhedron - prism
///\ru возвращает многогранник - призму
Polyhedron_3 *GetPrism(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base);
///\en returns polyhedron - pyramid
/// O is center of the polygon in the base, top is top of the pyramid, 
/// x and y are axis directions, polygon will be created relatively to these axis,
/// base is polygon , confinment is polyhedron which confines pyramid (if it is absent then pyramid is infinite in one direction)
///\ru возвращает многогранник - пирамиду
/// O - центр многоугольника, top - острие пирамиды, 
/// x и y - направления осей, относительно которых строится многоугольник,
/// base - многоугольник, confinment - многогранник, ограничивающий пирамиду (если его нет - пирамида бесконечна)
Polyhedron_3 *GetPyramida(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base);
///\en returns polyhedron which is intersection of two other polyhedra
///\ru возвращает многогранник, являющийся пересечением двух аргументов
Polyhedron_3 *GetConfinedPolyhedron(mngarg<Polyhedron_3> c1, mngarg<Polyhedron_3> c2);
///\en returns region which is intersection of some other region and polyhedron
template<class reg_t>
ConfinedRegion<reg_t> *GetConfinedRegion(mngarg<reg_t> c1,mngarg<Polyhedron_3> c2){
  return new ConfinedRegion<reg_t>(c1,c2);
}
///\en returns sphere of the radius R centered in the point center
Sphere *GetSphere(vec_type R, const Vector_3 &center);
///\en returns polyhedron for N sides which approximating sphere
///\ru возвращает многогранник из N граней, интерполирующий сферу
Plane_3* CreateSpherePlanes(const int &N, const vec_type &R, const Vector_3 &center);
///\en returns cylinder, origin is some point at cylinder axis, n is axis direction, R is radius
///\ru origin - произвольная точка на оси, n - направление оси, R - радиус циллиндра
Cylinder<Circle> *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R);
///\en height is height of finite cylinder, height is measured from origin in direction n
///\ru height - высота ограниченного цилиндра, отмеряемая от точки origin в направлении n
ConfinedRegion<Cylinder<Circle> > *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type height);
///\en height is height of cone, height is measured from origin in direction n
ConfinedRegion<Cone<Circle> > *GetCone(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type h, Vector_3 npl=0);
///\en R_top is upper radius (if R_top=0, then function works in a same way as GetCone)
ConfinedRegion<Cone<Circle> > *GetTruncatedCone(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type R_top, vec_type h);

///\en Volume integral of constant space density.
class volume_integral_const_t {
  vec_type value;
public:
  volume_integral_const_t(vec_type density = 1.):value(density){}
  ///\en returns an integral of some function over a volume enclosed between p1 and p2
  vec_type operator()(const Vector_3 &p1, const Vector_3 &p2) const{
    Box b(p1,p2);
    return b.Volume()*value;
  }
};

///\en Domain decomposition of the Box, abs(nproc) is number of domains.
/// nx, ny, nz - preferable number of sections per each dimension.
/// This function is used in parallel implementation of numerical algorithms requiring spatial domain decomposition: 
/// each processor simulates a subBox (domain) of the whole Box.
/// If any of the numbers is negative, then optimal number for the corresponding direction will be found.
/// Resulting subdomains will be collected in result, if nproc>0.
/// \return the number of domains that are generated (or will be generated if nproc<0)
/// \a weight_integral is a class which defines operator()(const Vector_3 &p1, const Vector_3 &p2) as a positive
/// function giving a weight proportional to CPU load inside a box limited by (p1,p2). 
/// The argument of type \ref volume_integral_const_t may be used when the CPU load is uniform inside \a box.
/// dx is space step used in domain decomposition algorithm (default value 0 means that dx will be chosen automatically)
template<class inp_it, class weight_integral_t>
int MakeBoxDomains(const Box &box, int nproc, int nx, int ny, int nz, inp_it result, const weight_integral_t &weight_integral, const Vector_3 &dx=0);

/*
class volume_increment{
public:
  /// rescale increments dl according to some weight function
  int operator()(const Vector_3 &p1, const Vector_3 &p2, int dir, vector<vec_type> &dl) const{
    return 1;
  }
};*/

template<class inp_it>
int MakeBoxDomains(const Box &box, int nproc, int nx, int ny, int nz, inp_it result){
  return MakeBoxDomains(box, nproc, nx, ny, nz, result, volume_integral_const_t(1.));
}

///\en gets the point having minimal projection on the direction k
/// infinite vectors are ignored
/// returns an infinite vector for void sequence
template<class point_it> 
Vector_3 get_min_point(const Vector_3 &k, point_it beg, point_it end){
  vec_type dmin=VEC_INFTY;
  Vector_3 res(VEC_INFTY);
  for(;beg!=end;++beg){
    Vector_3 v=*beg;
    if(!v.infinite()){
      vec_type d=k*v;
      if(d<dmin){
        dmin=d;
        res=v;
      }
    }
  }
  return res;
}

///\en in are some orthogonal vectors, innum is their number
/// out are 3-innum vectors of the unit length which are orthogonal to all vectors in and to each other 
int build_orth_basis(Vector_3 *in, int innum, Vector_3 *out);

///\en calculates 2 unit vectors x, y perpeindicular to unit vector z
int build_orth_basis(const Vector_3 &z, Vector_3 &x, Vector_3 &y);

///\en get arbitrary unit vector in 3D
Vector_3 random_direction();

///\en returns random unit vector uniformely distributed in 3D space (uses standard rand())
Vector_3 randdir();

///\en returns random unit vector uniformely distributed in 2D space (uses standard rand())
Vector_2 randdir2();

///\en calculates the intersection contour of the plane with the box, returns true if the intersection is nonvoid
/// *start is updated with the contour point having minimal projection on k direction 
bool get_first_corner(const Box &box,const Plane_3 &plane,const Vector_3 &k,Vector_3 *start);

///\en *start is updated with the box point having minimal projection on k direction 
/// @return true if box is nonvoid, false otherwise
bool get_first_corner(const Box &box,const Vector_3 &k,Vector_3 *start);

///\en Find the intersection of a convex region with the line l(t)=origin+dir*t.
/// @return true and the segment coordinates in the line system into t1 and t2 (t1<=t2) or
///         false if no intersection found
bool find_cross_segment(const Region_3 &reg, const Vector_3& orig,const Vector_3& dir,vec_type &t1, vec_type &t2);

///\en Fills the space inside convex region by the lattice given by 3 elementary translations (cell) starting from origin
/// (orig). Puts the filled positions into points by push_back operation (no clearing is performed).
/// Optionally (if ipoints is not NULL) fills integer lattice indices indicating the translation numbers for each point.
/// @return the number of points filled (on success)
///         -1 no intersection of the lattice with the region found
///         -2 the region is unbounded
int fill_lattice(const Region_3 &reg, const Vector_3& orig,const Basis_3 &cell,vector<Vector_3> &points, vector<iVector_3> *ipoints=NULL);

///\en until Basis_Nt and calculating inverse matrix are not implemented, 
/// works only for 2D rectangular case
Vector_3 * make_reciprocal_vectors(const vec_type a1, const vec_type a2, const vec_type max, int &num);

/**\en gets area fraction of the contour (if volume_step=0) or volume fraction of the control volume around contour (if volume_step>0),
  which is entering the region reg
  subcontour inside the region is recorded to subcont (if subcont!=NULL)
  normal to the region (in proximity to contour intersection with region) is recorded to normv (if norm!=NULL)
  if volume_step>0, calculate fraction of the control volume, considering area fractions of 2*volume_step+1 
  copies of the contour in this volume
  This function can be used, for example, to calculate averaged dielectric function in FDTD subpixel smoothing,
  see A. Deinega and I. Valuev, “Subpixel smoothing for conductive and dispersive media in the FDTD method”, Optics Letters 32, 3429 (2007)
  
  \param a0 contour area
  \param center contour center
  \param dir normal to the contour
*/
/** \ru возвращает площадь, отсекаемую телом от контура, согласно алгоритму contour path.
  Заполняет поля subcont и normv.
  При расчете normv делает усреднение нормали к телу по контрольному объему.
  \param reg тело
  \param cnt контур
  \param a0 площадь контура
  \param center центр контура
  \param dir направление, перпендикулярное контуру
  \param subcont подконтур, образуемый пересечением контура с регионом
  \param normv вектор нормали к телу в районе пересечения им контура
  */
template<class contour_t>
vec_type GetContourFraction(const Region_3 *reg, const contour_t &cnt, int volume_step,
vec_type a0, const Vector_3 &center, const Vector_3 &dir, VecContour<> *subcont, Vector_3 *normv);

template<class result_t>
class region3_function: public virt_unary_function<Vector_3,result_t>{
  Region_3 *reg;
  result_t y;
public:
  region3_function(Region_3 *reg_, result_t y_):reg(reg_),y(y_){}
  virtual result_t operator()(Vector_3 x){
    return reg->TestPoint(x) ? y : 0;
  }
};

template<class result_t>
class sphere_function: public virt_unary_function<Vector_3,result_t>{
  Vector_3 c;
  vec_type R;
  vec_type dec;
  result_t y;

public:
  sphere_function(const Vector_3 &c_, vec_type R_, vec_type dec_, result_t y_):c(c_),R(R_),dec(dec_),y(y_){}
  virtual result_t operator()(Vector_3 x){
    vec_type dist = R - (x-c).norm();
    return dist<0 ? 0 :exp(-dist/dec)*y;
  }
};

/// function that continuosly changes from 0 to 1 near surface of region reg
/// (used to specify alloys in Microvolt)
template<class reg_t>
class transition: public virt_unary_function<const Vector_3 &, vec_type>{
  mngptr<reg_t> reg;
  vec_type dist;
  int inv;

public:
  transition(mngarg<reg_t> reg_, vec_type dist_=0, int inv_=0): reg(reg_), dist(dist_),inv(inv_){}
  vec_type operator()(const Vector_3 &pos) {
    vec_type d=reg->SurfProject(pos); // distance to the surface (negative or positive according to the surface normal)
    if(d)
      d = dist ? d/dist : d/fabs(d);
    vec_type fr=fabs(d)>=1 ? d>=1 : (d+1)/2;
    if(inv)
      fr=1-fr;
    return fr;
  }
};

#endif
