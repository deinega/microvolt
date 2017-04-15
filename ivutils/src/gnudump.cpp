#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include "gnudump.hpp"

const int NUM_RADII=20; ///\en number of sides of regular polygon used to interpolate circle

int DumpDirections(FILE *f, int n, const Vector_3 *origin, const Vector_3 *k, const Box *B, int first){
  if(n<=0)
    return 0;

  vec_type mlt=0;
  if(B){
    Vector_3 sz=B->GetSize();
    for(int i=0; i<3; i++)
      mlt+=sz[i];
    mlt*=.5/3.;
  }
  else 
    mlt=1;
    
  if(first)
    fprintf(f,"#1-ox 2-oy 3-oz 4-kx 5-ky 6-kz\n");
  for(int i=0;i<n;i++){
    Vector_3 k2=k[i]*mlt;
    fprintf(f,"%g %g %g %g %g %g\n", origin[i][0], origin[i][1], origin[i][2], k2[0], k2[1], k2[2]);
  }
  return 1;
}

interp_t<2> InterpolateRegion(const Polygon_2 *poly){
  return *((VecContour<2> *)poly);
}

///\en returns contour which approximates given circle
///\ru Создать из 2D круга описывающий его набор точек (2D координат)
interp_t<2> InterpolateRegion(const Circle *cone){
  VecContour<2> cnt;
  Vector_2 center = cone->get_center();
  vec_type R = cone->get_radius();
  for (int i=0; i<NUM_RADII; i++){
    vec_type fi = 2*M_PI*vec_type(i)/vec_type(NUM_RADII);
    cnt.add( center+Vector_2(R*cos(fi),R*sin(fi)) );
  }
  return cnt;
}

interp_t<3> InterpolateRegion(const Box *B){
  return GetBox(B->get_p1(),B->get_p2());
}

interp_t<3> InterpolateRegion(const Polyhedron<Plane_3 *> *B){
  Polyhedron_3 *poly = new Polyhedron_3;
  poly->init<Plane_3 *>(B->planes_begin(),B->planes_end(),1);
  return poly;
}

interp_t<3> InterpolateRegion(const Sphere *B){
  vec_type R=B->get_radius();
  Vector_3 center=B->get_center();

  const int NUM_MERIDIAN=10;
  const int NUM_PARALLEL=10;

  Plane_3 *planes=new Plane_3[NUM_MERIDIAN*(NUM_PARALLEL+1)+1];
  int sz=0;
  
  // meridians
  for(int i=0;i<=NUM_PARALLEL;i++){
    vec_type theta=M_PI*vec_type(i)/vec_type(NUM_PARALLEL);
    vec_type sin_theta=R*sin(theta);
    vec_type cos_theta=R*cos(theta);
    for(int j=0;j<NUM_MERIDIAN;j++){
      vec_type fi=2*M_PI*vec_type(j)/vec_type(NUM_MERIDIAN);
      Vector_3 n=Vector_3(sin_theta*cos(fi),sin_theta*sin(fi),cos_theta);
      Vector_3 pos=center+n;
      planes[sz++].init(-n, pos);
      if(i==0 || i==NUM_PARALLEL)
        break;
    }
  }
  return new Polyhedron_3(planes, planes+sz);
   // not working - error in ProjectSimplex
//  Plane_3 *pl=CreateSpherePlanes(100, B->get_radius(), B->get_center());
//  Polyhedron_3 *pol = new Polyhedron_3(pl, pl+100, 1);
//  delete []pl;
//  return pol;
}


// Explicit instantiation of functions

// Explicit template specification is needed for MSVC
INSTANTIATE_GetRegionDumper(Region<2>);
INSTANTIATE_GetRegionDumper(Circle);
INSTANTIATE_GetRegionDumper(Region<3>);
INSTANTIATE_GetRegionDumper(Box);
INSTANTIATE_GetRegionDumper(Polyhedron<Box::plane_it>);
INSTANTIATE_GetRegionDumper(Polyhedron<Plane_3*>);
INSTANTIATE_GetRegionDumper(Sphere);
INSTANTIATE_GetRegionDumper(Cylinder<Circle>);
INSTANTIATE_GetRegionDumper(Cone<Circle>);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Region_3>);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Polyhedron_3>);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Cylinder<Circle> >);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Cone<Circle> >);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Sphere>);
INSTANTIATE_GetRegionDumper(StretchedRegion<Region<3> >);
INSTANTIATE_GetRegionDumper(Inverse<Region_3>);
INSTANTIATE_GetRegionDumper(Inverse<Polyhedron_3>);

//template interp_t<3> InterpolateRegion<Inverse<Region_3> >(const Inverse<Region_3> *B);
template interp_t<3> InterpolateRegion<Region_3>(const Inverse<Region_3> *B);
template interp_t<3> InterpolateRegion<Inverse<Polyhedron_3> >(const Inverse<Polyhedron_3> *B);

template int DumpContours(FILE *f, const vector<VecContour<2> > &cnt, int first);
template int DumpContours(FILE *f, const vector<VecContour<3> > &cnt, int first);  
template int DumpContoursVTK(FILE *f, const vector<VecContour<2> > &cnt, int first);
template int DumpContoursVTK(FILE *f, const vector<VecContour<3> > &cnt, int first);  

//template int GnuDumper::DumpContour(const VecContour<2> &cnt,const char *fname);

template int DumpVectorSet(FILE *f, const vector<Vector_3> *vs, const VecTransform *transform, int first);
template int DumpVectorSet(FILE *f, const SpaceVectorSet *vs, const VecTransform *transform, int first);
template int DumpVectorSet(FILE *f, const UniformGrid<Vector_3> *vs, const VecTransform *transform, int first);
template int DumpVectorSet(FILE *f, const BoxSurfaceSet *vs, const VecTransform *transform, int first);
template int DumpVectorSet(FILE *f, const CylinderSurfaceSet *vs, const VecTransform *transform, int first);
template int DumpVectorSet(FILE *f, const SphereSurfaceSet *vs, const VecTransform *transform, int first);
