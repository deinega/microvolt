/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *   $Revision: 1.9 $
 *   $Date: 2013/11/01 16:39:59 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/region_3.cpp,v 1.9 2013/11/01 16:39:59 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/region_3.cpp,v $
$Revision: 1.9 $
$Author: lesha $
$Date: 2013/11/01 16:39:59 $
*/
/*s****************************************************************************
 * $Log: region_3.cpp,v $
 * Revision 1.9  2013/11/01 16:39:59  lesha
 * GetCylinder is modified to handle negative heights
 *
 * Revision 1.8  2013/10/28 20:46:44  lesha
 * volume_integral_const_t is moved back
 *
 * Revision 1.7  2013/10/27 23:13:36  lesha
 * dx is introduced to MakeBoxDomain
 *
 * Revision 1.6  2013/10/26 23:06:07  lesha
 * Integrate is added to grid
 *
 * Revision 1.5  2013/10/26 14:13:06  valuev
 * modified MAkeBoxDomains
 *
 * Revision 1.4  2013/10/25 06:49:26  lesha
 * balanced domain decomposition interface
 *
 * Revision 1.3  2013/02/22 11:01:26  valuev
 * thin surface interface
 *
 * Revision 1.2  2012/12/06 23:23:26  lesha
 * randdir is moved to region
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.65  2012/12/05 21:28:37  lesha
 * GetSphere is added
 *
 * Revision 1.64  2012/11/20 01:36:22  lesha
 * GetPolyhedron functions are renamed
 *
 * Revision 1.63  2012/11/10 23:59:18  lesha
 * GetTruncatedCone is added
 *
 * Revision 1.62  2012/10/02 23:11:30  lesha
 * SpaceRegion is renamed to Region_3
 *
 * Revision 1.61  2012/09/27 21:29:18  lesha
 * *** empty log message ***
 *
 * Revision 1.60  2012/09/24 21:50:13  lesha
 * documentation
 *
 * Revision 1.59  2012/08/23 05:25:09  lesha
 * GetPolyhedronPrysma amd Pyramid are fixed
 *
*******************************************************************************/
/// \file \brief Non-template function definitions for  region.h

#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include <algorithm>
#include "region_2.hpp"
#include "region_3.hpp"


int Box::TestLine(const vector_t &p, const vector_t &dir, 
vec_type *frac, vector_t *surfp, vector_t *surfn, vec_type epsilon) const{
  Polyhedron<Box::plane_it> poly(planes_begin(),planes_end());
  return poly.TestLine(p,dir,frac,surfp,surfn,epsilon);
}

Polyhedron_3 *GetBox(const Vector_3 &p1, const Vector_3 &p2){
  Box B(p1,p2);
  Plane_3 *pl = new Plane_3[6];
  int i=0;
  for (Box::plane_it it=B.planes_begin(),e=B.planes_end();it!=e;++it){
    pl[i++]=*it;
  }
  return new Polyhedron_3(pl,pl+6);
}

Polyhedron_3 *GetHalfSpace(const Vector_3 &n, const Vector_3 &pos){
  if(n==0)return NULL;
  Plane_3 *pl=new Plane_3[2];
  pl[0].init(n, pos);
  return new Polyhedron_3(pl, pl+1);
}

Polyhedron_3 *GetPlate(const Vector_3 &n, const Vector_3 &pos, vec_type width){
  if(n==0)return NULL;
  Plane_3 *pl=new Plane_3[3];
  pl[0].init(n, pos);
  pl[1].init(n, pos);
  Vector_3 n2;
  vec_type d;
  pl[1].get_coeff(n2, d);
  n2=-n2;
  d=-d+width;
  pl[1].set_coeff(n2, d);
  return new Polyhedron_3(pl, pl+2);
}

Polyhedron_3 *GetPrism(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base){
  int sz=base->GetNPoints();

  Plane_3 *planes = new Plane_3[sz];
  sz=0;

  int sign=0;
  for (generic_edge_it<VecContour<2>::point_iterator> it(base->points_begin(), base->points_end()), e(base->points_end(), base->points_end(), 1); it!=e; ++it) {
    edge_t<VecContour<2>::point_iterator> ed=*it;
    Vector_2 a=ed.get_p1(), b=ed.get_p2();
    Vector_3 a3=O+a[0]*x+a[1]*y, b3=O+b[0]*x+b[1]*y;

    Vector_3 p1=top-O, p2=b3-a3;
    p1.normalize(), p2.normalize();
    Vector_3 n=p1%p2;
    n.normalize();
    if (sign==0){
      Vector_2 ct=base->GetCenter();
      Vector_3 ct3=O+ct[0]*x+ct[1]*y;
      Vector_3 c=ct3-a3;
      sign=(n*c>0)? 1 : -1;
    }

    planes[sz++].init(sign*n, a3);
  }
  
  return new Polyhedron_3(planes, planes+sz);
}

Polyhedron_3 *GetPyramida(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base){

  int sz=base->GetNPoints();

  Plane_3 *planes = new Plane_3[sz];
  sz=0;

  int sign=0;
  for (generic_edge_it<VecContour<2>::point_iterator> it(base->points_begin(), base->points_end()), e(base->points_end(), base->points_end(), 1); it!=e; ++it) {
    edge_t<VecContour<2>::point_iterator> ed=*it;
    Vector_2 a=ed.get_p1(), b=ed.get_p2();
    Vector_3 a3=O+a[0]*x+a[1]*y, b3=O+b[0]*x+b[1]*y;

    Vector_3 p1=top-a3, p2=b3-a3;
    p1.normalize(), p2.normalize();
    Vector_3 n=p1%p2;
    n.normalize();
    if (sign==0){
      Vector_2 ct=base->GetCenter();
      Vector_3 ct3=O+ct[0]*x+ct[1]*y;
      Vector_3 c=ct3-a3;
      sign=(n*c>0)? 1 : -1;
    }

    planes[sz++].init(sign*n, a3);
  }
  
  return new Polyhedron_3(planes, planes+sz);
}

Polyhedron_3 *GetConfinedPolyhedron(mngarg<Polyhedron_3> c1, mngarg<Polyhedron_3> c2){
  Polyhedron_3 *poly[2]={c1.first,c2.first};
  int sz=0;
  for(int i=0;i<2;i++){
    if(poly[i]){
      Polyhedron_3::plane_it it=poly[i]->planes_begin(), e=poly[i]->planes_end();
      for(; it!=e; ++it)
        sz++;
    }
  }

  Plane_3 *planes = new Plane_3[sz];
  sz=0;

  for(int i=0;i<2;i++){
    if(poly[i]){
      Polyhedron_3::plane_it it=poly[i]->planes_begin(), e=poly[i]->planes_end();
      for (; it!=e; ++it)
        planes[sz++]=*it;
    }
  }

  if(c1.second)
    delete c1.first;
  if(c2.second)
    delete c2.first;
  
  return new Polyhedron_3(planes, planes+sz);
}

Cylinder<Circle> *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R){
  n.normalize();
  Vector_3 x(n[2]-n[1],n[0]-n[2],n[1]-n[0]);
  if(x==0)
    x=Vector_3(1,-1,0);
  x.normalize();
  Vector_3 y=n%x;
  return new Cylinder<Circle>(origin,n,x,y,make_mngarg(new Circle(R,Vector_2())));
}

ConfinedRegion<Cylinder<Circle> > *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type height){
    if(height<0){
      height=-height;
      n=-n;
    }
  return new ConfinedRegion<Cylinder<Circle> >(
    make_mngarg(GetCylinder(origin,n,R)),
    make_mngarg(GetPlate(n,origin,height))
  );
}

ConfinedRegion<Cone<Circle> > *GetCone(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type h, Vector_3 npl){
  n.normalize();
  Vector_3 x(n[2]-n[1],n[0]-n[2],n[1]-n[0]);
  if(x==0)
    x=Vector_3(1,-1,0);
  x.normalize();
  Vector_3 y=n%x;

  return new ConfinedRegion<Cone<Circle> >(
    make_mngarg(new Cone<Circle>(origin,n,x,y,make_mngarg(new Circle(R,Vector_2())),h)),
    make_mngarg(GetHalfSpace(npl==0 ? n : npl,origin)));
}

ConfinedRegion<Cone<Circle> > *GetTruncatedCone(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type R_top, vec_type h){
  n.normalize();
  Vector_3 x(n[2]-n[1],n[0]-n[2],n[1]-n[0]);
  if(x==0)
    x=Vector_3(1,-1,0);
  x.normalize();
  Vector_3 y=n%x;

  return new ConfinedRegion<Cone<Circle> >(
    make_mngarg(new Cone<Circle>(origin,n,x,y,make_mngarg(new Circle(R,Vector_2())),h*R/(R-R_top))),
    make_mngarg(GetPlate(n,origin,h)));
}

int build_orth_basis(Vector_3 *in, int innum, Vector_3 *out){
  for(int i=0;i<innum;i++){
    if(in[i].norm()<=VEC_ZERO)
      return 0;
    in[i].normalize();
  }

  if(innum==0){
    for(int i=0;i<3;i++){
      out[i]=Vector_3();
      out[i][i]=1;
    }
  }
  if(innum==1){
/*    for(int i=0;i<3;i++){
      int i1=(i+1)%3;
      int i2=(i+2)%3;
      out[0][i] = in[0][i1]-in[0][i2];
    }
    out[0].normalize();
    out[1]=in[0]%out[0];*/
    build_orth_basis(in[0],out[0],out[1]);
  }
  if(innum==2){
    out[0]=in[0]%in[1];
    if(out[0].norm()<=VEC_ZERO)
      return 0;
  }
  return 1;
}

int build_orth_basis(const Vector_3 &z, Vector_3 &x, Vector_3 &y){
  int i=0;
  for(;i<3;i++){
    if(z[i])break;
  }
  if(i==3)return -1;
  int j=(i+1)%3;
  x[i]=z[j];
  x[j]=-z[i];
  x[(j+1)%3]=0;
  x.normalize();
  y=z%x;
  return 0;
}

Vector_3 random_direction(){
  vec_type fi=2*M_PI*(vec_type)rand()/(vec_type)RAND_MAX;
  vec_type tetta=acos(-1+2*(vec_type)rand()/(vec_type)RAND_MAX);
  vec_type sin_tetta=sin(tetta);
  return Vector_3(sin_tetta*cos(fi),sin_tetta*sin(fi),cos(tetta));
}

Vector_3 randdir(){
  vec_type xi1 = 2.*((vec_type)rand()) / RAND_MAX-1.;
  vec_type xi2 =   ((vec_type)rand()) / RAND_MAX;
  vec_type r1  = sqrt (1.-xi1*xi1);
  return  Vector_3(r1*cos(2.*M_PI*xi2),r1*sin(2.*M_PI*xi2),xi1);
}

Vector_2 randdir2(){
  vec_type xi2 =   ((vec_type)rand()) / RAND_MAX;
  return  Vector_2(cos(2.*M_PI*xi2),sin(2.*M_PI*xi2));
}

Sphere *GetSphere(vec_type R, const Vector_3 &center){
  return new Sphere(R,center);
}

Plane_3* CreateSpherePlanes(const int &N, const vec_type &R, const Vector_3 &center){
  Plane_3 *SpherePlanes=new Plane_3[N+1];
  for (int i=0;i<N;i++) {
    Vector_3 n=R*random_direction();
    Vector_3 pos=center+n;
    SpherePlanes[i].init(-n,pos);
  }
  return SpherePlanes;
}

bool get_first_corner(const Box &box,const Plane_3 &plane,const Vector_3 &k,Vector_3 *start){
  VecContour<> cnt;
  if(ProjectSimplex(plane, box.planes_begin(), box.planes_end(), cnt)<0)
    return false;
  *start=get_min_point(k,cnt.points_begin(),cnt.points_end());
  return true;
}

bool get_first_corner(const Box &box,const Vector_3 &k,Vector_3 *start){
  if(!box.valid())
    return false;
  *start=get_min_point(k,box.points_begin(),box.points_end());
  return true;
}

bool find_cross_segment(const Region_3 &reg, const Vector_3& orig,const Vector_3& dir,vec_type &t1, vec_type &t2){
  Vector_3 tdir=dir;
  if(reg.TestPoint(orig)){  // the point is inside
    t2=reg.TestRay(orig,tdir);
    if(t2<0) //strange, no intersection 
      return false; 
    t1=reg.TestRay(orig,-tdir);
    if(t1<0) //strange, no intersection 
      return false;
    t1=-t1;
    return true;
  }
  int inv1=1, inv2=1;
  t1=reg.TestRay(orig,tdir);
  if(t1<0){
    t1=reg.TestRay(orig,-tdir);
    if(t1<0)
      return false;
    tdir=-tdir;
    inv1=-1;
  }
  Vector_3 p1=orig+tdir*t1;
  Vector_3 p2=p1+tdir; // test point
  if(reg.TestPoint(p2)) // test point is inside, must go in the same direction
    t2=reg.TestRay(p2,tdir);
  else{ // opposite direction
    t2=reg.TestRay(p2,-tdir);
    inv2=-1;
  }
  if(t2<0)
    t2=t1;
  else
    t2=t1+1+inv2*t2;
 
  if(inv1<0){
    t1=-t1;
    t2=-t2;
    swap(t1,t2);
  }
  return true;
}

int fill_lattice(const Region_3 &reg, const Vector_3& orig,const Basis_3 &cell,vector<Vector_3> &points, vector<iVector_3> *ipoints){
  Vector_3 v1, v2, cnt;
  cnt=reg.GetBoundingBox(&v1,&v2);// always starting from region center
  if(v1.infinite() || v2.infinite())
    return -2; // can't fill unbounded region
  Vector_3 sh=cell.inv(cnt-orig);
  
  /*Box bb0(v1,v2);
  RegDumper *box_dmp=GetRegionDumper(&bb0);
  vector<VecContour> dmppoints;
  box_dmp->Dump(dmppoints);
  FILE *f1=fopen("lim.pol","wt");
  DumpVecContours(f1,dmppoints);
  fclose(f1);
  delete box_dmp;*/

  // finding limiting projection
  Box bb(v1-orig,v2-orig);
  int jmax=numeric_limits<int>::min();
  int jmin=numeric_limits<int>::max(); 
  int imin=jmin, imax=jmax;
  for(Box::point_it pit=bb.points_begin();pit!=bb.points_end();++pit){
    v1=cell.inv(*pit);
    imax=max((int)ceil(v1[0]),imax);
    imin=min((int)floor(v1[0]),imin);
    jmax=max((int)ceil(v1[1]),jmax);
    jmin=min((int)floor(v1[1]),jmin);
  } 

  /*jmax=20;
  jmin=-20;
  imax=20;
  imin=-20;*/
  /*f1=fopen("tst.d","wt");
  FILE *f2=fopen("tstc.d","wt");
  FILE *f3=fopen("tstn.d","wt");*/


  int i[2]={(int)floor(sh[0]),(int)ceil(sh[0])};
  if(i[0]==i[1])
    i[1]++;
  for(int id=0;id<2;id++){
    int ncrossi=0;
    do{
      int ncrossj=0;
      int j[2]={(int)floor(sh[1]),(int)ceil(sh[1])};
      if(j[0]==j[1])
        j[1]++;
      for(int jd=0;jd<2;jd++){
        do{
          Vector_3 cntj=orig+i[id]*cell[0]+j[jd]*cell[1]+sh[2]*cell[2];
          vec_type t1,t2;
          if(!find_cross_segment(reg,cntj,cell[2],t1,t2)){ //no intersection
            if(ncrossj)
              break;
          }
          else{
            ncrossj++;
            
            /*Vector_3 v1=orig+i[id]*cell[0]+j[jd]*cell[1]+(t1+sh[2])*cell[2];
            fprintf(f2,"%g %g %g \"%d %d\"\n",v1[0],v1[1], v1[2],i[id],j[jd]);
            v1=orig+i[id]*cell[0]+j[jd]*cell[1]+(t2+sh[2])*cell[2];
            fprintf(f2,"%g %g %g \"%d %d\"\n\n\n",v1[0],v1[1], v1[2],i[id],j[jd]);
            fflush(f2);*/

            int ks=(int)ceil(t1+sh[2]), ke=(int)floor(t2+sh[2]);
            for(int k=ks;k<=ke;k++){
              Vector_3 atompos=orig+i[id]*cell[0]+j[jd]*cell[1]+k*cell[2];
              points.push_back(atompos);
              if(ipoints)
                ipoints->push_back(iVector_3(i[id],j[jd],k));
            }
          }
          if(jd){
            j[jd]++;
            if(j[jd]>jmax)
              break;
          }
          else{
            j[jd]--;
            if(j[jd]<jmin)
              break;
          }
        }while(1);
      }// jd
      if(ncrossj==0 && ncrossi!=0)
        break;  // went out of intersection region
      ncrossi+=ncrossj;
      if(id){
        i[id]++;
        if(i[id]>imax)
          break;
      }
      else{
        i[id]--;
        if(i[id]<imin)
          break;
      }

    }while(1);
  }
  //fclose(f1);
  return (int)points.size();
}

class CompareNorms{
public:
  bool operator()(const Vector_3 a, const Vector_3 b){
    return a.norm2()<b.norm2();
  }
};

Vector_3 *make_reciprocal_vectors(const vec_type a1, const vec_type a2, const vec_type max, int &num){

  vec_type kmax=1/max;

  vec_type a[2]={a1,a2};
  int n[2];
  for(int i=0;i<2;i++){
    a[i]=1/a[i];
    n[i]=int(floor(kmax/a[i]));
  }
  num=(2*n[0]+1)*(2*n[1]+1);
  Vector_3 *v=new Vector_3[num];
  int k=0;
  for(int i=-n[0];i<=n[0];i++){
    for(int j=-n[1];j<=n[1];j++){
      v[k++]=Vector_3(a[0],0,0)*i+Vector_3(0,a[1],0)*j;
    }
  }
  sort(v,v+num,CompareNorms());

  return v;
}

template class Box_N<2>;
template class Box_N<3>;
template class Sphere_N<2>;
template class Sphere_N<3>;

template class Polyhedron<Box::plane_it>;
template class Polyhedron<Plane_3 *>;

template class Cylinder<Circle>;
template class Cone<Circle>;
template class StretchedRegion<Region_3>;
template class ConfinedRegion<Region_3>;

template vec_type Box::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
//template vec_type Box::TestContour(const Contour_N<modified_value_it<indexed_it<std::vector<Vector_3>::const_iterator, std::vector<size_t>::const_iterator>, plus_t<Vector_3> > > &, VecContour<> *, Vector_3 *) const;
template vec_type Sphere::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
template vec_type Polyhedron<Box::plane_it>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
template vec_type Polyhedron<Plane_3 *>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;


template vec_type Circle::TestContour(const PtrContour<2> &cnt, VecContour<2> *subcont, Vector_2 *subcenter) const;

template vec_type Cylinder<Circle>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
//template vec_type Cone<Circle>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
template vec_type StretchedRegion<Region<3> >::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;

template class ConfinedRegion<Sphere>;
template vec_type ConfinedRegion<Sphere>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;

template class ConfinedRegion<Cylinder<Circle> >;
template vec_type ConfinedRegion<Cylinder<Circle> >::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;

template class ConfinedRegion<Cone<Circle> >;

template int MonteCarloTestLine(const Region<3> *reg, const Vector_Nt<vec_type,3> &p, const Vector_Nt<vec_type,3> &dir, vec_type *frac, 
Vector_Nt<vec_type,3> *surfp, Vector_Nt<vec_type,3> *surfn, vec_type length, int nt);
template vec_type MonteCarloTestContour(const Region<2> *reg, const VecContour<2> &cnt, Vector_Nt<vec_type,2> *subcenter, int nt);
template vec_type MonteCarloTestContour(const Region<3> *reg, const VecContour<3> &cnt, Vector_Nt<vec_type,3> *subcenter, int nt);
