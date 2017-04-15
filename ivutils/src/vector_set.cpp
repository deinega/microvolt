/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 1995-2012        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   $Revision: 1.2 $
 *   $Date: 2013/06/29 18:49:21 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/vector_set.cpp,v 1.2 2013/06/29 18:49:21 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/vector_set.cpp,v $
$Revision: 1.2 $
$Author: lesha $
$Date: 2013/06/29 18:49:21 $
*/
/*s****************************************************************************
 * $Log: vector_set.cpp,v $
 * Revision 1.2  2013/06/29 18:49:21  lesha
 * bug is fixed
 *
 * Revision 1.1  2012/12/06 02:24:06  lesha
 * *** empty log message ***
 *
 * Revision 1.15  2012/05/18 13:45:15  valuev
 * fixed documentation, box sides by flux calculation
 *
 * Revision 1.14  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.13  2012/04/12 20:10:57  lesha
 * documentation
 *
 * Revision 1.12  2012/03/22 07:31:20  valuev
 * updated increments for iterators, prepared parallel/universal near-to-far (not working)
 *
*******************************************************************************/

#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include "vector_set.h"

Vector_3 get_iterator_dx(const SpaceVectorSet::iterator &it){
  return it.dx();
}

BoxSurfaceSet::iterator::iterator(const BoxSurfaceSet *sparent, int e): parent(sparent), gi(e ? 6 : 0) {
  if(!e){
    while (!(parent->used & 1<<gi) && gi<6)
      gi++; // find first used grid
    if (gi<6) it=parent->grd[gi].begin();
  }
}

Vector_3 BoxSurfaceSet::iterator::ds() const {
  int dir=gi%3; // direction perpendicular to the surface
  Vector_3 n;
  int sign=(gi/3) ? 1 : -1; // surface orientation
  n[dir]=sign*parent->ds[dir]; // oriented surface element
  return n;
}

BoxSurfaceSet::iterator& BoxSurfaceSet::iterator::operator++(){ //prefix
  ++it;
  if (!(it!=parent->grd[gi].end())) {
    do {
      gi++; // next grid
    } while (!(parent->used & 1<<gi) && gi<6);
    if (gi<6) it=parent->grd[gi].begin();
  }
  return *this;
}

void BoxSurfaceSet::init(const Vector_3 &v1, const Vector_3&v2, const int *sz, const int sused) {
  for (int dir=0; dir<3; dir++) {
    int dir1=(dir+1)%3, dir2=(dir+2)%3;
    Vector_3 dv;
    dv[dir1]=(v2[dir1]-v1[dir1])/sz[dir1];
    dv[dir2]=(v2[dir2]-v1[dir2])/sz[dir2];
    ds[dir]=dv[dir1]*dv[dir2]; // area of elementary surface at this side
    // constructing grids for opposite sides perpendicular to dir
    // grids size
    int s[3];
    s[dir]=1;
    s[dir1]=sz[dir1];
    s[dir2]=sz[dir2];
    // grids position
    Vector_3 shift;
    shift[dir1]=sz[dir1]>1 ? dv[dir1]/2 : 0;
    shift[dir2]=sz[dir2]>1 ? dv[dir2]/2 : 0;
    Vector_3 w1=v1, w2=v1;
    w1[dir]=v2[dir];
    w2[dir1]=v2[dir1];
    w2[dir2]=v2[dir2];
    grd[dir].init(v1+shift, w2-shift, -(dir+1), s);
    grd[dir+3].init(w1+shift, v2-shift, dir+1, s);
  }
  used=sused;
}

size_t BoxSurfaceSet::size() {
  size_t sz(0);
  for(int i=0;i<6;i++)
    if(used & (1<<i))
      sz+=grd[i].size();
  
  return sz;
}

Vector_3 BoxSurfaceSet::get_total_surface(){
  Vector_3 surface;
  for(int i=0; i<3; i++){
    if(used&(1<<i)){
      if(!(used&(1<<(i+3))))
        surface+=grd[i].get_total_surface();
    }
    else if(used&(1<<(i+3)))
      surface+=grd[i+3].get_total_surface(); // take opposite size!
  }
  return surface;
}

Vector_3 CylinderSurfaceSet::iterator::operator* ()const {
  int znum=parent->znum, finum=parent->finum;
  vec_type R=parent->R;
  Vector_3 origin=parent->origin, n=parent->n, x=parent->x, y=parent->y;
  vec_type phi=vec_type(fi)/finum*2*M_PI; // radial angle
  return origin+n*(zi+.5)/znum+R*(x*cos(phi)+y*sin(phi));
}

Vector_3 CylinderSurfaceSet::iterator::ds ()const {
  int finum=parent->finum;
  Vector_3 x=parent->x, y=parent->y;
  vec_type ds=parent->ds;
  vec_type phi=vec_type(fi)/finum*2*M_PI; // radial angle
  return ds*(x*cos(phi)+y*sin(phi));
}

void CylinderSurfaceSet::init(const Vector_3 &origin_, const Vector_3 &n_, const Vector_3 &x_, const Vector_3 &y_, 
vec_type R_, const int *sz, int surf_dir){
  origin=origin_, n=n_, x=x_, y=y_, R=R_;
  finum=sz[0], znum=sz[1];
  ds=(2*M_PI*R*n.norm())/(finum*znum);
  if(surf_dir<0)
    ds=-ds;
}

Vector_3 SphereSurfaceSet::iterator::operator*() const {
  Vector_3 Rtf=this->internal_coords();
  vec_type &R=Rtf[0];
  vec_type &theta=Rtf[1];
  vec_type &fi=Rtf[2];
  return parent->origin+Vector_3(R*sin(theta)*cos(fi), R*sin(theta)*sin(fi), R*cos(theta));
}

Vector_3 SphereSurfaceSet::iterator::internal_coords() const {
  vec_type theta=parent->theta0+parent->dtheta*th;
  vec_type fi=parent->fi0+parent->dfi*f;
  return Vector_3(parent->R, theta, fi);
}

SphereSurfaceSet::iterator& SphereSurfaceSet::iterator::operator++() {
  f++;
  if (f>=parent->fi_num) {
    f=0;
    th++;
  }
  return *this;
}

int SphereSurfaceSet::init(const Vector_3 &origin_, vec_type R_, int theta_num_, int fi_num_, vec_type theta0_, vec_type fi0_, vec_type theta1, vec_type fi1) {
  origin=origin_;
  R=R_;

  if (theta0_>theta1) return -1;
  if (fi0_>fi1) return -1;
  if (theta_num_<1) return -1;
  if (fi_num_<1) return -1;

  if (theta0_<0) theta0_=0;
  if (theta1>M_PI) theta1=M_PI;
  if (fi0_<0) fi0_=0;
  if (fi1>2*M_PI) fi1=2*M_PI;

  theta0=theta0_;
  fi0=fi0_;

  theta_num=theta_num_;
  fi_num=fi_num_;

  if (theta_num>1)
    dtheta=(theta1-theta0)/(theta_num-1);

  if (fi_num>1)
    dfi=(fi1-fi0)/(fi_num-1);

  return 1;
}

Vector_3 internal_coords(SphereSurfaceSet::iterator it) {
//  return it.internal_cooeds();
  Vector_3 coords=it.internal_coords();
  coords[1]*=180/M_PI;
  coords[2]*=180/M_PI;
  coords[1]=fabs(coords[1]);
  coords[2]=fabs(coords[2]);
  return coords;
};

Vector_3 internal_coords(const Vector_3 &F, SphereSurfaceSet::iterator it) {
  Vector_3 pos=it.internal_coords();
  vec_type &theta=pos[1];
  vec_type &fi=pos[2];
  return Vector_3(F[0]*sin(theta)*cos(fi)+F[1]*sin(theta)*sin(fi)+F[2]*cos(theta), 
                  F[0]*cos(theta)*cos(fi)+F[1]*cos(theta)*sin(fi)-F[2]*sin(theta), 
                 -F[0]*sin(fi)+F[1]*cos(fi));
};
