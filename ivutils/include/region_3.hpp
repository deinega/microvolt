/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *   $Revision: 1.13 $
 *   $Date: 2013/11/25 01:28:57 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/region_3.hpp,v 1.13 2013/11/25 01:28:57 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/region_3.hpp,v $
$Revision: 1.13 $
$Author: lesha $
$Date: 2013/11/25 01:28:57 $
*/
/*s****************************************************************************
 * $Log: region_3.hpp,v $
 * Revision 1.13  2013/11/25 01:28:57  lesha
 * MakeBoxDomain is fixed
 *
 * Revision 1.12  2013/11/24 18:15:49  lesha
 * MakeBoxDomain is fixed for float precision
 *
 * Revision 1.11  2013/11/23 22:38:44  lesha
 * float is supported
 *
 * Revision 1.10  2013/10/28 20:46:44  lesha
 * volume_integral_const_t is moved back
 *
 * Revision 1.9  2013/10/27 23:13:36  lesha
 * dx is introduced to MakeBoxDomain
 *
 * Revision 1.8  2013/10/26 20:25:19  lesha
 * Integrate is added to grid
 *
 * Revision 1.7  2013/10/26 16:57:40  lesha
 * bugs are fixed
 *
 * Revision 1.6  2013/10/26 14:13:06  valuev
 * modified MAkeBoxDomains
 *
 * Revision 1.5  2013/10/25 06:49:26  lesha
 * balanced domain decomposition interface
 *
 * Revision 1.4  2013/10/25 06:21:52  lesha
 * documentation
 *
 * Revision 1.3  2013/10/01 01:11:28  lesha
 * GetAreaFraction is simplified
 *
 * Revision 1.2  2013/09/30 22:18:13  lesha
 * GetAreaFraction is moved here
 *
 * Revision 1.1  2013/04/23 01:06:20  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2013/02/27 12:45:55  valuev
 * confined TestLine: fixed touching case
 *
 * Revision 1.3  2013/02/27 10:58:02  valuev
 * confined TestLine: added surface points for infinite solutions
 *
 * Revision 1.2  2013/02/27 10:48:54  valuev
 * TestLine for confined region covering all cases
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.49  2012/09/27 19:41:35  lesha
 * acccomp -> accomp
 *
 * Revision 1.48  2012/09/24 23:09:43  lesha
 * MonteCarlo implementation is moved from region_2 to region
 *
 * Revision 1.47  2012/09/24 21:51:10  lesha
 * documentation
 *
 * Revision 1.46  2012/09/24 13:44:16  lesha
 * pencil, seqpack and table_function are excluded from transport
 *
 * Revision 1.45  2012/08/23 05:24:47  lesha
 * small is intreduced to TestLine
 *
 * Revision 1.44  2012/04/12 20:10:56  lesha
 * documentation
 *
 * Revision 1.43  2012/03/21 17:23:21  lesha
 * documentation
 *
 * Revision 1.42  2012/01/10 01:07:12  lesha
 * SurfProject is added
 *
 * Revision 1.41  2011/09/29 07:12:02  lesha
 * ConfinedRegion::TestLine - other case is considered
 *
 * Revision 1.40  2011/03/18 03:34:14  lesha
 * TestLine is realized for ConfinedRegion (not fully)
 *
*******************************************************************************/

/// \file \brief Template function definitions for region.h

#include "region_3.h"

#if !defined(DUMP)
template<class reg_t>
RegDumper *GetRegionDumper(const reg_t *reg){
  return NULL;
}
#endif

// countour with functon point which return random point inside the contour
template<class contour_t,int N>
class MCContour{};

// specialization of MCContour for 2D
template<class contour_t>
class MCContour<contour_t,2>{
  Vector_2 p1,sz;
  contour_t const *cnt;
public:
  MCContour():cnt(NULL){}
  MCContour(const contour_t &cnt_){init(cnt_);}
  void init(const contour_t &cnt_){
    cnt=&cnt_;
    cnt->GetBoundingBox(&p1,&sz);
    sz-=p1;
  }
  // random point inside the contour
  Vector_2 point(){
    Vector_2 p;
    do{
      for(int i=0;i<2;i++){
        p[i]=p1[i]+((vec_type)rand())/RAND_MAX*sz[i];
      }
    }while(!cnt->TestPoint(p));
    return p;
  }
};

// specialization of MCContour for 3D
template<class contour_t>
class MCContour<contour_t,3>{
  Contour_3to2<contour_t> cnt2;
  MCContour<Contour_3to2<contour_t>,2> MC2;
  Vector_3 c,x,y;
public:
  MCContour(const contour_t &cnt_){
    c=*(cnt_.points_begin());
    Vector_3 n=cnt_.GetNormVect();
    n.normalize();
    build_orth_basis(n,x,y);
    cnt2=Contour_3to2<contour_t>(cnt_,c,x,y);
    MC2=MCContour<Contour_3to2<contour_t>,2>(cnt2);
  }
  // random point inside the contour
  Vector_3 point(){
    Vector_2 point2=MC2.point();
    Vector_3 p=Vector_2to3(point2,c,x,y);
    return p;
  }
};

template<int N,class contour_t>
vec_type MonteCarloTestContour(const Region<N> *reg, const contour_t &cnt, Vector_Nt<vec_type,N> *subcenter, int nt){
  MCContour<contour_t,N> c(cnt);
  int in=0;
  if(subcenter)
    *subcenter=0;
  for(int t=0;t<nt;t++){
    Vector_Nt<vec_type,N> point=c.point();
    bool test=reg->TestPoint(point);
    if(test){
      in++;
      if(subcenter)
        *subcenter+=point;
    }
  }
  if(subcenter)
    *subcenter/=in;
  return vec_type(in)/vec_type(nt)*cnt.Area();
}

template<int N>
int MonteCarloTestLine(const Region<N> *reg, const Vector_Nt<vec_type,N> &p, const Vector_Nt<vec_type,N> &dir, vec_type *frac, 
Vector_Nt<vec_type,N> *surfp, Vector_Nt<vec_type,N> *surfn, vec_type length, int nt){
  vec_type dir_norm=dir.norm();
  bool step=reg->TestPoint(p-length/dir_norm*dir); // left side
  int ind=0;
  vec_type dr=length/(vec_type(nt)*dir_norm);
  for(int i=-nt+1;i<=nt;i++){
    vec_type dist=i*dr;
    Vector_Nt<vec_type,N> pos=p+dist*dir;
    if(reg->TestPoint(pos)!=step){
      if(ind>=2)
        return ind;
//         throw -1;
      for(int i1=1;i1<=nt;i1++){
        vec_type dist=(i-1+vec_type(i1)/nt)*dr;
        pos=p+dist*dir;
        if(reg->TestPoint(pos)!=step)
          break;
      }
      if(surfp)
        surfp[ind]=pos;
      if(surfn){
        Vector_Nt<vec_type,N> bs[N], nbs[N];
        bs[0]=dir;
        if(!build_orth_basis(bs,1,bs+1)) // make basis, first eletemnt is dir
          throw -1;
        int cnum=100;
        vec_type drd=dr*dir_norm;
        for(int di=1;di<N;di++){
          bool was=reg->TestPoint(pos+drd*bs[0]);
          int dind=0;
          Vector_Nt<vec_type,N> nm[2];
          int cival[2]={-1,-1};
          for(int ci=1;ci<=cnum;ci++){
            vec_type phi=ci*2.*M_PI/cnum;
            Vector_Nt<vec_type,N> tp = bs[0]*cos(phi) + bs[di]*sin(phi);
            if(reg->TestPoint(pos+drd*tp)!=was){
              int cnum1=100;
              for(int ci1=1;ci1<=cnum1;ci1++){
                phi=(ci-1+vec_type(ci1)/cnum1)*2.*M_PI/cnum;
                tp = bs[0]*cos(phi) + bs[di]*sin(phi);
                if(reg->TestPoint(pos+drd*tp)!=was)
                  break;
              }
              was=!was;
              if(dind>=2)
                throw -1;
              cival[dind]=ci;
              nm[dind++]=tp;
            }
          }
          if(dind<=1)
            throw -1;
          nbs[di-1]=nm[1]-nm[0]; // tangent to surface
        }
        if(!build_orth_basis(nbs,N-1,nbs+N-1))
          throw -1;
        surfn[ind]=nbs[N-1]; // normal to all tangents is normal to surface
        surfn[ind].normalize();
        if(reg->TestPoint(pos+drd*surfn[ind]))
          surfn[ind]*=-1;
      }
      if(frac){
        frac[ind]=dist;
      }
      ind++;
      step=!step;
    }
  }
  return ind;
}

template <class contour_t>
vec_type Box::TestContour(const contour_t &cnt, VecContour<3> *subcont, vector_t *subcenter) const{
  Polyhedron<plane_it> poly(planes_begin(),planes_end());
  return poly.TestContour(cnt,subcont,subcenter);
}

template <class plane_itt>
void Polyhedron<plane_itt>::init(plane_it beg, plane_it end){
  b=beg;
  e=end;
  /// checking for coordinate planes
  Vector_3 n, p1(-VEC_INFTY), p2(VEC_INFTY);
  vec_type d;
  plane_it it;
  int i;
  for(it=b;it!=e;it++){
    const Plane_3 &p=*it;
    p.get_coeff(n,d);
    for(i=0;i<3;i++){
      if(fabs(fabs(n[i])-1.)<VEC_ZERO){ // coordinate direction
        break;
      }
    }
    if(i<3){ // adding as bound
      if(n[i]>0){ // left bound
        if(p1[i]<-d)p1[i]=-d;
      }
      else{  // right bound
        if(p2[i]>d)p2[i]=d;
      }
    } 
  }
  box.init(p1,p2);
}

template <class plane_itt>
bool Polyhedron<plane_itt>::TestPoint(const Vector_3 &p) const{
  // checking bounding box
  if(!box.TestPoint(p))return false;
  plane_it it;
  int res=1;
  for(it=b;it!=e;it++){
    if((*it).distance(p)<0){
      res=0;
      break;
    }
  }
  if(!res)return false;
  return true;
}

template <class plane_itt>
vec_type Polyhedron<plane_itt>::MinPlaneDist(const Vector_3 &pos, plane_it *mit){
  plane_it it=planes_begin(), e=planes_end();
  vec_type md=-1.;
  for(it=b;it!=e;it++){
    vec_type d=fabs((*it).distance(pos));
    if(md<0 || d<md){
      md=d;
      if(mit)*mit=it;
    }
  }
  return md;
}

template <class plane_itt>
template <class contour_t>
vec_type Polyhedron<plane_itt>::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{
  /// testing bounding box
  Vector_3 bc1, bc2;
  cnt.GetBoundingBox(&bc1,&bc2);
  Vector_3 p1=box.get_p1(), p2=box.get_p2();
  for(int i=0;i<3;i++){
    if(bc1[i]>p2[i] || bc2[i]<p1[i])return 0.;
  }
  aggregate_it<plane_it, typename contour_t::plane_it, Plane_3> a_b(b,e,cnt.planes_begin()),a_e(e,e,cnt.planes_end());
  // solving for intersection
  VecContour<> s_cont;
  int res=ProjectSimplex(cnt.GetPlane(),a_b,a_e,s_cont);
  if(res==1){
    vec_type a=s_cont.Area();
    if(subcenter)*subcenter=s_cont.GetCenter();
    if(subcont)subcont->GetPoints().swap(s_cont.GetPoints());
    return a; 
  }
  return 0.;
}

template <class plane_itt>
int Polyhedron<plane_itt>::TestLine(const vector_t &p, const vector_t &dir, 
vec_type *frac, vector_t *surfp, vector_t *surfn, vec_type epsilon) const{

  double small=VEC_ZERO;

  vec_type tleft=-numeric_limits<vec_type>::max(), tright=numeric_limits<vec_type>::max();
  Vector_3 nleft, nright;
  
  for(plane_it pi=b;pi!=e && tleft<=tright ;++pi){
    Vector_3 vnorm=(*pi).get_normal();
    vec_type dprod=-dir*vnorm; // '-' takes the opposite normal direction in ivutils into account
    // dprod>0 - я вхожу в плоскость изнутри, dprod<0 - снаружи
    vec_type dist=(*pi).distance(p); // distance to the plane

    if(fabs(dprod)<small){ // dir is parallel to the plane
      if(dist>=0)
        continue; // dir is all inside, OK
      else
        return 0; // dir is all outside, no solution
    }
    vec_type t=dist/dprod;
    if(dprod>0){ // right bound, searching minimum
      if(tright>t){
        tright=t;
        nright=vnorm;
      }
    }
    else{// left bound, searching maximum
      if(tleft<t){
        tleft=t;
        nleft=vnorm;
      }
    }
  }
  if(tleft>tright)
    return -1; // no solution on this line
  frac[0]=tleft;
  frac[1]=tright;

  if(surfp){
     surfp[0]=p+dir*tleft;
     surfp[1]=p+dir*tright;
  }
  if(surfn){
    surfn[0]=-nleft;
    surfn[1]=-nright;
  }
  return 1+(tleft!=tright);
}

template <class plane_itt>
vec_type Polyhedron<plane_itt>::SurfProject(const vector_t &p, vector_t *surfp, vector_t *surfn) const{

  bool in=TestPoint(p);

  vec_type t = in ? numeric_limits<vec_type>::max() : -numeric_limits<vec_type>::max();
  Vector_3 n;
  
  for(plane_it pi=b;pi!=e;++pi){
    Vector_3 vnorm=(*pi).get_normal();
    // dist>0 - я вхожу в плоскость изнутри, dist<0 - снаружи
    vec_type dist=(*pi).distance(p); // distance to the plane

    if(in && dist>=0){ // right bound, searching minimum
      if(t>dist){
        t=dist;
        n=vnorm;
      }
    }
    else if(!in && dist<0){// left bound, searching maximum
      if(t<dist){
        t=dist;
        n=vnorm;
      }
    }
  }

  if(surfp)
     *surfp=p+n*t;

  if(surfn)
    *surfn=-n;

  return t;
}

template <class contour_t>
vec_type Sphere::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{

  if(subcenter)
    *subcenter=0;

  Vector_3 nn=cnt.GetNormVect(); // normal to contour
  nn.normalize();
    
  vec_type dist=(*(cnt.points_begin())-center)*nn; // "distance" from center of sphere to contour plane (can be negative)
  vec_type r2=(R*R-dist*dist); // radius^2 of this circle forming by intersection between contour plane and sphere
  if (r2<=0) return 0; // contour plane does not intersect sphere

  Vector_3 c=center+dist*nn; // center of circle 
  Circle C(sqrt(r2),Vector_2());
  Vector_3 ax,ay;
  build_orth_basis(nn,ax,ay);
  Contour_3to2<contour_t> cnt2(cnt,c,ax,ay);
  Vector_2 subcenter2;
  vec_type area=C.TestContour(cnt2,NULL,subcenter ? &subcenter2 : NULL);
  if(subcenter && area)
    *subcenter=Vector_2to3(subcenter2,c,ax,ay);

  return area;
}

template<class base_tt>
template<class contour_t>
vec_type Cylinder<base_tt>::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{
  if(cnt.GetNPoints()==0)
    return 0;

  Vector_3 ncnt=cnt.GetNormVect();
  ncnt.normalize();
  vec_type coef=ncnt*n;

  if(accomp(coef,vec_type(0))){ // контур параллелен оси циллиндра, его нельзя спроектировать на основание
    Vector_3 k=ncnt%n,b=*(cnt.points_begin());
    Vector_2 k2=::Vector_3to2(k,Vector_3(),x,y),b2=Vector_3to2(b);
    vec_type frac[2];
    Vector_2 surfp[2];
    int tl=base->TestLine(b2,k2,frac,surfp);
    if(tl<=1)
      return 0; // all outside, no intersection
    Plane_3 planes[3];
    planes[0].init(k, Vector_2to3(surfp[0]));
    planes[1].init(-k, Vector_2to3(surfp[1]));
    Polyhedron_3 poly(planes, planes+2, 0);
    return poly.TestContour(cnt, NULL, subcenter);
  }

  // контур проецируется на основание
  Contour_3to2<contour_t> cnt2(cnt,origin,x,y);
  Vector_2 subcenter2;
  vec_type area=base->TestContour(cnt2, NULL, subcenter ? &subcenter2 : NULL);
  area/=fabs(coef);
  if(area!=0 && subcenter){
    Vector_3 b=*(cnt.points_begin()); // начальная вершина контура
    vec_type hc=(b-origin)*n; // ее высота
    Vector_3 c=Vector_2to3(subcenter2)+hc*n; // переводит в 3D на уровне точки b
    Plane_3 pl(ncnt,b); // плоскость контура
    pl.TestRay(c,n,subcenter); // поднимаем в плоскость контура
  }
  return area;
}

template<class base_tt>
int Cylinder<base_tt>::TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
vector_t *surfp, vector_t *surfn, vec_type epsilon) const{
  Vector_3 v=dir%n;
//  if(accomp(v.norm2()))
  if(accomp(v.norm(),vec_type(0)))
    return 0; // parallel to axis

  Vector_2 p2=Vector_3to2(p), dir2=::Vector_3to2(dir,0,x,y);
  Vector_2 surfn2[2];
  int tl=base->TestLine(p2,dir2,frac,NULL,surfn ? surfn2 : NULL);
  for(int i=0;i<tl;i++){
    if(surfp)
      surfp[i]=p+frac[i]*dir;
    if(surfn)
      surfn[i]=surfn2[i][0]*x+surfn2[i][1]*y;
  }
  return tl;
}

template<class base_tt>
int Cone<base_tt>::TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
vector_t *surfp, vector_t *surfn, vec_type epsilon) const{

  // we solve square equation for t:
  // (p + t*dir - origin, x)^2 + (p + t*dir - origin, y)^2 = (R * (1 - (p + t*dir - origin, n)/L)))^2

  return -1;
}

template<class reg_tt>
template<class contour_t>
vec_type StretchedRegion<reg_tt>::TestContour(const contour_t &cnt, VecContour<reg_tt::dimension> *subcont, typename reg_tt::vector_t *subcenter) const{
  VecContour<reg_tt::dimension> cnt1,subcont1;
  vector_t subcenter1;
  for(typename contour_t::point_iterator it=cnt.points_begin(),e=cnt.points_end();it!=e;++it)
    cnt1.add(basis_inv(*it-shift));
  vec_type area=reg->TestContour(cnt1,subcont ? &subcont1 : NULL,subcenter ? &subcenter1 : NULL);

  Vector_3 nn=cnt.GetNormVect(),x,y;
  build_orth_basis(nn,x,y);
  Vector_3 x1=basis(x),y1=basis(y);
  vec_type change=sqrt((x%y).norm2()/(x1%y1).norm2());
  area*=change;

  if(subcont){
    subcont->Clear();
    for(int i=0;i<subcont1.GetNPoints();i++)
      subcont->add(shift+basis(subcont1.GetPoint(i)));
  }
  if(subcenter){
    *subcenter=shift+basis(subcenter1);
  }

  return area;
}

template<class reg_tt>
template<class contour_t>
vec_type ConfinedRegion<reg_tt>::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{
  VecContour<> subcnt;
  vec_type area=poly->TestContour(cnt, &subcnt, NULL);
  return area ? reg->TestContour(subcnt,subcont,subcenter) : 0;
}

# if 0  // old version

template<class reg_tt>
int ConfinedRegion<reg_tt>::TestLine(const Vector_3 &p, const Vector_3 &dir, vec_type *frac,
Vector_3 *surfp, Vector_3 *surfn, vec_type epsilon) const{
  vec_type frac2[2][2];
  vector_t surfp2[2][2],surfn2[2][2];
  int res[2];
  bool in[2];
  Region<dimension> *r[2]={reg.ptr(),poly.ptr()};
  for(int i=0;i<2;i++){
    res[i]=r[i]->TestLine(p,dir,frac2[i],surfp2[i],surfn2[i]);
    in[i]=r[i]->TestPoint(p);
  }
  if(in[0] && in[1]){ // inside
    vec_type min=-VEC_INFTY, max=VEC_INFTY;
    int ii[2]={-1,-1,},ij[2]={-1,-1};
    for(int i=0;i<2;i++){
      for(int j=0;j<res[i];j++){
        if(frac2[i][j]>min && frac2[i][j]<0){
          min=frac2[i][j];
          ii[0]=i, ij[0]=j;
        }
        if(frac2[i][j]<max && frac2[i][j]>0){
          max=frac2[i][j];
          ii[1]=i, ij[1]=j;
        }
      }
    }
    int ind=0;
    for(int i=0;i<2;i++){
      if(ii[i]>=0){
        frac[ind]=frac2[ii[i]][ij[i]];
        if(surfp)
          surfp[ind]=surfp2[ii[i]][ij[i]];
        if(surfn)
          surfn[ind]=surfn2[ii[i]][ij[i]];
        ind++;
      }
    }
    return ind;
  }
  else if(!in[0] && !in[1]){ // outside both regions... not tested!
    for(int i=0;i<2;i++){
      if(res[i]<=0) // line does not intersect at least one region
        return 0;
    }
    if((frac2[0][0]>=0 && frac2[1][0]<=0) || (frac2[1][0]>=0 && frac2[0][0]<=0)) // regions are at different sides of line
      return 0;

    int open = fabs(frac2[1][0]) > fabs(frac2[0][0]); // far opening region
    vec_type open_val = fabs(frac2[open][0]);
    int close=-1;
    vec_type close_val=VEC_INFTY;
    for(int i=0;i<2;i++){
      if(res[i]==2){
        if(fabs(frac2[i][1]) < close_val){
          close_val=fabs(frac2[i][1]);
          close=i;
        }
      }
    }
    if(open_val>=close_val)
      return 0;
    frac[0]=frac2[open][0];
    if(surfp)
      surfp[0]=surfp2[open][0];
    if(surfn)
      surfn[0]=surfn2[open][0];

    if(close==-1)
      return 1;

    frac[1]=frac2[close][1];
    if(surfp)
      surfp[1]=surfp2[close][1];
    if(surfn)
      surfn[1]=surfn2[close][1];

    return 2;
  }
  // other cases still are not programmed
  return 0;
}

# endif

# if 1  // new version


template<class reg_tt>
int ConfinedRegion<reg_tt>::TestLine(const Vector_3 &p, const Vector_3 &dir, vec_type *frac,
Vector_3 *surfp, Vector_3 *surfn, vec_type epsilon) const{
  vec_type frac2[2][2];
  vector_t surfp2[2][2],surfn2[2][2];
  Region<dimension> *r[2]={reg.ptr(),poly.ptr()};
  for(int i=0;i<2;i++){
    int res=r[i]->TestLine(p,dir,frac2[i],surfp ? surfp2[i] : NULL, surfn ? surfn2[i] : NULL);
    if(res<=1){ // line not crossing region or touching it
      if(TestPoint(p)){  
        if(res<=0){// line fully inside, making infinite limits
          frac2[i][0]=-VEC_INFTY;
          frac2[i][1]=VEC_INFTY;
          if(surfp){  // should not be used anyway ?
            surfp2[i][0] = frac2[i][0]*dir;
            surfp2[i][1] = frac2[i][1]*dir;
          }
        }
        else{ // res =1, touching, limited by single point
          frac2[i][1] = frac2[i][0]; // copying data to right limit
          if(surfp)
            surfp2[i][1]=surfp2[i][0];
          if(surfn)
            surfn2[i][1]=surfn2[i][0];
        }
      }
      else // line fully outside one of the regions, no crossing with region intersection
        return 0;
    }    
    if(frac2[i][0]>frac2[i][1]){ // test that we have all intersections ordered
      swap(frac2[i][0],frac2[i][1]);
      if(surfp)
        swap(surfp2[i][0],surfp2[i][1]);
      if(surfn)
        swap(surfn2[i][0],surfn2[i][1]);
    }
  }
  // analysing intersections
  if(frac2[0][0]>frac2[1][1] || frac2[1][0]>frac2[0][1] ) // no region intersection along the line
    return 0;
  // selecting left and right limit on the line
  int ileft=0, iright=0;
  if(frac2[0][0]<frac2[1][0])
    ileft=1;
  if(frac2[0][1]>frac2[1][1])
    iright=1;
  frac[0]= frac2[ileft][0];
  if(surfp)
    surfp[0]= surfp2[ileft][0];
  if(surfn)
    surfn[0]= surfn2[ileft][0];
  frac[1]= frac2[iright][1];
  if(surfp)
    surfp[1]= surfp2[iright][1];
  if(surfn)
    surfn[1]= surfn2[iright][1];
  
  return 1+(frac[0]!=frac[1]);
}
 
# endif





/*
template<class reg_tt>
vec_type ConfinedRegion<reg_tt>::TestEdge(const Vector_3 &p1, const Vector_3 &p2, Vector_3 *surfp, Vector_3 *surfn) const{
  Vector_3 surfp2[2], surfn2[2];
  vec_type res[2];
  res[0]=reg->TestEdge(p1, p2, surfp2, surfn2);
  res[1]=poly->TestEdge(p1, p2, surfp2+1, surfn2+1);
  int ires;
  if (res[0]==1 && res[1]==1)
    return 1;
  if (res[0]==-1 || res[1]==-1)
    return -1;
  if (res[0]==1)
    ires=1;
  else if (res[1]==1)
    ires=0;
  else {
    if (reg->TestPoint(surfp2[1]))
      ires=1;
    else
      ires=0;
  }
  if (surfp)
    *surfp=surfp2[ires];
  if (surfn)
    *surfn=surfn2[ires];
  return res[ires];
}
*/

template<class inp_it, class weight_integral_t >
int MakeBoxDomains(const Box &box,int nproc, int nx, int ny, int nz, inp_it result, const weight_integral_t  &weight_integral, const Vector_3 &dx){
  // this function calls itself recursively for eash next subbox 
  // until minimal subbox is reached (minimal subbox is domain recorded to result)
  // 3 directions are considered in the order starting from z-direction
  // the box is split by a bisecting plane perpendicular to the current splitting direction, 
  // each of the formed subboxes then subdivided again recursively
  int ndir[3]={nx,ny,nz}; // specified number of subsections for each directiion
  int nspl[3]={nx,ny,nz}; // these arguments will be substituted to the next recursive call of MakeBoxDomains
  // if nproc is positive, then the found domains will be recorded to result
  // otherwise return the number of formed domains
  int np=(nproc>0 ? nproc : -nproc);
  int startdir=2; // we start with z-direction, assuming that this direction is most preferable to be cut by sections
  int curdir=-1; // direction which will be considered at this recursive call
  int nauto=0; // number of automatic directions (where the count of splitting planes is chosen automatically)
  int autod[3]={-1, -1, -1}; // directions where the number of splitting planes will be chosen automatically
  int ndef=0; // number of skipped directions at the current recursive call (ndir[dir]=0)
  for(int dd=0;dd<3;dd++){ // looking for working direction
    int dir=(dd+startdir)%3;
    if(ndir[dir]>0 && (curdir<0 || ndir[dir]>ndir[curdir]))
      curdir=dir; // current direction is found
    else if(ndir[dir]<0){ // indicates automatic direction
      autod[nauto++]=dir;
    }
    if(ndir[dir]==0)ndef++; // indicates skipped direction
  }
  if(curdir<0 && nauto>0){
    // working direction is not found,
    // and there are directions where splits count should be found automatically
    if(nauto==1){ // this is the only working direction
      ndir[autod[0]]=np; // number of sections per this direction is trivially equal to np
      curdir=autod[0];
    }
    else{
      // there are more than one possible working (automatic) direction,
      // we will chose working direction with the largest box side (if it can be divided yet)
      vec_type lprod=1.; // box volume (area) for automatic directions
      for(int i=0;i<nauto;i++){
        lprod*=box.GetSize()[autod[i]];
      }
      vec_type k=pow((vec_type)np/lprod,vec_type(1)/vec_type(nauto)); // mean 1D size per one subbox
      int mnp=0, mdir=0;
      for(int i=0;i<nauto;i++){ // now finding the direction that can produce maximal number of splits
        int tnp=(int)acfloor((box.GetSize()[autod[i]]*k)); // integer number of splits per direction
        if(mnp<tnp){
          mnp=tnp; // maximal number of splits
          mdir=i; // direction which contains maximal number of splits
        }
      }
      if(mnp>0){ // we choose this direction as a working direction
        ndir[autod[mdir]]=mnp;
        curdir=autod[mdir];
      }
    }
  }
  if(curdir>=0){
    int ndom=0; // number of resulted domains
    // determining sub-nps.
    // THIS IS NP SWITCH : !!!
    if(np<ndir[curdir])ndir[curdir]=np;
    if(np>ndir[curdir] && ndef==2){
      // the last direction (all other directions are considered already and skipped at this recursive call) 
      ndir[curdir]=np;
    }

    int snp=(np>ndir[curdir] ? np : ndir[curdir]);
    int mnp=snp/ndir[curdir]; // number of splits per orthogonal direction, which will be considered in a next recursive call
    int dntot=snp%ndir[curdir];
    vec_type c0=(vec_type)mnp/snp, c1=(vec_type)(mnp+1)/snp;
    Vector_3 p1=box.get_p1();
    Vector_3 p2=box.get_p2();
    
    vec_type l=box.GetSize()[curdir];

#if 1
    vec_type total = weight_integral(box.get_p1(),box.get_p2()); //  total integral for this volume
    vector< pair<vec_type, int> > fraction; // pair is {required_integral_at_split_pos, number_of_procs}
    vec_type sum = 0.;
    // distributing snp processors between ndir[curdir] slices
    for(int i=0;i<ndir[curdir];i++){ // cutting subboxes for othogonal direction
      // if dntot>0, cutting larger subboxes in order to make final subdomains equal
      vec_type dl=(dntot>0 ? total*c1 : total*c0);
      sum+= dl;
      int cnp=(dntot>0 ? mnp+1 : mnp);
      if(cnp<=0)break;
      fraction.push_back(make_pair(sum,cnp));
      dntot--;
    }
    int snum = dx[curdir] ? int(ceil(acdiv(l/ndir[curdir],dx[curdir]))) : 100;
    // step for analysing integral in current direction
    vec_type step = l/ndir[curdir]/snum; 
    sum =0.;
    size_t split = 0;
    Vector_3 v1 = box.get_p1(), v2 = box.get_p2(); // integration limits
    int assigned_np = 0; // total number of processors assigned to subvolumes 
    for(int si=0;si<ndir[curdir]*snum;si++){
      v1[curdir] = box.get_p1()[curdir]+si*step;
      v2[curdir] = v1[curdir]+step;
      sum+= weight_integral(v1,v2); // integral of weight between x-step and x
      int cur_np = 0;
      while(split<fraction.size() && (acless(fraction[split].first,sum,fraction[split].first) || si>=ndir[curdir]*snum-1)){ // integral reached required value 
        cur_np += fraction[split].second;
        split++;
      }
      if(cur_np>0){  // this slice has processors to assign
        if(split >= fraction.size()){ // this is the last split, make it exactly aligned
          p2[curdir]=box.get_p2()[curdir];
          cur_np = snp - assigned_np; //assign exactly the number of processors left
        }
        else
          p2[curdir]=v2[curdir];
        Box subdiv(p1,p2);
        nspl[curdir]=0; // this direction is frozen for the next recursive call, where orthogonal direction will be split
        // recursive call; subboxes will be recorded to result, only if nproc>0
        int ddom=MakeBoxDomains(subdiv,(nproc>0 ? cur_np: -cur_np),nspl[0],nspl[1],nspl[2],result,weight_integral,dx);
        assigned_np +=  cur_np;
        if(nproc>0)
          result+=ddom;
        ndom+=ddom;
        p1[curdir]=p2[curdir];
      }
    }
# else  // old version
    // distributing snp processors between ndir[curdir] slices
    for(int i=0;i<ndir[curdir];i++){ // cutting subboxes for othogonal direction
      // if dntot>0, cutting larger subboxes in order to make final subdomains equal
      vec_type dl=(dntot>0 ? l*c1 : l*c0);
      int cnp=(dntot>0 ? mnp+1 : mnp);
      if(cnp<=0)break;
      dntot--;
      //dl=incr(p1,p2,curdir,dl); // balanced domain decomposition: dl is recalculated taking into account volume load
      p2[curdir]=p1[curdir]+dl;
      Box subdiv(p1,p2);
      nspl[curdir]=0; // this direction is frozen for the next recursive call, where orthogonal direction will be split
      // recursive call; subboxes will be recorded to result, only if nproc>0
      int ddom=MakeBoxDomains(subdiv,(nproc>0 ? cnp: -cnp),nspl[0],nspl[1],nspl[2],result,weight_integral);
      if(nproc>0)
        result+=ddom;
      ndom+=ddom;
      p1[curdir]=p2[curdir];
    }
# endif
    return ndom;
  }
  // called with (0,0,0) split: put this box into the set
  if(nproc>0)
    *result=box;
  return 1;   
}

template<class contour_t>
vec_type GetContourFraction(const Region_3 *reg, const contour_t &cnt, int volume_step, 
vec_type a0, const Vector_3 &center, const Vector_3 &dir, VecContour<> *subcont, Vector_3 *normv){

  const vec_type eps=1e-5;
  Vector_3 insidep; // subcontour center
  vec_type a=reg->TestPtrContour(cnt,subcont,&insidep); // subcontour area
  bool isout=false, isin=false; // if there are control points outside / inside the region
  Vector_3 cout, cin; // point outside and inside the region

  if(normv){ // getting surface normal vector
    if (a/a0<eps) { // the whole contour is outside the region
      cout=center; // control point outside the region is center of the contour
      isout=true; //such point is found
      *normv=0.;
    }
    else if (a/a0>1-eps) { // the whole contour is inside the region
      cin=center; // control point inside the region is center of the contour
      isin=true; //such point is found
      *normv=0.;
    }
    else {
      cin=insidep; // control point inside the region is center of the subcontour
      cout=(center*a0-insidep*a)/(a0-a); // calculate control point outside the region
      isout=isin=true;
      reg->TestEdge(cout,cin,NULL,normv); // calculate contribution to normv
    }
  }
   
  if(volume_step){ // control volume
    vec_type v=cnt.GetControlVolume();
    vec_type h=v/a0; 
    Vector_3 shift=dir*h/2/volume_step;
    for(int i=-volume_step;i<=volume_step;i++){
      if(i==0)continue; // this is initial contour, it was already considered above
      VecContour<> vc(cnt.points_begin(), cnt.points_end(),1.,shift*i); // shifted contour
      vec_type ai=reg->TestPtrContour(vc,NULL,&insidep);
      a+=ai;
      if(normv){
        Vector_3 normvi;
        if (ai/a0<eps) {
          if (isin) {
            Vector_3 scenter=center+shift*i;
            reg->TestEdge(scenter,cin,NULL,&normvi);
          }
        }
        else if (ai/a0>1-eps) {
          if (isout) {
            Vector_3 scenter=center+shift*i;
            reg->TestEdge(cout,scenter,NULL,&normvi);
          }
        }
        else {
          if (isout) {
            reg->TestEdge(cout,insidep,NULL,&normvi);
          }
          else {
            Vector_3 tcenter=((center+shift*i)*a0-insidep*ai)/(a0-ai);
            reg->TestEdge(tcenter,cin,NULL,&normvi);
          }
        }
        *normv+=normvi;
      }
    }
    a/=(2*volume_step+1); //control volume fraction inside the region
  }
  return a;
}
