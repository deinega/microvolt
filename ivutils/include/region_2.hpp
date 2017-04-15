#include "region_2.h"

template<int N>
void Box_N<N>::init(const typename Region<N>::vector_t &sp1, const typename Region<N>::vector_t &sp2){
  p1=sp1;
  p2=sp2;
  for(int i=0;i<N;i++){
    if(p1[i]>p2[i]) { 
      vec_type tmp=p1[i];
      p1[i]=p2[i];
      p2[i]=tmp;
    }
    sz[i]=p2[i]-p1[i];
  }
}

template<int N>
int Sphere_N<N>::TestLine(const typename Region<N>::vector_t &p, const typename Region<N>::vector_t &dir, vec_type *frac, 
typename Region<N>::vector_t *surfp, typename Region<N>::vector_t *surfn, vec_type epsilon) const{

  vec_type r2=R*R;

  vector_t p0=p-center;
  vec_type p02=p0*p0;

  // |p0+k*dir|=R where k is some number
  // p0^2+k^2*dir^2+2k(p0*dir)=R^2
  vec_type a=dir*dir;
  vec_type b=2*(p0*dir);
  vec_type c=p02-r2;
  vec_type D=b*b-4*a*c;

  if(D<0)
    return 0; // there are no intersections
  if(D)
    D=sqrt(D);

  for(int i=0;i<2;i++){
    int sign=2*i-1;
    vec_type k=(-b+sign*D)/(2*a);

    vector_t touch=p0+k*dir;
    if(surfp)
      surfp[i]=center+touch;
    if(surfn){
      surfn[i]=touch;
      surfn[i].normalize();
    }
    frac[i]=k;
  }
  return 1+(D>0); // intersections number
}

template<int N>
typename Sphere_N<N>::vector_t Sphere_N<N>::GetBoundingBox(typename Region<N>::vector_t *v1, typename Region<N>::vector_t *v2) const{
  *v1=center-Vector_Nt<vec_type,N>(R);
  *v2=center+Vector_Nt<vec_type,N>(R);
  return center;
}

template <class contour_t>
vec_type Circle::TestContour(const contour_t &cnt, VecContour<2> *subcont, Vector_2 *subcenter) const{

  vec_type r2=R*R;
  // simple test

  int rotate(0);
  bool change_sign(false);

  Vector_2 a, b;
  vec_type a2, b2;

  typename contour_t::point_iterator it=cnt.points_begin(), e=cnt.points_end();
  a=*it-center;
  a2=a*a;
  bool ina = a2<r2 ? true : false;

  bool first(true);

  while (first) {
    ++it;
    if (!(it!=e)) {
      if (ina) {
        if (subcenter) 
          *subcenter=cnt.GetCenter();
        return cnt.Area(); // contour is inside circle
      }
      it=cnt.points_begin();
      first=false;
    }

    b=*it-center;
    b2=b*b;
    if (ina) {
      if (b2>r2) break; // there are points inside and outside circle      
    }
    else {
      if (b2<r2) break; // there are points inside and outside circle
      Vector_2 dr=b-a;
      vec_type dr2=dr*dr;
      vec_type mix=a*dr;
      vec_type k=-mix/dr2;
      vec_type d=a2+mix*k; // distance from center of circle to line containing a and b
      if (d<r2 && k<1 && k>0) {
//        if (!first)
//          first=false;
        first=true; // in order to passing 'if (!first) ...'
        break; // there are intersections
      }

      vec_type cur=a[0]*b[1]-a[1]*b[0];
      if (cur) {
        if (!rotate)
          rotate=cur>0 ? 1 : -1;
        else if (rotate!=(cur>0 ? 1 : -1)) 
          change_sign=true;
      }
      a=b;
      a2=b2;
    }
  }

  if (!first) {
    if (change_sign)
      return 0; // contour is outside circle
    else {
      if (subcenter) 
        *subcenter=center;
      return M_PI*r2; // circle is incide contour
    }
  }

  // contour intersects circle; more complicated analysis

  if (subcenter) 
    *subcenter=Vector_2();

  vec_type angle(0); // total sectors angle
  vec_type area(0); // total triangles area
    
  it=cnt.points_begin();
  a=*it-center;
  a2=a*a;
  ina = a2<r2 ? true : false;

  first=true;
    
  while (first) {
    ++it;
    if (!(it!=e)) {
      it=cnt.points_begin();
      first=false;
    }

    b=*it-center;
    b2=b*b;
    bool inb = b2<r2 ? true : false;

    do {
      Vector_2 inter;
      vec_type inter2;

      if (ina && inb) {
        inter=b; // a and b are inside circle
      }
      else {
        Vector_2 dr=b-a;
        vec_type a_=dr*dr;
        vec_type b_=2*(a*dr);
        vec_type c_=a2-r2;
        vec_type D=b_*b_-4*a_*c_;

        if (D<=0) {
          inter=b; // a and b are outside circle and connecting line does not intersect it (or touch it)
        }
        else {
          vec_type k=(-b_+sqrt(D)*(ina ? 1 : -1))/(2*a_); // if a is inside circle I find the next intersection toward to b; else I find the nearest intersection
          if (k<0 || k>=1)
            inter=b; // a and b are outside circle and connecting line does not intersect it or a is inside and b touches circle
          else
            inter=a+k*dr;
        }
      }
      inter2=inter*inter;        

      if (ina) {
        vec_type darea=(a[0]*inter[1]-a[1]*inter[0])/2; // triangle
        area+=darea;
        if (subcenter)
          *subcenter+=darea*(a+inter)/3; // medians intersection
      }
      else {
        vec_type sign=a[0]*inter[1]-a[1]*inter[0];
        if (sign) {
          sign = sign>0 ? 1 : -1;
          vec_type cos_dangle=a*inter/(sqrt(a2*inter2));
          vec_type dangle;
          if (cos_dangle>=1) {
            cos_dangle=1;
            dangle=0;
          } // posible error if a*inter/(sqrt(a2*inter2)) = 1 + small number
          else if (cos_dangle<=-1) {
            cos_dangle=-1;
            dangle=M_PI*sign;
          } // posible error if a*inter/(sqrt(a2*inter2)) = -1 - small number
          else
            dangle=acos(cos_dangle)*sign; // sector
          angle+=dangle;
          if (subcenter) {
            Vector_2 ax=a;
            Vector_2 ay(-a[1], a[0]);
            ax.normalize();
            ay.normalize();
            vec_type sin_dangle=sign*sqrt(1-cos_dangle*cos_dangle); // -pi/2 < dangle < pi/2
            Vector_2 seccent=(R*r2/3)*(ax*sin_dangle+ay*(1-cos_dangle)); // center*ds
            *subcenter+=seccent;
          }
        }
      }

      a=inter;
      a2=inter2;
      if (inter==b) {
        ina=inb;
        break;
      }
      ina=!ina;

    } while (true);
  }
  area+=angle*r2/2;

  if (subcenter) {
    *subcenter/=area;
    *subcenter+=center;
  }

  return fabs(area);
}

template<class reg_tt>
int StretchedRegion<reg_tt>::TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
vector_t *surfp, vector_t *surfn, vec_type epsilon) const{
  Vector_3 p1=p-shift;
  p1=basis_inv(p1);
  Vector_3 dir1=basis_inv(dir);
  Vector_3 surfp1[2];
  int res=reg->TestLine(p1,dir1,frac,surfp1,surfn,epsilon);
  vec_type dir_n=dir.norm2();
  vec_type dir_n1=dir1.norm2(); // for test
  for(int i=0;i<res;i++){
    if(surfp)
      surfp[i]=basis(surfp1[i])+shift;
    if(surfn){
      surfn[i]=basis(surfn[i]);
      surfn[i].normalize();
    }
    vec_type norm=(surfp[i]-p)*dir;
    vec_type norm1=(surfp1[i]-p1)*dir1; // for test
    frac[i]=norm/dir_n;
  }
  return res;
}
