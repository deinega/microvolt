#include <stdio.h>
#include <memory>
#include "gnudump.h"
#include "cpp11features.h"

#define INSTANTIATE_GetRegionDumper(type) \
	template RegDumper<type::dimension> *GetRegionDumper< type >(const type *);

template<int N>
int Dump(FILE *f, const pair<Vector_Nt<vec_type,N>, Vector_Nt<vec_type,N> > &p){
  for(int i=0;i<N;i++)fprintf(f,"%g%c",p.first[i],9);
  for(int i=0;i<N-1;i++)fprintf(f,"%g%c",p.second[i],9);
  fprintf(f,"%g\n",p.second[N-1]);
  return 1;
}


template<class contour_t>
int DumpContour(FILE *f, const contour_t &cnt, int first, bool deltas){
  typedef typename contour_t::point_iterator point_it;
  if(first){
    int i=1;
    fprintf(f,"#%d-x %d-y",i,i+1);
    i+=2;
    if(contour_t::vector_t::dimension==3)
      fprintf(f," %d-z",i++);
    if(deltas){
      fprintf(f," %d-delta_x %d-delta_y",i,i+1);
      i+=2;
      if(contour_t::vector_t::dimension==3)
        fprintf(f," %d-delta_z",i++);
    }
    fprintf(f,"\n\n");
  }
  //if(cnt.GetNPoints()){
    bool have_points=false;
    typename contour_t::vector_t pprev;
    for(point_it pit=cnt.points_begin(), pe=cnt.points_end();pit!=pe;++pit){
      typename contour_t::vector_t p=*pit;
      if(have_points && deltas)
        Dump(f,p-pprev, true); // duplicate prev to draw with lines
      Dump(f,p, !deltas);
      pprev = p;
      have_points=true;
    }
    if(have_points){
      typename contour_t::vector_t p=*(cnt.points_begin());
      if(deltas)
        Dump(f,p-pprev, true); // duplicate prev to draw with lines
      Dump(f,p,!deltas); // loop the point
      if(deltas)
        Dump(f,Vector_3(),true); // last delta is zero
      fprintf(f,"\n\n");
    }
  //}
  fflush(f);
  return 1;
}

template<class contour_t>
int DumpContours(FILE *f, const vector<contour_t> &cnt, int first){
  for(size_t i=0;i<cnt.size();i++)
    DumpContour(f,cnt[i],i==0 && first);
  return 1;
}

template<class contour_t>
int DumpContoursVTK(FILE *f, const vector<contour_t> &cnt, int first){
    
  int pcount=0; // total number of points
  for(size_t i=0;i<cnt.size();i++)
    pcount+=cnt[i].GetNPoints();

  string tailer="  </PolyData>\n</VTKFile>\n";
  if(!first){
    long tailer_sz=(long)tailer.size();
    fseek(f,0,SEEK_END);
    long pos=ftell(f);
    fseek(f,pos-tailer_sz,SEEK_SET); // overwrite the tailer if this is not the first write
  }
  else{ // the first write to this file
    fprintf(f,"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(f,"  <PolyData>\n");
  }
  fprintf(f,"    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
                 "NumberOfStrips=\"0\" NumberOfPolys=\"%d\">\n",pcount,(int)cnt.size());
  fprintf(f,"      <Points>\n");
  fprintf(f,"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(f,"          ");
  // writing coordinates
  int ipcount=0;
  for(size_t i=0;i<cnt.size();i++){
    int np=cnt[i].GetNPoints();
    for(int j=0;j<np;j++){
      fprintf(f,"%g %g %g ",cnt[i].GetPoints()[j][0],cnt[i].GetPoints()[j][1],cnt[i].GetPoints()[j][2]);
      ipcount++;
      if(!(ipcount%9) && pcount!=ipcount)
        fprintf(f,"\n          ");
    }
  }
  fprintf(f,"\n");
  fprintf(f,"        </DataArray>\n");
  fprintf(f,"      </Points>\n");
  fprintf(f,"      <Polys>\n");
  fprintf(f,"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  fprintf(f,"          ");
  // writing the order
  ipcount=0;
  for(size_t i=0;i<cnt.size();i++){
    int np=cnt[i].GetNPoints();
    for(int j=0;j<np;j++){
      fprintf(f,"%d ",ipcount);
      ipcount++;
      if(!(ipcount%9) && pcount!=ipcount)
        fprintf(f,"\n          ");
    }
  }
  fprintf(f,"\n");
  fprintf(f,"        </DataArray>\n");
  fprintf(f,"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  fprintf(f,"          ");
  // writing offsets
  ipcount=0;
  int offset=0;
  for(size_t i=0;i<cnt.size();i++){
    offset+=cnt[i].GetNPoints();
    fprintf(f,"%d ",offset);
    ipcount++;
    if(!(ipcount%20) && (int)cnt.size()!=ipcount){
      fprintf(f,"\n          ");
    }
  }
  fprintf(f,"\n");
  fprintf(f,"        </DataArray>\n");
  fprintf(f,"      </Polys>\n");
  fprintf(f,"    </Piece>\n");
  fprintf(f,"%s", tailer.c_str()); // closing the XML record anyway
  return 1;
}

template<class vset_t>
int DumpVectorSet(FILE *f, const vset_t *vs, const VecTransform *trans, int first){
  if(first)
    fprintf(f,"#1-x 2-y 3-z 4-vec_num\n");
  typename vset_t::const_iterator it=vs->begin(), e=vs->end();
  for(int k=0; it!=e; ++it, k++){
    Vector_3 v=trans ? (*trans)(*it) : *it;
    fprintf(f,"%g %g %g %d\n",v[0],v[1],v[2],k);
  }
  return 1;
}

template<class base_t>
interp_t<3> InterpolateRegion(const Cylinder<base_t> *B) {
  VecContour<2> cnt=InterpolateRegion(B->get_base()); // contour for base
  int sz=cnt.GetNPoints(); // number of vertices 
  if(sz<3)return interp_t<3>((Polyhedron_3 *)NULL,0);
  Plane_3 *planes = new Plane_3[sz+1];
  // consider all vertices of vertex
  VecContour<2>::point_iterator it=cnt.points_begin(), e=cnt.points_end();
  Vector_2 v2=*it; // first vertex of the base contour (2D coordinate in the base plane)
  Vector_3 b0=B->Vector_2to3(v2); // find 3D coordinates
  Vector_3 b1=b0; // remember this vertex
  ++it; // switch to second vertex
  for (int i=0; it!=e; ++it, i++) { // consider all other vertices
    Vector_2 v2=*it;
    Vector_3 b2=B->Vector_2to3(v2);
    Vector_3 n=-(b2-b1)%(B->n); // get normal for side plane
    planes[i].init(n,b1); // initialize side plane
    b1=b2;
  }
  Vector_3 n=-(b0-b1)%(B->n);
  planes[sz-1].init(n,b1); // last side plane
  // return polyhedron defined by planes
  return new Polyhedron_3(planes, planes+sz);
}

template<class base_t>
interp_t<3> InterpolateRegion(const Cone<base_t> *B) {
  Vector_3 top=B->origin+B->L*B->n;
  VecContour<2> cnt=InterpolateRegion(B->get_base()); // contour for base
  int sz=cnt.GetNPoints(); // number of vertices 
  if(sz<3)return interp_t<3>((Polyhedron_3 *)NULL,0);
  Plane_3 *planes = new Plane_3[sz+1];
  VecContour<2>::point_iterator it=cnt.points_begin(), e=cnt.points_end();
  Vector_2 v2=*it; // first vertex of the base contour (2D coordinate in the base plane)
  Vector_3 b0=B->Vector_2to3(v2); // find 3D coordinates
  Vector_3 b1=b0; // remember this vertex
  ++it; // switch to second vertex
  for (int i=0; it!=e; ++it, i++) { // consider all other vertices
    Vector_2 v2=*it;
    Vector_3 b2=B->Vector_2to3(v2);
    Vector_3 n=(b2-b1)%(b2-top); // get normal for side plane
    planes[i].init(n,b1); // initialize side plane
    b1=b2;
  }
  Vector_3 n=(b0-b1)%(b0-top);
  planes[sz-1].init(n,b1); // last side plane
  // return polyhedron defined by planes
  return new Polyhedron_3(planes, planes+sz);
}
/*
template<class reg_t>
Polyhedron_3 *InterpolateRegion(const StretchedRegion<reg_t> *B,int *poly_num) {
  Polyhedron_3 *poly=InterpolateRegion(B->reg.ptr());
  for(Polyhedron_3::plane_it it=poly->planes_begin(), e=poly->planes_end();it!=e;++it){
    Vector_3 n;
    vec_type d;
    it->get_coeff(n,d);
    Vector_3 pos=-d*n;
    it->init(B->basis(n),B->origin+B->basis(pos));
  }
  return poly;
}
*/
template<class reg_t>
interp_t<3> InterpolateRegion(const ConfinedRegion<reg_t> *B) {
  UNIQUE_OR_AUTO(RegDumper<3>) dmp(B->reg->CreateDumper());
  if(!dmp.get())
    return NULL;
  interp_t<3> plreg=dmp->InterpolateRegion();

  Polyhedron_3 *plc=B->poly.ptr();
  
  int psz=plreg.num;
  if(psz<1)
    return interp_t<3>((Polyhedron_3 *)NULL,0);

  Polyhedron_3 *pl = psz==1 ? new Polyhedron_3 : new Polyhedron_3[psz];

  for(int i=0;i<psz;i++){
    int sz=0;
    Polyhedron_3::plane_it it=(plreg.poly+i)->planes_begin(), e=(plreg.poly+i)->planes_end();
    for (; it!=e; ++it)
      sz++;
    it=plc->planes_begin(), e=plc->planes_end();
    for (; it!=e; ++it)
      sz++;

    Plane_3 *planes = new Plane_3[sz+1];
    sz=0;

    it=(plreg.poly+i)->planes_begin(), e=(plreg.poly+i)->planes_end();
    for (; it!=e; ++it)
      planes[sz++]=*it;
    it=plc->planes_begin(), e=plc->planes_end();
    for (; it!=e; ++it)
      planes[sz++]=*it;
    (pl+i)->init(planes, planes+sz);
  }
  return pl;
}

template<class reg_t>
interp_t<3> InterpolateRegion(const Inverse<reg_t> *B) {
  UNIQUE_OR_AUTO(RegDumper<3>) dmp(B->reg->CreateDumper());
  if(!dmp.get())
    return NULL;
  return dmp->InterpolateRegion();
}

template<class reg_t, class limiting_t>
bool Dump(const StretchedRegion<reg_t> *reg, vector<VecContour<> > &sides, limiting_t *lim){
  UNIQUE_OR_AUTO(RegDumper<3>) dmp(reg->reg->CreateDumper());
  if(!dmp.get())
    return false;
  dmp->Dump(sides, lim);

  for(size_t i=0;i<sides.size();i++){
    for(size_t j=0;j<sides[i].GetPoints().size();j++){
      sides[i].GetPoints()[j]=reg->shift+reg->basis(sides[i].GetPoints()[j]);
    }
  }
  return true;
}
