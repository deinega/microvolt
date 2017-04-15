#ifndef GNUDUMP_H
#define GNUDUMP_H

/// \en @file gnudump.h
/// Classes used to record information about vectors / contours / bodies 
/// to text files in gnuplot format (for visualization)
/// \ru @file gnudump.h
///  лассы, позвол€ющие записывать информацию о векторах / контурах / телах в файлы,
/// которые потом можно визуализировать с помощью программы gnuplot

#include <string>
#include "vector_set.h"
#include "region_3.h"

///\en records vector coordinates in text file
///\ru записывает координаты вектора в файл
template<int N>
int Dump(FILE *f, const Vector_Nt<vec_type,N> &p, bool eol=true){
  for(int i=0;i<N-1;i++)fprintf(f,"%g%c",p[i],9);
  fprintf(f,"%g",p[N-1]);
  if(eol)
    fprintf(f,"\n");
  else
    fprintf(f,"%c",9);

  return 1;
}

///\en records coordinates of two vectors in text file
/// used to plot vectors by gnuplot using option 'w vec'
///\ru записывает в одну строчку координаты пары векторов
/// »спользуетс€ тогда, когда рисуем вектора в гнуплоте (точка начала вектора и сам вектор) с помощью опции 'w vec'
template<int N>
int Dump(FILE *f, const pair<Vector_Nt<vec_type,N>, Vector_Nt<vec_type,N> > &p);

///\en records sequence contour vertices to text file
/// first vertex is recorded twice as it is required by gnuplot format
/// if first<>0 records line with text file column description
/// \a deltas means print delta column to draw with vectors
///\ru записывает в файл вершины контура 
/// координаты начальной вершины повтор€ютс€ в заключительной строчке, чтобы гнуплоту было €сно, что это замкнутный контур
/// если first<>0 в начале записываетс€ дополнительна€ строчка с указанием наименований столбцов
template<class contour_t>
int DumpContour(FILE *f, const contour_t &cnt, int first=0, bool deltas = false);

///\en calls DumpContour for each contour
/// if first<>0 records line with text file column description
///\ru по очереди вызывает DumpContour дл€ вектора контуров
template<class contour_t>
int DumpContours(FILE *f, const vector<contour_t> &cnt, int first=0);

template<class contour_t>
int DumpContoursVTK(FILE *f, const vector<contour_t> &cnt, int first=0);

///\en records sequence of vectors
template<class vset_t>
int DumpVectorSet(FILE *f, const vset_t *vset, const VecTransform *trans=NULL, int first=0);

///\en records sequence of unit vectors difined by their origins and directions k
/// Box size is used to normalize direction k (if B!=NULL)
/// in order to keep vector inside the box
/// if first<>0 records line with text file column description
int DumpDirections(FILE *f, int n, const Vector_3 *origin, const Vector_3 *k, const Box *B=NULL, int first=0);

///\en interpolation of num-dimensional body
/// in our model 2D bodies are interpolated by VecContour<2> 
/// and 3D bodies are interpolated by array of Polyhedron_3
template<int num>
struct interp_t{};

template<>
struct interp_t<2>: public VecContour<2>{
  interp_t<2>(){}
  interp_t<2>(const VecContour<2> &other):VecContour<2>(other){}
};

template<>
struct interp_t<3>{
  Polyhedron_3 *poly;
  int num;

  interp_t<3>(Polyhedron_3 *poly_=NULL, int num_=0):poly(poly_),num(poly_ ? (num_ ? num_ : 1) : 0){}

  ~interp_t<3>(){
    if(num==1)delete poly;
    else if(num>1)delete[]poly;
  }
};

///\en interpolates region for dumping
/// by default region is interpolated by nothing
template<class reg_t>
interp_t<reg_t::dimension> InterpolateRegion(const reg_t *reg){
  return interp_t<reg_t::dimension>();
}

interp_t<2> InterpolateRegion(const Polygon_2 *poly);
interp_t<2> InterpolateRegion(const Circle *base);

interp_t<3> InterpolateRegion(const Box *B);
interp_t<3> InterpolateRegion(const Polyhedron<Plane_3 *> *B);
interp_t<3> InterpolateRegion(const Sphere *B);
template<class base_t>
interp_t<3> InterpolateRegion(const Cylinder<base_t> *B);
template<class base_t>
interp_t<3> InterpolateRegion(const Cone<base_t> *B);
template<class reg_t>
interp_t<3> InterpolateRegion(const ConfinedRegion<reg_t> *B);
template<class reg_t>
interp_t<3> InterpolateRegion(const Inverse<reg_t> *B);

///\en records sides of given polyhedron
template<class poly_t, class limiting_t>
bool DumpPolyhedron(const poly_t *reg, vector<VecContour<> > &sides, limiting_t *lim){
  if(lim){
    typedef aggregate_it<typename poly_t::plane_it, typename limiting_t::plane_it> iter_t;
    iter_t beg(reg->planes_begin(),reg->planes_end(),lim->planes_begin());
    iter_t end(reg->planes_end(),reg->planes_end(),lim->planes_end());
    return GetFaces(sides,beg,end)>0;
  }
  else{
    return GetFaces(sides,reg->planes_begin(),reg->planes_end())>0;
  }
}

///\en interpolate 3d body as a number of polyhedra and records their sides
/// these sides can be plotted by gnuplot afterwards
/// if lim!=NULL consider intersection between body and lim
///\ru аппроксимирует тело многогранником и записывае его грани в sides
/// при не равном нулю указателе lim, рассматриваетс€ часть тела, помещающа€с€ в соответствующий объект lim
template<class reg_t, class limiting_t>
bool Dump(const reg_t *reg, vector<VecContour<> > &sides, limiting_t *lim){
  interp_t<3> poly=InterpolateRegion(reg);
  for(int i=0;i<poly.num;i++)
    DumpPolyhedron(poly.poly+i,sides,lim);
  return true;
}

///\en interpolate 2d body as a contour that can be plotted by gnuplot afterwards
/// if lim!=NULL consider intersection between body and lim
template<class reg_t, class limiting_t>
bool Dump(const reg_t *reg, vector<VecContour<2> > &sides, limiting_t *lim){
  sides.push_back(InterpolateRegion(reg));
  return true;
}

///\en specifications of Dump for some specific bodies

template<class plane_it, class limiting_t>
bool Dump(const Polyhedron<plane_it> *reg, vector<VecContour<> > &sides, limiting_t *lim){
  return DumpPolyhedron(reg,sides,lim);
}

template<class limiting_t>
bool Dump(const Polyhedron_3 *reg, vector<VecContour<> > &sides, limiting_t *lim){
  return DumpPolyhedron(reg,sides,lim);
}

template<class reg_t, class limiting_t>
bool Dump(const StretchedRegion<reg_t> *reg, vector<VecContour<> > &sides, limiting_t *lim);


///\en abstract class with Dump and InterpolateRegion functionality
template<int N>
class RegDumper{
public:
  ///\en see documentation to global function Dump
  ///\ru см. описание к глобальной функции Dump
  virtual bool Dump(vector<VecContour<N> > &sides, Box *lim=NULL) const=0;

  ///\en see documentation to global function InterpolateRegion
  virtual interp_t<N> InterpolateRegion() const=0;
};

///\en implementation of RegDumper for particular body
template<class reg_t>
class RegionDumper: public RegDumper<reg_t::dimension>{
  const reg_t *reg;
public:
  RegionDumper(const reg_t *reg_):reg(reg_){}
  ///\en implementation of Dump for some particular body type
  ///\ru реализаци€ функции Dump дл€ конкретного типа тела
  bool Dump(vector<VecContour<reg_t::vector_t::dimension> > &sides, Box *lim=NULL) const{
    return ::Dump/*<reg_t, Box>*/((const reg_t *)reg,sides,lim); // expicit template arguments are needed here for gcc 4.4
  }

  interp_t<reg_t::dimension> InterpolateRegion() const{
    return ::InterpolateRegion(reg);
  }
};

template<class reg_t>
RegDumper<reg_t::dimension> *GetRegionDumper(const reg_t *reg){
  return new RegionDumper<reg_t>(reg);
}

///\en test class to dump bodies
///\ru тестовый рисовальщик
class GnuDumper{
  mngptr<Box> lim;
public:
  GnuDumper(mngarg<Box> lim_=NULL):lim(lim_){}
  template<class contour_t>
  int DumpContour(const contour_t &cnt,const char *fname){
    FILE *f=fopen(fname,"w");
    ::DumpContour(f,cnt,0);
    fclose(f);
    return 1;
  }
  template<class reg_t>
  int Dump(const reg_t *reg,const char *fname, Box *box=NULL){
    RegDumper<reg_t::dimension> *dmp=reg->CreateDumper();
    if(!dmp)
      return -1;
    vector<VecContour<reg_t::dimension> > sides;
    dmp->Dump(sides,box ? box : lim.ptr());
    if(!sides.size())
      return 0;
    FILE *f=fopen(fname,"w");
    DumpContours(f,sides,0);
    fclose(f);
    return 1;
  }
};

#endif
