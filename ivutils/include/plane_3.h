/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.1 $
 *   $Date: 2012/12/06 02:24:05 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/plane_3.h,v 1.1 2012/12/06 02:24:05 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/plane_3.h,v $
$Revision: 1.1 $
$Author: lesha $
$Date: 2012/12/06 02:24:05 $
*/
/*s****************************************************************************
 * $Log: plane_3.h,v $
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.25  2012/10/02 23:11:29  lesha
 * SpaceRegion is renamed to Region_3
 *
 * Revision 1.24  2012/06/07 17:22:23  valuev
 * started adding scalable dipole
 *
 * Revision 1.23  2012/04/09 07:31:40  valuev
 * yet another fix for scatter_data
 *
 * Revision 1.22  2012/03/21 20:59:08  lesha
 * shift polyhedron is fixed
 *
 * Revision 1.21  2012/03/13 23:19:18  lesha
 * linear solvers are commented
 *
 * Revision 1.20  2010/03/04 16:12:35  valuev
 * corrected sphere radius in opt_module
 *
 * Revision 1.19  2010/02/06 23:36:04  valuev
 * updated documentation example
 *
 * Revision 1.18  2010/02/05 13:56:24  valuev
 * Tuned the Doxyfile for ivutils and added documentation for Plane_3
 *
 * Revision 1.17  2010/02/04 21:47:22  valuev
 * added doxyfile for ivutils
 *
 * Revision 1.16  2010/01/27 11:25:39  lesha
 * make SINGLE_PRECISION compilable
 *
 * Revision 1.15  2009/11/30 15:03:08  lesha
 * some comments are added
 *
 * Revision 1.14  2009/05/20 07:31:06  valuev
 * corrected normal direction in TestRay for plane_3
 *
 * Revision 1.13  2009/05/19 21:50:17  valuev
 * Added TestRay for plane
 *
 * Revision 1.12  2009/03/13 16:33:04  lesha
 * Ilya's changes concerning dipole
 *
 * Revision 1.11  2008/01/22 10:15:02  lesha
 * unnecessary include files are removed
 *
 * Revision 1.10  2008/01/06 03:28:43  lesha
 * #include "pencil.h" is removed
 *
 * Revision 1.9  2007/04/17 10:50:32  valuev
 * Added new contour path models: volume and tensor
 *
 * Revision 1.8  2007/04/03 14:13:39  lesha
 * set_coeffs is added
 *
 * Revision 1.7  2007/02/23 14:23:35  valuev
 * Added internal FFT to emSourceWave,
 * corrected transfers for PBC
 *
 * Revision 1.6  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.5  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/
# ifndef _PLANE_3_H
# define _PLANE_3_H

///     @file plane_3.h \brief Classes for planes in 3 dimensions 
/// \ru @file plane_3.h \brief Классы для работы с 3х-мерными плоскостями 

# include "vector_3.h"



/// \en Base plane class implementing plane
///     equation \f$\vec{n}\vec{x}+d=0\f$. 
/// \ru Класс, реализующий уравнение плоскости
///     \f$\vec{n}\vec{x}+d=0\f$.  
class Plane_3{
  Vector_3 nn; ///< \en plane vector \ru перпендикулярный вектор 
  /// \en constant \ru константа
  vec_type dd; 
public:
  

  /// \en Constructor with precomputed d. 
  /// \ru Конструктор с константой d, вычисленной заранее.
  /// \en \param[in] n plane vector (normalized automatically) 
  ///     \param[in] d constant 
  /// \ru \param[in] n вектор, перпендикулярный плоскости (нормализуется автоматически)
  ///     \param[in] d константа
  Plane_3(const Vector_3& n=Vector_3(1,1,1), vec_type d=0){
    nn=n;
    nn.normalize();
    dd=d;
  }


  /// \en Constructor from a given position. 
  /// \ru Конструктор плоскости, проходящей через заданную точку.
  /// \en \param[in] n plane normal vector (normalized automatically)  
  ///     \param[in] pos a point that belongs to the plane  
  /// \ru \param[in] n вектор, перпендикулярный плоскости (нормализуется автоматически)
  ///     \param[in] pos точка на плоскости
  Plane_3(const Vector_3&n, const Vector_3&pos){
    init(n,pos);
  }

  /// \en Initialization analogous to
  ///     Plane_3::Plane_3(const Vector_3& n, const Vector_3& pos). 
  /// \ru Инициализация, аналогичная 
  ///     Plane_3::Plane_3(const Vector_3& n, const Vector_3& pos).
  /// \en \retval <0 if \a n is \ref small_vectors "small"
  ///     \retval >0 otherwise
  /// \ru \retval <0 если \a n -- \ref small_vectors "малый вектор"
  ///     \retval >0 в противном случае
  int init(const Vector_3&n, const Vector_3&pos){
    nn=n;
    int res=1;
    if(nn.normalize()<VEC_ZERO)res=-1;
    set_pos(pos);
    return res;
  }

  /// \en Sets the position of a starting point not changing the normal vector.
  /// \ru Меняет плоскость, задавая положение точки, лежащей на ней, без изменения вектора нормали. 
  /// \en \param[in] newpos new plane position
  /// \ru \param[in] newpos новое пложение плоскости
  void set_pos(const Vector_3& newpos){
    dd=-(newpos*nn);
  }

  /// \en Shifts the plane by specified vector.
  void shift(const Vector_3& shift){
    dd-=(shift*nn);
  }

  /// \en Returns the distance \f$\vec{n}\vec{x}+d\f$ between 
  ///     a point \a x in space and the plane. 
  /// \ru Возвращает расстояние \f$\vec{n}\vec{x}+d\f$ от 
  ///     точки \a x в пространстве до плоскости.
  vec_type distance(const Vector_3 &x) const{
    return dd+x*nn;
  }

  /// \en Gets the plane vector \a n and the constant \a d.
  /// \ru Возварщает вектор \a n и константу \a d.
  void get_coeff(Vector_3& n, vec_type &d) const{
    n=nn;
    d=dd;
  }
  /// \en Sets the plane vector \a n and the constant \a d.
  ///     The plane vector is NOT normalized, for normalization use Plane_3::init()
  /// \ru Устанавливает вектор \a n и константу \a d.
  ///     Вектор \a n НЕ нормализуется, для нормализации см. Plane_3::init()
  void set_coeff(const Vector_3& n, const vec_type &d){
    nn=n;
    dd=d;
  }

  /// \en Gets the plane vector \a n.
  /// \ru Возварщает вектор \a n.
  Vector_3 get_normal() const {
    return nn;
  }

  /// \en Tests whether a segment \f$(\vec{v}_1,\vec{v}_2)\f$ crosses the plane.
  /// \en \param[in]  v1 segment starting position 
  ///     \param[in]  v2 segment end position 
  ///     \param[out] cross crossing position, returned if \a cross is not NULL
  /// \en \retval true if crossing point exists
  ///     \retval false otherwise
  /// \ru Проверяет, пересекает ли отрезок \f$(\vec{v}_1,\vec{v}_2)\f$ плоскость.
  /// \ru \param[in]  v1 начальная точка отрезка
  ///     \param[in]  v2 конечная точка отрезка
  ///     \param[out] точка пересечения с плоскостью, возвращается, если \a cross не равен NULL
  /// \ru \retval true если точка пересечения существует
  ///     \retval false в противном случае
  bool TestEdge(const Vector_3 &v1, const Vector_3 &v2, Vector_3 *cross=NULL) const {
    vec_type d1=distance(v1);
    vec_type d2=distance(v2);
    if(d1*d2>=0.)
      return false;
    if(cross){
      d1=fabs(d1);
      d2=fabs(d2);
      vec_type d=d2+d1;
      if(d<VEC_ZERO)
        return false;
      *cross=(v1*d2+v2*d1)/d;
    }
    return true;
  }

  /// \en Returns \f$t_0\f$ on the ray \f$\vec{p}=\vec{p}_1+\vec{k}t_0\f$ 
  ///     corresponding to the intersection with the plane in the positive direction.
  ///     Analogous to Region_3::TestRay() function.
  ///     \param[in] p1 ray starting point
  ///     \param[in] k ray direction
  ///     \param[out] surfp if not NULL, will contain the crossing point
  ///     \param[out] surfn if not NULL, will contain the plane normal \a n (normal at the crossing point)
  ///     \param[in]  epsilon intersections with \f$t_0<=\epsilon\f$ are ignored (-epsilon is reported as a solution)
  ///     \return intersection coordinate on the ray (\f$t_0\f$) 
  ///     \retval -epsilon means that intersection is found in epsilon vicinity of the starting point \a p1 
  ///                      (may be interpreted as \a p1 lying in the plane)
  ///     \retval <-epsilon means no intersection in the positive ray direction. Negative \f$t_0\f$
  ///                       value may however be used as a solution in the oposite direction
  ///     \retval -VEC_INFTY means that there is no intersection (ray is parallel to the plane with the accuracy of VEC_ZERO)
  /// \ru Возвращает значение \f$t_0\f$ на луче \f$\vec{p}=\vec{p}_1+\vec{k}t_0\f$, 
  ///     соответствующее пересечению луча с плоскостью в направленни возрастания \a t.
  ///     Аналогична функции Region_3::TestRay().
  ///     \param[in] p1 начальная точка луча
  ///     \param[in] k направление луча
  ///     \param[out] surfp если не NULL, содержит точку пересечения
  ///     \param[out] surfn если не  NULL, содержит нормаль \a n (нормаль в точке пересечения)
  ///     \param[in]  epsilon пересеченичя с \f$t_0<=\epsilon\f$ игнорируются (функция возвращает -epsilon)
  ///     \return координату \f$t_0\f$ точки пересечения вдоль луча
  ///     \retval -epsilon пересечение в эпсилон-окрестности начальной точки \a p1 
  ///                      (может интерпретироваться как принадлежность \a p1 плоскости)
  ///     \retval <-epsilon означает отсутствие пересечения в направлении возрастания координаты луча. 
  ///                       Отрицательное значение \f$t_0\f$
  ///                       может использоваться как решение, соответствующие обратному направлению луча 
  ///     \retval -VEC_INFTY означает отсутствие пересечения (луч параллелен плоскости с точностью VEC_ZERO)
  vec_type TestRay(const Vector_3 &p1, const Vector_3 &k, Vector_3 *surfp=NULL, Vector_3 *surfn=NULL, vec_type epsilon=0) const {
    Vector_3 vnorm=get_normal();
    vec_type dprod=k*vnorm;
    vec_type dist=distance(p1);
    if(fabs(dprod)<VEC_ZERO){ // all inside or all outside
      return -VEC_INFTY; // no solution
    }
    vec_type t=-dist/dprod;
    Vector_3 dv=t*k;
    if(epsilon!=0. && dv.norm2()<epsilon*epsilon)
      return -epsilon; // no solution
    if(surfp)
      *surfp=p1+dv;
    if(surfn){
      if(dist<0) // from outside     
        *surfn=-vnorm;
      else // from inside
        *surfn=vnorm;
    }  
    return t; // may be negative
  }

  //void get_syscoord(){}
};

/// \en Pointer to Plane_3. \ru Указатель на Plane_3.
typedef Plane_3 *Plane_3P; // не используется. убрать?

# endif
