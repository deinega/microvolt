/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *   $Revision: 1.4 $
 *   $Date: 2013/11/24 01:51:34 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/contour_t.cpp,v 1.4 2013/11/24 01:51:34 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/contour_t.cpp,v $
$Revision: 1.4 $
$Author: lesha $
$Date: 2013/11/24 01:51:34 $
*/
/*s****************************************************************************
 * $Log: contour_t.cpp,v $
 * Revision 1.4  2013/11/24 01:51:34  lesha
 * *** empty log message ***
 *
 * Revision 1.3  2013/02/22 11:01:26  valuev
 * thin surface interface
 *
 * Revision 1.2  2013/01/25 15:27:39  valuev
 * restructured contour: got rid of virtual functions and unnecessary storage
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2012/02/22 03:11:49  lesha
 * *** empty log message ***
 *
 * Revision 1.3  2012/02/22 02:52:52  lesha
 * *** empty log message ***
 *
 * Revision 1.2  2012/02/17 00:20:36  lesha
 * PtrContour is added
 *
 * Revision 1.1  2012/02/16 02:30:53  lesha
 * moving geometry installation to ivutils
 *
 * Revision 1.11  2012/02/16 02:03:34  lesha
 * *** empty log message ***
 *
 * Revision 1.10  2010/09/28 18:03:48  valuev
 * fixed Box:TestContour
 *
 * Revision 1.9  2010/06/08 22:43:23  lesha
 * make UNIX compilable
 *
 * Revision 1.8  2010/01/18 12:05:13  lesha
 * *** empty log message ***
 *
 * Revision 1.7  2010/01/01 20:39:43  lesha
 * Region, Box_N, Sphere_N, Basis_N, MonteCarloTestContour are added
 *
 * Revision 1.6  2009/12/08 21:13:30  lesha
 * *** empty log message ***
 *
 * Revision 1.5  2009/12/08 21:04:39  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2009/12/01 22:43:33  lesha
 * *** empty log message ***
 *
 * Revision 1.3  2009/12/01 18:38:26  lesha
 * contour is modified
 *
 * Revision 1.2  2009/06/21 10:40:32  lesha
 * region_2.h is added
 *
 * Revision 1.1  2009/01/30 14:02:16  valuev
 * restructured as a library
 *
*******************************************************************************/

/** @file contour_t.h @brief Template instantiations required for photonic emtl library from ivutils 
*/ 

#include "contour.hpp"
#include "region_3.h"

template Vector_3 GetBoundingBox(vector<Vector_3>::const_iterator, vector<Vector_3>::const_iterator,Vector_3 *, Vector_3 *);
template Vector_2 GetBoundingBox(vector<Vector_2>::const_iterator, vector<Vector_2>::const_iterator,Vector_2 *, Vector_2 *);
template Vector_2 GetBoundingBox(contour_projection_it<vector<Vector_3>::const_iterator>, contour_projection_it<vector<Vector_3>::const_iterator>,Vector_2 *, Vector_2 *);
template Vector_3 GetBoundingBox(Vector_3*, Vector_3*, Vector_3 *, Vector_3 *);
template Vector_2 GetBoundingBox(Vector_2*, Vector_2*, Vector_2 *, Vector_2 *);


template class Contour< Vector_3 * >;
template class Contour< Vector_3 *, array_store_t<Vector_3, 4> >;
template class Contour<vector<Vector_3>::const_iterator, vec_store_t<Vector_3> >;
template class Contour_N< Vector_3* >;
template class Contour_N< Vector_3*, Vector_3, array_store_t<Vector_3, 4> >;
template class Contour_N<vector<Vector_3>::const_iterator, Vector_3, vec_store_t<Vector_3> >;
//template class Contour_N<modified_value_it<indexed_it<std::vector<Vector_3>::const_iterator, std::vector<size_t>::const_iterator>, plus_t<Vector_3> > >;

template class Contour< Vector_2* >;
template class Contour< Vector_2 *, array_store_t<Vector_2, 4> >;
template class Contour<vector<Vector_2>::const_iterator, vec_store_t<Vector_2> >;
template class Contour_N< Vector_2* >;
template class Contour_N<vector<Vector_2>::const_iterator, Vector_2, vec_store_t<Vector_2> >;
template class Contour<contour_projection_it<Vector_3 * > >;
template class Contour<contour_projection_it<vector<Vector_3>::const_iterator> >;
template class Contour_N<contour_projection_it<vector<Vector_3>::const_iterator> >;

template int ProjectSimplex(const Plane_3 &,
aggregate_it<Plane_3 *, normplanes_it<generic_edge_it<Vector_3* > >, Plane_3>,
aggregate_it<Plane_3 *, normplanes_it<generic_edge_it<Vector_3* > >, Plane_3>,
VecContour<> &,int);

template int ProjectSimplex(const Plane_3 &,
aggregate_it<Plane_3 *, normplanes_it<generic_edge_it<vector<Vector_3>::const_iterator> >, Plane_3>,
aggregate_it<Plane_3 *, normplanes_it<generic_edge_it<vector<Vector_3>::const_iterator> >, Plane_3>,
VecContour<> &,int);

template int ProjectSimplex(const Plane_3 &,
aggregate_it<Box::plane_it, normplanes_it<generic_edge_it<Vector_3* > >, Plane_3>,
aggregate_it<Box::plane_it, normplanes_it<generic_edge_it<Vector_3* > >, Plane_3>,
VecContour<> &,int);

/*template int ProjectSimplex(const Plane_3 &,
aggregate_it<Box::plane_it, normplanes_it<generic_edge_it<vector<Vector_3>::const_iterator> >, Plane_3>,
aggregate_it<Box::plane_it, normplanes_it<generic_edge_it<vector<Vector_3>::const_iterator> >, Plane_3>,
VecContour<> &,int);*/

template int ProjectSimplex(const Plane_3 &,Box::plane_it,Box::plane_it,VecContour<> &,int);

//<aggregate_it<std::_Vector_iterator<Plane_3,std::allocator<Plane_3> >,
//normplanes_it<generic_edge_it<arr_point_it<Vector_3 > > >,Plane_3>,VecContour<3> >
/*template int ProjectSimplex(const Plane_3 &,
aggregate_it<std::_Vector_iterator<Plane_3,std::allocator<Plane_3> >,
normplanes_it<generic_edge_it<arr_point_it<Vector_3> > >,Plane_3>,
aggregate_it<std::_Vector_iterator<Plane_3,std::allocator<Plane_3> >,
normplanes_it<generic_edge_it<arr_point_it<Vector_3> > >,Plane_3>,
VecContour<> &,int);*/

template int GetFaces(vector<VecContour<> > &, Box::plane_it, Box::plane_it);
template int GetFaces(vector<VecContour<> > &, aggregate_it<Box::plane_it, Box::plane_it, Plane_3>, aggregate_it<Box::plane_it, Box::plane_it, Plane_3>);
template int GetFaces(vector<VecContour<> > &, Plane_3 *, Plane_3 *);
template int GetFaces(vector<VecContour<> > &, aggregate_it<Plane_3 *, Box::plane_it, Plane_3>, aggregate_it<Plane_3 *, Box::plane_it, Plane_3>);
template int GetFaces(vector<VecContour<> > &, vector<Plane_3>::iterator , vector<Plane_3>::iterator);
template int GetFaces(vector<VecContour<> > &, aggregate_it<vector<Plane_3>::iterator, Box::plane_it, Plane_3>, aggregate_it<vector<Plane_3>::iterator, Box::plane_it, Plane_3>);
