/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2013        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.11 $
 *   $Date: 2013/04/22 12:44:44 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/facetedsurf.h,v 1.11 2013/04/22 12:44:44 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/facetedsurf.h,v $
$Revision: 1.11 $
$Author: valuev $
$Date: 2013/04/22 12:44:44 $
*/
/*s****************************************************************************
 * $Log: facetedsurf.h,v $
 * Revision 1.11  2013/04/22 12:44:44  valuev
 * reverted to include contour hpp
 *
 * Revision 1.10  2013/04/20 06:20:27  lesha
 * *** empty log message ***
 *
 * Revision 1.9  2013/02/22 11:01:26  valuev
 * thin surface interface
 *
 * Revision 1.8  2013/02/13 19:14:31  valuev
 * added cube facets
 *
 * Revision 1.7  2013/02/12 10:10:40  valuev
 * fixed facet region crossing
 *
 * Revision 1.6  2013/02/11 17:00:26  valuev
 * fixes
 *
 * Revision 1.5  2013/02/11 07:38:33  valuev
 * faceted region crossing
 *
 * Revision 1.4  2013/02/02 19:56:29  valuev
 * fixed new vector product
 *
 * Revision 1.3  2013/02/01 18:01:12  valuev
 * added vector volume (needs testing)
 *
 * Revision 1.2  2013/01/29 16:41:32  valuev
 * added faceted surface and rotation quaternions
 *
 * Revision 1.1  2013/01/28 16:59:46  valuev
 * added faceted surface generator class
 *
 * PtrContour is added
 *
*******************************************************************************/
#ifndef FACETEDSURF_H
#define FACETEDSURF_H

/// \en @file facetedsurf.h \brief Class using contour.h to define a surface consisting of multiple facets (contours).

#include <vector>
#include "contour.h"
#include "utiltl.h"
#include "region_3.h"
#include "contour.hpp" // for contour implementations

class FacetedSurfaceData {
public:
  std::vector<Vector_3> points;
  std::vector<size_t> index;
  size_t num_contours;


  FacetedSurfaceData():num_contours(0){}

   ///\en Adds a point to the point array. The points 
  ///    are then referred to by their indices in the array and may
  ///    be used as contour vertices for facets.
  size_t AddPoint(const Vector_3& point){
    points.push_back(point);
    return points.size();
  }
  ///\en Adds a contour of \a num_vertices points. The points are given by
  ///    a sequence of size_t indices started and iterated by \a index_beg. The indices refer to
  ///    the point array.
  template <class forward_it>
  size_t AddContour(size_t num_vertices, forward_it index_beg){
    index.push_back(num_vertices);
    for(size_t i=0;i<num_vertices; i++)
      index.push_back((size_t)*index_beg++);
    num_contours++;
    return num_contours; 
  }

  ///\en Iterates over edge indicies for all facets
  class edge_index_it {
    friend class FacetedSurfaceData;
    const FacetedSurfaceData *parent;
    size_t pos;
    size_t iedge;
  public:
    typedef forward_iterator_tag iterator_category;
    typedef pair<size_t,size_t> value_type;
    

    edge_index_it(const FacetedSurfaceData *parent_ = NULL, size_t pos_=0):parent(parent_), pos(pos_), iedge(0){}    

    ///\en edge is returned in pair: first is the source vertex, second is the destination(arrow) vertex
    value_type operator*() const {
      if(iedge<parent->index[pos]-1)
        return make_pair(parent->index[pos+1+iedge],parent->index[pos+1+iedge+1]);
      else // last, first
        return make_pair(parent->index[pos+1+iedge],parent->index[pos+1]);
    }
    
    edge_index_it &operator++(){
      iedge++;
      if(iedge>=parent->index[pos]){ // next contour
        iedge=0;
        if(pos<parent->index.size())
          pos+=parent->index[pos]+1;
      }
      return *this;
    }

    edge_index_it operator++(int){ // postfix
      edge_index_it tmp=*this;
      ++*this;
      return tmp;
    }

   
    bool operator==(const edge_index_it &other) const {
      return (pos==other.pos && iedge==other.iedge);
    }
    bool operator!=(const edge_index_it &other) const {
      return !(*this==other);
    }
  };

  edge_index_it edge_index_begin() const {
    return edge_index_it(this,0);
  }

  edge_index_it edge_index_end() const {
    return edge_index_it(this,index.size());
  }
  
  ///\en Gets full size (number of edges) of the contour which iterated edge belongs to
  size_t get_cur_contour_size(const edge_index_it &edgeit) const {
    return index[edgeit.pos];
  }


};

template< class transform_t = noop_t<Vector_3> >
class FacetedSurfaceGenerator {
  const FacetedSurfaceData *data;
  transform_t tr;
public:
  typedef std::vector<Vector_3>::const_iterator src_point_it;
//  typedef modified_value_it<indexed_it<src_point_it, index_it>, transform_t > result_point_it;
  typedef std::vector<size_t>::const_iterator index_it;
  class contour_t: public Contour_N<modified_value_it<indexed_it<src_point_it, index_it>, transform_t > > {
  public:
    typedef indexed_it<src_point_it, index_it> src_indexed_it;
    typedef Contour_N<modified_value_it<indexed_it<src_point_it, index_it>, transform_t > > base_t;
    typedef typename base_t::vector_t vector_t;
    typedef typename base_t::point_iterator point_iterator;
    static const int dimension = base_t::dimension;
    typedef typename base_t::edge_iterator edge_iterator;
    typedef typename base_t::edge_type edge_type;

    
    contour_t(const FacetedSurfaceGenerator &parent, size_t start_ind, size_t num): 
      base_t(
             point_iterator(src_indexed_it(parent.data->points.begin(),parent.data->index.begin()+start_ind), parent.tr),
             point_iterator(src_indexed_it(parent.data->points.begin(),parent.data->index.begin()+start_ind+num), parent.tr)
             ) {};

    //int GetNPoints() const {
    //  return (int);
    //}
  };

  class contour_it {
    const FacetedSurfaceGenerator *parent;
    size_t pos;
  public:
    typedef forward_iterator_tag iterator_category;
    typedef contour_t value_type;
    

    contour_it(const FacetedSurfaceGenerator *parent_ = NULL, size_t pos_=0):parent(parent_), pos(pos_){}    

    
    value_type operator*() const {
      return contour_t(*parent,pos+1,parent->data->index[pos]);
    }
    
    contour_it &operator++(){
      if(pos<parent->data->index.size())
        pos+=parent->data->index[pos]+1;
      return *this;
    }

    contour_it operator++(int){ // postfix
      contour_it tmp=*this;
      ++*this;
      return tmp;
    }

   
    bool operator==(const contour_it &other) const {
      return (pos==other.pos);
    }
    bool operator!=(const contour_it &other) const {
      return !(*this==other);
    }
   
  };

  typedef modified_value_it<src_point_it, transform_t>  point_it; ///<\en iterator of modified data points
  typedef FacetedSurfaceData::edge_index_it edge_index_it; ///<\en iterator of edges as index pairs

  FacetedSurfaceGenerator(const FacetedSurfaceData *data_, const transform_t &tr_=transform_t()):data(data_), tr(tr_){}

  size_t contours_size() const {
    return data->num_contours;
  }

  contour_it contours_begin() const {
    return contour_it(this,0);
  }
  contour_it contours_end() const {
    return contour_it(this,data->index.size());
  }


  size_t points_size() const {
    return data->points.size();
  }

  point_it points_begin() const {
    return point_it(data->points.begin(), tr);
  }
  point_it points_end() const {
    return point_it(data->points.end(), tr);
  }

  edge_index_it edge_index_begin() const {
    return data->edge_index_begin();
  }

  edge_index_it edge_index_end() const {
    return data->edge_index_end();
  }

  ///\en Gets full size (number of edges) of the contour which iterated edge belongs to
  size_t get_cur_contour_size(const edge_index_it &edgeit) const {
    return data->get_cur_contour_size(edgeit);
  }

  ///\en Calculates the middle point
  Vector_3 GetCenter() const {
   point_it p_beg=points_begin(), p_end=points_end();
    size_t n=0;
    Vector_3 sum;
    for(;p_beg!=p_end;++p_beg, ++n)
      sum+=*p_beg;
    if(n)
      sum/=n;
    return sum;
  }

  ///\en Calculates the volume. The volume is positive for anti-clockwise 
  ///    (right hand rule) facet contours with normals pointing INSIDE the volume.
  vec_type Volume() const{
    Vector_3 cnt=GetCenter();
    contour_it c_beg=contours_begin(), c_end=contours_end();
    vec_type sum=0.;
    for(;c_beg!=c_end;++c_beg){
      Vector_3 base;
      Vector_3 varea=(*c_beg).VecArea(&base);
      sum+=(cnt-base)*varea;
    }
    return sum;
  }

  Vector_3 GetPoint(size_t ind) const {
    Vector_3 v(data->points[ind]);
    return tr(v);
  }

  ///\en Creates (new) data for the faceted surface part belonging to given Region_3. Optionally in \a points_inside the reg.TestPoint()
  ///    results for each of the faceted surface points (in their order in data->points) may be supplied if known in advance.
  ///    The region must be convex. To detect intersection the region must contain one or more vertices of faceted surface, otherwise
  ///    void intersection is always returned.
  FacetedSurfaceData * CreateRegionSection(const Region_3 &reg, const vector<bool> *points_inside = NULL) const ; 


};

///\en Creates (new) data whcih are afacets aligned on a 2D (x,y) grid with total size Lx x Ly and number of segments nx x ny.
///    If the bit flag loop is set: 0x1 -- means looping (connecting opposite) facets in the X direction, 0x2 - in the Y direction.
///    If cap is not zero, the nonlooped directions are capped: a facet connecting all boundary edges is added.
FacetedSurfaceData *CreateFacetGrid(vec_type Lx, vec_type Ly, size_t nx, size_t ny, int loop=0, int cap=0);

///\en Creates (new) data for a cylinder with the base plane intersecting start and norlal to L, radius \a R, \a n_segments is the number of angular segments.
///
FacetedSurfaceData *CreateFacetedCylinder(const Vector_3& start, const Vector_3 &L, vec_type R, size_t n_segments, int cap=0);

///\en Creates (new) data for a faceted cube with the side X
inline FacetedSurfaceData *CreateFacetedCube(vec_type side = 1.){
  return CreateFacetedCylinder(Vector_3(side/2,side/2,0),Vector_3(0,0,side),side/sqrt(2.),4,1);
}

const Vector_3 cube_vertices[] = { Vector_3(0,0,0),Vector_3(1,0,0), Vector_3(1,1,0), Vector_3(0,1,0), 
                                   Vector_3(0,0,1),Vector_3(1,0,1), Vector_3(1,1,1), Vector_3(0,1,1)};
const size_t cube_index[]= { 4,0,1,2,3,  4,0,4,5,1, 4,1,5,6,2,   4,2,6,7,3,   4,3,7,4,0,  4,4,7,6,5}; 


const struct cube_storage_t {
  FacetedSurfaceData data;
  cube_storage_t(){
    data.points.assign(cube_vertices,cube_vertices+8);
    data.index.assign(cube_index,cube_index+30);
    data.num_contours = 6;
  }
} cube_facets;

///\en Alternative to \ref CreateFacetedCube using predefined constant facet data  
inline const FacetedSurfaceData *GetFacetedCube() {
  return &cube_facets.data;
}


template< class transform_t >
FacetedSurfaceData * FacetedSurfaceGenerator<transform_t>::CreateRegionSection(const Region_3 &reg, const vector<bool> *points_inside) const {

  size_t res_ind=0; // first point to be replaced
    
  //1. Creating result data
  FacetedSurfaceData *result_data = new FacetedSurfaceData;
  bool found_external = false;
  bool found_internal = false;
  
  vector<bool> tvec;
  const vector<bool> &point_tests = points_inside ? *points_inside : tvec;
  
  if(!points_inside){
    for(size_t i=0;i<data->points.size(); i++){
      tvec.push_back(reg.TestPoint(GetPoint(i)));
    }
  }
    
  for(size_t i=0;i<point_tests.size(); i++){
    if(point_tests[i]){        
      result_data->AddPoint(GetPoint(i));
      found_internal = true;
    }
    else{
      if(!found_external){
        found_external = true;
        res_ind = i;
      }
      result_data->AddPoint(Vector_3());
    }
  }
  if(!found_internal){ // faceted surface is outside region
    result_data->points.clear();
    return result_data; // returning void surface
  }
  if(!found_external){ // faceted region is fully inside reg: copying all contours
    result_data->index = data->index;
    result_data->num_contours = data->num_contours;
    return result_data;
  }

  typedef map< pair<size_t, size_t>, size_t> edgemap_t;
  typedef typename edgemap_t::iterator edgemap_it;
  edgemap_t edgemap;
  map<size_t, size_t> capmap;
  vector<size_t> cnt_points;       
                   
  //2. inspecting edges
  edge_index_it eib = edge_index_begin(), eie = edge_index_end(); // all edges as index pairs 
  int loopid=0; 
  bool entered = false;
  int first_entry = -1;
  size_t cross_ind, entry_ind,  exit_ind;
  while(eib!=eie){
    bool finish_contour=false;
    bool have_inner = false;
    
    
    size_t nedges = get_cur_contour_size(eib);
    edge_index_it eib0 = eib;
    for(size_t i=0;i<nedges;i++){
      pair<size_t,size_t> iedge=*eib, aedge;
      bool add_edge = false ;
      bool extra_point = false;
      if(point_tests[iedge.first] && point_tests[iedge.second]){ // both in: will be added at second traversal
        if(entered || loopid==1){ // put this to result
          add_edge = true;
          aedge=iedge;      
        }
        have_inner = true;
      }
      // entry or exit from region -- find first entry and subsequent exit
      if((!point_tests[iedge.first] && point_tests[iedge.second]) || (point_tests[iedge.first] && !point_tests[iedge.second] && entered) ){ // enter or exiting previously enetered
        if(!entered && first_entry>=0 && iedge.second==(size_t)first_entry){
          finish_contour = true; 
          capmap[entry_ind]=exit_ind;
        }
        else{
          // inspecting database
          pair<size_t,size_t> iquery= iedge.first< iedge.second ? make_pair(iedge.first, iedge.second) : make_pair(iedge.second, iedge.first);
          edgemap_it it = edgemap.find(iquery);
          if(it==edgemap.end()){
            // finding intersection
            // should be cross position
            Vector_3 v1= point_tests[iedge.first] ? result_data->points[iedge.first] : GetPoint(iedge.first);
            Vector_3 v2= point_tests[iedge.second] ? result_data->points[iedge.second] : GetPoint(iedge.second);
            
            vec_type frac[2], cross_pos=-1;
            int ncross= reg.TestLine(v1,v2-v1,frac);
            // now finding out the actuall cross position
            for(int i=0;i<ncross;i++){
              if(frac[i]<0)
                continue;
              if(frac[i]>1.) // this and all others beyond second point
                break;
              if(cross_pos>0)  // entry and exit: skipping this configuration
                cross_pos=-1.; 
              else
                cross_pos=frac[i]; // entered
            }
            if(cross_pos<0) // TestEdge did not work, we must assign some intersection anyway
              cross_pos=numeric_limits<vec_type>::epsilon(); // assuming the inersection was close to the first point
            Vector_3 vcross = v1+(v2-v1)*cross_pos;
            
            cross_ind = res_ind;
            
            // registering new point
            if(res_ind < result_data->points.size() ){
              result_data->points[res_ind++] = vcross;
              for(;res_ind < point_tests.size();res_ind++){
                if(!point_tests[res_ind])
                  break;
              }
            }
            else{
              result_data->points.push_back(vcross);
              res_ind++;
            }
            // putting into database
            edgemap[iquery]= cross_ind;
          } // not found in database
          else{ // found
            cross_ind = it->second;
          }
          // adding
          add_edge = true;
          if(!entered) { //  not enetred yet
            aedge = make_pair(cross_ind,iedge.second);
            extra_point = true;
          }
          else // exiting
            aedge = make_pair(iedge.first,cross_ind);

          if(!entered){
            
            entered = true;
            if(first_entry<0){
              first_entry = (size_t) iedge.second;
              entry_ind = cross_ind;
            }
            else
              capmap[cross_ind]=exit_ind;
          }
          else{
            // refreshing map element at exit with known entry point
            exit_ind =cross_ind;
            entered = false;
            //finish_contour = true;  
          }
        } // first entry
      }

      if(add_edge){
        if(!cnt_points.size() || extra_point) // for inner contours skip the first point (duplicated by edge_index iterator)
          cnt_points.push_back(aedge.first); // first in 1st edge
        cnt_points.push_back(aedge.second);
      }
      
      if(finish_contour)
        for(;i<nedges;i++)++eib;
      else
        ++eib;

    } //i
    if(loopid>=1)
      finish_contour = true; // put all inner

    bool looped_edges = false;
    if(finish_contour){
      result_data->AddContour(cnt_points.size(), cnt_points.begin());
      cnt_points.clear();
    }
    else{
      //looping iterator
      if(first_entry>=0 || (have_inner && loopid<1)){
        eib=eib0;
        loopid++;
        looped_edges = true;
      }
    }
    if(!looped_edges){ // switching off looping
      loopid=0;
      entered = false;
      first_entry = -1;
    }
    
  }// all contours
  // making caps
  for(map<size_t, size_t>::iterator it = capmap.begin(); !capmap.empty() ; ){
    cnt_points.push_back(it->first);
    size_t next_point = it->second;
    capmap.erase(it);
    it = capmap.find(next_point);
    if(it ==capmap.end()){ // next cap
      result_data->AddContour(cnt_points.size(), cnt_points.begin()); 
      cnt_points.clear();
      it = capmap.begin();
    }
  }
  
  
  return result_data;
}


#endif
