#ifndef SC_FIX_H
#define SC_FIX_H

/* @file sc_fix.h Fix classes for getting value for detectors, interpolation between meshes, 
   boundary conditions. 
   Documentation to the general structure of the code can be found at doc/tutorial.doxc
*/

#include "region_3.h"
#include "component.h"
#include "transfer.h"
#include "sc_media_array.h"
#include "sc.h"

typedef Interpolation<valtype,indtype> scInterpolation;
typedef InterpPacker<valtype,valtype,indtype> scInterpPacker;


// calculate electron (car=0) or hole (car=1) concentration using Scharfetter-Gummel approximation
// n0, n1, psi0, psi1 are concentration and potential at some points
// coef identify position between these points (0 <= coef <= 1) where concentration should be calculated
valtype SG_approx(valtype coef,valtype n0,valtype n1,valtype psi0,valtype psi1,int car);

/* calculate potential, electron (car=0) or hole (car=1) quasi-fermi pontential and concentration 
using 1D or 2D Scharfetter-Gummel approximation, described in 
A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012),
see Fig.3 and right column of page 2131 there.
Point is interpolated by values of 2 other points or 4 other points which are vertices of a square
 */
template<class cont_t>
int SG_interp(const scInterpolation &form,int car,const local_t &med,
  const cont_t *cont,valtype *psi,valtype *fermi,valtype *conc);

/*
This class is used to collect data (potential, concentration, etc) in chosen points 
and store it in specified buffer(s)
which can be used by detectors (see detector.h) to produce output text files
To assign buffers use register_buffer and substitute buffer
To assign points in which data will be collected use put_full_request
To record data to buffer(s) use compute_local
*/
template<class container_tt>
class scFixInterp: public apComponent{
public:
  typedef container_tt container_t;
  typedef valtype val_t; ///<\en type of interpolated value (for example, double or complex)

protected:

  container_t *cont; // scBlockContainer which interpolates data at chosen points
  vector<val_t *> buff; // buffers to store data

  // information about assigned points

  vector<int> ibufs; // buffers numbers
  vector<int> stinds; // indices in buffers

  scInterpPacker forms; // packer of interpolations for scalar variables
  scInterpPacker j_forms; // packer of interpolations for current

  scMediaArray meds; // medium properties
  vector<valtype> gens; // generation rate values

//  valtype delta; // coord increment for E and J calculation
//  vector<pair<scInterpolation,int> > idelta; // for E-field

  // bit flags for argument and output data type (see argTYPE and scoutTYPE)
  // this flag defines which data will be stored in buffers
  int argtype,outtype;

  // Timer used to measure time consumed for collecting data from container.
  // Timers concept is implemented in the base class apComponent.
  // Here we just assign this timer by overriding function apComponent::usingAddDefaultTimers
  struct{
    int tr_lcl;
  } tid;

public:
  scFixInterp(container_t *cont_=NULL);

  void AddDefaultTimers(const vector<int> &p_id);

  // identify node which mantains given position (reserved for MPI-version)
  int test_rank(Vector_3 &pos) const{
    return cont->test_rank(pos);
  }

  // register a buffer to record data
  // the pointer is not managed (must be deleted externally)
  // return the buffer index
  int register_buffer(val_t *bufp){
    buff.push_back(bufp);
    return (int)(buff.size()-1);
  }

  // substitute buffer with some buffer index
  valtype *substitute_buffer(int bufnum,val_t *ptr){
    val_t *tmp=buff[bufnum];
    buff[bufnum]=ptr;
    return tmp;
  }

  /* record interpolation and all other neccessary information in given point pos,
   and associate it with index ind in buffer ibuf
   increment value of index ind for the next record
   if rec_out=0 then if position outside container return -1, otherwise record empty interpolation
   */
  int put_full_request(int ibuf, int &ind, const Vector_3 &pos, int rec_out=1);

  // collect data from container to chosen buffers
  int compute_local();

  // reserved for MPI-version
  int prepare_transfers(){return 1;}

  // reserved for MPI-version
  void start_transfers(){};

  // reserved for MPI-version
  void complete_transfers(){};

  size_t data_size() const{
    return forms.data_size() + j_forms.data_size() + meds.sizeInBytes();
  }

  size_t packed_size() const{return data_size();}

  // reserved for MPI-version
  int mpi_size(int out=0) const{return 0;}

  // return size of one record in bytes
  size_t get_full_data_size()const;

  void set_full_data(int outtype_){
    outtype=outtype_;
  }

  void clear(){
    buff.clear();

    ibufs.clear();
    stinds.clear();

    forms.clear();
    j_forms.clear();

    meds.clear();
    gens.clear();
  }
};

/*
This is collection of nodes in the area of intersection between lower and higher level meshes
which require special procedure to update corresponding elements of Jacobain,
see detailed documentation in scBlockContainer::calc_transfer_jacobian

To work with this class you should:
1. record external mesh (I) nodes and their interpolations by (X) nodes of other higher level mesh

2. call function init to
 - identify internal (O) nodes which have these external (I) nodes as adjacent
 - rearrange collection of all nodes in the order X nodes -> group of I nodes -> O nodes, 
 which is suitable for scBlockContainer::calc_transfer_jacobian
*/
template<class cont_t>
struct FixTransfer{

  cont_t *cont;

  // exteral mesh nodes which require interpolation from other mesh
  // (I nodes in documentation to scBlockContainer::calc_transfer_jacobian)
  vector<int> i_ind; // indices
  vector<int> i_mesh; // lower level meshes

  // groups of internal nodes which have these external mesh nodes as adjacent
  // (O in documentation to scBlockContainer::calc_transfer_jacobian)
  vector<int> o_gr_ind;
  vector<int> o_grsz; // number of nodes in each group

  // interpolation of external mesh nodes by higher level mesh nodes
  // (X in documentation to scBlockContainer::calc_transfer_jacobian)
  scInterpPacker forms;

  vector<int> x_ind; // X nodes packed in the vector

  static const int max_interp=100; 

  struct ind_t{
    int igroup[max_interp]; // group of indices in I vector out
    int ogroups[max_interp]; // group of starting indices of groups in O vector o_gr_ind
    int n; // size of the group
    ind_t():n(0){}
  };
  // structures ind_t arranged in the order corresponding to x_ind
  vector<ind_t> backs;

  // used to store current values in I nodes, while we rewrite them by new values 
  // after shifting values at X and making interpolation in I.
  vector<valtype> toval;
  int tind; // currently read index in toval (used in set_previous_value)

  // this is auxiliary function is used only in init
  // map: key - X node, value - its position in vector backs
  // shift is X node
  // oi is position in I vector out
  // pn is position of starting index of group in O vector o_gr_ind
  void fill_back_item(map<int,int> &imap,int shift,int in,int ogn);

public:

  FixTransfer():cont(NULL),tind(0){}

  void set_cont(cont_t *cont_){
    cont=cont_;
  }

  int record(int indi, int meshi, scInterpolation &form){
    i_mesh.push_back(meshi);
    i_ind.push_back(cont->blocks[meshi]->offset+indi);
    forms.next_record(form);

    cont->blocks[meshi]->block->fill_adjacent_border(indi,o_gr_ind,o_grsz);

    return 1;
  }

  // fill x_ind and backs vectors
  int init();

  void reset_previous_value(){
    toval.clear();
    tind=0;
  }

  // set interpolated value for scalar variable(s) specified by representation flag rep to vi element of i_ind
  // update value of scalar variable specified by srp (usually concentration)
  void set_value(int vi,int rep,int srp,int remember=1);

  // set previous value there (using toval)
  void set_previous_value(int vi,int rep,int srp);

  size_t size(){
    return i_ind.size();
  }

  int clear();

};

/*
Jacobian updates for all types of equations different from Poisson and continuity (for example, boundary conditions) 
are organized in mesh independent way given by Fix classes. 
Note that unlike the traditional approach of including the specialized updates into the main update loop 
and testing each time the equation type flags, these updates are performed as separate Fix loops. 
Unlike the basic update, the Fixes are implemented in a general mesh-independent way.

To work with any fix you should:
1. record neccesary information for applying fix (usually node where fix is applied, 
its adjacent nodes, and sometimes some extra information)

2. call y to calculate residue of corresponding discretized equation
*/

// Base class for fix
template<class cont_t>
class scFix{
protected:
  cont_t *cont;

public:

  vector<int> out; // nodes where fix is applied

  vector<int> in; // groups of nodes which are adjacent to nodes out
  vector<int> in_sz; // sizes of each group

  scFix():cont(NULL){}

  void set_cont(cont_t *cont_){
    cont=cont_;
  }

  // record nodes where fix is applied (oindi, obi - index and mesh number of fix node)
  // function to specify and record group of adjacent nodes should be defined in derived classes
  int record(int oindi,int obi){
    out.push_back(cont->blocks[obi]->offset+oindi);
    return 1;
  }

  // return 1, if fix works for discertized equation (out=1) or variable (out=0) of representation rep
  int working_repr(int rep, int out)const{
    return 1;
  }

  void clear(){
    out.clear();
    in.clear();
    in_sz.clear();
  }
};

// this fix implements equation psi = psi_0, where psi_0 is some constant value
// residue is calculated as psi - psi_0
template<class cont_t>
class FixPsi: public scFix<cont_t>{

  vector<valtype *> psis; // applied voltage
  vector<valtype> psi_builts;
  vector<valtype> wf_deltas; // workfunction differences at boundary (chi_metal-chi) for Schottky contact, =0 for ohmic contact
public:

  using scFix<cont_t>::cont;
  using scFix<cont_t>::out;
  using scFix<cont_t>::in;
  using scFix<cont_t>::in_sz;

  int record(int oindi,int obi,valtype *psi, valtype psi_built, valtype wf_delta);

  valtype y(int vi,int repr);

  // works only for equation for potential
  int working_repr(int rep, int out)const{
    return rep&rpPsi;
  }

  void clear(){
    psis.clear();
    psi_builts.clear();
    wf_deltas.clear();
    scFix<cont_t>::clear();
  }
};

/*
Fix for boundary condition which connects current towards the interface and concentration at the interface
with given surface recombination velocity
*/
template<class cont_t>
class FixSR: public scFix<cont_t>{

  scInterpPacker concs; // interpolation for node at the boundary
  vector<CarrierParams> concs0; // equilibrium concentration at the boudnary
  vector<valtype *> srs; // surface recombination velocity
  // workfunction differences at boundary (chi_metal-chi) for Schottky contact, =0 for ohmic contact
  vector<valtype> wf_deltas;
  vector<int> bound_type; // boundary type: 1 - metal, 2 - dielectric

  // flag specifying how to calculate current towards to the surface:
  // 0 - interpolate by values of the current at mesh nodes
  // 1 - use function SG_interp
  // It is better to use calc_cur=1 for curved interface, 
  // however, we still keep possibility to chose calc_cur=0 for test purposes
  int calc_cur;

  vector<int> ob; // mesh number of each boundary node
  scInterpPacker Js; // interpolation for current towards the surface
  vector<int> bjs; // index of J in bondary mesh array

  scInterpPacker nforms; // nearest point (for direct current calculation)
  vector<valtype> deltas; // space difference to nearest point

  scMediaArray meds; // media for direct current calculation

public:

  using scFix<cont_t>::cont;
  using scFix<cont_t>::out;
  using scFix<cont_t>::in;
  using scFix<cont_t>::in_sz;

  FixSR():calc_cur(1){}

  /* oindi, obi - index and mesh number of boundary node
   bj - index in boundary array
   eq - boundary type (1 - metal, 2 - dielectric)
   center, near, nvect, surfp - point B, point O, nornal vector to the surface point S,
   see documentation to scBlockContainer::find_normal_inside and 
   A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012), Fig.3 and right column of page 2131.
   material parameters at the surface:
   wf_delta - workfunction difference (chi_metal-chi) for Schottky contact, =0 for ohmic contact
   med - medium
   dop - doping
   sr - surface recombination for electrons and holes (array of two elements)
  */
  int record(int oindi,int obi,int bj,int eq,const Vector_3 &center, const Vector_3 &near, const Vector_3 &nvect,
    const Vector_3 &surfp, valtype wf_delta, const local_t &med, valtype *sr);

  valtype y(int vi,int repr);

  int working_repr(int rep, int out)const{
    return out==1 ? rep&rpC : 1;
  }

  void clear(){
    concs.clear(); 
    concs0.clear();
    srs.clear();
    wf_deltas.clear();
    bound_type.clear();

    ob.clear();
    Js.clear();
    bjs.clear();

    nforms.clear();
    deltas.clear();
    meds.clear();

    scFix<cont_t>::clear();
  }
};

/* fix for reflective boundary conditions \nabla var * n = 0
where var is scalar variable (potential or quasi-fermi potential),
and n is normal to the boundary surface.
reflective boundary conditions are used for dielectric interface with zero surface recombination.
*/
template<class cont_t>
class FixNeumann: public scFix<cont_t>{

  scInterpPacker forms;
  scMediaArray meds;
  int wrepr;

public:

  using scFix<cont_t>::cont;
  using scFix<cont_t>::out;
  using scFix<cont_t>::in;
  using scFix<cont_t>::in_sz;

  FixNeumann():wrepr(rpC|rpP){}

  void set_wrepr(int repr_){
    wrepr=repr_;
  }

  int working_repr(int rep, int out)const{
    return rep&wrepr;
  }

  int record(int oindi,int obi,scInterpolation &form, const local_t &med);

  valtype y(int vi,int repr);

  void clear(){
    forms.clear();
    meds.clear();
    scFix<cont_t>::clear();
  }
};

// fix for reflective boundary conditions \nabla var * n = 0
// this fix doesn't work well, wee keep it for test purposes only
template<class cont_t>
class FixNeumannDer: public scFix<cont_t>{

  scInterpPacker cforms;
  scInterpPacker nforms[3];
  vector<Vector_3> norms;
  valtype dr;
  int wrepr;

public:

  using scFix<cont_t>::cont;
  using scFix<cont_t>::out;
  using scFix<cont_t>::in;
  using scFix<cont_t>::in_sz;

  FixNeumannDer():dr(1e-10),wrepr(rpC|rpP){}

  void set_wrepr(int wr){wrepr=wr;}

  int record(int oindi,int obi,const Vector_3 &center,const Vector_3 &norm);

  valtype y(int vi,int repr);

  int working_repr(int rep, int out)const{
    return rep&wrepr;
  }

  void clear(){
    cforms.clear();
    for(int i=0;i<3;i++)
      nforms[i].clear();
    norms.clear();
    scFix<cont_t>::clear();
  }
};

// This class is used to represent Poisson and continuity equations as a fix
// Used for calculation of these equations at the axis nodes in radial scheme
// Can be used for all nodes for test purposes (this option is turned on by emBlockContainer::if_fix_bulk)
template<class cont_t>
class FixBlock: public scFix<cont_t>{

  vector<int> oind; // output mesh index
  vector<int> ob; // output mesh number

  vector<int> ginds;
  vector<int> bind,bindsz; // all neighbour points (for block action)

  // 1 if fix is applied for axis notes in radial scheme
  int axis;

public:

  using scFix<cont_t>::cont;
  using scFix<cont_t>::out;
  using scFix<cont_t>::in;
  using scFix<cont_t>::in_sz;

  FixBlock():axis(0){}

  void set_axis(){
    axis=1;
  }

  int record(int oindi,int obi);

  valtype y(int vi,int repr);

  void clear(){
    oind.clear();
    ob.clear();

    ginds.clear();
    bind.clear();
    bindsz.clear();
    scFix<cont_t>::clear();
  }
};

#endif
