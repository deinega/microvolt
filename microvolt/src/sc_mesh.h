#ifndef SC_MESH_H
#define SC_MESH_H

// @file sc_mesh.h Declaration of scBlockContainer which is the main access class of MICROVOLT model
// Documentation to the general structure of the code can be found at doc/tutorial.doxc

#include "region_3.h"
#include "refobj.h"
#include "sc_carrier_params.h"
#include "physconst.h"
#include "sc_generation.h"
#include "sc_fix.h"

class RectBlock;

/*
Simulation is performed on one or multiple meshes which are collected in scBlockContainer.
Container manages interpolation and connection between different meshes. 
It also has boundary regions for which boundary conditions are defined. 
Boundary conditions determine the behavior when a certain mesh node lies outside the container boundary.

To work with this class you should:
1. Define calculated space (SetRegion), add meshes (AddMesh), 
medium regions, boundary regions, generation rate profile(AddMediumRegion, AddBoundaryRegion, AddGenerationRegion) 
and defining some settings (otherwise default settings will be used).
2. Organize meshes memory layout and assign coefficients values in discretized equations (Prepare)
3. Solve resulted discretized equations by Newton iterations until converged solution will be found (Step).

Field in given point of calculated space can be obtained using GetField.
Interpolation form (set of mesh nodes and interpolation coefficients) is obtained using create_interpolation.
*/
template <class block_tt=RectBlock>
class scBlockContainer: public apComponent, public restorer{

protected:

  typedef block_tt block_t;
  typedef scBlockContainer<block_t> container_t;

  // Older compilers don't support template parameters as friends (it's a C++11 feature)
  // but some tricks may work in a similar way
  // Some compilers (MSVC, Intel) don't provide a complete support of this feature
  // in their older versions, but allow "friend T" syntax
#if defined(HAS_FRIENDEXT) || defined(_MSC_VER) || defined(__INTEL_COMPILER)
  friend block_t;
#else
  // This line was tested on GCC 4.1 and 4.4
  friend class std::pair<block_t, int>::first_type;
#endif
  friend struct meth_cont_t;

  friend struct FixTransfer<container_t>;
  friend class scFixInterp<container_t>;
  friend class scFix<container_t>;
  friend class FixPsi<container_t>;
  friend class FixSR<container_t>;
  friend class FixNeumann<container_t>;
  friend class FixNeumannDer<container_t>;
  friend class FixBlock<container_t>;

  template<class cont_t>
  friend valtype *get_pointer(cont_t &cont, int rep);

  template<class cont_t> 
  friend void set_repres(cont_t &cont, int rp, int i);

//  template<class cont_t>
//  friend int update_repres(cont_t &cont,int rp,int i);

  friend int SG_interp<container_t>(const scInterpolation &,int,const local_t &,
  const container_t *,valtype *,valtype *,valtype *);

  int dim; // dimensionality of container (1 <= dim <= 3)
  static const int dim_max=3; // maximal dimensionality
  // flags for each dimension, -1 - not working, 0 - rectangular, 
  // 1 - radial in cylindrical system with rotational symmetry
  int dim_type[3];

  // in the case of rotational symmetry one can reduce dimensionality of the problem
  // using radial coordinate, for example (r,z) instead of (x,y,z)
  int rad; // radial coordinate (-1, if not used)
  valtype rad_axis; // radial coordinate of axis
  bool rad_right; // if true, calculated space is at the right of axis, otherwise it is at the left

  int per_flags; // bit flag indicating periodicity in specific directions: 0x1 for x, 0x2 for y, 0x4 for z
  Vector_3 p1, cell; // specify periodic elementary cell [p1, p1+cell] for directions activated by per_flags

  /*
  We solve Poisson and continuity equations with electrostatic potential 
  and quasi-Fermi potentials for electrons and for holes as independent variables.
  However, we also store electrons and holes concentrations and currents, since 
   - it is more convinient to formulate discretized equation using them 
   - it helps to reduce calculation of exponents
   - it easier to get concentrations and currents values in a runtime if it is stored in memory

  quasi-Fermi potentials (and therefore quasi-Fermi levels) are zero 
  for the case of zero applied voltage and no photogeneration
  however, for the output quasi-Fermi levels are assumed to be equal to 
  built-in potential plus intrisic Fermi level for some chosen medium in equilibrium (global_ifermi)
  this medium is found in zero_pos position
  
  We assume that all scalar variables (psi, fermi_n, fermi_p, concentration_n, concentration_p) 
  are stored at the same mesh nodes (for each mesh).
  Therefore interpolation made for one scalar variable (interpolation is a set of mesh nodes 
  and interpolation coefficients) is applied for other scalar variable.

  Currents for electrons and for holes could be stored at some other mesh nodes,
  not necessarily the same as occupied by scalar variables.
  */

  int SZ; // number of mesh nodes (of all meshes) for scalar variables
  int SZJ; // number of mesh nodes (of all meshes) for currents

  /* Below we define mesh arrays for scalar variables (currents).
   Size of each array is SZ (SZJ). 
   Each array is a collection of subarrays for different meshes.
   Starting index of each subarray is shifted from starting index of the array 
   on some offset value which is stored in mesh descriptor
  */

  valtype *psi_ptr; // array for potential values
//  valtype *sbpsi; // exp(psi*q/kT)

  valtype *fermi[2]; // array for n and p quasi-fermi potential values
  valtype *conc[2]; // array for n and p concentration values 

  valtype *excitons[2]; // array for triplet and singlet excitons
//  valtype *slotboom[2]; // slotboom variable

  valtype *J[2]; // array n and p current values

  valtype *ch; // array of doping concentration normalized on n_0
  valtype *gen; // array of generation rate

  // true if the initial values are estimated using linear extrapolation of
  // distributions potential and QF levels obtained at previous voltage levels
  bool useLinearEstimate;

  // Noise level for extrapolated distributions
  // Exact linear extrapolation sometimes produces an estimate that isn't good
  // as a starting point (convergence may be poor). Slightly perturbed values
  // may be better
  valtype estimationNoise;

  valtype prevVoltage; // previous value of the applied voltage

  // previous distribution of the variables (psi, fermi)
  // inner index is the node index (SZ elements),
  // outer index is the variable number (0 - psi, 1 - fermiN, 2 - fermiP)
  std::vector<std::vector<valtype> > prevVars;

  scMediaArray medc; // array of media properties at mesh nodes for scalar variables. size of this array is SZ

  // array of the size SZ which specifies if mesh node is internal (0) or belongs to the boundary 
  // of some type. Type of the boundary is characterized by some integer number (see scBoundaryRegion::type).
  // This array is used only for dump purposes
  int *bcond; 

  /* array of the size 3*SZ which consists of three subarrays (of the size SZ) corresponding to
   Poisson equation, continuity equation for electrons or continuity equation for holes.
   Typical value of array element is 0.
   Value 1 is used if value at corresponding mesh node is fixed for some equation 
   (for example, potential value for ohmic contact).
   Value 2 is used if value at corresponding mesh mode is interpolated by some other mesh.
   Mesh nodes with positive values (1, 2) are excluded form calculation of Jacobian for corresponding 
   equation, because values th these mesh nodes are fixed or interpolated by other mesh nodes.
  */
  int *fixed;

  // mesh descriptor
  template <class block_t>
  struct sc_bdescr: public restorer{
    mngptr<block_t> block; // mesh
    int id; // number in vector of meshes
    int level; // mesh level which defines priority between meshes in the region where they intersect
    size_t offset; // mesh offset index in the array for scalar variables
    size_t offset_cur; // mesh offset index in the array for current
    sc_bdescr(mngarg<block_t> block_,int level_):block(block_),level(level_),offset(0),offset_cur(0){}
  };

  typedef sc_bdescr<block_t> sc_bdescr_t;
  // vector of specified meshes
  refvector<sc_bdescr_t> blocks;

  /* We use multidimensional Newton-Rhapson method to solve Poisson and continuity equations,
  see http://en.wikipedia.org/wiki/Newton-Rhapson_algorithm
  The idea of the method is to start with some initial guess x and iteratively move to solution:
  1. calculate dx as a solution of system of linear equations: dx * J = y,
    where J is Jacobian matrix delta_y/delta_x and y is residue
  2. update x: x = x + dx
  The iterative process is repeated until y will be small enough so we can consider x as a solution.
  */

  valtype *y_ar; // array for residue
  valtype *x_ar; // array for current value of x
  valtype *dx_ar; // array for calculated dx

  // scalar variables are shifted at this value to calculate Jacobian matrix delta_y/delta_x numerically
  valtype delta_x;

  /* Multidimensional Newton-Rhapson method can be optimized by "smoothening" 
   of dx value calculated from dx * J = y: if dx is too big, we make it smaller.
   We used smoothening procedure proposed in PC1D program: see PC1D Help -> Numerical Method, equations A.21 - A.22.
   Clamp is a parameter of this smoothening procedure, case clamp = 0 corresponds to no smoothening
  */
  valtype clamp;

  /* Linear solvers to store matrix elements of different Jacobians and solve system of linear equations: dx * J = y.
   p corresponds to Jacobian for Poisson equation, 
   c - continuity equation (for electrons or for holes),
   c2 - both continuity equations,
   pc - Poisson and some continuity equation,
   pc2 - Poisson and both continuity equations.
   pc2 is the biggest solver, other solvers use memory managed by pc2 solver */
  lin_solver ls_pc2, ls_p, ls_c2, ls_pc, ls_c;

  // this parameter is used to check if Newton iterations for all equations 
  // converged to some minimum (global, which means that we found solution, or local, which is undesirable).
  int jacstop;

	// residue for methods in chains 
  valtype rsd2_stored[10];

  int ch_it; // count how many times function Step was called with current method

  valtype residue_threshhold; // residue value which is treated as small (see stop_ratio and small_stop_ratio)

  // If the residue decreases by less than stop_ratio or small_stop_ratio
  // then we assume that further iterations will yield no improvement of accuracy
  // If the residue is greater than residue_threshhold, then we use stop_ratio
  // (iterations stall and don't converge to an accurate solution),
  // otherwise we use small_stop_ratio (iterations are likely to
  // converge to an accurate solution; rounding errors may significantly influence
  // the residue value, thus we set a stronger condition to the decrease of the residue
  // in order to avoid useless iterations when a good accuracy is already achieved)
  valtype stop_ratio, small_stop_ratio;

  // These are experimental parameters for changing sequences of methods to improve convergence
  // Normally vector of methods consists of only one method, and methi is always 0
  vector<meth_chain_t> methods; // methods chains
  int methi; // current methods chain

  /* 
  As initial guess we put potential equal to built-in potential.
  Quasi-Fermi potentials are put to be zero that
  corresponds to equilibrium concentration (gen_eq = 0).
  However, if generation rate is nonzero, sometimes (but sometimes not) it could be beter to have initial concentration
  equal to equilibrium concentration plus rate of generated carriers multiplied on their timelife (gen_eq = 1).
  0 < gen_eq < 1 corresponds to intermediate values of concentration.
  */
  valtype gen_eq;

  /* It was shown that using semiconductor of high quality in the depletion region 
  around the pn- or hetero- junction improves open-circuit voltage V_oc of the solar cell, 
  even if the silicon quality elsewhere is small
  (but it does not influence short-circuit current J_sc),
  see A. Deinega, S. John, J. Appl. Phys. 112, 074327 (2012).
  To demonstrate this we can simulate solar cell using a high diffusion length 
  in an area confined within some distance depl_decr from the junction.
  If depl_decr=0 (default value) diffusion length around the junction is the same as elsewhere.
  */
  Vector_3 depl_decr;

  const scRegionMediumMap* media_regions;
  
  // list of specified boundary regions, int number corresponds to level of the boundary region 
  multimap<int,mngptr<scBoundaryRegion> > boundaries;

  // list of specified generation regions, int number corresponds to level of the generation region 
  multimap<int, mngptr<scGenerationRegion> > gens;
  
  // if for ohmic contact 
  int use_drc[2];
  /* which fix is applied for reflective boundary conditions \nabla var * n = 0
   where var is scalar variable (potential or quasi-Fermi potential),
   and n is normal to the boundary surface.
   reflective boundary conditions are used for dielectric interface with zero surface recombination.
   case 0 corresponds to applying FixSR (works only for quasi-Fermi potentials), 
   1 - FixNeumann, 2 - FixNeumannDer.
   It is better to use FixSR for curved interface (FixNeumann can lead to unphysical solution) */
  int refl_psi, refl_fermi;

  // this is test option, if if_fix_bulk=1, then mesh bulk update 
  // for Poisson and continuity equations will be performed by fix_block
  int if_fix_bulk;

  // below there are objects to handle equations different from Poisson and continuity equations
  // these objects update their part of Jacobian
  // documentation could be found in sc_fix.h

  FixPsi<scBlockContainer<block_t> > fix_psi; // potential at metal contact
  // boundary condition for dielectric and metal with given surface recombination
  FixSR<scBlockContainer<block_t> > fix_sr;
  // different types of realization for reflective boundary conditions
  FixNeumann<scBlockContainer<block_t> > fix_neum_psi, fix_neum_fermi;
  FixNeumannDer<scBlockContainer<block_t> > fix_der_psi, fix_der_fermi;
  FixBlock<scBlockContainer<block_t> > fix_block; // bulk mesh cycle made as a fix
  FixBlock<scBlockContainer<block_t> > fix_block_cyl; // part of mesh block cycle for cylindrical axis

  // collection of external mesh nodes and their interpolation by internal nodes of other mesh.
  // used in calc_transfer_jacobian to calculate corresponding part of Jacobian
  FixTransfer<scBlockContainer<block_t> > fix_trans;

  /* In order to avoid large numbers in calculations, 
  variables in equations are scaled using coefficients
  taken from Vasileska, Computational Electronics, Table 3.1 (page 39):

  electrostatic and quasi-Fermi potentials and measured in J0 = kB/q*T
  length is measured in l0 = sqrt(eps0*kB*T/q/q/n_i)
  concentration is measured in n0 = n_i
  time is measured in t0 = l0*l0/D0
  (see 1,2,3 lines in Table),

  eps0, kB, q are defined in physconst.h and correspond to \varepsilon, k_B, q from Table
  T is temperature
  n_i is intrisic concentration of the filling medium
  D0 is some diffusion coefficient value, f.e. 1e-2 m2/s

  1. Rescaling for Poisson equation
  eps * d^2 psi / dr^2 = q * n
   or
  eps0/q * d^2 psi / dr^2 = n
  becomes
  eps0/q * d^2 (psi*J0) / (dr*l0)^2 = n*n0
  or 
  coefpsi * d^2 psi / dr^2 = n, 
  where coefpsi = eps0/q*J0/(l0*l0)/n0
  for chosen values J0, l0 and n0 (see formulae above), coefpsi = 1

  2. Rescaling for continuity equation (Scharfetter-Gummel discretization)
  D * d^2 n / dx^2 * exp(q/(kB*T) * psi) = G - R
  becomes
  D/t0 * d^2 (n*n0) / (dr*l0)^2 * exp(q/(kB*T) * psi*J0) = (G - R)*n0/t0,
  or after multiplying on t0^2
  coefc * D * d^2 n / dr^2 * exp(coefexp * psi) = coefr(G - R),
  where coefc = t0/(l0)^2/n0, coefr = t0/n0, coefexp = q/(kB*T)*J0
  for chosen values J0, l0 and n0 (see formulae above) coefexp = 1
  */
  valtype J0,l0,n0; // scaling units for potentials, length and concentrations
  valtype lu; // output length unit in meters (f.e., if lu=1e-6 output length unit is micron)
  valtype T; // temperature

  valtype coefpsi; // coef in the left handside of the Poisson equation
  valtype coefc; // coef in the left handside of the contiuity equation
  valtype coefr; // coef in the right handside of the contiuity equation
  valtype coefexp; // coef in exponenta before psi or fermi level in Scharfetter-Gummel discretization

  // used only for output:
  // if zero_pos!=VEC_INFTY, in the case of zero applied voltage and no photogeneration 
  // electrostatic potential = 0 and fermi level = global_ifermi in zero_pos 
  Vector_3 zero_pos; 
  valtype global_ifermi;

  typedef std::map<int,virt_unary_function<const Vector_3 &, vec_type>* > ini_funcmap_t;
  // Map that stores functions loading initial values at Prepare stage of calculation. The map keys are repesentations.
  ini_funcmap_t ini_funcmap;
  /// initial guess value for fermi levels
  mngptr<virt_unary_function<const Vector_3 &, vec_type> > fermi_guess[2];

  scInterpolation tmpi; // temporary interpolation
  static const int max_interp=100; // maximal size of the temporary interpolation
  vector<int> pshifts; // index array for temporary interpolation
  vector<valtype> pcoeffs; // coefficients array for temporary interpolation

  /* We use timers to measure time consumed by some chosen procedures.
  Using timers allow to find more time consuming parts of the code and optimize them.
  Timers concept is implemented in the base class apComponent.
  Here we just assign some timers by overriding function apComponent::usingAddDefaultTimers
  */
  struct tid_t{
    int specify_nodes,analyze_geometry,prepare_transfers,dump,linsys,
      block_jacobian,transfer_jacobian,jacobian_linsys,jacobian_step,jacobian_transfer_values,matrix_prepare;
  } tid; // numbers of timers for each chosen procedure

  void AddDefaultTimers(const vector<int> &p_id);

  // clear vector of specified medium / boundary / generation regions
  int clean_geometry();

  // return mesh nodes classification bit flag for chosen mesh node
  int classify_node(const Vector_3 &c,int lev);

  // this is auxiliary function which is called from specify_nodes.
  // check adjacent nodes of the node pit
  // and mark those of them which are internal / boundary / transfer in map memmap.
  // adjacent nodes which should be analyzed by specify_adjacent_nodes will be recorded to vector vp.
  // bflag is used only for mFixBoundSurround nodes (see documentation for mFixBoundSurround)
  int specify_adjacent_nodes(sc_bdescr_t *bdescr,typename block_t::mesh_it &pit,int bflag,
    vector<pair<typename block_t::mesh_it,int> > &vp, map<int,pair<int,int> > &memmap);

  /*
  Each mesh has iterator mesh_it on all its nodes.
  In this function container iterate mesh_it and analyze position of each mesh node 
  compare to regions where fixes should be applied (boundaries) and other meshes of heigher level.

  Those nodes which lie inside mesh are marked as internal nodes (mInt).
  They participate in bulk Jacobian update for Poisson and continuity equations.

  Those nodes which are adjacent to internal nodes and lie in the area managed by higher level meshes
  are marked as transfer nodes (mTrans).
  To update corresponding part of Jacobian FixTransfer is used.

  Those nodes where fix should be applied (equation different from Poisson or contiuity) 
  are marked as fix nodes (mFix).
  These nodes include:
   - boundary fix nodes: mFixBound and mFixBoundSurround (boundary fix node 
       which surround other boundary fix nodes, not used now)
   - internal fix nodes: mFixInt

  All mesh nodes which do not fall to one of this category, should not be stored in the mesh.

  Each mesh should have map 
    with key - index which identify node location at the mesh, 
    value - pair of node memory index and flag mFlag

  This map is returned by mesh get_memmap and filled by this function scBlockContainer::specify_nodes

  After this map will be used by mesh function organize_memory to specify its memory layout.
  */
  int specify_nodes(sc_bdescr_t *bdescr);

  // finds which block manages interpolation of given point (if there is no such a block, returns -1)
  // block with id-number will be prefered among blocks with equal level
  // if strong=1 and interpolation cannot be done within one block only, returns -1
  int test_in_blocks(const Vector_3 &pos, int id=0, int strong=0) const;

  int test_in_blocks(pair<Vector_3, Vector_3> arg, int id=0, int strong=0) const{
    return test_in_blocks(arg.first,id,strong);
  }

  size_t sc_bdescr_t::*get_offset(const Vector_3 &pos){
    return &sc_bdescr_t::offset;
  }

  size_t sc_bdescr_t::*get_offset(const pair<Vector_3, Vector_3> &cur){
    return &sc_bdescr_t::offset_cur;
  }

  // this function is recursivelly called from create_interpolation
  // it fills elements of pshifts and pcoeffs arrays starting from index cur_ind
  // (elements before are assumed to be assigned already)
  template<class arg_t>
  int fill_interpolation_array(int cur_ind, const arg_t &arg, size_t sc_bdescr_t::*offset, int id=0,int forse_id=0,int incompl=0,int non_nonlocal=0);

  // create pools for the local media data and return references to their elements
  // The size is SZ, indexing is the same as for bcond
  int construct_local_media();

  // returns 0 (is point is not at the boundary) or boundary type (otherwise)
  // records boundary region parameters
  int test_bond(const Vector_3 &p,valtype **psi=NULL,valtype **sr=NULL, valtype *workf=NULL, Region_3 **reg=NULL)const;

  // return generation rate at given point
  valtype test_gen(const Vector_3 &p) const;

  /* 
  This function is used for numerical implementation of reflective boundary conditions,
  described in A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012),
  see Fig.3 and right column of page 2131 there.
  Input: point p (B at Fig.3) which belongs to sc_bdescr and lies outside region reg.
  Output: point surfp at the region surface (S at Fig.3) which is closed to p,
  normal nvect to this surface at point surfp,
  point inside region that is aligned with mesh and lies at the ray B - S (O at Fig.3)
  */
  int find_normal_inside(Region_3 *reg, sc_bdescr<block_t> *sc_bdescr, 
    typename block_t::node_t &p, Vector_3 &point, Vector_3 &nvect, Vector_3 *surfpp=NULL);

  // analyze mesh node position relatively to specified media and generation rate gen at this point
  // and send this information to mesh which will assign coefficients values in discretized equations
  // put some initial guess value for scalar variables
  // flag bord is used only in test regime when if_fix_bulk=1
  // ini_repr - representation where initial values are recorded
  int record_internal(sc_bdescr<block_t> *sc_bdescr,const typename block_t::node_t &p,int bord,valtype gen,int ini_repr=rpAll);

  // analyze mesh nodes position relatively to specified boundaries
  // and make a record in corresponding fix
  int record_bound(sc_bdescr<block_t> *sc_bdescr,const typename block_t::node_it &pit,typename block_t::node_t &p,int aflag);

  // this bit flag is used as a parameter of function analyze to specify which instructions should be done
  enum{
    scIniVal=0x1, // intialize mesh points by default values
    scIniVoltage=0x2, // apply voltage
    scIniGen=0x4, // apply generation
    scIni=0x8, // all other instructions
    scIniAll=0xfffff
  };

  /* analyze mesh nodes position relatively to specified media and boundaries. 
    This information is used to assign coefficients values in discretized equations
    for mesh cycle and fixes.
    aflag allows to turn on/off some steps of analyze (see corresponding enumeration)
    if GR!=NULL, record there generation rate integrated along calculated space
    ini_repr - representation where initial values are recorded
  */
  int analyze(int aflag=scIniAll, valtype *GR=NULL, int ini_repr=rpAll);

  int clear_geometry();

  /* 
  update external mesh nodes with values interpolated by internal nodes of other meshes
  flag = 0: just update
  flag = 1: update and remember current values at external mesh nodes
  flag = -1: restore these values if they were overwritten due to calling this function with flag=1
  */
  int transfer_values(FixTransfer<scBlockContainer<block_t> > &trans,int repr,int srp,int flag=0);

  // calculate residue for internal nodes of all meshes, record it to y_ptr and return sum y*y
  valtype calc_mesh_y(const meth_cont_t &mc,valtype *y_ptr);

  // calculate residue for chosen fix, record it to y_ptr and return sum y*y
  template<class fix_t>
  valtype calc_fix_y(fix_t &fix,const meth_cont_t &mc,valtype *y_ptr);

  // calculate total residue (internal nodes of all meshes and all fixes) and return sum y*y
  valtype calc_y(const meth_cont_t &mc,valtype *y_ptr);

  /*
  Discretized equation at each mesh node (lets call it central O node) is contributed by this node 
  and some adjacent mesh nodes (lets call them I nodes).
  Variables (potential, quasi-fermi potentials) values in O node and adjacent I nodes influence 
  residue value of disretized equation at O node.
  To update Jacobian delta_y/delta_x we shift scalar variable at O node and adjacent I nodes at some
  small value delta_x and calculate resulted change of residue delta_y in O node

  Node could be considered as O node or I node depending on in which node we solve discretized equation.

  Some of the mesh nodes could lie in the priority area of higher level mesh (lets call it mesh 2).
  These nodes can not be O nodes, because we do not solve discretized equation there.
  Value in this nodes should be interpolated by values of some nodes 
  of this other higher level mesh (lets call them X nodes of mesh 2).
  However, these nodes (of the lower level mesh 1) could be adjacent I nodes for some other O nodes of lower level mesh.
  Therefore to update Jacobian elements for O nodes we should
   - shift scalar variable in nodes X (of the higher level mesh 2) on some small value delta_x
   - calculate shifted interpolated value at nodes I (of the mesh 1)
   - calcuate shifted residue delta_y at nodes O (of the mesh 1)

  Jacobian update procedure for such type of nodes is implemented in function calc_transfer_jacobian
  */
  int calc_transfer_jacobian(FixTransfer<scBlockContainer<block_t> > &trans,const meth_cont_t &mc);

  // calculate jacobian for chosen fix
  template<class fix_t>
  int calc_fix_jacobian(fix_t &fix,const meth_cont_t &mc);

  // calculate Jacobian
  int calc_jacobian(const meth_cont_t &mc);

  /*
  calculate Jacobian matrix J = delta_y/delta_x and residue y,
  find dx as a solution of system of linear equations: dx * J = y,
  and update current values of scalar variables x: x = x + dx
  meth - Jacobian structure
  rsd - current total residue, rsd2 - new total residue (for updated x)
  return 1 if everything ok, 0 if Newton iteration is not moving, -1 if Jacobian matrix is well conditioned
  */
  int Jacobian(const meth_cont_t &mc,valtype &rsd,valtype &rsd2);

public:

  scBlockContainer(int *dim_type_, valtype T_);

  ~scBlockContainer();

  int get_dim()const{return dim;}
  
  /* Dump geometry of specified regions (medium, boundary, generation)
   to text files which could be plotted by gnuplot program.
   Dump settings are specified in the base class apComponent.
   Here we just assign regions to be dumped by overriding function apComponent::Dump
  */
//  void Dump(bool lim=false);

  /* create linear interpolation for some argument which could be position 
  (this case corresponds to scalar variable: potential, fermi-level or concentration),
  or pair (position, direction) (this case corresponds to current).
  block with id-number will be prefered among blocks with equal level
  if forse_id, then block with id-number will be chosen for sure
  if interpolation cannot be done, zero interpolation will be returned
  if incompl=1 then even interpolation is incomplete (not enough mesh nodes to make interpolation 
  for some point close to the boundary), it will be created (otherwise zeroth interpolation will be returned)
  if non_nonlocal=1 then if interpolation cannot be done within one mesh, zero interpolation will be returned
  if persistent=0 subsequent calls of create_interpolation may destoy the retuned value
   */
  template<class arg_t>
  scInterpolation create_interpolation(const arg_t &arg,int id=0,int forse_id=0,
    int incompl=0,int non_nonlocal=0,int persistent=0);

  // return interpolated value for variable corresponding to representation flag
  valtype get_value(const scInterpolation &form,int rep) const{
    return form.get_value(get_pointer(*this,rep));
  }

  valtype get_value(const scInterpolation &form,int rep,int is) const{
    return get_value(form,rep|(is ? rpH : rpE)|(is ? rpT : rpS)); // added excitons
  }

  // tests if value at this point can be interpolated
  bool TestPoint(const Vector_3 &pos){
    scInterpolation form=create_interpolation(pos,0,0,1); // weak condition for detectors
    return form.valid();
//    return test_bond(pos)==0;
//    return test_in_blocks(pos)>=0; // weak condition for detectors
  }

  // identify node which mantains given position (reserved for MPI-version)
  int test_rank(Vector_3 &pos) const;

  // rescaling discretized equations using length for potential (eV), length (lu), time and concentration ()
  // lu is user length unit in micron
  int calibrate(valtype j=1, valtype l=1, valtype t=1, valtype n=1, valtype lu=1);

  // calibration according to Vasileska for chosen temperature, concentration unit and user length unit
  int calibrate_Tn(valtype T, valtype n, valtype lu);

  void set_method(const meth_chain_t &meth,int clear=1){
    if(clear)
      methods.clear();
    methods.push_back(meth);
  }
  
  /// Set loader function for initial values for given representation. The func pointer is unmanaged.
  /// The values are loaded by record_internal() at Prepare stage.
  void SetIniValueFunction(int repr, virt_unary_function<const Vector_3 &, vec_type> *func){
    ini_funcmap[repr] = func;
  }

  /// adding initial guess value for fermi levels
  void AddFermiGuess(virt_unary_function<const Vector_3 &, vec_type> *g0, virt_unary_function<const Vector_3 &, vec_type> *g1){
    fermi_guess[0].reset(g0,0);
    fermi_guess[1].reset(g1,0);
  }
  
  /// load values from some table to array corresponding to chosen representation
  int load_ini_values(virt_unary_function<const Vector_3 &, vec_type> *table, int repr){
    for(size_t i=0;i<blocks.size();i++)
      blocks[i]->block->load_ini_values(table,repr);
    return 1;
  }

  int FillIniValues(int ini_repr){
    return analyze(scIniVal,NULL,ini_repr);
  }

  void set_convergence_check(valtype ratio=0.05, valtype threshold=-1, valtype small_ratio=0.5){
    stop_ratio=ratio;
    residue_threshhold = threshold;
    small_stop_ratio = small_ratio;
  }

  void set_clamp(valtype cl){
    clamp=cl;
  }

  void set_gen_eq(valtype gen_eq_){
    gen_eq=gen_eq_;
  }

  void set_depl_decr(Vector_3 dd){
    depl_decr=dd;
  }

  void SetZeroPosition(Vector_3 zero_pos_){
    zero_pos=zero_pos_;
  }

  void EnableLinearEstimate(bool val = true, valtype noiseLevel = 0) {
    useLinearEstimate = val;
    estimationNoise = noiseLevel;
  }

  // add mesh and specify its level
  int AddMesh(mngarg<block_t> ptr,int level=0);

  void ShareMediumRegions(const scRegionMediumMap& meds);

  // add boundary region and specify its level
  int AddBoundaryRegion(mngarg<scBoundaryRegion> bond, int level);

  // set radial coordinate, position of radial axis (r=0)
  // if rad_right = true (false) then positive (negative) values of radial coordinate are considered
  void SetRadialCoordinate(int rad_, valtype rad_axis_, bool rad_right_){
    dim_type[rad_]=1;
    rad=rad_;
    rad_axis=rad_axis_;
    rad_right=rad_right_;
  }

  // set bit flag indicating periodicity in specific directions: 0x1 for x, 0x2 for y, 0x4 for z
  // for periodic elementary cell [p1, p1+cell]
  void SetPeriodicBoundaries(const Vector_3 &p1_, const Vector_3 &cell_, int per_flags_){
    p1=p1_;
    cell=cell_;
    per_flags=per_flags_;

    for(int i=0;i<3;i++){
      if(dim_type[i]<0)
        per_flags&=~(0x1<<i);
    }
  }

  // add generation region and specify its level
  int AddGenerationRegion(mngarg<scGenerationRegion> gen, int level=0);

  // apply voltage to metal contacts where it should be applied
  int ApplyVoltage(valtype psi);

  // allocate memory and organize memory layout
  // dmp is bit flag regulating mesh nodes dumping (used for test purposes), first bit - dump mesh nodes, 
  // second bit - create different files for different mesh nodes (internal, boundary, transfer)
  int AllocateMemory(int dmp=0);

  // calculate and assign coefficients in discretized equations
  // if GR!=NULL, record there generation rate integrated along calculated space
  int BuildGeometry(valtype *GR=NULL, int ini_repr=rpAll);

  int Prepare(valtype *GR=NULL, int dmp=0){
    if(AllocateMemory(dmp)<0)
      return -1;
    return BuildGeometry(GR);
  }


  // make step for current chain methi of Jacobian structures
  // return 1 if everything ok, 0 if Newton iteration for all methods are not moving, 
  // -1 if Jacobian matrix is well conditioned
  int Step();

  // return total residue and print residues for Poisson and continuity equations (if print=1)
  // residue is normalized on number of mesh nodes (of all meshes) for scalar variables
  valtype get_residue(int print=0){
    valtype rsdp[3];
    rsdp[0]=calc_y(meth_cont_t(*this,method_t("p")),y_ar);
    rsdp[1]=calc_y(meth_cont_t(*this,method_t("e")),y_ar);
    rsdp[2]=calc_y(meth_cont_t(*this,method_t("h")),y_ar);
    if(print)
      message(vblMESS1,0,"\t%g %g %g ",rsdp[0]/SZ,rsdp[1]/SZ,rsdp[2]/SZ);

    return (rsdp[0]+rsdp[1]+rsdp[2])/SZ;
  }

  int DumpMeshes(const string &suf, bool dif, int outtype);

};

#endif
