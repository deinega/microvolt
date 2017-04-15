#ifndef FD_H
#define FD_H

// @file sc_fd.h Declaration of RectBlock 
// Documentation to the general structure of the code can be found at doc/tutorial.doxc

#include <algorithm>
#include "grid.h"
#include "sc.h"
#include "sc_media_array.h"

/*
This class implements finite difference scheme in 1D, 2D and 3D 
and radial finite difference scheme for structures with cylindrical symmetry
described in A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012).
It supports all functionality required for mesh to be included in scBlockContainer

In finite difference scheme each central node has 2*dim adjacent nodes, where dim is dimensionality of the mesh
*/
class RectBlock: public apComponent{

  template<class cont_t>
  friend valtype *get_pointer(cont_t &cont, int rep);

  template<class cont_t> 
  friend void set_repres(cont_t &cont, int rp, int i);

public:

  static const int dim_max=3;
  // maximal number of (1) interpolation nodes, (2) central node plus adjacent nodes 
  // used to define size of some static arrays
  // 9 = pow(2,dim_max)+1
  static const int vert_max=9;
  // the same as vert_max, but in one dimension
  // 3 = pow(2,1)+1
  static const int vert_dir_max=3;

  class node_t;

  // Mesh array is considered as a consequence of slices.
  // These slices are used as arguments of functions with suffix _slice
  // that perform a loop along the slice.

  struct slice_t{
    int oind; // memory index of the start central node (central nodes are always internal)
    int iind[6]; // memory indices of its adjacent nodes
    int ic[3]; // position on a grid of start cental node
    int gind; // geometry index of start central node (packed ic[3])
    int sz; // size of the slice
    slice_t(int oind_=-1, int gind_=1):oind(oind_),gind(gind_),sz(-1){
      for(int i=0;i<6;i++)
        iind[i]=-1;
    }
    bool operator <(const slice_t &other) const{
      return oind<other.oind;
    }
  };

protected:

  // ----------------

  // meaning of variables in the block below is the same as in scBlockContainer, see documentation there
  // these variables are assigned by scBlockContainer which calls function reset
  // arrays will be shifted on some offset value compare to corersponding global arrays in scBlockContainer

  int dim;
  int dim_type[3];

  int extsz; // scBlockContainer::SZ

  int offset; // offset of scalar variable arrays from corresponding global container arrays

  valtype *psi_ptr;

  valtype *fermi[2];
  valtype *conc[2];

  valtype *J[2];

  valtype *excitons[2]; // I.V. added

  scMediaArray medc;

  valtype *ch;

  valtype *gen;

  int *bcond;
  int *fixed;

  valtype *y_ar; // not shifted from global container array

  lin_solver *ls_pc2, *ls_p, *ls_c2, *ls_pc, *ls_c;

  valtype delta_x;

  valtype T;
  valtype J0,l0,n0;
  valtype lu;
  valtype coefpsi;
  valtype coefexp;
  valtype coefc;
  valtype coefr;

  valtype global_ifermi;

  // ----------------

  int SZ; // number of mesh nodes for scalar variables
  int N[3]; // mesh nodes number along each dimension
  int fd; // fastest dimension in array indexing

  scMediaArray medh; // additional array for media properties at mid-nodes (where current is calculated)

  valtype *recomb_ptr; // recombination rate value at mesh nodes (used in Jacobian calculation)

  typedef NonUniformGrid<valtype,1,indtype> grid_t;
  // grid which is used to obtain interpolation coefficients for arbitrary points inside mesh
  // and pack / unpack grid indices to single integer numbers
  grid_t grd;

  Box rbox; /// bounding box of the mesh

  // distance between mesh point and axis (used in radial scheme for structures with cylindrical symmetry)
  valtype *r_arr;
  valtype *dx_arr[3]; // mesh step for each dimension

  mngptr<Region_3> conf; // if conf!=NULL, grid will be confined within this region

  // map for mesh nodes classification
  // with key - index which identify node location at the mesh, 
  // value - pair of node memory index and flag mFlag
  // see detailed documentation in scBlockContainer::specify_nodes
  typedef map<int,pair<int,int> > memmap_t;
  memmap_t memmap;

  vector<slice_t> slices; // slices for internal nodes

  vector<slice_t> b_slices; // slices for boundary nodes (their size is always 1)

  // these arrays are used in packing / unpacking directions
  // unpacked direction is two numbers:
  // di, that takes values 0,1,2 (x,y,z), and ni, that takes values 0,1 (back and forward)
  int dnf[3],dnb[3];

  ///\en multipliers for gradient components in each direction, default is 1,1,1
  /// using other values can improve convergence, however, this question is still not answered
  Vector_3 grad_split;

  // pack direction
  inline int dnum(int di,int ni) const {
    return 2*dnf[di]+ni;
  }

  // unpack direction
  inline void dnum(int ind, int &di,int &ni) const {
    ni=ind%2;
    di=dnb[ind/2];
  }


  // allocate memory for internal arrays and intialize them, 
  // rad_axis is radial coordiante of axis which is used in the case of radial finite difference scheme
  int reset(valtype rad_axis);

  // determine context parameters (electric field, concentration et al.) at the given point
  scPointContext getContext(const slice_t& slice, int shift) const;

  // calculate recombination and generation rate at the specified mesh point
  void recomb_gen(const slice_t& slice, int shift, valtype &R, valtype &G);

  // calculate recombination and generation rate at mesh point with node with memory index i
  void recomb_gen(int i,valtype &R,valtype &G);

  // fill recomb_ptr array with current recombination values
  void fill_recomb();

  /* calculate Jacobian element
   \param oind - start node;
   \param nind - nodes which influence residue in start node
   \nindz - their number
   \param shift - node where residue is calculated is shifted on this value from start node
   the same is about nodes which influence residue in this node
   \mc - Jacobian structure
  */
  int Jacobian_node(int oind, const int *nind, int nind_sz, int shift, const meth_cont_t &mc);

  /**\en Calculate electron and hole current or its derivatives according to Scharfetter-Gummel approximation:
   \param slice - slice;
   \param shift - node shift from the beginning start of the slice
   \param schemeCoeffs - difference scheme coefficients, e.g.:
     {(0.5,0.5,0.5,0.5,0.5,0.5)} - for current (half sum of the adjacent current values),
     {(-1/h,-1/h,-1/h,1/h,1/h,1/h)} - for current derivatives on a homogeneous mesh;
   \param J - result values (J[0] for electrons, J[1] for holes).
   \param boundary - if current is calculated for boundary point
   */
  void calcJ(const slice_t &slice, int shift, const valtype* schemeCoeffs, Vector_3* J, int boundary) const;

  // record header (columns titles) to dump file
  int Dump_header(FILE *f, int outtype);

  // record node coordinates and data at this node to dump file
  int Dump_element(const slice_t &slice, int shift, FILE *f, int outtype);

public:

  RectBlock(const Box &B, const int *ndim);

  ~RectBlock();

  // initialize nonuniform grid direction
  // beg, end - iterator on nodes coordinates in this direction
  template<class inp_it>
  int SetMeshSteps(size_t aind, inp_it beg, inp_it end){
    Vector_3 p1=rbox.get_p1(),p2=rbox.get_p2();
    p1[aind]=*beg;
    p2[aind]=*(end-1);
    rbox.init(p1,p2);
    grd.SetMeshSteps(aind,beg,end);
    grd.GetInterpolationRange(N);
    return 1;
  }

  int SetConfinement(mngarg<Region_3> conf_){
    conf.reset(conf_);
    return 1;
  }

  // number of mesh nodes for scalar variables
  size_t GetSize() const {
    return SZ;
  }

  // number of mesh nodes for current
  size_t GetCurrentSize() const {
    return dim*SZ;
  }

  // test if point is inside the grid
  // if only=true, check if point belongs to confining region only
  bool TestPoint(const Vector_3 &pos, bool only=false) const {
    if(conf.ptr() && !conf->TestPoint(pos))
      return false;
    else if(only)
      return true;
    for(int i=2;i>=0;i--){
      if(dim_type[i]<0)
        continue;
      if(!acless(rbox.get_p1()[i],pos[i]) || !acless(pos[i],rbox.get_p2()[i]))
        return false;
    }
    return true;
  }

  // test if point can be interpolated with this mesh nodes only
  bool TestInterpolation(const Vector_3 &pos){
    int gind[vert_max], mind[vert_max];
    valtype coeff[vert_max];
    vector<Vector_3> nonlocal;
    int n=create_interpolation(pos,gind,mind,coeff,&nonlocal);
    return n>0 && nonlocal.size()==0;
  }

  void SetDirSplit(const Vector_3 &grad_split_){
    grad_split=grad_split_;
  }

//  memmap_t &get_memmap(){ // used in container
//    return memmap;
//  }

  void set_memmap(const memmap_t &memmap_){
    memmap=memmap_;
  }

  /* Gets interpolation coefficients (gind - geometry indices, mind - memory indices, coeff - coefficients in interpolation).
   the arrays must be at least 8 elements long.
   the negative gind and mind indicies correspond to absent grid points,
   the space location of kth absent point is pushed_back to nonlocal[-int_ind[k]-1].
   if incompl=1 then even interpolation is incomplete (not enough mesh nodes to make interpolation
   for some point close to the boundary), it will be created (otherwise -1 will be returned).
   returns number of elements in interpolation
   */
  int create_interpolation(const Vector_3 &pos, int *gind, int *mind, valtype *coeff, vector<Vector_3> *nonlocal, int incompl=0);

  // Get interpolation coefficients for current defined as a pair position-direction
  int create_interpolation(const pair<const Vector_3,const Vector_3> &arg,int *gind, int *mind, valtype *coeff, 
    vector<pair<Vector_3,Vector_3> > *nonlocal, int incompl=0);

  // see detailed documentation in scBlockContainer::find_normal_inside.
  // for given point p and direction surfn finds point which is aligned with mesh
  int find_point(const node_t &p,const Vector_3 &nvect, Vector_3 &point);

  // assign all parameters that are defined directly by scBlockContainer
  // allocate memory for internal arrays and intialize them
  template<class cont_t>
  int reset(cont_t &cont, int offset_, int offset_cur){
    if(!SZ)
      return 0;

    if(cont.rad>=0)
      dim_type[cont.rad]=1;

    offset=offset_;
    extsz=cont.SZ;

    psi_ptr = cont.psi_ptr + offset;

    for(int car=0;car<2;car++){
      conc[car] = cont.conc[car] + offset;
      fermi[car] = cont.fermi[car] + offset;
      J[car] = cont.J[car] + offset_cur;
      excitons[car] = cont.excitons[car] + offset; // I.V.
    }

    ch=cont.ch + offset, gen=cont.gen + offset;
 //   medc.map_to(cont.medc, offset, SZ);
    bcond=cont.bcond + offset, fixed=cont.fixed + offset;
    y_ar=cont.y_ar;

    ls_pc2=&cont.ls_pc2, ls_p=&cont.ls_p, ls_c2=&cont.ls_c2, ls_pc=&cont.ls_pc, ls_c=&cont.ls_c;

    delta_x=cont.delta_x;
    T=cont.T;
    J0=cont.J0, l0=cont.l0, n0=cont.n0, lu=cont.lu;
    coefpsi=cont.coefpsi, coefexp=cont.coefexp, coefc=cont.coefc, coefr=cont.coefr;
    return reset(cont.rad_axis);
  }

  template<class cont_t>
  int reset_media(cont_t &cont, int offset){
    global_ifermi=cont.global_ifermi;
    if(!SZ)
      return 0;
    medc.map_to(cont.medc, offset, SZ);
    return 1;
  }

  int clear_geometry(){
    medc.clear();
    medh.clear();
    return 1;
  }

  // assign memory indices to nodes in memmap, define internal slices and collect them in vector slices,
  // the rest boundary nodes are collected in b_slices
  int organize_memory();

  int set_local_media(const scRegionMediumMap* mr_map);

    // get node position on a grid
  Vector_3 get_position(int ix, int iy, int iz) const{
    return grd.Position(ix,iy,iz);
  }

  // get node position on a grid
  Vector_3 get_position(const node_t &p) const{
    int ic[3];
    grd.unpack_ind(p.gind,ic[0],ic[1],ic[2]);
    return grd.Position(ic[0],ic[1],ic[2]);
  }


  /// load values from some table to array corresponding to chosen representation
  int load_ini_values(virt_unary_function<const Vector_3 &, vec_type> *table, int repr);

  // returns grid node for given geometry index
  inline void unpack_ind(int ind, int &ix, int &iy, int &iz) const{
    grd.unpack_ind(ind, ix, iy, iz);
  }

  /* for given mesh node with memory index ind and geometry index gind, 
   put memory indices of adjacent nodes to vector n and their amount to vector nsz
   flag:
   first bit (1) - put given node as well to n
   second bit (2) - add offset to all indices
   third bit (4) - do not sort indices (otherwise they will be sorted in the order of their memory indices)
   fourth bit (8) - return -1, if some node is stored in memory but do not belong to mflag filter
   fifth bit (10) - put nodes that are not stored in memory (as -1)
   mflag: bit flag filter (see mesh nodes classification), 
   only nodes which which belong to chosen classification will be put
   */
  int fill_adjacent(int ind,int gind,vector<int> &n,vector<int> &nsz,int flag,int mflag);

  // call fill_adjacent for internal node
  // record its geometry index gind
  int fill_adjacent_internal(int ind,vector<int> &n,vector<int> &nsz,int flag=0,int mflag=mInt,int *gind_=NULL);

  // call fill_adjacent for boundary node
  // record its geometry index gind
  int fill_adjacent_border(int ind,vector<int> &n,vector<int> &nsz,int flag=0,int mflag=mInt,int *gind_=NULL);

  // calculate residue of discretized Poisson equation in node 
  // defined by shift position whithin some slice 
  // record calculated residue to *y
  // return residue squared
  valtype Poisson(const slice_t &slice, int shift, valtype *y, int bf=0, int carsz=0);

  // calculate residue of discretized continuity equation(s) in node
  // defined by shift position whithin some slice 
  // record calculated residue to y[0] and / or y[carsz] 
  // depending on set up bits in bf flag, see documentation to get_bf
  // return residue squared
  valtype Continuity(const slice_t &slice, int shift, valtype *y, int bf, int carsz);

  typedef valtype (RectBlock::*bulk_fun_t)(const slice_t &slice, int shift, valtype *y, int bf, int carsz);

  bulk_fun_t get_bulk_function(int rep){
    if(rep&rpPsi)return &RectBlock::Poisson;
    if(rep&rpFermi)return &RectBlock::Continuity;
    return NULL;
  }

  // record residue to array y for all internal nodes
  // return sum_y^2
  valtype internal_nodes(int repr, valtype *y){
    bulk_fun_t bulk_fun = get_bulk_function(repr);
    int bf = get_bf(repr);
    valtype val=0;
    for(node_it it=begin_int(),e=end_int();it!=e;++it){
      int loc;
      slice_t &slice=it.get_slice_position(loc);
      val+=(this->*bulk_fun)(slice,loc,y+slice.oind+loc,bf,extsz);
    }
    return val;
  }

  // record residue to y for szind nodes from sind
  // skipe those nodes from szind which are not internal
  // if node is not skipped set touch[ind]=true
  void chosen_nodes(int repr, int *sind, bool *touch, int szind, int shift, valtype *y){
    bulk_fun_t bulk_fun = get_bulk_function(repr);
    int bf = get_bf(repr);
    for(chosen_it it=begin_chosen(slices,sind,szind,shift),e=end_chosen(slices);it!=e;++it){
      int loc;
      slice_t &slice=it.get_slice_position(loc);
      int cind=it.get_cind();
      (this->*bulk_fun)(slice,loc,y+cind,bf,2*dim+1);
      touch[cind]=true;
    }
  }

  // update Jacobian
  int Jacobian(const meth_cont_t &mc);

  // update J array with current values
  int update_J();

  // update current value for boundary point
  int update_bJ(int bj);

  /* Record mesh node coordinates and data of chosen type outtype at these nodes to text file
     if outtype=0 then only mesh coordinates will be recorded  
     if dif=false then text file name is msuf.d, where m is mesh name, and suf is argument of this function
     if dif=true then different files will be created for different nodes:
       msuf.d: for internal mesh nodes,
       mxsuf.d: for boundary mesh nodes, where x is c for contact, d for dielectric, r for radial axis
       mtsuf.d: for transfer mesh nodes
  */
  int DumpMesh(const string &suf, bool dif, int outtype);

  // mesh node and its adjacent nodes
  class node_t{
    friend class mesh_it;
    friend class node_it;
    friend class RectBlock;

    Vector_3 center; // node position
    valtype v; //control volume, can be used for integration of some variables (f. e. generation rate)
    int ind; // memory index (or grid index if node_t is returned by mesh_it::operator *)
    int gind; // grid index

    Vector_3 adj[2*dim_max]; // adjacent nodes
    // memory indices (or grid indices if node_t is returned by mesh_it::operator *) of adjacent nodes
    int adj_ind[2*dim_max];

    Vector_3 cur[2*dim_max]; // mid-step adhacent nodes (for current calculation)

    int dim; // mesh dimension (used in _end functions)

  public:

    Vector_3 GetCenter()const{return center; }
    valtype GetControlVolume()const{return v;}
    Vector_3 *adj_begin()const{return (Vector_3 *)adj; }
    Vector_3 *adj_end()const{return (Vector_3 *)(adj+2*dim); }
    Vector_3 *cur_begin()const{return (Vector_3 *)cur; }
    Vector_3 *cur_end()const{return (Vector_3 *)(cur+2*dim); }
    int get_ind()const{return ind; }
  };

  // make object node_t for node specified by its position ind_shift inside some slice
  // when make_node is called from mesh_it::operator *, slice.oind = -1, slice.iind = NULL
  // since at this moment memory indices are not assigned.
  // in this case we assume that all nodes are internal, and slices are set correspondingly
  node_t make_node(const slice_t &slice, int ind_shift=0);

  // mesh nodes iterating within some slices (used in function scBlockContainer::analyze)
  class node_it{
    friend class RectBlock;

  protected:
    RectBlock *parent;
   
    vector<slice_t>::iterator slice_it; // current slice
    int i; // position within a slice

  public:

    node_it(RectBlock *sparent, vector<slice_t> &slices, int e=0):parent(sparent){
      slice_it = e ? slices.end() : slice_it=slices.begin();
      i = e ? -1 : 0;
    }

    // should be calles only for boundary slices
    // return slice number
    int bound_index() const{
      return slice_it-parent->b_slices.begin();
    }


    node_t operator*(){
      return parent->make_node(*slice_it,i);
    }

    // return current slice and position within a slice
    slice_t &get_slice_position(int &shift){
      shift=i;
      return *slice_it;
    }

    node_it& operator++(){ //prefix
      i++;
      if(i>=slice_it->sz){ // go to the next slice
        ++slice_it;
        i=0;
      }
      return *this;
    }
     
    node_it operator++(int){ // postfix
      node_it tmp=*this;
      ++*this;
      return tmp;
    }

    /// tests for the end only
    bool operator!=(const node_it &other) const{
      return slice_it != other.slice_it;
    }

  };

  node_it begin_int(){
    return node_it(this,slices);
  }

  node_it end_int(){
    return node_it(this,slices,1);
  }

  node_it begin_bound(){
    return node_it(this,b_slices);
  }

  node_it end_bound(){
    return node_it(this,b_slices,1);
  }


  // iterator on chosen nodes which are filtered if they do not belong to some slices
  class chosen_it{
    friend class RectBlock;

  protected:
    RectBlock *parent;
   
    vector<slice_t>::iterator slice_it, slice_e; // iterator on current and end slices
    int *sind; // chosen nodes
    int szind; // size of sind
    int cind; // current element of sind
    int shift; // if shift!=0, current node is shifted: sind[cind]+shift

  public:

    chosen_it(RectBlock *sparent, vector<slice_t> &slices, int *sind_, int szind_, int shift_, int e=0):
    parent(sparent), sind(sind_), szind(szind_), shift(shift_){
      cind = e ? szind : -1;
      slice_e = slices.end();
      slice_it = e ? slice_e : lower_bound(slices.begin(),slices.end(),slice_t(sind[0]+shift));
      if(!e){
        if(slice_it!=slices.begin())
          slice_it--;
        operator++();
      }
    }

    node_t operator*(){
      return parent->make_node(*slice_it,sind[cind]+shift);
    }

    slice_t &get_slice_position(int &pos){
      pos = (sind[cind]+shift)-slice_it->oind;
      return *slice_it;
    }

    int get_cind(){
      return cind;
    }

    chosen_it& operator++(){ //prefix

      cind++;
      if(cind>=szind){
        slice_it=slice_e;
        return *this;
      }

      while(((sind[cind]+shift) < slice_it->oind) || ((sind[cind]+shift) >= slice_it->sz + slice_it->oind)){
        // execute this loop untill we get into some slice

        while((sind[cind]+shift) < slice_it->oind){ // go to the next cind
          cind++;
          if(cind>=szind){ // no more chosen nodes
            slice_it=slice_e;
            return *this;
          }
        }

        while((sind[cind]+shift) >= slice_it->sz + slice_it->oind){ // go to the next slice
          ++slice_it;
          if(!(slice_it!=slice_e)) // no more slices
            return *this;
        }
      }
      return *this;
    }
     
    chosen_it operator++(int){ // postfix
      chosen_it tmp=*this;
      ++*this;
      return tmp;
    }

    /// tests for the end only
    bool operator!=(const chosen_it &other) const{
      return (slice_it != other.slice_it);
    }

  };

  chosen_it begin_chosen(vector<slice_t> &slices, int *sind, int szind, int shift){
    return chosen_it(this,slices,sind,szind,shift);
  }

  chosen_it end_chosen(vector<slice_t> &slices){
    return chosen_it(this,slices,NULL,0,0,1);
  }


  // iterator on all mesh points (including unused ones), used in scBlockContainer::specify_nodes only
  class mesh_it{
    friend class RectBlock;
    typedef grid_t::iterator grid_it;
  protected:
    RectBlock *parent;
    grid_it it;

    mesh_it(RectBlock *sparent):parent(sparent){
      it=parent->grd.begin();
    }

  public:

    /// default constructor, constructs sequence end
    mesh_it(){}

    // total number of adjacent nodes (some of them might be absent if this is border node)
    int get_adjacent_num()const{
      return 2*parent->dim;
    }

    // return iterator i-th adjacent node
    // if the case of bondary node, some adjacent nodes could be absent, in this case return end iterator
    mesh_it get_adjacent(int i)const{
      int sh[3]={0,0,0};
      int di;
      int ni;
      parent->dnum(i,di,ni);
      sh[di] = ni ? 1 : -1;
      int ic=it.get_ind(di); // index in shifted dimension
      int exc=ic+sh[di]; // shift this index
      if(exc<0 || exc>=parent->N[di])
        return mesh_it();
      mesh_it bpit=*this;
      bpit.it=bpit.it.shift(sh);
      return bpit;
    };

    node_t operator*(){
      int ind[3];
      for(int i=0;i<3;i++)
        ind[i]=it.get_ind(i);
      return parent->make_node(slice_t(-1,parent->grd.pack_ind(ind[0],ind[1],ind[2])));
    }

    mesh_it& operator++(){ //prefix
      ++it;
      return *this;
    }
     
    mesh_it operator++(int){ // postfix
      mesh_it tmp=*this;
      ++*this;
      return tmp;
    }

    /// tests for the end only
    bool operator!=(const mesh_it &other) const{
      return it!=other.it;
    }
  };

  mesh_it begin(){
    return mesh_it(this);
  }

  mesh_it end(){
    return mesh_it();
  }
};

#endif
