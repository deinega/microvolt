#ifndef SC_UI_H
#define SC_UI_H

/** @file sc_ui.h
This is user interface to MICROVOLT
*/

#include <memory>

#include "region_3.h"
#include "cpp11features.h"
#include "sc_carrier_params.h"
#include "sc_mediumfactory.h"
#include "sc_alloyfactory.h"
#include "sc_generation.h"

class scAbstractMedium;
struct scBoundaryRegion;

/// class that connects user interface and MICROVOLT implementation
struct sc_storage;

/// MICROVOLT initialization
/// argc, argv - arguments of function main (command line parameters)
/// outdir - output directory
int scInit(int argc, char **argv, string outdir="");

/// this bit flag is used to specify detectors type
enum DET_TYPES{
  DET_VAR=0,  /// detectors for data in chosen point (potential, concentrations, fermi-levels etc.)
  DET_FLUX=1, /// detectors for current flux through chosen surface
};

/** User interface class to specify numerical experiments
It manages MICROVOLT container (scBlockContainer) and detecors objects which implement numerical model
*/
class scExperiment: public apComponent{

protected:

  sc_storage *store;

  /// Parameters below dublicate parameters in MICROVOLT container with the same name
  /// They are directly transfered to scBlockContainer after intialization

  int dim_type[3];

  int rad;
  valtype rad_axis;
  bool rad_right;

  int per_flags;

  valtype clamp;

  valtype residue_threshhold;

  valtype stop_ratio, small_stop_ratio;

  valtype gen_eq;

  Vector_3 depl_decr;

  valtype n0;
  valtype lu;
  valtype T;

  Vector_3 zero_pos;

  Vector_3 grad_split;

  valtype sys_t0; /// program start time

  int tid; /// main timer id

  /// execution phases
  enum EXE_PHASES {
    PH_GEOMETRY = 0x1, 
    PH_CALC = 0x2
  }; 

  int phase; /// current phase
  int all_phases; /// phases to be executed

  Vector_3 ip1,ip2; /// calculated box size

  /// if residue of converged solution is bigger than min_residue, then this solution is considered as incorrect.
  /// in the case of solar cell simulation, calculated values for Isc (Voc) 
  /// will not be presented in corresponding output files (see documentation to Calculate)
  valtype min_residue;

  string flux_det; /// detector used to calculate Isc, Voc and efficiency of solar cell
  /// area of the surface of this detector. 
  /// Isc is calculated as a flux through this surface divided on its area
  valtype flux_norm;

  valtype VStart; /// starting value of the voltage

  /// test flags to specify detectors and mesh output
  int dump_iter; /// detectors (and meshes) will produce output for each calculation step
  /// bit flag regulating meshes dumping (used for test purposes), first bit - dump meshes, 
  /// second bit - create different files for different mesh nodes (internal, boundary, transfer)
  int dm;

  /// true if the initial values are estimated using linear extrapolation of
  /// distributions potential and QF levels obtained at previous voltage levels
  bool useLinearEstimate;

  /// Noise level for extrapolated distributions
  /// Exact linear extrapolation sometimes produces an estimate that isn't good
  /// as a starting point (convergence may be poor). Slightly perturbed values
  /// may be better
  valtype estimationNoise;

  /// mesh descriptor
  struct block_t{
    Vector_3 dr; /// mesh step along each dimension
    iVector_3 sz; /// mesh steps numbers along each dimension
    /// mesh grids coordinates along each dimension
    /// if axis[i]==NULL, then mesh coordinates are defined from sz or dr
    mngptr<valtype> axis[3];
    Region_3 *confinement; /// if not NULL, mesh is confined within this region
    block_t(const iVector_3 &sz_,Region_3 *conf):sz(sz_),confinement(conf),dr(0){}
    block_t(const Vector_3 &dr_,Region_3 *conf):dr(dr_),confinement(conf),sz(0){}
  };

  /// boundary regions
  typedef scBoundaryRegion bond_t;

  /// generation rate profiles
  typedef scGenerationRegion gen_t;

  /// detectors descriptor
  struct det_t{
    int flag; /// detector type, see DET_TYPES

    string name; /// detector name (used in output file names)
    int gr; /// vector set type: 0 - UniformGrid, 1 - BoxSurfaceSet, 2 - CylinderSurfaceSet, 3 - vector<Vector_3>
    int ind; /// index of vector set in corresponding refvector from storage

    int output; /// output type (see scoutTYPE)
    string format; /// double values format used in text output files

    Region_3 *conf; /// if conf!=NULL, detectors will be confined within this region 

    det_t(const string &name_,int gr_,int ind_,int flag_=0,Region_3 *conf_=NULL):
      name(name_),gr(gr_),ind(ind_),flag(flag_),conf(conf_),output(0xffff&(~outJtot)),format("%g"){}
  };

  /// vectors of specified objects, second element of pair specifies priority level
  vector<pair<block_t,int> > blocks; /// meshes

  vector<pair<bond_t *,int > > bonds; /// boundary regions
  vector<pair<gen_t *,int > > gens; /// generation rate profiles
  refmap<string,det_t> dets; // detectors

  typedef std::map<std::string,virt_unary_function<const Vector_3 &, vec_type>* > ini_funcmap_t;
  // Map that stores functions loading initial values at Prepare stage of calculation. The map keys are repesentations.
  ini_funcmap_t ini_funcmap;
  /// initial guess value for fermi levels
  mngptr<virt_unary_function<const Vector_3 &, vec_type> > fermi_guess[2];
  
  // used to load initial values from file
//  string ini_name;
//  int ini_s, ini_sa;
//  Vector_3 ini_sym_center;

  /// create container, add meshes, media, boundaries, generation rate profile
  int Init();

  /// allocate memory and organize memory layout, calculate and assign coefficients in discretized equations,
  /// create detectors and allocate their memory buffers
  int Prepare();

public:

  /// Sets the position of calculation space 
  /// if some for some dimension p1[i]==p2[i], then this dimension is excluded
  scExperiment(const Vector_3 &p1,const Vector_3 &p2);

  ~scExperiment();

  /// specify executed phases: "g" - geometry initialization only (used for test purposes),
  /// "c" - initialization and calculation
  void SetPhases(const string &s){
    if(s=="g")
      all_phases=PH_GEOMETRY;
    else if(s=="c")
      all_phases=PH_GEOMETRY|PH_CALC;
  }

  void SetClamp(valtype cl){
    clamp=cl;
  }

  void set_convergence_check(valtype ratio=0.05, valtype threshold=-1, valtype small_ratio=0.5){
    stop_ratio=ratio;
    residue_threshhold = threshold;
    small_stop_ratio = small_ratio;
  }


  void SetMinResidue(valtype res){
    min_residue=res;
  }


  void SetGenEq(valtype gen_eq_){
    gen_eq=gen_eq_;
  }

  void SetDepletionDecreasing(Vector_3 dd){
    depl_decr=dd;
  }

  void SetTemperature(valtype T_=300){
    T=T_;
  }

  void SetConcentrationUnit(valtype n0_=1e10){
    n0=n0_;
  }

  void SetConcentrationUnit(const scMediumFactory &med);

  void SetLengthUnit(valtype lu_){
    lu=lu_;
  }

  void SetZeroPosition(Vector_3 zero_pos_){
    zero_pos=zero_pos_;
  }

  void SetDirSplit(valtype cx=1., valtype cy=1., valtype cz=1.){
    grad_split=Vector_3(cx,cy,cz);
  }

  void DumpEachIteration(){
    dump_iter=1;
  }

  void DumpMeshes(int dm_=1){
    dm=dm_;
  }

  void SetStartVoltage(valtype vStart) {
    VStart = vStart;
  }

  void EnableLinearEstimate(valtype noiseLevel = 0) {
    useLinearEstimate = true;
    estimationNoise = noiseLevel;
  }

  /// Set the numbers of mesh steps for each direction for default mesh 
  int SetResolutionN(iVector_3 sz){
    if(blocks.size())
      return message(vblMESS1,-1,"SetResolutionN: default mesh is already defined\n");
    int dim=0;
    for(int i=0;i<3;i++){
      if(sz[i]<0)
        return message(vblMESS1,-1,"Mesh step number is negative in %d direction...\n",i);
      if(sz[i])
        dim++;
    }
    if(!dim)
      return message(vblMESS1,-1,"Zero mesh step number along all directions...\n");
    AddMeshN(sz);
    return 1;
  }

  /// Set mesh step for each direction for default mesh 
  int SetResolution(Vector_3 dr){
    if(blocks.size())
      return message(vblMESS1,-1,"SetResolution: default mesh is already defined\n");
    int dim=0;
    for(int i=0;i<3;i++){
      if(dr[i]<0)
        return message(vblMESS1,-1,"Mesh step is negative in %d direction...\n",i);
      if(dr[i])
        dim++;
    }
    if(!dim)
      return message(vblMESS1,-1,"Zero mesh step number along all directions...\n");
    AddMesh(dr);
    return 1;
  }

  /// Add mesh of some level with some numbers of mesh steps sz for each direction 
  /// if conf!=NULL, mesh will be confined within this region
  void AddMeshN(iVector_3 sz, Region_3 *conf=NULL, int level=0){
    blocks.push_back(make_pair(block_t(sz,conf),level));
  }

  /// Add mesh of some level with some mesh step dr for each direction 
  /// if conf!=NULL, mesh will be confined within this region
  void AddMesh(Vector_3 dr, Region_3 *conf=NULL, int level=0){
    blocks.push_back(make_pair(block_t(dr,conf),level));
  }

  /// define nonuniform grid for di-th dimension of i-th mesh (if i==-1, then last mesh)
  /// nodes coordinates in this dimension are stored in some container cont_t (could be vector)
  template<class cont_t>
  int SetMeshSteps(int di, cont_t cont, int i=-1){
    if(!blocks.size())
      return message(vblMESS1,-1,"SetMeshSteps: default mesh is not defined\n");
    if(i>=(int)blocks.size())
      return message(vblMESS1,-1,"SetMeshSteps: there is no such a mesh (check the last argument)\n");

    if(i<0)
      i=blocks.size()-1;

    int &n=blocks[i].first.sz[di];
    n=cont.size();
    if(!n)
      return 0;

    blocks[i].first.axis[di].reset(new valtype[n],0x8);
    valtype *ptr = blocks[i].first.axis[di].first;

    for(int ii=0;ii<n;ii++)
      ptr[ii]=cont[ii];

    return 1;
  }

  /// define nonuniform grid for di-th dimension of i-th mesh (if i==-1, then last specified mesh)
  /// beg, end - iterator on nodes coordinates in this dimension
  template<class inp_it>
  int SetMeshSteps(int di, inp_it beg, inp_it end, int i=-1){
    if(!blocks.size())
      return message(vblMESS1,-1,"SetMeshSteps: default mesh is not defined\n");
    if(i>=(int)blocks.size())
      return message(vblMESS1,-1,"SetMeshSteps: there is no such a mesh (check the last argument)\n");

    if(i<0)
      i=blocks.size()-1;

    int &n=blocks[i].first.sz[di];
    n=0;
    for(inp_it it=beg;it!=end;++it)
      n++;
    if(!n)
      return 0;

    blocks[i].first.axis[di].reset(new valtype[n],0x8);
    valtype *ptr = blocks[i].first.axis[di].first;

    int ii=0;
    for(inp_it it=beg;it!=end;++it)
      ptr[ii++]=*it;

    return 1;
  }

  /// set confining region for i-th mesh (if i==-1, then last specified mesh)
  int SetConfinement(Region_3 *conf,int i=-1){
    if(!blocks.size())
      return message(vblMESS1,-1,"SetConfinement: default mesh is not defined\n");
    if(i>=(int)blocks.size())
      return message(vblMESS1,-1,"SetConfinement: there is no such a mesh (check the last argument)\n");

    if(i<0)
      i=blocks.size()-1;
    blocks[i].first.confinement=conf;

    return 1;
  }

  /// set radial coordinate, position of radial axis (r=0)
  /// if rad_right = true (false) then positive (negative) values of radial coordinate are considered
  void SetRadialCoordinate(int rad_, valtype rad_axis_, bool rad_right_){
    rad=rad_;
    rad_axis=rad_axis_;
    rad_right=rad_right_;
  }

  // set bit flag indicating periodicity in specific directions: 0x1 for x, 0x2 for y, 0x4 for z
  void SetPeriodicBoundaries(int per_flags_){
    per_flags=per_flags_;
  }

  /**
  see documentation to method_t and meth_chain_t
  
  set of Jacobians to be solved in sequence
  Jacobians are coded in string format (see documentation to method_t)
  sequence of Jacobians are defined eq as a sequence of these coding strings separated by commas, f.e.
  p,eh correspond to Jacobian for Poisson equation and Jacobian for both continuity equations

  an the each calculation step specified Jacobians will be solved n times, 
  where n is corresponding number in string num (or 1 if this number is not specified), 
  f.e. eq="eh,p" num="2,1" means that Jacobian for continuity equations will be solved twice 
  and Jaobian for Poisson equation will be solved once.
  this is the same as "eh,p" "2", since unspecified number corresponds to one

  sequence of specified Jacobians will be used until solution converged to some (local) minimum, 
  but not more than tot_num times (if tot_num>0)

  snum_gs, it_gs are gauss seidel method parameters */
  void SetMethod(const string &eq, const string &num="", int tot_num=50, int snum_gs=1,int it_gs=1);

#ifdef HAS_SMART_PTR
  /// add a medium region
  /// level is region priority
  void AddCreatedMediumRegion(const std::shared_ptr<const scAbstractMedium>& med,
      const std::shared_ptr<const Region_3>& reg = std::shared_ptr<const Region_3>(), int level = 0);

  /// add a medium region
  /// level is region priority
  /// legacy version with a raw pointer for the region
  void AddCreatedMediumRegion(const std::shared_ptr<const scAbstractMedium>& med,
      const Region_3 *reg, int level = 0);
#else
  /// add a medium region
  /// level is region priority
  // auto_ptr_ref is needed instead of auto_ptr in order to allow
  // calls like AddMediumRegion(medFactory.create(), reg, level);
  void AddCreatedMediumRegion(std::auto_ptr_ref<const scAbstractMedium> med,
    const Region_3 *reg = 0, int level = 0);

#endif

  template<class factory_t>
  void AddMediumRegion(const factory_t& med, const Region_3 *reg, int level = 0){
    AddCreatedMediumRegion(med.create(),reg,level);
  }

  /// delete all specified medium regions
  void ClearMediumRegions();

  /// initialize container and detectors
  /// ini_repr is flag for variables which will be reinitialized with default initial values
  /// other variables will be kept as they are (from previous solution obtained after function Calculate)
  int BuildMediumRegions(int ini_repr);

  /// add Ohmic contact
  /// level specifies corresponding boundary region priority
  void AddContact(Region_3 *reg, int level=0);

  /// add Shottky contact with specified electron and hole surface recombination (sn, sp), cm/s
  /// and contact workfunction metal_workf, eV (default VEC_INFTY sets simple Ohmic contact)
  void AddContact(Region_3 *reg, valtype sn, valtype sp, valtype metal_workf= VEC_INFTY,int level=0);

  /// see AddContact
  /// voltage will be gradually applied to this contact
  void AddVoltageContact(Region_3 *reg, int level=0);

  /// see AddContact
  /// voltage will be gradually applied to this contact
  void AddVoltageContact(Region_3 *reg, valtype sn, valtype sp, valtype metal_workf= VEC_INFTY,int level=0);

  /// add reflective boundary of some level
  void AddReflectiveBoundary(Region_3 *reg, int level=0);

  /// add insulator interface with specified electron and hole surface recombination (sn, sp), cm/s
  /// level specifies corresponding boundary region priority
  void AddReflectiveBoundary(Region_3 *reg, valtype sn, valtype sp, int level=0);

  /// add generation profile 
  /// if reg!=NULL then region is specified within some region
  /// level specifies gereation region priority
  void AddGeneration(scGeneration *gen_,Region_3 *reg=NULL,int level=0){
    gens.push_back(make_pair(new gen_t(gen_,reg),level));
  }
  
  /// Set loader function for initial values for given representation. The func pointer is unmanaged.
  /// \a type is a string specifying function type. Avaulable types:\n
  /// "potential" -- sets initial potential;
  /// "fermi_n" -- sets initial fermi level for electrons
  /// "fermi_p" -- sets initial fermi level for holes
  /// The values are loaded by record_internal() at Prepare stage.
  void SetIniValueFunction(const string& type, virt_unary_function<const Vector_3 &, vec_type> *func);

  /// adding initial guess value for fermi levels
  void AddFermiGuess(virt_unary_function<const Vector_3 &, vec_type> *g0, virt_unary_function<const Vector_3 &, vec_type> *g1){
    fermi_guess[0].reset(g0,1);
    fermi_guess[1].reset(g1,1);
    SetIniValueFunction("fermi_n",g0);
    SetIniValueFunction("fermi_p",g0);
  }

  /// add a 3D detector grid set (p1, p2 are opposite grid vertices, 
  /// n is amount of grid steps for each dimension).
  /// The name is used in filename of output files.
  /// If n[i]=INT_INFTY, then detector grid resolution for this direction will be the same as default mesh resolution
  /// out is a bit flag which defines output data type (see scoutTYPE)
  int AddDetector(const string &name, Vector_3 p1, Vector_3 p2, iVector_3 n,
    int out=0xffff^(outJtot|outResidue));

  /// add a 3D detector grid set (p1, p2 are opposite grid vertices, 
  /// rsl is resolution for each dimension corresponding to main mesh resolution).
  /// The name is used in filename of output files.
  /// If n[i]=INT_INFTY, then detector grid resolution for this direction will be the same as default mesh resolution
  /// out is a bit flag which defines output data type (see scoutTYPE)
  int AddDetector(const string &name, Vector_3 p1, Vector_3 p2, Vector_3 rsl,
    int out=0xffff^(outJtot|outResidue));

  /// add box for calculation of current flux through its chosen sides (see BOX_SIDES)
  /// (p1, p2 are opposite box vertices, n is amount of grid steps for each dimension).
  /// current flux is calculated as integral{J * dS}, where J is measured in A/cm2, 
  /// and dS is measured in user length units squared
  int AddBoxDetector(const string &name, Vector_3 p1, Vector_3 p2, iVector_3 n, int sides);

  /// add plane for calculation of corresponding current flux
  /// plane is perpendicular to direction dir and passes through pos coordinate corresponding to this direction
  /// dr defines spacing between detectors at the plane (if dr=0 then this spacing will be close to default mesh step)
  void AddFluxPlane(const string &name, size_t dir, valtype pos, Vector_3 dr=0);

  /// add cylindrical grid set of detectors
  /// cylindrical grid is defined by cylinder origin, height direction norm, radius R,
  /// number of mesh steps n along radial and axis directions
  /// surf_dir defines direction of the normal to cylindrical surface (inside or outside)
  /// which is used to define the sign of the current flux
  /// flag is detector type (see DET_TYPES)
  void AddCylinderDetector(const string &name, const Vector_3 &origin, const Vector_3 &norm,
    valtype R, const iVector_2 &n, int surf_dir, int flag=0);

  /// specifies detector to calculate Isc in the case of solar cell simulation.
  /// if area!=0, specifies the area used for Isc normalization 
  /// (otherwise, it would be calculated automatically)
  void SetFluxIV(const string &name, valtype area=0){
    flux_det=name;
    flux_norm=area;
  }

//  void SetIniValue(const string &name, int s,int sa, const Vector_3 &sym_center);
//  int IniValue(const string &name, int s,int sa, const Vector_3 &sym_center);

  int FillIniValues(int ini_repr);

  /// create and initialize container and detectors, record dump files
  int BuildGeometry(){
    if(!(all_phases&PH_GEOMETRY))
      return 0;
    if(Init()<0)
      return -1;
    if(Prepare()<0)
      return -1;
    return 1;
  }

  /**
  Make a calculation for set of applied voltages
  if dV=0 then only zero applied voltage will be calculated
   if dV<>0 then:
     applied voltage will be increased by step dV, and calculation will be performed 
       until voltage will get value V
     if there is generation profile is specified and SetFluxIV is called, voltage will be increased
       until current value for flux_det detector will change sign.
       At the end of calculation you could find the file sc.dat with calculated Isc, Voc and efficiency of solar cell
       which are measured using this detector specified ib SetFluxIV.
       You alse can find file sc_line.dat with the same data printed in one line
       (if you have optimzation loop, different files sc_line.dat could be collected 
       in a table format to plot the results)
  */
  int Calculate(valtype dV=0, valtype V=0);

  /**
  Make a calculation for a set of applied voltages
  If V is empty, then only zero applied voltage will be simulated (or is it better to do nothing in this case?),
  otherwise the calculation will be run for the specified voltages. 

  If a generation profile is specified and SetFluxIV is called, then the voltage will be increased
    until current value for flux_det detector changes sign. If the sign was the same for all specified
    voltages and V contains 2 or more elements, then the calculation will also be performed for voltages
    V.back() + (V.back() - V[V.size()-2])*n, where n = 1,2,... until the current changes sign.
  */
  int Calculate(std::vector<valtype> V);

  // Dump ath their current state meshes using filename suffix
  int DoDumpMeshes(const std::string &suffix);

private:
  // Dump the current state of the calculation
  // Returns 1 on success, -1 if an error occured during dumping
  int dumpCurrentStep(int nIteration, valtype voltage);
};

#endif
