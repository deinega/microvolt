#include "sc_ui.h"
#include "sc_mesh.h"
#include "sc_fd.h"
#include "sc_dump.h"
#include "detector.h"
#include "string_utils.h"

#include "sc_medium_hqs.h"

// this is storage for all MICROVOLT objects which are hidden from user interfalce declared in sc_ui.h
struct sc_storage{

  mngptr<scBlockContainer<RectBlock> > cont;
  vector<BaseDetectorSet *> dets;

  refvector<UniformGrid<Vector_3> > gsets;
  refvector<BoxSurfaceSet> bsets;
  refvector<CylinderSurfaceSet> csets;
  refvector<vector<Vector_3> > vsets;

  vector<meth_chain_t> methods;

  scRegionMediumMap meds;
};

string scInfo(){
  string info;
  info+="Microvolt, (c) University of Toronto\n";
# ifdef USE_MPI
  info+="Parallel version. ";
# else
  info+="Serial version. ";
# endif
# ifndef SINGLE_PRECISION
  info+="Double precision. ";
# else
  info+="Single precision. ";
# endif
  info+="\n";
  return info;
}

int scInit(int argc, char **argv, string outdir){
  // substituting info function
  theConfig->info=&scInfo;
  // registering dumpers for components
  scAdd2Dump();
  theDump->SetFormat(GNUPLOT);
  return apInit(argc,argv,outdir,outdir);
}

scExperiment::scExperiment(const Vector_3 &p1,const Vector_3 &p2): 
  rad(-1),rad_axis(0),rad_right(true),per_flags(0),
  clamp(1./5),gen_eq(0),/*n0(9978706067256290),*/lu(1),T(300),grad_split(1.,1.,1.),
  phase(0),all_phases(PH_GEOMETRY|PH_CALC),
  min_residue(VEC_INFTY),flux_norm(0),VStart(0),dump_iter(0),dm(0),
  useLinearEstimate(false), estimationNoise(0) {

  ut=1;
  dump=1;

  set_convergence_check();

  ip1=p1;
  ip2=p2;

  SetConcentrationUnit(getscSi());

  store=new sc_storage;

  sys_t0=apTimer::gettime(false); // start time
  tid=AddTimersField("all",-1);
}

scExperiment::~scExperiment(){

  if(ut) // printing timers data
    log->TimeUsage(fmt("timers_%03d.d",theComm.get_myrank()),TABLE_TAB);

  for(size_t i=0;i<bonds.size();i++)
    delete bonds[i].first;
  for(size_t i=0;i<gens.size();i++)
    delete gens[i].first;
  for(refmap<string,det_t>::iterator it=dets.begin(),e=dets.end();it!=e;++it)
    delete it->second->conf;
  for(vector<BaseDetectorSet *>::iterator it = store->dets.begin(); it != store->dets.end(); it++)
    delete *it;
  delete store;

  valtype sec=apTimer::gettime(false)-sys_t0; // calculating simulation time
  string st=stime(int(sec));
  message(vblMESS1,0,"Elapsed time is %s\n",st.c_str());
}

void scExperiment::SetConcentrationUnit(const scMediumFactory &med){
  n0=med.createSpecific()->getLocalMedium()->get_ni();
}

int scExperiment::AddDetector(const string &name, Vector_3 p1, Vector_3 p2, iVector_3 n, int out){

  for(int i=0;i<3;i++){
    if(n[i]==INT_INFTY && ip1[i]!=ip2[i]){
      if(!blocks.size())
        return message(vblMESS1,-1,"AddDetector: default mesh resolution is not defined, change the third parameter to finite value\n");
      valtype dr=blocks[0].first.dr[i];
      if(!dr)
        dr = acdiv(ip2[i]-ip1[i],valtype(blocks[0].first.sz[i]));
      n[i]=int(ceil(acdiv(p2[i]-p1[i],dr)))+1;
    }
    if(ip1[i]==ip2[i] || n[i]<1){
      n[i]=1;
    }
  }

  store->gsets.push_back(new UniformGrid<Vector_3>(p1,p2,(int *)&n));
  det_t *d = new det_t(name,0,store->gsets.size()-1);
  d->output=out;
  dets[name]=d;

  return 1;
}

int scExperiment::AddDetector(const string &name, Vector_3 p1, Vector_3 p2, Vector_3 rsl, int out){
  if(!blocks.size())
    return message(vblMESS1,-1,"AddDetector: default mesh resolution is not defined\n");
  iVector_3 n;
  for(int i=0;i<3;i++){
    if(ip1[i]!=ip2[i]){
      valtype dr=blocks[0].first.dr[i];
      if(!dr)
        dr = acdiv(ip2[i]-ip1[i],valtype(blocks[0].first.sz[i]));
      n[i]=int(rsl[i]*ceil(acdiv(p2[i]-p1[i],dr)))+1;
    }
  }
  return AddDetector(name,p1,p2,n,out);
}

int scExperiment::AddBoxDetector(const string &name, Vector_3 p1, Vector_3 p2, iVector_3 n, int sides){

  int sd=BOX_BACK_X|BOX_FRONT_X;
  for(int i=0;i<3;i++,sd<<=1){
    if(n[i]==INT_INFTY && ip1[i]!=ip2[i]){
      if(!blocks.size())
        return message(vblMESS1,-1,"AddDetector: default mesh resolution is not defined, change the third parameter to finite value\n");
      valtype dr=blocks[0].first.dr[i];
      if(!dr)
        dr = (ceil(acdiv(ip2[i]-ip1[i],valtype(blocks[0].first.sz[i]))));
      n[i]=int(ceil(acdiv(p2[i]-p1[i],dr)));
    }
    if(n[i]<1)
      n[i]=1;
    if(ip1[i]==ip2[i] && !(sides&sd)){
      p1[i]=ip1[i]-.5, p2[i]=ip1[i]+.5;
      n[i]=1;
    }
  }

  store->bsets.push_back(new BoxSurfaceSet(p1,p2,(int *)&n,sides));
  dets[name] = new det_t(name,1,store->bsets.size()-1,DET_FLUX);

  return 1;
}

void scExperiment::AddFluxPlane(const string &name,size_t dir,valtype pos,Vector_3 dr){
  Vector_3 p1=ip1, p2=ip2;
  p1[dir]=p2[dir]=pos;
  iVector_3 n;
  n[dir]=1;
  for(int j=1;j<3;j++){
    int di=(dir+j)%3;
    if(ip1[di]==ip2[di])
      n[di]=1;
    else if(dr[di])
      n[di]=int(ceil(acdiv(ip2[di]-ip1[di],dr[di])));
    else{
      if(blocks[0].first.sz[di])
        n[di]=blocks[0].first.sz[di];
      else
        n[di]=int(ceil(acdiv(p2[di]-p1[di],blocks[0].first.dr[di])));
    }
  }
  int sides=BOX_FRONT_X<<dir;
  AddBoxDetector(name,p1,p2,n,sides);
}

void scExperiment::AddCylinderDetector(const string &name,const Vector_3 &origin,const Vector_3 &norm,valtype R,const iVector_2 &n, int surf_dir, int flag){
  Vector_3 x,y;
  build_orth_basis(norm,x,y);
  x=Vector_3(0,1,0);
  y=Vector_3(0,0,1);
  store->csets.push_back(new CylinderSurfaceSet(origin,norm,x,y,R,(int *)&n,surf_dir));
  dets[name] = new det_t(name,2,store->csets.size()-1,flag);
}

void scExperiment::SetMethod(const string &eq, const string &num, int tot_num, int snum_gs,int it_gs){
  vector<string> eqs=separate_list(eq,"_,",1);
  vector<string> nums=separate_list(num,"_,",1);

  meth_chain_t ch;
  for(size_t i=0;i<eqs.size();i++){
    int n=1;
    if(nums.size()>i){
      int n2=atoi(nums[i].c_str());
      if(n2>0)
        n=n2;
    }
    ch.add(method_t(eqs[i],snum_gs,it_gs),n);
  }
  ch.set_tot_num(tot_num);
  store->methods.push_back(ch);
}

/// add a medium region
/// level is region priority
#ifdef HAS_SMART_PTR
void scExperiment::AddCreatedMediumRegion(const std::shared_ptr<const scAbstractMedium>& med,
    const std::shared_ptr<const Region_3>& reg, int level) {
#else
void scExperiment::AddCreatedMediumRegion(std::auto_ptr_ref<const scAbstractMedium> med,
    const Region_3 *reg, int level) {
#endif
  Vector_3 dd = depl_decr;
  for(int i = 0; i < dd.dimension; i++) {
    if (dim_type[i] < 0)
      dd[i] = 0;
  }

  SHARED_OR_AUTO(const scAbstractMedium) pmed(med);
  if (dd != 0) {
    const scMedium* mmed = dynamic_cast<const scMedium*>(pmed.get());
    if(mmed) {
      pmed = SHARED_OR_AUTO(const scMediumWithHQSurface)(
        new scMediumWithHQSurface(*mmed, &store->meds, depl_decr));
    }
  }

  // null region causes crashes at dumping, creating a default Region_3 is a workaround
  store->meds.add(pmed, (reg ? reg : scRegionMediumMap::RegionPtr(new Region_3)), level);
}

#ifdef HAS_SMART_PTR
void scExperiment::AddCreatedMediumRegion(const std::shared_ptr<const scAbstractMedium>& med,
    const Region_3 *reg, int level) {
  AddCreatedMediumRegion(med, std::shared_ptr<const Region_3>(reg), level);
}
#endif

void scExperiment::ClearMediumRegions(){
  store->meds.clear();
}

int scExperiment::BuildMediumRegions(int ini_repr){
  if(phase<PH_GEOMETRY)
    return 0;
  store->cont->ShareMediumRegions(store->meds);
  store->cont->BuildGeometry(NULL,ini_repr);

  vector<BaseDetectorSet *> &sdets=store->dets;
  for(size_t di=0;di<sdets.size();di++){
    sdets[di]->clear();
    if(sdets[di]->Init()<0)
      return -1;
  }
  return 1;
}

void scExperiment::AddContact(Region_3 *reg, int level){
  bonds.push_back(make_pair(new bond_t(make_mngarg(reg),1,VEC_INFTY,VEC_INFTY,VEC_INFTY),level));
}

void scExperiment::AddContact(Region_3 *reg, valtype sn, valtype sp, valtype metal_workf,int level){
  bonds.push_back(make_pair(new bond_t(make_mngarg(reg),1,sn,sp,metal_workf),level));
}

void scExperiment::AddVoltageContact(Region_3 *reg, int level){
  bonds.push_back(make_pair(new bond_t(make_mngarg(reg),1,VEC_INFTY,VEC_INFTY,VEC_INFTY,1),level));
}

void scExperiment::AddVoltageContact(Region_3 *reg, valtype sn, valtype sp, valtype metal_workf,int level){
  bonds.push_back(make_pair(new bond_t(make_mngarg(reg),1,sn,sp,metal_workf,1),level));
}

void scExperiment::AddReflectiveBoundary(Region_3 *reg, int level){
  bonds.push_back(make_pair(new bond_t(make_mngarg(reg),2,0,0,VEC_INFTY),level));
}

void scExperiment::AddReflectiveBoundary(Region_3 *reg, valtype sn, valtype sp, int level){
  bonds.push_back(make_pair(new bond_t(make_mngarg(reg),2,sn,sp,VEC_INFTY),level));
}

int scExperiment::Init(){

  for(int i=0;i<3;i++){
    if(ip1[i]==ip2[i]){
      dim_type[i]=-1;
      // artificially make nonzero mesh size in this direction, 
      // otherwise scBlockContainer::TestPoint will always return false
      ip1[i]-=1e-10, ip2[i]+=1e-10;
    }
    else{
      if(blocks[0].first.dr==0){
        if(blocks[0].first.sz[i]<=1)
          return message(vblMESS1,-1,"Zero mesh size in %d dimension\n",i);
      }
      else if(blocks[0].first.dr[i]<=0)
        return message(vblMESS1,-1,"Zero mesh step in %d dimension\n",i);
      dim_type[i]=0;
    }
  }
  if(rad>=0){
    dim_type[rad]=1;
  }
  Box box(ip1,ip2);

  for(int i=0;i<3;i++){
    if(dim_type[i]<0)
      per_flags&=~(0x1<<i);
  }

  store->cont.reset(new scBlockContainer<RectBlock>(dim_type,T),1);
  store->cont->SetLimitingBox(make_mngarg(new Box(box)));
  store->cont->InitComponent("",0xffff,1,makevec<int>(1,tid));
  store->cont->set_convergence_check(stop_ratio, residue_threshhold, small_stop_ratio);
  store->cont->AddFermiGuess(fermi_guess[0].ptr(),fermi_guess[1].ptr());
  store->cont->set_clamp(clamp);
  store->cont->set_gen_eq(gen_eq);
  store->cont->set_depl_decr(depl_decr);
  store->cont->SetZeroPosition(zero_pos);
  store->cont->EnableLinearEstimate(useLinearEstimate, estimationNoise);

  for(size_t i=0;i<blocks.size();i++){
    RectBlock *block;

    Vector_3 p1=ip1,p2=ip2;
    iVector_3 csz;

    for(int di=0;di<3;di++){
      if(dim_type[di]>=0){
        if(blocks[i].first.sz[di])
          csz[di]=blocks[i].first.sz[di];
        else
          csz[di]=int(ceil(acdiv(ip2[di]-ip1[di],blocks[i].first.dr[di])));
      }
      else
        csz[di]=1;

      if(per_flags&(0x1<<di)){
        if(di==rad)
          return message(vblMESS1,-1,"Radial direction is assigned as periodic\n");
        // extending uniform mesh outside periodic boundaries on a one step
        valtype dr=(ip2[di]-ip1[di])/csz[di];
        csz[di]+=2;
        p1[di]-=dr;
        p2[di]+=dr;
      }

      if(di==rad){ // shortenning the mesh box 
        valtype old_width = p2[rad]-p1[rad];

        if(rad_right) p1[rad]=rad_axis;
        else p2[rad]=rad_axis;
        valtype new_width = p2[rad]-p1[rad];

        csz[rad]=int(csz[rad]*new_width/old_width + .5);
      }
    }

    block = new RectBlock(Box(p1,p2),(int *)&csz);

    if(blocks[i].first.confinement)
      block->SetConfinement(make_mngarg(blocks[i].first.confinement));
    for(int di=0;di<3;di++){
      if(blocks[i].first.axis[di].ptr()){
        valtype *ax=blocks[i].first.axis[di].ptr();
        int sz=blocks[i].first.sz[di];
        if(per_flags&(0x1<<di)){
          // extending nonuniform mesh outside periodic boundaries on a one step
          ax = new valtype[sz+2];
          for(int j=0;j<sz;j++)
            ax[j+1]=blocks[i].first.axis[di][j];
          ax[0]=ip1[di]-(ip2[di]-ax[sz-1]);
          ax[sz+1]=ip2[di]+(ax[2]-ip1[di]);
          sz+=2;
        }
        if(block->SetMeshSteps(di,ax,ax+sz)<0)
          return message(vblMESS1,-1,"Some step is equal to zero in mesh %d, direction %d...\n",i,di);
        if(per_flags&(0x1<<di))
          delete[]ax;
      }
    }
    block->InitComponent(i ? fmt("m%d",i+1) : "m",0,0); // this name will be used in names of output files produced by mesh
    block->SetDirSplit(grad_split);
    store->cont->AddMesh(make_mngarg(block),blocks[i].second);
  }

  if(rad>=0)
    store->cont->SetRadialCoordinate(rad,rad_axis,rad_right);

  store->cont->SetPeriodicBoundaries(ip1,ip2-ip1,per_flags);

  store->cont->ShareMediumRegions(store->meds);

  for(size_t i=0;i<bonds.size();i++)
    store->cont->AddBoundaryRegion(bonds[i].first,bonds[i].second);

  for(size_t i=0;i<gens.size();i++){
    gens[i].first->gen->SetLengthUnit(lu);
    store->cont->AddGenerationRegion(gens[i].first,gens[i].second);
  }

  return 1;
}

int scExperiment::Prepare(){

  if(store->methods.size()==0)
    SetMethod("p_eh","1_10"); // default method

  for(size_t i=0;i<store->methods.size();i++)
    store->cont->set_method(store->methods[i],!i);

  valtype GR=0; // intergared generation rate

  store->cont->calibrate_Tn(T,n0,lu);

  // loading initial functions
  for(ini_funcmap_t::iterator it = ini_funcmap.begin(); it!=ini_funcmap.end(); ++it){
     if(it->first=="potential")
       store->cont->SetIniValueFunction(rpPsi,it->second);
     else if(it->first=="fermi_n")
       store->cont->SetIniValueFunction(rpFermi|rpE,it->second);
     else if(it->first=="fermi_p")
       store->cont->SetIniValueFunction(rpFermi|rpP,it->second);
  }

  if(store->cont->Prepare(&GR,dm)<0)
    return -1;

  Region_3 *dconf[3]={NULL,NULL,NULL}; // filtering detectors which are outside calculated space for 1D and 2D cases
  for(int i=0;i<3;i++){
    if(dim_type[i]<0){
      Vector_3 axis;
      axis[i]=1;
      valtype c=(ip1[i]+ip2[i])/2;
      valtype dc=1e-9;
      dconf[i] = GetPlate(axis,Vector_3(0,c-dc,0),2*dc);
    }
  }

  for(refmap<string,det_t>::iterator it=dets.begin(),e=dets.end();it!=e;++it){

    det_t &deti=*(it->second);

    RegionUnion<3> *un=NULL;
    if(dconf[0] || dconf[1] || dconf[2] || deti.conf){
      un = new RegionUnion<3>();
      for(int i=0;i<3;i++){
        if(dconf[i])
          un->AddRegion(dconf[i]);
      }
      if(deti.conf)
        un->AddRegion(deti.conf);
    }

    SpaceVectorSet *vect = new SpaceVectorSet; // here we will put detectors inside calculated space
    size_t vsz=0;
    int ind=deti.ind;
    if(deti.gr==0){
      filter_vset(store->cont.ptr(),un,*(store->gsets[ind]),vect->pos,NULL);
      vsz=store->gsets[ind]->size();
    }
    else if(deti.gr==1){
//      if(it->first==flux_det)
//        flux_norm=(store->bsets[ind]->get_total_surface()).norm();
      filter_vset(store->cont.ptr(),un,*(store->bsets[ind]),vect->pos,&(vect->dpos));
      vsz=store->bsets[ind]->size();
    }
    else if(deti.gr==2){
      filter_vset(store->cont.ptr(),un,*(store->csets[ind]),vect->pos,&(vect->dpos));
      vsz=store->csets[ind]->size();
    }

    if(rad>=0){ // for radial dimension dS should be renormalized
      for(size_t i=0;i<vect->dpos.size();i++){
        valtype R=(vect->pos[i])[rad]-rad_axis;
        valtype len=fabs(2*M_PI*R);
        vect->dpos[i]*=len;
      }
    }

    if(it->first==flux_det && !flux_norm)
      flux_norm=(vect->get_total_surface()).norm();

    TextRecorder<SpaceVectorSet> *rec = deti.flag ?
      new TextFluxRecorder<SpaceVectorSet>() :
      new TextRecorder<SpaceVectorSet>();

    rec->SetFormat(deti.format);

    int argtype=0, outtype;
    for(int i=0;i<3;i++){
      if(dim_type[i]>=0)
        argtype|=argx<<i;
    }
    if(!deti.flag){
      outtype=deti.output;
    }
    else // flux detectors measure flux only
      outtype=outJtot;

    DetectorSet<SpaceVectorSet,scFixInterp<scBlockContainer<> > > *det = 
      new DetectorSet<SpaceVectorSet,scFixInterp<scBlockContainer<> > >
      (make_mngarg(rec),store->cont.ptr(),1,make_mngarg(vect));

    vector<string> columns=out_names(deti.flag ? 0 : argtype, outtype);
    det->SetTypes(argtype,outtype,&columns);

    det->InitComponent(deti.name,0xffff,1,makevec<int>(2,tid,tid));
    det->Dump();
    store->dets.push_back(det);

    delete un;
  }

  if(GR){
    GR/=1e4/(lu*1e-6);
    if(flux_norm)
      message(vblMESS1,0,"integrated generation rate = %g, contact area = %g\n", GR, flux_norm);
    else
      message(vblMESS1,0,"integrated generation rate = %g\n", GR);
  }

  for(int i=0;i<3;i++)
    delete dconf[i];

  vector<BaseDetectorSet *> &sdets=store->dets;

  for(size_t di=0;di<sdets.size();di++)
    if(sdets[di]->Init()<0)
      return -1;

  message(vblMESS1,0,"\n");

  log->Usage();
  message(vblMESS1,0,"\n");

  phase=PH_GEOMETRY;

  return 1;
}


void scExperiment::SetIniValueFunction(const string& type, virt_unary_function<const Vector_3 &, vec_type> *func){
  ini_funcmap[type] = func;
}

int scExperiment::FillIniValues(int ini_repr){
  return store->cont->FillIniValues(ini_repr);
}

int scExperiment::Calculate(valtype dV, valtype V){
  vector<valtype> voltages(1, VStart);
  if (dV == 0) {
    if (!accomp(V, VStart))
      voltages.push_back(V);
  }
  else if(flux_det != "") // solar cell efficiency simulation
    voltages.push_back(VStart + dV);
  else {
    // fill voltages with values from VStart to V with step dV
    // More precisely, the first value in the vector should be VStart;
    // the last value should be greater than or equal to V if V>VStart
    // and less than or equal to V if V<VStart;
    // the other elements should be between VStart and V

    if ((V - VStart)*dV < 0)
      dV = -dV; // dV should be positive if V>VStart and negative otherwise

    for(valtype vCurrent = VStart+dV; !acless((V - voltages.back())*dV, 0.); vCurrent += dV)
      voltages.push_back(vCurrent);
  }

  return Calculate(voltages);
}

int scExperiment::Calculate(std::vector<valtype> V) {
  if (V.empty())
    V.push_back(0.);

  if(phase<PH_GEOMETRY){
    if(BuildGeometry()<0)
      return -1;
  }
  if(!(all_phases&PH_CALC))
    return 0;

  scBlockContainer<> &cont=*store->cont;
  vector<BaseDetectorSet *> &sdets=store->dets;

  if(!dump_iter){
    refmap<string,det_t>::iterator it=dets.begin();
    for(size_t di=0;di<sdets.size();di++,++it){
      int flag=it->second->flag;
      if(flag&DET_FLUX){
        string dname=sdets[di]->get_name();
        sdets[di]->SetFileName(dname + ".d");
        sdets[di]->SetExtra("V");
        if(sdets[di]->StartRecord()<0)return -1;
      }
    }
  }

  valtype rsd_first=-1;

  bool lastIteration = false;
  for(vector<valtype>::size_type vi = 0; !lastIteration && vi < V.size(); vi++) { // inside this loop we gradually increase the voltage and find solution for this new voltage value

    valtype v = V[vi];

 //   if(v){
      int res=cont.ApplyVoltage(v);
      if(res==0)
        return message(vblMESS1,-1,"No contacts where voltage is applied are specified\n");
      if(res<0) // something wrong at the analysis stage
        return -1;
 //   }

    message(vblMESS1,0,"Voltage %g.\n",v);
    cont.get_residue(1); // print residue of guess solution
    message(vblMESS1,0,"\n");

    int stepResult = 1, nStep = 0;
    for(; stepResult > 0; nStep++) { // calculation steps
      if(dump_iter) {
        if(dumpCurrentStep(nStep, v) < 0)
          return -1;
      }

      stepResult = cont.Step();
    }

    if(stepResult < 0)
      return -1; // something wrong with Newton iteration

    if(dumpCurrentStep(nStep, v) < 0)
      return -1;

    if(!vi)
      rsd_first = cont.get_residue();

    if(flux_det!="" && gens.size()) { // for solar cell efficiency simulation only
      string sj=theConfig->GetOutputDir()+flux_det+".d";
      valtype *sval;
      int sn1,sn2;
      int tres=read_table(sj.c_str(),&sval,sn1,sn2,true,true);
      if(tres<0)
        return -1;
      if(sn2!=2)
        return message(vblMESS1,-1,"File %s is not current flux detector file\n",sj.c_str());
      int se=(sn1-1)*sn2;
      valtype Isc=sval[1];
      valtype Ilast=sval[se+1];
      lastIteration = (Ilast<0 && Isc>0) || (Ilast>0 && Isc<0);

      // The sign of the current didn't change, append an element to V if we are at its end
      if (!lastIteration && V.size() > 1 && vi == V.size()-1)
        V.push_back(2*V.back() - V[V.size()-2]);

      if(lastIteration || vi == 0) {
//        valtype Voc=sval[se];
        valtype Voc=0;
        if(se){
          valtype Voc1=sval[se], Voc0=sval[se-sn2];
          valtype Ilast0=sval[se-sn2+1];
          Voc=(Ilast0*Voc1 - Ilast*Voc0)/(Ilast0-Ilast);
        }
        valtype ff=0;
        valtype Imax=0,Vmax=0;
        valtype Pmax=0;
        // interpolation should be applied
        for(int ffi=0;ffi<sn1;ffi++){
          int si=ffi*sn2;
          valtype Ic=sval[si+1];
          valtype Vc=sval[si];
          valtype Pc=Ic*Vc;
          if(Pc>Pmax){ // interpolation should be applied
            Imax=Ic;
            Vmax=Vc;
            Pmax=Pc;
          }
        }
        ff=Pmax/(Isc*Voc);
        FILE *f=fopen((theConfig->GetOutputDir()+"sc.dat").c_str(),"w");
        fprintf(f,"area\t%g\n",flux_norm);
        fprintf(f,"Isc,A\t%g\n",Isc);
        fprintf(f,"Isc,A/cm2\t%g\n",Isc/flux_norm);
        const valtype AM15=.097; // AM1.5 solar intensity, Wt/cm2
        if(Ilast<0){
          fprintf(f,"Voc,V\t%g\n",Voc);
          fprintf(f,"ff\t%g\n",ff);
          fprintf(f,"Pmax,W/cm2\t%g\n",Pmax/flux_norm);
          fprintf(f,"efficiency\t%g\n",Pmax/AM15/flux_norm);
        }
        fclose(f);
        if(name!=""){
          f=fopen((theConfig->GetOutputDir()+"sc_line.dat").c_str(),"w");
            
          fprintf(f,"%s\t%s",name.c_str(),rsd_first<=min_residue ? fmt("%g",Isc/flux_norm) : "-");
          if(Ilast>=0)
            fprintf(f,"\n");
          else{
            valtype rsd_last=cont.get_residue();
            fprintf(f,"\t%s\t%g",rsd_last<=min_residue ? fmt("%g",Voc) : "-",Pmax/AM15/flux_norm);
          }
          fclose(f);
        }
      }

      delete[] sval;
    }

  }

  if(!dump_iter){
    for(size_t di=0;di<sdets.size();di++)
      sdets[di]->EndRecord();
  }

  return 1;
}

int scExperiment::dumpCurrentStep(int nIteration, valtype voltage) {
  vector<BaseDetectorSet *> &sdets = store->dets;

  refmap<string,det_t>::iterator it=dets.begin();
  for(size_t di=0;di<sdets.size();di++,++it){
    int flag=it->second->flag;
    if(dump_iter || !(flag&DET_FLUX)){
      string dname=sdets[di]->get_name();
      if(dump_iter)
        dname += fmt("_%d_", nIteration);
      if(voltage)
        dname += fmt("%g", voltage);
      
      sdets[di]->SetFileName(dname + ".d");
      if(sdets[di]->StartRecord()<0)
        return -1;
    }
    else
      sdets[di]->SetExtra(fmt("%g",voltage));

    sdets[di]->StartNextRecord();
    sdets[di]->CompleteNextRecord();

    if(dump_iter || !(flag&DET_FLUX))
      sdets[di]->EndRecord();
  }
  if(dm) {
    string suf;
    if(voltage)
      suf=fmt(".%g",voltage);
    if(dump_iter)
      suf+=fmt("_%d_",nIteration); 
    store->cont->DumpMeshes(suf,dm&2 ? 1 : 0, 0xffff^outJtot);
  }

  return 1;
}


int scExperiment::DoDumpMeshes(const std::string &suffix){
  if(store->cont.ptr())
    return store->cont->DumpMeshes(suffix,dm&2 ? 1 : 0, 0xffff^outJtot);
  else
    return LOGERR(-1, "scExperiment::DoDumpMeshes: container is not initialized!\n",0);
}