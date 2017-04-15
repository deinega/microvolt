#include "sc_mesh.h"
#include "sc_fd.h"
#include "linsolv.hpp"
#include "vector_set.h"

template<class block_t>
scBlockContainer<block_t>::scBlockContainer(int *dim_type_, valtype T_):
  T(T_),rad(-1),rad_right(true),rad_axis(0),per_flags(0),SZ(0),SZJ(0),
  psi_ptr(NULL),ch(NULL),gen(NULL),useLinearEstimate(false),estimationNoise(0),medc(),
  bcond(NULL),fixed(NULL),dx_ar(NULL),x_ar(NULL),y_ar(NULL),
  jacstop(0), ch_it(0), methi(0), gen_eq(0), global_ifermi(0), pshifts(max_interp), pcoeffs(max_interp) {

  tmpi=scInterpolation(0,max_interp,pshifts.begin(),pcoeffs.begin());

  dim=0;

  for(int i=2;i>=0;i--){
    dim_type[i] = dim_type_[i];
    if(dim_type[i]>=0)
      dim++;
  }

  for(int car=0;car<2;car++){
    conc[car]=NULL;
    fermi[car]=NULL;
    J[car]=NULL;
  }

  clamp=1./5; // value from PC1D
  delta_x= 1e-5;

  calibrate_Tn(T,1e16,1);
  set_convergence_check();

  use_drc[0]=1;
  use_drc[1]=1;
  refl_psi=1;
  refl_fermi=1;
  if_fix_bulk=0;

  fix_trans.set_cont(this);

  fix_psi.set_cont(this);
  fix_sr.set_cont(this);

  fix_neum_psi.set_cont(this);
  fix_neum_psi.set_wrepr(rpPsi);
  fix_neum_fermi.set_cont(this);
  fix_neum_fermi.set_wrepr(rpFermi);

  fix_der_psi.set_cont(this);
  fix_der_psi.set_wrepr(rpPsi);
  fix_der_fermi.set_cont(this);
  fix_der_fermi.set_wrepr(rpFermi);

  fix_block.set_cont(this);
  fix_block_cyl.set_cont(this);

  fix_block_cyl.set_axis();
}

template<class block_t>
scBlockContainer<block_t>::~scBlockContainer(){
  delete[]psi_ptr;
  delete[]bcond;
  delete[]fixed;
  delete[]dx_ar;
  clean_geometry();
}

template<class block_t>
void scBlockContainer<block_t>::AddDefaultTimers(const vector<int> &p_id){
  int calc=p_id[0];
  tid.specify_nodes=AddTimersField("specify_nodes", calc);
  tid.analyze_geometry=AddTimersField("analyze_geometry", calc);
  tid.prepare_transfers=AddTimersField("prepare_transfers", calc);
  tid.dump=AddTimersField("dump", calc);
  tid.linsys=AddTimersField("linsys", calc);
  tid.block_jacobian=AddTimersField("block_jacobian", calc);
  tid.transfer_jacobian=AddTimersField("transfer_jacobian", calc);
  tid.jacobian_linsys=AddTimersField("jacobian_linsys", calc);
  tid.jacobian_step=AddTimersField("jacobian_step", calc);
  tid.jacobian_transfer_values=AddTimersField("jacobian_transfer_values", tid.jacobian_step);
  tid.matrix_prepare=AddTimersField("matrix_prepare", calc);
}


template<class block_t>
int scBlockContainer<block_t>::calibrate(valtype j, valtype l, valtype t, valtype n, valtype lu_){

  J0=j,l0=l,n0=n;
  lu=lu_*1e-6;

  valtype t0=t;

  // unscaled values
  coefexp=phys::q/(phys::kB*T);
  coefpsi=phys::eps0/phys::q;
  coefc=coefr=1;

  // scaling
  coefexp*=J0; // psi*k/q/T
  coefpsi*=J0/(l0*l0)/n0; // eps*psi/dx/dx = q*n
  coefc*=t0/(l0*l0);
  coefr*=t0/n0;

  if(accomp(coefpsi,1.))
    coefpsi=1;
  if(accomp(coefexp,1.))
    coefexp=1;
  if(accomp(coefc,1.))
    coefc=1;

  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::calibrate_Tn(valtype T, valtype n, valtype lu){
  valtype J=phys::kB/phys::q*T; // V
  valtype l=sqrt(phys::eps0*phys::kB*T/phys::q/phys::q/n); // m
  valtype D=1e-2; // m2/s
  valtype t=l*l/D; // s
  return calibrate(J,l,t,n,lu);
}

/*template<class block_t>
void scBlockContainer<block_t>::Dump(bool lim){
  if(dump)
    DumpOther(this);
}*/

template<class block_t>
int scBlockContainer<block_t>::clean_geometry(){
  boundaries.clear();
  gens.clear();
  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::test_rank(Vector_3 &pos) const{
  if(test_in_blocks(pos)<0)
    return -1;
  return 0;
}

template<class block_t>
template<class arg_t>
int scBlockContainer<block_t>::fill_interpolation_array(int cur_ind, const arg_t &arg,
  size_t sc_bdescr_t::*off, int id,int forse_id,int incompl,int non_nonlocal){

  int b;
  if(!forse_id){
    b=test_in_blocks(arg,id);
    if(b<0)
      return -1;
  }
  else
    b=id;
  int gind[3*block_t::vert_max];
  int mind[3*block_t::vert_max];
  valtype coef[3*block_t::vert_max];

  vector<arg_t> nonlocal;

  int n=blocks[b]->block->create_interpolation(arg,gind,mind,coef,&nonlocal,incompl);
  if(n<=0)
    return n;

  int i1=0;
  for(int i=0;i<n;i++){
    if(mind[i]<0) // interpolation on other meshes should be used
      continue;
    if(i1>=max_interp-1)
      return -2;
    pcoeffs[cur_ind+i1]=coef[i];
    pshifts[cur_ind+i1]=mind[i]+blocks[b]->*off;
    i1++;
  }
  int cur_n=i1;

  i1=0;
  for(int i=0;i<n;i++){
    if(mind[i]<0){ // interpolation on other meshes
      if(non_nonlocal)
        return -1;
      int n2=fill_interpolation_array(cur_ind+cur_n,nonlocal[i1++],off,0,0,incompl);
      if(n2<=0)
        return n2;
      for(int j=0;j<n2;j++){
        pcoeffs[cur_ind+cur_n+j]*=coef[i];
      }
      cur_n+=n2;
    }
  }
  return cur_n;
}

template<class block_t>
template<class arg_t>
scInterpolation scBlockContainer<block_t>::create_interpolation(const arg_t &arg,int id,int forse_id,
  int incompl,int non_nonlocal,int persistent){

  int n=fill_interpolation_array(0,arg,get_offset(arg),id,forse_id,incompl,non_nonlocal);

  if(n<=0)
    return scInterpolation();

  if(!persistent){ // make a reference
    tmpi.set_size(n);
    return tmpi;
  }
  else // making a copy
    return scInterpolation(0,n,pshifts.begin(),pcoeffs.begin(),2);
}

template scInterpolation scBlockContainer<RectBlock>::create_interpolation(const Vector_3 &arg,int id,int forse_id,int incompl,int non_nonlocal,int persistent);
template scInterpolation scBlockContainer<RectBlock>::create_interpolation(const pair<Vector_3,Vector_3> &arg,int id,int forse_id,int incompl,int non_nonlocal,int persistent);

template<class block_t>
int scBlockContainer<block_t>::classify_node(const Vector_3 &c,int lev){

  int type=test_bond(c);

  if(type==1)
//    return mFixBoundSurround; // contact boundary points surround insulator boundary points
    return mFixBound;
  else if(type==2)
    return mFixBound;
  else if(type==3)
    return mFixInt;
  else if(type==4)
    return mFixBound;
  else if(type<0)
    return type;

  for(size_t ii=0;ii<blocks.size();ii++){
    int ilev=blocks[ii]->level;
    if(ilev>lev && blocks[ii]->block->TestPoint(c)){ // belongs to other mesh with higher level
      return mTrans;
    }
  }

  return 0;
}

template<class block_t>
int scBlockContainer<block_t>::specify_adjacent_nodes(sc_bdescr_t *bdescr,typename block_t::mesh_it &pit,
int bflag,vector<pair<typename block_t::mesh_it,int> > &vp, map<int,pair<int,int> > &memmap){

  block_t *block = bdescr->block.ptr();

//  map<int,pair<int,int> > &memmap = block->get_memmap();

  for(int dn=0;dn<pit.get_adjacent_num();dn++){ // checking adjacent nodes

    // create iterator for dn-th adjacent node of pit or return -1 if this adjacent node does exist
    typename block_t::mesh_it bpit = pit.get_adjacent(dn);
    if(!(bpit!=block->end()))
      continue;
    typename block_t::node_t bp=*bpit;
    Vector_3 b=bp.GetCenter(); // center of the adjacent node

    int mflag = classify_node(b,bdescr->level);
    if(mflag<0)
      return mflag;

    if(!mflag && block->TestPoint(b,true)) // if it is inside, than it should be already stored in specify_nodes
      continue;

    int b_pind=bp.get_ind();

    if(memmap.find(b_pind)!=memmap.end()) // it is already stored
      continue;

    if(bflag==0){ // bflag - second element of pair in vp. it is used only for mFixBoundSurround nodes
      if(mflag&(mFixBound|mFixBoundSurround)){
        memmap[b_pind]=make_pair(0,mflag); // this is boundary node
        // this non-mFixBoundSurround boundary node could be surround by mFixBoundSurround boundary nodes
        if(!(mflag&mFixBoundSurround))
          vp.push_back(make_pair(bpit,1));
      }
      else if(mflag==mFixInt){ // this is not boundary node
        memmap[b_pind]=make_pair(0,mflag);
        vp.push_back(make_pair(bpit,0));
      }
      else if(mflag==mTrans)
        memmap[b_pind]=make_pair(0,mflag);
      else{ // external node which does not belong to boundary and higher level mesh (mflag=0)
//        scInterpolation form=create_interpolation(b);
//        if(!form.valid()){ // cannot be interpolated using other meshes
        if(test_in_blocks(b,0,1)<0){ // cannot be interpolated using other meshes
          memmap[b_pind]=make_pair(0,mInt); // making internal
          vp.push_back(make_pair(bpit,0));
        }
        else
          memmap[b_pind]=make_pair(0,mTrans);
      }
    }
    else{
      if(mflag&mFixBoundSurround)
        // surrounding non-mFixBoundSurround boundary node by mFixBoundSurround boundary nodes
        memmap[b_pind]=make_pair(0,mflag);
    }
  }
  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::specify_nodes(sc_bdescr_t *bdescr){

  block_t *block = bdescr->block.ptr();

  // key - grid index, value - memory index and flag (see mesh nodes classification)
//  map<int,pair<int,int> > &memmap = block->get_memmap();
  map<int,pair<int,int> > memmap;

  vector<pair<typename block_t::mesh_it,int> > vp[2]; // nodes requiring analyzing their adjacent nodes
  int ivp=0; // working index of vp

  // analyze all mesh nodes
  for(typename block_t::mesh_it pit=block->begin(),pe=block->end();pit!=pe;++pit){

    typename block_t::node_t p=*pit;
    Vector_3 c=p.GetCenter();

    int mflag = classify_node(c,bdescr->level);
    if(mflag<0)
      return mflag;

    if(mflag || !block->TestPoint(c,true)) // node is outside
      continue;

    for(int dn=0;dn<pit.get_adjacent_num();dn++){ // checking adjacent nodes
      if(!(pit.get_adjacent(dn)!=block->end()))
        // this is mesh boundary mesh node which is marked as internal
        // because it does not lie inside any boundary or higher level mesh.
        // possibly container is not entirely surrounded by boundaries
        return message(vblMESS1,-1,"Calculated space is not entirely surrounded by boundaries\n");
    }


    int gind=p.get_ind();
    memmap[gind]=make_pair(0,mInt); // store node as internal 
    vp[ivp].push_back(make_pair(pit,0)); // we need to analyze adjacent its nodes
  }

  while(vp[ivp].size()){
    // analyze nodes in vp[ivp]
    for(typename vector<pair<typename block_t::mesh_it,int> >::iterator vit=vp[ivp].begin(),ve=vp[ivp].end();vit!=ve;++vit)
      if(specify_adjacent_nodes(bdescr,vit->first,vit->second,vp[(ivp+1)%2],memmap)<0)
        return -1;
    // switch working index ivp
    vp[ivp].clear();
    ivp=(ivp+1)%2;
  }

  block->set_memmap(memmap);

  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::test_in_blocks(const Vector_3 &pos, int id, int strong) const{
//  id=0;
  int lev=-1;
  int index=-1;
  int n=blocks.size();
  for(int i=0;i<n;i++){
    int ii=(i+id)%n;
    int res = strong ? blocks[ii]->block->TestInterpolation(pos) : blocks[ii]->block->TestPoint(pos);
    if(res){
      if((index>=0 && lev<blocks[ii]->level) || index<0){
        index=ii;
        lev=blocks[ii]->level;
      }
    }
  }
  return index;
}

// TODO: eliminate duplicate code from RectBlock::set_local_media<>()
template<class block_t>
int scBlockContainer<block_t>::construct_local_media() {
  
  vector<Vector_3> points(SZ); // space position of mesh nodes

  for(size_t i = 0; i < blocks.size(); i++) {
    block_t *block = blocks[i]->block.ptr();

    // allocate memory for array with media properties, and assign values for its elements
    if (block->set_local_media(media_regions) < 0)
      return -1;

    typedef typename block_t::node_it node_it;
    for(node_it pit = block->begin_bound(), pe = block->end_bound(); pit != pe; ++pit) {
      int pind = blocks[i]->offset + (*pit).get_ind();
      points[pind] = (*pit).GetCenter();
    }

    for(node_it pit = block->begin_int(), pe = block->end_int(); pit != pe; ++pit) {
      int pind = blocks[i]->offset + (*pit).get_ind();
      points[pind] = (*pit).GetCenter();
    }
  }

  if(medc.build(*media_regions, points.begin(), points.end())<0)
    return -1;

  for(size_t i=0;i<blocks.size();i++){
    block_t *block=blocks[i]->block.ptr();
    block->reset_media(*this,blocks[i]->offset);
  }

  return 1;
}

template<class block_t>
valtype scBlockContainer<block_t>::test_gen(const Vector_3 &p) const{

  for (multimap<int, mngptr<scGenerationRegion> >::const_reverse_iterator git=gens.rbegin(), ge=gens.rend(); git!=ge; ++git){
    const scGenerationRegion &gr=*(git->second);
    if(!gr.reg.ptr() || gr.reg->TestPoint(p))
      return (*gr.gen)(p);
  }
  return 0;
}

template<class block_t>
int scBlockContainer<block_t>::test_bond(const Vector_3 &p,valtype **psi,valtype **sr, valtype *workf, Region_3 **reg)const{

  for (multimap<int, mngptr<scBoundaryRegion> >::const_reverse_iterator mit=boundaries.rbegin(), me=boundaries.rend(); mit!=me; ++mit){
    scBoundaryRegion &mr=*(mit->second);
    if(mr.reg->TestPoint(p)){
      if(psi)
        *psi=&(mr.psi);
      if(sr)
        *sr=mr.sr;
      if(workf)
        *workf=mr.workf;
      if(reg)
        *reg=mr.reg.ptr();
      return mr.type;
    }
  }
  if(rad>=0){ // check if this is axial node in radial scheme
    if(accomp(p[rad],rad_axis)) // axial node
      return 3;
    if(rad_right){
      if(p[rad]<rad_axis)
        return message(vblMESS1,-1,"Some mesh point is inconsistent with used radial system (%g, %g, %g)\n", p[0],p[1],p[2]);
    }
    else if(p[rad]>rad_axis)
      return message(vblMESS1,-1,"Some mesh point is inconsistent with used radial system (%g, %g, %g)\n", p[0],p[1],p[2]);
  }

  Vector_3 im=p.rpcell(p1,cell,per_flags); // periodic image
  if(im!=p) // periodic boundary conditions should be applied
    return 4;

  return 0;
}

template<class block_t>
int scBlockContainer<block_t>::find_normal_inside(Region_3 *reg, sc_bdescr<block_t> *sc_bdescr, 
typename block_t::node_t &p, Vector_3 &point, Vector_3 &nvect, Vector_3 *surfpp){

  block_t *block=sc_bdescr->block.ptr();
  Vector_3 c=p.GetCenter();

  Vector_3 norm; // averaged normal
  int nm=0; // number of nodes bordering with boundary node
  for(Vector_3 *v=p.adj_begin(), *ve=p.adj_end();v!=ve;v++){
    Region_3 *reg2;
    int flag2=test_bond(*v,NULL,NULL,NULL,&reg2);
    if(flag2==0 || flag2==3 || flag2==4){
      Vector_3 surfp,surfn;
      if(reg->TestEdge(c,*v,&surfp,&surfn)<0)
        return -1;
      norm+=surfn;
      nm++;
    }
  }
  if(nm){
    norm/=nm;
    Vector_3 surfp;
    if(reg->TestRay(c,norm,&surfp,&nvect)<0) // correct normal to the surface (norm -> nvect)
      return -1;
    if(surfpp)
      *surfpp=surfp;
    if(block->find_point(p,nvect,point)<0)
      return -1;
    for(int i=0;i<3;i++){
      // removing components of nonused dimension that could appear
      // as a result of find_normal_inside 
      if(dim_type[i]<0)
        point[i]=c[i];
    }
    return 1;
  }
  return 0;
}

template<class block_t>
int scBlockContainer<block_t>::record_internal(sc_bdescr<block_t> *sc_bdescr,const typename block_t::node_t &p,int bord,valtype gen, int ini_repr){

  int pind = sc_bdescr->offset + p.get_ind(); // memory index
 
  Vector_3 pos = sc_bdescr->block->get_position(p);
//  Vector_3 pos=p.GetCenter();

  const local_t& med = medc[pind];

  valtype psi;
  CarrierParams c;
  CarrierParams f;
  
  if(med.isOrganic()){
    
    psi = 0.;

    ini_funcmap_t::mapped_type func = ini_funcmap[rpPsi]; // locating ini function for potential
    if(func) // if found, initialize from it
      psi=(*func)(pos);
    //psi = -5.05;
    c = med.calcConcentration(CarrierParams(0.),psi); // Fermi value for contacts
    //c = med.get_N()/2;
    //c = med.calcConcentration(CarrierParams(5.05),psi); // Fermi value for contacts
    //c = med.get_N()*1e-5/n0; // I.V.:  this is for organics  
    psi = 0.;
  }
  else{
    ini_funcmap_t::mapped_type func = ini_funcmap[rpPsi]; // locating ini function for potential
    if(func) // if found, initialize from it
      psi=(*func)(pos);
    else
      psi = med.get_workfunc();
    c = med.get_conc0(); // equilibrium concentration
  }

  ch[pind] = med.get_dop()/n0;
  if(ini_repr&rpPsi)
    psi_ptr[pind] = psi/J0; // initial guess value for potential

  if (gen) {
    int minor = c[0] > c[1]; // minority carriers
    valtype tau = med.getRecombParams().tau[minor]; // tau is defined by minority carriers
    CarrierParams gen_eq_car(gen_eq);
    c += gen*tau*gen_eq_car; // add concentration hold by generation (G = \Delta n/ \tau)
  }
  
  f = med.calcQFermi(c, psi);

  if(med.isOrganic()){
    psi = -f[1];
    if(ini_repr&rpPsi)
      psi_ptr[pind] = psi/J0; // rescaling for organics
    f[0]=f[1]=0.;
  }
  
  //CarrierParams f2 = med.calcQFermi(c2, psi);
  CarrierParams c_post = med.calcConcentration(f,psi); // just to ensure that c = c_post 
  //CarrierParams c2_post = med.calcConcentration(f2,psi);

  for(CarrierParams::size_type car = 0; car < c.size(); car++){
    if(ini_repr&(rpFermi|(car ? rpH : rpE))){
      if(fermi_guess[car].ptr()) // substituting initial guess value
        fermi[car][pind] = (-fermi_guess[car]->operator()(pos)+global_ifermi)/J0;
      else
        fermi[car][pind] = f[car]/J0;
    }
  }

  //valtype cor2=(car ? 1 : -1)*med.get_Eg()/2/J0*gen_eq_car[car]; // fix is defined by bandgap

  if(!bord && if_fix_bulk)
    fix_block.record(p.get_ind(),sc_bdescr->id);

  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::record_bound(sc_bdescr<block_t> *sc_bdescr,
const typename block_t::node_it &pit,typename block_t::node_t &p,int aflag){

//  block_t *block=sc_bdescr->block.ptr();

  Vector_3 c=p.GetCenter();

  valtype *psi;
  valtype *sr;
  Region_3 *reg;
  valtype metal_workf; // contact metal workfunction for Schottky contact
  int flag=test_bond(c,&psi,&sr,&metal_workf,&reg);

  if(flag!=1 && !(aflag&scIni))
    return 0;

  int pind=sc_bdescr->offset+p.get_ind();
//  bcond[pind]=flag;

  const local_t& med = medc[pind];

  valtype wf_delta=0.;
  if(metal_workf!=VEC_INFTY){ // wf difference for Schottky contact, VEC_INFTY indicates simple Ohmic contact
    wf_delta = -med.get_workfunc() - metal_workf;
  }
 
  if(flag==1){ // metal or Schottky metal
    valtype pb = med.get_workfunc();
    psi_ptr[pind]=(pb+*psi+wf_delta)/J0;
    
    //if(med.isOrganic())
     // psi_ptr[pind] = *psi;

    if(use_drc[0]){
      fixed[pind]=1;
    }
    else if(aflag&scIni){
      fix_psi.record(p.get_ind(),sc_bdescr->id,psi,pb,wf_delta);
    }
    if(sr[0]==VEC_INFTY && sr[1]==VEC_INFTY && use_drc[1]){
      for(int car=0;car<2;car++)
        fermi[car][pind]=*psi/J0;
      fixed[pind+SZ]=fixed[pind+2*SZ]=1;
    }
    else if(aflag&scIni){ // finite surface recombination
      Vector_3 point,nvect,surfp;
      if(find_normal_inside(reg,sc_bdescr,p,point,nvect,&surfp)<=0){
//        find_normal_inside(reg,sc_bdescr,p,point,nvect);
        return message(vblMESS1,-1,"Cannot apply boundary condition for meal contact (%g, %g, %g)\n", c[0],c[1],c[2]);
      }
      if(fix_sr.record(p.get_ind(),sc_bdescr->id,pit.bound_index(),flag,c,point,nvect,surfp,wf_delta,med,sr)<0)
        return message(vblMESS1,-1,"More than one mesh are close to the same metal contact\n");

      if(med.isOrganic()){
        for(int car=0;car<2;car++) // A.D. - does it work with semiconductors?
          fermi[car][pind]=*psi/J0; //-(pb+wf_delta)/J0; // I.V.: initial value is equal to the metal Fermi level
      }
    }
  }
  
  if(flag==2){ // dielectric
    Vector_3 point,nvect,surfp;   
    if(find_normal_inside(reg,sc_bdescr,p,point,nvect,&surfp)<=0)
      return message(vblMESS1,-1,"Cannot apply boundary condition for dielectric interface (%g, %g, %g)\n", c[0],c[1],c[2]);

    scInterpolation form=create_interpolation(point,sc_bdescr->id,0,1); // incomplete interpolation is allowed
    if(!form.valid()){
      return message(vblMESS1,-1,"Can not apply neumann boundary conditions (%g, %g, %g):\n"
        " dielectric interface is too sharp or\n"
        " more than one mesh are close to the same dielectric interface\n", c[0],c[1],c[2]);
    }
    if(refl_psi==1)
      fix_neum_psi.record(p.get_ind(),sc_bdescr->id,form,med);
    else if(refl_psi==2)
      fix_der_psi.record(p.get_ind(),sc_bdescr->id,c,nvect);

    if(sr[0]==0 && sr[1]==0 && refl_fermi==1)
      fix_neum_fermi.record(p.get_ind(),sc_bdescr->id,form,med);
    else if(sr[0]==0 && sr[1]==0 && refl_fermi==2)
      fix_der_fermi.record(p.get_ind(),sc_bdescr->id,c,nvect);
    else{
      if(fix_sr.record(p.get_ind(),sc_bdescr->id,pit.bound_index(),flag,c,point,nvect,surfp,wf_delta,med,sr)<0){
        return message(vblMESS1,-1,"More than one mesh are close to the same dielectric interface\n");
      }
    }
  }

  if(flag==3) // cylindrical axis
    fix_block_cyl.record(p.get_ind(),sc_bdescr->id);

  if(flag==4){
    Vector_3 im=c.rpcell(p1,cell,per_flags);
    scInterpolation form=create_interpolation(im,sc_bdescr->id);
    if(!form.valid()){
      return ::pmessage(vblMESS1,-1,"Periodic boundary conditions are not consistent with mesh layout\n");
    }
    fix_trans.record(p.get_ind(),sc_bdescr->id,form);
    fixed[pind]=fixed[pind+SZ]=fixed[pind+2*SZ]=2;
  }

  if(flag==0){ // mesh transfer
    int next=(sc_bdescr->id+1)%blocks.size();
    scInterpolation form=create_interpolation(c,next);
    if(!form.valid()){
      return ::pmessage(vblMESS1,-1,"BlockContainer is not covered by meshes properly\n");
    }
    fix_trans.record(p.get_ind(),sc_bdescr->id,form);
    fixed[pind]=fixed[pind+SZ]=fixed[pind+2*SZ]=2;
  }

  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::analyze(int aflag, valtype *GR, int ini_repr){

  if(aflag&scIni)
    fix_trans.clear();

  valtype V=0; // total volume
  valtype G=0; // integrated generation rate

  for(size_t i=0;i<blocks.size();i++){

    block_t *block=blocks[i]->block.ptr();

    // boudary points
    for(typename block_t::node_it pit=block->begin_bound(),pe=block->end_bound();pit!=pe;++pit){
      typename block_t::node_t p=*pit;
      valtype gr=test_gen(p.GetCenter());
      valtype dv=p.GetControlVolume();
      if(rad>=0){
        Vector_3 c=p.GetCenter();
        valtype R=c[rad]-rad_axis;
        dv*=fabs(2*M_PI*R);
      }
      V+=dv;
      G+=dv*gr;

      if(aflag&scIniVal){
        if(record_internal(blocks[i],p,1,gr,ini_repr)<0)
          return -1;
      }
      
      if(aflag&scIniGen){
        int pind=blocks[i]->offset + p.get_ind(); // memory index
        gen[pind]=gr;
      }

      if(aflag&scIniVoltage){
        if(record_bound(blocks[i],pit,p,aflag)<0)
          return -1;
      }
    }

    // internal points
    for(typename block_t::node_it pit=block->begin_int(),pe=block->end_int();pit!=pe;++pit){
      typename block_t::node_t p=*pit;
      valtype gr=test_gen(p.GetCenter());
      valtype dv=p.GetControlVolume();
      if(rad>=0){
        Vector_3 c=p.GetCenter();
        valtype R=c[rad]-rad_axis;
        dv*=fabs(2*M_PI*R);
      }
      V+=dv;
      G+=dv*gr;

      if(aflag&scIniVal){
        if(record_internal(blocks[i],p,0,gr,ini_repr)<0)
          return -1;
      }

      if(aflag&scIniGen){
        int pind=blocks[i]->offset + p.get_ind(); // memory index
        gen[pind]=gr;
      }

    }
  }
  if(GR)
    *GR=G*phys::q;
  // update external mesh nodes
  transfer_values(fix_trans,rpPsi|rpFermi|rpE|rpH,rpConc);

  return 1;
};

template<class block_t>
int scBlockContainer<block_t>::clear_geometry(){
  medc.clear();
  fix_psi.clear();
  fix_sr.clear();
  fix_neum_psi.clear();
  fix_neum_fermi.clear();
  fix_der_psi.clear();
  fix_der_fermi.clear();
  fix_block.clear();
  fix_block_cyl.clear();
  fix_trans.clear();

  for(size_t i=0;i<blocks.size();i++){
    block_t *block=blocks[i]->block.ptr();
    block->clear_geometry();
  }
  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::AllocateMemory(int dmp){

  start(tid.specify_nodes);

  message(vblMESS1,0,"Organizing meshes memory layout...\n");

  int min_lev=INT_INFTY, max_lev=-INT_INFTY+1;

  for(size_t i=0;i<blocks.size();i++){ // find minimal and maximal level for specified meshes
    int lev=blocks[i]->level;
    if(lev>max_lev)
      max_lev=lev;
    if(lev<min_lev)
      min_lev=lev;
  }

  for(int li=min_lev;li<=max_lev;li++){
    for(size_t i=0;i<blocks.size();i++){
      int lev=blocks[i]->level;
      if(li!=lev)
        continue;

      if(specify_nodes(blocks[i])<0) // specify nodes which will be used by this mesh
        return -1;

      block_t *block=blocks[i]->block.ptr();
      block->organize_memory(); // mesh organizes memory layout for nodes which will be used

      blocks[i]->offset=SZ;
      SZ+=block->GetSize();
      blocks[i]->offset_cur=SZJ;
      SZJ+=block->GetCurrentSize();
    }
  }

  if(!SZ){
    ::pmessage(vblMESS1,0,"Zero size of all meshes\n");
    comm->set_exit();
  }
  if(comm->check_exit())return -1;

  int max_matr=0; // maximal amount of blocks in Jacobian
  for(size_t i=0;i<methods.size();i++){
    for(int j=0;j<methods[i].n;j++){
      if(methods[i].meth[j].nr>max_matr)
        max_matr=methods[i].meth[j].nr;
    }
  }

  // it should be max_matr, not 3. Why did Rusakov put 3?
  const int x_ar_count = 3; // number of equation groups (p, e, h)

  try{
    psi_ptr = new valtype[9*SZ+2*SZJ];

    bcond = new int[SZ];
    fixed = new int[3*SZ];
    dx_ar = new valtype[3*x_ar_count*SZ];
#ifdef USE_PARDISO
    if(blocks.size()==1)
      ls_pc2.set_sparce(4*get_dim()*(block_t::vert_dir_max));
    else // for possible interpolation between nodes
      ls_pc2.set_sparce(2*4*get_dim()*(block_t::vert_dir_max));
#endif
    ls_pc2.init(max_matr*SZ,3);
  }
  catch(...){
    ::pmessage(vblMESS1,0,"Not enough memory for specified mesh(es)\n");
    comm->set_exit();
  }
  if(comm->check_exit())return -1;

  for(int i=0;i<9*SZ+2*SZJ;i++)
    psi_ptr[i]=0;
  for(int i=0;i<SZ;i++)
    bcond[i]=0;
  for(int i=0;i<3*SZ;i++)
    fixed[i]=0;
  for(int i=0;i<3*x_ar_count*SZ;i++)
    dx_ar[i]=0;
  
  ls_p.init(SZ,0);
  ls_p.share(ls_pc2);

  for(size_t mi=0;mi<methods.size();mi++){
    for(int j=0;j<methods[mi].n;j++){
      if(methods[mi].meth[j].nr==2){
        ls_c2.init(2*SZ,0);
        ls_c2.share(ls_pc2);
        ls_pc.init(2*SZ,0);
        ls_pc.share(ls_pc2);
        break;
      }
    }
  }
  ls_c.init(SZ,0);
//    ls_c[car].share(ls30,(1+car)*SZ);
  ls_c.share(ls_pc2);

  for(int car=0;car<2;car++){
    conc[car] = psi_ptr+(1+car)*SZ;
    fermi[car] = psi_ptr+(3+car)*SZ;
    J[car]=psi_ptr+7*SZ+car*SZJ;
    excitons[car]=psi_ptr+7*SZ+2*SZJ+car*SZ;
  }
  ch=psi_ptr+5*SZ;
  gen=psi_ptr+6*SZ;
  y_ar=dx_ar+x_ar_count*SZ;
  x_ar=dx_ar+2*x_ar_count*SZ;

  stop(tid.specify_nodes);

  int local_size=0; // local meshes arrays size
  for(size_t i=0;i<blocks.size();i++){
    block_t *block=blocks[i]->block.ptr();
    int lsz=block->reset(*this,blocks[i]->offset,blocks[i]->offset_cur);
    if(lsz<0)
      return -1;
    local_size+=lsz;
    // filling bcond for DumpMeshPoints
    for(typename block_t::node_it pit=block->begin_bound(),pe=block->end_bound();pit!=pe;++pit){
      typename block_t::node_t p=*pit;
      Vector_3 c=p.GetCenter();
      int pind=blocks[i]->offset+p.get_ind();
      bcond[pind]=test_bond(c);
    }
    if(dmp&1)
      block->DumpMesh("",(dmp&2)!=0,0);
  }

  log->add_memusage("global arrays",SZ*((9+3*x_ar_count)*sizeof(valtype)+4*sizeof(int)+sizeof(local_t))+SZJ*2*sizeof(valtype));
  log->add_memusage("local mesh arrays",local_size);
  log->add_memusage("matrix",ls_pc2.data_size(7));

  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::BuildGeometry(valtype *GR, int ini_repr){
  clear_geometry();
  try{
    // allocate memory for array with media properties, and assign values for its elements
    construct_local_media();
  }
  catch(...){
    ::pmessage(vblMESS1,0,"Not enough memory for specified mesh(es)\n");
    comm->set_exit();
  }
  if(comm->check_exit())return -1;

  if(zero_pos != VEC_INFTY){
    if(test_bond(zero_pos, NULL, NULL, &global_ifermi) == 1 && global_ifermi != VEC_INFTY)
      global_ifermi = -global_ifermi;
    else {
      scRegionMediumMap::MediumPtr med = media_regions->at(zero_pos);
      if(!med)
        return message(vblMESS1,-1,"The zero potential point is not covered by media\n");;

      scAbstractMedium::LocalPtr medLocal = med->createLocalMedium(zero_pos);
      global_ifermi = medLocal->get_workfunc();
    }
  }

  start(tid.analyze_geometry);
  message(vblMESS1,0,"Analyzing geometry...\n");
  if(analyze(scIniAll, GR, ini_repr)<0)
    comm->set_exit();
  if(comm->check_exit())return -1;

  stop(tid.analyze_geometry);

  start(tid.prepare_transfers);
  fix_trans.init();
  stop(tid.prepare_transfers);

  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::transfer_values(FixTransfer<scBlockContainer<block_t> > &trans,int repr,int srp,int flag){

  for(size_t oi=0;oi<trans.i_ind.size();oi++){
    if(flag>=0)
      trans.set_value(oi,repr,srp,flag);
    else
      trans.set_previous_value(oi,repr,srp);
  }
  return 1;
}

template<class block_t>
valtype scBlockContainer<block_t>::calc_mesh_y(const meth_cont_t &mc,valtype *y_ptr){

  valtype tot_residue=0;

  int szi=0;
  for(int j0=0;j0<mc.snum;j0++){
    for(size_t i=0;i<blocks.size();i++){
      int offset=blocks[i]->offset;
      tot_residue+=blocks[i]->block->internal_nodes(mc.gyrepr[j0],y_ptr+szi+offset);
    }
    szi+=mc.fstep[j0]*mc.sz[j0];
  }

  return tot_residue;
}

template<class block_t>
template<class fix_t>
valtype scBlockContainer<block_t>::calc_fix_y(fix_t &fix,const meth_cont_t &mc,valtype *y_ptr){

  valtype residue=0;

  int szj=0;
  int blocki=0;
  for(int j0=0;j0<mc.snum;j0++){
    for(int js=0;js<mc.fstep[j0];js++,blocki++){
      int repy=mc.meth.yrepr[blocki]; // representation of equation
      if(!fix.working_repr(repy,1))
        continue;

      for(size_t oi=0;oi<fix.out.size();oi++){
        int out=fix.out[oi];

        if(mc.fx[j0][out+js*mc.sz[j0]]) // fixed variable or external mesh point, corresponding y = 0
          y_ptr[out+js*mc.sz[j0]+szj]=0;
        else{
          valtype y=fix.y(oi,repy);
          residue+=y*y;
          y_ptr[out+js*mc.sz[j0]+szj]=y;
        }
      }
    }
    szj+=mc.fstep[j0]*mc.sz[j0];
  }

  return residue;
}

template<class block_t>
valtype scBlockContainer<block_t>::calc_y(const meth_cont_t &mc,valtype *y_ptr){

  for(int i0=0,szi=0;i0<mc.snum;i0++){
    for(int i=0;i<mc.fstep[i0]*mc.sz[i0];i++){
      y_ptr[szi+i]=0;
    }
    szi+=mc.fstep[i0]*mc.sz[i0];
  }

  for(int i=0;i<SZ;i++)
    set_repres(*this,mc.srp,i);

  valtype residue=0;
  if(!if_fix_bulk)
    residue+=calc_mesh_y(mc,y_ptr);

  residue+=calc_fix_y(fix_sr,mc,y_ptr);
  residue+=calc_fix_y(fix_psi,mc,y_ptr);
  residue+=calc_fix_y(fix_neum_psi,mc,y_ptr);
  residue+=calc_fix_y(fix_neum_fermi,mc,y_ptr);
  residue+=calc_fix_y(fix_der_psi,mc,y_ptr);
  residue+=calc_fix_y(fix_der_fermi,mc,y_ptr);
  residue+=calc_fix_y(fix_block,mc,y_ptr);
  residue+=calc_fix_y(fix_block_cyl,mc,y_ptr);

/*  residue=0; // just to check again
  for(int i0=0,szi=0;i0<snum;i0++){
    for(int i=0;i<fstep[i0]*sz[i0];i++)
      residue+=y_ptr[szi+i]*y_ptr[szi+i];
    szi+=fstep[i0]*sz[i0];
  }*/
  return residue;
}

template<class block_t>
int scBlockContainer<block_t>::calc_transfer_jacobian(FixTransfer<scBlockContainer<block_t> > &trans,const meth_cont_t &mc){

  valtype ny[3*(2*dim_max+1)];
  int nind[2*dim_max+1];
  bool touch[2*dim_max+1];

  int ii=0;
  int nindi=0;

  for(size_t oi=0;oi<trans.x_ind.size();oi++){ // internal points of next mesh
    int x_ind=trans.x_ind[oi];

    if(fixed[x_ind]==1)
      continue;

    int szi=0;
    int blocki=0;
    for(int i0=0;i0<mc.snum;i0++){
      for(int is=0;is<mc.fstep[i0];is++,blocki++){

        int rep_set=mc.meth.xrepr[blocki]; // used in transfer.set_value

        valtype iprev=mc.x[i0][mc.sz[i0]*is+x_ind];
        mc.x[i0][mc.sz[i0]*is+x_ind]+=delta_x; // shift value in internal point of next mesh
        set_repres(*this,mc.srp,x_ind);

        int sz=trans.backs[oi].n;

        for(int ii=0;ii<sz;ii++){ // calculate value in external point of mesh
          int toi=trans.backs[oi].igroup[ii];
          trans.set_value(toi,rep_set,mc.srp);
        }

        for(int ii=0;ii<sz;ii++){

          int toi=trans.backs[oi].igroup[ii];
          int ntoi=trans.backs[oi].ogroups[ii];
          int o_grsz=trans.o_grsz[toi];
          for(int nindi2=0;nindi2<o_grsz;nindi2++){ // filling indices of influenced equations in internal mesh points
            nind[nindi2]=trans.o_gr_ind[ntoi+nindi2];
            touch[nindi2]=0;
          }

          int id=trans.i_mesh[toi]; // to what block

          for(int j=0;j<3*(2*dim+1);j++)
            ny[j]=0;

          int szj=0;
          int szjf=0;
          for(int j0=0;j0<mc.snum;j0++){
            blocks[id]->block->chosen_nodes(mc.gyrepr[j0],nind,touch,o_grsz,0,ny+(2*dim+1)*szjf);
            for(int js=0;js<mc.fstep[j0];js++){
              for(int j=0;j<o_grsz;j++){ // cycle for influenced equations in internal mesh points
                valtype delta_y=ny[(szjf+js)*(2*dim+1)+j]-mc.ls->v(blocks[id]->offset+js*mc.sz[j0]+szj+nind[j]);
                // set delta_y/delta_x
                if(mc.ls->set_m(mc.sz[i0]*is+szi+x_ind,blocks[id]->offset+js*mc.sz[j0]+szj+nind[j],delta_y/delta_x)<0)
                  return -1;
              }
            }
            szj+=mc.fstep[j0]*mc.sz[j0];
            szjf+=mc.fstep[j0];
          }
        }

        for(int ii=0;ii<sz;ii++){ // restore value in external point
          int toi=trans.backs[oi].igroup[ii];
          trans.set_previous_value(toi,rep_set,mc.srp);
        }
        trans.reset_previous_value();

        mc.x[i0][mc.sz[i0]*is+x_ind]=iprev; // restore value in internal point of next mesh (shift back)
        set_repres(*this,mc.srp,x_ind);
      }
      szi+=mc.fstep[i0]*mc.sz[i0];
    }
  }

  return 1;
}

template<class block_t>
template<class fix_t>
int scBlockContainer<block_t>::calc_fix_jacobian(fix_t &fix,const meth_cont_t &mc){

  for(size_t oi=0,isz=0;oi<fix.out.size();oi++){ // cycle for points recorded in fix
    int out=fix.out[oi];

    int sz=fix.in_sz[oi];

    for(size_t ii=isz;ii<isz+sz;ii++){ // cycle for points influencing point oi
      int in=fix.in[ii];

      int szi=0;
      int blocki=0;
      for(int i0=0;i0<mc.snum;i0++){

        for(int is=0;is<mc.fstep[i0];is++,blocki++){

          int rep_x=mc.meth.xrepr[blocki]; // representation to change

          if(!fix.working_repr(rep_x,0))
            continue;
          if(mc.fx[i0][mc.sz[i0]*is+in])
            continue; // constant or transfer variable x

          valtype iprev=mc.x[i0][mc.sz[i0]*is+in];
          mc.x[i0][mc.sz[i0]*is+in]+=delta_x;
          set_repres(*this,mc.srp,in);

          int szj=0;
          int blockj=0;
          for(int j0=0;j0<mc.snum;j0++){

            for(int js=0;js<mc.fstep[j0];js++,blockj++){

              int repy=mc.meth.yrepr[blockj]; // equation representation
              int const_y=0;
              if(!fix.working_repr(repy,1))
                continue;

              valtype y=fix.y(oi,repy);
              valtype delta_y=y-mc.ls->v(out+js*mc.sz[j0]+szj);

              if(mc.ls->set_m(mc.sz[i0]*is+szi+in,out+js*mc.sz[j0]+szj,delta_y/delta_x)<0)
                return -1;
            }
            szj+=mc.fstep[j0]*mc.sz[j0];
          }

          mc.x[i0][mc.sz[i0]*is+in]=iprev;
          set_repres(*this,mc.srp,in);

        }
        szi+=mc.fstep[i0]*mc.sz[i0];
      }
    }
    isz+=fix.in_sz[oi];
  }
  
  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::calc_jacobian(const meth_cont_t &meth){

  // each block calculates Jacobian for its internal points 
  // which are influenced by other internal points only

  if(!if_fix_bulk){
    start(tid.block_jacobian);
    for(size_t i=0;i<blocks.size();i++){
      if(!blocks[i]->block->GetSize())
        continue;
      if(blocks[i]->block->Jacobian(meth)<0)
        return message(vblMESS1,-1,"Jacobian cannot be updated by mesh(es)\n");
    }
    stop(tid.block_jacobian);
  }

  // filling Jacobian for fixes (metal or dielectric boundary conditions, etc.)
  calc_fix_jacobian(fix_sr,meth);
  calc_fix_jacobian(fix_psi,meth);
  calc_fix_jacobian(fix_neum_psi,meth);
  calc_fix_jacobian(fix_neum_fermi,meth);
  calc_fix_jacobian(fix_der_psi,meth);
  calc_fix_jacobian(fix_der_fermi,meth);
  calc_fix_jacobian(fix_block,meth);
  calc_fix_jacobian(fix_block_cyl,meth);

  // filling Jacobian for internal points of each block
  // which are influenced not only by internal points, but by external points as well
  start(tid.transfer_jacobian);
  if(calc_transfer_jacobian(fix_trans,meth)<0)
    return message(vblMESS1,-1,"Jacobian cannot be updated for transfer nodes\n");
  stop(tid.transfer_jacobian);

  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::Step(){

  meth_chain_t &ch = methods[methi];
  int moving=1; // 0 if local minimum is achieved
	
  for(int i=0;i<ch.n;i++){

    meth_cont_t mc(*this,ch.meth[i]);
    string estr=ch.meth[i].code;
//    if(ch.num[i]>1)
//      estr+=str_sprintf("%d",ch.num[i]);

    message(vblMESS1,0,"%s ",estr.c_str());

    valtype rsd=-1,rsd2=-1; // initial and next residue
    int res=0;
    for(int j=0;j<ch.num[i];j++){
      valtype rsd_temp=-1;
      // if j==0 then initial residue should be recorded to rsd, otherwise it is recorded to dummy rsd_temp
      int res_temp=Jacobian(mc,j==0 ? rsd : rsd_temp,rsd2);
      if(j==0 || res_temp<0)
        res=res_temp;
      if(res_temp<=0)
        break;
    }
    if(res<0){
      message(vblMESS1,0,"System is not well conditioned...\n");
      return -1;
    }

    get_residue(1);
    if(rsd>=0){
			rsd2_stored[i] = rsd2; // storing the final residue for the method i in chain
      if(res==0){
        jacstop++;
        message(vblMESS1,0," %g <->",rsd2/SZ);
      }
      else{
        message(vblMESS1,0," %g -> %g",rsd/SZ,rsd2/SZ);
        jacstop=0;
      }
      message(vblMESS1,0,"\n");

      if(jacstop>=ch.n){
				// checking that residue is not spoiled for other methods in chain
				int res_other = 0;
				for(int ii =0; ii< ch.n ; ii++){
					valtype rsd2 = calc_y(meth_cont_t(*this,ch.meth[ii]),y_ar);

					if(rsd2> rsd2_stored[ii] && rsd2 > residue_threshhold)
						res_other = 1;  

				  //res_other = (1.-rsd2/rsd2_stored[ii] > (rsd2 > residue_threshhold ? stop_ratio : small_stop_ratio)) ? 1 : 0; // same as end of Jacobian, make a function
					if(res_other)
						break;
				}
				if(res_other == 0){ // all methods not moving
          message(vblMESS1,0,"Residue minimum is achieved\n");
          moving=0;
				}
				else
					jacstop = 0;
      }
    }
    if(!moving)break;
  }
  ch_it++;

  if((!moving || (ch.tot_num>0 && ch_it>=ch.tot_num))){
//    methi=(methi+1)%methods.size();
    if(methi>=int(methods.size())-1)
      moving=0;
    //else{
    else if(moving){ // appply next method only when needed - is it correct ? The idea is to move to the next method when previous method does not move (moving=0)
      methi++;
      moving=1;
      jacstop=0;
      ch_it=0;
    }
  }

  for(size_t i=0;i<blocks.size();i++){
    block_t *block=blocks[i]->block.ptr();
    block->update_J();
  }


#ifdef USE_PARDISO
  //if(methi==methods.size()-1)
  if(methods.size()==1)
    ls_pc2.delete_first_record_arays();
#endif
  return moving;
}

template<class block_t>
int scBlockContainer<block_t>::DumpMeshes(const string &suf, bool dif, int outtype){

  calc_y(meth_cont_t(*this,method_t("peh")),y_ar);

  start(tid.dump);
  for(size_t i=0;i<blocks.size();i++)
    blocks[i]->block->DumpMesh(suf,dif,outtype);
  stop(tid.dump);
  return 1;
}

template<class block_t>
int scBlockContainer<block_t>::Jacobian(const meth_cont_t &mc,valtype &rsd,valtype &rsd2){

  mc.ls->zero_matrix();
  mc.ls->zero_vector();

  // calculate vector y and record its squared norm (sum y*y) to rsd, rsd2
  rsd=rsd2=calc_y(mc,mc.ls->get_v());
#ifdef LIMIT_RESIDUE_GROWTH
  valtype rsdFull = calc_y(meth_cont_t(*this, method_t("peh")), y_ar);
#endif

//  valtype rsd_psi, rsd_y, rsd_n, rsd_p, rsd_sr, rsd_out;

  // calculate Jacobian
  if(calc_jacobian(mc)<0)
    return -1;

//  ls.prepare_center();
//  prepare_val(ls);
//  ls.prepare_max();

  /* Jacobian is not calculated for external nodes of each block (let denote their indeces as i)
  Corresponding columns and raws of Jacobian J, and elements of residue vector dy are zero:
  J_{i,j} = J_{j,i} = dy_i = 0 (i - external nodes, j - any nodes)
  However equation J*dx = dy can be solved if J is nondegenerate (otherwise solver cannot work)
  Therefore wee need to assign elements J_{i,i} with some nonzero value
  Since dy_i = 0 it will lead to dx_i = 0
  */
  int deg=mc.ls->make_nondegenerate();
  start(tid.matrix_prepare);
  mc.ls->end_record();
  stop(tid.matrix_prepare);

  int szi=0;
  for(int j0=0;j0<mc.snum;j0++)
    szi+=mc.fstep[j0]*mc.sz[j0];

  start(tid.jacobian_linsys);
  if(mc.meth.snum_gs>1){
    // solving by iterative Gauss–Seidel method
    // see http://en.wikipedia.org/wiki/Gauss_seidel
    if(gauss_seidel(*mc.ls,szi,dx_ar,mc.meth.snum_gs,NULL,mc.meth.it_gs)<0)
      return -1;
  }
  else{
    // direct solving
    if(mc.ls->operator()(dx_ar)<0)
      return -1;
  }
  stop(tid.jacobian_linsys);


  mc.ls->start_record();

/*  valtype ymax=0;
  for(int i=0;i<szi;i++){
    valtype ycur=fabs(ls.get_v()[i]);
    if(ycur>ymax)
      ymax=ycur;
  }*/

  valtype clamp_scaled = clamp*coefexp*J0;
  if(1){
    for(int i=0;i<szi;i++){
			if(test_NaN(dx_ar[i])){
        return message(vblMESS1,0,"NaN value in solution\n");
			  //message(vblMESS1,0,"NaN value in solution\n");
				//dx_ar[i]= 0.;
			}
      // see PC1D Help -> Numerical Method, eq. A.21
      valtype denum=1+fabs(dx_ar[i]*clamp_scaled);
//      if(meth.&rpPsi)
      dx_ar[i]/=denum;
//      valtype ycur=fabs(ls.get_v()[i]);
    }
  }
  else if(mc.meth.nr==3 && mc.snum==2){
    for(int i=0;i<mc.sz[0];i++){
      if(test_NaN(dx_ar[i]))
        return message(vblMESS1,0,"NaN value in solution\n");
      // see PC1D Help -> Numerical Method, eq. A.21
      valtype denum=1+fabs(dx_ar[i]*clamp_scaled);
      dx_ar[i]/=denum;
      for(int car=0;car<mc.fstep[1];car++){
        if(test_NaN(dx_ar[mc.sz[0]+i+car*mc.sz[1]]))
          return message(vblMESS1,0,"NaN value in solution\n");
        valtype denum=1+fabs(dx_ar[mc.sz[0]+i+car*mc.sz[1]]*clamp_scaled);
        dx_ar[mc.sz[0]+i+car*mc.sz[1]]/=denum;
/*        // see PC1D Help -> Numerical Method, eq. A.22
        valtype dfi=dx_ar[sz[0]+i+car*sz[1]];
        valtype dif=dfi-dx_ar[i];
        denum=1+fabs(dif*clamp_1);
        dx_ar[sz[0]+i+car*sz[1]]=dx_ar[i]+dif/denum;*/
      }
    }
  }

  valtype dx_mult=1;
  int rot=0;
//  if(repr==rpPsi)
//    dx_mult=.1;
  
  int res_ind=0;
  valtype residue2[3]={0,0,0}; // residue corresponding to x - dx_mult*dx (calculating below)
  valtype dx_factors[3] = {0,1,0}; // values for dx_mult
  bool chooseSide = false;

  szi=0;
  for(int i0=0;i0<mc.snum;i0++){
    for(int i=0;i<mc.fstep[i0]*mc.sz[i0];i++){
      x_ar[szi+i]=mc.x[i0][i]; // memorizing current value of x
    }
    szi+=mc.fstep[i0]*mc.sz[i0];
  }

  start(tid.jacobian_step);

  while(1){
    
    szi=0;
    for(int i0=0;i0<mc.snum;i0++){
      for(int i=0;i<mc.fstep[i0]*mc.sz[i0];i++){
        /* x = x - dx_mult*dx, where dx is calculated from J*dx = dy.
        dx_mult is initially 1 which corresponds to Newton's methos,
        see http://en.wikipedia.org/wiki/Newton-Rhapson_algorithm
        Sometimes residue corresponding to x - dx is higher then residue corresponding to x.
        In this case we try to find such a value 0 < dx_mult < 1, that residue corresponding to x - dx_mult*dx 
        will be lower than residue corresponding to x
        */
        mc.x[i0][i]=x_ar[szi+i]-dx_mult*dx_ar[szi+i];
      }
      szi+=mc.fstep[i0]*mc.sz[i0];
    }

/*
    FILE *fxx30=fopen("xx30.txt","w");
    FILE *fx30=fopen("x30.txt","w");
    FILE *fx=fopen("x.txt","w");
    FILE *fy=fopen("y.txt","w");

    szi=0;
    for(int i0=0;i0<snum;i0++){
      for(int i=0;i<fstep[i0]*sz[i0];i++){
        fprintf(fxx30,"%g\n",xx30[szi+i]);
        fprintf(fx30,"%g\n",x30[szi+i]); // delta x
        fprintf(fx,"%g\n",x[i0][i]);
        fprintf(fy,"%g\n",ls.get_v()[szi+i]);
      }
      szi+=fstep[i0]*sz[i0];
    }

    fclose(fxx30);
    fclose(fx30);
    fclose(fx);
    fclose(fy);
*/
    // filling external nodes of each block with interpolated values from bordering blocks
    // current values at external nodes are memorized (transfer_values is called with '1' flag value)
    start(tid.jacobian_transfer_values);
    transfer_values(fix_trans,mc.meth.xreprs,mc.srp,1);
    stop(tid.jacobian_transfer_values);

#ifdef LIMIT_RESIDUE_GROWTH
    valtype rsdFull2 = calc_y(meth_cont_t(*this, method_t("peh")), y_ar);
#endif
    residue2[res_ind] = rsd2 = calc_y(mc,y_ar);

    if(test_NaN(rsd2))
      return message(vblMESS1,-1,"NaN value while calculating residue\n");
#ifdef CALC_OPTIMAL_DAMPING
    if (rot == 0) {
      if (rsd2 < rsd/2) {
        // If the residue decreases at an acceptable rate without using the damping
        // factor, then we won't spend time on finding an optimal value for the factor
        fix_trans.reset_previous_value(); // forget values at external nodes that were memorized during transfer_values
        break;
      }
      else {
        // Heuristically estimate the first value for the damping factor
        // assuming that the residue has a parabolic dependency on the factor
        // and that the minimal residue equals zero
        residue2[1] = rsd2;
        dx_mult = max(0.1, 1-1/(1+sqrt(residue2[0]/residue2[1])));
      }
    }
    else {
      if (chooseSide) {
        chooseSide = false;
        if (rsd2 > residue2[1]) {
          swap(residue2[1], rsd2);
          swap(dx_factors[1], dx_mult);
        }
        if (dx_mult > dx_factors[1]) {
          residue2[0] = residue2[1];
          dx_factors[0] = dx_factors[1];
          residue2[1] = residue2[2];
          dx_factors[1] = dx_factors[2];
        }
      }

      if (dx_mult < 1e-20 || (dx_mult - dx_factors[0])/dx_mult < 1e-3 ||
          (dx_factors[1] - dx_mult)/dx_factors[1] < 1e-3) {
        // The last guess for the damping factor didn't differ significantly from
        // the previous one; we'll assume it to be optimal
        fix_trans.reset_previous_value(); // forget values at external nodes that were memorized during transfer_values
        break;
      }
      valtype r1 = min(residue2[0], residue2[1]);
      valtype r2 = max(residue2[0], residue2[1]);
      if (rsd2 >= r1 + (r2 - r1)*pow(
          (dx_factors[residue2[0] > residue2[1] ? 1 : 0] - dx_mult) /
          (dx_factors[1] - dx_factors[0]), 2)) {
        // The residue depends on the factor non-parabolically, thus we choose
        // the next point in a way similar to the first iteration (see above).
        // We always choose the interval closer to 0
        residue2[1] = rsd2;
        dx_factors[1] = dx_mult;
        dx_mult = dx_factors[0] + (dx_factors[1] - dx_factors[0]) *
          max(0.1, min(0.9, 1-1/(1+sqrt(residue2[0]/residue2[1]))));
      }
      else {
        // The dependency is likely to be parabolic
        valtype a1 = dx_factors[0]*(residue2[1] - rsd2);
        valtype a2 = dx_factors[1]*(rsd2 - residue2[0]);
        valtype a3 = dx_mult*(residue2[0] - residue2[1]);
        valtype d = (a1*dx_factors[0] + a2*dx_factors[1] + a3*dx_mult)/(2*(a1 + a2 + a3));
        residue2[2] = residue2[1];
        dx_factors[2] = dx_factors[1];
        residue2[1] = rsd2;
        dx_factors[1] = dx_mult;
        dx_mult = d;
        chooseSide = true;
      }
    }
#else
    if(!res_ind){
#ifdef LIMIT_RESIDUE_GROWTH
      if(residue2[res_ind] <= rsd && rsdFull2 <= 10*rsdFull){ // we found such a dx_mult that residue for x - dx_mult*dx is lower then residue for x
#else
      if(residue2[res_ind]<=rsd){// && rsdFull2 <= 10*rsdFull){ // we found such a dx_mult that residue for x - dx_mult*dx is lower then residue for x
#endif
        fix_trans.reset_previous_value(); // forget values at external nodes that were memorized during transfer_values
        break; // leave 'dx_mult-decreasing' cycle
        // if break is commented then we try to find if there is some smaller dx_mult that leads to smaller residue
        res_ind++;
      }
    }
    else if(res_ind==1){
      if(residue2[1]>=residue2[0]){
        res_ind++;
        dx_mult*=4;
      }
      else
        residue2[0]=residue2[1];
    }
    else if(res_ind==2){
      fix_trans.reset_previous_value();
      break;
    }
#endif
    szi=0; // restoring x from recorded initial value
    for(int i0=0;i0<mc.snum;i0++){
      for(int i=0;i<mc.fstep[i0]*mc.sz[i0];i++){
        mc.x[i0][i]=x_ar[szi+i];
      }
      szi+=mc.fstep[i0]*mc.sz[i0];
    }

    transfer_values(fix_trans,mc.meth.xreprs,mc.srp,-1); // restore memorized values at external points ('-1' flag value)
    fix_trans.reset_previous_value(); // forget these memorized values

    rot++;

#ifndef CALC_OPTIMAL_DAMPING
    dx_mult/=2; // decreasing dx_mult
#endif

    //if(dx_mult==0){
    if(dx_mult<1e-4){
      residue2[0]=calc_y(mc,y_ar);
      break;
    }
    
  }
  stop(tid.jacobian_step);

  if(rot){
    message(vblMESS1,0,"%d ",rot);
  }

  //rsd2=residue2[0];

  // if the new residue differs just a little from the old one, then we assume that we already reach minimum
  return (1-rsd2/rsd > (rsd2 > residue_threshhold ? stop_ratio : small_stop_ratio)) ? 1 : 0;
}


template<class block_t>
int scBlockContainer<block_t>::AddMesh(mngarg<block_t> ptr, int lev){
  sc_bdescr_t *bd=new sc_bdescr_t(ptr,lev);
  bd->id=blocks.size();
  blocks.push_back(bd);
  return (int)blocks.size()-1;
}

template<class block_t>
void scBlockContainer<block_t>::ShareMediumRegions(const scRegionMediumMap& meds) {
  media_regions = &meds;

  if(dump)
    DumpOther(media_regions, true);

  vector<const scAbstractMedium*> media = media_regions->get_media();
  for(vector<const scAbstractMedium*>::const_iterator it = media.begin(), e = media.end(); it != e; ++it) {
    // If medium parameters depend on variables (carrier concentration etc.)
    // then we have to use general formulas for boundary conditions
    if(!(*it)->getDefaultContext().empty())
      use_drc[0] = use_drc[1] = 0;
  }
}

template<class block_t>
int scBlockContainer<block_t>::AddBoundaryRegion(mngarg<scBoundaryRegion> bond, int level){
   // to m/s
  for(int i=0;i<2;i++)
    if(bond->sr[i]!=VEC_INFTY) bond->sr[i]*=1e-2;
  multimap<int,mngptr<scBoundaryRegion> >::iterator it = boundaries.insert(make_pair(level, bond.ptr()));
  it->second.second=bond.second;
  if(dump)
    DumpOther(bond.ptr(),true);
  return (int)boundaries.size()-1;
}

template<class block_t>
int scBlockContainer<block_t>::ApplyVoltage(valtype psi){
  methi=jacstop=ch_it=0;
  valtype currVoltage = 0;
  int num=0;
  for (multimap<int, mngptr<scBoundaryRegion> >::const_reverse_iterator mit=boundaries.rbegin(), me=boundaries.rend(); mit!=me; ++mit){
    scBoundaryRegion &mr=*(mit->second);
    if(mr.appl_volt){
      // assume that the voltage is the same on all boundaries where it is applied
      currVoltage = mr.psi;

      mr.psi=psi;
      num++;
    }
  }
  if(!num)return 0;

  if (useLinearEstimate) {
    const size_t NVars = 3;
    valtype* currVars[NVars] = {psi_ptr, fermi[0], fermi[1]};
    if (prevVars.empty()) {
      prevVars.resize(NVars); // 3 variables: psi, fermiN, fermiP
      for(vector<vector<valtype> >::size_type nvar = 0; nvar < prevVars.size(); nvar++)
        prevVars[nvar].assign(currVars[nvar], currVars[nvar]+SZ);
    }
    else {
      for(vector<vector<valtype> >::size_type nvar = 0; nvar < prevVars.size(); nvar++) {
        for(vector<valtype>::size_type i = 0; i < prevVars[nvar].size(); i++) {
          valtype prev = currVars[nvar][i];
          currVars[nvar][i] += (1 + estimationNoise*(2.*rand()/RAND_MAX - 1)) *
            (currVars[nvar][i] - prevVars[nvar][i]) *
            (psi - currVoltage) / (currVoltage - prevVoltage);
          prevVars[nvar][i] = prev;
        }
      }
      prevVoltage = currVoltage;
    }
  }

  int upd_flag=scIniVoltage;
  for (multimap<int, mngptr<scGenerationRegion> >::const_reverse_iterator git=gens.rbegin(), ge=gens.rend(); git!=ge; ++git){
    int res=git->second->gen->change_voltage(psi);
    if(res>0) // we update generation array if it depends on voltage
      upd_flag|=scIniGen;
    else if(res<0)
      return -1;
  }

  return analyze(upd_flag);
}

template<class block_t>
int scBlockContainer<block_t>::AddGenerationRegion(mngarg<scGenerationRegion> gen,int level){
  multimap<int,mngptr<scGenerationRegion> >::iterator it = gens.insert(make_pair(level, gen.ptr()));
  it->second.second=gen.second;
  return (int)gens.size()-1;
}

template class scBlockContainer<>;
