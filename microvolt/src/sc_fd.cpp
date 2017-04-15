#include <algorithm>

#include "sc_fd.h"
#include "grid.hpp"
#include "data_flags.h"
#include "string_utils.h"

#include "physconst.h"


RectBlock::RectBlock(const Box &B, const int *ndim):
rbox(B),extsz(0),offset(0), 
psi_ptr(NULL),medc(),ch(NULL),gen(NULL),bcond(NULL),fixed(NULL),y_ar(NULL),
ls_pc2(NULL),ls_p(NULL),ls_c2(NULL),ls_pc(NULL),ls_c(NULL),
recomb_ptr(NULL),medh(),r_arr(NULL),
delta_x(0),T(0),J0(0),l0(0),n0(0),lu(0),coefpsi(0),coefexp(0),coefc(0),coefr(0),
SZ(0),grad_split(1.,1.,1.){

  for(int car=0;car<2;car++){
    conc[car]=NULL;
    fermi[car]=NULL;
    J[car]=NULL;
  }

  int dnc=0;
  fd=-1;
  for(int i=2;i>=0;i--){
    dx_arr[i]=NULL;

    if(ndim[i]>1){
      if(fd<0)
        fd=i;
      N[i]=ndim[i]+1;
      dim_type[i]=0;
      dnf[i]=dnc;
      dnb[dnc]=i;
      dnc++;
    }
    else{
      N[i]=1;
      dim_type[i]=-1;
      dnf[i]=-1;
    }
  }

  dim=dnc;

  for(;dnc<3;dnc++)
    dnb[dnc]=-1;

  grd.init(rbox.get_p1(),rbox.get_p2(),0,N);
}

RectBlock::~RectBlock(){
  for(int i=0;i<3;i++)
    delete[] dx_arr[i];

  delete[] r_arr;
  delete[] recomb_ptr;
}

int RectBlock::reset(valtype rad_axis){

  try{
    recomb_ptr = new valtype[SZ];

    for(int i=2;i>=0;i--){
      if(dim_type[i]>=0)
        dx_arr[i] = new valtype[N[i]+1];
      if(dim_type[i]==1)
        r_arr = new valtype[N[i]];
    }
  }
  catch(...){
    ::pmessage(vblMESS1,0,"Not enough memory for specified mesh(es)\n");
    comm->set_exit();
  }
  if(comm->check_exit())return -1;

  for(int i=0;i<SZ;i++){
    recomb_ptr[i]=0;
  }

  for(int di=2;di>=0;di--){
    if(dim_type[di]<0)
      continue;
    int ic[3]={0,0,0};
    valtype x1=get_position(ic[0],ic[1],ic[2])[di];
    if(dim_type[di]==1)
      r_arr[0]=(x1-rad_axis)/(l0/lu);
    for(int i=1;i<N[di];i++){
      ic[di]++;
      valtype x2=get_position(ic[0],ic[1],ic[2])[di];
      dx_arr[di][i]=(x2-x1)/(l0/lu);
      if(dx_arr[di][i]==0) // two mesh nodes has the same position
        return LOGERR(-1,"RectBlock::reset: two mesh nodes have the same position\n",0);
      if(dim_type[di]==1)
        r_arr[i]=(x2-rad_axis)/(l0/lu);
      x1=x2;
    }
    dx_arr[di][0]=0;
    dx_arr[di][N[di]]=0;
//    dx_arr[di][0]=dx_arr[di][1];
//    dx_arr[di][N[di]]=dx_arr[di][N[di]-1];
  }

  return SZ*(sizeof(local_t)*2*dim+sizeof(valtype));
}

int RectBlock::find_point(const node_t &p,const Vector_3 &nvect, Vector_3 &point){
  Vector_3 c=p.GetCenter();
  int found=-1;
  valtype frp=VEC_INFTY;
  for(Vector_3 *v=p.adj_begin(), *ve=p.adj_end();v!=ve;v++){
    Vector_3 nbr=*v;
    Plane_3 pl(nbr-c,nbr); // mesh surface
    Vector_3 s2; // intersection between normal to boundary and mesh surface
    valtype fr2=pl.TestRay(c,nvect,&s2);
    if(fr2>=0 && fr2<frp){ // find closest intersection between normal to boundary and mesh surface
      frp=fr2;
      point=s2;
      found=1;
    }
  }
  return found;
}

int RectBlock::organize_memory(){

  // 1.assign memory index to each point
  SZ=0;
  for(memmap_t::iterator it=memmap.begin(),e=memmap.end();it!=e;++it)
    it->second.first=SZ++; // memory index

  // 2. filling slices and b_slices
  int mind=-1; // for memory indeces
  int memszi=0;
  slice_t slice;

  for(memmap_t::iterator it=memmap.begin(),e=memmap.end();it!=e;++it){
    int mind1=it->second.first; // memory index
    int fl=it->second.second; // mesh point classification
    int pind1=it->first; // geometry index

    if(!(fl&mInt)){ // this is boundary point, put to b_slices
      slice_t b_slice;
      b_slice.oind=mind1;
      b_slice.gind=pind1;

      grd.unpack_ind(pind1,b_slice.ic[0],b_slice.ic[1],b_slice.ic[2]);
      int sh[3]={0,0,0};
      for(int di=2;di>=0;di--){ // find adjacent nodes memory indices
        if(dim_type[di]<0)
          continue;
        for(int ni=0;ni<2;ni++){
          int mdind = -1; // index of adjacent node
          int dn=dnum(di,ni);
          sh[di]=ni?1:-1;
          int exc=b_slice.ic[di]+sh[di];
          if(exc>=0 && exc<N[di]){
            int dind=grd.pack_ind(b_slice.ic[0]+sh[0],b_slice.ic[1]+sh[1],b_slice.ic[2]+sh[2]);
            memmap_t::iterator mdit=memmap.find(dind);
            if(mdit!=memmap.end()){
                mdind = memmap[dind].first; // adjacent node is stored
            }
          }
          b_slice.iind[dn]=mdind;
        }
        sh[di]=0;
      }
      b_slice.sz=1; // boundary slice size is 1
      b_slices.push_back(b_slice);
      continue;
    }
    
    if(mind1-mind>1){

      if(mind>=0){ // previous slice is finished, put it to slices
        slice.sz=memszi;
        slices.push_back(slice);
      }
      memszi=1; // start new slice

      slice.oind=mind1;
      slice.gind=pind1;

      grd.unpack_ind(pind1,slice.ic[0],slice.ic[1],slice.ic[2]);
      int sh[3]={0,0,0};
      for(int di=2;di>=0;di--){ // find adjacent nodes memory indices
        if(dim_type[di]<0)
          continue;
        for(int ni=0;ni<2;ni++){
          int dn=dnum(di,ni);
          sh[di]=ni?1:-1;
          int dind=grd.pack_ind(slice.ic[0]+sh[0],slice.ic[1]+sh[1],slice.ic[2]+sh[2]);
          int mdind=memmap[dind].first;
          slice.iind[dn]=mdind;
        }
        sh[di]=0;
      }
    }
    else
      memszi++;
    mind=mind1;
  }

  if(memszi){ // put last slice to slices
    slice.sz=memszi;
    slices.push_back(slice);
  }
  return 1;
}

// TODO: eliminate duplicate code from scBlockContainer<>::construct_local_media()
int RectBlock::set_local_media(const scRegionMediumMap* mr_map) {
  vector<Vector_3> points(2*dim*SZ); // space position of half-step mesh nodes
  for(node_it pit = begin_bound(), pe = end_bound(); pit != pe; ++pit) {
    const node_t& p = *pit;
    int pi = 2*dim*p.get_ind();
    for(Vector_3 *v = p.cur_begin(), *ve = p.cur_end(); v != ve; v++, pi++)
//      points[pi] = (*pit).GetCenter(); // this is bug !
      points[pi] = *v; // this is correct and agrees with dnf[]
  }

  for(node_it pit = begin_int(), pe = end_int(); pit != pe; ++pit) {
    const node_t& p = *pit;
    int pi = 2*dim*p.get_ind();
    for(Vector_3 *v = p.cur_begin(), *ve = p.cur_end(); v != ve; v++, pi++)
//      points[pi] = (*pit).GetCenter(); // this is bug !
      points[pi] = *v;
  }
  return medh.build(*mr_map, points.begin(), points.end());
}

int RectBlock::load_ini_values(virt_unary_function<const Vector_3 &, vec_type> *table, int repr){
  valtype *ptr = get_pointer(*this,repr);
  for(memmap_t::iterator it=memmap.begin(),e=memmap.end();it!=e;++it){
    int gind=it->first;
    int mind=it->second.first;
//    int flag=it->second.second;
    int ic[3];
    grd.unpack_ind(gind,ic[0],ic[1],ic[2]);
    Vector_3 pos=get_position(ic[0],ic[1],ic[2]);
    valtype val=(*table)(pos);
    ptr[mind]=val/J0;
  }
  for(int i=0;i<SZ;i++)
    set_repres(*this,rpConc,i);
  update_J();
  return 1;
}

int RectBlock::fill_adjacent(int mind,int gind,vector<int> &n,vector<int> &nsz,int flag,int mflag){

  int ic[3];
  grd.unpack_ind(gind,ic[0],ic[1],ic[2]);

  int fsz=0;
  int nsorted[vert_max];

  int sh[3]={0,0,0};
  for(int di=2;di>=0;di--){
    if(dim_type[di]<0)
      continue;
    for(int ni=0;ni<2;ni++){
      sh[di]=ni?1:-1;
      int exc=ic[di]+sh[di];
      if(exc>=0 && exc<N[di]){
        int iind=grd.pack_ind(ic[0]+sh[0],ic[1]+sh[1],ic[2]+sh[2]);
        if(memmap.find(iind)!=memmap.end()){
          pair<int,int> finded=memmap[iind];
          if(finded.second&mflag)
            nsorted[fsz++]=finded.first;
          else if(flag&8)
            return -1; // node does not belong to mflag filter
//            nsorted[fsz++]=-1;
        }
        else if(flag&0x10)
          nsorted[fsz++]=-1; // there is no such node in memory array
      }
      else if(flag&0x10)
        nsorted[fsz++]=-1; // there is no such node in memory array
    }
    sh[di]=0;
  }
  if(flag&1)
    nsorted[fsz++]=mind;
  if(!(flag&4))
    sort(nsorted,nsorted+fsz);
  for(int nsi=0;nsi<fsz;nsi++){
    int curind=nsorted[nsi];
    if(curind>=0 && flag&2)
      curind+=offset;
    n.push_back(curind);
  }
  nsz.push_back(fsz);

  return fsz;
}

int RectBlock::fill_adjacent_internal(int ind,vector<int> &n,vector<int> &nsz,int flag,int mflag,int *gind_){

  pair<vector<int> *,vector<int> *> p(&n,&nsz);
  for(size_t i=0;i<slices.size();i++){
    int mno=slices[i].oind;
    int msz=slices[i].sz;
    if(ind>=mno && ind<mno+msz){
      int gind=slices[i].gind+ind-mno;
      int sz=fill_adjacent(ind,gind,n,nsz,flag,mflag);
      if(gind_)
        *gind_=gind;
      return sz;
    }
  }
  return 0;
}

int RectBlock::fill_adjacent_border(int ind,vector<int> &n,vector<int> &nsz,int flag,int mflag,int *gind_){

  pair<vector<int> *,vector<int> *> p(&n,&nsz);
  for(size_t i=0;i<b_slices.size();i++){
    if(ind==b_slices[i].oind){
      int gind=b_slices[i].gind;
      int sz=fill_adjacent(ind,gind,n,nsz,flag,mflag);
      if(gind_)
        *gind_=gind;
      return sz;
    }
  }
  return 0;
}

int RectBlock::create_interpolation(const Vector_3 &pos, int *gind, int *mind, valtype *coeff, vector<Vector_3> *nonlocal, int incompl){
  int n=grd.GetCoeff(pos,gind,coeff,0,nonlocal);
  int num=0;
  valtype coef_div=0;
  for(int i=0;i<n;i++){
    if(gind[i]>=0){
      memmap_t::iterator mit=memmap.find(gind[i]);
      if(mit!=memmap.end()){
        int fl=mit->second.second;
        if(fl&mTrans){ // covered by other higher level mesh
          int ic[3];
          grd.unpack_ind(gind[i],ic[0],ic[1],ic[2]);
          Vector_3 nl=grd.Position(ic[0],ic[1],ic[2]);
          nonlocal->push_back(nl);
          mind[num]=gind[num]=-int((nonlocal->size()));
        }
        else{
          mind[num]=mit->second.first;
          gind[num]=gind[i];
        }
        coeff[num]=coeff[i];
        coef_div+=coeff[i];
        num++;
      }
      else{
        if(incompl){
/*          int ic[3];
          grd.unpack_ind(gind[i],ic[0],ic[1],ic[2]);
          Vector_3 nl=grd.Position(ic[0],ic[1],ic[2]);
          nonlocal->push_back(nl);
          mind[num]=gind[num]=-int((nonlocal->size()));
          coeff[num]=coeff[i];
          coef_div+=coeff[i];
          num++;*/
        }
        // in the organizing mesh memory layout stage there are cases when lower level mesh 
        // tries to interpolate boundary node of higher level mesh. 
        // because lower level mesh resolves boundary surface less accurate, 
        // it could not have enough nodes for this interpolation. 
        // therefore using nonlocal vector here can lead to infinte recursion.
        // to avoid this problem, we accept that currently this interpolation is unavailable,
        // and boundary node of higher level mesh will be marked as internal.
        else
          return -1;
      }
    }
    else{ // nonlocal is filled in grd.GetCoeff... however, Iam not sure that code will be executed here
      mind[num]=gind[num]=gind[i];
      coeff[num]=coeff[i];
      coef_div+=coeff[i];
      num++;
    }
  }
  for(int i=0;i<num;i++){
    coeff[i]/=coef_div;
  }
  return num;
}

int RectBlock::create_interpolation(const pair<const Vector_3,const Vector_3> &arg, 
int *gind, int *mind, valtype *coeff, vector<pair<Vector_3,Vector_3> > *nonlocal, int incompl){

  const Vector_3 &pos = arg.first;
  const Vector_3 &dir = arg.second;

  int gind0[vert_max], mind0[vert_max];
  valtype coeff0[vert_max];
  vector<Vector_3> nonlocal0;
  int n=create_interpolation(pos,gind0,mind0,coeff0,&nonlocal0,incompl);
  if(n<=0)
    return n;
  int nj=0;
  for(int di=2;di>=0;di--){
    if(dim_type[di]<0)
      continue;
    if(!dir[di])
      continue;
    for(int i=0;i<n;i++){
      if(gind0[i]>=0){
        gind[nj]=gind0[i];
        mind[nj]=mind0[i] + dnf[di]*SZ;
        coeff[nj]=coeff0[i]*dir[di];
        nj++;
      }
    }
  }
  int i1=0;
  for(int i=0;i<n;i++){
    if(gind0[i]<0){
      nonlocal->push_back(make_pair(nonlocal0[i1],dir));
      mind[nj]=gind[nj]=-int((nonlocal->size()));
      coeff[nj]=coeff0[i];
      nj++;
      i1++;
    }
  }

  return nj;
}

//#define OLD_CALCJ
#ifdef OLD_CALCJ
void RectBlock::calcJ(const slice_t &slice, int shift, const valtype* schemeCoeffs, Vector_3* J, int boundary) const {

  // in boundary case always sz = 1, index = 0

  int oind_i = slice.oind + shift;

  for(int car = 0; car < 2; car++)
    J[car] = 0;

  for(int di = 2; di >= 0; di--) {
    if(dim_type[di] < 0)
      continue;

    int ic_i=slice.ic[di];
    if(di==fd)
      ic_i+=shift;

    for(int ni = 0; ni < 2; ni++){
      int dn = dnum(di, ni);
      if(!schemeCoeffs[dn] || (boundary && slice.iind[dn] < 0)) // CHECK the last condition
        continue;

      int idn = dnum(di, 1-ni);

      int sign = ni ? 1: -1;
      int iind_i = slice.iind[dn] + shift;

      int grad = medc[oind_i].if_gradual(medc[iind_i]);

      valtype cfermi = medc[oind_i].get_ifermi()/J0;
      valtype c1fermi = medc[iind_i].get_ifermi()/J0;
      valtype hfermi  = medh[2*dim*oind_i + dn].get_ifermi()/J0;

      valtype cni = medc[oind_i].get_ni()/n0;
      valtype c1ni = medc[iind_i].get_ni()/n0;
      valtype hni  = medh[2*dim*oind_i + dn].get_ni()/n0;

      valtype psidif = psi_ptr[iind_i] - psi_ptr[oind_i];
      valtype ifermidif = c1fermi - cfermi;
      valtype logndif = c1ni!=cni ? (std::log(c1ni) - std::log(cni))/coefexp : 0;

      if (ifermidif == 0 && logndif == 0)
        grad = 1;

      for(int car = 0; car < 2; car++) {
        int csn = car ? 1 : -1;

        valtype dif = psidif - (grad ? ifermidif + csn*logndif : 0);
        // TODO: additional array to store exp(psi) to exlude calculation exponenta from this function
        valtype difexp = exp(dif*coefexp);
        valtype B[2];
        B[0]=accomp(difexp,1.) ? 1 :dif/(difexp-1); // Bernoulli function, B(x) = x/(e^x - 1)
        B[1]=B[0]*difexp; // B(-x) = B(x)*e^x

        valtype val;

        if(grad) // gradual
          val = B[car]*conc[car][iind_i]-B[1-car]*conc[car][oind_i]; // works only with gradual heterostructures
        else{ // abrupt
          /* in order to work with abrupt heterostructures, 
          coefficients r1, r2 corresponding to different fermi levels and intrisic concentrations
          are introduces to Scharfetter-Gummel scheme */
          valtype r1=(hni/c1ni)*exp(-csn*(c1fermi-hfermi)*coefexp);
          valtype r2=(hni/cni)*exp(-csn*(cfermi-hfermi)*coefexp);
          val = B[car]*conc[car][iind_i]*r1-B[1-car]*conc[car][oind_i]*r2;
        }

        valtype step = dx_arr[di][ic_i + ni]*l0; // mesh step in direction (di,ni)

        // Medium for diffusion and migration coefficients
        const local_t& cmed = medh[boundary ? 2*dim*iind_i + idn : 2*dim*oind_i + dn];

        // constant coefficient
        // We don't provide support for variable diffusion coefficients here since
        // using the optimized code for calcJ (below) is highly encouraged
        val *= cmed.get_D(car); 

        J[car][di] += -schemeCoeffs[dn]*csn*sign*val/step;
      }
    }
  }
}

#else

namespace {
  // Bernoulli function, B(x) = x/(e^x - 1)
  // B(-x) = B(x)*e^x
  inline valtype bernoulli(valtype dif, valtype difexp) {
    return accomp(difexp, 1.) ? 1 : dif/(difexp-1); // old
    //return abs(dif) < 0.003 ? 1 - dif/2*(1-dif/6) : dif/(difexp-1);
  }
}

# define ORGANIC_SCHG // use organic Scharfetter-Gummel scheme

void RectBlock::calcJ(const slice_t &slice, int shift, const valtype* schemeCoeffs, Vector_3* J, int boundary) const {

  // in boundary case always sz = 1, index = 0

  int oind_i = slice.oind + shift;

  for(int car = 0; car < 2; car++)
    J[car] = 0;

  valtype cfermi = medc[oind_i].get_ifermi()/J0;
  valtype cni = medc[oind_i].get_ni()/n0;

  // Preload default medium context.
  // We assume the set of variables is the same for the center and all the neighboring points
  scPointContext ctx = medc[oind_i].getDefaultContext();

  for(int di = 2; di >= 0; di--) { // dimensions
    if(dim_type[di] < 0)
      continue;

    int ic_i=slice.ic[di];
    if(di==fd) // fast dimension
      ic_i+=shift; // slice shift is in fast dimension

    for(int ni = 0; ni < 2; ni++){ // left - right
      int dn = dnum(di, ni);
      if(!schemeCoeffs[dn] || (boundary && slice.iind[dn] < 0)) // CHECK the last condition
        continue;

      int idn = dnum(di, 1-ni);

      int sign = ni ? 1: -1;
      int iind_i = slice.iind[dn] + shift;

      // Medium for diffusion and migration coefficients

      const local_t& cmed = medh[boundary ? 2*dim*iind_i + idn : 2*dim*oind_i + dn]; // half-step
      valtype c1fermi = medc[iind_i].get_ifermi()/J0;
      valtype c1ni = medc[iind_i].get_ni()/n0;

      CarrierParams oconc(conc[0][oind_i], conc[1][oind_i]);
      CarrierParams iconc(conc[0][iind_i], conc[1][iind_i]);

#ifndef ORGANIC_SCHG 

      valtype psidif = psi_ptr[iind_i] - psi_ptr[oind_i];
      valtype ifermidif = c1fermi - cfermi;
      valtype logndif = c1ni == cni ? 0 : (std::log(c1ni) - std::log(cni))/coefexp;

      int grad = (ifermidif == 0 && c1ni == cni) ? 1 :
        medc[oind_i].if_gradual(medc[iind_i]);

      valtype dif = psidif - (grad ? ifermidif - logndif : 0);
      valtype difexp = exp(dif*coefexp);
      valtype B[4];
      // for car=0
      B[0] = bernoulli(dif, difexp); // Bernoulli function, B(x) = x/(e^x - 1)
      B[2] = B[0]*difexp; // B(-x) = B(x)*e^x
      // for car=1
      B[3] = -B[0];
      B[1] = -B[2];

      if (!grad) {
#if !defined(OLD_J)
        CarrierParams ofermi(fermi[0][oind_i], fermi[1][oind_i]);
        CarrierParams ifermi(fermi[0][iind_i], fermi[1][iind_i]);
        oconc = cmed.calcConcentration(ofermi, psi_ptr[oind_i]);
        iconc = cmed.calcConcentration(ifermi, psi_ptr[iind_i]);
#else
        
//        dif -= 2*logndif; // same as dif = psidif - (grad ? ifermidif + logndif : 0);
//        difexp *= pow(cni/c1ni, 2);
//        B[3] = -bernoulli(dif, difexp); // Bernoulli function, B(x) = x/(e^x - 1)
//        B[1] = B[3]*difexp; // B(-x) = B(x)*e^x

        valtype hfermi  = medh[2*dim*oind_i + dn].get_ifermi()/J0;
        valtype hni  = medh[2*dim*oind_i + dn].get_ni()/n0;

        valtype r1c = hni/c1ni;
        valtype r2c = hni/cni;
        valtype r1e = exp(-(c1fermi-hfermi)*coefexp);
        valtype r2e = exp(-(cfermi-hfermi)*coefexp);
        B[0] *= r1c/r1e;
        B[1] *= r1c*r1e;
        B[2] *= r2c/r2e;
        B[3] *= r2c*r2e;
#endif
      }
      else if (c1ni != cni) {
        dif -= 2*logndif; // same as dif = psidif - (grad ? ifermidif + logndif : 0);
        difexp *= pow(cni/c1ni, 2);
        B[3] = -bernoulli(dif, difexp); // Bernoulli function, B(x) = x/(e^x - 1)
        B[1] = B[3]*difexp; // B(-x) = B(x)*e^x
      }

      valtype step = dx_arr[di][ic_i + ni]*l0; // mesh step in direction (di,ni)

#if !defined(OLD_J)
      CarrierParams mu, D;
      if (ctx.empty()) {
        // Use faster calls without the context parameter
        mu = cmed.get_mu();
        D = cmed.get_D();
      }
      else {
#ifndef NDEBUG
        std::set<scPointContext::VarID> vars = ctx.getVarList();
#endif

        if (ctx.contains(scPointContext::varPsi)) {
          valtype psi2 = (psi_ptr[oind_i] + psidif/2)*J0;
          ctx.set(scPointContext::varPsi, psi2);
#ifndef NDEBUG
          vars.erase(scPointContext::varPsi);
#endif
        }

        if (ctx.contains(scPointContext::varEField)) {
          // This is for 1D only!
          // TODO: calculate field in 2D and 3D here
          valtype field = abs(psidif*J0/step);
          ctx.set(scPointContext::varEField, field);
#ifndef NDEBUG
          vars.erase(scPointContext::varEField);
#endif
        }

        if (ctx.contains(scPointContext::varConcN)) {
          valtype n = conc[0][oind_i] + (conc[0][iind_i] - conc[0][oind_i])/(1 + sqrt(B[2]/B[0]));
          n *= n0;
          ctx.set(scPointContext::varConcN, n);
#ifndef NDEBUG
          vars.erase(scPointContext::varConcN);
#endif
        }

        if (ctx.contains(scPointContext::varConcP)) {
          valtype p = conc[1][oind_i] + (conc[1][iind_i] - conc[1][oind_i])/(1 + sqrt(B[1]/B[3]));
          p *= n0;
          ctx.set(scPointContext::varConcP, p);
#ifndef NDEBUG
          vars.erase(scPointContext::varConcP);
#endif
        }

        if (ctx.contains(scPointContext::varFermiN)) {
          valtype fermin = (fermi[0][oind_i] + fermi[0][iind_i])/2;
          fermin *= J0;
          ctx.set(scPointContext::varFermiN, fermin);
#ifndef NDEBUG
          vars.erase(scPointContext::varFermiN);
#endif
        }

        if (ctx.contains(scPointContext::varFermiP)) {
          valtype fermip = (fermi[1][oind_i] + fermi[1][iind_i])/2;
          fermip *= J0;
          ctx.set(scPointContext::varFermiP, fermip);
#ifndef NDEBUG
          vars.erase(scPointContext::varFermiP);
#endif
        }

        assert(vars.empty()); // assertion fails if the medium needs something we didn't provide

        mu = cmed.get_mu(ctx);
        D = cmed.get_D(ctx);
      }

      //if (mu[0] < 0 || mu[1] < 0 || D[0] < 0 || D[1] < 0) {
      //  if (mu[0] < 0 || mu[1] < 0)
      //    message(vblERR, -1, "Mobility coefficient is negative");
      //  if (D[0] < 0 || D[1] < 0)
      //    message(vblERR, -1, "Diffusion coefficient is negative");
      //  return;
      //}

      CarrierParams mukT = mu*phys::kB_eV*T;
      CarrierParams ratio = D/(D + mukT);
//      if (false) {
      if (!accomp(ratio[0], 0.5) || !accomp(ratio[1], 0.5)) { // accomp is time consuming
        // D and mu are not related by the classic Einstein formula
        // thus we need to recalculate the coefficients

        CarrierParams dif = psidif - ifermidif - negate_first(logndif);
        CarrierParams difmu = dif*mukT*coefexp;
        CarrierParams difc = difmu/D;
        CarrierParams difexp = exp(difc);
        const valtype minratio = 1/std::log(1e50);
        B[0] = ratio[0] < minratio ? -min(difmu[0], 0.) : D[0]*bernoulli(difc[0], difexp[0]);
        B[2] = ratio[0] < minratio ?  max(difmu[0], 0.) : B[0]*difexp[0];
        B[3] = ratio[1] < minratio ?  min(difmu[1], 0.) : -D[1]*bernoulli(difc[1], difexp[1]);
        B[1] = ratio[1] < minratio ? -max(difmu[1], 0.) : B[3]*difexp[1];
      }
      else {
        for(int i = 0; i < 4; i++)
          B[i] *= D[i%2];
      }

      CarrierParams val = CarrierParams(B[0], B[1])*iconc - CarrierParams(B[2], B[3])*oconc;
      for(int car = 0; car < 2; car++)
        J[car][di] += schemeCoeffs[dn]*sign*val[car]/step;


#else // old medium
      for(int car = 0; car < 2; car++) {
        valtype val = B[car]*conc[car][iind_i]-B[2+car]*conc[car][oind_i];
        val *= cmed.get_D(car); // constant coefficient
        J[car][di] += schemeCoeffs[dn]*sign*val/step;
      }
#endif

#else // organic Scharfetter-Gummel

      int grad = medc[oind_i].if_gradual(medc[iind_i]);

      valtype coef = J0*coefexp;
      
      valtype psidif = psi_ptr[iind_i] - psi_ptr[oind_i];

      CarrierParams ofermi(fermi[0][oind_i], fermi[1][oind_i]);
      CarrierParams ifermi(fermi[0][iind_i], fermi[1][iind_i]);

      oconc = cmed.calcConcentration(coef*ofermi, coef*psi_ptr[oind_i])/n0;
      iconc = cmed.calcConcentration(coef*ifermi, coef*psi_ptr[iind_i])/n0;


      CarrierParams ig3 = medc[iind_i].calcG3(coef*ifermi,coef*psi_ptr[iind_i]);
      CarrierParams og3 = medc[oind_i].calcG3(coef*ofermi,coef*psi_ptr[oind_i]);

      
      valtype step = dx_arr[di][ic_i + ni]*l0; // mesh step in direction (di,ni)

      CarrierParams iofermi[2] = { ifermi, ofermi};
      CarrierParams ioconc[2] = { iconc, oconc};
      valtype iopsi[2] = { psi_ptr[iind_i] , psi_ptr[oind_i] }; 
      
      CarrierParams iomu[2]; // to be calculated
#ifndef NDEBUG
      std::set<scPointContext::VarID> vars = ctx.getVarList(); // control set for debugging
#endif
      // field is always in the middlepoint:
      if (ctx.contains(scPointContext::varEField)) { // field is calculated in the middle
        // This is for 1D only!
        // TODO: calculate field in 2D and 3D here
        valtype field = abs(psidif*J0/step);
        ctx.set(scPointContext::varEField, field);
#ifndef NDEBUG
        vars.erase(scPointContext::varEField);
#endif
      }

      for(size_t i=0;i<2;i++){
        
        if (ctx.contains(scPointContext::varPsi)) {
          ctx.set(scPointContext::varPsi, iopsi[i]);
#ifndef NDEBUG
          vars.erase(scPointContext::varPsi);
#endif
        }

        if (ctx.contains(scPointContext::varConcN)) {
          ctx.set(scPointContext::varConcN, ioconc[i][0]*n0);
#ifndef NDEBUG
          vars.erase(scPointContext::varConcN);
#endif
        }

        if (ctx.contains(scPointContext::varConcP)) {
          ctx.set(scPointContext::varConcP, ioconc[i][1]*n0);
#ifndef NDEBUG
          vars.erase(scPointContext::varConcP);
#endif
        }

        if (ctx.contains(scPointContext::varFermiN)) {
          ctx.set(scPointContext::varFermiN, iofermi[i][0]);
#ifndef NDEBUG
          vars.erase(scPointContext::varFermiN);
#endif
        }

        if (ctx.contains(scPointContext::varFermiP)) {
          ctx.set(scPointContext::varFermiP, iofermi[i][1]);
#ifndef NDEBUG
          vars.erase(scPointContext::varFermiP);
#endif
        }

        assert(vars.empty()); // assertion fails if the medium needs something we didn't provide

        iomu[i] = cmed.get_mu(ctx);
        
      }

      CarrierParams muh(sqrt(iomu[0][0]*iomu[1][0]),sqrt(iomu[0][1]*iomu[1][1])); // geometric average for mue
      CarrierParams g3h(sqrt(ig3[0]*og3[0]),sqrt(ig3[1]*og3[1])); // geometric average for g3
# if 0 // new diff 
      CarrierParams difexp, val;
      for(int i =0 ; i<2 ;i++){
        difexp[i] = exp(psidif*coefexp/g3h[i]/2.);
        val[i] = muh[i]*phys::kB_eV*T;
      }
      CarrierParams conch;
      conch[0]= ioconc[0][0]*(1. - 1./(1.+difexp[0])) + ioconc[1][0]/(1.+difexp[0]);
      conch[1]= ioconc[1][1]*(1. - 1./(1.+difexp[1])) + ioconc[0][1]/(1.+difexp[1]);
      CarrierParams fermidif = iofermi[1]-iofermi[0];
      for(int car = 0; car < 2; car++)
        J[car][di] += schemeCoeffs[dn]*sign*val[car]*conch[car]*fermidif[car]/step;

# else

      CarrierParams val, difexp;
      for(int i =0 ; i<2 ;i++){
        difexp[i] = exp(psidif*coefexp/g3h[i]);
        val[i] = muh[i]*phys::kB_eV*T;
      }
      // electrons
      if(difexp[0]>1.){
        valtype  k = 1./difexp[0];
        val[0] *= psidif*(iconc[0]*k - oconc[0])/(1.-k); 
      }
      else{
        val[0] *= (iconc[0] - oconc[0]*difexp[0])*bernoulli(psidif,difexp[0]); 
      }
       
      // holes
      if(difexp[1]>1.){
        valtype  k = 1./difexp[1];
        val[1] *= psidif*(-iconc[1] + oconc[1]*k)/(1.-k); 
      }
      else{
        val[1] *= (-iconc[1]*difexp[1] + oconc[1])*bernoulli(psidif,difexp[1]); 
      }

      for(int car = 0; car < 2; car++)
        J[car][di] += schemeCoeffs[dn]*sign*val[car]/step;
# endif

#endif // organic Scharfetter-Gummel
    }
  }
}

#endif

int RectBlock::update_J(){

  for(int i=0;i<SZ;i++)
    set_repres(*this,rpConc,i);

  valtype coeff[2*dim_max];
  for(int dn = 0; dn < 2*dim_max; dn++)
    coeff[dn] = 0.5;
  Vector_3 Jlocal[2];

  for(node_it it=begin_int(),e=end_int();it!=e;++it){
    int loc;
    slice_t &slice=it.get_slice_position(loc);
    calcJ(slice,loc,coeff,Jlocal,0);

    for(int di = 2; di >= 0; di--) {
      if (dim_type[di] >= 0) {
        for(int car = 0; car < 2; car++)
          J[car][dnf[di]*SZ+slice.oind+loc] = Jlocal[car][di];
      }
    }
  }

  for(size_t i=0;i<b_slices.size();i++)
    update_bJ(i);

  return 1;
}

int RectBlock::update_bJ(int bj){

  valtype coeff[2*dim_max];
  Vector_3 Jlocal[2];

  for(int di = 2; di >= 0; di--) {
    if (dim_type[di] < 0)
      continue;

    int dvd = 0;
    for(int ni = 0; ni < 2; ni++)
      dvd += b_slices[bj].iind[dnum(di, ni)] >= 0;

    for(int ni = 0; ni < 2; ni++)
      coeff[dnum(di, ni)] = 1./dvd;
  }

  slice_t &slice=b_slices[bj];

  calcJ(slice, 0, coeff, Jlocal, 1);

  for(int di = 2; di >= 0; di--) {
    if (dim_type[di] >= 0) {
      for(int car = 0; car < 2; car++)
        J[car][dnf[di]*SZ + b_slices[bj].oind] = Jlocal[car][di];
    }
  }

  return 1;
}

scPointContext RectBlock::getContext(const slice_t& slice, int shift) const {
  int oind_i = slice.oind + shift;
  scPointContext ctx = medc[oind_i].getDefaultContext();

#ifndef NDEBUG
  std::set<scPointContext::VarID> vars = ctx.getVarList();
#endif

  if (ctx.contains(scPointContext::varPsi)) {
    ctx.set(scPointContext::varPsi, psi_ptr[oind_i]*J0);
#ifndef NDEBUG
    vars.erase(scPointContext::varPsi);
#endif
  }

  if (ctx.contains(scPointContext::varConcN)) {
    ctx.set(scPointContext::varConcN, conc[0][oind_i]*n0);
#ifndef NDEBUG
    vars.erase(scPointContext::varConcN);
#endif
  }

  if (ctx.contains(scPointContext::varConcP)) {
    ctx.set(scPointContext::varConcP, conc[1][oind_i]*n0);
#ifndef NDEBUG
    vars.erase(scPointContext::varConcP);
#endif
  }

  if (ctx.contains(scPointContext::varFermiN)) {
    ctx.set(scPointContext::varFermiN, fermi[0][oind_i]*J0);
#ifndef NDEBUG
    vars.erase(scPointContext::varFermiN);
#endif
  }

  if (ctx.contains(scPointContext::varFermiP)) {
    ctx.set(scPointContext::varFermiP, fermi[1][oind_i]*J0);
#ifndef NDEBUG
    vars.erase(scPointContext::varFermiP);
#endif
  }

  if (ctx.contains(scPointContext::varEField)) {
    Vector_3 field;

    for(int di = 2; di >= 0; di--) {
      if(dim_type[di] < 0)
        continue;

      int ic_i = slice.ic[di];
      if(di == fd)
        ic_i += shift;

      valtype stepWeight = 0;
      for(int ni = 0; ni < 2; ni++) {
        int dn = dnum(di, ni);
        if(slice.iind[dn] < 0)
          continue;

        int idn = dnum(di, 1-ni);

        int sign = ni ? 1 : -1;
        int iind_i = slice.iind[dn] + shift;

        valtype psidif = (psi_ptr[iind_i] - psi_ptr[oind_i])*J0;
        valtype invStep = 1/(dx_arr[di][ic_i + ni]*l0); // mesh step in direction (di,ni)

        // The field is calculated as a weighted average of the values on both sides of the point
        field[di] += (psidif*invStep - field[di]) * (invStep / (stepWeight += invStep));
      }
    }

    ctx.set(scPointContext::varEField, field.norm());
#ifndef NDEBUG
    vars.erase(scPointContext::varEField);
#endif
  }

  assert(vars.empty()); // assertion fails if the medium needs something we didn't provide

  return ctx;
}

void RectBlock::recomb_gen(const slice_t& slice, int shift, valtype &R, valtype &G) {
  int oind_i = slice.oind + shift;
  const local_t& med = medc[oind_i];
  if (med.getDefaultContext().empty())
    R = med.get_recomb(n0*conc[0][oind_i], n0*conc[1][oind_i]);
  else
    R = med.get_recomb(n0*conc[0][oind_i], n0*conc[1][oind_i], getContext(slice, shift));

  G = gen ? gen[oind_i] : 0;
}

void RectBlock::recomb_gen(int i,valtype &R,valtype &G){
  R = medc[i].get_recomb(n0*conc[0][i], n0*conc[1][i]);
  G = gen ? gen[i] : 0;
}

void RectBlock::fill_recomb(){
  for(int i=0;i<SZ;i++){
    valtype R,G;
    recomb_gen(i,R,G);
    recomb_ptr[i]=R;
  }
}

valtype RectBlock::Poisson(const slice_t &slice, int i, valtype *y, int bf, int carsz){

  int oind=slice.oind;
  const int *iind=slice.iind;
  const int *ic=slice.ic;

  // free charge concentration
  valtype free_charge=-conc[1][oind+i]+conc[0][oind+i];
  valtype total_charge=free_charge-ch[oind+i]; // total charge = free charge - dopants (N_D and N_A)
  valtype der=0; // \nabla (\varepsilon \nabla \psi)

  for(int di=2;di>=0;di--){

    if(dim_type[di]<0)
      continue;

    int ic_i=ic[di];
    if(di==fd)
      ic_i+=i;

    valtype r = dim_type[di]==1 ? r_arr[ic_i] : -1; // used only for radial coordinade

    if(bcond[oind+i]==3 && dim_type[di]==1){ // cylindrical coordinates, r = 0
      // A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012), formula (33)
      int lr = ic_i==N[di]-1; // if radial axis is located at the left from axis
      valtype hsum=dx_arr[di][ic_i+1-lr];
      for(int ni=0;ni<2;ni++){
        if(ni==lr)
          continue;
        int dn=dnum(di,ni);
//        int grad=medc[oind+i].if_gradual(medc[iind[dn]+i]);
        valtype epsilon=medh[2*dim*(oind+i)+dn].get_eps();
        // 2 directions, back and forward
        der+=4*epsilon*(psi_ptr[iind[dn]+i]-psi_ptr[oind+i])/(hsum*dx_arr[di][ic_i+ni]);
      }
      continue;
    }

//    int i0=dnum(di,0),i1=dnum(di,1);
//    int grad=medc[oind+i].if_gradual(medc[iind[i0]+i]) && medc[oind+i].if_gradual(medc[iind[i1]+i]);
    // code works wrong with abrupt change of epsilon (grad=0), 
    // for this case applying boundary conditions is necessary

    valtype hsum=(dx_arr[di][ic_i]+dx_arr[di][ic_i+1])/2;

    for(int ni=0;ni<2;ni++){
      int dn=dnum(di,ni);
      valtype epsilon=medh[2*dim*(oind+i)+dn].get_eps();
      der+=epsilon*(psi_ptr[iind[dn]+i]-psi_ptr[oind+i])/(hsum*dx_arr[di][ic_i+ni]);
    }
    if(dim_type[di]==1){ // cylindrical coordinates, adding 1/r * dpsi/dr
      // A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012), formula (31)
      valtype epsilon=medc[oind+i].get_eps();
      der+=epsilon*(psi_ptr[iind[dnum(di,1)]+i]-psi_ptr[iind[dnum(di,0)]+i])/(2*hsum*fabs(r));
    }

  }
  *y=total_charge-coefpsi*der;

  return (*y)*(*y);
}

valtype RectBlock::Continuity(const slice_t &slice, int i, valtype *y, int bf, int carsz){

  int oind=slice.oind;
//  const int *iind=slice.iind;
  const int *ic=slice.ic;

  valtype coeff[2*dim_max];
  Vector_3 Jdiff[2];

  valtype R,G; // recombination and generation
  recomb_gen(slice, i, R, G);
  //coefr = 0.;
  valtype RG = (R-G)*coefr; // rescaling

  valtype delta[2]={0,0};

  for(int di = 2; di >= 0; di--) {
    if(dim_type[di]<0)
      continue;

    int ic_i=ic[di];
    if(di==fd)
      ic_i+=i;

    // \Delta x_i^{avg} = (\Delta x_i + \Delta x_{i-1})/2
    valtype hsum;
    if(bcond[oind+i]==3 && dim_type[di]==1){ // cylindrical axis, r = 0
      int lr = ic_i==N[di]-1; // if radial axis is located at the left from axis
      hsum = dx_arr[di][ic_i+1-lr];
    }
    else
      hsum = (dx_arr[di][ic_i]+dx_arr[di][ic_i+1])/2;

    for(int ni = 0; ni < 2; ni++) {

      int dn = dnum(di, ni);
      int sign = ni ? 1 : -1;
      coeff[dn] = -sign/hsum; // current divergence (J_{1/2) - J_{-1/2))/hsum

      if (dim_type[di] == 1) { // radial finite-difference scheme
        if (bcond[oind+i]==3){ // cylindrical axis, r = 0
          // radius = 0, apply rectangular finite difference scheme
          // A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012), formula (33)
          int lr = ic_i==N[di]-1; // if radial axis is located at the left from axis
          if(ni==lr) 
            coeff[dn]=0;
          else 
            coeff[dn] *= 4;
        }
        else
          // radius != 0, add extra term proportional to 1/(2*radius)
          // A. Deinega and S. John, Comp. Phys. Comm. 183, 2128 (2012), formulae (32), (27)
          coeff[dn] -= 1/(2*r_arr[ic_i]);
      }
    }
  }

  calcJ(slice, i, coeff, Jdiff, 0);

  for(int di = 2; di >= 0; di--) {
    if (dim_type[di] >= 0) {
      for(int car = 0; car < 2; car++)
        delta[car] += (car ? 1 : -1) * Jdiff[car][di];
    }
  }

  valtype val=0;

  for(int car=0,cd=0;car<2;car++){
    int bt=0x1<<car;
    if(!(bf&bt))
      continue;
    delta[car]*=coefc*l0;

    valtype d1=delta[car]-RG;
    //valtype d2=log_reduce(delta[car]); //-log_reduce(R);

    y[cd*carsz]=d1;
    val+=d1*d1;
    cd++;
  }

  return val;
}

int RectBlock::Jacobian_node(int oind, const int *nind, int nind_sz, int shift, const meth_cont_t &mc){

  valtype ny[3*(2*dim_max+1)];
  bool tch[2*dim_max+1];

  for(int j=0;j<3*(2*dim+1);j++)
    ny[j]=0;

  for(int j=0;j<=2*dim;j++)
    tch[j]=0;

  for(int i0=0,szi=0;i0<mc.snum;i0++){

    for(int is=0;is<mc.fstep[i0];is++){

      if(mc.fx[i0][offset+mc.sz[i0]*is+oind+shift])
        continue;

      valtype prev=mc.x[i0][offset+mc.sz[i0]*is+oind+shift]; // shifting variable
      mc.x[i0][offset+mc.sz[i0]*is+oind+shift]+=delta_x;
      set_repres(*this,mc.srp,oind+shift);
      
      int szj=0;
      int szjf=0;
      for(int j0=0;j0<mc.snum;j0++){

        chosen_nodes(mc.gyrepr[j0],(int *)nind,tch,nind_sz,shift,(valtype *)(ny+(2*dim+1)*szjf));

        for(int js=0;js<mc.fstep[j0];js++){
          for(int j=0;j<nind_sz;j++){
            if(!tch[j])
              continue;
            valtype y2=ny[(szjf+js)*(2*dim+1)+j];
            valtype y1=mc.ls->v(offset+js*mc.sz[j0]+szj+nind[j]+shift);
            valtype delta_y=y2-y1;
            int ii=offset+mc.sz[i0]*is+szi+oind+shift;
            int jj=offset+js*mc.sz[j0]+szj+nind[j]+shift;
            if(mc.ls->set_m(ii,jj,delta_y/delta_x)<0)
              return -1;
          }
        }
        szj+=mc.fstep[j0]*mc.sz[j0];
        szjf+=mc.fstep[j0];
      }
      mc.x[i0][offset+mc.sz[i0]*is+oind+shift]=prev;
      set_repres(*this,mc.srp,oind+shift);
    }
    szi+=mc.fstep[i0]*mc.sz[i0];
  }
  return 1;
}

int RectBlock::Jacobian(const meth_cont_t &mc){

  for(int i=0;i<SZ;i++)
    set_repres(*this,rpConc,i);

  /*
  During Jacobian delta_y/delta_x calculation, we can add to the discretized equation 
  some function y_0 which does not depend on scalar varialbes x,
  sice this procedure doesn't influence on delta_y:
  delta_y = [y(x+delta_x) + y_0] - [y(x) + y_0] = y(x+delta_x) - y(x) 

  This trick might improve acuracy of calculated delta_y, if y_0 is chosen to be close to y(x).
  We can implement this trick using different generation term while calculating Jacobian
  All you need to do is just temporarily assign pointer gen to address of some specific array
  */

  valtype *gen_ptr=gen; // store pointer on actual generation array

//  gen=NULL; // turn off generation when calculating Jacobian

//  fill_recomb(); // fill recombination array
//  gen=recomb_ptr; // use generation equal to recombination when calculating Jacobian

  if(gen_ptr!=gen){ // calculate new y for internal mesh nodes
    for(int j0=0,szi=0;j0<mc.snum;j0++){
      if(mc.gyrepr[j0]&rpFermi){
        internal_nodes(mc.gyrepr[j0],mc.ls->get_v()+offset+szi);
      }
      szi+=mc.fstep[j0]*mc.sz[j0];
    }
  }


//  aggregate_it<node_it,bound_it> it;

  int loc0=INT_INFTY;
  int nind[2*dim_max+1];
  int nind_sz=0; // size of nind
  for(node_it it=begin_int(),e=end_int();it!=e;++it){

    int loc;
    slice_t &slice=it.get_slice_position(loc);

    if(loc<=loc0){
      loc0=loc;
      // we fill nind array in a ordered way
      nind_sz=0; // size of nind
      for(int i=dim-1;i>=0;i--){
        if(slice.iind[2*i]<0)
          continue;
        nind[nind_sz++]=slice.iind[2*i];
      }
      nind[nind_sz++]=slice.oind;
      for(int i=0;i<dim;i++){
        if(slice.iind[2*i+1]<0)
          continue;
        nind[nind_sz++]=slice.iind[2*i+1];
      }
    }

    if(Jacobian_node(slice.oind,nind,nind_sz,loc,mc)<0)
      return -1;
  }

  for(node_it it=begin_bound(),e=end_bound();it!=e;++it){
    int loc;
    slice_t &slice=it.get_slice_position(loc);

    if(loc<=loc0){
      loc0=loc;
      // we fill nind array in a ordered way
      nind_sz=0; // size of nind
      for(int i=dim-1;i>=0;i--){
        if(slice.iind[2*i]<0)
          continue;
        nind[nind_sz++]=slice.iind[2*i];
      }
      nind[nind_sz++]=slice.oind;
      for(int i=0;i<dim;i++){
        if(slice.iind[2*i+1]<0)
          continue;
        nind[nind_sz++]=slice.iind[2*i+1];
      }
    }

    if(Jacobian_node(slice.oind,nind,nind_sz,loc,mc)<0)
      return -1;
  }

  if(gen_ptr!=gen){
    gen=gen_ptr; // swithing back to generation
    for(int j0=0,szi=0;j0<mc.snum;j0++){ // swithcing back to old y
      if(mc.gyrepr[j0]&rpFermi){
        internal_nodes(mc.gyrepr[j0],mc.ls->get_v()+offset+szi);
      }
      szi+=mc.fstep[j0]*mc.sz[j0];
    }
  }

  return 1;
}

int RectBlock::DumpMesh(const string &suf, bool dif, int outtype){
  // If dif is true, then the values for the boundaries of each type
  // are dumped to separate files
  const int bcondList[] = { 0,   1,   2,   3,   4,   0};
  const string suffix[] = {"", "c", "d", "r", "p", "t"};

  string fnameBase = theConfig->GetOutputDir() + (name.empty() ? "m" : name);
  FILE *f = NULL;

  for(size_t i = 0; i < sizeof(bcondList)/sizeof(bcondList[0]); ++i) {
    if (i == 0 || dif) {
      if (f)
        fclose(f);
      string fname = fnameBase + suffix[i] + suf + ".d";
      f = fopen(fname.c_str(), "w");
      if (!f)
        return message(vblERR, -1, ("Failed to create output file " + fname).c_str());

      int result = Dump_header(f, outtype);
      if (result < 0)
        return result;
    }

    for(node_it it = i == 0 ? begin_int() : begin_bound(),
        e = i == 0 ? end_int() : end_bound();
        it != e; ++it) {
      int loc;
      const slice_t& slice = it.get_slice_position(loc);
      if(bcond[slice.oind + loc] == bcondList[i]) {
        int result = Dump_element(slice, loc, f, outtype);
        if (result < 0)
          return result;
      }
    }
    fprintf(f,"\n\n");
  }

  if (f)
    fclose(f);

  return 1;
}

int RectBlock::Dump_header(FILE *f, int outtype){

  fprintf(f,"# ");

  int argtype=0;
  for(int i=0;i<3;i++){
    if(dim_type[i]>=0)
      argtype|=argx<<i;
  }
  vector<string> columns=out_names(argtype, outtype);
  for(size_t i=0;i<columns.size();i++)
    fprintf(f,"%d-%s\t",int(i+1),columns[i].c_str());

  fprintf(f,"\n");

  return 1;
}

int RectBlock::Dump_element(const slice_t &slice, int shift, FILE *f, int outtype){
  int oind_i = slice.oind + shift;

  int icc[3];
  for(int i = 0; i < 3; i++)
    icc[i] = slice.ic[i];
  icc[fd] += shift;

  Vector_3 pos=get_position(icc[0],icc[1],icc[2]);

  for(int j=0;j<3;j++)
    if(dim_type[j]>=0)fprintf(f,"%g\t",pos[j]);

  if(outtype&outPsi){
    fprintf(f,"%g\t",psi_ptr[oind_i]*J0-global_ifermi);
  }
  if(outtype&outFermi){
    for(int car=0;car<2;car++){
      fprintf(f,"%g\t",-fermi[car][oind_i]*J0+global_ifermi);
    }      
  }
  if(outtype&outBands){
    valtype c=-psi_ptr[oind_i]*J0-medc[oind_i].get_chi()+global_ifermi;
    valtype v=c-medc[oind_i].get_Eg();
    fprintf(f,"%g\t",c);
    fprintf(f,"%g\t",v);
  }
  if(outtype&outConc){
    for(int car=0;car<2;car++)
      fprintf(f,"%g\t",conc[car][oind_i]*n0*1e-6);
  }
  if(outtype&outDop)fprintf(f,"%g\t",ch[oind_i]*n0*1e-6);

  valtype R=0,G=0;
  if(outtype&(outR_sc|outG_sc))
    recomb_gen(slice, shift, R, G);
  if(outtype&outR_sc)fprintf(f,"%g\t",R*1e-6);
  if(outtype&outG_sc)fprintf(f,"%g\t",G*1e-6);

  if(outtype&outJ){
    for(int car=0;car<2;car++){
      for(int di=0;di<3;di++){
        if(dim_type[di]<0)
          continue;
        fprintf(f,"%g\t",J[car][dnf[di]*SZ+oind_i]*n0*phys::q*1e-4);
      }
    }
  }
  if(outtype&outE_sc){
    const int *iind=slice.iind;
    for(int di=0;di<3;di++){
      if(dim_type[di]<0)
        continue;
      valtype E = slice.sz>1 ? J0*(psi_ptr[iind[dnum(di,1)]+shift]-psi_ptr[iind[dnum(di,0)]+shift])/(2*l0*dx_arr[di][icc[di]]) : 0;
      fprintf(f,"%g\t",E);
    }
  }
  if(outtype&outExcitons){
    for(int car=0;car<2;car++){
      // recalc from recombination
      valtype R=0,G=0;
      recomb_gen(slice, shift, R, G);
      CarrierParams exc_ini = medc[oind_i].get_excIniState(R); 
      fprintf(f,"%g\t",exc_ini[car]);

      //fprintf(f,"%g\t",excitons[car][oind_i]*n0*1e-6);
    }
  }
  if(outtype&outResidue){
    fprintf(f,"%g\t",y_ar[oind_i]);
    for(int car=0;car<2;car++)
      fprintf(f,"%g\t",y_ar[(1+car)*extsz+oind_i]);
  }
  fprintf(f,"\n");

  return 1;
}

RectBlock::node_t RectBlock::make_node(const slice_t &slice, int ind_shift){

  // if make_node is called from mesh_it, then oind = -1, iind = NULL
  // since at this moment memory indices are not assigned

  int oind=slice.oind;
  const int *iind=slice.iind;
  const int gind=slice.gind;

  node_t p;
  p.dim=dim;

  p.ind = oind<0 ? gind+ind_shift : oind+ind_shift;

  p.gind=gind+ind_shift;

  int ind[3];
  grd.unpack_ind(gind+ind_shift,ind[0],ind[1],ind[2]);
  p.center=get_position(ind[0],ind[1],ind[2]);
  p.v=grd.GetControlVolume(ind[0],ind[1],ind[2]);

  int sh[3]={0,0,0};
  for(int di=2;di>=0;di--){
    if(dim_type[di]<0)
      continue;
    for(int ni=0;ni<2;ni++){
      int sign=ni?1:-1;
      sh[di]=sign;

      int dn=dnum(di,ni); // pack direction
      if(iind)
        p.adj_ind[dn]=iind[dn]+ind_shift;
      else{ // make_node is called from mesh_it
        int iind=grd.pack_ind(ind[0]+sh[0],ind[1]+sh[1],ind[2]+sh[2]);
        p.adj_ind[dn]=iind;
      }

      if((ni==0 && ind[di]<=0) || (ni==1 && ind[di]>=N[di]-1)){
        p.adj_ind[dn]=-1; // this node does not exist
        continue;
      }

      p.adj[dn]=get_position(ind[0]+sh[0],ind[1]+sh[1],ind[2]+sh[2]);
      p.cur[dn]=(p.center+p.adj[dn])/2;
    }
    sh[di]=0;
  }

  return p;
}

template class NonUniformGrid<valtype,1,indtype>;
