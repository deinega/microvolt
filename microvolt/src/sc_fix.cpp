#include "sc_fix.h"
#include "sc_fd.h"
#include "sc_mesh.h"
#include "data_flags.h"

#include "physconst.h"

vector<string> out_names(int argtype, int outtype){

  vector<string> str;

  char x[3]={'x','y','z'}, cnc[2]={'n','p'}, exc[2]={'S','T'};

  if(argtype&argx)str.push_back("x");
  if(argtype&argy)str.push_back("y");
  if(argtype&argz)str.push_back("z");

  if(outtype&outPsi)str.push_back("psi[eV]");
  if(outtype&outFermi){
    for(int car=0;car<2;car++)
      str.push_back(string("fermi") + cnc[car]+"[eV]");
  }
  if(outtype&outBands){
    str.push_back("cond_band[eV]");
    str.push_back("val_band[eV]");
  }
  if(outtype&outConc){
    for(int car=0;car<2;car++)
      str.push_back(string(1, cnc[car]) + "[1/cm3]");
  }
  if(outtype&outDop)str.push_back("dop[1/cm3]");
  if(outtype&outR_sc)str.push_back("R[1/sec/cm3]");
  if(outtype&outG_sc)str.push_back("G[1/sec/cm3]");

  if(outtype&outJ){
    for(int car=0;car<2;car++){
      if(argtype){
        for(int i=0;i<3;i++){
          if(argtype&(argx<<i))
            str.push_back(string("J") + cnc[car] + x[i] + "[A/cm2]");
        }
      }
      else
        str.push_back(string("J") + cnc[car] + "[A/cm2]");
    }
  }
  if(outtype&outJtot){
    if(argtype){
      for(int i=0;i<3;i++){
        if(argtype&(argx<<i))
          str.push_back(string("J") + x[i] + "[A/cm2]");
      }
    }
    else
      str.push_back(string("J") + "[A/cm2]");
  }

  if(outtype&outE_sc){
    if(argtype){
      for(int i=0;i<3;i++){
        if(argtype&(argx<<i))
          str.push_back(string("E") + x[i] + "[eV/m]");
      }
    }
    else
      str.push_back(string("E[eV/m]"));
  }

  if(outtype&outExcitons){
    for(int car=0;car<2;car++)
      str.push_back(string(1, exc[car]));
  }
  if(outtype&outResidue){
    str.push_back(string("residue_psi"));
    for(int car=0;car<2;car++)
      str.push_back(string("residue_") + cnc[car]);
  }


  return str;
}

valtype SG_approx(valtype coef,valtype n0,valtype n1,valtype psi0,valtype psi1,int car){

  int sign = car ? -1: 1;

  valtype psidif=psi1-psi0;
  psidif*=sign;
  valtype e_ax=exp(psidif*coef);
  valtype e_ah=exp(psidif);

//  valtype psi_interp = psi0*(1-coef)+psi1*coef;
  valtype n = psidif ? 
    (n1*(e_ax-1)/(e_ah-1) + n0*(e_ah-e_ax)/(e_ah-1)) : 
    (n1*coef              + n0*(1-coef));
  return n;
}

template<class cont_t>
int SG_interp(const scInterpolation &form,int car,const local_t &med,
  const cont_t *cont,valtype *psi_val,valtype *fermi_val,valtype *conc_val){

  valtype coef[4]; // interpolation coefficients
  valtype psi[4],fermi[4],c[4]; // values at the square vertices

  for(int i=0;i<4;i++){
    if(i>=form.size()){
      psi[i]=fermi[i]=c[i]=coef[i]=0;
      continue;
    }

    psi[i]=cont->psi_ptr[form.shifts[i]];
    fermi[i]=cont->fermi[car][form.shifts[i]];
    c[i]=cont->conc[car][form.shifts[i]];
    coef[i]=form.coeffs[i];

    if(form.size()==1){
      if(conc_val)
        *conc_val=c[i];
      if(psi_val)
        *psi_val=psi[i];
      if(fermi_val)
        *fermi_val=fermi[i];
      return 1;
    }
  }

  valtype psi_line[2],c_line[2],sum[2];
  for(int i=0;i<2;i++){
    sum[i]=coef[2*i]+coef[2*i+1];
    if(sum[i]==0){
      psi_line[i]=psi_line[(i+i)%2]; // is it correct?
      c_line[i]=c_line[(i+i)%2]; // is it correct?
    }
    else{
      psi_line[i]=(psi[2*i]*coef[2*i]+psi[2*i+1]*coef[2*i+1])/sum[i];
      c_line[i] = SG_approx(coef[2*i+1]/sum[i],c[2*i],c[2*i+1],psi[2*i],psi[2*i+1],car);
    }
  }

//  valtype conc_interp = SG_approx(form.coef[1],c[0],c[1],psi[0],psi[1],sign);
//  valtype psi_interp = psi[0]*form.coef[0]+psi[1]*form.coef[1];

  valtype conc_interp = SG_approx(sum[1],c_line[0],c_line[1],psi_line[0],psi_line[1],car);
  if(conc_val)
    *conc_val=conc_interp;

  valtype psi_interp = psi_line[0]*sum[0]+psi_line[1]*sum[1];
  if(psi_val)
    *psi_val=psi_interp;

  valtype coefexp=cont->coefexp;
  valtype J0=cont->J0;
  valtype n0=cont->n0;

  valtype ni=med.get_ni();
  valtype ifermi=med.get_ifermi()/J0;

  psi_interp=psi_interp-ifermi;
  valtype cnc=std::log(conc_interp/(ni/n0))/coefexp;
  valtype fermi_interp = car ? psi_interp+cnc : psi_interp-cnc;
  if(fermi_val)
    *fermi_val=fermi_interp;

  return 1;
}

template<class container_t>
scFixInterp<container_t>::scFixInterp(container_t *cont_):cont(cont_),//delta(1e-9),
argtype(argx|argy|argz),outtype(0xffff&(~outJtot)){
  const int *dim_type=cont->dim_type;
  if(dim_type[0]<0)argtype&=~argx;
  if(dim_type[1]<0)argtype&=~argy;
  if(dim_type[2]<0)argtype&=~argz;
}

template<class container_t>
void scFixInterp<container_t>::AddDefaultTimers(const vector<int> &p_id){
  if(p_id.size()!=2)return;
  //int build=p_id[0];
  int calc=p_id[1];
  tid.tr_lcl=AddTimersField("tr_lcl",calc);
}

template<class container_t>
size_t scFixInterp<container_t>::get_full_data_size()const{
  int dim=cont->get_dim();
  int fs=0;
  if(outtype&outPsi)fs+=1;
  if(outtype&outFermi)fs+=2;
  if(outtype&outBands)fs+=2;
  if(outtype&outConc)fs+=2;
  if(outtype&outDop)fs+=1;
  if(outtype&outR_sc)fs+=1;
  if(outtype&outG_sc)fs+=1;
  if(outtype&outJ)fs+=2*dim;
  if(outtype&outJtot)fs+=dim;
  if(outtype&outExcitons)fs+=2;
  fs*=sizeof(valtype);
  return fs;
}

template<class container_t>
int scFixInterp<container_t>::put_full_request(int ibuf, int &ind, const Vector_3 &pos, int rec_out){

  scInterpolation interp=cont->create_interpolation(pos,0,0,1); // incomplete interpolation is allowed
  if(!interp.valid() && !rec_out)
    return -1;

  int stind=ind*get_full_data_size()/sizeof(valtype);
  ibufs.push_back(ibuf);
  stinds.push_back(stind);
  ind++;

  forms.next_record(interp);
  if(outtype&(outJ|outJtot)){
    for(int di=0;di<3;di++){
      if(cont->dim_type[di]<0)
        continue;
      Vector_3 dir;
      dir[di]=1;
      scInterpolation interpj=cont->create_interpolation(make_pair(pos,dir),0,0,1); // incomplete interpolation is allowed
      j_forms.next_record(interpj);
    }
  }

  meds.append(*(cont->media_regions), pos);

  valtype gen=cont->test_gen(pos);
  gens.push_back(gen);

/*  if(outtype&outE_sc){
    for(int di=2;di>=0;di--){
      if(cont->dim_type[di]<0)
        continue;
      for(int ni=1;ni>=0;ni--){
        Vector_3 sh;
        sh[di] = (ni ? 1: -1)*delta;
        scInterpolation interp=cont->create_interpolation(pos+sh);
        idelta.push_back(make_pair(interp,ni));
      }
    }
  }*/

  return 1;
}

template<class container_t>
int scFixInterp<container_t>::compute_local(){

  start(tid.tr_lcl);

//  valtype coefexp=cont->coefexp;
//  valtype J0=cont->J0;
  valtype n0=cont->n0;

//  int nb=(int)buff.size();
  int dim=cont->get_dim();
  const int *dim_type=cont->dim_type;
  
  for(size_t i=0;i<forms.size();i++){

    int ibuf=ibufs[i]; // buffer number
    int ind=stinds[i]; // indez in buffer

    const scInterpolation form=forms[i];
    val_t psi=cont->get_value(form,rpPsi);
    
    val_t fermi[2],conc[2], excitons[2];
    for(int car=0;car<2;car++){
      fermi[car]=cont->get_value(form,rpFermi,car);
      conc[car]=cont->get_value(form,rpConc,car);
      excitons[car]=cont->get_value(form,rpExcitons,car);
    }

    if(outtype&outPsi){
      buff[ibuf][ind++]=psi*cont->J0-cont->global_ifermi;
    }

    if(outtype&outFermi){
      for(int car=0;car<2;car++){
        buff[ibuf][ind++]=-fermi[car]*cont->J0+cont->global_ifermi;
      }
    }

    if(outtype&outBands){
      valtype c=-psi*cont->J0-meds[i].get_chi()+cont->global_ifermi;
      buff[ibuf][ind++]=c;
      buff[ibuf][ind++]=c-meds[i].get_Eg();
    }

    if(outtype&outConc){
      for(int car=0;car<2;car++)
        buff[ibuf][ind++]=conc[car]*n0*1e-6;
    }

    if(outtype&outDop)
      buff[ibuf][ind++] = meds[i].get_dop()*1e-6;

    valtype R = meds[i].get_recomb(cont->n0*conc[0], cont->n0*conc[1]);
    valtype G = gens[i];

    if(outtype&outR_sc)
      buff[ibuf][ind++]=R*1e-6;
    if(outtype&outG_sc)
      buff[ibuf][ind++]=G*1e-6;

    valtype J[2][3];

    for(int di=0,dic=0;di<3;di++){

      if(dim_type[di]<0)
        continue;

/*      valtype psiexp=exp(psidif*coefexp);
      valtype B[2];
      B[0]=accomp(psiexp,1.) ? 1 :psidif/(psiexp-1);
      B[1]=B[0]*psiexp;*/

      if(outtype&(outJ|outJtot)){
        for(int car=0;car<2;car++){
/*        valtype D=med.D[car];
        int sign=-nsign*(car?1:-1);

        // work unpredictibly
        valtype iconc=cont->get_value(*iform[0],rpConc,car);
        valtype val=B[car]*iconc-B[1-car]*conc[car]; // for electrons psi_der-psi is with c_der
        val*=D/delta;
        J[car][di]=sign*val;

        // work inaccurate at the borders and pn-junction
        valtype ifermi[2];
        for(int ni=0;ni<2;ni++)
          ifermi[(ni+1)%2]=cont->get_value(*iform[ni],rpFermi,car);
        valtype mu=med.mu[car];
        J[car][di]=-mu*conc[car]*J0*(ifermi[1]-ifermi[0])/(delta*dvd);
        J[car][di]=sign*D*exp(sign*coefexp*psi)*(exp(-sign*coefexp*ifermi[1])-exp(-sign*coefexp*ifermi[0]))/(delta*dvd);*/

          J[car][di]=cont->get_value(j_forms[dim*i+dic],rpJ,car)*n0*phys::q*pow(units::cm, 2);
        }
      }
      dic++;
    }

    if(outtype&outJ){
      for(int car=0;car<2;car++){
        for(int di=0;di<3;di++){
          if(dim_type[di]<0)
            continue;
          buff[ibuf][ind++]=J[car][di];
        }
      }
    }
    if(outtype&outJtot){
      for(int di=0;di<3;di++){
        if(dim_type[di]<0)
          continue;
        buff[ibuf][ind++]=(J[0][di]+J[1][di]);
      }
    }

    if(outtype&outExcitons){
      for(int car=0;car<2;car++)
        buff[ibuf][ind++]=excitons[car]*n0*1e-6;
    }

  }

  stop(tid.tr_lcl);

  return 1;
}

template<class cont_t>
void FixTransfer<cont_t>::fill_back_item(map<int,int> &imap,int shift,int in,int ogn){
  map<int,int>::iterator mip=imap.find(shift);
  int vind;
  if(mip!=imap.end()){
    vind=mip->second;
  }
  else{
    x_ind.push_back(shift);
    backs.push_back(ind_t());
    vind=x_ind.size()-1;
    imap[shift]=vind;
  }
  int fn=backs[vind].n;
  backs[vind].igroup[fn]=in; // index in transfer oind
  backs[vind].ogroups[fn]=ogn; // index in transfer nind
  backs[vind].n++;
  if(backs[vind].n>=max_interp)
    throw 1;
}

template<class cont_t>
int FixTransfer<cont_t>::init(){
  // key - memory index of X node
  // value - index in back packer
  map<int,int> imap;
  int pn=0;
  
  for(size_t oi=0;oi<forms.size();oi++){
    const scInterpolation form=forms[oi];
    int fon=form.size();
    for(int i=0;i<fon;i++)
      fill_back_item(imap,form.shifts[i],oi,pn);
    pn+=o_grsz[oi];
  }

  return 1;
}

template<class cont_t>
void FixTransfer<cont_t>::set_value(int vi,int rep,int srp,int remember){

  for(int mask=1;mask<=rpJ;mask<<=1){
    if(!(rep&mask))
      continue;
    for(int j=0;j<2;j++){
      valtype *ptr=NULL;
      if(mask&rpP){
        if(j==0)
          ptr=get_pointer(*cont,mask);
        else
          break;
      }
      else{
        int maskEH=rpE<<j;
        if(!(rep&maskEH))
          continue;
        ptr=get_pointer(*cont,mask|maskEH);
      }

      if(remember)
        toval.push_back(ptr[i_ind[vi]]); // remembering previous value

      const scInterpolation form=forms[vi];
      valtype sum=0;

      for(int i=0;i<form.size();i++){ // calculation interpolation
        sum+=ptr[form.shifts[i]]*form.coeffs[i];
      }
      ptr[i_ind[vi]]=sum;
    }
  }

  set_repres(*cont,srp,i_ind[vi]);
}

template<class cont_t>
void FixTransfer<cont_t>::set_previous_value(int vi,int rep,int srp){

  for(int mask=1;mask<=rpJ;mask<<=1){
    if(!(rep&mask))
      continue;
    for(int j=0;j<2;j++){
      valtype *ptr=NULL;
      if(mask&rpP){
        if(j==0)
          ptr=get_pointer(*cont,mask);
        else
          break;
      }
      else{
        int maskEH=rpE<<j;
        if(!(rep&maskEH))
          continue;
        ptr=get_pointer(*cont,mask|maskEH);
      }

      ptr[i_ind[vi]]=toval[tind++];
    }
  }

  set_repres(*cont,srp,i_ind[vi]);
}

template<class cont_t>
int FixTransfer<cont_t>::clear(){

  i_mesh.clear();
  i_ind.clear();

  forms.clear();

  o_gr_ind.clear();
  o_grsz.clear();

  x_ind.clear();
  backs.clear();

  toval.clear();

  return 1;
}

template<class cont_t>
int FixPsi<cont_t>::record(int oindi,int obi,valtype *psi,valtype psi_built, valtype wf_delta){

  scFix<cont_t>::record(oindi,obi);

  in.push_back(cont->blocks[obi]->offset+oindi);
  in_sz.push_back(1);

  psis.push_back(psi);
  psi_builts.push_back(psi_built);
  wf_deltas.push_back(wf_delta);
  return 1;
}

template<class cont_t>
valtype FixPsi<cont_t>::y(int vi,int repr){

  StaticInterpolation<valtype,100> st_form;
  st_form.shifts[0]=out[vi];
  st_form.coeffs[0]=1;
  st_form.n=1;

  scInterpolation form(0,1,st_form.shifts,st_form.coeffs,2);


  valtype psi0=(*psis[vi] + psi_builts[vi] + wf_deltas[vi]);
  psi0/= cont->J0;
  valtype psic=cont->get_value(form,rpPsi,0);
  return psic-psi0;
}

template<class cont_t>
int FixSR<cont_t>::record(int oindi,int obi,int bj, int eq,const Vector_3 &c, const Vector_3 &near, const Vector_3 &nvect, 
  const Vector_3 &surfp, valtype wf_delta,  const local_t &med,  valtype *sr){

  scFix<cont_t>::record(oindi,obi);
  ob.push_back(obi);

  if(cont->blocks[obi]->block->fill_adjacent_border(oindi,in,in_sz,1|2|8,mInt|mFix)<0)
    return -1;

  scInterpolation conc=cont->create_interpolation(c,obi,1);
  concs.next_record(conc);

  scInterpolation J=cont->create_interpolation(make_pair(c,nvect),obi,1);
  Js.next_record(J);

  if(!conc.valid() || !J.valid())
    return ::pmessage(vblMESS1,-1,"BlockContainer is not covered by meshes properly (metal contact initialization)\n");

  bjs.push_back(bj);

  concs0.push_back(med.get_conc0()/cont->n0);

  srs.push_back(sr);

  wf_deltas.push_back(wf_delta);

  bound_type.push_back(eq);

//  scInterpolation nform=cont->create_interpolation(near,obi);
  scInterpolation nform=cont->create_interpolation(near,obi,0,1,1);
  if(!nform.valid()){
    return ::pmessage(vblMESS1,-1,"Too sharp interface or two meshes border close to interface (%g, %g, %g)\n",near[0],near[1],near[2]);
  }
  nforms.next_record(nform);

  deltas.push_back((c-near).norm()*cont->lu);

  meds.append_ref_to(&med);

  return 1;
}

template<class cont_t>
valtype FixSR<cont_t>::y(int vi,int repr){
  using namespace phys;

  if(repr&rpPsi)
    return 0;

  int is = get_is(repr);
  int sign = is ? 1 : -1;
//  int conc_repr = rpConc | (repr&rpEH);
  int cur_repr = rpJ | (repr&rpEH);

  valtype res = 0;

  cont->blocks[ob[vi]]->block->update_bJ(bjs[vi]);
  valtype conc[2];
  for(int i = 0; i<2; i++)
    conc[i] = cont->get_value(concs[vi], rpConc, i);
  valtype cur = conc[is];

  // interpolation current
  valtype J_interp = sign*cont->get_value(Js[vi],cur_repr);

  // calculating current using SG_interp
  valtype delta = deltas[vi];

  valtype coefexp = cont->coefexp;
  valtype J0 = cont->J0;
  valtype n0 = cont->n0;

  valtype ni = meds[vi].get_ni();
  valtype ifermi = meds[vi].get_ifermi()/J0;

  valtype D = meds[vi].get_D()[is];
//  valtype mu = meds[vi].get_mu()[is];

  scInterpolation form_der = nforms[vi];
  valtype fermi = cont->get_value(concs[vi], rpFermi, is);
  valtype fermi_der = cont->get_value(form_der, rpFermi, is); // linear interpolation - leads to error
  SG_interp(form_der, is, meds[vi], cont, NULL, &fermi_der, NULL); // using SG_interp
  valtype psi = cont->get_value(concs[vi], rpPsi) - ifermi;
  valtype psi_der = cont->get_value(form_der, rpPsi) - ifermi;

  valtype psidif = psi_der - psi;
  valtype psiexp = exp(psidif*coefexp);
  valtype B[2];
  B[0] = accomp(psiexp, 1.) ? 1 : psidif/(psiexp - 1);
  B[1] = B[0]*psiexp;

  //valtype c = cont->get_value(concs[vi], rpConc, is);
  //valtype c_der = cont->get_value(form_der, rpConc, is);

  valtype c = ni/n0*exp(sign*(fermi - psi)*coefexp);
  valtype c_der = ni/n0*exp(sign*(fermi_der - psi_der)*coefexp);

  valtype val = B[is]*c_der - B[1-is]*c;
  val *= D/delta;
  valtype J_SG = -val;

//  valtype Jtest_cf=-sign*mu*c*J0*(fermi-fermi_der)/delta;
//  valtype Jtest_fp=-D*exp(-sign*coefexp*psi)*(exp(sign*coefexp*fermi_der)-exp(sign*coefexp*fermi))/delta;

  calc_cur=0;
  valtype J = calc_cur ? J_SG : J_interp;

  // depletion factor before concentrations for Schottky contact (wf_delta!=0)
  // (-) for electrons, (+) for holes 
  valtype deplet = exp(-sign*wf_deltas[vi]/(kB_eV*cont->T));

  valtype *sr = srs[vi];

  CarrierParams conc0 = concs0[vi];

  // The condition below is inadequate: it should choose between organic
  // and inorganic media, though it works in the current applications.
  // TODO: revise the code; branching via polymorphism may be a better choice
  if (!meds[vi].getDefaultContext().empty()) {
    CarrierParams barrier(wf_deltas[vi] + meds[vi].get_workfunc());

    // Adjust the barrier due to the reflected charge effect
    if (meds[vi].get_useReflectedCharge()) {
      int psiSign = psidif < 0 ? 1 : -1;
      valtype refChargeCorrection =
        sqrt(q/(4*M_PI*eps0)*abs(psidif)*(kB_eV*cont->T)/delta/meds[vi].get_eps());
      barrier[psidif < 0] -= psiSign*refChargeCorrection;
    }

    conc0 = meds[vi].calcConcentration(-barrier/*,vi? 0.: 0.05*/);
    //CarrierParams conc1 = meds[vi].calcConcentration(-barrier,1e-8)/n0;
    
    // Add surface recombination that makes n*p==ni^2 on the boundary
    /*
    valtype concDiff = conc0[0] - conc0[1];
    valtype np0 = meds[vi].get_np();
    valtype r = 2*(conc0[0]*conc0[1] - np0) /
      (conc0[0] + conc0[1] + sqrt(concDiff*concDiff + np0));
    int nMajor = conc0[0] < conc0[1];
    conc0[nMajor] -= r;
    conc0[1-nMajor] = np0 / conc0[nMajor]; */
    

    conc0 /= n0;
    //conc1-=conc0;
    /*
    valtype psi = cont->get_value(concs[vi], rpPsi)*J0;
    CarrierParams f = meds[vi].calcQFermi(conc0, psi);
    valtype fermi = cont->get_value(concs[vi], rpFermi, is)*J0;
    res = (fermi+f[is])*conc[is];
    res = fermi;
    if(is==1)
      printf("%g\n",res);
    */
    deplet = 1;
    

    // TODO: consider surface recombination according to Malliaras paper
    //sr[is] = 16*M_PI*meds[vi].get_eps()*phys::eps0*D*phys::kB_eV*cont->T/phys::q;
  }

  if(sr[is] == VEC_INFTY){ // inifinite surface recombination, equilibrium concentration at the surface
    res += cur - conc0[is]*deplet;
  }
  else if(sr[is] > 0){
    if(bound_type[vi] == 1){ // metal contact
      res += cur - conc0[is]*deplet + J/sr[is];
    }
    else{ // dielectric interface
      valtype num = conc[0]*conc[1] - conc0[0]*conc0[1];
//      valtype denum=(conc[0]+conc0[0])/sr[1]+(conc[1]+conc0[1])/sr[0];
      valtype denum = (conc0[0])/sr[1] + (conc0[1])/sr[0];
//      valtype denum=(conc[0])/sr[1]+(conc[1])/sr[0];
      res += num*deplet/denum + J;
/*      int min = conc0[1]<conc0[0];
      res+=conc[min]-conc0[min]+J/sr[min];
//      res+=sr[min]*(conc[min]-conc0[min])+J;*/
    }
  }
  else
    res += J;

  return res;
}

template<class cont_t>
int FixNeumann<cont_t>::record(int oindi,int obi,scInterpolation &form, const local_t &med){

  scFix<cont_t>::record(oindi,obi);
  for(int i=0;i<form.size();i++){
    in.push_back(form.shifts[i]); // influencing indices
  }
  in.push_back(cont->blocks[obi]->offset+oindi); // index influences itself (CHECK IF IT IS ALREADY RECORDED)
  in_sz.push_back(form.size()+1);
  forms.next_record(form);
  meds.append_ref_to(&med);
  return 1;
}

template<class cont_t>
valtype FixNeumann<cont_t>::y(int vi,int repr){

  const scInterpolation form=forms[vi];

  valtype val=cont->get_value(form,repr);
  valtype *ptr=get_pointer(*cont,repr);
  valtype cval=ptr[out[vi]];

  if(form.size()>1 && repr&rpFermi){ // Scharfetter-Gummel approximation
    int car=get_is(repr);
    SG_interp(form,car,meds[vi],cont,NULL,&val,NULL);
  }
  valtype res=val-cval;

  return res;
}

template<class cont_t>
int FixNeumannDer<cont_t>::record(int oindi,int obi,const Vector_3 &center,const Vector_3 &norm){

  scFix<cont_t>::record(oindi,obi);

  if(cont->blocks[obi]->block->fill_adjacent_border(oindi,in,in_sz,1|2|8,mInt|mFix)<0)
    return -1;

  scInterpolation form=cont->create_interpolation(center);
  cforms.next_record(form);
  for(int i=0;i<3;i++){
    Vector_3 pos=center;
    pos[i]+=dr*norm[i];
    form=cont->create_interpolation(pos);
    nforms[i].next_record(form);
  }
  norms.push_back(norm);

  return 1;
}

template<class cont_t>
valtype FixNeumannDer<cont_t>::y(int vi,int repr){

  if(!(repr&wrepr))
    return 0;

  valtype val=cont->get_value(cforms[vi],repr);
  valtype dval[3];
  for(int i=0;i<3;i++)
    dval[i]=cont->get_value(nforms[i][vi],repr);

  valtype res=0;
  Vector_3 norm=norms[vi];
  for(int i=0;i<3;i++)
    res+=(dval[i]-val)*fabs(norm[i]);

  return res/dr;
}

template<class cont_t>
int FixBlock<cont_t>::record(int oindi,int obi){
  scFix<cont_t>::record(oindi,obi);
  oind.push_back(oindi);
  ob.push_back(obi);

  !axis ?
    cont->blocks[obi]->block->fill_adjacent_internal(oindi,in,in_sz,1|2,mInt|mFix) : 
    cont->blocks[obi]->block->fill_adjacent_border(oindi,in,in_sz,1|2,mInt|mFix);

  int gind;
  // put every boundary point (including absent in memory array) and not sort
  // in order to properly use block_t::bulk_fun
  // put -1 for points not in memory array
  !axis ?
    cont->blocks[obi]->block->fill_adjacent_internal(oindi,bind,bindsz,4|0x10,mAll,&gind) :
    cont->blocks[obi]->block->fill_adjacent_border(oindi,bind,bindsz,4|0x10,mAll,&gind);

  ginds.push_back(gind);
  return 1;
}

template<class cont_t>
valtype FixBlock<cont_t>::y(int vi,int repr){

  typename cont_t::block_t::slice_t slice;
  slice.oind=oind[vi];
  int bsz=2*cont->dim*vi;
  for(int i=0;i<2*cont->dim;i++)
    slice.iind[i]=bind[bsz+i];
  slice.gind=ginds[vi];
  cont->blocks[ob[vi]]->block->unpack_ind(ginds[vi],slice.ic[0],slice.ic[1],slice.ic[2]);
  valtype val;

  if(repr&rpPsi)
    cont->blocks[ob[vi]]->block->Poisson(slice,0,&val);
  else
    cont->blocks[ob[vi]]->block->Continuity(slice,0,&val,get_bf(repr),0);

  return val;
}

//template int SG_interp(const scInterpolation &form,int car,const local_t &med, const scBlockContainer<> *cont,valtype *psi,valtype *fermi,valtype *n);

template class scFixInterp<scBlockContainer<> >;
template struct FixTransfer<scBlockContainer<> >;
template class FixPsi<scBlockContainer<> >;
template class FixSR<scBlockContainer<> >;
template class FixNeumannDer<scBlockContainer<> >;
template class FixNeumann<scBlockContainer<> >;
template class FixBlock<scBlockContainer<> >;
