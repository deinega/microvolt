#ifndef SC_H
#define SC_H

/// @file sc.h  Some auxiliary classes used to manage memory layout, Jacobian construction, etc
/// Documentation to the general structure of the code can be found at doc/tutorial.doxc

#include "linsolv.h"
#include "refobj.h"
#include "component.h"
#include "sc_carrier_params.h"

typedef int indtype;  // this is for 32-bit architecture
//typedef ptrdiff_t indtype; // universal
//typedef long long indtype; // test

// linear solver for Jacobian
#ifdef USE_PARDISO
typedef pardiso_solver<valtype> lin_solver;
#else
typedef linear_solver<valtype> lin_solver;
#endif

// return bit flag (first bit correspond to electrons, second bit correspond to holes) 
// from representation flag which is assumed to have rpE or rpH bit (some of them or both)
inline int get_bf(int repr){
  int bf=0;
  if(repr&rpE) bf|=1;
  if(repr&rpH) bf|=2;
  return bf;
}

// return 0 (for electrons) or 1 (for holes) from representation flag 
// which is assumed to have rpE or rpH bit (not both)
inline int get_is(int repr){
  if(repr&rpH) return 1;
  else return 0;
}

enum mFlags{ // mesh nodes classification, see documentation to scBlockContainer::specify_nodes
  mInt=0x1, // internal node
  mFixBound=0x2, // boundary fix node
  mFixBoundSurround=0x4, // boundary fix node which surround other boundary fix nodes (not used)
  mFixInt=0x8, // internal fix node

  mFix=mFixBound|mFixBoundSurround|mFixInt,

  mTrans=0x10, // node covered by other higher level mesh

  mAll=0xfff
};

// boundary (metal contact, insulator interface)
struct scBoundaryRegion: public apComponent{
  mngptr<Region_3> reg; // boundary position
  valtype psi; // applied voltage (used for metal contact)
  valtype sr[2]; // surface recombination velocities for electrons and holes, VEC_INFTY corresponds to infinity
  valtype workf; // metal workfunction for Schottky contact, VEC_INFTY indicates Ohmic contact (used for metal contact)
  // type of boundary
  // 1 - metal contact
  // 2 - insulator interface
  // 3 - cylindrical axis for simulation in radial coordinates
  // 4 - periodic boundary conditions
  int type;
  int appl_volt; // if voltage is applied (only for metal contacts)

  scBoundaryRegion(mngarg<Region_3> sreg, int type_, valtype sn, valtype sp, valtype workf_=VEC_INFTY, int appl_volt_=0): 
  type(type_), reg(sreg), psi(0), workf(workf_), appl_volt(appl_volt_){
    sr[0]=sn;
    sr[1]=sp;
  }
};

// convert binary flag for argument (argTYPE) and output (scoutTYPE) types into vector of names
std::vector<std::string> out_names(int argtype, int outtype);

/*
Size of Jacobian is SZ X SZ, 2*SZ X 2*SZ or 3*SZ X 3*SZ
where SZ is number of mesh nodes (of all meshes) for scalar variables

We can make different Jacobians delta_y_i / delta_x_j,
where index i corresponds to Poisson equation, continuity equation for ehectons or continuity equation for holes,
and index j corresponds to electrostatic potential, 
quasi-fermi potential for electrons or quasi-fermi potential for holes

We can denote these indices by following char symbols:
 p - Poisson equation
 e - continuity equation for electrons
 h - continuity equation for holes

 P - electrostatic potential
 E - quasi-fermi level for electrons
 H - quasi-fermi level for holes

These are examples of possible Jacobians:
 pPeEhH - Jacobian for all eqiations and variables
 pP - Jacobian for Poisson equation and electrostatic potential (as an argument of equation)
 pE - Jacobian for Poisson equation and quasi-fermi potential (as an argument of equation)

In the following if equation and variable are denoted by the same symbol, we skip capital letter, for example,
first and second examples of Jacobian could be labeled as peh and p.
*/

// this class determines Jacobian structure
struct method_t{
  string code; // Jacobian string label
  int nr; // number of equations (variables) in Jacobian
  int yrepr[3],xrepr[3]; // representation flag for each equation and variable
  int snum_gs,it_gs; // gauss seidel method parameters
  int yreprs, xreprs; // representation flag for all equations and variables in Jacobian

  // return representation flag from its char code
  // if char is not capital, sm will be true (otherwise false)
  int chr(char let,bool &sm){
    sm = let=='p' || let=='e' || let=='h';
    if(let=='p' || let=='P')
      return rpPsi;
    if(let=='e' || let=='E')
      return rpFermi|rpE;
    if(let=='h' || let=='H')
      return rpFermi|rpH;
    return 0;
  }

  // specify Jacobian structure from its char representation
  method_t(const string &code_="",int snum_gs_=1,int it_gs_=1): 
  code(code_),snum_gs(snum_gs_),it_gs(it_gs_),yreprs(0),xreprs(0){

    for(int i=0;i<3;i++){
      yrepr[i]=xrepr[i]=0;
    }

    nr=0;
    // wasi is 0 if we should fill xrepr for previous index
    for(size_t i=0,wasi=1;i<code.length();i++){
      bool if_y;
      int repr=chr(code[i],if_y);
      if(if_y){
        if(!wasi)
          xrepr[nr-1]=yrepr[nr-1];
        yrepr[nr]=repr;
        if(i==code.length()-1)
          xrepr[nr]=yrepr[nr];
        nr++;
        wasi=0;
      }
      else{
        xrepr[nr-1]=repr;
        wasi=1;
      }
    }

    for(int i=0;i<nr;i++){
      yreprs|=yrepr[i];
      xreprs|=xrepr[i];
    }
  }
};

/*
To find a solution we apply Newton iterations for sequences of different Jacobians.
For example, we can solve all equations by peh, or apply 10 iterations for eh
and 1 iteration for p (in literature solving continuity and Poisson equations sperately is called as Gummel method)

Structure meth_chain_t is used to specify sequence (chain) of Jacobians
These Jacobians are consequently solved in scBlockContainer::Step
*/
struct meth_chain_t{
  method_t meth[10]; // Jacobians
  int num[10]; // number of Newton's iterations for each Jacobian
  int n; // number of Jacobians in chain (maximal is 10)
  // chain of specified Jacobians will be used until solution converged to some (local) minimum, 
  // but not more than tot_num times (if tot_num>0)
  int tot_num;

  meth_chain_t():tot_num(50),n(0){}

  void add(method_t meth_,int num_=1){
    meth[n]=meth_;
    num[n]=num_;
    n++;
  }

  void set_tot_num(int tot_num_){
    tot_num=tot_num_;
  }
};

/*
In this class we specify a solver and arrays of variables 
for specific Jacobian defined in method_t

Some functions can process several variables at the same time, for example,
function Continuity can process fermi potential for electrons and for holes.

These types of variables can be grouped togerher in blocks
*/
struct meth_cont_t{

  method_t meth; // Jacobian structure

  int snum; // number of blocks

  // blocks:
  int gyrepr[3]; // representation flag for function
  valtype *x[3]; // array for variables
  int gxrepr[3]; // representation flag for these variables
  // array of flags which are turned on if corresponding variable is fixed,
  // see documentation to scBlockContainer::fixed
  int *fx[3];
  // number of variables for this block (more than one if variables are grouped to be processed by the same function)
  int fstep[3];
  // size of the array for each variable in the group 
  // currently it is always equal to SZ, but we keep variable for this for possible extension of the program
  int sz[3];
  int rp[3]; // representation that should be updated after variable is changed
  int srp; // sum of all representation flags rp

  lin_solver *ls; // solver for Jacobian

  // if gr_cont, then group variables that can be processed by the same function
  template<class cont_t>
  meth_cont_t(cont_t &cont, const method_t &meth_, int gr_cont=1): meth(meth_),snum(0){

    for(int i=0;i<3;i++){
      x[i]=NULL;
      gyrepr[i]=gxrepr[i]=0;
      fx[i]=NULL;
      fstep[i]=sz[i]=rp[i]=0;
    }

    for(int i=0,cf=0;i<meth.nr;i++){

      if(gr_cont && i<meth.nr-1){
        if(meth.yrepr[i]==meth.xrepr[i] && meth.yrepr[i+1]==meth.xrepr[i+1] &&
        (meth.yrepr[i]|meth.yrepr[i+1]) == (rpFermi|rpE|rpH)){ // grouping eh to one block
          fstep[snum]+=2;
          gyrepr[snum]=gxrepr[snum]=meth.yrepr[i]|meth.yrepr[i+1];
          snum++;
          i++;
          continue;
        }
      }
      fstep[snum]++;
      gyrepr[snum]=meth.yrepr[i];
      gxrepr[snum]=meth.xrepr[i];
      snum++;

/*      if(cf){
        if(meth.yrepr[i]==meth.xrepr[i] && meth.yrepr[i]&rpH){ // grouping eh to one block
          fstep[snum]++;
          gyrepr[snum]|=meth.yrepr[i];
          gxrepr[snum]|=meth.xrepr[i];
          snum++;
        }
        else{
          snum++;
          fstep[snum]++;
          gyrepr[snum]=meth.yrepr[i];
          gxrepr[snum]=meth.xrepr[i];
          snum++;
        }
        cf=0;
      }
      else if(gr_cont && meth.yrepr[i]==meth.xrepr[i] && meth.xrepr[i]&rpE){ // at the next step we check if this is rpH
        cf=1;
        fstep[snum]++;
        gyrepr[snum]=meth.yrepr[i];
        gxrepr[snum]=meth.xrepr[i];
        if(i==meth.nr-1)
          snum++;
      }
      else{
        fstep[snum]++;
        gyrepr[snum]=meth.yrepr[i];
        gxrepr[snum]=meth.xrepr[i];
        snum++;
      }*/
    }

    for(int i=0;i<snum;i++){

      if(gxrepr[i]&rpP){
        x[i]=cont.psi_ptr;
        fx[i]=cont.fixed;
      }
      if(gxrepr[i]&rpFermi){
        x[i]=cont.fermi[0];
        fx[i]=cont.fixed+cont.SZ;
        if(!(gxrepr[i]&rpE)){
          x[i]+=cont.SZ;
          fx[i]+=cont.SZ;
        }
      }

      sz[i]=cont.SZ;
      rp[i]=rpConc;
    }

    srp=0;
    for(int j0=0;j0<snum;j0++)
      srp|=rp[j0];

    if(meth.nr==1){
      if(meth.yrepr[0]&rpP)
        ls=&cont.ls_p;
      else
        ls=&cont.ls_c;
    }
    else if(meth.nr==2){
      if(meth.yreprs==(rpFermi|rpE|rpH)){
        ls=&cont.ls_c2;
      }
      else{
        ls=&cont.ls_pc;
      }
    }
    else
      ls=&cont.ls_pc2;
  }
};

// return pointer at array of variables corresponding to representation flag
template<class cont_t>
valtype *get_pointer(cont_t &cont, int rep){
  if(rep&rpPsi)return cont.psi_ptr;
  int ind=(rep&rpH)!=0;
  if(rep&rpFermi)return cont.fermi[ind];
  if(rep&rpConc)return cont.conc[ind];
  if(rep&rpJ)return cont.J[ind];
  ind=(rep&rpT)!=0;
  if(rep&rpExcitons)return cont.excitons[ind];
  return NULL;
}

// update value of scalar variable which is expressed from electrostatic and quasi-fermi potentials
template<class cont_t>
void set_repres(cont_t &cont, int rp, int i){

  valtype ifermi=cont.medc[i].get_ifermi()/cont.J0;
  valtype ni=cont.medc[i].get_ni()/cont.n0;
  valtype psi=cont.psi_ptr[i]-ifermi;

  if(rp&rpConc){
    if (!cont.medc[i].getDefaultContext().empty()) {
      valtype coef = cont.J0*cont.coefexp;
      CarrierParams frm(cont.fermi[0][i], cont.fermi[1][i]);
      CarrierParams conc = cont.medc[i].calcConcentration(coef*frm, coef*cont.psi_ptr[i])/cont.n0;

      for(CarrierParams::size_type car = 0; car < conc.size(); car++)
        cont.conc[car][i] = conc[car];
    }
    else { // this somehow influences convergence process
      CarrierParams conc_0;;
      for(int car=0;car<2;car++){
        valtype frm=cont.fermi[car][i];
        valtype expdif = car ? frm-psi : psi-frm;
        conc_0[car]=ni*exp(expdif*cont.coefexp);  // calculating concentrations from quasi fermi levels
        cont.conc[car][i]=ni*exp(expdif*cont.coefexp);  // calculating concentrations from quasi fermi levels
      }
    }
  }

/*  if(rp&rpFermi){
    for(int car=0;car<2;car++){
      valtype cnc=std::log(cont.conc[car][i]/ni)/cont.coefexp;
      cont.fermi[car][i] = car ? psi+cnc : psi-cnc;
    }
  }*/
}

#endif
