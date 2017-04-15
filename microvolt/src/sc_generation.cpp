#include "sc_generation.h"
#include "string_utils.h"
#include "physconst.h"

using namespace phys;

scLambert::scLambert(Plane_3 pl_, valtype I_, valtype wv, const char *dir, const char *fabs):pl(pl_),I(I_){

  string cdir=string(dir);
  correct_directory_name(cdir);

  char dfabs[1000];
  strcpy(dfabs,cdir.c_str());
  strcat(dfabs,fabs);

  valtype *aval;
  int an1,an2;
  if(read_table(dfabs,&aval,an1,an2,false,true)<0)
    throw 1;

  if(an2<3)
    throw 1;

  valtype k=get_table_value(wv,aval,an1,an2,0,2,false);
  wv*=1e-9; // translate to meters
  I/=phys::h*phys::c/wv; // J/sm2 / (Jsm/s/m ) = 1/sm2 (sm means sec*meter)
  valtype omega=2.*M_PI/wv;
  alpha=2*omega*k;

  delete[]aval;
}

scLambertSpectrum::scLambertSpectrum(Plane_3 pl_, const char *dir, const char *flight,
const char *fabs, valtype wv1, valtype wv2, valtype Imult_, const char *filter, int jfilter):
I(NULL),alpha(NULL),Imult(Imult_),pl(pl_){

  string cdir=string(dir);
  correct_directory_name(cdir);

  char dflight[1000];
  char dfabs[1000];
  strcpy(dflight,cdir.c_str());
  strcat(dflight,flight);
  strcpy(dfabs,cdir.c_str());
  strcat(dfabs,fabs);

  valtype *lval;
  int ln1,ln2;
  if(read_table(dflight,&lval,ln1,ln2,false,true)<0)
    throw 1;

  valtype *aval;
  int an1,an2;
  if(read_table(dfabs,&aval,an1,an2,false,true)<0)
    throw 1;

  valtype *fval=NULL;
  int fn1=0,fn2=0;
  if(filter && read_table(filter,&fval,fn1,fn2,false,true)<0)
    throw 1;

  if(ln2<2 || an2<3)
    throw 1;

  int i=0, j=-1;
  for(;i<ln1;i++){
    if(lval[i]<wv1)
      continue;
    else if(j<0)
      j=i;
    if(wv2<lval[i])
      break;
  }
  if(j<0 || i==j)
    throw 1;
  num=i-j;
  I = new valtype[num];
  alpha = new valtype[num];
  valtype Itot=0; // total number of photons

  valtype *I_=lval+1*ln1;
  valtype *wv=lval;
  for(i=j;i<num+j;i++){
    valtype k=get_table_value(wv[i],aval,an1,an2,0,2,false);

    I[i-j]=I_[i];
    I[i-j]*=(wv[i+1]-wv[i]);
    valtype fabs = filter ? get_table_value(1/wv[i],fval,fn1,fn2,0,jfilter,false) : 1;
    wv[i]*=1e-9; // translate to meters
    I[i-j]/=phys::h*phys::c/wv[i]; // J/sm2 / (Jsm/s/m ) = 1/sm2 (sm means sec*meter)
    Itot+=I[i-j]*fabs;

    valtype omega=2.*M_PI/wv[i];
    alpha[i-j]=2*omega*k;
  }
  valtype max=Itot*phys::q;
  delete[]lval;
  delete[]aval;
}

int scGenerationTable::pack(const int *N, const int *ord, const int *ic, int argnum){
  int mult=1, ind=0;
  for(int i=0;i<argnum;i++){
    int dim = ord ? ord[i] : i;
    ind+=mult*ic[dim];
    mult*=N[dim];
  }
  return ind;
}

scGenerationTable::scGenerationTable(const char *dir, const char *fgen, 
  const Vector_3 &p1, const Vector_3 &p2, const int *N, valtype Imult_): 
  gen(NULL),Imult(Imult_){

  grd.init(p1,p2,2,N);

  string cdir=string(dir);
  correct_directory_name(cdir);

  char dfgen[1000];
  strcpy(dfgen,cdir.c_str());
  strcat(dfgen,fgen);

  valtype *val;
  int n1,n2;
  if(read_table(dfgen,&val,n1,n2,true)<0)
    throw 1;

  if(N[0]*N[1]*N[2] != n1*n2)
    throw 1;

  gen = new valtype[N[0]*N[1]*N[2]];
  int ic[3];
  for(ic[0]=0;ic[0]<N[0];ic[0]++){
    for(ic[1]=0;ic[1]<N[1];ic[1]++){
      for(ic[2]=0;ic[2]<N[2];ic[2]++){
        int ind=grd.pack_ind(ic[0],ic[1],ic[2]);
        int vind=pack(N,NULL,ic,3);
        gen[ind]=val[vind];
//        gen[ind]=1;
      }
    }
  }
}

int scGenerationTable::init(const char *dir, const char *fgen, 
  const iVector_3 &arg, int jval, int arg_num, 
  const Basis_3 &bs_, const Vector_3 &shift_, valtype Imult_, int format){

  bs=bs_;
  shift=shift_;
  Imult=Imult_;

  string cdir=string(dir);
  correct_directory_name(cdir);

  char dfgen[1000];
  strcpy(dfgen,cdir.c_str());
  strcat(dfgen,fgen);

  double *ptr;
  int N1,N2;

  if(read_table(dfgen,&ptr,N1,N2,!(format&TABLE_HORIZONTAL),format&TABLE_HEADER)<0)
    return message(vblMESS1,-1,"Cannot read generation table %s\n",fgen);

  Vector_3 p1,p2;
  int sz[3],ord[3];
  if(table_grid(ptr,N1,N2,arg_num,arg.v,p1.v,p2.v,sz,ord)<0){ // read dimension of the table (p1, p2, sz)
    return message(vblMESS1,-1,"Cannot read generation table %s\n",fgen);
  }
  for(int i=2;i>=arg_num;i--){ // fill dimension that is not present in the table
    p1[i]=p2[i]=0;
    sz[i]=1;  
  }
  grd.init(p1,p2,2,sz);

  if(gen)delete[] gen;
  gen = new valtype[N1];
  int ic[3];
  for(ic[0]=0;ic[0]<sz[0];ic[0]++){ // fastest
    for(ic[1]=0;ic[1]<sz[1];ic[1]++){
      for(ic[2]=0;ic[2]<sz[2];ic[2]++){ // slowest

        int ind=grd.pack_ind(ic[0],ic[1],ic[2]);
        int vind=pack(sz,ord,ic,arg_num);
        gen[ind]=ptr[jval+vind*N2]*Imult;
//        gen[ind]=1e10;
/*        if(gen[ind])
          gen[ind]=0;
        else
          gen[ind]=1e21;*/
      }
    }
  }
  delete[]ptr;

  return 1;
}

valtype scGenerationTable::operator()(const Vector_3 &p) {

  Vector_3 pos=bs(p+shift); // project to the generation mesh

  int gind[8];
  valtype coef[8];
  int n=grd.GetCoeff(pos,gind,coef,1,NULL); // getting interpolation coefficient and indices

  // interpolation of generation value using values on closest grid points
  valtype coef_nz=0;
  for(int i=0;i<n;i++){
    int ind=gind[i];
    if(gen[ind])
      coef_nz+=coef[i];
  }
  if(!coef_nz)
    return 0;
  coef_nz=1;
  valtype val=0;
  for(int i=0;i<n;i++){
    int ind=gind[i];
    val+=gen[ind]*coef[i]/coef_nz;
  }
  return val;
}

valtype scGenerationCylinder::operator()(const Vector_3 &pos) {
  int sym1=(sym+1)%3;
  int sym2=(sym+2)%3;

  // test parameters
  valtype phist=0; // starting angle
  int st=0; // if phist will be used

  Vector_3 p=pos;
  p[sym]=center[sym];
  valtype R=(p-center).norm();
  valtype phi0 = st ? asin((p-center)[sym1]/R) : phist;

  p[sym]=pos[sym];
  valtype res=0;

  for(int i=0;i<num;i++){
    valtype phi=valtype(i)/num;
    phi*=2*M_PI;
    p[sym1] = center[sym1]+R*sin(phi0+phi);
    p[sym2] = center[sym2]+R*cos(phi0+phi);
    res+=(*gen)(p);
  }
  res/=num;
  return res;
}
