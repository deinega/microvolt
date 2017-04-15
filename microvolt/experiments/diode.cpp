#include "sc_ui.h"

/*
This is 1d diode test code. This code is used for comparison with PC1D and AFORS-HET.

1,.5 -1e-7 .01 - comparison with PC1D and AFORS-HET (rad=2 should be used),
see pc1d_si_diode.prm and afors_het_si_diode.het in transport/benchmark
there is exact agreement with AFORS-HET and good agreement with PC1D
in AFORS-HET photon flux 2.515e17 1/(cm^2 s) should be used

1,.01 1,1e5 .001 - comparison with Kayes paper
(B. Kayes et al.,'Comparison of the device physics principles of planar and radial p-n junction nanorod solar cells')
if diffusion length increase then surface recombination influence on Isc and Voc increases

1,.5 1,1e5 .01 - test to check correct work after code changed (depl should be turned on)
*/

int main(int argc,char **argv){

  string dir,spectra_dir;
  dir="../../../../runs/diode"; // output directory
  spectra_dir="../../spectra";

  if(scInit(argc,argv,dir)<0)
    return -1;

  /* Length is d
  left medium at (0, left_length) and right medium at (left_length, length) 
  grad_left_right defines intermediate alloy region for heterojunction case
  if 0, than there is not intermediate region (abrupt heterojunction)
  */

  valtype length=1; // diode width, micron
  valtype left_length=.5; // pn-junction position, if pn=-1 than pn=d/2
  valtype grad_left_right=.2; //activated only when using corresponding medium-factory
  valtype left_dop=1e18; // doping concentration of the left medium, cm-3
  valtype right_dop=-1e18; // doping concentration of the right and back medium. if right_dop=-1 than right_dpop=-left_dop

  scMediumFactory left_med, right_med;

  valtype L=-1e-7; // if negative - lifetime, sec; if positive - diffusion length, micron
  valtype sr=VEC_INFTY; // surface recombination, cm/s
  valtype depl=1*0.05; // large timelife and small recombination at depl around junction
  int rad=3; // 0 - no, 1 - uniform, 2 - 500nm, 3 - AM1.5
  valtype dr=0.01; // mesh size
  int wd=2; // working dimension (z as a default)

  int ind=1;
  int ac=theConfig->argv.size();
  if(ac>ind)extract_list(theConfig->argv[ind++],&length,&left_length);
  if(ac>ind)extract_list(theConfig->argv[ind++],&L,&sr);
  if(ac>ind)extract_list(theConfig->argv[ind++],&dr);

  string tname=fmt("%g\t%g\t%g",length,L,sr);

  left_med=getscSi();
  left_med.dop=left_dop;
//  left_med.N[0]=left.N[1]=0;
//  left_med.T=270;
//  left_med.ni=1e11;
//  left_med.mu[0]=270, left_med.mu[1]=95; // this are mobilities from Kayes paper
//  left_med=getscGe(); // germanium
  right_med=left_med; // silicon
  right_med.dop=right_dop;
//  right_med.Eg+=.2; // test when I was deleting old medium
//  right_med.chi+=.1; // test when I was deleting old medium
//  right_med.chi+=0.2;
//  right_med.Eg+=0.2;

//  left_med=getscCdS();
//  right_med=getscCdTe();


  if(L>0)
    left_med.L[0]=right_med.L[0]=L;
  else
    left_med.rcm.tau[0]=right_med.rcm.tau[0]=-L;
  Vector_3 sz(0);
  sz[wd]=length;
  scExperiment task(0,sz);
//  task.SetLengthUnit(1e6);

  task.set_name(tname);
  task.DumpMeshes(); // meshes will be dump their values
//  task.DumpEachIteration();

  task.SetResolution(dr);
//  task.SetMethod("p");
//  task.SetMethod("p_eh","1,10");
//  task.SetMethod("eh,p","10,1");
//  task.SetMethod("p_e_h_pe_ph");
//  task.SetMethod("pe,ph,hP,eP");
//  task.SetIfl(1);
  task.SetDepletionDecreasing(depl);
//  task.SetGenEq(1);

  // mesh generation
  vector<Vector_3> inc;
//  left_length=.5;
//  inc.push_back(Vector_3(left_length-10*dr,left_length+10*dr,dr/10)); // more resolution at pn-junction
//  inc.push_back(Vector_3(0,3*dr,dr/10)); // more resolution at contacts
//  inc.push_back(Vector_3(length-3*dr,length,dr/10)); // more resolution at contacts
  inc.push_back(Vector_3(0,2*left_length,dr/10)); // more resolution at pn-junction
  vector<valtype> mesh;
  if(generate_1d_mesh(Vector_3(0,length,dr),inc,mesh)<0)
    return -1;
//  task.SetMeshSteps(wd,mesh);

/*  valtype *zval;
  int zn1,zn2;
  if(read_table("D:\\FDTD_Work\\transport\\experiments\\zval_1mkm.txt",&zval,zn1,zn2)<0)
//  if(read_table("D:\\FDTD_Work\\transport\\experiments\\zval_kayes_1mkm.txt",&zval,zn1,zn2)<0)
    return -1;
//  task.SetMeshSteps(2,zval,zval+zn1); // using non-uniform rectangular mesh
  delete[]zval;*/

  // using 2 meshes
  // default mesh is confined by left side
//  task.SetConfinement(0,GetHalfSpace(Vector_3(0,0,-1),Vector_3(0,0,d/2+4*d/nz)));
  // right side mesh is added
//  task.AddBlock(iVector_3(1,1,2*nz),GetHalfSpace(Vector_3(0,0,1),Vector_3(0,0,d/2-0*d/nz)),1);

  // adding dense mesh at the depletion region
//  int sd=5;
//  task.AddBlock(iVector_3(1,1,10*nz),GetHalfSpace(Vector_3(0,0,1),Vector_3(0,0,pn-sd*d/nz),2*sd*d/nz),1);

  Vector_3 n;
  n[wd]=1;
  // abrupt junction
  task.AddMediumRegion(left_med,GetHalfSpace(-n,left_length*n));
  task.AddMediumRegion(right_med,GetHalfSpace(n,left_length*n));
//#ifndef NEW_MEDIUM
//  task.AddMediumRegion(left_med.create(),GetHalfSpace(-n,left_length*n));
//  task.AddMediumRegion(right_med,GetHalfSpace(n,left_length*n));
//#endif
  // gradual junction
  //scParametricMediumFactory::PositionFunction *fun(new transition<Polyhedron_3>(GetHalfSpace(-n,(left_length-dr/2)*n),.5*grad_left_right));
  //scVarDopFactory cont_med(left_med,right_dop,fun); // varying doping
  //scAlloyFactory<> cont_med(left_med,right_med,
    //SHARED_OR_AUTO(scParametricMediumFactory::PositionFunction)(fun)); // alloy
  //task.AddMediumRegion(cont_med,new Region_3());

  task.AddContact(GetHalfSpace(-n,VEC_ZERO*n),sr,sr);
  // voltage will be applied to this boundary 
  valtype wf=VEC_INFTY;
//  valtype wf=5.42;
  task.AddVoltageContact(GetHalfSpace(n,(length-VEC_ZERO*dr)*n),sr,sr,wf);

  if(rad==1){ // uniform absorption (unphysical)
    task.AddGeneration(new scGenerationUniform(1e20));
/*    valtype c=.9e-6;
    valtype w=.1e-6;
    Box *reg = new Box(Vector_3(0,0,c-w/2),Vector_3(1,1,c+w/2));
    task.AddGeneration(new scGenerationUniform(1e20*10),reg);*/
//    task.AddGeneration(new transition<Polyhedron_3>(GetPolyhedronPlane(Vector_3(0,0,1),Vector_3(0,0,d/2)),.1e-6,1e20));
  }
  else if(rad==2) // monochromatic absorption
    task.AddGeneration(new scLambert(Plane_3(n,Vector_3()),1e3,500,spectra_dir.c_str(),"si.nk"));
  else if(rad==3){ // solar spectrum absorption
    scGeneration *genAM = new scLambertSpectrum(Plane_3(n,Vector_3()),spectra_dir.c_str(),
//      "am1.5_afors_het.dat","si.nk",350,1240.8/si.Eg);
      "am1.5.in","si.nk",350,1240.8/left_med.Eg);
    task.AddGeneration(genAM);
  }
  else if(rad==4){ // resubstituting absorption from previous Microvolt simulation
    Basis_3 bs(Vector_3(0,0,1),Vector_3(0,1,0),Vector_3(1,0,0));
    // last argument is 1e6, since in Microvolt output generation is in [1/sec/cm3],
    // while scGenerationTable read generation in [1/sec/m3]
    scGeneration *gen=new scGenerationTable("..\\..\\..\\..\\..\\runs\\diode5","d.d",0,10,1,bs,0,1e6);
    task.AddGeneration(gen);
  }

  iVector_3 dsz(1);
//  dsz[wd]=int(acdiv(d,dr))+1;
//  task.AddDetector("d",0,sz,dsz);
  task.AddDetector("d",0,sz,iVector_3(INT_INFTY));
//  task.AddBoxDetector("j",0,sz,iVector_3(1,1,1),BOX_FRONT_Y);
  task.AddFluxPlane("j",wd,length);
  task.SetFluxIV("j"); // I-V calculation

  task.Calculate(.02);

  return 1;
}
