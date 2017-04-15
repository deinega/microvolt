#include "sc_ui.h"

/*
This is 1d silicon diode test in 3D.
Diod can be rotated in 3D.
*/

int main(int argc,char **argv){

  string dir,spectra_dir;

  dir="../../../../runs/diode_3d"; // output directory
  spectra_dir="../../spectra";

  if(scInit(argc,argv,dir)<0)
    return -1;

  valtype d=1; // diode length, mkm
  valtype pn=0.5;
  valtype left_dop=1e18; // doping concentration of the left medium, cm-3
  valtype right_dop=-1e18; // doping concentration of the right medium.
  scMediumFactory med=getscSi();
  valtype L=1;
  valtype sr=1e5; // surface recombination, cm/s
//  valtype sr=VEC_INFTY; // surface recombination, cm/s
  valtype sri=0;
  valtype depl=0;
  int rad=3; // 0 - no, 1 - uniform, 2 - 500nm, 3 - AM1.5

  valtype dr=0.01;
  Vector_3 bulk(d,10*dr,0*dr);

  valtype theta=M_PI/3, phi=0*M_PI/3;
  Vector_3 n_vect[3];
  n_vect[0]=Vector_3(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta)); // diode axis
//  n_vect[1]=Vector_3(1,0,0);
//  n_vect[2]=Vector_3(0,1,0); // perpendicular directions
  build_orth_basis(n_vect[0],n_vect[1],n_vect[2]);

  med.L[0]=L; 

  // calculating bounding box
  Vector_3 p1(VEC_INFTY),p2(-VEC_INFTY);
  valtype add=0*dr;
  for(int i=0;i<0x10;i++){
    Vector_3 p;
    for(int j=0;j<3;j++){
      if(bulk[j]){
        if(i&(0x1<<j))
          p+=(bulk[j]+add)*n_vect[j];
        else
          p-=add*n_vect[j];
      }
    }
    for(int j=0;j<3;j++){
      if(p[j]<p1[j])
        p1[j]=p[j];
      if(p[j]>p2[j])
        p2[j]=p[j];
    }
  }

  scExperiment task(p1,p2);
  task.DumpMeshes(1|2);

  Vector_3 c=(p1+p2)/2;

  task.SetResolution(dr);
//  task.SetPhases("g");

//  task.SetMethod("eh","1",0);
//  task.SetMethod("eh,p","10,1",0);
//  task.SetMethod("e,h,p","1,1,1");
//  task.SetMethod("eh,p","10,1");
//  task.SetMethod("peh,eh,p","10,10,1");
//  task.SetMethod("peh","1",0);
//  task.SetGenEq(1);
  task.SetDepletionDecreasing(depl);

  task.AddMesh(dr/2,GetHalfSpace(Vector_3(0,0,1),.3),1);
//  task.AddMesh(dr/2,GetHalfSpace(Vector_3(0,1,0),.91),1);
//  task.AddMesh(dr*2,GetHalfSpace(Vector_3(0,0,1),.3),-1);

  med.dop=left_dop;
  task.AddMediumRegion(med,GetHalfSpace(-n_vect[0],c));
  med.dop=right_dop;
  task.AddMediumRegion(med,GetHalfSpace(n_vect[0],c));

  for(int i=1;i<3;i++){
    if(bulk[i]){
//      task.AddContact(GetHalfSpace(-n_vect[i],c-n_vect[i]*(bulk[i]/2-VEC_ZERO)),-1);
//      task.AddContact(GetHalfSpace(n_vect[i],c+n_vect[i]*(bulk[i]/2-VEC_ZERO)),-1);
      task.AddReflectiveBoundary(GetHalfSpace(-n_vect[i],c-n_vect[i]*(bulk[i]/2-VEC_ZERO)),sri,sri,-1);
      task.AddReflectiveBoundary(GetHalfSpace(n_vect[i],c+n_vect[i]*(bulk[i]/2-VEC_ZERO)),sri,sri,-1);
    }
  }

  task.AddContact(GetHalfSpace(-n_vect[0],c-(d/2-VEC_ZERO)*n_vect[0]),sr,sr);
  task.AddVoltageContact(GetHalfSpace(n_vect[0],c+(d/2-VEC_ZERO)*n_vect[0]),sr,sr);

  if(rad==2){
    task.AddGeneration(new scLambert(Plane_3(n_vect[0],0),1e3,500,spectra_dir.c_str(),"si.nk"));
  }
  else if(rad==3){
    scGeneration *genAM = new scLambertSpectrum(Plane_3(n_vect[0],0),spectra_dir.c_str(),
      "am1.5.in","si.nk",350,1240.8/med.Eg);
    task.AddGeneration(genAM);
  }

  task.Calculate(.02);

  return 1;
}
