#include "sc_ui.h"

/*
Multijunction solar cells.

1,.2,1,.2 .2,1e5 .02
*/

int main(int argc,char **argv){

  string dir,spectra_dir;
  dir="../../../../runs/diode"; // output directory
  spectra_dir="../../spectra";

  if(scInit(argc,argv,dir)<0)
    return -1;

  valtype length=2; // diode vertical length
  valtype w=.2; // pn-junction position along diode length
  valtype a=1; // period
  valtype left_length=.4; // width of the layer
  valtype r=0; // internal radius
  int hp=2; // number of half-periods in calculated box
  int sym=0; // use radial coordinates

  valtype ndop=1e18; // n doping concentration, cm-3
  valtype pdop=-1e18; // p doping concentration
  valtype L=1; // diffusion length, mkm
  valtype sr=1e5; // surface recombination, cm/sec
  valtype sri=1e2; // insulator surface recombination, cm/sec
  valtype depl=0; // large timelife and small recombination at depl around junction
  int rad=3; // 0 - no, 1 - uniform, 3 - AM1.5
  valtype dz=.02;
  int nun=1;
//  valtype dz=.05;
//  valtype nun=2.5; // if >1, smaller mesh step is using around p-n junction

  int ind=1;
  int ac=theConfig->argv.size();
  if(ac>ind)extract_list(theConfig->argv[ind++],&length,&w,&a,&left_length,&r);
  if(ac>ind)extract_list(theConfig->argv[ind++],&L,&sr,&sri);
  if(ac>ind)extract_list(theConfig->argv[ind++],&dz,&nun);

  if(sr<0)
    sr=VEC_INFTY;

  scMediumFactory left=getscSi(); // silicon
  scMediumFactory right=left;
  left.dop=ndop;
  right.dop=pdop;
  left.L[0]=right.L[0]=L;

  Vector_3 sz(0,hp*a/2,length);

  scExperiment task(0,sz);
  valtype dy=1*dz;
  task.SetResolution(Vector_3(0,dy,dz));

  task.set_name(fmt("%g\t%g",length,L));
  task.DumpMeshes(1|3); // to see values of calculated variables at mesh nodes

  task.SetDepletionDecreasing(depl);

  if(nun>1){

    RegionUnion<3> *reg = new RegionUnion<3>();

    reg->AddRegion(GetBox(Vector_3(-1,-left_length/2,length-w+dz),Vector_3(1,left_length/2,length-w-dz)));
    for(int i=0;i<hp;i+=2){
      reg->AddRegion(GetBox(Vector_3(-1,i*a/2+left_length/2-dy,length-w),Vector_3(1,i*a/2+left_length/2+dy,w)));
      reg->AddRegion(GetBox(Vector_3(-1,i*a/2+left_length/2,w+dz),Vector_3(1,i*a/2+a-left_length/2,w-dz)));
    }
    for(int i=1;i<hp;i+=2){
      reg->AddRegion(GetBox(Vector_3(-1,(i+1)*a/2-left_length/2-dy,length-w),Vector_3(1,(i+1)*a/2-left_length/2+dy,w)));
      reg->AddRegion(GetBox(Vector_3(-1,(i+1)*a/2-left_length/2,length-w+dz),Vector_3(1,(i+1)*a/2+left_length/2,length-w-dz)));
    }

    // put local mesh with higher resolution and priority 1 inside region reg
    task.AddMesh(Vector_3(0,dy,dz)/nun,reg,1);
//    task.AddMesh(Vector_3(0,dy,dz)/4,GetHalfSpace(Vector_3(0,0,-1),w/2),1);

    vector<Vector_3> inc;
    inc.push_back(Vector_3(left_length/2-dy < 0 ? 0 : left_length/2-dy,left_length/2+dy,dy/nun)); // more resolution at pn-junction
    vector<valtype> mesh;
    if(!accomp(w,length-w)){
      if(generate_1d_mesh(Vector_3(0,sz[1],dy),inc,mesh)<0)
        return -1;
//      task.SetMeshSteps(1,mesh.begin(),mesh.end());
    }

    vector<Vector_3> incx;
    incx.push_back(Vector_3(length-w-dz,length-w+dz > length ? length : length-w+dz,dz/nun)); // more resolution at pn-junction
    if(!accomp(w,length-w))
      incx.push_back(Vector_3(w-dz < 0 ? 0 : w-dz,w+dz,dz/nun)); // more resolution at pn-junction
//    incx.push_back(Vector_3(0,dy,dy/2/nun)); // more resolution at contact
//    incx.push_back(Vector_3(length-dy,length,dy/2/nun)); // more resolution at contact
    vector<valtype> meshx;
    if(generate_1d_mesh(Vector_3(0,length,dz),incx,meshx)<0)
      return -1;
//    task.SetMeshSteps(2,meshx.begin(),meshx.end());
  }

  task.AddMediumRegion(right,new Region_3(),-1); // default medium
  if(accomp(w,length-w)) // planar diode case
    task.AddMediumRegion(left,GetHalfSpace(Vector_3(0,0,-1),w));
  else{
    for(int i=0;i<=hp;i+=2)
      task.AddMediumRegion(left,GetPlate(Vector_3(0,1,0),-.5*left_length+i*a/2,left_length));
    task.AddMediumRegion(left,GetHalfSpace(Vector_3(0,0,-1),w),1);
    task.AddMediumRegion(right,GetHalfSpace(Vector_3(0,0,1),length-w),1);
  }

  Region_3 *LR = GetHalfSpace(Vector_3(0,-1,0),VEC_ZERO);
  Region_3 *RR = GetHalfSpace(Vector_3(0,1,0),sz[1]-VEC_ZERO);
  Region_3 *BR = GetHalfSpace(Vector_3(0,0,-1),VEC_ZERO);
  Region_3 *UR = GetHalfSpace(Vector_3(0,0,1),length-VEC_ZERO);

  if(sym){
    // y is a radial coordinate, hp*a/2+r is coordinate of the radial axis,
    // we consider negative radial values
    task.SetRadialCoordinate(1,hp*a/2+r,false);
    task.AddReflectiveBoundary(LR,sri,sri,-1);
    if(r)
      task.AddReflectiveBoundary(RR,sri,sri,-1);
  }
  else{
//    task.AddReflectiveBoundary(LR,0,0,-1);
//    task.AddReflectiveBoundary(RR,0,0,-1);
    task.SetPeriodicBoundaries(0x2); // set periodicity for y-dimension
  }

  task.AddContact(BR,sr,sr,VEC_INFTY);
  task.AddVoltageContact(UR,sr,sr,VEC_INFTY);

/*  valtype cwidth = .2; // dot contact width
  Polyhedron_3 *conf = GetPlate(Vector_3(0,1,0),(sz[1]-cwidth)/2,cwidth);
  task.AddVoltageContact(new ConfinedRegion<Region_3>(BR,conf),sr,sr,VEC_INFTY,1);
  task.AddReflectiveBoundary(BR,sri,sri);
  task.AddContact(new ConfinedRegion<Region_3>(UR,conf),sr,sr,VEC_INFTY,1);
  task.AddReflectiveBoundary(UR,sri,sri);*/


/*
// playing with contact position and size
//  valtype sdz=x/2;
  valtype sddz=.05;
  sddz*=1e-6;
  valtype sdz=sddz;
  Polyhedron_3 *SD = GetBox(Vector_3(sdz,0,1e-15),Vector_3(sdz-sddz,1,-1e20));
//  task.AddContact(SD,sr,sr,2);
*/

  task.AddDetector("d",0,sz,iVector_3(1, hp*20, 40));

  task.AddFluxPlane("j",2,sz[2]/2);

  if(rad==1){ // uniform absorption (nonphysical)
    task.AddGeneration(new scGenerationUniform(1e20));
  }
  else if(rad==3){
    task.AddGeneration(new scLambertSpectrum(Plane_3(Vector_3(0,0,1),0),spectra_dir.c_str(),
      "am1.5.in","si.nk",350,1240.8/left.Eg));
  }

  task.SetFluxIV("j");
  task.Calculate(.02);

  return 1;
}
