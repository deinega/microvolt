#include "sc_ui.h"

/*
This is 1d silicon diode test code.
Diod can be rotated in 2D space.
Additional higher resolution mesh can be added.
*/

int main(int argc,char **argv){

  string dir="../../../../runs/diode2";
  if(scInit(argc,argv,dir)<0)
    return -1;

  valtype z=10; // length
  int nz=50;
  valtype dz=z/nz;
  int ny=10;
  valtype y=ny*dz; // width
  valtype alpha=0*(M_PI/2.);
//  valtype alpha=1.*(M_PI/2.);
//  valtype alpha=1./2*(M_PI/2.);
//  valtype alpha=1./3.5*(M_PI/2.);
  valtype y1=z*sin(alpha)+y*cos(alpha);
  valtype z1=z*cos(alpha)+y*sin(alpha);

  int ny1=int(ceil(accdiv(y1,dz)))+2;
  int nz1=int(ceil(accdiv(z1,dz)))+2;
//  int ny1=ny; // for exact compare with 1D
//  int nz1=nz; 

  y1=ny1*dz;
  z1=nz1*dz;

  Vector_3 sz(2,y1,z1);

  scMedium si=getscSi();
  si.L[0]=1;
  valtype dop=1e18;

  scExperiment task;
  task.DumpMeshes();
  task.SetInternalSpace(0,sz);
  task.SetResolutionN(iVector_3(1,ny1,nz1));
//  task.SetResolutionN(iVector_3(1,3*ny1,3*nz1));
  task.SetMemoryFlag(1);

  SpaceRegion *conf;
//  conf=GetPolyhedronPlane(Vector_3(0,0,-1),Vector_3(0,0,z1/2));
//  conf=new Inverse<SpaceRegion>(make_mngarg(new Sphere(2*dz,Vector_3(1,y1/2,z1/4)))); // alpha 0
//  conf=(new Inverse<SpaceRegion>(make_mngarg(new Sphere(1*dz,Vector_3(1,2.6*y1/4,z1/4)))));
//  task.SetConfinement(conf);

  conf=GetPolyhedronPlane(Vector_3(0,0,-1),Vector_3(0,0,z1/2-2*dz));
//  conf=new Sphere(2*dz,Vector_3(1,1*y1/2,z1/4));
//  conf=new Sphere(2*dz,Vector_3(1,2.6*y1/4,z1/4));
//  task.AddBlock(iVector_3(1,2*ny1,2*nz1),NULL,1);
//  task.AddBlock(iVector_3(1,7*ny1,7*nz1),NULL,1);
//  task.SetConfinement(conf);

  Vector_3 nvect(0,sin(alpha),-cos(alpha));
  Vector_3 pvect(0,-cos(alpha),-sin(alpha));
  Vector_3 sep(0,dz+z/2*sin(alpha),dz+z/2*cos(alpha));
//  Vector_3 sep(0,z/2*sin(alpha),z/2*cos(alpha)); // for exact compare with 1D

  task.AddMediumRegion(si,dop,GetPolyhedronPlane(nvect,sep));
  task.AddMediumRegion(si,-dop,GetPolyhedronPlane(-nvect,sep));

  task.AddReflectiveBoundary(GetPolyhedronPlane(pvect,sep-pvect*(dz*.5)),-1);
  task.AddReflectiveBoundary(GetPolyhedronPlane(-pvect,sep-pvect*(y-dz*.5)),-1);

//  valtype sr=VEC_INFTY;
  valtype sr=1e5;
  task.AddVoltageContact(GetPolyhedronPlane(-nvect,sep-nvect*z/2),sr,sr,VEC_INFTY,-2);
  task.AddContact(GetPolyhedronPlane(nvect,sep+nvect*z/2),sr,sr,VEC_INFTY,-2);

  task.AddDetector("d",0,sz,iVector_3(1,ny1+1,nz1+1));

  task.Calculate(.02);

  return 1;
}
