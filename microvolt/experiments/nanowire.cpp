#include "sc_ui.h"

int main(int argc,char **argv){

  string dir,spectra_dir,fdtd_dir;
  dir="../../../../runs/diode"; // output directory
  spectra_dir="../../spectra";
  fdtd_dir="../../../../runs"; // directory with fdtd generation file

  if(scInit(argc,argv,dir)<0)
    return -1;

  valtype R=.2; // nanowire radius
  valtype h=2; // nanowire height
  valtype uw=.1, dw=.1; // uw,dw - pn-junction distance from contact, if dw is negative, then dw=uw
  valtype pn=R/2; // pn-junction distance from nanowire surface (if negative, then pn will be equal to R/2)
  valtype ndop=1e18; // n-doping
  valtype pdop=1e18; // p-doping (if negative, then pdop will be the same as ndop)
  valtype L=.2; // diffusion length, um
  valtype sr=1e5; // contact surface recombination, cm/sec
  valtype sri=0; // insulator surface recombination, cm/sec
  valtype dx=.004; // mesh step along x-axis
  valtype nz_mult=1; // ratio dz/dx
  int nun=0; // if nonuniform mesh step will be used
  int mg=0; // if extra-grid wull be used
  int rad=1; // use radial coordinates

  int nx; // number of mesh points per 2*R
  if(dx>1){
    nx=int(dx);
    dx=2*R/nx;
  }
  else{
    nx=int(acdiv(2*R,dx));
    dx=acdiv(2*R,valtype(nx));
  }

  int nz=int(acdiv(h,nz_mult*dx)); // number of mesh points per h

  scMediumFactory top=getscSi(); // silicon
  // mobilities from Kayes paper
  // top.mu[0]=270, top.mu[1]=95;
  top.L[0]=L; // um
  scMediumFactory down=top;
  top.dop=ndop; down.dop=-pdop;

  // calculated box.
  Vector_3 p1(-R-dx,-R-dx,0), p2(R+dx,R+dx,h);
  if(rad)
    p1[1]=p2[1]=0; // y-axis is not working

  scExperiment task(p1,p2);
  task.SetResolutionN(iVector_3(nx+2,nx+2,nz));

  if(rad){
    task.SetRadialCoordinate(0,0,0); // x is radial axis
    // reflective boundary will be applied instead of axis (nanowire -> vertical multijunction)
//    task.AddReflectiveBoundary(GetCylinder(R,Vector_3(0,1,0),1e-9),-9);
  }

  if(mg && nun>1){ // put extra-mesh of the higher resolution at the pn-junction area
    task.AddMeshN(iVector_3(nun*(nx+2),nun*(nx+2),nz),NULL,1);
    RegionIntersection<3> *reg_dr = new RegionIntersection<3>;
    reg_dr->AddRegion(GetCylinder(Vector_3(0,0,0),Vector_3(0,0,1),R-pn+dx));
    reg_dr->AddRegion(new Inverse<Region_3>(make_mngarg(GetCylinder(Vector_3(0,0,0),Vector_3(0,0,1),R-pn-dx))));
    task.SetConfinement(reg_dr);
//  task.SetConfinement(GetCylinder(Vector_3(0,0,0),Vector_3(0,0,1),R-1*dr));
  }

  if(nun>1){ // use smaller mesh step at the pn-junction (but mesh is the same)
    if(!mg){
      vector<Vector_3> inc;
      // more resolution at the pn-junction
      inc.push_back(Vector_3(-pn-2*R/nx/2 < 0 ? 0 : -pn-2*R/nx/2,-pn+2*R/nx/2,2*R/nx/nun));
      vector<valtype> mesh;
      if(!accomp(dw,h-uw)){
        if(generate_1d_mesh(Vector_3(-R-dx,dx,2*R/nx),inc,mesh)<0)
          return -1;
        task.SetMeshSteps(0,mesh.begin(),mesh.end());
      }
    }
    vector<Vector_3> incz;
    incz.push_back(Vector_3(h-uw-h/nz,h-uw+h/nz > h ? h : h-uw+h/nz,h/nz/nun)); // more resolution at the contact
    if(!accomp(dw,h-uw))
      incz.push_back(Vector_3(dw-h/nz < 0 ? 0 : dw-h/nz,dw+h/nz,h/nz/nun)); // more resolution at the contact
    vector<valtype> meshz;
    if(generate_1d_mesh(Vector_3(0,h,h/nz),incz,meshz)<0)
      return -1;
/*    task.SetMeshSteps(2,meshz.begin(),meshz.end(),0);
    if(mg)
      task.SetMeshSteps(2,meshz.begin(),meshz.end(),1);*/
  }

  Cylinder<Circle> *OC = GetCylinder(0,Vector_3(0,0,1),R-1e-10);
  Inverse<Region_3> *IOC = new Inverse<Region_3>(make_mngarg((Region_3 *)OC)); // outer cylinder

  Polyhedron_3 *DP=GetHalfSpace(Vector_3(0,0,-1),Vector_3(0,0,1e-10)); // bottom plane
  Polyhedron_3 *UP=GetHalfSpace(Vector_3(0,0,1),Vector_3(0,0,h-1e-10)); // top plane

  task.AddReflectiveBoundary(IOC,sri,sri,-3); // outer cylinder is insulator boundary
  task.AddVoltageContact(DP,sr,sr,VEC_INFTY,-1); // bottom plane is contact
  task.AddContact(UP,sr,sr,VEC_INFTY,-1); // top plane is contact

  ConfinedRegion<Cylinder<Circle> > *pnC = GetCylinder(0,Vector_3(0,0,1),R-pn,h-uw);
  task.AddMediumRegion(down,pnC); // p-region

  task.AddMediumRegion(down,GetHalfSpace(Vector_3(0,0,-1),Vector_3(0,0,dw)),1); // p-region
  if(accomp(dw,h-uw))
    task.AddMediumRegion(top,GetHalfSpace(Vector_3(0,0,1),Vector_3(0,0,dw)),1); // n-region

  Inverse<Region_3> *IpnC = new Inverse<Region_3>((Region_3 *)pnC);
  task.AddMediumRegion(top,IpnC); // n-region

//  task.AddGeneration(new scGenerationUniform(1e27));

/*
  scGeneration *gen = new scLambertSpectrum(Plane_3(Vector_3(0,0,-1),Vector_3(0,0,h)),spectra_dir.c_str(),
    "am1.5.in","si.nk",350,1240.8/top.Eg);
  task.AddGeneration(gen);
  */
  
  const char *fname="abs_3d.d"; // file with fdtd results
  // linear transformation between FDTD and Microvolt coordinate systems
  Basis_3 bs(Vector_3(1,0,0),Vector_3(0,1,0),Vector_3(0,0,-1));
  valtype a=.8; // period of the nanowires array in FDTD
  Vector_3 shift(a/2,a/2,-h); // shift between FDTD and Microvolt coordinate systems
  // coordinates are in 0,1,2 columns; generation is in 3 column
  // 3 coordinates are used (3D).
  // FDTD output generation is in 1/sec/m/m/um; while Microvolt input generation is in 1/sec/m3.
  // to compensate for this, we use 1e6 factor
  scGeneration *gen=new scGenerationTable(fdtd_dir.c_str(),fname,iVector_3(0,1,2),3,3,bs,shift,1e6);
  gen = new scGenerationCylinder(make_mngarg(gen),2,0,64);
  task.AddGeneration(gen);
  
  // detectors resolution
  valtype detx=.2, detz=.05;
  task.AddDetector("d",p1,p2,Vector_3(detx,detx,detz));

//  task.AddDetector("yz",p1,p2,Vector_3(0,detx,detz));
//  task.AddDetector("xz",p1,p2,Vector_3(detx,0,detz));

  Vector_3 j0(-R-dx,-R-dx,h/2), j1(R+dx,R+dx,h/2);
  task.AddBoxDetector("j",j0,j1,iVector_3(nx+2,nx+2,1),BOX_BACK_Z);
  task.SetFluxIV("j");

//  for(int i=0;i<8;i++){
//    j0[2]=j1[2]=i*h/8;
//    task.AddBoxDetector(fmt("%s%d","jh",i),j0,j1,iVector_3(nx+2,nx+2,1),BOX_BACK_Z);
//  }

  task.Calculate(.025);

  return 1;
}
