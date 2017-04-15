#include "sc_ui.h"

/*
Planar solar cell with textured front surface

1.1,0 100,-1 .04
*/

int main(int argc,char **argv){

  string dir,spectra_dir,fdtd_dir;

  dir="../../../../runs/diode"; // output directory
  spectra_dir="../../spectra";
  fdtd_dir="../../../../runs"; // directory with fdtd generation file

  if(scInit(argc,argv,dir)<0)
    return -1;

 // example from web-page
  valtype length=1; // film thickness
  valtype a=.6; // period of square lattice
  valtype R=.25; // cones radius
  valtype h=.8; // cones height
  valtype pn=.5; // pn-junction position
  valtype dz=.02; // mesh step

/*  valtype length=1.6; // film thickness
  valtype a=.85; // period of square lattice
  valtype R=.5; // cones radius
  valtype h=1*length; // cones height
  valtype pn=1; // pn-junction position
  valtype dz=.04; // mesh step
  */
  valtype shift=0*a/2; // tip of the cone is shifted from the center of the square cell 
  valtype top=.0; // top contact is extending inside the cone until this value
  valtype dpn=0; // if dpn>0, curved junction will be used (like in multijunction geometry)
  valtype bsf=0; // thickness of the back surface field region
  valtype rcyl=0; // if nonzero, bottom contact is a dot of rcyl radius 
  valtype left_dop=1e18; // doping concentration of the left medium
  valtype right_dop=-1e18; // doping concentration of the right medium

  scMediumFactory med_left=getscSi(); // medium
  scMediumFactory med_right=getscSi(); // medium
  valtype L=100; // diffusion length
  valtype sr=VEC_INFTY; // contact surface recombination
//  valtype sr=1e5; // surface recombination
  valtype sri=0; // insulator surface recombination
  int rad=3; // 3 - AM1.5, 4 - FDTD results

  int nun=0; // for nonuniform mesh

  int ind=1;
  int ac=theConfig->argv.size();
  if(ac>ind)extract_list(theConfig->argv[ind++],&pn,&top,&dpn,&bsf,&rcyl);
  if(ac>ind)extract_list(theConfig->argv[ind++],&L,&sr,&sri);
  if(ac>ind)extract_list(theConfig->argv[ind++],&dz,&nun);

  if(sr<=0)
    sr=VEC_INFTY;

  string tname=fmt("%g",L);

  med_left.L[0]=L;
  med_right.L[0]=L;
  med_left.dop=left_dop;
  med_right.dop=right_dop;

  Vector_3 sz(a,a,length);
//  Vector_3 sz(0*a,a,length);

  scExperiment task(0,sz);
  task.SetResolution(dz);
  task.SetMinResidue(1e5);

  task.set_name(tname);
  task.DumpMeshes();

  if(nun>0){
    // mesh generation
    vector<Vector_3> inc;
    inc.push_back(Vector_3(pn-dz,pn+dz,dz/nun)); // more resolution at pn-junction
//    inc.push_back(Vector_3(0,3*dz,dz/10)); // more resolution at contacts
//    inc.push_back(Vector_3(d-3*dz,d,dz/10)); // more resolution at contacts
    vector<valtype> mesh;
    if(generate_1d_mesh(Vector_3(0,length,dz),inc,mesh)<0)
      return -1;
    task.SetMeshSteps(2,mesh);
  }

//  task.SetMethod("eh_p","10,1");

//  task.SetDepletionDecreasing(.05);

  task.AddMediumRegion(med_left,GetHalfSpace(-Vector_3(0,0,1),pn));
  task.AddMediumRegion(med_right,GetHalfSpace(Vector_3(0,0,1),pn));
//  task.AddMediumRegion(med,right_dop,new Region_3(),0);
  if(bsf){ // surface field
    med_left.dop=5*left_dop;
    med_right.dop=5*right_dop;
    task.AddMediumRegion(med_left,GetHalfSpace(Vector_3(0,0,1),bsf),2);
    task.AddMediumRegion(med_right,GetHalfSpace(Vector_3(0,0,1),length-bsf),2);
  }

  task.SetPeriodicBoundaries(0x1|0x2);

  vector<Vector_3> pos; // centers of central and adjacent square cells
  pos.push_back(Vector_3(sz[0]/2,sz[1]/2,0));
  for(size_t i=0;i<2;i++){
    if(sz[i]>0){
      Vector_3 im=pos[0];
      im[i]+=sz[i];
      pos.push_back(im);
      im=pos[0];
      im[i]-=sz[i];
      pos.push_back(im);
    }
  }

  if(R>0){
    for(size_t i=0;i<pos.size();i++){

      Region_3 *cone=GetPyramida(Vector_3(pos[i][0],pos[i][1],0),Vector_3(pos[i][0],pos[i][1]+shift,h),
        Vector_3(1,0,0),Vector_3(0,1,0),GetRegularPolygon(0,R,100));

      task.AddReflectiveBoundary(cone,sri,sri,1);

      if(top){ // top contact is extending inside the hole
        RegionIntersection<3> *cone_top = new RegionIntersection<3>();
        cone_top->AddRegion(cone);
        cone_top->AddRegion(GetHalfSpace(Vector_3(0,0,-1),top));
        task.AddContact(cone_top,sr,sr,VEC_INFTY,2);
      }

      if(rcyl){ // dot bottom contact
        task.AddReflectiveBoundary(GetHalfSpace(Vector_3(0,0,1),length-VEC_ZERO),sri,sri);
        Region_3 *cyl = GetCylinder(Vector_3(pos[i][0],pos[i][1],length-VEC_ZERO),Vector_3(0,0,1),rcyl,VEC_INFTY);
        task.AddVoltageContact(cyl,sr,sr,VEC_INFTY,1);
      }

      if(dpn){ // curved junction
        valtype dr=.1;
        valtype dh=dr/R*h;
        Polyhedron_3 *curved=GetPyramida(Vector_3(pos[i][0],pos[i][1],0),Vector_3(pos[i][0],pos[i][1]+shift,h+dh),
          Vector_3(1,0,0),Vector_3(0,1,0),GetRegularPolygon(0,R+dr,100));

        task.AddMediumRegion(med_left,
          GetConfinedPolyhedron(make_mngarg(curved),make_mngarg(GetHalfSpace(Vector_3(0,0,-1),pn+dpn))),1);
      }

    }
  }

  task.AddContact(GetHalfSpace(-Vector_3(0,0,1),VEC_ZERO),sr,sr); // top contact

  if(!rcyl) // planar bottom contact
    task.AddVoltageContact(GetHalfSpace(Vector_3(0,0,1),length-VEC_ZERO),sr,sr);

  if(rad==3){
    scGeneration *genAM = new scLambertSpectrum(Plane_3(Vector_3(0,0,1),0),spectra_dir.c_str(),
      "am1.5.in","si.nk",350,1240.8/med_left.Eg);
    task.AddGeneration(genAM);
  }
  else if(rad==4){
    const char *fname="abs_texture.d";
    Basis_3 bs(Vector_3(1,0,0),Vector_3(0,1,0),Vector_3(0,0,1));
    Vector_3 shift;
    scGeneration *gen=new scGenerationTable(fdtd_dir.c_str(),fname,iVector_3(0,1,2),3,3,bs,shift,1e6);
    gen = new scGenerationAverage(make_mngarg(gen),.2,10);
    task.AddGeneration(gen);
  }

//  task.AddDetector("d",Vector_3(sz[0]/2,0,0),Vector_3(sz[0]/2,sz[1],sz[2]),INT_INFTY);
  task.AddDetector("d",Vector_3(sz[0]/2,0,0),Vector_3(sz[0]/2,sz[1],sz[2]),iVector_3(1,20,20));

  task.AddDetector("d3d",Vector_3(0,0,0),Vector_3(sz[0],sz[1],sz[2]),iVector_3(6,6,10));

  if(sz[0]){
//    task.AddDetector("dl",Vector_3(dz,0,0),Vector_3(dz,sz[1],sz[2]),INT_INFTY);
    task.AddDetector("dl",Vector_3(dz,0,0),Vector_3(dz,sz[1],sz[2]),iVector_3(1,20,20));
    task.AddDetector("dh",Vector_3(sz[0]/4,0,0),Vector_3(sz[0]/4,sz[1],sz[2]),iVector_3(1,20,20));
  }

  int jn=4;
  for(int i=0;i<jn;i++){
    valtype jh=(i+.5)*length/jn;
    task.AddFluxPlane(fmt("j%d",i),2,jh);
  }
  task.AddFluxPlane("j",2,length/2);

  task.SetFluxIV("j"); // I-V calculation

  task.Calculate(.02);

  return 1;
}
