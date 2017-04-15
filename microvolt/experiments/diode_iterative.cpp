#include "sc_ui.h"

/*
This is 1d diode test code for iterative method to improve convergence
*/

int main(int argc,char **argv){

  string dir,spectra_dir;
  dir="../../../../runs/diode"; // output directory
  spectra_dir="../../spectra";

  if(scInit(argc,argv,dir)<0)
    return -1;

  /* Length is d
  left medium at (0, left_length) and right medium at (left_length, length) 
  */

  valtype length=1; // diode width, micron
  valtype left_length=.5; // pn-junction position, if pn=-1 than pn=d/2

  scMediumFactory left_med, right_med;

  valtype L=1e5; // diffusion length, micron
  valtype sr=1e5; // surface recombination, cm/s
  valtype dr=0.01; // mesh size
  int wd=2; // working dimension (z as a default)

  left_med=getscGaAs(recombRad);
  right_med=left_med;
  left_med.L[0]=right_med.L[0]=L;

  Vector_3 sz(0);
  sz[wd]=length;
  scExperiment task(0,sz);

  task.DumpMeshes(); // meshes will be dump their values

  task.SetResolution(dr);

  Vector_3 n;
  n[wd]=1;

  task.AddContact(GetHalfSpace(-n,VEC_ZERO*n),sr,sr);
  // voltage will be applied to this boundary 
  valtype wf=VEC_INFTY;
  task.AddVoltageContact(GetHalfSpace(n,(length-VEC_ZERO*dr)*n),sr,sr,wf);

  // solar spectrum absorption
   scGeneration *genAM = new scLambertSpectrum(Plane_3(n,Vector_3()),spectra_dir.c_str(),
     "am1.5.in","gaas.nk",350,1240.8/left_med.Eg);
  task.AddGeneration(genAM);

  iVector_3 dsz(1);
  task.AddDetector("d",0,sz,iVector_3(INT_INFTY));
  task.AddFluxPlane("j",wd,length);
  task.SetFluxIV("j"); // I-V calculation

//  const char *fname="d_17n_concentration_1000sun.d";  
//  Basis_3 bs(Vector_3(0,0,1),Vector_3(0,1,0),Vector_3(1,0,0));
//  dir+="\\..";
//  scGeneration *gen=new scGenerationTable(dir.c_str(),fname,0,10,1,bs,0,1e6);
//  task.AddGeneration(gen);

//  const char *fname_np="d_17n_concentration_1000sun.d";
//  scGeneration *fer_n=new scGenerationTable(dir.c_str(),fname_np,0,2,1,bs,0,1); 
//  scGeneration *fer_p=new scGenerationTable(dir.c_str(),fname_np,0,3,1,bs,0,1);
//  task.AddFermiGuess(fer_n,fer_p);

  task.SetConcentrationUnit(left_med);

  valtype dop=1e20; // doping concentration, cm-3

  valtype dop_it=1e14; // start dopping concentration, cm-3

  while(1){ // iterative cycle

    left_med.dop=dop_it;
    right_med.dop=-dop_it;
    task.AddMediumRegion(left_med,GetHalfSpace(-n,left_length*n));
    task.AddMediumRegion(right_med,GetHalfSpace(n,left_length*n));
    task.BuildMediumRegions(rpPsi); // submitting updated medium regions and keeping fermi levels from previous solution (but not potential psi)
    if(dop_it>=dop)
      break;
    if(task.Calculate(.0)<0)
      return -1;
    dop_it*=100;
    if(dop_it>dop)
      dop_it=dop;
    task.ClearMediumRegions();
  };

  task.Calculate(.0); // final solution

  return 1;
}
