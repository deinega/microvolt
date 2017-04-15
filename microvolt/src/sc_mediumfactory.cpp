
#include <cassert>
#include <cmath>

#include "sc_mediumfactory.h"
#include "sc_mediumbase.h"

using namespace std;
using namespace phys;

scMediumFactory::InstancePtr scMediumFactory::createSpecific() const {
  SHARED_OR_AUTO(scMedium) med(MAKE_SHARED_OR_AUTO(scMedium, ()));

  valtype coef2 = pow(lengthUnit, 2);
  valtype coef3 = pow(lengthUnit, 3);

  med->T = T;
  med->eps = eps;
  med->dop = dop/coef3;
  med->Eg = Eg;
  med->chi = chi;

  for(int i = 0; i < 2; i++) {
    med->N[i] = N[i] != 0 ? N[i]/coef3 : NCV*pow(m[i]*T, 1.5);
    med->m[i] = N[i] != 0 ? pow(N[i]/coef3/NCV, 1/1.5)/T : m[i];

    med->mu[i] = mu[i]*coef2;
    // Einstein relation
    med->D[i] = med->mu[i]*kB_eV*T;
  }

  // Nc/Nv
  valtype Ncv = med->N[0]/med->N[1]; // 1
  // this is necessary for heterostructures
  med->ifermi_bg = kB_eV*T*log(1/Ncv)/2; // Neamen, (4.25)
  med->ifermi = - (chi + Eg/2 - med->ifermi_bg);

  if (ni == 0) { // ni is not defined by user
    valtype np = med->N[0]*med->N[1] * exp(-med->Eg/(kB_eV*T));
    med->ni = sqrt(np);
  }
  else {
    med->ni = ni/coef3;
    med->N[0] = ni*exp(med->Eg/(2*kB_eV*T))*sqrt(Ncv);
    med->N[1] = med->N[0]/Ncv;
    for(int i = 0; i < 2; i++)
      med->m[i] = pow(med->N[i]/NCV, 1/1.5)/T;
  }

  rcmInit(&*med);

#ifdef HAS_SMART_PTR
  return med;
#else
  // under GCC, std::auto_ptr needs an explicit conversion here
  return InstancePtr(med);
#endif
}

scAbstractMediumPtr scMediumFactory::create() const {
#ifdef HAS_SMART_PTR
  return createSpecific();
#else
  // auto_ptr needs an additional conversion to auto_ptr_ref, so we do an explicit cast
  // in order to avoid compiler warning
  return scAbstractMediumPtr(createSpecific());
#endif
}

void scMediumFactory::rcmInit(scMedium* med) const {
  scMedium::Recomb& mrcm = med->rcm;

  valtype coef3 = pow(lengthUnit, 3);
  CarrierParams L_sc;

  mrcm.B = rcm.B * coef3;
  mrcm.trap = rcm.trap;

  for(int i = 0; i < 2; i++) {
    mrcm.C[i] = rcm.C[i] * coef3 * coef3;
    mrcm.tau[i] = rcm.tau[i];
    L_sc[i] = L[i] * units::um;

    // SRH recombination
    if(mrcm.tau[i]) // tau specified, find diffusion length
      L_sc[i] = sqrt(med->D[i]*mrcm.tau[i]);
    else if(L_sc[i]) // diffusion length specified, find tau
      mrcm.tau[i] = L_sc[i]*L_sc[i]/med->D[i];
    else { // tau and diffusion length are not specified, use tau of other carrier
      mrcm.tau[i] = rcm.tau[1-i] ?
        rcm.tau[1-i] : L_sc[1-i]*L_sc[1-i]/med->D[1-i];

      L_sc[i] = sqrt(med->D[i]*mrcm.tau[i]);
    }
  }

  med->rcm.ntrap[0] = med->ni*exp(rcm.trap/(kB_eV*T));
  med->rcm.ntrap[1] = med->ni*exp(-rcm.trap/(kB_eV*T));
//  med->rcm.ntrap[0]=med->N[0]*exp(Eg/(2*kB_eV*T)); // in Zvyagin notation
 // med->rcm.ntrap[1]=med->N[1]*exp(-Eg/(2*kB_eV*T));
}

scMediumFactory getscSi(int recomb){
  scMediumFactory med;
  med.eps=11.9;
  med.Eg=1.124;
  med.chi=4.05;

  med.mu[0]=1107;
  med.mu[1]=424.6;
  med.m[0]=1.08;
  med.m[1]=0.56;

  med.N[0]=2.86e19;
  med.N[1]=2.66e19;

  if(recomb&recombRad)
    med.rcm.B=9.5e-15;
  if(recomb&recombAuger){
    med.rcm.C[0]=2.2e-31;
    med.rcm.C[1]=9.9e-32;
  }
  // SRH recombination, default value
  if(recomb&recombSRH)
    med.rcm.tau[0]=1e-7;

  return med;
};

scMediumFactory getscaSi(int recomb){
  scMediumFactory med;
  med.eps=11.7;
  med.Eg=1.6;
  med.chi=4;//

  med.mu[0]=1;
  med.mu[1]=0.5;
  med.m[0]=1.08;//
  med.m[1]=0.56;//

  med.N[0]=2.86e19;//
  med.N[1]=2.66e19;//

  if(recomb&recombRad)
    med.rcm.B=9.5e-15;//
  if(recomb&recombAuger){
    med.rcm.C[0]=2.2e-31;//
    med.rcm.C[1]=9.9e-32;//
  }
  // SRH recombination, default value
  if(recomb&recombSRH)
    med.rcm.tau[0]=1e-7;

  return med;
}

scMediumFactory getscGe(int recomb){
  scMediumFactory med;
  med.eps=16.2;
  med.Eg=0.66;
  med.chi=4.13;

  med.mu[0]=3900;
//  med.mu[1]=1900;
  med.mu[1]=2500;
//  med.mu[0]=1107;
//  med.mu[1]=424.6;
  med.m[0]=0.55;
  med.m[1]=0.37;
  med.N[0]=1.04e19;
  med.N[1]=.6e19;

  if(recomb&recombAuger){
    med.rcm.C[0]=1e-30;
    med.rcm.C[1]=1e-30;
  }
  // SRH recombination, default value
  if(recomb&recombSRH)
    med.rcm.tau[0]=1e-7;

  return med;
}

scMediumFactory getscCdTe(){
  scMediumFactory med;
  med.eps=9.4;
  //med.Eg=1.45; // Sites
  // med.Eg=1.49; // Spath
  med.Eg=1.5; // Zvyagin
  //med.chi=4.4; // Zvyagin
  med.chi=4.28; // Sites
  //med.chi=4.5; // Spath

 // med.mu[0]=320;
 // med.mu[1]=40;
    med.mu[0]=1000;
  med.mu[1]=1000;
  med.m[0]=0.11;
  med.m[1]=0.6;
// med.N[0]=2.86e19;
// med.N[1]=2.66e19;
//  med.B=9.5e-15;
//  med.C[0]=2.2e-31;
//  med.C[1]=9.9e-32;
  // SRH recombination, default value
 // med.ni=4e5;
  med.rcm.tau[0]=2e-9;
  med.rcm.tau[1]=2e-9;
  return med;
};

scMediumFactory getscCdS(){
  scMediumFactory med;
  med.eps=10;
  med.Eg=2.42; // Sites
 // med.Eg=2.4; // Zvyagin
  //med.chi=4.5; // Zvyagin
  med.chi=4.4; // Sites

  med.mu[0]=100;
  med.mu[1]=25;
  //med.ni=0.1; // sm-3
  med.m[0]=0.2;
  med.m[1]=0.8;
// med.N[0]=2.86e19;
// med.N[1]=2.66e19;

//  med.B=9.5e-15;
//  med.C[0]=2.2e-31;
//  med.C[1]=9.9e-32;
  // SRH recombination, default value
  med.rcm.tau[0]=1e-9;
  med.rcm.tau[1]=1e-9;
  return med;
};

scMediumFactory getscGaAs(int recomb){
  scMediumFactory med;
 //dielectric constant
  med.eps=12.95;
  // electron bandgap
  med.Eg=1.425; 
  //affinity
  med.chi=4.07; // 
  //mobility
  med.mu[0]=8500; //electron
  med.mu[1]=400;  //hole 
 
  // effective electron (hole) mass
  med.m[0]=0.067;
  med.m[1]=0.45; 

  //effective density of states in conduction and valence band
  med.N[0]=4.35e17;
  med.N[1]=7.57e18;

  //The radiative recombination coefficient of GaAs from laser delay
  //measurements and effective nonradiative carrier lifetimes
  //G. W. Hooft Citation: Appl. Phys. Lett. 39, 389 (1981)
  if(recomb&recombRad)
    med.rcm.B=1.3e-10;
  //if(recomb&recombAuger){
   // med.rcm.C[0]=2.2e-31;
   // med.rcm.C[1]=9.9e-32;
  //}
  // SRH recombination, default value
 //if(recomb&recombSRH)
 // med.rcm.tau[0]=1e-9;
 // med.rcm.tau[1]=1e-9;

  return med;
};

scMediumFactory getscZnO(){
  scMediumFactory med;
  med.eps=9;
  med.Eg=3.3; 
  med.chi=4.6;

  med.mu[0]=100;
  med.mu[1]=25;
  med.m[0]=0.2;
  med.m[1]=1.2;

  // SRH recombination, default value
  med.rcm.tau[0]=1e-9;
  med.rcm.tau[1]=1e-9;
  return med;
};

scMediumFactory getscCIGS(){
  scMediumFactory med;
  med.eps=13.6;
  med.Eg=1.15;
  med.chi=4.6;

  med.mu[0]=100;  
  med.mu[1]=12.5;  
  med.m[0]=0.09;  
  med.m[1]=0.72;

  // SRH recombination, default values
  med.rcm.tau[0]=5e-9;
  med.rcm.tau[1]=5e-9;
  
  return med;
};
