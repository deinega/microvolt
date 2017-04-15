
#include <cmath>

#include "sc_mediumbase.h"

#include "physconst.h"

using namespace std;

CarrierParams scAbstractMedium::Local::get_conc0() const {
  valtype dop = get_dop();
  int major = dop < 0; // if doping is positive then electrons are major carriers
  valtype d2 = fabs(dop)/2;
  valtype np = get_np();
  CarrierParams conc;
  conc[major] = d2 + sqrt(d2*d2 + np); // Neaman, Eq. 4.60
  conc[1-major] = np/conc[major];
//  valtype test=-d2+sqrt(d2*d2+np); // another way for conc[1-major] calculating (this way is less accurate)
  return conc;
}

valtype scAbstractMedium::Local::get_builtin() const {
  return phys::kB_eV*get_T()*log(get_conc0()[0]/get_ni());
}

CarrierParams scAbstractMedium::Local::get_N() const {
  return phys::NCV*pow(get_m()*get_T(), 1.5);
}

// Calculation of carrier concentrations (in 1/m3) from quasi-Fermi levels (in eV)
CarrierParams scAbstractMedium::Local::calcConcentration(CarrierParams fermi, valtype psi) const {
  return get_ni()*altexp((fermi - (psi - get_ifermi()))/(phys::kB_eV*get_T()));
}

// Calculation of quasi-Fermi levels (in eV) from carrier concentrations (in 1/m3)
CarrierParams scAbstractMedium::Local::calcQFermi(CarrierParams concentr, valtype psi) const {
  return phys::kB_eV*get_T()*altlog(concentr/get_ni()) + (psi - get_ifermi());
}

valtype scMedium::Recomb::get_recomb(valtype n, valtype p, valtype np) const {
  valtype difnp = n*p - np; // n*p-n0*p0

  valtype R_rad = B*difnp; // radiative recombination

  valtype R_Auger = (C[0]*n + C[1]*p)*difnp; // Auger recombination

  valtype R_SRH = 0;
  if (tau[0] > 0 || tau[1] > 0) {
    valtype denom_SRH = tau[1]*(n + ntrap[0]) + tau[0]*(p + ntrap[1]);
    R_SRH = difnp/denom_SRH; // SRH recombination
  }

  valtype R = R_rad + R_Auger + R_SRH;

  return R;
}

