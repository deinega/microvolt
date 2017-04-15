#ifndef PHYSCONST_H
#define PHYSCONST_H

/// @file physconst.h  Useful physical constants and scaling factors for different units of measurement

#include <cmath>

// some useful physical constants (from NIST 2010 list)
namespace phys {

const double c = 299792458; // speed of light, m / s
const double q = 1.602176565e-19; // electron charge C
const double kB_eV = 8.6173324e-5; // coefficient from Kelvin to eV (1eV=11604,505K), eV/K
const double kB = q*kB_eV; // Boltzman constant, J/K
const double h_eV = 4.135667516e-15; // Planck constant, eV s (Energy = _H * _C / wavelength)
const double h = q*h_eV; // Planck constant, J s (this is Plank constant without bar)
const double m_e = 9.10938291e-31; // elecron mass, kg
const double mu0 = 4*M_PI*1e-7; // magnetic constant, H / m
const double eps0 = 1/(mu0*c*c); // C^2 / J m

// effective density of the conduction band states (see eq. (5.5) Sheng book), 1/(K^1.5*m^3)
const double NCV = 2*std::pow(2*M_PI*m_e*kB/(h*h), 1.5);

}

// Scaling factors for different units of measurement
// TODO: replace with a more reliable system like Boost.Units
namespace units {

const double m = 1;
const double cm = 0.01;
const double um = 1e-6;
const double nm = 1e-9;

}

#endif
