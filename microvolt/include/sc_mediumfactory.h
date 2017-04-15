#ifndef SC_MEDIUM_FACTORY_H
#define SC_MEDIUM_FACTORY_H

/// @file sc_medium_factory.h \brief Medium factory classes for specification of parameters

#include <memory>

#include "cpp11features.h"
#include "physconst.h"
#include "sc_carrier_params.h"

#ifndef HAS_SMART_PTR
// In order to delete std::auto_ptr<scMedium>, we need to know how to delete scMedium
// The same is for std::unique_ptr<scMedium>, but not std::shared_ptr<scMedium> -
// the latter contains a pointer to the deleter
#include "sc_mediumbase.h"
#endif

class scAbstractMedium;
typedef SHARED_OR_AUTO(const scAbstractMedium) scAbstractMediumPtr;


class scMedium;

class scMediumFactory {
  void rcmInit(scMedium*) const;

public:

  typedef SHARED_OR_AUTO(const scMedium) InstancePtr;

  struct scRecomb {
    // coefficient before the recombination 
    valtype rate;
    // radiative recombination
    valtype B; // cm3/s
    // Auger recombination
    CarrierParams C; // cm6/s

    // SRH recombination
    // difference between trap level and intrinsic fermi level
    valtype trap; // eV
    // lifetime, see Sheng book, eq (6.24)
    CarrierParams tau; // s

    // Nonzero depl_decr vector specifies an area where recombination rate is chosen to be artificially low.
    // this option is used to simulate high quality medium (low recombination) in depletion region, as in paper
    // B. Kayes et al.,'Comparison of the device physics principles of planar and radial p-n junction nanorod solar cells'
    Vector_3 depl_decr;

    scRecomb() : B(0), C(0,0), trap(0), tau(0,0), rate(1.) { }
  };


  // parameters that should be specified

  // Default length units (see namespace units in physconst.h)
  valtype lengthUnit;

  // temperature, K
  valtype T;

  // dielectric permittivity
  valtype eps;
  // concentration of dopants, 1/cm3
  valtype dop;
  // bandgap and electron affinity
  valtype Eg, chi; // eV
  // mobility
  CarrierParams mu; // cm2/Vs

  // either m or N shoud be defined (if m is defined, than N is expresses via m, and vice versa)

  // effective electron (hole) mass relative to electron mass in vacuum
  CarrierParams m; // 1
  // effective density of states in conduction and valence band,
  // see coef in eq.(5.9) from Sheng book
  CarrierParams N; // 1/cm3 (Nc and Nv)
  // intrinsic carriers concentration
  valtype ni; // 1/cm3

  // diffusion length, sqrt(D*tau), it is in micrometers
  // unlike the other length quantities (lengthUnit is ignored)
  CarrierParams L; // um

  // recombination parameters
  scRecomb rcm;

  
  scMediumFactory()
    : lengthUnit(units::cm), T(300), eps(1), dop(0), Eg(0), chi(0),
    mu(0,0), m(0,0), N(0,0), ni(0), L(0,0) { }

  // We need two functions which do the same, but return different types of pointers
  // in order to avoid inclusion of heavy headers where media types are defined

  // Creates an instance of scMedium and returns a pointer to it as an abstract medium
  scAbstractMediumPtr create() const;

  // Creates an instance of scMedium and returns a pointer to it as an specific medium
  InstancePtr createSpecific() const;
};


// functions to construct scMedium for some semiconductors

// recombination types
enum{
  recombRad=0x1, // radiative
  recombAuger=0x2, // Auger
  recombSRH=0x4 // Shockley-Reed-Hall, SRH
};

// cSi parameters taken from PC1D and AFORS-HET
// recomb is recombination bit flag that turns on / off different recombination types
scMediumFactory getscSi(int recomb=0);

scMediumFactory getscaSi(int recomb=0);

scMediumFactory getscGe(int recomb=0);

scMediumFactory getscCdTe();

scMediumFactory getscCdS();

// GaAs parameters from http://www.ioffe.ru/SVA/NSM/Semicond/GaAs
scMediumFactory getscGaAs(int recomb=0);

scMediumFactory getscZnO();

scMediumFactory getscCIGS();

#endif
