#ifndef SC_MEDIUMBASE_H
#define SC_MEDIUMBASE_H

/// @file sc_mediumbase.h Classes to specify medium properties
/// scAbstractMedium is a base class for description of medium properties
/// scMedium is derived class for semiconductors with constant doping
/// Other derived class are possible, for example, for extension of Microvolt for organic solar cells.

#include <utility>
#include <vector>
#include <memory>
#include <functional>
#include <cassert>

#include "vector_3.h"
#include "cpp11features.h"
#include "sc_context.h"
#include "sc_carrier_params.h"


/*
Abstract class to specify properties of arbitrary medium in space
It has cass Local which is optimized to store medium property in spefic point,
and class LocalPool which is optimized to store medium prorties in a set of points 
(typically, mesh nodes).
If medium properties do not vary in space, class LocalPool is trivial
*/
class scAbstractMedium {
protected:
#ifdef _MSC_VER
  // The compiler shouldn't allow anyone to create instances of abstract classes, but MSVC does
  // (e. g., it allows constructions like
  // base b[1] = {derived()};
  // where base is abstract, confirmed for MSVC10),
  // so we explicitly put the constructors to the protected section
  scAbstractMedium() { }
  scAbstractMedium(const scAbstractMedium&) { }
#endif

public:
  // parameters to specify different types of recombination
  // (radiative, Auger and Shockley-Read-Hall)
  // all length variables are in m
  // 
  // The structure may be irrelevant for some types of media
  // but the main code significantly depends on it, so it'll
  // be here until all medium-specific code will be moved to
  // specialized classes
  struct Recomb {
    // radiative recombination
    valtype B; // m3/s
    // Auger recombination
    CarrierParams C; // m6/s

    // SRH recombination
    // difference between trap level and intrinsic fermi level
    valtype trap; // eV
    // lifetime, see Sheng book, eq (6.24)
    CarrierParams tau; // s
    // carriers density n1 and p1, when the Fermi level coincides with the trap level, 
    // see equation (6.14) from the Sheng book
    CarrierParams ntrap; // 1/ m3

    Recomb() : B(0), C(0,0),  trap(0), tau(0,0), ntrap(0,0) { }

    valtype get_recomb(valtype n, valtype p, valtype np0) const;
  };

  // Local parameters block
  class Local {
    const scAbstractMedium* parent_;

  protected:
    typedef Local base;

    const scAbstractMedium* parent() const { return parent_; }

  public:
    Local(const scAbstractMedium* parentMedium) : parent_(parentMedium) { }
    virtual ~Local() { }

    // assign medium properties corresponding to given space position
    // by default, medium properties are constant, and this function is empty
    virtual void setPoint(const Vector_3&) { }

    // Return a default context for calculation of parameters. The context should
    // contain default values for all variables needed.
    const scPointContext& getDefaultContext() const { return parent()->getDefaultContext(); }

    virtual valtype get_eps() const = 0;
    virtual valtype get_dop() const = 0;
    virtual CarrierParams get_mu() const = 0;
    virtual CarrierParams get_D() const = 0;
    virtual valtype get_ni() const = 0;
    virtual valtype get_ifermi() const = 0;
    virtual valtype get_T() const = 0;
    virtual valtype get_Eg() const = 0;
    virtual valtype get_chi() const = 0;
    virtual CarrierParams get_m() const = 0;
    virtual Recomb getRecombParams() const = 0;

    // returns 1 if medium properties change continuosly at this point.
    // this function can be used in some dischretization schemes.
    virtual bool if_gradual(const Local& other, valtype fr_lim = 1) const {
      return parent() == other.parent();
    }

    virtual valtype get_recomb(valtype n, valtype p) const {
      return getRecombParams().get_recomb(n, p, get_np());
    }
    virtual CarrierParams get_conc0() const;
    virtual valtype get_builtin() const;
    virtual CarrierParams get_N() const;

    valtype get_np() const { valtype ni = get_ni(); return ni*ni; }
    valtype get_workfunc() const { return get_builtin() + get_ifermi(); }

    // Calculation of carrier concentrations (in 1/m3) from quasi-Fermi levels (in eV)
    virtual CarrierParams calcConcentration(CarrierParams fermi, valtype psi = 0) const;

    // Calculation of quasi-Fermi levels (in eV) from carrier concentrations (in 1/m3)
    virtual CarrierParams calcQFermi(CarrierParams concentr, valtype psi = 0) const;

    // Calculation of enchancement factor g3 = (conc/kT)/(dconc/dfermi) as a function of fermi levels
    virtual CarrierParams calcG3(CarrierParams fermi, valtype psi = 0) const {
      return CarrierParams(1.,1.);
    }

    // Context-dependent parameters
    virtual valtype get_eps(const scPointContext&) const { assert(getDefaultContext().empty()); return get_eps(); }
    virtual CarrierParams get_mu(const scPointContext&) const { assert(getDefaultContext().empty()); return get_mu(); }
    virtual CarrierParams get_D(const scPointContext&) const { assert(getDefaultContext().empty()); return get_D(); }
    virtual Recomb getRecombParams(const scPointContext&) const { assert(getDefaultContext().empty()); return getRecombParams(); }

    virtual valtype get_recomb(valtype n, valtype p, const scPointContext& ctx) const {
      return getRecombParams(ctx).get_recomb(n, p, get_np());
    }

    virtual CarrierParams get_excIniState(valtype recomb) const { return CarrierParams();}

    // TODO: revise boundary conditions so that this function won't be needed here
    virtual bool get_useReflectedCharge() const { return false; }

    virtual bool isOrganic() const {return false;}
  };

  // Pool of Local's descendants for set of points (enumerated by some index)
  class LocalPool {
  private:
    // The class is non-copyable
    LocalPool(const LocalPool&) DELETED;
    const LocalPool& operator=(const LocalPool&) DELETED;

  protected:
    LocalPool() {  }

  public:
    virtual ~LocalPool() {  }
    virtual Local& operator[](int n) = 0; // gets Local for some spefic (mesh) index
    virtual const Local& operator[](int n) const = 0;
    virtual size_t sizeInBytes() const { return sizeof(*this); }
  };

  typedef UNIQUE_OR_AUTO(LocalPool) LocalPoolPtr;
  typedef UNIQUE_OR_AUTO(Local) LocalPtr;

  virtual ~scAbstractMedium() { }

  virtual LocalPoolPtr getLocalMediumPool(int nElements) const = 0;
  virtual LocalPtr getLocalMedium() const = 0;

  // Create a pool of local medium data for n points with
  // coordinates returned by pointFunc
#if defined(HAS_LAMBDAS) // && HAS_LAMBDAS >= 11 // Not sure if any implementation will work with that
  LocalPoolPtr createLocalMediumPool(int nElements,
    const std::function<Vector_3(int)>& pointFunc) const {
#else
  template<class T>
  LocalPoolPtr createLocalMediumPool(int nElements, T pointFunc) const {
#endif
    LocalPoolPtr p = getLocalMediumPool(nElements);
    for(int i = 0; i < nElements; ++i)
      (*p)[i].setPoint(pointFunc(i));

    return p;
  }

  // Create a pool of locals for an array of points
  // specified with an iterator range
  template<class InputIterator> LocalPoolPtr createLocalMediumPool(
    InputIterator positionBegin, InputIterator positionEnd) const {
    int nElems = 0;
    for(InputIterator it = positionBegin; it != positionEnd; ++it)
      nElems++;

    LocalPoolPtr p = getLocalMediumPool(nElems);
    int i = 0;
    for(InputIterator it = positionBegin; it != positionEnd; ++it)
      (*p)[i++].setPoint(*it);

    return p;
  }

  LocalPtr createLocalMedium(const Vector_3& p = Vector_3()) const {
    LocalPtr med = getLocalMedium();
    med->setPoint(p);
    return med;
  }

  // Return a default context for calculation of parameters. The context should
  // contain default values for all variables needed.
  const scPointContext& getDefaultContext() const { return defaultContext; }

protected:
  scPointContext defaultContext;
};

// semiconductor with constant properties in space
class scMedium : public scAbstractMedium {
  friend class scMediumFactory;

  // temperature, K
  valtype T;

  // dielectric permittivity
  valtype eps;
  // concentration of dopants, 1/m3
  valtype dop;
  // bandgap and electron affinity
  valtype Eg, chi; // eV
  // mobility
  CarrierParams mu; // m2/Vs

  // either m or N shoud be defined (if m is defined, than N is expresses via m, and vice versa)

  // effective electron (hole) mass relative to electron mass in vacuum
  CarrierParams m; // 1
  // effective density of states in conduction and valence band,
  // see coef in eq.(5.9) from Sheng book
  CarrierParams N; // 1/m3 (Nc and Nv)

  // recombination parameters
  Recomb rcm;   
  
  // parameters that are calculated using already parameters

  // intrinsic fermi level, difference between ifermi and middle of the bandgap
  valtype ifermi,ifermi_bg; // eV
  // diffusion coefficient
  CarrierParams D; //m2/s
  // intrinsic carriers concentration
  valtype ni; // 1/m3
  
public:

  class LocalPool;

  // Local parameters block
  class Local : public scAbstractMedium::Local {
  protected:
    const scMedium* parent() const {
      return static_cast<const scMedium*>(base::parent());
    }

  public:
    explicit Local(const scMedium* constParams = 0)
      : base(constParams) { }

    virtual ~Local() { }

    virtual valtype get_eps() const { return parent()->eps; }
    virtual valtype get_dop() const { return parent()->dop; }
    virtual CarrierParams get_mu() const { return parent()->mu; }
    virtual CarrierParams get_D() const { return parent()->D; }
    virtual valtype get_ni() const { return parent()->ni; }
    virtual valtype get_ifermi() const { return parent()->ifermi; }
    virtual valtype get_T() const { return parent()->T; }
    virtual valtype get_Eg() const { return parent()->Eg; }
    virtual valtype get_chi() const { return parent()->chi; }
    virtual CarrierParams get_m() const { return parent()->m; }
    virtual CarrierParams get_N() const { return parent()->N; }

    //virtual CarrierParams get_excIniState(valtype recomb) const { return CarrierParams();}

    virtual Recomb getRecombParams() const { return parent()->rcm; }

    // Called frequently - faster version
    virtual valtype get_recomb(valtype n, valtype p) const {
      return Local::getRecombParams().get_recomb(n, p, Local::get_ni()*Local::get_ni());
    }

    virtual bool isOrganic() const {return false;};
  };

  class LocalPool : public scAbstractMedium::LocalPool {
    Local elem; // store only one element, since properties in space are constant

  public:
    explicit LocalPool(const scMedium* params) : elem(params) {  }
    virtual Local& operator[](int n) { return elem; } // properties in space are constant
    virtual const Local& operator[](int n) const { return elem; }

    virtual size_t sizeInBytes() const { return sizeof(*this); }
  };

  virtual ~scMedium() { }

  virtual LocalPoolPtr getLocalMediumPool(int /*nElements*/) const {
    return LocalPoolPtr(new LocalPool(this));
  }

  virtual LocalPtr getLocalMedium() const {
    return LocalPtr(new Local(this));
  }
};

#endif
