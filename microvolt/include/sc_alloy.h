#ifndef SC_ALLOY_H
#define SC_ALLOY_H

/// @file sc_alloy.h  Classes to specify properties of alloys (mixtures of two media)

#include <memory>
#include <utility>

#include "vector_3.h"
#include "sc_mediumbase.h"
#include "physconst.h"

// Base class for media with coordinate-dependent parameters
// We assume that all local values of parameters can be expressed
// through one position-dependent quantity and constant "global" parameters.
class scParametricMedium : public scAbstractMedium {
public:
  // It's better to use std::function<vec_type(const Vector_3&)> if supported
  // TODO: implement switching between virt_unary_function<> and std::function<>
  //   depending on availability of the latter
  typedef virt_unary_function<const Vector_3 &, vec_type> PositionFunction;
  typedef scAbstractMedium::Local MLocal;
  typedef SHARED_OR_AUTO(const PositionFunction) PosFuncPtr;

#ifndef HAS_SMART_PTR
  // the class is non-copyable
  scParametricMedium(const scParametricMedium&);
  const scParametricMedium& operator=(const scParametricMedium&);
#endif

private:
  // Position-dependent function. Returns characteristic value
  // (e.g. dopant fraction) of a position
  PosFuncPtr fraction;

protected:
#ifdef HAS_SMART_PTR
  scParametricMedium(const PosFuncPtr& fractionFunction)
#else
  scParametricMedium(PosFuncPtr fractionFunction)
#endif
    : fraction(fractionFunction) { }

public:
  class Local : public scAbstractMedium::Local {
    PositionFunction::result_type fraction;

  protected:
    const scParametricMedium* parent() const {
      return static_cast<const scParametricMedium*>(base::parent());
    }

    Local(const scParametricMedium* constParams = 0)
      : base(constParams) { }

    PositionFunction::result_type get_fr() const { return fraction; }

    template<class T>
    T lin_interp(T val1, T val2) const
    { return get_fr()*val1 + (1-get_fr())*val2; }

  public:
    virtual ~Local() { }

    // stores function value at the given point
    virtual void setPoint(const Vector_3& p) {
      fraction = parent()->fractionFunction(p);
    }

    // By default, the parameters are calculated from the respective parameters
    // of the media using getParam(); descendants may do it in a different way
    virtual valtype get_eps() const { return getParam(&MLocal::get_eps); }
    virtual valtype get_dop() const { return getParam(&MLocal::get_dop); }
    virtual CarrierParams get_mu() const { return getParam(&MLocal::get_mu); }
    virtual CarrierParams get_D() const { return getParam(&MLocal::get_D); }
    virtual valtype get_ni() const { return getParam(&MLocal::get_ni); }
    virtual valtype get_ifermi() const { return getParam(&MLocal::get_ifermi); }
    virtual valtype get_T() const { return getParam(&MLocal::get_T); }
    virtual valtype get_Eg() const { return getParam(&MLocal::get_Eg); }
    virtual valtype get_chi() const { return getParam(&MLocal::get_chi); }
    virtual CarrierParams get_m() const { return getParam(&MLocal::get_m); }

    virtual bool if_gradual(const scAbstractMedium::Local& other, valtype fr_lim = 1) const {
      return scAbstractMedium::Local::if_gradual(other, fr_lim) &&
        std::abs(fraction - static_cast<const Local&>(other).fraction) < fr_lim;
    }

  private:
    // Calculate parameter from the values of the respective parameters
    // of the media
    virtual valtype getParam(valtype (MLocal::*)() const) const = 0;
    virtual CarrierParams getParam(CarrierParams (MLocal::*)() const) const = 0;
  };

  template<class T> class LocalPool : public scAbstractMedium::LocalPool {
    typedef typename T::Local TLocal;

    std::vector<TLocal> elems;

  public:
    LocalPool(int nElements, const T* alloy)
      : elems(nElements, TLocal(alloy)) { }

    virtual TLocal& operator[](int n) { return elems[n]; }
    virtual const TLocal& operator[](int n) const { return elems[n]; }

    virtual size_t sizeInBytes() const {
      return sizeof(*this) + elems.size()*sizeof(TLocal);
    }
  };

  virtual ~scParametricMedium() { }

  PositionFunction::result_type fractionFunction(const Vector_3& p) const {
    // The base class, virt_unary_function, was incorrectly defined with operator() having no
    // const qualifier, though the functions are immutable in fact (at least the ones in this
    // library). In order to preserve compatibility with user code which may contain classes
    // derived from virt_unary_function we don't add const to the base class, but use
    // const_cast here instead
    return (const_cast<PositionFunction&>(*fraction))(p);
  }
};

/* Medium with variable doping.
Doping is linear combination of scMediumFactory::dop and dop2 with
weights defined by fractionFunction.
*/
class scVarDopMedium : public scParametricMedium {
public:
  typedef SHARED_OR_AUTO(const scMedium) MediumPtr; // semiconductor medium

#ifdef HAS_SMART_PTR
  scVarDopMedium(const MediumPtr& med_, valtype dop2_, const PosFuncPtr& fractionFunction)
#else
  scVarDopMedium(MediumPtr med_, valtype dop2_, PosFuncPtr fractionFunction)
#endif
    : scParametricMedium(fractionFunction), med(med_), dop2(dop2_) { }

  // Local parameters block
  class Local : public scParametricMedium::Local {
  private:
    scMedium::Local med; // keeps pointer on scMedium
    valtype dop;

    virtual valtype getParam(valtype (MLocal::*param)() const) const {
      return (med.*param)();
    }

    virtual CarrierParams getParam(CarrierParams (MLocal::*param)() const) const {
      return (med.*param)();
    }

  protected:
    const scVarDopMedium* parent() const {
      return static_cast<const scVarDopMedium*>(base::parent());
    }

  public:
    // doping is initialized with the second doping constant instead of
    // local doping since we don't know the local value of fraction yet
    explicit Local(const scVarDopMedium* params = 0)
      : scParametricMedium::Local(params),
      med(&*(params->med)),
      dop(params ? params->dop2 : 0) { }

    virtual ~Local() { }

    virtual void setPoint(const Vector_3& p) {
      scParametricMedium::Local::setPoint(p);
      med.setPoint(p);
      dop = lin_interp(med.get_dop(), parent()->dop2);
    }

    // These ones are called frequently, so we provide optimal ways to determine the parameters
    virtual CarrierParams get_D() const { return med.scMedium::Local::get_D(); }
    virtual valtype get_ni() const { return med.scMedium::Local::get_ni(); }
    virtual valtype get_ifermi() const { return med.scMedium::Local::get_ifermi(); }
    virtual valtype get_recomb(valtype n, valtype p) const { return med.scMedium::Local::get_recomb(n, p); }

    virtual valtype get_dop() const { return dop; }
    virtual Recomb getRecombParams() const { return med.getRecombParams(); }
  };

  virtual ~scVarDopMedium() { }

  virtual LocalPoolPtr getLocalMediumPool(int nElements) const {
    return LocalPoolPtr(new LocalPool<scVarDopMedium>(nElements, this));
  }

  virtual LocalPtr getLocalMedium() const {
    return LocalPtr(new Local(this));
  }

private:
  // doping is varied between med.dop and dop2
  MediumPtr med;
  valtype dop2;
};

template <class M1, class M2> class scAlloyFactory;
/* 
 Properties of alloy med1_{fr} med2_{1-fr} based on properties of its components
 med1, med2 and their fill fractions fr, 1-fr (0<=fr<=1).
 most of parameters are linearly interpolated and don't need to be stored.
 some parameters (ni, ifermi) which are calculated using more complicated expressions, 
 are stored in order to improve code efficiency
*/
template<class M1 = scMedium, class M2 = M1>
class scAlloy : public scParametricMedium {
public:
  typedef M1 Medium1;
  typedef M2 Medium2;

  typedef SHARED_OR_AUTO(const Medium1) Medium1Ptr;
  typedef SHARED_OR_AUTO(const Medium2) Medium2Ptr;

#ifdef HAS_SMART_PTR
  scAlloy(const Medium1Ptr& med1_, const Medium2Ptr& med2_, const PosFuncPtr& fraction)
#else
  scAlloy(Medium1Ptr med1_, Medium2Ptr med2_, PosFuncPtr fraction)
#endif
    : scParametricMedium(fraction), med1(med1_), med2(med2_) { }

  // Local parameters block
  class Local : public scParametricMedium::Local {
    typedef scAbstractMedium::Recomb Recomb;

    typename Medium1::Local med1;
    typename Medium2::Local med2;

    // the same meaning as for scMedium
    // stored in order not to recalculate all the time
    valtype ni, ifermi;

    // Calculate interpolated parameters of recombination
    Recomb interpRecomb() const {
      const Recomb& r1 = med1.getRecombParams();
      const Recomb& r2 = med2.getRecombParams();
      Recomb result;
      result.B = lin_interp(r1.B, r2.B);
      result.C = lin_interp(r1.C, r2.C);
      result.tau = lin_interp(r1.tau, r2.tau);
      result.ntrap = lin_interp(r1.ntrap, r2.ntrap);
      return result;
    }

    virtual valtype getParam(valtype (MLocal::*param)() const) const
    { return lin_interp((med1.*param)(), (med2.*param)()); }

    virtual CarrierParams getParam(CarrierParams (MLocal::*param)() const) const 
    { return lin_interp((med1.*param)(), (med2.*param)()); }

  public:
    explicit Local(const scAlloy* params = 0)
      : scParametricMedium::Local(params),
      med1(&*(params->med1)),
      med2(&*(params->med2)) { }

    virtual ~Local() { }

    virtual void setPoint(const Vector_3& p) {
      scParametricMedium::Local::setPoint(p);
      med1.setPoint(p);
      med2.setPoint(p);

      CarrierParams N = get_N();
      valtype Eg = get_Eg(), chi = get_chi();
      
      valtype np = N[0]*N[1]*std::exp(-Eg/(phys::kB_eV*get_T()));
      ni = std::sqrt(np);

      valtype ifermi_bg = phys::kB_eV*get_T()*std::log(N[1]/N[0])/2; // Neamen, (4.25)
      ifermi = - (chi + Eg/2 - ifermi_bg);
    }

    virtual valtype get_ni() const { return ni; }
    virtual valtype get_ifermi() const { return ifermi; }
    virtual Recomb getRecombParams() const { return interpRecomb(); }
  };

  virtual ~scAlloy() { }

  virtual LocalPoolPtr getLocalMediumPool(int nElements) const {
    return LocalPoolPtr(new LocalPool<scAlloy>(nElements, this));
  }

  virtual LocalPtr getLocalMedium() const {
    return LocalPtr(new Local(this));
  }

private:
  Medium1Ptr med1;
  Medium2Ptr med2;
};

#endif
