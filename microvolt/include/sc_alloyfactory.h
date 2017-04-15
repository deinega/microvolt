#ifndef SC_ALLOY_FACTORY_H
#define SC_ALLOY_FACTORY_H

/// @file sc_alloyfactory.h  Alloy factory classes for specification of parameters

#include <memory>

#include "sc_mediumfactory.h"

#ifndef HAS_SMART_PTR
// In order to delete std::auto_ptr<scMedium>, we need to know how to delete scMedium
// The same is for std::unique_ptr<scMedium>, but not std::shared_ptr<scMedium> -
// the latter contains a pointer to the deleter
#include "sc_alloy.h"
#endif

class scParametricMediumFactory {
public:
#ifndef HAS_SMART_PTR
  // the class is non-copyable due to presence of std::auto_ptr
  scParametricMediumFactory(const scParametricMediumFactory&);
  const scParametricMediumFactory& operator=(const scParametricMediumFactory&);
#endif

  typedef virt_unary_function<const Vector_3 &, vec_type> PositionFunction;
  typedef SHARED_OR_AUTO(const PositionFunction) PosFuncPtr;

#ifdef HAS_SMART_PTR
  PosFuncPtr fraction;
#else
  // Normally a factory don't modify itself on object creation
  // but we can make a duplicate neither of the fraction pointer
  // (since it's std::auto_ptr), nor of the function object itself
  // (it has children)
  mutable PosFuncPtr fraction;
#endif

protected:
  scParametricMediumFactory() { }
  
#ifdef HAS_SMART_PTR
  explicit scParametricMediumFactory(const PosFuncPtr& fractionFunction)
#else
  explicit scParametricMediumFactory(PosFuncPtr fractionFunction)
#endif
    : fraction(fractionFunction) { }
};


class scVarDopMedium;

/* Medium with variable doping.
Doping is linear combination of scMediumFactory::dop and dop2 with
weights defined by fractionFunction.
*/
class scVarDopFactory : public scMediumFactory, public scParametricMediumFactory {
public:
  typedef SHARED_OR_AUTO(const scVarDopMedium) InstancePtr;

  valtype dop2;

  scVarDopFactory() { }
  scVarDopFactory(const scMediumFactory& med_, valtype dop2_,
#ifdef HAS_SMART_PTR
    const PosFuncPtr& fractionFunction)
#else
    PosFuncPtr fractionFunction)
#endif
    : scMediumFactory(med_), scParametricMediumFactory(fractionFunction), dop2(dop2_) { }

  // We need two functions which do the same, but return different types of pointers
  // in order to avoid inclusion of heavy headers where media types are defined

  // Creates an instance of scMedium and returns a pointer to it as an abstract medium
  scAbstractMediumPtr create() const;

  // Creates an instance of scMedium and returns a pointer to it as an specific medium
  InstancePtr createSpecific() const;
};


template<class Medium1, class Medium2>
class scAlloy;

// Mixture of two media
// Allows to obtain mixture of more than one media, 
// if one of the template arguments is scAlloyFactory
template<class Medium1Factory = scMediumFactory, class Medium2Factory = Medium1Factory>
class scAlloyFactory : public scParametricMediumFactory {
public:
  typedef typename Medium1Factory::InstancePtr::element_type Medium1;
  typedef typename Medium2Factory::InstancePtr::element_type Medium2;

  typedef scAlloy<Medium1, Medium2> Alloy;

  typedef SHARED_OR_AUTO(const Alloy) InstancePtr;

  Medium1Factory med1;
  Medium2Factory med2;

  scAlloyFactory() { }

  scAlloyFactory(const Medium1Factory& med1_, const Medium2Factory& med2_,
#ifdef HAS_SMART_PTR
    const PosFuncPtr& fractionFunction)
#else
    PosFuncPtr fractionFunction)
#endif
    : scParametricMediumFactory(fractionFunction), med1(med1_), med2(med2_) { }

  // We need two functions which do the same, but return different types of pointers
  // in order to avoid inclusion of heavy headers where media types are defined

  // Creates an instance of scMedium and returns a pointer to it as an abstract medium
  scAbstractMediumPtr create() const;

  // Creates an instance of scMedium and returns a pointer to it as an specific medium
  InstancePtr createSpecific() const;
};

#endif
