
#include "sc_alloyfactory.h"
#include "sc_alloy.h"

scVarDopFactory::InstancePtr scVarDopFactory::createSpecific() const {
  valtype coef3 = lengthUnit*lengthUnit*lengthUnit;
  scVarDopMedium::MediumPtr med(scMediumFactory::createSpecific());
  return MAKE_SHARED_OR_AUTO(InstancePtr::element_type, (med, dop2/coef3, fraction));
}

scAbstractMediumPtr scVarDopFactory::create() const {
#ifdef HAS_SMART_PTR
  return createSpecific();
#else
  // auto_ptr needs an additional conversion to auto_ptr_ref, so we do an explicit cast
  // in order to avoid compiler warning
  return scAbstractMediumPtr(createSpecific());
#endif
}

template<class Medium1Factory, class Medium2Factory>
typename scAlloyFactory<Medium1Factory, Medium2Factory>::InstancePtr
    scAlloyFactory<Medium1Factory, Medium2Factory>::createSpecific() const {
  return MAKE_SHARED_OR_AUTO(typename InstancePtr::element_type,
    (med1.createSpecific(), med2.createSpecific(), fraction));
}

template<class Medium1Factory, class Medium2Factory>
scAbstractMediumPtr scAlloyFactory<Medium1Factory, Medium2Factory>::create() const {
#ifdef HAS_SMART_PTR
  return createSpecific();
#else
  // auto_ptr needs an additional conversion to auto_ptr_ref, so we do an explicit cast
  // in order to avoid compiler warning
  return scAbstractMediumPtr(createSpecific());
#endif
}

template class scAlloyFactory<>;
