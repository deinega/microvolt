#ifndef SC_MEDIUM_HQS_H
#define SC_MEDIUM_HQS_H

/// @file sc_medium_hqs.h  Specification of medium with low recombination level near the junction 
/// (like in paper B. Kayes et al.,'Comparison of the device physics principles of planar and radial p-n junction nanorod solar cells')

#include "sc_mediumbase.h"

class scMediumWithHQSurface : public scMedium {
  const scRegionMediumMap* rm_map;

  // Nonzero depl_decr vector specifies an area where recombination rate is chosen to be artificially low.
  // this option is used to simulate high quality medium (low recombination) in depletion region, as in paper
  // B. Kayes et al.,'Comparison of the device physics principles of planar and radial p-n junction nanorod solar cells'
  Vector_3 depl_decr;
  
public:
  scMediumWithHQSurface(const scMedium& baseMedium,
      const scRegionMediumMap* rmMap, const Vector_3& deplDecr)
    : scMedium(baseMedium), rm_map(rmMap), depl_decr(deplDecr) { }

  class Local : public scMedium::Local {
    bool low_recomb; // if recombination is low (position close to the junction)

  protected:
    const scMediumWithHQSurface* parent() const {
      return static_cast<const scMediumWithHQSurface*>(scMedium::Local::parent());
    }

  public:
    explicit Local(const scMediumWithHQSurface* constParams = 0)
      : scMedium::Local(constParams), low_recomb(false) { }

    virtual ~Local() { }

    virtual void setPoint(const Vector_3& p) {
      low_recomb = false;
      const scMediumWithHQSurface* par = parent();

      if(par->depl_decr == 0)
        return;

      scRegionMediumMap::MediumPtr med = par->rm_map->at(p);
      if (!med)
        return;

      for(int di = 2; di >= 0; di--) { // check if node is close to the junction
        if(par->depl_decr[di] == 0)
          continue;
        Vector_3 shift;
        for(int ni = 0; ni < 2; ni++) {
          shift[di] = ni ? par->depl_decr[di] : -par->depl_decr[di];
          scRegionMediumMap::MediumPtr med2 = par->rm_map->at(p + shift);
          // this is possible bug (if *med and *med2 are different objects
          // with the same content)
          if(med2 && med2 != med) {
            low_recomb = true;
            return;
          }
        }
      }
    }

    virtual Recomb getRecombParams() const {
      Recomb rcm = scMedium::Local::getRecombParams();
      if (low_recomb)
        rcm.tau.fill(100);

      return rcm;
    }

    virtual valtype get_recomb(valtype n, valtype p) const {
      return Local::getRecombParams().get_recomb(n, p, Local::get_ni()*Local::get_ni());
    }
  };

  class LocalPool : public scAbstractMedium::LocalPool {
    std::vector<Local> elems;

  public:
    LocalPool(int nElements, const scMediumWithHQSurface* parent)
      : elems(nElements, Local(parent)) { }

    virtual Local& operator[](int n) { return elems[n]; }
    virtual const Local& operator[](int n) const { return elems[n]; }

    virtual size_t sizeInBytes() const {
      return sizeof(*this) + elems.size()*sizeof(Local);
    }
  };

  virtual ~scMediumWithHQSurface() { }

  virtual LocalPoolPtr getLocalMediumPool(int nElements) const {
    return LocalPoolPtr(new LocalPool(nElements, this));
  }

  virtual LocalPtr getLocalMedium() const {
    return LocalPtr(new Local(this));
  }
};

#endif
