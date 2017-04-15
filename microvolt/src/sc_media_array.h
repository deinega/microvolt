#ifndef SC_MEDIA_ARRAY_H
#define SC_MEDIA_ARRAY_H

/// @file sc_media_array.h Classes for storage of media

#include <cmath>
#include <cassert>
#include <map>

#include "sc_mediumbase.h"
#include "component.h"

typedef scAbstractMedium::Local local_t;

// In order to produce coinciding results we use the same
// underlying level->MediumRegion map (with reverse iteration
// instead of reverse ordering). For better performance
// forward iteration is preferred.
// The difference between different mappings is seen at points
// covered by two or more regions of the same level (e.g. on
// boundaries). For the forward-iterated map, the algorithm
// chooses the medium which was added first, for the backward-iterated
// map it chooses the last one 
#define USE_REVERSE_MEDIA_MAP

class scRegionMediumMap : public apComponent {
public:
  // we can't use std::auto_ptr - the pointers are put into containers
  typedef SHARED_OR_RAWPTR(const scAbstractMedium) MediumPtr;
  typedef SHARED_OR_RAWPTR(const Region_3) RegionPtr;

#ifdef HAS_SMART_PTR
  /// add a medium region
  /// level is region priority
  void add(const MediumPtr& med, const RegionPtr& reg, int level = 0) {
    meds.insert(make_pair(level, MediumRegion(med, reg)));
  }
#endif

  /// add a medium region
  /// level is region priority
  /// legacy version with a raw pointer for the region
  void add(SHARED_OR_AUTO(const scAbstractMedium) med, const Region_3 *reg, int level = 0) {
#ifdef HAS_SMART_PTR
    add(med, RegionPtr(reg), level);
#else
    MediumRegion mr(med.release(), reg);
    medStorage.add(mr);
    meds.insert(make_pair(level, mr));
#endif
  }

  /// return a medium at the specified point
  MediumPtr at(const Vector_3& point) const {
#ifdef USE_REVERSE_MEDIA_MAP
    base_map::const_reverse_iterator mit = meds.rbegin(), me = meds.rend();
#else
    base_map::const_iterator mit = meds.begin(), me = meds.end();
#endif
    while (mit != me && mit->second.reg && !mit->second.reg->TestPoint(point))
      ++mit;

    return mit == me ? MediumPtr() : mit->second.med;
  }

  /// Return the list of regions used
  std::vector<const Region_3*> get_regions() const {
    std::vector<const Region_3*> result;
    for(base_map::const_iterator mit = meds.begin(), me = meds.end(); mit != me; ++mit)
#ifdef HAS_SMART_PTR
      result.push_back(mit->second.reg.get());
#else
      result.push_back(mit->second.reg);
#endif
    return result;
  }

  /// Return the list of media used
  std::vector<const scAbstractMedium*> get_media() const {
    std::set<const scAbstractMedium*> result;
    for(base_map::const_iterator mit = meds.begin(), me = meds.end(); mit != me; ++mit)
#ifdef HAS_SMART_PTR
      result.insert(mit->second.med.get());
#else
      result.insert(mit->second.med);
#endif
    return std::vector<const scAbstractMedium*>(result.begin(), result.end());
  }

  void clear(){
    meds.clear();
#ifndef HAS_SMART_PTR
    medStorage.clear();
#endif
  }

private:
  struct MediumRegion {
    MediumPtr med;
    RegionPtr reg;

    explicit MediumRegion(const MediumPtr& med_ = MediumPtr(), const RegionPtr& sreg = RegionPtr())
      : med(med_), reg(sreg) { }
  };

  typedef std::multimap<int, MediumRegion
#ifndef USE_REVERSE_MEDIA_MAP
    , std::greater<int>
#endif
  > base_map;

  base_map meds;

#ifndef HAS_SMART_PTR
  class Storage {
    // Disable copying
    Storage(const Storage&);
    const Storage& operator=(const Storage&);

    std::vector<const scAbstractMedium*> meds;
    std::vector<const Region_3*> regs;
  public:
    Storage() {}

    void add(const MediumRegion& mr) {
      meds.push_back(mr.med);
      regs.push_back(mr.reg);
    }

    void clear(){
      for(std::vector<const scAbstractMedium*>::size_type i = 0; i < meds.size(); i++)
        delete meds[i];
      for(std::vector<const Region_3*>::size_type i = 0; i < regs.size(); i++)
        delete regs[i];
      meds.clear();
      regs.clear();
    }

    ~Storage() {
      clear();
    }

  } medStorage;
#endif
};

// this class allows obtaining MLocal by index number (in mesh)
// it organizes optimal storages of different media and their indexation
class scMediaArray {
  typedef scAbstractMedium::LocalPoolPtr PoolPtr;
  typedef scAbstractMedium::Local MLocal;

#ifdef HAS_SMART_PTR
  std::vector<PoolPtr> pools; // why we need to store it here?
#else
  PoolPtr* pools;
  size_t nPools;
#endif
  // pointers on MLocals enumareted according to corresponding mesh indices
  std::vector<const MLocal*> med;

  // The class is non-copyable (it can be made moveable, do we need it?)
  scMediaArray(const scMediaArray&) DELETED;
  const scMediaArray& operator=(const scMediaArray&) DELETED;

public:
  typedef std::vector<const MLocal*>::const_iterator const_iterator;

#ifdef HAS_SMART_PTR
  scMediaArray() { }
#else
  scMediaArray() : pools(), nPools(0) { }
#endif

  const MLocal& operator[](int n) const { 
    return *(med[n]); 
  }

  const_iterator begin() const { return med.begin(); }
  const_iterator end() const { return med.end(); }

  size_t sizeInBytes() const {
#ifdef HAS_SMART_PTR
    size_t nPools = pools.size();
#endif

    size_t size = sizeof(*this) +
      sizeof(const MLocal*)*med.size() +
      sizeof(PoolPtr)*nPools;

    for(size_t i = 0; i < nPools; i++)
      size += pools[i]->sizeInBytes();

    return size;
  }

  // Map the array to a part of another media array without copying
  void map_to(const_iterator medBegin, const_iterator medEnd) {
    clear();
    med.assign(medBegin, medEnd);
  }

  void map_to(const scMediaArray& m, size_t offset, size_t count) {
    map_to(m.begin() + offset, m.begin() + offset + count);
  }

  // Append without copying
  void append_refs_to(const_iterator medBegin, const_iterator medEnd) {
    med.insert(med.end(), medBegin, medEnd);
  }

  // Append without copying
  void append_ref_to(const MLocal* medLocal) {
    med.push_back(medLocal);
  }

  // Create a pool of local medium data for points specified with an iterator range
  // Media regions allocation is taken from media
  template<class InputIterator> 
  int build(const scRegionMediumMap& media, 
      InputIterator positionBegin, InputIterator positionEnd) {
    clear();
    return append(media, positionBegin, positionEnd);
  }

  // Append one point to the pool of local medium data for points specified with an iterator range
  int append(const scRegionMediumMap& media, const Vector_3& p) {
    return append(media, &p, &p + 1);
  }

  // Append points to the pool of local medium data for points specified with an iterator range
  template<class InputIterator> 
  int append(const scRegionMediumMap& media, 
      InputIterator positionBegin, InputIterator positionEnd) {
    typedef scRegionMediumMap::MediumPtr MediumPtr;

    // Lists of points and their indices for each medium
    typedef std::map<MediumPtr, pair<vector<Vector_3>, vector<int> > > PointsMap;
    PointsMap medPoints;

    std::vector<const MLocal*>::size_type newMedSize = med.size();
    MediumPtr mprev = MediumPtr();
    for(InputIterator it = positionBegin; it != positionEnd; ++it, ++newMedSize) {
      const Vector_3& p = *it;
      MediumPtr m = media.at(p); //get medium at position p

      if (!m){
        message(vblMESS1,-1,"Calculated space is not entirely covered by media, using previous point\n");
        m = mprev;
        if(!m)
          return LOGERR(-1,"Previous point not found!\n",0);
      }

      PointsMap::mapped_type& medPoint = medPoints[m];
      medPoint.first.push_back(p); // space position
      medPoint.second.push_back(newMedSize); // index
      mprev = m;
    }

    med.resize(newMedSize);
    
#ifndef HAS_SMART_PTR
    PoolPtr* oldPools = pools;
    pools = new PoolPtr[nPools + medPoints.size()];
    for(size_t i = 0; i < nPools; i++)
      pools[i] = oldPools[i];
    delete[] oldPools;
#endif
    for(PointsMap::const_iterator it = medPoints.begin(), e = medPoints.end(); it != e; ++it) {
      const std::vector<Vector_3>& points = it->second.first;
      scAbstractMedium::LocalPoolPtr pool =
        it->first->createLocalMediumPool(points.begin(), points.end());

      const std::vector<int>& index = it->second.second;
      for(std::vector<int>::size_type i = 0; i < index.size(); i++)
        med[index[i]] = &((*pool)[i]);

#ifdef HAS_SMART_PTR
      pools.push_back(std::move(pool));
#else
      pools[nPools++] = pool;
#endif
    }
    
    return 1;
  }

#ifdef HAS_SMART_PTR
  void clear() {
    med.clear();
    pools.clear();
  }
#else
  void clear() {
    med.clear();
    if(pools)
      delete[] pools;
    pools=NULL;
    nPools=0;
  }

  ~scMediaArray() {
    clear();
  }
#endif
};

#endif
