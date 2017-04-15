#include "sc_dump.h"
#include "sc_fd.h"
#include "sc_mesh.h"
#include "detector.h"

void scAdd2Dump(){

  AddDumpComponent<apDump, scRegionMediumMap>();
  AddDumpComponent<apDump, scBoundaryRegion>();
//  AddDumpComponent<apDump, scBlockContainer<RectBlock> >();

  AddDumpComponent<apDump,DetectorSet<vector<Vector_3>, scFixInterp<scBlockContainer<> > > >();
  AddDumpComponent<apDump,DetectorSet<SpaceVectorSet, scFixInterp<scBlockContainer<> > > >();
  AddDumpComponent<apDump,DetectorSet<UniformGrid<Vector_3>, scFixInterp<scBlockContainer<> > > >();
  AddDumpComponent<apDump,DetectorSet<BoxSurfaceSet, scFixInterp<scBlockContainer<> > > >();
  AddDumpComponent<apDump,DetectorSet<CylinderSurfaceSet, scFixInterp<scBlockContainer<> > > >();
}

void DumpComponent(const scRegionMediumMap *comp, apDump *dump, bool lim){
  vector<const Region_3*> regions = comp->get_regions();
  for(vector<const Region_3*>::const_iterator
      it = regions.begin(), e = regions.end(); it != e; ++it) {
    dump->DumpRegion(*it, "med", lim);
  }
}

void DumpComponent(const scBoundaryRegion *comp,apDump *dump, bool lim){
  dump->DumpRegion(comp->reg.ptr(), comp->type==1 ? "cont" : "di" ,lim);
}

/*template<class block_t>
void DumpComponent(const scBlockContainer<block_t> *comp,apDump *dump, bool lim){
  dump->DumpRegion(&(comp->GetRegion()),"cell",lim);
}*/

//template void DumpComponent(const scBlockContainer<> *comp,apDump *dump, bool lim);

template<class vset_tt, class xfer_t>
void DumpComponent(const DetectorSet<vset_tt,xfer_t> *comp,apDump *dump, bool lim){
  dump->DumpVectorSet(comp->GetVectorSet(),comp->GetVecTransform(),comp,"vset.d",lim);
}

template void DumpComponent(const DetectorSet<vector<Vector_3>, scFixInterp<scBlockContainer<> > > *comp,apDump *dump, bool lim);
template void DumpComponent(const DetectorSet<SpaceVectorSet, scFixInterp<scBlockContainer<> > > *comp,apDump *dump, bool lim);
template void DumpComponent(const DetectorSet<UniformGrid<Vector_3>, scFixInterp<scBlockContainer<> > > *comp,apDump *dump, bool lim);
template void DumpComponent(const DetectorSet<BoxSurfaceSet, scFixInterp<scBlockContainer<> > > *comp,apDump *dump, bool lim);
template void DumpComponent(const DetectorSet<CylinderSurfaceSet, scFixInterp<scBlockContainer<> > > *comp,apDump *dump, bool lim);
