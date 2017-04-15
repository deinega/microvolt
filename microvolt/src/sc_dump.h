#ifndef SC_DUMP_H
#define SC_DUMP_H

/** @file sc_dump.h Functions used to record information about all necessary components
 to text files in gnuplot and VTK format (for visualization)
 (see also component.h and gnudump.h) */

#include "component.h"

/// calls AddDumpComponent (see component.h) for all necessary components
void scAdd2Dump();

/// emty definition 
template<class comp_t>
void DumpComponent(const comp_t *comp,apDump *dump=theDump, bool lim=false){}

/// below there are instantiations of DumpComponent for all necessary components

template<class vset_tt,class xfer_t>
class DetectorSet;

template<class vset_tt, class xfer_t>
void DumpComponent(const DetectorSet<vset_tt,xfer_t> *comp,apDump *dump=theDump, bool lim=false);

class scRegionMediumMap;
void DumpComponent(const scRegionMediumMap *comp, apDump *dump = theDump, bool lim = false);

struct scBoundaryRegion;
void DumpComponent(const scBoundaryRegion *comp,apDump *dump=theDump, bool lim=false);

template<class block_t>
class scBlockContainer;

template<class block_t>
void DumpComponent(const scBlockContainer<block_t> *comp,apDump *dump=theDump, bool lim=false);

#endif
