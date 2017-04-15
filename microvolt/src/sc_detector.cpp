#include "detector.hpp"
#include "sc_fd.h"
#include "sc_mesh.h"

/// instantiation of used detectors

template class TextRecorder<vector<Vector_3> >;
template class TextFluxRecorder<SpaceVectorSet>;
template class TextFluxRecorder<CylinderSurfaceSet>;

template class DetectorSet<vector<Vector_3>, scFixInterp<scBlockContainer<> > >;
template class DetectorSet<SpaceVectorSet, scFixInterp<scBlockContainer<> > >;
template class DetectorSet<UniformGrid<Vector_3>, scFixInterp<scBlockContainer<> > >;
template class DetectorSet<BoxSurfaceSet, scFixInterp<scBlockContainer<> > >;
template class DetectorSet<CylinderSurfaceSet, scFixInterp<scBlockContainer<> > >;
