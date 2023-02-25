#ifndef VERTEXFITTERACTS_ANALYZERS_H
#define VERTEXFITTERACTS_ANALYZERS_H

#include <cmath>
#include <vector>

#include "EVENT/Track.h"
#include "EVENT/Vertex.h"
#include <IMPL/TrackImpl.h>

#include "ROOT/RVec.hxx"
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include "VertexingUtils.h"

namespace VertexFitterActs
{
  VertexingUtils::FCCAnalysesVertex VertexFitterFullBilloir(ROOT::VecOps::RVec<lcio::TrackImpl> tracks);
}

#endif
