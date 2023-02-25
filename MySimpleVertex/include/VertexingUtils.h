#ifndef VERTEXINGUTILS_ANALYZERS_H
#define VERTEXINGUTILS_ANALYZERS_H

#include <cmath>
#include <vector>

#include "lcio.h"

#include "ROOT/RVec.hxx"
#include "EVENT/Track.h"
#include "IMPL/TrackImpl.h"
#include "EVENT/Vertex.h"
#include "IMPL/VertexImpl.h"

#include "EVENT/ReconstructedParticle.h"

#include "TVectorD.h"
#include "TVector3.h"
#include "TMatrixDSym.h"

/** Vertexing utilities 
*/

namespace VertexingUtils
{

  /// Structure to keep useful track information that is related to the vertex
  struct FCCAnalysesVertex
  {
    lcio::VertexImpl vertex;
    int ntracks;
    int mc_ind; ///index in the MC vertex collection if any
    ROOT::VecOps::RVec<int> reco_ind;
    ROOT::VecOps::RVec<float> reco_chi2;
    ROOT::VecOps::RVec<TVector3> updated_track_momentum_at_vertex;
    ROOT::VecOps::RVec<TVectorD> updated_track_parameters;
    ROOT::VecOps::RVec<float> final_track_phases;
  };

  /// Retrieve the number of reconstructed vertices from the collection of vertex object
  int get_Nvertex(ROOT::VecOps::RVec<FCCAnalysesVertex> TheVertexColl);

  /// Retrieve a single FCCAnalyses vertex from the collection of vertex object
  FCCAnalysesVertex get_FCCAnalysesVertex(ROOT::VecOps::RVec<FCCAnalysesVertex> TheVertexColl, int index);

  /// Retrieve the lcio::VertexImpl from the vertex object
  lcio::VertexImpl get_VertexData(FCCAnalysesVertex TheVertex);

  /// Retrieve a vector of lcio::VertexImpl from the collection of vertex object
  ROOT::VecOps::RVec<lcio::VertexImpl> get_VertexData(ROOT::VecOps::RVec<FCCAnalysesVertex> TheVertexColl);

  /// Retrieve a lcio::VertexImpl from the collection of vertex object at a given index
  lcio::VertexImpl get_VertexData(ROOT::VecOps::RVec<FCCAnalysesVertex> TheVertexColl, int index);

  /// Retrieve the number of tracks from FCCAnalysesVertex
  int get_VertexNtrk(FCCAnalysesVertex TheVertex);

  /// Retrieve the tracks indices from FCCAnalysesVertex
  ROOT::VecOps::RVec<int> get_VertexRecoInd(FCCAnalysesVertex TheVertex);

  /// Return the number of tracks in a given track collection
  int get_nTracks(ROOT::VecOps::RVec<lcio::TrackImpl> tracks);

  // --- Internal methods needed by the code of  Franco B :
  TVectorD get_trackParam(lcio::Track &atrack);
  TMatrixDSym get_trackCov(lcio::Track &atrack);

  TVectorD ParToACTS(TVectorD Par);
  TMatrixDSym CovToACTS(TMatrixDSym Cov, TVectorD Par);

}
#endif
