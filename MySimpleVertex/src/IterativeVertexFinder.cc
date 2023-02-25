#include "MyVertexFinder.h"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "Acts/Vertexing/ZScanVertexFinder.hpp"

ROOT::VecOps::RVec<Acts::Vertex<Acts::BoundTrackParameters>> MyVertexFinder::IterativeVertexFinder(const ROOT::VecOps::RVec<lcio::TrackImpl> &inParticles) const
{
  // First thing, convert tracks to ACTS FORMAT
  int Ntr = (int)(inParticles.size());
  std::vector<Acts::BoundTrackParameters> allTracks;
  allTracks.reserve(Ntr);

  for (Int_t i = 0; i < Ntr; i++)
  {
    // get the LCIO track
    lcio::TrackImpl ti = inParticles[i];

    // convert LCIO track to FB track
    TVectorD trackFB = VertexingUtils::get_trackParam(ti);
    TMatrixDSym covFB = VertexingUtils::get_trackCov(ti);

    // use FB tool to convert to ACTS track
    TVectorD trackACTS = VertexingUtils::ParToACTS(trackFB);
    TMatrixDSym covACTS = VertexingUtils::CovToACTS(covFB, trackFB);

    // Acts::BoundTrackParameters::ParametersVector newTrackParams;
    Acts::BoundVector newTrackParams;
    newTrackParams << trackACTS[0], trackACTS[1], trackACTS[2], trackACTS[3], trackACTS[4], trackACTS[5];

    for (int ii = 0; ii < 6; ii++)
    {
      for (int jj = 0; jj < 6; jj++)
      {
        covACTS(jj, ii) = double(covACTS(jj, ii));
      }
    }

    // Get track covariance vector
    using Covariance = Acts::BoundSymMatrix;
    Covariance covMat;
    covMat << covACTS(0, 0), covACTS(1, 0), covACTS(2, 0), covACTS(3, 0), covACTS(4, 0), covACTS(5, 0),
        covACTS(0, 1), covACTS(1, 1), covACTS(2, 1), covACTS(3, 1), covACTS(4, 1), covACTS(5, 1),
        covACTS(0, 2), covACTS(1, 2), covACTS(2, 2), covACTS(3, 2), covACTS(4, 2), covACTS(5, 2),
        covACTS(0, 3), covACTS(1, 3), covACTS(2, 3), covACTS(3, 3), covACTS(4, 3), covACTS(5, 3),
        covACTS(0, 4), covACTS(1, 4), covACTS(2, 4), covACTS(3, 4), covACTS(4, 4), covACTS(5, 4),
        covACTS(0, 5), covACTS(1, 5), covACTS(2, 5), covACTS(3, 5), covACTS(4, 5), covACTS(5, 5);

    // Create track parameters and add to track list
    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface;
    Acts::Vector3 beamspotPos;
    beamspotPos << 0.0, 0.0, 0.0;
    perigeeSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(beamspotPos);

    allTracks.emplace_back(perigeeSurface, newTrackParams, std::move(covMat));
  }

  // retrieve input tracks and convert into the expected format
  std::vector<const Acts::BoundTrackParameters *> trackParametersPointers;
  trackParametersPointers.reserve(allTracks.size());
  for (const auto &trackParam : allTracks)
  {
    trackParametersPointers.push_back(&trackParam);
  }

  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  using ImpactPointEstimator =
      Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;
  using VertexSeeder = Acts::ZScanVertexFinder<VertexFitter>;
  using VertexFinder = Acts::IterativeVertexFinder<VertexFitter, VertexSeeder>;
  using VertexFinderOptions = Acts::VertexingOptions<Acts::BoundTrackParameters>;

  // Set up EigenStepper
  Acts::EigenStepper<> stepper(magneticField());

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(std::move(vertexFitterCfg));
  // Setup the track linearizer
  Linearizer::Config linearizerCfg(magneticField(), propagator);
  Linearizer linearizer(std::move(linearizerCfg));
  // Setup the seed finder
  ImpactPointEstimator::Config ipEstCfg(magneticField(), propagator);
  ImpactPointEstimator ipEst(std::move(ipEstCfg));
  VertexSeeder::Config seederCfg(ipEst);
  VertexSeeder seeder(std::move(seederCfg));
  // Set up the actual vertex finder
  VertexFinder::Config finderCfg(std::move(vertexFitter), std::move(linearizer),
                                 std::move(seeder), ipEst);
  finderCfg.maxVertices = 200;
  finderCfg.reassignTracksAfterFirstFit = true;
  VertexFinder finder(finderCfg);
  VertexFinder::State state(*magneticField(), magneticFieldContext());
  VertexFinderOptions finderOpts(geometryContext(), magneticFieldContext());

  // find vertices and measure elapsed time
  auto result = finder.find(trackParametersPointers, finderOpts, state);
  ROOT::VecOps::RVec<Acts::Vertex<Acts::BoundTrackParameters>> myPVs;

  if (not result.ok())
  {
    Acts::Vertex<Acts::BoundTrackParameters> myPV;
    myPV.setFullPosition({0., 0., 0., 0.});
    myPVs.push_back(myPV);

    return myPVs;
  }
  auto vertices = *result;

  // show some debug output
  streamlog_out(DEBUG2) << "Found " << vertices.size() << " vertices in event" << std::endl;

  for (const auto &vtx : vertices)
  {
    streamlog_out(DEBUG2) << "Found vertex at " << vtx.fullPosition().transpose() << " with "
                          << vtx.tracks().size() << " tracks." << std::endl;

    myPVs.push_back(vtx);
  }

  if (myPVs.size() == 0)
  {
    // If no PVs are found, take the beamspot
    Acts::Vertex<Acts::BoundTrackParameters> myPV;
    myPV.setFullPosition({0., 0., 0., 0.});

    myPVs.push_back(myPV);
  }

  return myPVs;
}