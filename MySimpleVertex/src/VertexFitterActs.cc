#include "VertexFitterActs.h"

// ACTS
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/Vertex.hpp"

//V5.0
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

#include "TMath.h"

using namespace VertexFitterActs;
using namespace Acts::UnitLiterals;

VertexingUtils::FCCAnalysesVertex VertexFitterActs::VertexFitterFullBilloir(ROOT::VecOps::RVec<lcio::TrackImpl> tracks)
{

  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;

  // Create a test context
  //Acts::GeometryContext geoContext = Acts::GeometryContext();
  //Acts::MagneticFieldContext magFieldContext = Acts::MagneticFieldContext();
  const auto &geoContext = Acts::GeometryContext();
  const auto &magFieldContext = Acts::MagneticFieldContext();

  // Set up EigenStepper
  //Acts::ConstantBField bField(Acts::Vector3(0., 0., 2_T));
  //Acts::EigenStepper<Acts::ConstantBField> stepper(bField);
  auto bField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0., 0., 3.56_T));
  Acts::EigenStepper<> stepper(bField);

  // Set up the propagator
  //using Propagator = Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>>;
  auto propagator = std::make_shared<Propagator>(stepper);

  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Set up Billoir Vertex Fitter
  using VertexFitter = Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  VertexFitter::Config vertexFitterCfg;
  VertexFitter billoirFitter(vertexFitterCfg);
  //VertexFitter::State state(magFieldContext);
  VertexFitter::State state(bField->makeCache(magFieldContext));

  // Constraint for vertex fit
  Acts::Vertex<Acts::BoundTrackParameters> myConstraint;
  // Some abitrary values
  Acts::SymMatrix4 myCovMat = Acts::SymMatrix4::Zero();
  myCovMat(0, 0) = 30.;
  myCovMat(1, 1) = 30.;
  myCovMat(2, 2) = 30.;
  myCovMat(3, 3) = 30.;
  myConstraint.setFullCovariance(std::move(myCovMat));
  myConstraint.setFullPosition(Acts::Vector4(0, 0, 0, 0));

  Acts::VertexingOptions<Acts::BoundTrackParameters> vfOptions(geoContext, magFieldContext);
  Acts::VertexingOptions<Acts::BoundTrackParameters> vfOptionsConstr(geoContext, magFieldContext, myConstraint);

  int Ntr = tracks.size();

  std::vector<Acts::BoundTrackParameters> allTracks;
  allTracks.reserve(Ntr);

  for (Int_t i = 0; i < Ntr; i++)
  {
    //get the LCIO track
    lcio::TrackImpl ti = tracks[i];

    //convert LCIO track to FB track
    TVectorD trackFB = VertexingUtils::get_trackParam(ti);
    TMatrixDSym covFB = VertexingUtils::get_trackCov(ti);

    //use FB tool to convert to ACTS track
    TVectorD trackACTS = VertexingUtils::ParToACTS(trackFB);
    TMatrixDSym covACTS = VertexingUtils::CovToACTS(covFB, trackFB);

    //Acts::BoundTrackParameters::ParametersVector newTrackParams;
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

  //std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
  //int count = 0;
  //for (const auto& trk : allTracks) {
  //  std::cout << count << ". track: " << std::endl;
  //  std::cout << "params: " << trk << std::endl;
  //  count++;
  // }

  std::vector<const Acts::BoundTrackParameters *> tracksPtr;
  for (const auto &trk : allTracks)
  {
    tracksPtr.push_back(&trk);
  }

  VertexingUtils::FCCAnalysesVertex TheVertex;
  lcio::VertexImpl edm4hep_vertex;
  ROOT::VecOps::RVec<float> reco_chi2;
  ROOT::VecOps::RVec<TVectorD> updated_track_parameters;
  ROOT::VecOps::RVec<int> reco_ind;
  ROOT::VecOps::RVec<float> final_track_phases;
  ROOT::VecOps::RVec<TVector3> updated_track_momentum_at_vertex;

  TheVertex.vertex = edm4hep_vertex;
  TheVertex.reco_chi2 = reco_chi2;
  TheVertex.updated_track_parameters = updated_track_parameters;
  TheVertex.reco_ind = reco_ind;
  TheVertex.final_track_phases = final_track_phases;
  TheVertex.updated_track_momentum_at_vertex = updated_track_momentum_at_vertex;

  TheVertex.ntracks = Ntr;
  if (Ntr <= 1)
  {
    return TheVertex; // can not reconstruct a vertex with only one track...
  }

  Acts::Vertex<Acts::BoundTrackParameters> fittedVertex =
      billoirFitter.fit(tracksPtr, linearizer, vfOptions, state).value();
  //Acts::Vertex<Acts::BoundTrackParameters> fittedVertexConstraint =
  //  billoirFitter.fit(tracksPtr, linearizer, vfOptionsConstr, state).value();

  //std::cout << "Fitting nTracks: " << Ntr << std::endl;
  //  std::cout << "Fitted Vertex: " << fittedVertex.position() << std::endl;
  //  std::cout << "Fitted constraint Vertex: "
  //             << fittedVertexConstraint.position() << std::endl;

  TheVertex.ntracks = Ntr;
  edm4hep_vertex.setChi2(fittedVertex.fitQuality().first / fittedVertex.fitQuality().second);
  //std::cout << "Fitted chi2: " << edm4hep_vertex.chi2 << std::endl;

  float vpos[3] = {fittedVertex.position()[0], fittedVertex.position()[1], fittedVertex.position()[2]};
  edm4hep_vertex.setPosition(vpos);

  auto vtxCov = fittedVertex.covariance();
  std::vector<float> covMatrix; // covMat in lcio is a LOWER-triangle matrix.
  covMatrix.push_back(vtxCov(0, 0));
  covMatrix.push_back(vtxCov(1, 0));
  covMatrix.push_back(vtxCov(1, 1));
  covMatrix.push_back(vtxCov(2, 0));
  covMatrix.push_back(vtxCov(2, 1));
  covMatrix.push_back(vtxCov(2, 2));

  edm4hep_vertex.setAlgorithmType("3");
  edm4hep_vertex.setCovMatrix(covMatrix);

  TheVertex.vertex = edm4hep_vertex;

  return TheVertex;
}
