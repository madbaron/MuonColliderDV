#include "MyVertexFinder.h"

#include "TMath.h"
#include "TH1.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "IMPL/TrackImpl.h"

// ACTS
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include <Acts/Vertexing/ImpactPointEstimator.hpp>
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//--------------------------------------------------------
//   Template routine for 2track secondary vertices selection
//

using namespace Acts;
using namespace Acts::UnitLiterals;

void MyVertexFinder::select2TrVrt(ROOT::VecOps::RVec<Acts::BoundTrackParameters> &allTracks,
                                  const ROOT::VecOps::RVec<float> track_innermostRadius,
                                  const Acts::Vertex<Acts::BoundTrackParameters> &primVrt) const
{
   using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
   using Estimator = Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;

   // Set up EigenStepper
   Acts::EigenStepper<> stepper(magneticField());

   // Create the IP Estimator
   Estimator::Config cfg(magneticField(), std::make_shared<Propagator>(std::move(stepper)));
   Estimator ipEstimator = Estimator(cfg);

   // Set up EigenStepper
   Acts::EigenStepper<> vstepper(magneticField());

   using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
   Linearizer::Config ltConfig(magneticField(), std::make_shared<Propagator>(std::move(vstepper)));
   Linearizer linearizer(ltConfig);

   int Ntr = (int)(allTracks.size());

   // Constraint for vertex fit
   Acts::Vertex<Acts::BoundTrackParameters> myPV = primVrt;
   streamlog_out(DEBUG0) << "PV at " << myPV.position()[0] << " " << myPV.position()[1] << " " << myPV.position()[2] << std::endl;

   double signifR = 0., signifZ = 0.;
   std::vector<double> trackSignif(Ntr), dRdZratio(Ntr);
   std::vector<double> trackd0(Ntr);
   std::vector<double> trackd0sig(Ntr);

   for (int i = 0; i < Ntr; i++)
   {
      // check that computed impact parameters are meaningful
      auto result = ipEstimator.estimateImpactParameters(allTracks[i], myPV, geometryContext(), magneticFieldContext());
      if (not result.ok())
      {
         trackSignif[i] = 0.;
         dRdZratio[i] = 0.;
         continue;
      }

      ImpactParametersAndSigma impact = result.value();

      // Acts::ImpactParametersAndSigma impact = ipEstimator.estimateImpactParameters(trackACTS, myPV, geometryContext(), magneticFieldContext());
      // signifR = impact.IPd0 / impact.PVsigmad0;
      // signifZ = impact.IPz0 / impact.PVsigmaz0;
      signifR = impact.IPd0 / impact.sigmad0;
      signifZ = impact.IPz0 / impact.sigmaz0;
      trackSignif[i] = sqrt(signifR * signifR + signifZ * signifZ);
      dRdZratio[i] = std::abs(signifR / signifZ);
      trackd0[i] = impact.IPd0;
      trackd0sig[i] = signifR;

      // if (trackd0[i] > m_trkMind0 && signifR > m_trkSigCut)
      // if (trackd0[i] > m_trkMind0 && trackSignif[i] > m_trkSigCut && dRdZratio[i] > m_dRdZRatioCut)
      //{
      streamlog_out(DEBUG0) << "Track pT " << allTracks[i].transverseMomentum() << " IPd0 " << impact.IPd0 << " sd0 " << impact.sigmad0 << " pvsd0 " << impact.PVsigmad0 << " dRdz " << dRdZratio[i] << std::endl;
      //}
   }

   // Set up Billoir Vertex Fitter
   using VertexFitter = Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
   VertexFitter::Config vertexFitterCfg;
   vertexFitterCfg.maxIterations = 10;

   VertexFitter billoirFitter(vertexFitterCfg);

   // Build initial graph
   for (int i = 0; i < Ntr - 1; i++)
   {
      if (fabs(trackd0[i]) < m_trkMind0)
         // if (fabs(trackd0[i]) < m_trkMind0 || trackSignif[i] < m_trkSigCut || dRdZratio[i] < m_dRdZRatioCut)
         // if (fabs(trackd0[i]) < m_trkMind0 || trackd0sig[i] < m_trkSigCut)
         continue;
      for (int j = i + 1; j < Ntr; j++)
      {
         if (fabs(trackd0[j]) < m_trkMind0)
            // if (fabs(trackd0[j]) < m_trkMind0 || trackSignif[j] < m_trkSigCut || dRdZratio[j] < m_dRdZRatioCut)
            //  if (fabs(trackd0[j]) < m_trkMind0 || trackd0sig[j] < m_trkSigCut)
            continue;

         if (std::abs(allTracks[i].transverseMomentum() - allTracks[j].transverseMomentum()) == 0)
            continue; // remove duplicated tracks

         TVector3 p3_i = {allTracks[i].momentum()[0], allTracks[i].momentum()[1], allTracks[i].momentum()[2]};
         TLorentzVector p4_i;
         p4_i.SetPtEtaPhiM(allTracks[i].transverseMomentum(), p3_i.PseudoRapidity(), p3_i.Phi(), m_massPi);

         TVector3 p3_j = {allTracks[j].momentum()[0], allTracks[j].momentum()[1], allTracks[j].momentum()[2]};
         TLorentzVector p4_j;
         p4_j.SetPtEtaPhiM(allTracks[j].transverseMomentum(), p3_j.PseudoRapidity(), p3_j.Phi(), m_massPi);
         TLorentzVector sumTrkP = p4_i + p4_j;

         float ihitR = track_innermostRadius[i];
         float jhitR = track_innermostRadius[j];
         if (std::abs(ihitR - jhitR) > m_deltaRadiusInnermostHits2VrtCut)
         {
            streamlog_out(DEBUG0) << "Fail innermost hit pair distance " << std::abs(ihitR - jhitR) << " R1 " << ihitR << " R2 " << jhitR << std::endl;
            continue; //- FMPs are in very different layers
         }

         std::vector<const Acts::BoundTrackParameters *>
             tracksForFit;
         tracksForFit.push_back(&allTracks[i]);
         tracksForFit.push_back(&allTracks[j]);

         VertexFitter::State state(magneticField()->makeCache(magneticFieldContext()));
         Acts::VertexingOptions<Acts::BoundTrackParameters> vfOptions(geometryContext(), magneticFieldContext());
         // Acts::VertexingOptions<Acts::BoundTrackParameters> vfOptionsConstr(geometryContext(), magneticFieldContext(), myPV);

         auto result = billoirFitter.fit(tracksForFit, linearizer, vfOptions, state);
         if (not result.ok())
         {
            streamlog_out(DEBUG0) << "Can't fit these tracks " << std::endl;
            continue;
         }

         Acts::Vertex<Acts::BoundTrackParameters> fittedVertex = result.value();

         double Prob2v = TMath::Prob(fittedVertex.fitQuality().first, fittedVertex.fitQuality().second);
         if (Prob2v < m_sel2VrtProbCut)
         {
            streamlog_out(DEBUG0) << "Prob too low " << Prob2v << " (" << fittedVertex.fitQuality().first << " " << fittedVertex.fitQuality().second << ")" << std::endl;
            continue;
         }

         double vertexRadius = sqrt(fittedVertex.position()[0] * fittedVertex.position()[0] + fittedVertex.position()[1] * fittedVertex.position()[1]);
         if (vertexRadius > m_maxSVRadiusCut)
         {
            streamlog_out(DEBUG0) << "Vertex radius in IT volume " << vertexRadius << std::endl;
            continue; // Too far from interaction point
         }

         double cosSVPV = projSV_PV({fittedVertex.position()[0], fittedVertex.position()[1], fittedVertex.position()[2]}, {myPV.position()[0], myPV.position()[1], myPV.position()[2]}, sumTrkP);
         if (cosSVPV < m_cosSVPVCut)
         {
            streamlog_out(DEBUG0) << "DV decays backwards " << cosSVPV << std::endl;
            continue;
         }

         if (sumTrkP.Pt() < 1.)
         {
            streamlog_out(DEBUG0) << "DV pT " << sumTrkP.Pt() << std::endl;
            continue;
         } // Add cuts to reject material regions

         // Save good candidate for multi-vertex fit
         add_edge(i, j, *m_compatibilityGraph);
      }
   }

   return;
}
