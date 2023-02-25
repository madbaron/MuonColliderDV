#include "MyVertexFinder.h"
#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include "IMPL/TrackImpl.h"

#include "VertexFitterSimple.h"
#include "VertexFitterActs.h"

#include "VertexingUtils.h"

#include "TMath.h"
#include "TVector3.h"

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

#include "boost/graph/bron_kerbosch_all_cliques.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

ROOT::VecOps::RVec<std::pair<Acts::Vertex<Acts::BoundTrackParameters>, std::deque<long int>>> MyVertexFinder::getAllDV(const ROOT::VecOps::RVec<lcio::TrackImpl> &inParticles,
                                                                                                                       const ROOT::VecOps::RVec<float> track_innermostRadius,
                                                                                                                       const ROOT::VecOps::RVec<int> track_mother,
                                                                                                                       const Acts::Vertex<Acts::BoundTrackParameters> &primaryVx) const
{
   // First thing, convert tracks to ACTS FORMAT
   int Ntr = (int)(inParticles.size());
   ROOT::VecOps::RVec<Acts::BoundTrackParameters> allTracks;
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

   // Now convert vertex
   Acts::Vertex<Acts::BoundTrackParameters> myPV = primaryVx;

   const double probVrtMergeLimit = 0.01;

   int i, j;
   int inNPart = inParticles.size();

   ROOT::VecOps::RVec<std::pair<Acts::Vertex<Acts::BoundTrackParameters>, std::deque<long int>>> finalVertices(0);

   if (inNPart < 2)
   {
      return finalVertices;
   }

   // TODO implement
   // selParticles = selGoodTrkParticle(inParticles, primaryVx);
   ROOT::VecOps::RVec<Acts::BoundTrackParameters> selParticles = allTracks;
   int nTracks = selParticles.size();

   if (nTracks < 2)
   {
      return finalVertices;
   }

   //  inpTrk[]           - input track list
   //  listSelTracks[]    - list of good tracks in jet for vertex search
   //------------------------------------------------------------
   //                     Initial track list ready
   //                     Find 2track vertices

   select2TrVrt(selParticles, track_innermostRadius, primaryVx);

   streamlog_out(DEBUG0) << " Defined edges in the graph: " << num_edges(*m_compatibilityGraph) << std::endl;

   //  m_Incomp[]           -  main vector of pointers for multivertex search
   //-----------------------------------------------------------------------------------------------------
   //            Secondary track list is ready
   //            Creation of initial vertex set

   std::unique_ptr<std::vector<WrkVrt>> wrkVrtSet = std::make_unique<std::vector<WrkVrt>>();
   WrkVrt newvrt;
   newvrt.Good = true;
   long int NPTR = 0, nth = 2; // VK nth=2 to speed up PGRAPH when it's used

   using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
   using Estimator = Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;

   // Set up EigenStepper
   Acts::EigenStepper<> vstepper(magneticField());

   using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
   Linearizer::Config ltConfig(magneticField(), std::make_shared<Propagator>(std::move(vstepper)));
   Linearizer linearizer(ltConfig);

   // Set up Billoir Vertex Fitter
   using VertexFitter = Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
   VertexFitter::Config vertexFitterCfg;
   VertexFitter billoirFitter(vertexFitterCfg);
   VertexFitter::State state(magneticField()->makeCache(magneticFieldContext()));
   Acts::VertexingOptions<Acts::BoundTrackParameters> vfOptions(geometryContext(), magneticFieldContext());

   std::vector<std::vector<int>> allCliques;
   bron_kerbosch_all_cliques(*m_compatibilityGraph, clique_visitor(allCliques));
   for (int cq = 0; cq < (int)allCliques.size(); cq++)
   {
      newvrt.selTrk.clear();
      NPTR = allCliques[cq].size();
      for (i = 0; i < NPTR; i++)
      {
         newvrt.selTrk.push_back(allCliques[cq][i]);
      }

      std::vector<const Acts::BoundTrackParameters *> tmpListTracks;
      std::vector<double> chi2PerTrk;
      TLorentzVector sum_p4;
      for (i = 0; i < NPTR; i++)
      {
         tmpListTracks.push_back(&selParticles.at(newvrt.selTrk[i]));
         TVector3 p3_i = {selParticles[newvrt.selTrk[i]].momentum()[0], selParticles[newvrt.selTrk[i]].momentum()[1], selParticles[newvrt.selTrk[i]].momentum()[2]};
         TLorentzVector p4_i;
         p4_i.SetPtEtaPhiM(selParticles[newvrt.selTrk[i]].transverseMomentum(), p3_i.PseudoRapidity(), p3_i.Phi(), m_massPi);
         sum_p4 += p4_i;
         chi2PerTrk.push_back(inParticles[newvrt.selTrk[i]].getChi2());
      }

      auto result = billoirFitter.fit(tmpListTracks, linearizer, vfOptions, state);
      if (not result.ok())
         continue;

      Acts::Vertex<Acts::BoundTrackParameters> fittedVertex = result.value();

      streamlog_out(DEBUG0) << "Fit DV at " << fittedVertex.position()[0] << " " << fittedVertex.position()[1] << " " << fittedVertex.position()[2] << std::endl;
      if (NPTR == 2 && fittedVertex.fitQuality().first > 10.)
      {
         streamlog_out(DEBUG0) << " bad fit quality " << fittedVertex.fitQuality().first << " " << NPTR << std::endl;
         continue;
      }
      newvrt.vertex = fittedVertex;
      newvrt.Good = true;
      newvrt.vertexMom = sum_p4;
      newvrt.chi2 = fittedVertex.fitQuality().first;
      newvrt.projectedVrt = MomProjDist({fittedVertex.position()[0], fittedVertex.position()[1], fittedVertex.position()[2]}, {primaryVx.position()[0], primaryVx.position()[1], primaryVx.position()[2]}, newvrt.vertexMom); // 3D SV-PV distance
      newvrt.chi2PerTrk = chi2PerTrk;

      wrkVrtSet->push_back(newvrt);
   }

   if ((*wrkVrtSet).size() == 0)
      return finalVertices;

   printWrkSet(wrkVrtSet.get(), "Initial Vertices");

   // Count track participation in different vertices
   std::vector<int> trkNPairs(nTracks, 0);
   for (auto &vrt : (*wrkVrtSet))
   {
      int ntInV = vrt.selTrk.size() - 1;
      for (auto &trk : vrt.selTrk)
         trkNPairs.at(trk) += ntInV;
   }

   // Resolve all overlapped vertices
   std::multimap<double, std::pair<int, int>> vrtWithCommonTrk;
   while (true)
   {
      int nSoluI = (*wrkVrtSet).size();
      vrtWithCommonTrk.clear();
      unsigned int nTComMax = 0;
      for (int iv = 0; iv < nSoluI - 1; iv++)
      {
         if (!(*wrkVrtSet)[iv].Good)
            continue;
         if ((*wrkVrtSet)[iv].selTrk.size() < nTComMax)
            continue; // Optimisation. Only biggest overlap matters
         for (int jv = iv + 1; jv < nSoluI; jv++)
         {
            if (!(*wrkVrtSet)[jv].Good)
               continue;
            if ((*wrkVrtSet)[jv].selTrk.size() < nTComMax)
               continue; // Optimisation. Only biggest overlap matters
            unsigned int nTCom = nTrkCommon(wrkVrtSet.get(), iv, jv);
            if (!nTCom)
               continue;
            if (nTCom < nTComMax)
               continue;
            double sumChi2 = (*wrkVrtSet)[iv].chi2 + (*wrkVrtSet)[jv].chi2;
            sumChi2 = std::min(sumChi2, 999.) * 1.e-3;
            vrtWithCommonTrk.emplace(nTCom + sumChi2, std::make_pair(iv, jv));
            nTComMax = std::max(nTComMax, nTCom);
         }
      }
      if (vrtWithCommonTrk.size() == 0)
         break;

      unsigned int nTCom = (*vrtWithCommonTrk.rbegin()).first;
      WrkVrt &v1 = (*wrkVrtSet)[(*vrtWithCommonTrk.rbegin()).second.first];
      WrkVrt &v2 = (*wrkVrtSet)[(*vrtWithCommonTrk.rbegin()).second.second];
      //--First check if one vertex is fully contained in another
      if (nTCom == v1.selTrk.size() || nTCom == v2.selTrk.size())
      {
         if (nTCom == v1.selTrk.size())
         {
            v1.Good = false;
            continue;
         }
         if (nTCom == v2.selTrk.size())
         {
            v2.Good = false;
            continue;
         }
      }
      //--Then check if 2 vertices with common tracks can be simply merged
      if (nTCom > 1 && TMath::Prob(v1.chi2, 2 * v1.selTrk.size() - 3) > probVrtMergeLimit && TMath::Prob(v2.chi2, 2 * v2.selTrk.size() - 3) > probVrtMergeLimit)
      {
         double prbV = mergeAndRefitVertices(v1, v2, newvrt, selParticles);
         if (prbV > probVrtMergeLimit)
         {
            streamlog_out(DEBUG0) << "Merging vertices!" << std::endl;
            v1.Good = false;
            v2.Good = false;
            newvrt.Good = true;
            newvrt.projectedVrt = MomProjDist({newvrt.vertex.position()[0], newvrt.vertex.position()[1], newvrt.vertex.position()[2]}, {primaryVx.position()[0], primaryVx.position()[1], primaryVx.position()[2]}, newvrt.vertexMom); // 3D SV-PV distance
            newvrt.chi2 = newvrt.vertex.fitQuality().first;

            std::swap(v1, newvrt); // Replace v1 by new vertex
            continue;
         }
      }
      //--If not mergeable - refine them
      refineVerticesWithCommonTracks(v1, v2, selParticles);
   }

   /*
   // Clean duplicated 1track vertices if they exist
   if (m_multiWithOneTrkVrt)
   {
      for (auto &v1t : (*wrkVrtSet))
      {
         if (v1t.selTrk.size() != 1 || !v1t.Good)
            continue;
         int ind_t = v1t.selTrk[0];
         if (trkNPairs[ind_t] < 2)
         {
            v1t.Good = false;
            continue;
         } // Remove 1tr-vertex if track crosses only one other track
         if (selParticles[ind_t]->pt() < m_cutPt * 2)
         {
            v1t.Good = false;
            continue;
         }; // Tighten track_pt cut for 1-track vertex
         for (auto &vrt : (*wrkVrtSet))
         { // Check if the track is present in another vertex, including other 1-track ones
            if (!vrt.Good || &v1t == &vrt)
               continue;
            if (std::find(vrt.selTrk.begin(), vrt.selTrk.end(), ind_t) != vrt.selTrk.end())
            {
               v1t.Good = false;
               break;
            }
         }
      }
   }
   */

   // Remove all bad vertices from the working set
   int tmpV = 0;
   while (tmpV < (int)(*wrkVrtSet).size())
      if (!(*wrkVrtSet)[tmpV].Good)
      {
         (*wrkVrtSet).erase((*wrkVrtSet).begin() + tmpV);
      }
      else
      {
         tmpV++;
      }
   if ((*wrkVrtSet).size() == 0)
      return finalVertices;

   printWrkSet(wrkVrtSet.get(), "Intermediate Vertices");

   for (auto &tmpV : (*wrkVrtSet))
      tmpV.projectedVrt = MomProjDist({tmpV.vertex.position()[0], tmpV.vertex.position()[1], tmpV.vertex.position()[2]}, {primaryVx.position()[0], primaryVx.position()[1], primaryVx.position()[2]}, tmpV.vertexMom); // Setup ProjectedVrt
                                                                                                                                                                                                                       //----------------------------------------------------------------------------

   //----------------------------------------------------------------------------
   // Final check/merge for close vertices

   if ((*wrkVrtSet).size() > 1)
   {

      int foundV1 = -1, foundV2 = -1;
      std::vector<double> checkedDst(0);
      double minDistVV = minVrtVrtDist(wrkVrtSet.get(), foundV1, foundV2, checkedDst);

      // recalculate VV distances
      while (minDistVV < m_vertexMergeCut)
      {
         if (foundV1 < foundV2)
         {
            int tmp = foundV1;
            foundV1 = foundV2;
            foundV2 = tmp;
         }
         double probV = mergeAndRefitVertices((*wrkVrtSet)[foundV1], (*wrkVrtSet)[foundV2], newvrt, selParticles);
         streamlog_out(DEBUG0) << "Merged vertex prob=" << probV << " Vrt1=" << foundV1 << " Vrt2=" << foundV2 << " dst=" << minDistVV << std::endl;
         if (probV < probVrtMergeLimit)
         { //--- If merged vertex is bad - try to remove the worst track
            int pos = std::max_element(newvrt.chi2PerTrk.begin(), newvrt.chi2PerTrk.end()) - newvrt.chi2PerTrk.begin();
            streamlog_out(DEBUG0) << "Removing track " << pos << " with " << newvrt.chi2PerTrk[pos] << std::endl;

            newvrt.detachedTrack = newvrt.selTrk[pos];
            newvrt.selTrk.erase(newvrt.selTrk.begin() + pos);
            probV = refitVertex(newvrt, selParticles);
            streamlog_out(DEBUG0) << "Attempt to improve prob=" << probV << std::endl;
         }

         if (probV > probVrtMergeLimit)
         {
            //  Good merged vertex found
            newvrt.projectedVrt = MomProjDist({newvrt.vertex.position()[0], newvrt.vertex.position()[1], newvrt.vertex.position()[2]}, {primaryVx.position()[0], primaryVx.position()[1], primaryVx.position()[2]}, newvrt.vertexMom);
            std::swap((*wrkVrtSet)[foundV1], newvrt);
            (*wrkVrtSet)[foundV2].Good = false;   // Drop vertex
            (*wrkVrtSet)[foundV2].selTrk.clear(); // Clean dropped vertex
         }
         else
         {
            checkedDst.push_back(minDistVV);
         }

         minDistVV = minVrtVrtDist(wrkVrtSet.get(), foundV1, foundV2, checkedDst);
      }
   }

   /*
   // Try to improve vertices with big Chi2 if something went wrong. Just precaution.
   for (int iv = 0; iv < (int)wrkVrtSet->size(); iv++)
   {
      if (!(*wrkVrtSet)[iv].Good)
         continue; // don't work on vertex which is already bad
      if ((*wrkVrtSet)[iv].selTrk.size() < 3)
         continue;
      double tmpProb = TMath::Prob((*wrkVrtSet)[iv].chi2, 2 * (*wrkVrtSet)[iv].selTrk.size() - 3); // Chi2 of the original vertex
      if (tmpProb < m_globVrtProbCut)
      {
         ATH_MSG_DEBUG("BAD vertex found prob=" << tmpProb);
         tmpProb = improveVertexChi2((*wrkVrtSet)[iv], inParticles->listSelTracks, *state, false);
         (*wrkVrtSet)[iv].projectedVrt = MomProjDist((*wrkVrtSet)[iv].vertex, primaryVx, (*wrkVrtSet)[iv].vertexMom);
      }
   }
   */

   // Final vertex selection/cleaning
   double signif3D = 0., signif2D = 0.;

   //-----  Vertices with >1 tracks
   for (int iv = 0; iv < (int)wrkVrtSet->size(); iv++)
   {
      WrkVrt &curVrt = (*wrkVrtSet)[iv];
      nth = (*wrkVrtSet)[iv].selTrk.size();
      if (nth == 1)
         continue; // 1track vertices for later...
      if (!curVrt.Good)
         continue; // don't work on vertex which is already bad
      (*wrkVrtSet)[iv].Good = false;
      if (nth < 1)
         continue;
      if ((*wrkVrtSet)[iv].projectedVrt < 0.)
         continue;
      if (TMath::Prob(curVrt.chi2, 2 * nth - 3) < m_globVrtProbCut)
         continue;
      //-----------------------------------------------------------------------------------------
      signif3D = vrtVrtDist(myPV, curVrt.vertex);

      //---
      if (signif3D < m_selVrtSigCut)
         continue; // Main PV-SV distance quality cut

      double vertexRadius = sqrt(curVrt.vertex.position()[0] * curVrt.vertex.position()[0] + curVrt.vertex.position()[1] * curVrt.vertex.position()[1]);
      if (vertexRadius > m_maxSVRadiusCut)
         continue; // Too far from interaction point
      curVrt.Good = true;
   }

   printWrkSet(wrkVrtSet.get(), "Final Vertices");

   // Preparing the output RVector
   int nNtrVrt = 0;
   for (auto &iv : (*wrkVrtSet))
   {
      nth = iv.selTrk.size();
      if (iv.Good && nth > 0)
      {
         std::pair<Acts::Vertex<Acts::BoundTrackParameters>, std::deque<long int>> myPair;
         myPair.first = iv.vertex;
         myPair.second = iv.selTrk;
         finalVertices.push_back(myPair);
      }
   }

   //  Saving of results
   return finalVertices;
}
