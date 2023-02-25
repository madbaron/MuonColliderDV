#ifndef MyVertexFinder_h
#define MyVertexFinder_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <map>
#include <vector>

#include "TLorentzVector.h"

#include "ROOT/RVec.hxx"

#include "EVENT/Track.h"
#include "EVENT/Vertex.h"
#include <IMPL/TrackImpl.h>

// ACTS
#include "Acts/Vertexing/Vertex.hpp"
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#define BOOST_ALLOW_DEPRECATED_HEADERS
#include "boost/graph/adjacency_list.hpp"

#include <EVENT/LCCollection.h>

#include "VertexingUtils.h"

using namespace lcio;
using namespace marlin;

/**  Hit filtering processor for marlin.
 *
 * @param TrackCollectionName Name of the input track collection
 * @param PVCollection Name of the output PV collection
 * @param SVCollection Name of the output SV collection
 *
 * @author F. Meloni, DESY
 * @version $Id: MyVertexFinder.h,v 0.1 2021-11-04 09:27:21 fmeloni Exp $
 */

class MyVertexFinder : public Processor
{

public:
  virtual Processor *newProcessor() { return new MyVertexFinder; }

  MyVertexFinder();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  // Call to get collections
  void getCollection(LCCollection *&, std::string, LCEvent *);

protected:
  const Acts::MagneticFieldContext &magneticFieldContext() const;
  std::shared_ptr<Acts::MagneticFieldProvider> magneticField() const;
  const Acts::GeometryContext &geometryContext() const;

  // Collection names for (in/out)put
  std::string m_inputTrackCollection = "";
  std::string m_outputPVCollection = "";
  std::string m_outputSVCollection = "";
  std::string m_outputTrack2VertexCollection = "";

  std::string m_inColMCP = "";
  std::vector<std::string> m_inColMC2T;

  /** Output collection names.
   */
  std::string _outColMC2T;
  float m_trkSigCut{};
  float m_trkMind0{};

  float m_dRdZRatioCut{};
  float m_sel2VrtProbCut{};
  float m_maxSVRadiusCut{};
  float m_cosSVPVCut{};
  float m_vertexMergeCut{};
  float m_selVrtSigCut{};
  float m_globVrtProbCut{};
  float m_deltaRadiusInnermostHits2VrtCut{};

  std::unique_ptr<boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>>
      m_compatibilityGraph{};

  int _nRun{};
  int _nEvt{};

private:
  double m_massPi{};
  bool m_multiWithOneTrkVrt{};

  ROOT::VecOps::RVec<Acts::Vertex<Acts::BoundTrackParameters>> IterativeVertexFinder(const ROOT::VecOps::RVec<lcio::TrackImpl> &inParticles) const;

  ROOT::VecOps::RVec<std::pair<Acts::Vertex<Acts::BoundTrackParameters>, std::deque<long int>>> getAllDV(const ROOT::VecOps::RVec<lcio::TrackImpl> &inParticles,
                                                                                                         const ROOT::VecOps::RVec<float> track_innermostRadius,
                                                                                                         const ROOT::VecOps::RVec<int> track_mother,
                                                                                                         const Acts::Vertex<Acts::BoundTrackParameters> &primaryVx) const;
  // void selGoodTrkParticle(workVectorArr *inParticles,
  //                         const VertexingUtils::FCCAnalysesVertex &primaryVx);
  void select2TrVrt(ROOT::VecOps::RVec<Acts::BoundTrackParameters> &allTracks,
                    const ROOT::VecOps::RVec<float> track_innermostRadius,
                    const Acts::Vertex<Acts::BoundTrackParameters> &primVrt) const;

  double projSV_PV(const TVector3 &SV, const TVector3 &PV, const TLorentzVector &Direction) const;
  double MomProjDist(const TVector3 &SecVrt, const TVector3 &primVrt, const TLorentzVector &Mom) const;

  struct WrkVrt
  {
    bool Good = true;
    std::deque<long int> selTrk;
    Acts::Vertex<Acts::BoundTrackParameters> vertex;
    TLorentzVector vertexMom;
    long int vertexCharge{};
    std::vector<double> vertexCov;
    std::vector<double> chi2PerTrk;
    std::vector<std::vector<double>> trkAtVrt;
    double chi2{};
    double projectedVrt = 0.;
    int detachedTrack = -1;
  };

  void printWrkSet(const std::vector<WrkVrt> *WrkSet, const std::string &name) const;
  int nTrkCommon(std::vector<WrkVrt> *WrkVrtSet, int indexV1, int indexV2) const;

  double refitVertex(WrkVrt &Vrt, ROOT::VecOps::RVec<Acts::BoundTrackParameters> &SelectedTracks) const;

  double mergeAndRefitVertices(WrkVrt &v1, WrkVrt &v2, WrkVrt &newvrt,
                               ROOT::VecOps::RVec<Acts::BoundTrackParameters> &AllTrackList) const;

  double refineVerticesWithCommonTracks(WrkVrt &v1, WrkVrt &v2, ROOT::VecOps::RVec<Acts::BoundTrackParameters> &allTrackList) const;
  double minVrtVrtDist(std::vector<WrkVrt> *WrkVrtSet, int &indexV1, int &indexV2, std::vector<double> &check) const;
  double vrtVrtDist(const Acts::Vertex<Acts::BoundTrackParameters> &vrt1, const Acts::Vertex<Acts::BoundTrackParameters> &vrt2) const;

  Acts::MagneticFieldContext m_magneticFieldContext;
  Acts::GeometryContext m_geometryContext;
  std::shared_ptr<Acts::MagneticFieldProvider> m_magneticField;
};

struct clique_visitor
{
  clique_visitor(std::vector<std::vector<int>> &input) : m_allCliques(input) { input.clear(); }

  template <typename Clique, typename Graph>
  void clique(const Clique &clq, Graph &)
  {
    std::vector<int> new_clique(0);
    for (auto i = clq.begin(); i != clq.end(); ++i)
      new_clique.push_back(*i);
    m_allCliques.push_back(new_clique);
  }

  std::vector<std::vector<int>> &m_allCliques;
};

#endif
