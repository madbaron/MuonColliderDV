#include "MyVertexFinder.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include "IMPL/TrackImpl.h"
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCRelationNavigator.h>

#include "VertexFitterSimple.h"
#include "VertexFitterActs.h"
#include "VertexingUtils.h"

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "Acts/MagneticField/ConstantBField.hpp"

#include "boost/graph/bron_kerbosch_all_cliques.hpp"

#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

MyVertexFinder aMyVertexFinder;

MyVertexFinder::MyVertexFinder() : Processor("MyVertexFinder")
{

    // Modify processor description
    _description = "MyVertexFinder finds displaced vertices";

    // Input collections
    registerProcessorParameter("TrackCollectionName",
                               "Name of the Track input collection",
                               m_inputTrackCollection,
                               std::string("InputCollection"));

    registerInputCollection(LCIO::MCPARTICLE,
                            "MCParticleCollectionName",
                            "Name of the MCParticle input collection",
                            m_inColMCP,
                            std::string("MCParticle"));

    registerInputCollections(LCIO::LCRELATION,
                             "Particle2TrackRelationName",
                             "Map from MC particle to reconstructed track.",
                             m_inColMC2T,
                             {});

    // Output collection
    registerProcessorParameter("PVCollectionName",
                               "Reconstructed primary vertices",
                               m_outputPVCollection,
                               std::string("VertexCollection"));

    // Output collection
    registerProcessorParameter("SVCollectionName",
                               "Reconstructed displaced vertices",
                               m_outputSVCollection,
                               std::string("VertexCollection"));

    registerOutputCollection(LCIO::LCRELATION,
                             "Track2VertexRelationName",
                             "Map from Vertex to associated tracks.",
                             m_outputTrack2VertexCollection,
                             std::string("Track2VertexRelationName"));
}

void MyVertexFinder::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    m_compatibilityGraph = std::make_unique<boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>>();

    _nRun = 0;
    _nEvt = 0;

    m_trkSigCut = 2.0;
    m_trkMind0 = 0.5;
    m_dRdZRatioCut = 0.25;
    m_sel2VrtProbCut = 0.02;
    m_maxSVRadiusCut = 400.;
    m_cosSVPVCut = 0.;
    m_vertexMergeCut = 4.;
    m_selVrtSigCut = 3.0;
    m_globVrtProbCut = 0.005;
    m_deltaRadiusInnermostHits2VrtCut = 50.;

    m_massPi = 0.139570; // GeV
    m_multiWithOneTrkVrt = false;

    // Get the magnetic field
    dd4hep::Detector &lcdd = dd4hep::Detector::getInstance();
    const double position[3] = {0, 0, 0};                      // position to calculate magnetic field at (the origin in this case)
    double magneticFieldVector[3] = {0, 0, 0};                 // initialise object to hold magnetic field
    lcdd.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep

    // Build ACTS representation of field
    // Note:
    //  magneticFieldVector[2] = 3.57e-13
    //  dd4hep::tesla = 1e-13
    //  Acts::UnitConstants::T = 0.000299792
    m_magneticField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(
        magneticFieldVector[0] / dd4hep::tesla * Acts::UnitConstants::T,
        magneticFieldVector[1] / dd4hep::tesla * Acts::UnitConstants::T,
        magneticFieldVector[2] / dd4hep::tesla * Acts::UnitConstants::T));
}

void MyVertexFinder::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void MyVertexFinder::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the collection of tracker hits
    LCCollection *trackCollection = 0;
    getCollection(trackCollection, m_inputTrackCollection, evt);

    // Load relations
    std::vector<std::shared_ptr<LCRelationNavigator>> trk2mcs;
    for (const std::string &name : m_inColMC2T)
    {
        // Get the collection of tracker hit relations
        LCCollection *MCtoTrackRelationCollection = 0;
        getCollection(MCtoTrackRelationCollection, name, evt);
        std::shared_ptr<LCRelationNavigator> trk2mc = std::make_shared<LCRelationNavigator>(MCtoTrackRelationCollection);
        trk2mcs.push_back(trk2mc);
    }

    // Load MC particles
    LCCollection *particleCollection = 0;
    getCollection(particleCollection, m_inColMCP, evt);

    // Make the output collections
    LCCollection *pvCollection = new LCCollectionVec("Vertex");
    LCCollection *svCollection = new LCCollectionVec("Vertex");
    LCRelationNavigator relV2T(LCIO::VERTEX, LCIO::TRACK);

    size_t nTracks = trackCollection->getNumberOfElements();
    streamlog_out(DEBUG2) << "Event with " << nTracks << " tracks" << std::endl;

    // Make vector of tracks
    ROOT::VecOps::RVec<lcio::TrackImpl> allVec;
    ROOT::VecOps::RVec<int> track_mother;
    ROOT::VecOps::RVec<float> track_innermostRadius;
    for (size_t itTrk = 0; itTrk < nTracks; itTrk++)
    {
        lcio::TrackImpl *trk = static_cast<lcio::TrackImpl *>(trackCollection->getElementAt(itTrk));

        int origin_PDGid = 0;
        // Find the related particle
        EVENT::MCParticle *mcp = nullptr;
        for (std::shared_ptr<LCRelationNavigator> trk2mc : trk2mcs)
        {
            const LCObjectVec &relMCP = trk2mc->getRelatedFromObjects(trk);
            if (!relMCP.empty())
            { // Found the related MC particle
                mcp = dynamic_cast<MCParticle *>(relMCP.at(0));

                // Look for sbottom mothers
                EVENT::MCParticle *mc_mother = nullptr;
                MCParticleVec momVec = mcp->getParents();
                while (!momVec.empty() && fabs(origin_PDGid) != 1000005)
                {
                    mc_mother = dynamic_cast<MCParticle *>(momVec.at(0));
                    origin_PDGid = mc_mother->getPDG();
                    momVec = mc_mother->getParents();
                }
            }
        }
        allVec.push_back(*trk);
        track_mother.push_back(origin_PDGid);

        // Find position of innermost hit
        float innerHitRadius = 99999999999.;
        const EVENT::TrackerHitVec &hits = trk->getTrackerHits();
        for (const EVENT::TrackerHit *hit : hits)
        {
            float hit_R = sqrt(hit->getPosition()[0] * hit->getPosition()[0] + hit->getPosition()[1] * hit->getPosition()[1]);
            if (hit_R < innerHitRadius)
            {
                innerHitRadius = hit_R;
            }
        }

        track_innermostRadius.push_back(innerHitRadius);
    }

    // Find PVs
    ROOT::VecOps::RVec<Acts::Vertex<Acts::BoundTrackParameters>> listPV = IterativeVertexFinder(allVec);

    if (listPV.size() > 0)
    {
        // Find highest sumPT2 vertex
        int i_highest = 0;
        double maxSumpT = 0.;
        for (int i = 0; i < listPV.size(); i++)
        {
            double sumpT = 0.;
            for (const auto &trk : listPV[i].tracks())
            {
                sumpT += trk.fittedParams.transverseMomentum() * trk.fittedParams.transverseMomentum();
            }

            if (sumpT > maxSumpT)
            {
                maxSumpT = sumpT;
                i_highest = i;
            }
        }

        // Make the LCIO vertex
        lcio::VertexImpl *myPV = new VertexImpl(); // Empty Vertex
        myPV->setChi2(listPV[i_highest].fitQuality().first / listPV[i_highest].fitQuality().second);

        float vpos[3] = {listPV[i_highest].position()[i_highest], listPV[i_highest].position()[1], listPV[i_highest].position()[2]};
        myPV->setPosition(vpos);

        auto vtxCov = listPV[i_highest].covariance();
        std::vector<float> covMatrix; // covMat in lcio is a LOWER-triangle matrix.
        covMatrix.push_back(vtxCov(0, 0));
        covMatrix.push_back(vtxCov(1, 0));
        covMatrix.push_back(vtxCov(1, 1));
        covMatrix.push_back(vtxCov(2, 0));
        covMatrix.push_back(vtxCov(2, 1));
        covMatrix.push_back(vtxCov(2, 2));

        myPV->setAlgorithmType("3");
        myPV->setCovMatrix(covMatrix);
        myPV->setPrimary(1);

        pvCollection->addElement(myPV);
        streamlog_out(DEBUG2) << "Found PV (" << listPV[i_highest].position()[0] << ", " << listPV[i_highest].position()[1] << ", " << listPV[i_highest].position()[2] << ") nTracks " << listPV[i_highest].tracks().size() << std::endl;

        ROOT::VecOps::RVec<std::pair<Acts::Vertex<Acts::BoundTrackParameters>, std::deque<long int>>> listDV = getAllDV(allVec, track_innermostRadius, track_mother, listPV[i_highest]);

        m_compatibilityGraph->clear();

        // Loop over vertices and add them to the output
        for (size_t itDV = 0; itDV < listDV.size(); itDV++)
        {
            streamlog_out(DEBUG2) << "Found SV (" << listDV[itDV].first.position()[0] << ", " << listDV[itDV].first.position()[1] << ", " << listDV[itDV].first.position()[2] << ") nTracks " << listDV[itDV].first.tracks().size() << std::endl;

            // Make the LCIO vertex
            lcio::VertexImpl *myDV = new VertexImpl(); // Empty Vertex
            myDV->setChi2(listDV[itDV].first.fitQuality().first / listDV[itDV].first.fitQuality().second);

            float vpos[3] = {listDV[itDV].first.position()[0], listDV[itDV].first.position()[1], listDV[itDV].first.position()[2]};
            myDV->setPosition(vpos);

            auto vtxCov = listDV[itDV].first.covariance();
            std::vector<float> covMatrix; // covMat in lcio is a LOWER-triangle matrix.
            covMatrix.push_back(vtxCov(0, 0));
            covMatrix.push_back(vtxCov(1, 0));
            covMatrix.push_back(vtxCov(1, 1));
            covMatrix.push_back(vtxCov(2, 0));
            covMatrix.push_back(vtxCov(2, 1));
            covMatrix.push_back(vtxCov(2, 2));

            myDV->setAlgorithmType("3");
            myDV->setCovMatrix(covMatrix);
            myDV->setPrimary(false);

            svCollection->addElement(myDV);

            for (int jt = 0; jt < listDV[itDV].second.size(); jt++)
            {
                EVENT::Track *track = static_cast<EVENT::Track *>(trackCollection->getElementAt(listDV[itDV].second[jt]));
                relV2T.addRelation(myDV, track);
            }
        }
    }

    LCCollectionVec *outColT2V = (LCCollectionVec *)relV2T.createLCCollection();

    // Store the filtered hit collections
    evt->addCollection(pvCollection, m_outputPVCollection);
    evt->addCollection(svCollection, m_outputSVCollection);
    evt->addCollection(outColT2V, m_outputTrack2VertexCollection);

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG4) << "   done processing event: " << evt->getEventNumber()
                          << "   in run:  " << evt->getRunNumber() << std::endl;

    _nEvt++;
}

void MyVertexFinder::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void MyVertexFinder::end()
{

    //   std::cout << "MyVertexFinder::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void MyVertexFinder::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
{
    try
    {
        collection = evt->getCollection(collectionName);
    }
    catch (DataNotAvailableException &e)
    {
        streamlog_out(DEBUG5) << "- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
        return;
    }
    return;
}
