#include "TightTrackSelector.h"

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>

#include <IMPL/TrackImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/SimTrackerHitImpl.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>

#include <marlin/AIDAProcessor.h>
#include <marlinutil/GeometryUtil.h>

#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

TightTrackSelector aTightTrackSelector;

TightTrackSelector::TightTrackSelector() : Processor("TightTrackSelector")
{

    // Modify processor description
    _description = "Select good tracks";

    // --- Processor parameters:
    registerProcessorParameter("TrackCollectionName",
                               "Name of reconstructed track input collection",
                               m_inputTrackCollection,
                               std::string("Tracks"));

    registerProcessorParameter("HitsCollectionName",
                               "Name of hit input collection",
                               m_inputHitCollection,
                               std::string("Hits"));

    registerProcessorParameter("TrackerHitInputRelations",
                               "Name of the tracker hit relation collections",
                               m_inputTrackerHitRelNames,
                               {});

    registerProcessorParameter("SelectedTracksCollectionName",
                               "Name of the selected tracks output collection",
                               m_outputTrackCollection,
                               std::string("SelectedTracks"));

    registerProcessorParameter("UsedHitsCollectionName",
                               "Name of the used hits output collection",
                               m_outputHitsCollection,
                               std::string("UsedHits"));

    registerProcessorParameter("UsedSimHitOutputCollections",
                               "Name of the used simhit output collections",
                               m_outputUsedSimHitsCollection,
                               std::string("UsedSimHits"));

    registerProcessorParameter("UsedHitOutputRelations",
                               "Name of the used hit relation collections",
                               m_outputUsedHitRelCollection,
                               std::string("UsedHitsRelations"));

    registerProcessorParameter("FillHistograms",
                               "Flag to fill the diagnostic histograms",
                               m_fillHistos,
                               false);
}

void TightTrackSelector::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;

    // --- Get the value of the magnetic field
    m_magneticField = MarlinUtil::getBzAtOrigin();

    // --- Initialize the AIDAProcessor and book the diagnostic histograms:
    AIDAProcessor::histogramFactory(this);

    m_d0 = new TH1F("d0", "track d0 [mm]", 200, -150., 150.);
    m_z0 = new TH1F("z0", "track z0 [mm]", 200, -150., 150.);
    m_pT = new TH1F("pT", "track pT [GeV]", 100, 0., 100.);
    m_theta = new TH1F("theta", "track theta", 100, 0., M_PI);
    m_nhits = new TH1F("nhits", "track nhits", 20, 0., 20.);
    m_chi2ndf = new TH1F("chi2ndf", "track chi2ndf", 100, 0., 10.);
}

void TightTrackSelector::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void TightTrackSelector::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the input collections
    LCCollection *trackCollection = 0;
    getCollection(trackCollection, m_inputTrackCollection, evt);

    LCCollection *trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputHitCollection, evt);
    std::string encoderString = trackerHitCollection->getParameters().getStringVal("CellIDEncoding");
    UTIL::CellIDDecoder<TrackerHitPlane> myCellIDEncoding(encoderString);

    // --- Get the input hit relation collections:
    const unsigned int nTrackerHitCol = m_inputTrackerHitRelNames.size();
    std::vector<std::shared_ptr<LCRelationNavigator>> hit2simhits;
    for (unsigned int icol = 0; icol < nTrackerHitCol; ++icol)
    {
        // Get the collection of tracker hit relations
        LCCollection *trackerHitRelationCollection = 0;
        getCollection(trackerHitRelationCollection, m_inputTrackerHitRelNames[icol], evt);
        std::shared_ptr<LCRelationNavigator> hit2simhit = std::make_shared<LCRelationNavigator>(trackerHitRelationCollection);
        hit2simhits.push_back(hit2simhit);
    }

    // Make the output collections
    LCCollectionVec *SelectedTracksCollection = new LCCollectionVec(LCIO::TRACK);
    // Enable the track collection to point back to hits
    LCFlagImpl trkFlag(0);
    trkFlag.setBit(LCIO::TRBIT_HITS);
    SelectedTracksCollection->setFlag(trkFlag.getFlag());

    // Clone encoding and flags from other input collections
    LCCollectionVec *UsedHitsCollection = new LCCollectionVec(trackerHitCollection->getTypeName());
    UsedHitsCollection->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *UsedSimHitsCollection = new LCCollectionVec(LCIO::SIMTRACKERHIT);
    UsedSimHitsCollection->parameters().setValue("CellIDEncoding", encoderString);
    // LCFlagImpl lcFlag_sim(inputSimHitColls[0]->getFlag());
    // UsedSimHitsCollection->setFlag(lcFlag_sim.getFlag());
    LCCollectionVec *UsedHitsRelationCollection = new LCCollectionVec(LCIO::LCRELATION);
    // LCFlagImpl lcFlag_rel(hit2simhits[0]->getFlag());
    // UsedHitsRelationCollection->setFlag(lcFlag_rel.getFlag());

    size_t nTracks = trackCollection->getNumberOfElements();
    streamlog_out(MESSAGE) << "Total tracks: " << nTracks << std::endl;

    // Loop over tracker hits
    for (size_t itTrack = 0; itTrack < nTracks; itTrack++)
    {
        // Get the track
        EVENT::Track *track = static_cast<EVENT::Track *>(trackCollection->getElementAt(itTrack));

        if (m_fillHistos)
        {
            m_d0->Fill(track->getD0());
            m_z0->Fill(track->getZ0());
            m_pT->Fill(0.3 * m_magneticField / fabs(track->getOmega() * 1000.));
            m_theta->Fill(M_PI / 2. - atan(track->getTanLambda()));
            m_nhits->Fill(track->getTrackerHits().size());
            m_chi2ndf->Fill(track->getChi2() / track->getNdf());
        }

        if (track->getTrackerHits().size() < 6)
        {
            continue;
        }

        if (track->getChi2() / track->getNdf() > 3)
        {
            continue;
        }

        if (fabs(track->getD0()) > 150)
        {
            continue;
        }

        if (fabs(track->getZ0()) > 75)
        {
            continue;
        }

        // Make track object
        IMPL::TrackImpl *new_track = new IMPL::TrackImpl;
        // Fit state
        new_track->setChi2(track->getChi2());
        new_track->setNdf(track->getNdf());
        // Track states
        const EVENT::TrackStateVec &trackstates = track->getTrackStates();
        for (const EVENT::TrackState *state : trackstates)
        {
            IMPL::TrackStateImpl *trackState = new IMPL::TrackStateImpl();

            // Basic properties
            trackState->setLocation(state->getLocation());

            trackState->setPhi(state->getPhi());
            trackState->setTanLambda(state->getTanLambda());
            trackState->setOmega(state->getOmega());
            trackState->setD0(state->getD0());
            trackState->setZ0(state->getZ0());

            trackState->setCovMatrix(state->getCovMatrix());

            new_track->trackStates().push_back(trackState);
        }

        // Hits
        const EVENT::TrackerHitVec &hits = track->getTrackerHits();
        for (EVENT::TrackerHit *hiti : hits)
        {
            TrackerHitPlane *hit = dynamic_cast<TrackerHitPlane *>(hiti);
            IMPL::TrackerHitPlaneImpl *hit_new = new IMPL::TrackerHitPlaneImpl();

            hit_new->setCellID0(hit->getCellID0());
            hit_new->setCellID1(hit->getCellID1());
            hit_new->setType(hit->getType());
            hit_new->setPosition(hit->getPosition());
            hit_new->setU(hit->getU());
            hit_new->setV(hit->getV());
            hit_new->setdU(hit->getdU());
            hit_new->setdV(hit->getdV());
            hit_new->setEDep(hit->getEDep());
            hit_new->setEDepError(hit->getEDepError());
            hit_new->setTime(hit->getTime());
            hit_new->setQuality(hit->getQuality());

            UsedHitsCollection->addElement(hit_new);

            // Find the sim hit
            for (std::shared_ptr<LCRelationNavigator> hit2simhit : hit2simhits)
            {
                const LCObjectVec &simHitVector = hit2simhit->getRelatedToObjects(hiti);
                const EVENT::FloatVec &rel_weights = hit2simhit->getRelatedToWeights(hiti);
                if (!simHitVector.empty())
                { // Found the sim hit
                    EVENT::SimTrackerHit *simhit = dynamic_cast<SimTrackerHit *>(simHitVector.at(0));
                    float rel_wei = rel_weights.at(0);

                    SimTrackerHitImpl *simhit_new = new SimTrackerHitImpl();

                    simhit_new->setCellID0(simhit->getCellID0());
                    simhit_new->setCellID1(simhit->getCellID1());
                    simhit_new->setPosition(simhit->getPosition());
                    simhit_new->setEDep(simhit->getEDep());
                    simhit_new->setTime(simhit->getTime());
                    simhit_new->setMCParticle(simhit->getMCParticle());
                    simhit_new->setMomentum(simhit->getMomentum());
                    simhit_new->setPathLength(simhit->getPathLength());
                    simhit_new->setQuality(simhit->getQuality());
                    simhit_new->setOverlay(simhit->isOverlay());
                    simhit_new->setProducedBySecondary(simhit->isProducedBySecondary());

                    UsedSimHitsCollection->addElement(simhit_new);

                    LCRelationImpl *rel_new = new LCRelationImpl();

                    rel_new->setFrom(hit_new);
                    rel_new->setTo(simhit_new);
                    rel_new->setWeight(rel_wei);

                    UsedHitsRelationCollection->addElement(rel_new);

                    break;
                }
            }

            new_track->addHit(hit_new);
        }
        // Other stuff
        new_track->setTypeBit(track->getType());
        new_track->setdEdx(track->getdEdx());
        new_track->setdEdxError(track->getdEdxError());
        new_track->setRadiusOfInnermostHit(track->getRadiusOfInnermostHit());

        SelectedTracksCollection->addElement(new_track);
    }

    streamlog_out(MESSAGE) << "Selected tracks: " << SelectedTracksCollection->getNumberOfElements() << std::endl;

    // Store the filtered hit collections
    evt->addCollection(SelectedTracksCollection, m_outputTrackCollection);
    evt->addCollection(UsedHitsCollection, m_outputHitsCollection);
    evt->addCollection(UsedSimHitsCollection, m_outputUsedSimHitsCollection);
    evt->addCollection(UsedHitsRelationCollection, m_outputUsedHitRelCollection);

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG4) << "   done processing event: " << evt->getEventNumber()
                          << "   in run:  " << evt->getRunNumber() << std::endl;

    _nEvt++;
}

void TightTrackSelector::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void TightTrackSelector::end()
{
    streamlog_out(DEBUG8) << "TightTrackSelector::end()  " << name()
                          << " processed " << _nEvt << " events in " << _nRun << " runs "
                          << std::endl;
}

void TightTrackSelector::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
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
