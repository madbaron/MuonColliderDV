#include "HitSelectorLayer.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>

#include "TMath.h"
#include "TVector3.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

HitSelectorLayer aHitSelectorLayer;

HitSelectorLayer::HitSelectorLayer() : Processor("HitSelectorLayer")
{

    // Modify processor description
    _description = "HitSelectorLayer applies space selections to reduce the BIB";

    // Input collection
    registerProcessorParameter("TrackerHitCollectionName",
                               "Name of the TrackerHit input collection",
                               m_inputHitCollection,
                               std::string("VertexBarrelCollection"));

    // Output collection
    registerProcessorParameter("GoodHitCollection",
                               "Good hits from tracker",
                               m_outputHitCollection,
                               std::string("VertexBarrelGoodCollection"));
}

void HitSelectorLayer::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void HitSelectorLayer::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void HitSelectorLayer::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the collection of tracker hits
    LCCollection *trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputHitCollection, evt);

    std::string encoderString = trackerHitCollection->getParameters().getStringVal("CellIDEncoding");
    UTIL::CellIDDecoder<TrackerHitPlane> myCellIDEncoding(encoderString);

    // Make the output collections
    LCCollectionVec *GoodHitsCollection = new LCCollectionVec(trackerHitCollection->getTypeName());
    GoodHitsCollection->setSubset(true);
    GoodHitsCollection->parameters().setValue("CellIDEncoding", encoderString);

    int nHits = trackerHitCollection->getNumberOfElements();

    // Loop over tracker hits
    for (int itHit = 0; itHit < nHits; itHit++)
    {

        // Get the hit
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(itHit));
        unsigned int layer = myCellIDEncoding(hit)["layer"];
        unsigned int detector = myCellIDEncoding(hit)["system"];
        int side = myCellIDEncoding(hit)["side"];

        if (detector == 2)
        {
            if (layer > 3)
            {
                streamlog_out(DEBUG0) << " --> rejected EndCap hit. Layer = " << layer << " " << detector << " " << side << std::endl;
                continue;
            }
        }

        if (detector == 3)
        {
            if (layer > 0)
            {
                streamlog_out(DEBUG0) << " --> rejected IT hit. Layer = " << layer << " " << detector << " " << side << std::endl;
                continue;
            }
        }

        // Add hits to output
        GoodHitsCollection->addElement(hit);
    }

    // Store the filtered hit collections
    evt->addCollection(GoodHitsCollection, m_outputHitCollection);

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG) << "   done processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;
    
    _nEvt++;
}

void HitSelectorLayer::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void HitSelectorLayer::end()
{

    //   std::cout << "HitSelectorLayer::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void HitSelectorLayer::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
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
