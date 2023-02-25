#include "MyMCParticleFilter.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

MyMCParticleFilter aMyMCParticleFilter;

MyMCParticleFilter::MyMCParticleFilter() : Processor("MyMCParticleFilter")
{

    // Modify processor description
    _description = "Select MCParticles based on pdgID";

    // Input collection
    registerProcessorParameter("MCParticleCollectionName",
                               "Name of MCParticle input collection",
                               m_inputMCParticleCollection,
                               std::string("MCParticles"));

    registerProcessorParameter("MCPdgId",
                               "PdgID for filtering",
                               m_pdgId,
                               m_pdgId);

    // Output collection
    registerProcessorParameter("SelectedMCParticlesCollectionName",
                               "Name of the selected MCParticles output collection",
                               m_outputMCParticleCollection,
                               std::string("SelectedMCParticles"));
}

void MyMCParticleFilter::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void MyMCParticleFilter::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void MyMCParticleFilter::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the collection of MCParticles
    LCCollection *MCParticleCollection = 0;
    getCollection(MCParticleCollection, m_inputMCParticleCollection, evt);

    // Make the output collections
    LCCollectionVec *SelectedMCParticlesCollection = new LCCollectionVec(MCParticleCollection->getTypeName());
    SelectedMCParticlesCollection->setSubset(true);

    size_t nMCParticles = MCParticleCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "Total MCParticles: " << nMCParticles << std::endl;

    // Loop over MCParticle
    for (size_t itMCParticle = 0; itMCParticle < nMCParticles; itMCParticle++)
    {
        // Get the MCParticle
        EVENT::MCParticle *MCParticle = static_cast<EVENT::MCParticle *>(MCParticleCollection->getElementAt(itMCParticle));
        if (abs(MCParticle->getPDG()) == m_pdgId)
        {
            bool good = false;

            // Look for daughters with the same pdgid
            EVENT::MCParticle *mc_daughter = nullptr;
            MCParticleVec dauVec = MCParticle->getDaughters();
            for (size_t itDaughter = 0; itDaughter < dauVec.size(); itDaughter++)
            {
                mc_daughter = static_cast<EVENT::MCParticle *>(dauVec.at(itDaughter));
                if (abs(mc_daughter->getPDG()) == 1000022)
                {
                    good = true;
                }
            }

            if (good)
            {
                SelectedMCParticlesCollection->addElement(MCParticle);
                streamlog_out(DEBUG0) << "Found candidate at: " << MCParticle->getEndpoint()[0] << " " << MCParticle->getEndpoint()[1] << " " << MCParticle->getEndpoint()[2] << std::endl;
            }
        }
    }

    streamlog_out(DEBUG4) << "Selected MCParticles: " << SelectedMCParticlesCollection->getNumberOfElements() << std::endl;

    // Store the filtered hit collections
    evt->addCollection(SelectedMCParticlesCollection, m_outputMCParticleCollection);

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG4) << "   done processing event: " << evt->getEventNumber()
                          << "   in run:  " << evt->getRunNumber() << std::endl;

    _nEvt++;
}

void MyMCParticleFilter::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void MyMCParticleFilter::end()
{

    //   std::cout << "MCParticleFilter::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void MyMCParticleFilter::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
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
