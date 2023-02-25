#ifndef MyMCParticleFilter_h
#define MyMCParticleFilter_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <EVENT/LCCollection.h>


using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param MCParticleCollectionName Name of the input MCParticle collection
 * @param SelectedMCParticleCollection Base name of the output MCParticle collections
 * 
 * @author F. Meloni, DESY
 * @version $Id: MyMCParticleFilter.h,v 0.1 2022-05-19 12:57:39 fmeloni Exp $ 
 */

class MyMCParticleFilter : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new MyMCParticleFilter ; }
  
  
  MyMCParticleFilter() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  // Call to get collections
  void getCollection(LCCollection *&, std::string, LCEvent *);

  
protected:
  // Collection names for (in/out)put
  std::string m_inputMCParticleCollection = "";
  int32_t m_pdgId=13; 
  std::string m_outputMCParticleCollection = "";

  int _nRun{} ;
  int _nEvt{} ;
} ;

#endif



