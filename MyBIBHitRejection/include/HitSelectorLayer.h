#ifndef HitSelectorLayer_h
#define HitSelectorLayer_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <map>
#include <vector>

#include <EVENT/LCCollection.h>

using namespace lcio;
using namespace marlin;

/**  Hit filtering processor for marlin.
 * 
 * @param TrackerHitCollectionName Name of the input hit collection
 * @param GoodHitCollection Base name of the output hit collections
 * 
 * @author F. Meloni, DESY; R. Simoniello, CERN
 * @version $Id: HitSelectorLayer.h,v 0.1 2020-09-27 11:24:21 fmeloni Exp $ 
 */

class HitSelectorLayer : public Processor
{

protected:
  static const size_t MAX_NHITS = 10000000;

public:
  virtual Processor *newProcessor() { return new HitSelectorLayer; }

  HitSelectorLayer();

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
  // Collection names for (in/out)put
  std::string m_inputHitCollection = "";
  std::string m_outputHitCollection = "";

  // map hits in the detector by layer
  std::map < unsigned int, std::vector<size_t> > m_hitsMap;

  // hit decisions
  bool m_accepted[MAX_NHITS];

  int _nRun{};
  int _nEvt{};
};

#endif
