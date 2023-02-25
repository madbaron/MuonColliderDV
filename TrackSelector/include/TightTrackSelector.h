#ifndef TightTrackSelector_h
#define TightTrackSelector_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include <TH1F.h>

#include <EVENT/LCCollection.h>

using namespace lcio;
using namespace marlin;

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
 * @param TrackCollectionName Name of the input track collection
 * @param SelectedTrackCollection Base name of the output track collections
 *
 * @author F. Meloni, DESY
 * @version $Id: TightTrackSelector.h,v 0.1 2022-05-19 12:57:39 fmeloni Exp $
 */

class TightTrackSelector : public Processor
{

public:
  virtual Processor *newProcessor() { return new TightTrackSelector; }

  TightTrackSelector();

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
  std::string m_inputTrackCollection = "";
  std::string m_outputTrackCollection = "";
  std::string m_inputHitCollection = "";
  std::string m_outputHitsCollection = "";
  std::string m_outputUsedSimHitsCollection = "";
  std::string m_outputUsedHitRelCollection = "";
  std::vector<std::string> m_inputTrackerHitRelNames{};

  int _nRun{};
  int _nEvt{};

  double m_magneticField{};

  // --- Processor parameters:
  bool m_fillHistos{};

  // --- Diagnostic histograms:
  TH1F *m_d0 = nullptr;
  TH1F *m_z0 = nullptr;
  TH1F *m_pT = nullptr;
  TH1F *m_theta = nullptr;
  TH1F *m_nhits = nullptr;
  TH1F *m_chi2ndf = nullptr;
};

#endif
