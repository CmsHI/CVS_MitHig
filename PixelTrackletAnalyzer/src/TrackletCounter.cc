// -*- C++ -*-
//
// Package:    TrackletCounter
// Class:      TrackletCounter
// 
/**\class TrackletCounter TrackletCounter.cc MitHig/TrackletCounter/src/TrackletCounter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yilmaz Yetkin
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "MitHig/PixelTracklet/interface/TrackletCorrections.h"
#include "MitHig/PixelTrackletAnalyzer/interface/TrackletFinder.h"

using namespace std;

//
// class decleration
//

class TrackletCounter : public edm::EDAnalyzer {
   public:
      explicit TrackletCounter(const edm::ParameterSet&);
      ~TrackletCounter();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  const char* betafile;
  TrackletFinder* finder_;
  TrackletCorrections* corrections_;
  const TrackerGeometry* trGeo;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackletCounter::TrackletCounter(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


TrackletCounter::~TrackletCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrackletCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   finder_->setEvent(iEvent);

   vector<Tracklet> tracklets = finder_->getTracklets();

   cout<<"Number of tracklets : "<<tracklets.size()<<endl;


}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackletCounter::beginJob(const edm::EventSetup& iSetup)
{
  
  TFile* infile = new TFile(betafile,"read");
  corrections_  = new TrackletCorrections(infile);

  edm::ESHandle<TrackerGeometry> tGeo;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
  trGeo = tGeo.product();


  finder_ = new TrackletFinder(corrections_,trGeo,true);
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackletCounter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackletCounter);
