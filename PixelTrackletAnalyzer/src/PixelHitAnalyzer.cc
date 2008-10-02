// -*- C++ -*-
//
// Package:    PixelHitAnalyzer
// Class:      PixelHitAnalyzer
// 
/**\class PixelHitAnalyzer PixelHitAnalyzer.cc MitHig/PixelHitAnalyzer/src/PixelHitAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yilmaz Yetkin
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id: PixelHitAnalyzer.cc,v 1.1 2008/09/30 13:56:05 yilmaz Exp $
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TNtuple.h"

using namespace std;

//
// class decleration
//

struct PixelEvent{
  //  int nhits;
  double z;
};


class PixelHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PixelHitAnalyzer(const edm::ParameterSet&);
      ~PixelHitAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  const char* betafile;
  TrackletFinder* finder_;
  TrackletCorrections* corrections_;
  const TrackerGeometry* trGeo;
  edm::Service<TFileService> fs;           

  TTree* pixelTree_;
  TBranch* branch0_;
  TBranch* branch1_;
  TBranch* branch2_;

  TNtuple* hits1_;
  TNtuple* hits2_;

  PixelEvent evt_;

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
PixelHitAnalyzer::PixelHitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


PixelHitAnalyzer::~PixelHitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PixelHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   hits1_->Reset();
   hits2_->Reset();

   finder_->setEvent(iEvent);

   //   vector<Tracklet> tracklets = finder_->getTracklets();

   finder_->fillPixelEvent(hits1_,hits2_);

   pixelTree_->Fill();

   //   cout<<"Number of tracklets : "<<tracklets.size()<<endl;


}


// ------------ method called once each job just before starting event loop  ------------
void 
PixelHitAnalyzer::beginJob(const edm::EventSetup& iSetup)
{
  
  //  TFile* infile = new TFile(betafile,"read");
  //  corrections_  = new TrackletCorrections(infile);
  corrections_  = new TrackletCorrections(1,1,1);

  edm::ESHandle<TrackerGeometry> tGeo;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
  trGeo = tGeo.product();

  hits1_ = new TNtuple("layer1","Layer 1 Hits","x:y:z");
  hits2_ = new TNtuple("layer1","Layer 2 Hits","x:y:z");

  finder_ = new TrackletFinder(corrections_,trGeo,true);

  pixelTree_ = fs->make<TTree>("PixelTree","Tree of Pixel Hits");
  branch1_ = pixelTree_->Branch("Layer1","TNtuple",&hits1_,64000,1);
  branch1_ = pixelTree_->Branch("Layer2","TNtuple",&hits2_,64000,1);
  //  branch0_ = pixelTree_->Branch("Event",&evt_,"nhits/I:z/F");
  branch0_ = pixelTree_->Branch("Event",&evt_,"z/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelHitAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelHitAnalyzer);
