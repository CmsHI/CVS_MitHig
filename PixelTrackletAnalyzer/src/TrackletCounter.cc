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
// $Id: TrackletCounter.cc,v 1.1 2008/09/30 13:56:05 yilmaz Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <string>

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

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "TNtuple.h"

using namespace std;
using namespace reco;
using namespace edm;

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

  string betafile;
  string vertexSrc_;

  TrackletFinder* finder_;
  TrackletCorrections* corrections_;
  const TrackerGeometry* trGeo;

  TNtuple* nt;

  edm::Service<TFileService> fs;
  edm::ESHandle < ParticleDataTable > pdt;

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
  betafile = iConfig.getParameter<string>("correctionFile");
  vertexSrc_ = iConfig.getUntrackedParameter<string>("vertexSrc","pixelVertices");

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

   int nbins = corrections_->size();

   vector<double> nchgd;
   vector<double> nlep;
   vector<double> nneut;
   vector<double> ntrc;
   vector<double> ntru;

   nchgd.reserve(nbins);
   nlep.reserve(nbins);
   nneut.reserve(nbins);
   ntrc.reserve(nbins);
   ntru.reserve(nbins);

   for(int ib = 0; ib <nbins; ++ib){
     nchgd[ib] = 0;
     nlep[ib] = 0;
     nneut[ib] = 0;
     ntrc[ib] = 0;
     ntru[ib] = 0;
   }

   // Get reconstructed vertices
   const reco::VertexCollection * recoVertices;
   edm::Handle<reco::VertexCollection> vertexCollection;
   iEvent.getByLabel(vertexSrc_,vertexCollection);
   recoVertices = vertexCollection.product();
   math::XYZVector vertex(0,0,0);
   int greatestvtx = 0; 
   for (unsigned int iv = 0 ; iv< recoVertices->size(); ++iv){
     int daughter = (*recoVertices)[iv].tracksSize();
     if( daughter > (*recoVertices)[greatestvtx].tracksSize())
       {
	 greatestvtx = iv;
       }
   }
   if(recoVertices->size()>0)
     {
       vertex = math::XYZVector((*recoVertices)[greatestvtx].position().x(),
				(*recoVertices)[greatestvtx].position().y(),
				(*recoVertices)[greatestvtx].position().z());
     }

   double z = vertex.z();
   
   finder_->setEvent(iEvent);
   finder_->setVertex(vertex);
   finder_->sortLayers();
   vector<Tracklet> tracklets = finder_->getTracklets();
   //   cout<<"Number of tracklets : "<<tracklets.size()<<endl;
   double nhits = finder_->getNHits();

   // Apply DeltaR cut and beta correction
   for(int i = 0; i < tracklets.size(); ++i){
     double eta = tracklets[i].eta1();
     int bin = corrections_->findBin(nhits,eta,z);
     if(tracklets[i].dR()>corrections_->getDeltaRCut()) continue;
     ntru[bin] += 1;
     ntrc[bin] += 1-corrections_->beta(bin);
   }

   //Read MC Info

   Handle<HepMCProduct> mc;
   iEvent.getByLabel("source",mc);
   const HepMC::GenEvent * evt = mc->GetEvent();
   
   int all = evt->particles_size();
   HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
   HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
   for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
     if((*it)->status() == 1){
       float pdg_id = (*it)->pdg_id();
       float eta = (*it)->momentum().eta();
       int bin = corrections_->findBin(nhits,eta,z);
       float pt = (*it)->momentum().perp();
       const ParticleData * part = pdt->particle(pdg_id );
       float charge = part->charge();
       if(charge == 0 )nneut[bin] += 1;
       else{
	 nchgd[bin] += 1;
	 int id = fabs(pdg_id);
	 if (id==11 || id ==13 || id==15 ){
	   nchgd[bin] -= 1;
	   nlep[bin] += 1;
	 }
       }
     }
   }

   for(int j = 0; j <nbins; ++j){
   nt->Fill(j,nchgd[j],nlep[j],nneut[j],ntrc[j],ntru[j]);
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackletCounter::beginJob(const edm::EventSetup& iSetup)
{
  
  TFile* infile = new TFile(betafile.data(),"read");
  corrections_  = new TrackletCorrections(infile);

  edm::ESHandle<TrackerGeometry> tGeo;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
  trGeo = tGeo.product();


  finder_ = new TrackletFinder(corrections_,trGeo,true);

  nt = fs->make<TNtuple>("nt","Counts in Events","bin:chgd:lep:neut:trc:tru");
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackletCounter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackletCounter);
