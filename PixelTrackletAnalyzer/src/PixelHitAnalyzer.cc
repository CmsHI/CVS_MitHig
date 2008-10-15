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
// $Id: PixelHitAnalyzer.cc,v 1.1 2008/10/02 00:14:59 yilmaz Exp $
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class decleration
//

#define MAXHITS 1000
#define MAXVTX 10

struct PixelEvent{
   int nhits1;
   int nhits2;

   int mult;
   //   int mult2;

   int nv;
   float vz[MAXVTX];

   float eta1[MAXHITS];
   float phi1[MAXHITS];
   float r1[MAXHITS];
   int id1[MAXHITS];
   float cs1[MAXHITS];
   float cg1[MAXHITS];

   float eta2[MAXHITS];
   float phi2[MAXHITS];
   float r2[MAXHITS];
   int id2[MAXHITS];
   float cs2[MAXHITS];
   float cg2[MAXHITS];

};


class PixelHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PixelHitAnalyzer(const edm::ParameterSet&);
      ~PixelHitAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
   void fillVertices(const edm::Event& iEvent);
   void fillHits(const edm::Event& iEvent);
   int associateSimhitToTrackingparticle(unsigned int trid );
   bool checkprimaryparticle(TrackingParticleRef tp);

      // ----------member data ---------------------------

   //  const char* betafile;
   //  TrackletFinder* finder_;
   //  TrackletCorrections* corrections_;

   bool doMC_;
   vector<string> vertexSrc_;
   double etaMult_;

  const TrackerGeometry* geo_;
  edm::Service<TFileService> fs;           
   edm::Handle<TrackingParticleCollection> trackingParticles;

  TTree* pixelTree_;
  PixelEvent pev_;

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
   doMC_             = iConfig.getUntrackedParameter<bool>  ("doMC",true);
   vertexSrc_ = iConfig.getParameter<vector<string> >("vertexSrc");
   etaMult_ = iConfig.getUntrackedParameter<double>  ("nHitsRegion",1.);

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
   pev_.nhits1 = 0;
   pev_.nhits2 = 0;
   pev_.mult = 0;

   pev_.nv = 0;

   fillVertices(iEvent);
   fillHits(iEvent);

   pixelTree_->Fill();

}

void
PixelHitAnalyzer::fillVertices(const edm::Event& iEvent){

   if(doMC_){
      unsigned int daughter = 0;
      unsigned int nVertex = 0;
      unsigned int greatestvtx = -1;
      Handle<TrackingVertexCollection> vertices;
      iEvent.getByLabel("mergedtruth","MergedTrackTruth", vertices);
      nVertex = vertices->size();
      for (unsigned int i = 0 ; i< vertices->size(); ++i){
	 daughter = (*vertices)[i].nDaughterTracks();
	 if( daughter > (*vertices)[greatestvtx].nDaughterTracks()) greatestvtx = i;
      }
      
      if(vertices->size()>0){
	 pev_.vz[pev_.nv] = (*vertices)[greatestvtx].position().z();
      }else{
	 pev_.vz[pev_.nv] =  -99; 
      }
      pev_.nv++;
   }
   
   for(int iv = 0; iv < vertexSrc_.size(); ++iv){
      const reco::VertexCollection * recoVertices;
      edm::Handle<reco::VertexCollection> vertexCollection;
      iEvent.getByLabel(vertexSrc_[iv],vertexCollection);
      recoVertices = vertexCollection.product();
      int daughter = 0;
      int nVertex = 0;
      int greatestvtx = -1;
      
      nVertex = recoVertices->size();
      for (unsigned int i = 0 ; i< recoVertices->size(); ++i){
	 daughter = (*recoVertices)[i].tracksSize();
	 if( daughter > (*recoVertices)[greatestvtx].tracksSize()) greatestvtx = i;
      }
      
      if(recoVertices->size()>0){
	 pev_.vz[pev_.nv] = (*recoVertices)[greatestvtx].position().z();
      }else{
	 pev_.vz[pev_.nv] =  -99;
      }
      pev_.nv++;
   }

}

void
PixelHitAnalyzer::fillHits(const edm::Event& iEvent){

   TrackerHitAssociator theHitAssociator(iEvent);
   if(doMC_)iEvent.getByLabel("mergedtruth","MergedTrackTruth",trackingParticles);
   
   const SiPixelRecHitCollection* rechits;
   Handle<SiPixelRecHitCollection> rchts;
   iEvent.getByLabel("siPixelRecHits",rchts);
   rechits = rchts.product();

   for(SiPixelRecHitCollection::id_iterator id = rechits->id_begin(); id!= rechits->id_end(); id++){
      if((*id).subdetId() == int(PixelSubdetector::PixelBarrel)){
	 PXBDetId pid(*id);
	 SiPixelRecHitCollection::range range;
	 int layer = pid.layer();
	 if(layer == 1 || layer == 2) range = rechits->get(*id);
	 for(SiPixelRecHitCollection::const_iterator recHit = range.first; recHit!= range.second; recHit++){
	    
	    const SiPixelRecHit* recHit1 = &*recHit;

	    // SIM INFO
	    int trid = -9999;

	    if (doMC_) {
	       vector<PSimHit> simHits1 = theHitAssociator.associateHit(*recHit1);
	       const PSimHit * bestSimHit1 = 0;
	       int simIdx =0;

	       //gets the primary simhit and its specifications for the rechit   	     
	       for(vector<PSimHit>::const_iterator simHit1 = simHits1.begin(); simHit1!= simHits1.end(); simHit1++){  
		  simIdx++;
		  unsigned int associatedTPID = associateSimhitToTrackingparticle((&(*simHit1))->trackId());
		  if (associatedTPID == -1) continue;    // doesn't match to any Trackingparticle
		  
		  TrackingParticleRef associatedTP(trackingParticles, associatedTPID);
		  int ptype = (&(*simHit1))->processType();
		  
		  bool isprimary = checkprimaryparticle(associatedTP);
		  
		  if (isprimary && bestSimHit1==0){ 
		     bestSimHit1 = &(*simHit1);
		     break;
		  }
	       } 
	       
	       if(bestSimHit1!=0){
		  trid = bestSimHit1->trackId();  
	       }
	    }
	    
	    // GEOMETRY INFO
	    const PixelGeomDetUnit* pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit1->geographicalId()));
	    GlobalPoint gpos = pixelLayer->toGlobal(recHit1->localPosition());
	    
	    if(layer == 1){ 
	       pev_.eta1[pev_.nhits1] = gpos.eta();
	       pev_.phi1[pev_.nhits1] = gpos.phi();
	       pev_.r1[pev_.nhits1] = gpos.perp();
	       pev_.id1[pev_.nhits1] = trid;
	       pev_.cs1[pev_.nhits1] = recHit1->cluster()->size(); //Cluster Size
               pev_.cg1[pev_.nhits1] = recHit1->cluster()->charge(); //Cluster Charge
	       pev_.nhits1++;
	       if(fabs(gpos.eta()) < etaMult_ ) pev_.mult++;
	    }
	    if(layer == 2){
	       pev_.eta2[pev_.nhits2] = gpos.eta();
	       pev_.phi2[pev_.nhits2] = gpos.phi();
	       pev_.r2[pev_.nhits2] = gpos.perp();
               pev_.id2[pev_.nhits2] = trid;
	       pev_.cs2[pev_.nhits2] = recHit1->cluster()->size(); //Cluster Size
               pev_.cg2[pev_.nhits2] = recHit1->cluster()->charge(); //Cluster Charge
	       pev_.nhits2++;
	    } 
	    
	 }
      }
   }

}

int PixelHitAnalyzer::associateSimhitToTrackingparticle(unsigned int trid )
{
   int ref=-1;

   const TrackingParticleCollection TPCProd = *(trackingParticles.product());
   for (TrackingParticleCollection::size_type i=0; i<TPCProd.size(); i++){
      TrackingParticleRef tp(trackingParticles, i);
      vector <PSimHit> particlesimhits = tp->trackPSimHit();
      for(vector<PSimHit>::const_iterator simhit = particlesimhits.begin(); simhit != particlesimhits.end(); ++simhit)
	 {
	    //cout <<"       matching TP: "<<i<<" TPsimhitid: "<<simhit->trackId()<<" simhitId: "<<trid<<endl;
	    if(simhit->trackId()==trid)
	       {
		  ref=i;
		  break;
	       }
	 }
      if (ref!=-1) break;
   }

   return ref;
}
   

bool PixelHitAnalyzer::checkprimaryparticle(TrackingParticleRef tp)
{
   int primarycheck=2;
   if(((tp->charge()==1)||(tp->charge()==-1))&&(tp->vertex().Rho()<0.2))
      {
	 primarycheck=1;
      } else {
	 primarycheck=0;
      }
   return primarycheck;
}       



// ------------ method called once each job just before starting event loop  ------------
void 
PixelHitAnalyzer::beginJob(const edm::EventSetup& iSetup)
{
  
  //  TFile* infile = new TFile(betafile,"read");
  //  corrections_  = new TrackletCorrections(infile);
   //  corrections_  = new TrackletCorrections(1,1,1);

  edm::ESHandle<TrackerGeometry> tGeo;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
  geo_ = tGeo.product();

  //  finder_ = new TrackletFinder(corrections_,trGeo,true);

  pixelTree_ = fs->make<TTree>("PixelTree","Tree of Pixel Hits");
  pixelTree_->Branch("nhits1",&pev_.nhits1,"nhits1/I");
  pixelTree_->Branch("nhits2",&pev_.nhits2,"nhits2/I");
  pixelTree_->Branch("mult",&pev_.mult,"mult/I");
  //  pixelTree_->Branch("mult2",&pev_.mult2,"mult2/I");
  pixelTree_->Branch("nv",&pev_.nv,"nv/I");
  pixelTree_->Branch("vz",pev_.vz,"vz[nv]/F");
  pixelTree_->Branch("eta1",pev_.eta1,"eta1[nhits1]/F");
  pixelTree_->Branch("phi1",pev_.phi1,"phi1[nhits1]/F");
  pixelTree_->Branch("r1",pev_.r1,"r1[nhits1]/F");
  pixelTree_->Branch("id1",pev_.id1,"id1[nhits1]/I");
  pixelTree_->Branch("cs1",pev_.cs1,"cs1[nhits1]/F");
  pixelTree_->Branch("cg1",pev_.cg1,"cg1[nhits1]/F");
  pixelTree_->Branch("eta2",pev_.eta2,"eta2[nhits2]/F");
  pixelTree_->Branch("phi2",pev_.phi2,"phi2[nhits2]/F");
  pixelTree_->Branch("r2",pev_.r2,"r2[nhits2]/F");
  pixelTree_->Branch("id2",pev_.id2,"id2[nhits2]/I");
  pixelTree_->Branch("cs2",pev_.cs2,"cs2[nhits2]/F");
  pixelTree_->Branch("cg2",pev_.cg2,"cg2[nhits2]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelHitAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelHitAnalyzer);
