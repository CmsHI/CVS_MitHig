// -*- C++ -*-
//
// Package:    TrackAnalyzer
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc MitHig/TrackAnalyzer/src/TrackAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     Prepare the Treack Tree for analysis
*/
//
// Original Author:  Yilmaz Yetkin, Yen-Jie Lee
// Updated: Frank Ma
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id: TrackAnalyzer.cc,v 1.6 2011/04/01 13:35:45 frankma Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <map>

// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"


// Root include files
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class decleration
//

#define PI 3.14159265358979

#define MAXTRACKS 50000
#define MAXVTX 100

struct TrackEvent{

   // event information
   int nRun;
   int nEv;
   int nLumi;
   int nBX;

   // Vertex information
   int nv;
   float vx[MAXVTX];
   float vy[MAXVTX];
   float vz[MAXVTX];
   float vxError[MAXVTX];
   float vyError[MAXVTX];
   float vzError[MAXVTX];


   // track
   int nTrk;
   float trkEta[MAXTRACKS];
   float trkPhi[MAXTRACKS];
   float trkPt[MAXTRACKS];
   float trkPtError[MAXTRACKS];
   int trkNHit[MAXTRACKS];
   int trkNlayer[MAXTRACKS];
   int trkNlayer3D[MAXTRACKS];
   int trkQual[MAXTRACKS];
   float trkChi2[MAXTRACKS];
   float trkChi2hit1D[MAXTRACKS];
   float trkNdof[MAXTRACKS];
   float trkDz[MAXTRACKS];
   float trkDz1[MAXTRACKS];               // dZ to the first vertex
   float trkDz2[MAXTRACKS];               // dZ to the second vertex
   float trkDzError[MAXTRACKS];
   float trkDzError1[MAXTRACKS];          
   float trkDzError2[MAXTRACKS];          
   float trkDxy[MAXTRACKS];
   float trkDxyBS[MAXTRACKS];
   float trkDxy1[MAXTRACKS];              // d0 to the first vertex
   float trkDxy2[MAXTRACKS];              // d0 to the second vertex
   float trkDxyError[MAXTRACKS];          
   float trkDxyErrorBS[MAXTRACKS];          
   float trkDxyError1[MAXTRACKS];          
   float trkDxyError2[MAXTRACKS];          
   float trkVx[MAXTRACKS];
   float trkVy[MAXTRACKS];
   float trkVz[MAXTRACKS];
   bool  trkFake[MAXTRACKS];

   float trkExpHit1Eta[MAXTRACKS];
   float trkExpHit2Eta[MAXTRACKS];
   float trkExpHit3Eta[MAXTRACKS];
};

class TrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackAnalyzer(const edm::ParameterSet&);
      ~TrackAnalyzer();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

   void fillVertices(const edm::Event& iEvent);
   void fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);
   bool hitDeadPXF(const reco::Track& tr);
   
   template <typename TYPE>
   void                          getProduct(const std::string name, edm::Handle<TYPE> &prod,
                                            const edm::Event &event) const;    
   template <typename TYPE>
   bool                          getProductSafe(const std::string name, edm::Handle<TYPE> &prod,
                                                const edm::Event &event) const;

   int associateSimhitToTrackingparticle(unsigned int trid );
   bool checkprimaryparticle(const TrackingParticle* tp);

      // ----------member data ---------------------------

   bool doTrack_;
   bool doTrackExtra_;
   bool doSimTrack_;

   double trackPtMin_;
   double genTrackPtMin_;
   bool fiducialCut_;
   edm::InputTag trackSrc_;
   edm::InputTag tpFakeSrc_;

   vector<string> vertexSrc_;

   const TrackerGeometry* geo_;
   edm::Service<TFileService> fs;           
   edm::ESHandle < ParticleDataTable > pdt;
   edm::Handle<TrackingParticleCollection> trackingParticles;
   edm::InputTag beamSpotProducer_;

   // Root object
   TTree* trackTree_;

   TrackEvent pev_;

};

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)

{
   doTrack_             = iConfig.getUntrackedParameter<bool>  ("doTrack",true);
   doTrackExtra_             = iConfig.getUntrackedParameter<bool>  ("doTrackExtra",false);
   doSimTrack_             = iConfig.getUntrackedParameter<bool>  ("doSimTrack",false);
   trackPtMin_             = iConfig.getUntrackedParameter<double>  ("trackPtMin",0.4);
   fiducialCut_ = (iConfig.getUntrackedParameter<bool>("fiducialCut",false));
   trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
   tpFakeSrc_ =  iConfig.getUntrackedParameter<edm::InputTag>("tpFakeSrc",edm::InputTag("cutsTPForFak"));
   vertexSrc_ = iConfig.getParameter<vector<string> >("vertexSrc");
   beamSpotProducer_  = iConfig.getUntrackedParameter<edm::InputTag>("towersSrc",edm::InputTag("offlineBeamSpot"));   

}

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::~TrackAnalyzer()
{
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get tracker geometry
   edm::ESHandle<TrackerGeometry> tGeo;
   iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
   geo_ = tGeo.product();
   iSetup.getData(pdt);

   pev_.nEv = (int)iEvent.id().event();
   pev_.nRun = (int)iEvent.id().run();
   pev_.nLumi = (int)iEvent.luminosityBlock();
   pev_.nBX = (int)iEvent.bunchCrossing();

   pev_.nv = 0;

   //cout <<"Fill Vtx"<<endl;
   fillVertices(iEvent);

   //cout <<"Fill Tracks"<<endl;
   if (doTrack_) fillTracks(iEvent, iSetup);
   trackTree_->Fill();
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillVertices(const edm::Event& iEvent){

   // Vertex 0 : pev_vz[0] MC information from TrackingVertexCollection
   // Vertex 1 - n : Reconstructed Vertex from various of algorithms
   pev_.vx[0]=0;
   pev_.vy[0]=0;
   pev_.vz[0]=0;
   pev_.vxError[0]=0;
   pev_.vyError[0]=0;
   pev_.vzError[0]=0;
   pev_.nv++;
   // Fill reconstructed vertices.   
   for(unsigned int iv = 0; iv < vertexSrc_.size(); ++iv){
      const reco::VertexCollection * recoVertices;
      edm::Handle<reco::VertexCollection> vertexCollection;
      cout <<vertexSrc_[iv]<<endl;
      iEvent.getByLabel(vertexSrc_[iv],vertexCollection);
      recoVertices = vertexCollection.product();
      unsigned int daughter = 0;
      int nVertex = 0;
      int greatestvtx = 0;
      
      nVertex = recoVertices->size();
      for (unsigned int i = 0 ; i< recoVertices->size(); ++i){
	 daughter = (*recoVertices)[i].tracksSize();
	 if( daughter > (*recoVertices)[greatestvtx].tracksSize()) greatestvtx = i;
	 //         cout <<"Vertex: "<< (*recoVertices)[i].position().z()<<" "<<daughter<<endl;
      }
      
      if(recoVertices->size()>0){
	 pev_.vx[pev_.nv] = (*recoVertices)[greatestvtx].position().x();
	 pev_.vy[pev_.nv] = (*recoVertices)[greatestvtx].position().y();
	 pev_.vz[pev_.nv] = (*recoVertices)[greatestvtx].position().z();
	 pev_.vxError[pev_.nv] = (*recoVertices)[greatestvtx].xError();
	 pev_.vyError[pev_.nv] = (*recoVertices)[greatestvtx].yError();
	 pev_.vzError[pev_.nv] = (*recoVertices)[greatestvtx].zError();

      }else{
	 pev_.vx[pev_.nv] =  -99;
	 pev_.vy[pev_.nv] =  -99;
	 pev_.vz[pev_.nv] =  -99;
	 pev_.vxError[pev_.nv] =  -99;
	 pev_.vyError[pev_.nv] =  -99;
	 pev_.vzError[pev_.nv] =  -99;
      }
      pev_.nv++;
      cout <<pev_.nv<<endl;
   }

}



//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup){
      Handle<vector<Track> > etracks;
      iEvent.getByLabel(trackSrc_, etracks);
      const string qualityString = "highPurity";
      reco::BeamSpot beamSpot;
      edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
      iEvent.getByLabel(beamSpotProducer_,recoBeamSpotHandle);
      beamSpot = *recoBeamSpotHandle;

      // do reco-to-sim association
      Handle<TrackingParticleCollection>  TPCollectionHfake;
      Handle<edm::View<reco::Track> >  trackCollection;
      iEvent.getByLabel(trackSrc_, trackCollection);
      ESHandle<TrackAssociatorBase> theAssociator;
      const TrackAssociatorByHits *theAssociatorByHits;
      reco::RecoToSimCollection recSimColl;
  
      if(doSimTrack_) {
        iEvent.getByLabel(tpFakeSrc_,TPCollectionHfake);
        iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
        theAssociatorByHits = (const TrackAssociatorByHits*) theAssociator.product();
      
 //       simRecColl= theAssociatorByHits->associateSimToReco(trackCollection,TPCollectionHeff,&iEvent);
        recSimColl= theAssociatorByHits->associateRecoToSim(trackCollection,TPCollectionHfake,&iEvent);
      }

      pev_.nTrk=0;
      for(unsigned it=0; it<etracks->size(); ++it){
	 const reco::Track & etrk = (*etracks)[it];
         if (etrk.pt()<trackPtMin_) continue;
	 if(fiducialCut_ && hitDeadPXF(etrk)) continue; // if track hits the dead region, igonore it;

         pev_.trkQual[pev_.nTrk]=0;
	 if(etrk.quality(reco::TrackBase::qualityByName(qualityString))) pev_.trkQual[pev_.nTrk]=1;

         trackingRecHit_iterator edh = etrk.recHitsEnd();
         int count1dhits=0;
         for (trackingRecHit_iterator ith = etrk.recHitsBegin(); ith != edh; ++ith) {
           const TrackingRecHit * hit = ith->get();
           DetId detid = hit->geographicalId();
           if (hit->isValid()) {
            if (typeid(*hit) == typeid(SiStripRecHit1D)) ++count1dhits;
           }
         }
         pev_.trkEta[pev_.nTrk]=etrk.eta();
         pev_.trkPhi[pev_.nTrk]=etrk.phi();
         pev_.trkPt[pev_.nTrk]=etrk.pt();
         pev_.trkPtError[pev_.nTrk]=etrk.ptError();
         pev_.trkNHit[pev_.nTrk]=etrk.numberOfValidHits();
         pev_.trkDxy[pev_.nTrk]=etrk.dxy();
         pev_.trkDxyError[pev_.nTrk]=etrk.dxyError();
         pev_.trkDz[pev_.nTrk]=etrk.dz();
         pev_.trkDzError[pev_.nTrk]=etrk.dzError();
         pev_.trkChi2[pev_.nTrk]=etrk.chi2();
         pev_.trkNdof[pev_.nTrk]=etrk.ndof();
         pev_.trkChi2hit1D[pev_.nTrk]=(etrk.chi2()+count1dhits)/double(etrk.ndof()+count1dhits);
         pev_.trkVx[pev_.nTrk]=etrk.vx();
         pev_.trkVy[pev_.nTrk]=etrk.vy();
         pev_.trkVz[pev_.nTrk]=etrk.vz();

         math::XYZPoint v1(pev_.vx[1],pev_.vy[1], pev_.vz[1]);
         pev_.trkDz1[pev_.nTrk]=etrk.dz(v1);
         pev_.trkDzError1[pev_.nTrk]=sqrt(etrk.dzError()*etrk.dzError()+pev_.vzError[1]*pev_.vzError[1]);
         pev_.trkDxy1[pev_.nTrk]=etrk.dxy(v1);
         pev_.trkDxyError1[pev_.nTrk]=sqrt(etrk.dxyError()*etrk.dxyError()+pev_.vxError[1]*pev_.vyError[1]);

         math::XYZPoint v2(pev_.vx[2],pev_.vy[2], pev_.vz[2]);
         pev_.trkDz2[pev_.nTrk]=etrk.dz(v2);
         pev_.trkDzError2[pev_.nTrk]=sqrt(etrk.dzError()*etrk.dzError()+pev_.vzError[2]*pev_.vzError[2]);
         pev_.trkDxy2[pev_.nTrk]=etrk.dxy(v2);
         pev_.trkDxyError2[pev_.nTrk]=sqrt(etrk.dxyError()*etrk.dxyError()+pev_.vxError[2]*pev_.vyError[2]);

         pev_.trkDxyBS[pev_.nTrk]=etrk.dxy(beamSpot.position());
         pev_.trkDxyErrorBS[pev_.nTrk]=sqrt(etrk.dxyError()*etrk.dxyError()+beamSpot.BeamWidthX()*beamSpot.BeamWidthY());
 
         pev_.trkNlayer[pev_.nTrk] = etrk.hitPattern().trackerLayersWithMeasurement();
         pev_.trkNlayer3D[pev_.nTrk] = etrk.hitPattern().pixelLayersWithMeasurement() + etrk.hitPattern().numberOfValidStripLayersWithMonoAndStereo();

         //
         if (doSimTrack_) {
             reco::TrackRef trackRef=reco::TrackRef(etracks,it);
             if(recSimColl.find(edm::RefToBase<reco::Track>(trackRef)) == recSimColl.end())
	     pev_.trkFake[pev_.nTrk]=1;
         }

         // Very rough estimation of the expected position of the Pixel Hit
         double r = 4.4; // averaged first layer rho
         double x = r*cos(etrk.phi())+etrk.vx();
         double y = r*sin(etrk.eta())+etrk.vy();
         double z = r/tan(atan(exp(-etrk.eta()))*2)+etrk.vz();
         ROOT::Math::XYZVector tmpVector(x-pev_.vx[1],y-pev_.vy[1],z-pev_.vz[1]);
         double eta1 = tmpVector.eta();
         double phi1 = etrk.phi();

         double r2 = 7.29; // averaged 2nd layer rho
         x = r2*cos(etrk.phi())+etrk.vx();
         y = r2*sin(etrk.eta())+etrk.vy();
         z = r2/tan(atan(exp(-etrk.eta()))*2)+etrk.vz();
         ROOT::Math::XYZVector tmpVector2(x-pev_.vx[1],y-pev_.vy[1],z-pev_.vz[1]);
         double eta2 = tmpVector2.eta();


         double r3 = 10.16; // averaged 3rd layer rho
         x = r3*cos(etrk.phi())+etrk.vx();
         y = r3*sin(etrk.eta())+etrk.vy();
         z = r3/tan(atan(exp(-etrk.eta()))*2)+etrk.vz();
         ROOT::Math::XYZVector tmpVector3(x-pev_.vx[1],y-pev_.vy[1],z-pev_.vz[1]);
         double eta3 = tmpVector3.eta();

         if (doTrackExtra_) {
            pev_.trkExpHit1Eta[pev_.nTrk]=eta1;
            pev_.trkExpHit2Eta[pev_.nTrk]=eta2;
            pev_.trkExpHit3Eta[pev_.nTrk]=eta3;
         }
         //pev_.trkNhit[pev_.nTrk]=tr.numberOfValidHits();
         pev_.nTrk++;
      }
      
}

// ---------------
bool
TrackAnalyzer::hitDeadPXF(const reco::Track& tr){

   //-----------------------------------------------
   // For a given track, check whether this contains 
   // hits on the dead region in the forward pixel 
   //-----------------------------------------------

   bool hitDeadRegion = false;

   for(trackingRecHit_iterator recHit = tr.recHitsBegin();recHit!= tr.recHitsEnd(); recHit++){

      if((*recHit)->isValid()){

	 DetId detId = (*recHit)->geographicalId();
	 if(!geo_->idToDet(detId)) continue;

	 Int_t diskLayerNum=0, bladeLayerNum=0, hcylLayerNum=0;
	 
	 unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());

	 if (subdetId == PixelSubdetector::PixelEndcap){
	    
	    PixelEndcapName pxfname(detId.rawId());
	    diskLayerNum = pxfname.diskName();
	    bladeLayerNum = pxfname.bladeName();
	    hcylLayerNum = pxfname.halfCylinder();
	    
	    // hard-coded now based on /UserCode/Appeltel/PixelFiducialRemover/pixelfiducialremover_cfg.py
	    if((bladeLayerNum==4 || bladeLayerNum==5 || bladeLayerNum==6) &&
	       (diskLayerNum==2) && (hcylLayerNum==4)) hitDeadRegion = true;
	 }
	 
      }// end of isValid
   }

   return hitDeadRegion;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrackAnalyzer::beginJob()
{

  trackTree_ = fs->make<TTree>("trackTree","Tree of Pixel Hits");

  // event
  trackTree_->Branch("nEv",&pev_.nEv,"nEv/I");
  trackTree_->Branch("nLumi",&pev_.nLumi,"nLumi/I");
  trackTree_->Branch("nBX",&pev_.nBX,"nBX/I");
  trackTree_->Branch("nRun",&pev_.nRun,"nRun/I");
  
  // vertex
  trackTree_->Branch("nv",&pev_.nv,"nv/I");
  trackTree_->Branch("vx",pev_.vx,"vx[nv]/F");
  trackTree_->Branch("vy",pev_.vy,"vy[nv]/F");
  trackTree_->Branch("vz",pev_.vz,"vz[nv]/F");

  // Tracks
  trackTree_->Branch("nTrk",&pev_.nTrk,"nTrk/I");
  trackTree_->Branch("trkPt",&pev_.trkPt,"trkPt[nTrk]/F");
  trackTree_->Branch("trkPtError",&pev_.trkPtError,"trkPtError[nTrk]/F");
  trackTree_->Branch("trkNHit",&pev_.trkNHit,"trkNHit[nTrk]/I");
  trackTree_->Branch("trkNlayer",&pev_.trkNlayer,"trkNlayer[nTrk]/I");
  trackTree_->Branch("trkNlayer3D",&pev_.trkNlayer3D,"trkNlayer3D[nTrk]/I");
  trackTree_->Branch("trkEta",&pev_.trkEta,"trkEta[nTrk]/F");
  trackTree_->Branch("trkPhi",&pev_.trkPhi,"trkPhi[nTrk]/F");
  trackTree_->Branch("trkQual",&pev_.trkQual,"trkQual[nTrk]/I");
  trackTree_->Branch("trkChi2",&pev_.trkChi2,"trkChi2[nTrk]/F");
  trackTree_->Branch("trkChi2hit1D",&pev_.trkChi2hit1D,"trkChi2hit1D[nTrk]/F");
  trackTree_->Branch("trkNdof",&pev_.trkNdof,"trkNdof[nTrk]/F");
  trackTree_->Branch("trkDz",&pev_.trkDz,"trkDz[nTrk]/F");
  trackTree_->Branch("trkDzError",&pev_.trkDzError,"trkDzError[nTrk]/F");
  trackTree_->Branch("trkDzError1",&pev_.trkDzError1,"trkDzError1[nTrk]/F");
  trackTree_->Branch("trkDzError2",&pev_.trkDzError2,"trkDzError2[nTrk]/F");
  trackTree_->Branch("trkDxy",&pev_.trkDxy,"trkDxy[nTrk]/F");
  trackTree_->Branch("trkDxyBS",&pev_.trkDxyBS,"trkDxyBS[nTrk]/F");
  trackTree_->Branch("trkDxy1",&pev_.trkDxy1,"trkDxy1[nTrk]/F");
  trackTree_->Branch("trkDxy2",&pev_.trkDxy2,"trkDxy2[nTrk]/F");
  trackTree_->Branch("trkDxyError",&pev_.trkDxyError,"trkDxyError[nTrk]/F");
  trackTree_->Branch("trkDxyErrorBS",&pev_.trkDxyErrorBS,"trkDxyErrorBS[nTrk]/F");
  trackTree_->Branch("trkDxyError1",&pev_.trkDxyError1,"trkDxyError1[nTrk]/F");
  trackTree_->Branch("trkDxyError2",&pev_.trkDxyError2,"trkDxyError2[nTrk]/F");
  trackTree_->Branch("trkDz1",&pev_.trkDz1,"trkDz1[nTrk]/F");
  trackTree_->Branch("trkDz2",&pev_.trkDz2,"trkDz2[nTrk]/F");
  trackTree_->Branch("trkVx",&pev_.trkVx,"trkVx[nTrk]/F");
  trackTree_->Branch("trkVy",&pev_.trkVy,"trkVy[nTrk]/F");
  trackTree_->Branch("trkVz",&pev_.trkVz,"trkVz[nTrk]/F");
  trackTree_->Branch("trkFake",&pev_.trkFake,"trkFake[nTrk]/O");

  // Track Extra
  if (doTrackExtra_) {
     trackTree_->Branch("trkExpHit1Eta",&pev_.trkExpHit1Eta,"trkExpHit1Eta[nTrk]/F");
     trackTree_->Branch("trkExpHit2Eta",&pev_.trkExpHit2Eta,"trkExpHit2Eta[nTrk]/F");
     trackTree_->Branch("trkExpHit3Eta",&pev_.trkExpHit3Eta,"trkExpHit3Eta[nTrk]/F");
  }

  
}

//--------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void TrackAnalyzer::getProduct(const std::string name, edm::Handle<TYPE> &prod,
                                    const edm::Event &event) const
{
  // Try to access data collection from EDM file. We check if we really get just one
  // product with the given name. If not we throw an exception.

  event.getByLabel(edm::InputTag(name),prod);
  if (!prod.isValid()) 
    throw edm::Exception(edm::errors::Configuration, "TrackAnalyzer::GetProduct()\n")
      << "Collection with label '" << name << "' is not valid" <<  std::endl;
}

//--------------------------------------------------------------------------------------------------
template <typename TYPE>
inline bool TrackAnalyzer::getProductSafe(const std::string name, edm::Handle<TYPE> &prod,
                                        const edm::Event &event) const
{
  // Try to safely access data collection from EDM file. We check if we really get just one
  // product with the given name. If not, we return false.

  if (name.size()==0)
    return false;

  try {
    event.getByLabel(edm::InputTag(name),prod);
    if (!prod.isValid()) 
      return false;
  } catch (...) {
    return false;
  }
  return true;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
