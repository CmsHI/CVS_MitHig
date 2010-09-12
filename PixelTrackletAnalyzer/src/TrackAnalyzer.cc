// -*- C++ -*-
//
// Package:    TrackAnalyzer
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc MitHig/TrackAnalyzer/src/TrackAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     Prepare the Hit Tree for analysis
*/
//
// Original Author:  Yilmaz Yetkin, Yen-Jie 
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id: TrackAnalyzer.cc,v 1.25 2010/09/07 13:08:52 yjlee Exp $
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

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

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

#define MAXPARTICLES 500000
#define MAXHITS 50000
#define MAXVTX 100
#define MAXHLTBITS 100

struct PixelEvent{

   int nRun;
   int nEv;
   int nLumi;
   int nBX;
   
   int nhits1;
   int nhits2;
   int nhits3;
   int nhitsF1;
   int nhitsF2;
   int ntrks;
   int ntrksCut;
   
   int mult;
   
   // vertex
   int nv;
   float vz[MAXVTX];
   float vzMinDeltaR;

   // pixel hit
   float eta1[MAXHITS];
   float phi1[MAXHITS];
   float r1[MAXHITS];
   int id1[MAXHITS];
   float cs1[MAXHITS];
   float ch1[MAXHITS];
   int gp1[MAXHITS];
   int type1[MAXHITS];

   float eta2[MAXHITS];
   float phi2[MAXHITS];
   float r2[MAXHITS];
   int id2[MAXHITS];
   float cs2[MAXHITS];
   float ch2[MAXHITS];
   int gp2[MAXHITS];
   int type2[MAXHITS];

   float eta3[MAXHITS];
   float phi3[MAXHITS];
   float r3[MAXHITS];
   int id3[MAXHITS];
   float cs3[MAXHITS];
   float ch3[MAXHITS];
   int gp3[MAXHITS];
   int type3[MAXHITS];

   float etaF1[MAXHITS];
   float phiF1[MAXHITS];
   float rF1[MAXHITS];
   int idF1[MAXHITS];
   float csF1[MAXHITS];
   float chF1[MAXHITS];
   int gpF1[MAXHITS];
   int typeF1[MAXHITS];

   float etaF2[MAXHITS];
   float phiF2[MAXHITS];
   float rF2[MAXHITS];
   int idF2[MAXHITS];
   float csF2[MAXHITS];
   float chF2[MAXHITS];
   int gpF2[MAXHITS];
   int typeF2[MAXHITS];

   // genparticle
   int nparticle;
   float pt[MAXPARTICLES];
   float eta[MAXPARTICLES];
   float phi[MAXPARTICLES];
   int pdg[MAXPARTICLES];
   int chg[MAXPARTICLES];
   float x[MAXPARTICLES];
   float y[MAXPARTICLES];
   float z[MAXPARTICLES];
   int evtType;

   // track
   int nTrk;
   float trkEta[MAXHITS];
   float trkPhi[MAXHITS];
   float trkPt[MAXHITS];
   int trkNHit[MAXHITS];
   int trkQual[MAXHITS];
   float trkChi2[MAXHITS];
   float trkNdof[MAXHITS];
   float trkD0[MAXHITS];
   float trkDz[MAXHITS];
   float trkVx[MAXHITS];
   float trkVy[MAXHITS];
   float trkVz[MAXHITS];

   // hlt
   int nHLTBit;
   bool hltBit[MAXHLTBITS];

   // l1
   int nL1TBit;
   bool l1TBit[MAXHLTBITS];
   int nL1ABit;
   bool l1ABit[MAXHLTBITS];

   // HI
   int cBin;
   int nbins;
   int binsize;
   float hf;
   float hftp;
   float hftm;
   float eb;
   float eep;
   float eem;
   float nparti;
   float npartiSigma;
   float ncoll;
   float ncollSigma;
   float nhard;
   float nhardSigma;
   float b;
   float bSigma;
   float pixel;
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
   void fillHits(const edm::Event& iEvent);
   void fillTracks(const edm::Event& iEvent);
   void fillParticles(const edm::Event& iEvent);
   void fillPixelTracks(const edm::Event& iEvent);
   void fillL1Bits(const edm::Event& iEvent);
   void fillHLTBits(const edm::Event& iEvent);
   void fillCentrality(const edm::Event& iEvent, const edm::EventSetup& iSetup);
   
   template <typename TYPE>
   void                          getProduct(const std::string name, edm::Handle<TYPE> &prod,
                                            const edm::Event &event) const;    
   template <typename TYPE>
   bool                          getProductSafe(const std::string name, edm::Handle<TYPE> &prod,
                                                const edm::Event &event) const;

   int associateSimhitToTrackingparticle(unsigned int trid );
   bool checkprimaryparticle(const TrackingParticle* tp);

      // ----------member data ---------------------------

   bool doMC_;
   bool doCentrality_;
   bool doTrackingParticle_;
   bool doPixel_;

   vector<string> vertexSrc_;
   edm::InputTag trackSrc_;
   edm::InputTag L1gtReadout_; 
   double etaMult_;

   const TrackerGeometry* geo_;
   edm::Service<TFileService> fs;           
   edm::ESHandle < ParticleDataTable > pdt;
   edm::Handle<TrackingParticleCollection> trackingParticles;

   map<int,int> tpmap_;

   std::string                   hltResName_;         //HLT trigger results name
   std::vector<std::string>      hltProcNames_;       //HLT process name(s)
   std::vector<std::string>      hltTrgNames_;        //HLT trigger name(s)

   std::vector<int>              hltTrgBits_;         //HLT trigger bit(s)
   std::vector<bool>             hltTrgDeci_;         //HLT trigger descision(s)
   std::vector<std::string>      hltTrgUsedNames_;    //HLT used trigger name(s)
   std::string                   hltUsedResName_;     //used HLT trigger results name

   // Root object
   TTree* pixelTree_;

   PixelEvent pev_;

};

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)

{
   doMC_             = iConfig.getUntrackedParameter<bool>  ("doMC",true);
   doCentrality_             = iConfig.getUntrackedParameter<bool>  ("doCentrality",true);
   doTrackingParticle_             = iConfig.getUntrackedParameter<bool>  ("doTrackingParticle",false);
   doPixel_             = iConfig.getUntrackedParameter<bool>  ("doPixel",false);
   vertexSrc_ = iConfig.getParameter<vector<string> >("vertexSrc");
   etaMult_ = iConfig.getUntrackedParameter<double>  ("nHitsRegion",1.);
   trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
   L1gtReadout_ = iConfig.getParameter<edm::InputTag>("L1gtReadout");
   hltResName_ = iConfig.getUntrackedParameter<string>("hltTrgResults","TriggerResults");
   
   // if it's not MC, don't do TrackingParticle
   if (doMC_ == false) doTrackingParticle_ = false;
   if (iConfig.exists("hltTrgNames"))
    hltTrgNames_ = iConfig.getUntrackedParameter<vector<string> >("hltTrgNames");

   if (iConfig.exists("hltProcNames"))
      hltProcNames_ = iConfig.getUntrackedParameter<vector<string> >("hltProcNames");
   else {
      hltProcNames_.push_back("FU");
      hltProcNames_.push_back("HLT");
   }

   
}

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::~TrackAnalyzer()
{
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::ESHandle<TrackerGeometry> tGeo;
   iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
   geo_ = tGeo.product();
   iSetup.getData(pdt);

   tpmap_.clear();
   pev_.nhits1 = 0;
   pev_.nhits2 = 0;
   pev_.nhits3 = 0;
   pev_.nhitsF1 = 0;
   pev_.nhitsF2 = 0;
   pev_.ntrks = 0;
   pev_.ntrksCut = 0;
   pev_.mult = 0;
   pev_.nparticle = 0;
   pev_.nEv = (int)iEvent.id().event();
   pev_.nRun = (int)iEvent.id().run();
   pev_.nLumi = (int)iEvent.luminosityBlock();
   pev_.nBX = (int)iEvent.bunchCrossing();

   pev_.nv = 0;
   cout <<"Fill MC"<<endl;
   if (doMC_) fillParticles(iEvent);
   cout <<"Fill Vtx"<<endl;
   fillVertices(iEvent);
   cout <<"Fill Hits"<<endl;
   if (doPixel_) fillHits(iEvent);
   cout <<"Fill Tracks"<<endl;
   fillTracks(iEvent);
//   fillPixelTracks(iEvent);
   cout <<"Fill L1"<<endl;
   fillL1Bits(iEvent);
   cout <<"Fill HLT"<<endl;
   //fillHLTBits(iEvent);
   cout <<"Fill Centrality"<<endl;
   if (doCentrality_) fillCentrality(iEvent, iSetup);
   map<int,int>::iterator begin = tpmap_.begin();
   map<int,int>::iterator end = tpmap_.end();

   pixelTree_->Fill();
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillVertices(const edm::Event& iEvent){

   // Vertex 0 : pev_vz[0] MC information from TrackingVertexCollection
   // Vertex 1 - n : Reconstructed Vertex from various of algorithms
   if(doMC_){
      unsigned int daughter = 0;
      int nVertex = 0;
      int greatestvtx = 0;
      if (doTrackingParticle_) {
         Handle<TrackingVertexCollection> vertices;
         iEvent.getByLabel("mergedtruth","MergedTrackTruth", vertices);
         nVertex = vertices->size();
         for (unsigned int i = 0 ; i< vertices->size(); ++i){
   	    daughter = (*vertices)[i].nDaughterTracks();
   	    if( daughter >(*vertices)[greatestvtx].nDaughterTracks()&&fabs((*vertices)[i].position().z())<30000) greatestvtx = i;
         }
      
         if(vertices->size()>0&&fabs((*vertices)[greatestvtx].position().z())<30000){
   	    pev_.vz[pev_.nv] = (*vertices)[greatestvtx].position().z();
         }else{
	    pev_.vz[pev_.nv] =  -99; 
         }
      } else {
            pev_.vz[pev_.nv] = -99;
      }
      pev_.nv++;
   } else {
      // Fill a dummy MC information
      pev_.vz[pev_.nv] = -99;
      pev_.nv++;
   }
   
   // Fill reconstructed vertices.   
   for(unsigned int iv = 0; iv < vertexSrc_.size(); ++iv){
      const reco::VertexCollection * recoVertices;
      edm::Handle<reco::VertexCollection> vertexCollection;
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
	 pev_.vz[pev_.nv] = (*recoVertices)[greatestvtx].position().z();
      }else{
	 pev_.vz[pev_.nv] =  -99;
      }
      pev_.nv++;
   }

}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillTracks(const edm::Event& iEvent){
      Handle<vector<Track> > etracks;
      iEvent.getByLabel(trackSrc_, etracks);
      const string qualityString = "highPurity";
      pev_.nTrk=0;
      for(unsigned it=0; it<etracks->size(); ++it){
	 const reco::Track & etrk = (*etracks)[it];
         pev_.trkQual[pev_.nTrk]=0;
	 if(etrk.quality(reco::TrackBase::qualityByName(qualityString))) pev_.trkQual[pev_.nTrk]=1;
	 //if(fabs(etrk.eta())<etaCut_evtSel && etrk.pt()>ptMin_) mult++;
         pev_.trkEta[pev_.nTrk]=etrk.eta();
         pev_.trkPhi[pev_.nTrk]=etrk.phi();
         pev_.trkPt[pev_.nTrk]=etrk.pt();
         pev_.trkNHit[pev_.nTrk]=etrk.numberOfValidHits();
         pev_.trkD0[pev_.nTrk]=etrk.d0();
         pev_.trkDz[pev_.nTrk]=etrk.dz();
         pev_.trkChi2[pev_.nTrk]=etrk.chi2();
         pev_.trkNdof[pev_.nTrk]=etrk.ndof();
         pev_.trkVx[pev_.nTrk]=etrk.vx();
         pev_.trkVy[pev_.nTrk]=etrk.vy();
         pev_.trkVz[pev_.nTrk]=etrk.vz();

         //pev_.trkNhit[pev_.nTrk]=tr.numberOfValidHits();
         pev_.nTrk++;
      }
      
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillHits(const edm::Event& iEvent){

   double matchEtaMax = 0.005;
   double matchPhiMax = 0.01;

   if(doMC_&&doTrackingParticle_) iEvent.getByLabel("mergedtruth","MergedTrackTruth",trackingParticles);
   
   const SiPixelRecHitCollection* rechits;
   Handle<SiPixelRecHitCollection> rchts;
   iEvent.getByLabel("siPixelRecHits",rchts);
   rechits = rchts.product();

   for (SiPixelRecHitCollection::const_iterator it = rechits->begin(); it!=rechits->end();it++)
   {
      SiPixelRecHitCollection::DetSet hits = *it;
      DetId detId = DetId(hits.detId());
      SiPixelRecHitCollection::const_iterator recHitMatch = rechits->find(detId);
      const SiPixelRecHitCollection::DetSet recHitRange = *recHitMatch;
      unsigned int detType=detId.det();    // det type, tracker=1
      unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
      if (detType!=1) continue;

      if (subid==1) {
      PXBDetId pdetId= PXBDetId(detId);
      unsigned int layer=0;
      layer=pdetId.layer();
      for ( SiPixelRecHitCollection::DetSet::const_iterator recHitIterator = recHitRange.begin(); 
	 recHitIterator != recHitRange.end(); ++recHitIterator) {
         const SiPixelRecHit * recHit = &(*recHitIterator);

         // SIM INFO
         bool isprimary    = false;
         bool issecondary  = false;
         bool isbackground = false;
         int ptype = -99;
         int gpid = -9999;
         int trid = -9999;

         const PixelGeomDetUnit* pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit->geographicalId()));
         GlobalPoint gpos = pixelLayer->toGlobal(recHit->localPosition());
         math::XYZVector rechitPos(gpos.x(),gpos.y(),gpos.z()-pev_.vz[1]);

         // position
         double eta = rechitPos.eta();
         double phi = rechitPos.phi();
         double r   = rechitPos.rho();

         int type = -99;
	 if(isbackground) type = 0;
	 if(isprimary) type = 1;
	 if(ptype != 2) type = 2;
         if(issecondary) type = 3;

	
	 if(layer == 1){ 
	    pev_.eta1[pev_.nhits1] = eta;
	    pev_.phi1[pev_.nhits1] = phi;
	    pev_.r1[pev_.nhits1] = r;
	    pev_.id1[pev_.nhits1] = trid;
	    pev_.cs1[pev_.nhits1] = recHit->cluster()->size(); //Cluster Size
            pev_.ch1[pev_.nhits1] = recHit->cluster()->charge(); //Cluster Charge
	    pev_.gp1[pev_.nhits1] = gpid;
	    pev_.type1[pev_.nhits1] = type;
	    pev_.nhits1++;
	    if(fabs(gpos.eta()) < etaMult_ ) pev_.mult++;
	 }
	 
	 if(layer == 2){
	    pev_.eta2[pev_.nhits2] = eta;
	    pev_.phi2[pev_.nhits2] = phi;
	    pev_.r2[pev_.nhits2] = r;
            pev_.id2[pev_.nhits2] = trid;
	    pev_.cs2[pev_.nhits2] = recHit->cluster()->size(); //Cluster Size
            pev_.ch2[pev_.nhits2] = recHit->cluster()->charge(); //Cluster Charge
	    pev_.gp2[pev_.nhits2] = gpid;
	    pev_.type2[pev_.nhits2] = type;
	    pev_.nhits2++;
	 } 

	 if(layer == 3){
	    pev_.eta3[pev_.nhits3] = eta;
	    pev_.phi3[pev_.nhits3] = phi;
	    pev_.r3[pev_.nhits3] = r;
            pev_.id3[pev_.nhits3] = trid;
	    pev_.cs3[pev_.nhits3] = recHit->cluster()->size(); //Cluster Size
            pev_.ch3[pev_.nhits3] = recHit->cluster()->charge(); //Cluster Charge
	    pev_.gp3[pev_.nhits3] = gpid;
	    pev_.type3[pev_.nhits3] = type;
	    pev_.nhits3++;
	 } 
      }
      } else {
      PXFDetId pdetId= PXFDetId(detId);
      unsigned int layer=0;
      layer=pdetId.disk();
      for ( SiPixelRecHitCollection::DetSet::const_iterator recHitIterator = recHitRange.begin(); 
	 recHitIterator != recHitRange.end(); ++recHitIterator) {
         const SiPixelRecHit * recHit = &(*recHitIterator);

         // SIM INFO
         bool isprimary    = false;
         bool issecondary  = false;
         bool isbackground = false;
         int ptype = -99;
         int gpid = -9999;
         int trid = -9999;

         const PixelGeomDetUnit* pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit->geographicalId()));
         GlobalPoint gpos = pixelLayer->toGlobal(recHit->localPosition());
         math::XYZVector rechitPos(gpos.x(),gpos.y(),gpos.z()-pev_.vz[1]);

         // position
         double eta = rechitPos.eta();
         double phi = rechitPos.phi();
         double r   = rechitPos.rho();

         int type = -99;
	 if(isbackground) type = 0;
	 if(isprimary) type = 1;
	 if(ptype != 2) type = 2;
         if(issecondary) type = 3;

	
	 if(layer == 1){ 
	    pev_.etaF1[pev_.nhitsF1] = eta;
	    pev_.phiF1[pev_.nhitsF1] = phi;
	    pev_.rF1[pev_.nhitsF1] = r;
	    pev_.idF1[pev_.nhitsF1] = trid;
	    pev_.csF1[pev_.nhitsF1] = recHit->cluster()->size(); //Cluster Size
            pev_.chF1[pev_.nhitsF1] = recHit->cluster()->charge(); //Cluster Charge
	    pev_.gpF1[pev_.nhitsF1] = gpid;
	    pev_.typeF1[pev_.nhitsF1] = type;
	    pev_.nhitsF1++;
	 }
	 
	 if(layer == 2){
	    pev_.etaF2[pev_.nhitsF2] = eta;
	    pev_.phiF2[pev_.nhitsF2] = phi;
	    pev_.rF2[pev_.nhitsF2] = r;
            pev_.idF2[pev_.nhitsF2] = trid;
	    pev_.csF2[pev_.nhitsF2] = recHit->cluster()->size(); //Cluster Size
            pev_.chF2[pev_.nhitsF2] = recHit->cluster()->charge(); //Cluster Charge
	    pev_.gpF2[pev_.nhitsF2] = gpid;
	    pev_.typeF2[pev_.nhitsF2] = type;
	    pev_.nhitsF2++;
	 } 
      }
      }
   }
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillParticles(const edm::Event& iEvent)
{
   Handle<HepMCProduct> mc;
   iEvent.getByLabel("generator",mc);
   const HepMC::GenEvent* evt = mc->GetEvent();

   int evtType = evt->signal_process_id();
   
   pev_.evtType = evtType;
   
   HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
   HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
   for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
      if((*it)->status() != 1) continue;
	 tpmap_[(*it)->barcode()] = pev_.nparticle;
	 pev_.pdg[pev_.nparticle] = (*it)->pdg_id();
	 pev_.eta[pev_.nparticle] = (*it)->momentum().eta();
         pev_.phi[pev_.nparticle] = (*it)->momentum().phi();
	 pev_.pt[pev_.nparticle] = (*it)->momentum().perp();
	 const ParticleData * part = pdt->particle(pev_.pdg[pev_.nparticle]);
	 pev_.chg[pev_.nparticle] = (int)part->charge();
         pev_.x[pev_.nparticle] = (*it)->production_vertex()->position().x();
         pev_.y[pev_.nparticle] = (*it)->production_vertex()->position().y();
         pev_.z[pev_.nparticle] = (*it)->production_vertex()->position().z();
         if (pev_.chg[pev_.nparticle]==0) continue;
         if (fabs(pev_.eta[pev_.nparticle])>3) continue;
	 pev_.nparticle++;
   }
}

//--------------------------------------------------------------------------------------------------
int TrackAnalyzer::associateSimhitToTrackingparticle(unsigned int trid )
{
   int ref=-1;
   const TrackingParticleCollection* TPCProd = trackingParticles.product();
   for (TrackingParticleCollection::size_type i=0; i<TPCProd->size(); i++){
      const TrackingParticle* tp = &(*TPCProd)[i];
      vector <PSimHit> particlesimhits = tp->trackPSimHit();
      for(vector<PSimHit>::const_iterator simhit = particlesimhits.begin(); simhit != particlesimhits.end(); ++simhit)
	 {
	    //cout <<"       matching TP: "<<i<<" TPsimhitid: "<<simhit->trackId()<<" simhitId: "<<trid<<endl;
	    if(simhit->trackId()==trid)//  checkprimaryparticle(tp))
	       {
		  ref=i;
		  break;
	       }
	 }
      if (ref!=-1) break;
   }

   return ref;
}

//--------------------------------------------------------------------------------------------------   
bool TrackAnalyzer::checkprimaryparticle(const TrackingParticle* tp)
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
TrackAnalyzer::beginJob()
{

  pixelTree_ = fs->make<TTree>("PixelTree","Tree of Pixel Hits");
  pixelTree_->Branch("nEv",&pev_.nEv,"nEv/I");
  pixelTree_->Branch("nLumi",&pev_.nLumi,"nLumi/I");
  pixelTree_->Branch("nBX",&pev_.nBX,"nBX/I");
  pixelTree_->Branch("nRun",&pev_.nRun,"nRun/I");
  
  pixelTree_->Branch("nhits1",&pev_.nhits1,"nhits1/I");
  pixelTree_->Branch("nhits2",&pev_.nhits2,"nhits2/I");
  pixelTree_->Branch("nhits3",&pev_.nhits3,"nhits3/I");
  pixelTree_->Branch("nhitsF1",&pev_.nhitsF1,"nhitsF1/I");
  pixelTree_->Branch("nhitsF2",&pev_.nhitsF2,"nhitsF2/I");
  pixelTree_->Branch("ntrks",&pev_.ntrks,"ntrks/I");
  pixelTree_->Branch("ntrksCut",&pev_.ntrksCut,"ntrksCut/I");
  pixelTree_->Branch("mult",&pev_.mult,"mult/I");
  pixelTree_->Branch("nv",&pev_.nv,"nv/I");
  pixelTree_->Branch("vz",pev_.vz,"vz[nv]/F");
  pixelTree_->Branch("vzMinDeltaR",&pev_.vzMinDeltaR,"vzMinDeltaR/F");
  pixelTree_->Branch("eta1",pev_.eta1,"eta1[nhits1]/F");
  pixelTree_->Branch("phi1",pev_.phi1,"phi1[nhits1]/F");
  pixelTree_->Branch("r1",pev_.r1,"r1[nhits1]/F");
  pixelTree_->Branch("id1",pev_.id1,"id1[nhits1]/I");
  pixelTree_->Branch("cs1",pev_.cs1,"cs1[nhits1]/F");
  pixelTree_->Branch("ch1",pev_.ch1,"ch1[nhits1]/F");
  pixelTree_->Branch("gp1",pev_.gp1,"gp1[nhits1]/I");
  pixelTree_->Branch("type1",pev_.type1,"type1[nhits1]/I");

  pixelTree_->Branch("eta2",pev_.eta2,"eta2[nhits2]/F");
  pixelTree_->Branch("phi2",pev_.phi2,"phi2[nhits2]/F");
  pixelTree_->Branch("r2",pev_.r2,"r2[nhits2]/F");
  pixelTree_->Branch("id2",pev_.id2,"id2[nhits2]/I");
  pixelTree_->Branch("cs2",pev_.cs2,"cs2[nhits2]/F");
  pixelTree_->Branch("ch2",pev_.ch2,"ch2[nhits2]/F");
  pixelTree_->Branch("gp2",pev_.gp2,"gp2[nhits2]/I");
  pixelTree_->Branch("type2",pev_.type2,"type2[nhits2]/I");

  pixelTree_->Branch("eta3",pev_.eta3,"eta3[nhits3]/F");
  pixelTree_->Branch("phi3",pev_.phi3,"phi3[nhits3]/F");
  pixelTree_->Branch("r3",pev_.r3,"r3[nhits3]/F");
  pixelTree_->Branch("id3",pev_.id3,"id3[nhits3]/I");
  pixelTree_->Branch("cs3",pev_.cs3,"cs3[nhits3]/F");
  pixelTree_->Branch("ch3",pev_.ch3,"ch3[nhits3]/F");
  pixelTree_->Branch("gp3",pev_.gp3,"gp3[nhits3]/I");
  pixelTree_->Branch("type3",pev_.type3,"type3[nhits3]/I");

  pixelTree_->Branch("etaF1",pev_.etaF1,"etaF1[nhitsF1]/F");
  pixelTree_->Branch("phiF1",pev_.phiF1,"phiF1[nhitsF1]/F");
  pixelTree_->Branch("rF1",pev_.rF1,"rF1[nhitsF1]/F");
  pixelTree_->Branch("idF1",pev_.idF1,"idF1[nhitsF1]/I");
  pixelTree_->Branch("csF1",pev_.csF1,"csF1[nhitsF1]/F");
  pixelTree_->Branch("chF1",pev_.chF1,"chF1[nhitsF1]/F");
  pixelTree_->Branch("gpF1",pev_.gpF1,"gpF1[nhitsF1]/I");
  pixelTree_->Branch("typeF1",pev_.typeF1,"typeF1[nhitsF1]/I");

  pixelTree_->Branch("etaF2",pev_.etaF2,"etaF2[nhitsF2]/F");
  pixelTree_->Branch("phiF2",pev_.phiF2,"phiF2[nhitsF2]/F");
  pixelTree_->Branch("rF2",pev_.rF2,"rF2[nhitsF2]/F");
  pixelTree_->Branch("idF2",pev_.idF2,"idF2[nhitsF2]/I");
  pixelTree_->Branch("csF2",pev_.csF2,"csF2[nhitsF2]/F");
  pixelTree_->Branch("chF2",pev_.chF2,"chF2[nhitsF2]/F");
  pixelTree_->Branch("gpF2",pev_.gpF2,"gpF2[nhitsF2]/I");
  pixelTree_->Branch("typeF2",pev_.typeF2,"type2[nhits2]/I");
  pixelTree_->Branch("evtType",&pev_.evtType,"evtType/I");
  pixelTree_->Branch("npart",&pev_.nparticle,"npart/I");
  pixelTree_->Branch("pt",pev_.pt,"pt[npart]/F");
  pixelTree_->Branch("eta",pev_.eta,"eta[npart]/F");
  pixelTree_->Branch("phi",pev_.phi,"phi[npart]/F");
  pixelTree_->Branch("pdg",pev_.pdg,"pdg[npart]/I");
  pixelTree_->Branch("chg",pev_.chg,"chg[npart]/I");

  /* Not needed anymore 
  pixelTree_->Branch("x",pev_.x,"x[npart]/F");
  pixelTree_->Branch("y",pev_.y,"y[npart]/F");
  pixelTree_->Branch("z",pev_.z,"z[npart]/F");
  */

  pixelTree_->Branch("nHLTBit",&pev_.nHLTBit,"nHLTBit/I");
  pixelTree_->Branch("hltBit",pev_.hltBit,"hltBit[nHLTBit]/O");

  pixelTree_->Branch("nL1T",&pev_.nL1TBit,"nL1T/I");
  pixelTree_->Branch("L1T",pev_.l1TBit,"L1T[nL1T]/O");

  pixelTree_->Branch("nL1A",&pev_.nL1ABit,"nL1A/I");
  pixelTree_->Branch("L1A",pev_.l1ABit,"L1A[nL1A]/O");

  // Tracks
  pixelTree_->Branch("nTrk",&pev_.nTrk,"nTrk/I");
  pixelTree_->Branch("trkPt",&pev_.trkPt,"trkPt[nTrk]/F");
  pixelTree_->Branch("trkNHit",&pev_.trkNHit,"trkNHit[nTrk]/I");
  pixelTree_->Branch("trkEta",&pev_.trkEta,"trkEta[nTrk]/F");
  pixelTree_->Branch("trkPhi",&pev_.trkPhi,"trkPhi[nTrk]/F");
  pixelTree_->Branch("trkQual",&pev_.trkQual,"trkQual[nTrk]/I");
  pixelTree_->Branch("trkChi2",&pev_.trkChi2,"trkChi2[nTrk]/F");
  pixelTree_->Branch("trkNdof",&pev_.trkNdof,"trkNdof[nTrk]/F");
  pixelTree_->Branch("trkD0",&pev_.trkD0,"trkD0[nTrk]/F");
  pixelTree_->Branch("trkDz",&pev_.trkDz,"trkDz[nTrk]/F");
  pixelTree_->Branch("trkVx",&pev_.trkVx,"trkVx[nTrk]/F");
  pixelTree_->Branch("trkVy",&pev_.trkVy,"trkVy[nTrk]/F");
  pixelTree_->Branch("trkVz",&pev_.trkVz,"trkVz[nTrk]/F");

  // HI related
  pixelTree_->Branch("hf",&pev_.hf,"hf/F");
  pixelTree_->Branch("hftp",&pev_.hftp,"hftp/F");
  pixelTree_->Branch("hftm",&pev_.hftm,"hftm/F");
  pixelTree_->Branch("eb",&pev_.eb,"eb/F");
  pixelTree_->Branch("eep",&pev_.eep,"eep/F");
  pixelTree_->Branch("eem",&pev_.eem,"eem/F");
  pixelTree_->Branch("cBin",&pev_.cBin,"cBin/I");
  pixelTree_->Branch("nbins",&pev_.nbins,"nbins/I"); 
  pixelTree_->Branch("binsize",&pev_.binsize,"binsize/I");
  pixelTree_->Branch("nparti",&pev_.nparti,"nparti/F");
  pixelTree_->Branch("npartiSigma",&pev_.npartiSigma,"npartiSigma/F");
  pixelTree_->Branch("ncoll",&pev_.ncoll,"ncoll/F");
  pixelTree_->Branch("ncollSigma",&pev_.ncollSigma,"ncollSigma/F");
  pixelTree_->Branch("nhard",&pev_.nhard,"nhard/F");
  pixelTree_->Branch("nhardSigma",&pev_.nhardSigma,"nhardSigma/F");
  pixelTree_->Branch("b",&pev_.b,"b/F");
  pixelTree_->Branch("bSigma",&pev_.bSigma,"bSigma/F");
  pixelTree_->Branch("pixel",&pev_.pixel,"pixel/F");

  HLTConfigProvider hltConfig;

  cout <<"Configure hlt"<<endl;
  bool isinit = false;
  string teststr;

  // setup "Any" bit
  hltTrgBits_.clear();
  hltTrgBits_.push_back(-1);
  hltTrgDeci_.clear();
  hltTrgDeci_.push_back(true);
  hltTrgUsedNames_.clear();
  hltTrgUsedNames_.push_back("Any");

  // figure out relation of trigger name to trigger bit and store used trigger names/bits
  for(size_t i=0;i<hltTrgNames_.size();++i) {
    const string &n1(hltTrgNames_.at(i));
    bool found = 0;
    for(size_t j=0;j<hltConfig.size();++j) {
      const string &n2(hltConfig.triggerName(j));
      if (n2==n1) {
        hltTrgBits_.push_back(j);
        hltTrgUsedNames_.push_back(n1);
        hltTrgDeci_.push_back(false);
        cout <<Form("Added trigger %d with name %s for bit %d", 
                     hltTrgBits_.size()-1, n1.c_str(), j)<<endl;
        found = 1;
        break;
      }
    }      
    if (!found) {
      cout <<Form("Could not find trigger bit for %s", n1.c_str())<<endl;
    }
  }

  // ensure that trigger collections are of same size
  if (hltTrgBits_.size()!=hltTrgUsedNames_.size())
    cout <<Form("Size of trigger bits not equal used names: %d %d",
                 hltTrgBits_.size(), hltTrgUsedNames_.size())<<endl;
  if (hltTrgDeci_.size()!=hltTrgUsedNames_.size())
    cout <<Form("Size of decision bits not equal names: %d %d",
                 hltTrgDeci_.size(), hltTrgUsedNames_.size())<<endl;

  
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillPixelTracks(const edm::Event& iEvent){
  // First fish the pixel tracks out of the event
  edm::Handle<reco::TrackCollection> trackCollection;
  edm::InputTag trackCollName = trackSrc_; 
  iEvent.getByLabel(trackCollName,trackCollection);
  const reco::TrackCollection tracks = *(trackCollection.product());
  reco::TrackRefVector trks;
  for (unsigned int i=0; i<tracks.size(); i++) {
    if (tracks[i].pt() > 0.15)     
      trks.push_back( reco::TrackRef(trackCollection, i) );
  }
  pev_.ntrksCut = trks.size();  
  pev_.ntrks = tracks.size();  
}

void TrackAnalyzer::fillL1Bits(const edm::Event &iEvent)
{
  edm::Handle< L1GlobalTriggerReadoutRecord >  L1GlobalTrigger;

  iEvent.getByLabel(L1gtReadout_, L1GlobalTrigger); 
  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = L1GlobalTrigger->technicalTriggerWord();

  for (int i=0; i<64;i++)
  {
    pev_.l1TBit[i] = technicalTriggerWordBeforeMask.at(i);
  }
  pev_.nL1TBit = 64;
}

void TrackAnalyzer::fillHLTBits(const edm::Event &iEvent)
{
  // Fill HLT trigger bits.

  Handle<TriggerResults> triggerResultsHLT;
  
  getProduct(hltUsedResName_, triggerResultsHLT, iEvent);

  for(size_t i=0;i<hltTrgBits_.size();++i) {
    if (hltTrgBits_.at(i)<0) 
      continue; //ignore unknown trigger 
    size_t tbit = hltTrgBits_.at(i);
    if (tbit<triggerResultsHLT->size()) {
      hltTrgDeci_[i] = triggerResultsHLT->accept(tbit);
    } else {
      cout << Form("Problem slot %i for bit %i for %s",
                   i, tbit, triggerResultsHLT->size(), hltTrgUsedNames_.at(i).c_str())<<endl;
    }
  }
  
  pev_.nHLTBit = hltTrgBits_.size();
  
  // fill correlation histogram
  for(size_t i=0;i<hltTrgBits_.size();++i) {
    pev_.hltBit[i]=false;
    if (hltTrgDeci_.at(i)) pev_.hltBit[i]=true;
  }
}

//--------------------------------------------------------------------------------------------------
void TrackAnalyzer::fillCentrality(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const CentralityBins * cbins_;
  cbins_ = getCentralityBinsFromDB(iSetup);

  edm::Handle<reco::Centrality> cent;
  iEvent.getByLabel(edm::InputTag("hiCentrality"),cent);
 

  double hf = cent->EtHFhitSum();
  cout <<hf<<endl;
  pev_.hf = hf;
  pev_.hftp = (double)cent->EtHFtowerSumPlus();
  pev_.hftm = (double)cent->EtHFtowerSumMinus();
  pev_.eb = (double)cent->EtEBSum();
  pev_.eep = (double)cent->EtEESumPlus();
  pev_.eem = (double)cent->EtEESumMinus();
  pev_.cBin = (int)cbins_->getBin(hf);
  pev_.nbins = (int)cbins_->getNbins(); 
  pev_.binsize = (int)(100/cbins_->getNbins() );
  pev_.nparti = (double)cbins_->NpartMean(hf);
  pev_.npartiSigma = (double)cbins_->NpartSigma(hf);
  pev_.ncoll = (double)cbins_->NcollMean(hf);
  pev_.ncollSigma = (double)cbins_->NcollSigma(hf);
  pev_.nhard = (double)cbins_->NhardMean(hf);
  pev_.nhardSigma = (double)cbins_->NhardSigma(hf);
  pev_.b = (double)cbins_->bMean(hf);
  pev_.bSigma = (double)cbins_->bSigma(hf);
  pev_.pixel = (double)cent->multiplicityPixel();
  
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
