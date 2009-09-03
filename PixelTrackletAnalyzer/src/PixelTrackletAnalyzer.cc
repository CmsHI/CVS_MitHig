// -*- C++ -*-
//
// Package:    PixelTrackletAnalyzer
// Class:      PixelTrackletAnalyzer
// 
/**\class PixelTrackletAnalyzer PixelTrackletAnalyzer.cc PixelTrackletAnalyzer/PixelTrackletAnalyzer/src/PixelTrackletAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Arya Tafvizi, Yen-Jie Lee
//         Created:  Tue Jul 22 07:59:06 EDT 2008
// $Id: PixelTrackletAnalyzer.cc,v 1.17 2008/10/15 08:38:20 yilmaz Exp $
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
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"


#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TNtuple.h"
#include "TH1F.h"
#include "TFile.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "MitHig/PixelTracklet/interface/Tracklet.h"
#include "MitHig/PixelTracklet/interface/TrackletCorrections.h"
#include "MitHig/PixelTrackletAnalyzer/interface/TrackletFinder.h"


// ROOT includes
#include <Math/VectorUtil.h>

using namespace std;
using namespace reco;
using namespace edm;

//
// class decleration
//

namespace {
   bool compareDeltaR(Tracklet a,Tracklet b) { return fabs(a.dR2())<fabs(b.dR2());}
   bool compareDeltaEta(Tracklet a,Tracklet b) {return fabs(a.deta())<fabs(b.deta());}
}

class PixelTrackletAnalyzer : public edm::EDAnalyzer {
   public:
     explicit PixelTrackletAnalyzer(const edm::ParameterSet&);
      ~PixelTrackletAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
 
      void fillGeneratorInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      vector<Tracklet> makeTracklets(const edm::Event& iEvent, vector<const SiPixelRecHit*> layer1, vector<const SiPixelRecHit*> layer2,math::XYZVector vertex, bool invert);
      vector<Tracklet> cleanTracklets(vector<Tracklet> input,int matchNumber = 0);
      void analyzeTracklets(vector<Tracklet> input, vector<Tracklet> invertedInput, vector<Tracklet> misMatchedInput);
      int associateSimhitToTrackingparticle(unsigned int trid );
      bool checkprimaryparticle(TrackingParticleRef tp);

   // ----------member data -------------------------------------------------------------

   double beta_;
   double alpha_;
   double etaMax_;
   double deltaCut_;

   int etaBins_;
   int eventCounter_;

   int skipBest_;
   bool trySecondVtx_;
   bool doMC_;
   bool checkSecondLayer_;
   bool verbose_;
   bool useRecoVertex_;
   bool useDeltaPhi_;
   string vertexSrc_;
   string vertexSrc2_;
   edm::ParameterSet pSet_;

   TNtuple* ntevent;
   TNtuple* ntgen;
   TNtuple* ntrechits;
   TNtuple* ntmatched;
   TNtuple* ntInvMatched;
   TNtuple* ntMisMatched;
   TNtuple* ntsim;
   TNtuple* ntvertex;

   edm::Service<TFileService> fs;           
   const TrackerGeometry* trGeo;
   const CaloGeometry *caloGeo;
   const PixelGeomDetUnit* pixelLayer;
   edm::Handle<TrackingParticleCollection> trackingParticles ;
   
   float particles[72];
   float tracklets[72];
   float layer1Hits[72];
   float signalTracklets[72];
   float layer1HitInEta1_;

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
PixelTrackletAnalyzer::PixelTrackletAnalyzer(const edm::ParameterSet& iConfig)

{
   beta_             = iConfig.getUntrackedParameter<double>("inputBeta",1);
   doMC_             = iConfig.getUntrackedParameter<bool>  ("doMC",true);
   useRecoVertex_             = iConfig.getUntrackedParameter<bool>  ("useRecoVertex",true);
   checkSecondLayer_ = iConfig.getUntrackedParameter<bool>  ("checkSecondLayer", true);
   verbose_          = iConfig.getUntrackedParameter<bool>  ("verbose",true);
   vertexSrc_ = iConfig.getUntrackedParameter<string>("vertexSrc","pixelVertices");
   trySecondVtx_ = iConfig.getUntrackedParameter<bool>  ("trySecondVertex", false);
   vertexSrc2_ = iConfig.getUntrackedParameter<string>("vertexSrc2","pixelVertexFromClusters");
   pSet_ = iConfig.getParameter<edm::ParameterSet>("AssociatorParameters");
   skipBest_             = iConfig.getUntrackedParameter<int>("skipBestMatch",1);
   useDeltaPhi_ = iConfig.getUntrackedParameter<bool>  ("useDeltaPhi", true);

   etaMax_ = 3.;
   etaBins_ = 12;

   deltaCut_ = iConfig.getUntrackedParameter<double>("deltaCut",0.2);
   eventCounter_ = 0;
}


PixelTrackletAnalyzer::~PixelTrackletAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PixelTrackletAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  math::XYZVector vertex1(-99,-99,-99);
  math::XYZVector vertex2(-99,-99,-99);
  math::XYZVector vertexsim(-99,-99,-99);
  math::XYZVector vertex(0,0,0);

  int greatestvtx = 0;
  int greatestvtx2 = 0;
  int greatestvtxsim = 0;

  vector<const SiPixelRecHit*> layer1;
  vector<const SiPixelRecHit*> layer2;

  eventCounter_++;

  // Get reconstructed vertices
  const reco::VertexCollection * recoVertices;
  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByLabel(vertexSrc_,vertexCollection);
  recoVertices = vertexCollection.product();

  const reco::VertexCollection * recoVertices2;
  edm::Handle<reco::VertexCollection> vertexCollection2;

  if(trySecondVtx_){
     iEvent.getByLabel(vertexSrc2_,vertexCollection2);
     recoVertices2 = vertexCollection2.product();
  }


  //Get MonteCarlo information
  Handle<TrackingVertexCollection> vertices;

  if (doMC_) 
  {
     iEvent.getByLabel("mergedtruth","MergedTrackTruth", vertices);
     iEvent.getByLabel("mergedtruth","MergedTrackTruth",trackingParticles);
  }

  // reset counters
  for (int i = 0; i < 72; i++) 
  {
     particles[i]       = 0;
     tracklets[i]       = 0;
     signalTracklets[i] = 0;
     layer1Hits[i]      = 0;
  }

  //Get reconstructed hits and geometry  
  const SiPixelRecHitCollection* rechits;
  Handle<SiPixelRecHitCollection> rchts;
  iEvent.getByLabel("siPixelRecHits",rchts);
  rechits = rchts.product();
      
  //Fill generator information

  if (doMC_) fillGeneratorInfo(iEvent, iSetup);

  // Prepare the primary vertex coordinates

  unsigned int daughter = 0;
  unsigned int nVertex = 0;
  unsigned int daughter2 = 0;
  unsigned int nVertex2 = 0;
  unsigned int daughtersim = 0;
  unsigned int nVertexsim = 0;

  double chi1 = -99;
  double chi2 = -99;

  if (doMC_) {
     nVertexsim = vertices->size();

     for (unsigned int i = 0 ; i< vertices->size(); ++i)
     {
        daughtersim = (*vertices)[i].nDaughterTracks();
        if( daughtersim > (*vertices)[greatestvtxsim].nDaughterTracks())
        {
           greatestvtxsim = i;
        }
     }

     if(nVertexsim>0)
     {
        vertexsim = math::XYZVector((*vertices)[greatestvtxsim].position().x(),
                                 (*vertices)[greatestvtxsim].position().y(),
                                 (*vertices)[greatestvtxsim].position().z());
     }
  }

  if (useRecoVertex_) {
     nVertex = recoVertices->size();
     for (unsigned int i = 0 ; i< recoVertices->size(); ++i)
	{
        daughter = (*recoVertices)[i].tracksSize();
        if( daughter > (*recoVertices)[greatestvtx].tracksSize())
	   {
	      greatestvtx = i;
	   }
	}
     
     if(trySecondVtx_){
	nVertex2 = recoVertices2->size();
	for (unsigned int i = 0 ; i< nVertex2; ++i)
	   {
	      daughter2 = (*recoVertices2)[i].tracksSize();
	      if( daughter2 > (*recoVertices2)[greatestvtx2].tracksSize())
		 {
		    greatestvtx2 = i;
		 }
	   }
     }
     if(nVertex>0)
	{
	   vertex1 = math::XYZVector((*recoVertices)[greatestvtx].position().x(),
				     (*recoVertices)[greatestvtx].position().y(),
				     (*recoVertices)[greatestvtx].position().z());
	   chi1 = (*recoVertices)[greatestvtx].chi2();
	}
     
     if(trySecondVtx_ && nVertex2>0){
	vertex2 = math::XYZVector((*recoVertices2)[greatestvtx2].position().x(),
				  (*recoVertices2)[greatestvtx2].position().y(),
				  (*recoVertices2)[greatestvtx2].position().z());
	chi2 = (*recoVertices2)[greatestvtx2].chi2();
     }
  }
  
  ntvertex->Fill(vertex1.x(),vertex1.y(),vertex1.z(),vertex2.z(),vertexsim.z(),daughter,daughter2,daughtersim,nVertex,nVertex2,nVertexsim,chi1,chi2);
  
  if(!useRecoVertex_ && doMC_ && nVertexsim > 0) vertex = vertexsim;
  if(useRecoVertex_){
     if(trySecondVtx_ && nVertex2>0) vertex = vertex2;
     if(nVertex>0) vertex = vertex1;
  }  
  if (verbose_) cout <<"vertex: "<<vertex<<endl;
 
  // Prepare the reconstructed hits
  
  /*
  for(SiPixelRecHitCollection::const_iterator recHitIdIterator = rechits->begin(); recHitIdIterator!= rechits->end(); recHitIdIterator++)
  {
     SiPixelRecHitCollection::DetSet hits = *recHitIdIterator;
     DetId detId = DetId(hits.detId()); // Get the Detid object
     
     if((detId).subdetId() == int(PixelSubdetector::PixelBarrel))
     {
	PXBDetId pid(detId);
	SiPixelRecHitCollection::range range;
	int layer = pid.layer();
	if(layer == 1 || layer == 2)
        {
	   range = rechits->get(detId);
	   pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(detId));
	}
	
	for(SiPixelRecHitCollection::const_iterator recHit = range.first; recHit!= range.second; recHit++)
	{
	   if(layer == 1) layer1.push_back(&(*recHit));
	   if(layer == 2) layer2.push_back(&(*recHit));
	}
     }
  }
   */
   for (SiPixelRecHitCollection::const_iterator it = rechits->begin(); it!=rechits->end();it++)
   {
      SiPixelRecHitCollection::DetSet hits = *it;
      DetId detId = DetId(hits.detId());
      SiPixelRecHitCollection::const_iterator recHitMatch = rechits->find(detId);
      const SiPixelRecHitCollection::DetSet recHitRange = *recHitMatch;
      unsigned int detType=detId.det();    // det type, tracker=1
      unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
      if (detType!=1||subid!=1) continue;

      PXBDetId pdetId = PXBDetId(detId);
      unsigned int layer=0;
      layer=pdetId.layer();
      if (layer == 1 || layer == 2 ) {
         for ( SiPixelRecHitCollection::DetSet::const_iterator recHitIterator = recHitRange.begin(); 
	    recHitIterator != recHitRange.end(); ++recHitIterator) {
            const SiPixelRecHit * recHit = &(*recHitIterator);
    	    if(layer == 1) layer1.push_back(&(*recHit));
	    if(layer == 2) layer2.push_back(&(*recHit));        
         }
      }
   }

  layer1HitInEta1_ = 0;

  for(unsigned int i1 = 0; i1 < layer1.size(); ++i1)          //loops over and gets spatial information and associated simhits for each rechit
  {
     const SiPixelRecHit* recHit1 = layer1[i1];
    
     pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(recHit1->geographicalId()));
    
     GlobalPoint gpos1 = pixelLayer->toGlobal(recHit1->localPosition());
    
     // Calculate the rechit position with respect to the vertex
     math::XYZVector rechitPos(gpos1.x(),gpos1.y(),gpos1.z()-vertex.z());
     double eta1 = rechitPos.eta();
     for(int ietat = 0 ; ietat < etaBins_; ++ietat)
     {
	double etaBin = ietat * (2.*etaMax_/etaBins_) - etaMax_;
	if(eta1<etaBin || eta1>=etaBin+(2*etaMax_/etaBins_)) continue;
        ++layer1Hits[ietat];
     }
     if (fabs(eta1)<1) layer1HitInEta1_++;
      
  }

  vector<Tracklet> protoTracklets;
  vector<Tracklet> recoTracklets;

  vector<Tracklet> protoInvertedTracklets;
  vector<Tracklet> recoInvertedTracklets;

  vector<Tracklet> recoMisMatchedTracklets;

  protoTracklets = makeTracklets(iEvent,layer1,layer2,vertex,0);
  recoTracklets  = cleanTracklets(protoTracklets,0);

  protoInvertedTracklets = makeTracklets(iEvent,layer1,layer2,vertex,1);
  recoInvertedTracklets  = cleanTracklets(protoInvertedTracklets);

  recoMisMatchedTracklets  = cleanTracklets(protoTracklets,skipBest_);

  analyzeTracklets(recoTracklets,recoInvertedTracklets,recoMisMatchedTracklets);

  if (verbose_) cout <<"number of reconstructed Tracklets: "<<recoTracklets.size()<<endl;

  Float_t tmpvar[48];

  tmpvar[0] = eventCounter_;

  for (int it = 0 ; it <12; it++) 
  {
     tmpvar[it+1] = tracklets[it];
     tmpvar[it+13] = signalTracklets[it];
     tmpvar[it+25] = layer1Hits[it];
  }
  ntevent->Fill(tmpvar);

}

// ------------ method called once each job just before starting event loop  ------------                                                                                                                                       
void
PixelTrackletAnalyzer::beginJob(const edm::EventSetup& iSetup){

   edm::ESHandle<CaloGeometry> pGeo;
   iSetup.get<CaloGeometryRecord>().get(pGeo);
   caloGeo = pGeo.product();

   edm::ESHandle<TrackerGeometry> tGeo;
   iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
   trGeo = tGeo.product();

   ntevent =  fs->make<TNtuple>("ntevent","","evtid:trt1:trt2:trt3:trt4:trt5:trt6:trt7:trt8:trt9:trt10:trt11:trt12:strt1:strt2:strt3:strt4:strt5:strt6:strt7:strt8:strt9:strt10:strt11:strt12:hit1:hit2:hit3:hit4:hit5:hit6:hit7:hit8:hit9:hit10:hit11:hit12");
   ntmatched = fs->make<TNtuple>("ntmatched","","eta1:matchedeta:phi1:matchedphi:deta:dphi:signalCheck:tid:r1id:r2id:evtid:nhit1:sid:ptype");
   ntInvMatched = fs->make<TNtuple>("ntInvMatched","","eta1:matchedeta:phi1:matchedphi:deta:dphi:evtid:nhit1");
   ntMisMatched = fs->make<TNtuple>("ntMisMatched","","eta1:matchedeta:phi1:matchedphi:deta:dphi:signalCheck:tid:r1id:r2id:evtid:nhit1:sid:ptype");
   ntrechits =  fs->make<TNtuple>("ntrechits","","eta1:eta2:phi1:phi2");
   ntsim = fs->make<TNtuple>("ntsim","","eta1:eta2:phi1:phi2:pabs:pt:pid:ptype:energyloss:isprimary");
   ntgen = fs->make<TNtuple>("ntgen","","had1:had2:had3:had4:had5:had6:had7:had8:had9:had10:had11:had12:lep1:lep2:lep3:lep4:lep5:lep6:lep7:lep8:lep9:lep10:lep11:lep12");
   ntvertex = fs->make<TNtuple>("ntvertex","","x:y:z1:z2:zsim:ntrk1:ntrk2:ntrksim:nvtx1:nvtx2:nvtxsim:chi1:chi2");

}

// ------------ method called once each job just after ending the event loop  ------------  

void
PixelTrackletAnalyzer::endJob() {

}

// Make Tracklets from hits
vector<Tracklet> PixelTrackletAnalyzer::makeTracklets(const edm::Event& iEvent,vector<const SiPixelRecHit*> layer1, vector<const SiPixelRecHit*> layer2, math::XYZVector vertex, bool invert)
{
  vector<Tracklet> recoTracklets;
  TrackerHitAssociator theHitAssociator(iEvent,pSet_);
  
  for(unsigned int i1 = 0; i1 < layer1.size(); ++i1)          //loops over and gets spatial information and associated simhits for each rechit
  {
     // Ids
     int rechit1Type = 0;
     int rechit2Type = 0;
     int trackletType = -1;
     int signalExistCheck = 0;
     int signalCheck = 0;
    
     const SiPixelRecHit* recHit1 = layer1[i1];
    
     pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(recHit1->geographicalId()));
    
     GlobalPoint gpos1 = pixelLayer->toGlobal(recHit1->localPosition());
    
     // Calculate the rechit position with respect to the vertex
     math::XYZVector rechitPos(gpos1.x(),gpos1.y(),gpos1.z()-vertex.z());
     double phi1 = rechitPos.phi();
     double eta1 = rechitPos.eta();

     // Get the associated simHit
     int trid = -9999; 

     if (doMC_) {
        vector<PSimHit> simHits1 = theHitAssociator.associateHit(*recHit1);
        const PSimHit * bestSimHit1 = 0;
        if (verbose_) cout <<"Rechit "<<i1<<" "<<eta1<<" "<<phi1<<" "<<endl;
        if (verbose_) cout <<"Number of matched simHits:"<<simHits1.size()<<endl;
        int simIdx =0;
        for(vector<PSimHit>::const_iterator simHit1 = simHits1.begin(); simHit1!= simHits1.end(); simHit1++)   //gets the primary simhit and its specifications for the rechit 
        {  
           simIdx++;
           int associatedTPID = associateSimhitToTrackingparticle((&(*simHit1))->trackId());
           if (verbose_) cout <<"AssociatedPID: "<<associatedTPID<<" "<<endl;
           if (associatedTPID == -1) continue;    // doesn't match to any Trackingparticle

           TrackingParticleRef associatedTP(trackingParticles, associatedTPID);

           int detid = (&(*simHit1))->detUnitId();
           int ptype = (&(*simHit1))->processType();

           if (ptype != 2) continue;

           const PixelGeomDetUnit* pxlayer;
           PXBDetId pid(detid);
           pxlayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(pid));
           if(!pxlayer) continue;
           GlobalPoint gpos=pxlayer->toGlobal((&(*simHit1))->localPosition());
           GlobalPoint vgpos(gpos.x(),gpos.y(),gpos.z()-vertex.z());
           bool isprimary = checkprimaryparticle(associatedTP);

           if (verbose_) cout <<"Matched "<<simIdx<<" : "<<vgpos.eta()<<" "<<vgpos.phi()<<" , "<<"TP: "<<associatedTP->eta()<<" "<<associatedTP->phi()<<" Primary? "<<isprimary<<" Processtype:"<<ptype<<endl; 

           if (isprimary && bestSimHit1==0)
   	   { 
	     bestSimHit1 = &(*simHit1);
             break;
	   }           
        } 

        if(bestSimHit1!=0)
        {
           trid = bestSimHit1->trackId();  
        }
     }

     // Match with second layer reconstructed hits
     for(unsigned int i2 = 0; i2 < layer2.size(); ++i2)
     {
        const SiPixelRecHit* recHit2 = layer2[i2];
        pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(recHit2->geographicalId()));
        GlobalPoint gpos2 = pixelLayer->toGlobal(recHit2->localPosition());

        // Calculate the rechit position with respect to the vertex
        math::XYZVector rechit2Pos(gpos2.x(),gpos2.y(),gpos2.z()-vertex.z());
        double phi2 = rechit2Pos.phi();
        double eta2 = rechit2Pos.eta();
        if (invert==1) eta2*=-1;
        ntrechits->Fill(eta1,eta2,phi1,phi2);

        double trid2=-9999;
        int trackletProcessType=0;
        if(doMC_) {
 
           vector<PSimHit> simHits2 = theHitAssociator.associateHit(*recHit2);
           const PSimHit * bestSimHit2 = 0;
           for(vector<PSimHit>::const_iterator simHit2 = simHits2.begin(); simHit2!= simHits2.end(); simHit2++)
           {

             int ptype = (&(*simHit2))->processType();
             if (ptype != 2) continue;
              int associatedTPID2 = associateSimhitToTrackingparticle((&(*simHit2))->trackId());

              if (associatedTPID2 == -1) continue;    // doesn't match to any Trackingparticle

              TrackingParticleRef associatedTP2(trackingParticles, associatedTPID2);
              bool isprimary = checkprimaryparticle(associatedTP2);
              if (isprimary)
	      {
	         bestSimHit2 = &(*simHit2);
	         break;
	      }
           }
      
           if(bestSimHit2!=0)
           {  //for each simhit on the first layer, finds the simhits on the second layer with the same trackid. 
  	      if(bestSimHit2->trackId() == (unsigned int)trid && trid !=-9999)
              {
	         double pabs = bestSimHit2->momentumAtEntry().mag();
	         int particletype = bestSimHit2->particleType();
	         trackletProcessType = bestSimHit2->processType();
	         double energyloss = bestSimHit2->energyLoss();
	         double pt = bestSimHit2->momentumAtEntry().perp();
                 int associatedTPID = associateSimhitToTrackingparticle(trid);

                 if (associatedTPID == -1) continue;    // doesn't match to any Trackingparticle

                 TrackingParticleRef associatedTP2(trackingParticles, associatedTPID);

                 bool isprimary = checkprimaryparticle(associatedTP2);
                 float var[20];
                 var[0]=eta1;
                 var[1]=eta2;
                 var[2]=phi1;
                 var[3]=phi2;
                 var[4]=pabs;
                 var[5]=pt;
                 var[6]=particletype;
                 var[7]=trackletProcessType;
                 var[8]=energyloss;
      	         var[9]=isprimary;
   	         if (!invert) ntsim->Fill(var);                         
	      }
              trid2 = bestSimHit2->trackId();
           }          
        }
        
        // Tracklet id
        if(trid!=-9999)
        {
	   if(trid2 != -9999)
           {
	      if(trid2 == trid)
              {
	         signalCheck = 1;
	      }
           }
        }

        if(trid!=-9999)
        {
           rechit1Type = 1;
        }
           
        if(trid2!=-9999)
        {
           rechit2Type = 1;
        }

        if (rechit1Type == 1) 
        {
           if (rechit2Type == 1) 
           {
              if (signalCheck == 1 ) trackletType = 1; 
              else                   trackletType = 2;
           } else {
              trackletType = 3;
           }
         } else {
           if (rechit2Type == 1)
           {
              trackletType = 4;
           } else {
              trackletType = 5;
           }
        }

        Tracklet mytracklet(eta1,eta2,phi1,phi2);

        mytracklet.setIt1(i1);
        mytracklet.setIt2(i2);
        mytracklet.setId1(rechit1Type);
        mytracklet.setId2(rechit2Type);
        mytracklet.setId(trackletType);
        mytracklet.setSId(signalExistCheck);
        mytracklet.setType(trackletProcessType);
        recoTracklets.push_back(mytracklet);
     }
  }


  return recoTracklets;
}


// Clean the Tracklet with multiple use
vector<Tracklet> PixelTrackletAnalyzer::cleanTracklets(vector<Tracklet> input, int matchNumber)
{
   vector<Tracklet> output;

   if(useDeltaPhi_)
      sort( input.begin() , input.end() , compareDeltaR);
   else
      sort( input.begin() , input.end() , compareDeltaEta);

   if (verbose_) {
      for (unsigned int i = 0; i < input.size(); i++)
      {
         cout <<input[i].deta()<<" "<<input[i].getIt1()<<" "<<input[i].getIt2()<<endl;
      }
   }

   int used1[1000];
   int used2[1000];

   for (int i=0;i<1000;i++) { 
      used1[i]=0;
      used2[i]=0;
   } 

   cout<<"Printing Hits"<<endl;
   
   for (unsigned int i = 0; i < input.size(); i++){
      

      if(useDeltaPhi_)
	 cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta R : "<<input[i].dR()<<endl;
      else
	 cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta Eta : "<<input[i].deta()<<endl; 
      
      int i1=input[i].getIt1();
      int i2=input[i].getIt2();
      
      if (used1[i1]==0&&used2[i2]==matchNumber) {
	 Tracklet tmp = input[i];
	 output.push_back(tmp);
	 used1[i1]++;
	 if (checkSecondLayer_) used2[i2]++;
      }
      if (used1[i1]==0&&used2[i2]<matchNumber) {
	 if (checkSecondLayer_) used2[i2]++;
      }
   }
   if (verbose_) {
      cout <<"Output:"<<endl;
      for (unsigned int i = 0; i < output.size(); i++)
      {
         cout <<output[i].deta()<<" "<<output[i].getIt1()<<" "<<output[i].getIt2()<<endl;
      }
   }
   
   return output;
}

// Make Ntuple
void PixelTrackletAnalyzer::analyzeTracklets(vector<Tracklet> input, vector<Tracklet> invertedInput, vector<Tracklet> mismatchedInput)
{

  for (unsigned int i = 0; i < input.size(); i++) 
  {
     float var[100];
     int signalCheck=0;
     if (input[i].getId()==1) signalCheck=1; 

     var[0]=input[i].eta1();
     var[1]=input[i].eta2();
     var[2]=input[i].phi1();
     var[3]=input[i].phi2();
     var[4]=input[i].deta();
     var[5]=input[i].dphi();
     var[6]=signalCheck; 
     var[7]=input[i].getId();
     var[8]=input[i].getId1();
     var[9]=input[i].getId2();
     var[10]=eventCounter_;
     var[11]=layer1HitInEta1_;
     var[12]=input[i].getSId();
     var[13]=input[i].getType();
     ntmatched->Fill(var);

	cout<<"Correctly matched Entries :: "<<"Eta 1 : "<<var[0]<<" ; Eta 2 :  "<<var[1]<<" DeltaR : "<<sqrt(var[4]*var[4]+var[5]*var[5])<<endl;
    
     if(input[i].deta() < deltaCut_)
     {
        for(int ietat = 0 ; ietat < etaBins_; ++ietat)
        {
	   double etaBin = ietat * (2.*etaMax_/etaBins_) - etaMax_;
	   if(input[i].eta1()<etaBin || input[i].eta1()>=etaBin+(2*etaMax_/etaBins_)) continue;

	   ++tracklets[ietat];
	   if (signalCheck) ++signalTracklets[ietat];
        }
     }
  }

  for (unsigned int i = 0; i < invertedInput.size(); i++) 
  {
     float var[100];

     var[0]=invertedInput[i].eta1();
     var[1]=invertedInput[i].eta2();
     var[2]=invertedInput[i].phi1();
     var[3]=invertedInput[i].phi2();
     var[4]=invertedInput[i].deta();
     var[5]=invertedInput[i].dphi();
     var[6]=eventCounter_;
     var[7]=layer1HitInEta1_;

     ntInvMatched->Fill(var);
  }

  for (unsigned int i = 0; i < mismatchedInput.size(); i++)
     {
	float var[100];
	int signalCheck=0;
	if (input[i].getId()==1) signalCheck=1;
	var[0]=mismatchedInput[i].eta1();
	var[1]=mismatchedInput[i].eta2();
	var[2]=mismatchedInput[i].phi1();
	var[3]=mismatchedInput[i].phi2();
	var[4]=mismatchedInput[i].deta();
	var[5]=mismatchedInput[i].dphi();
	var[6]=signalCheck;
	var[7]=mismatchedInput[i].getId();
	var[8]=mismatchedInput[i].getId1();
	var[9]=mismatchedInput[i].getId2();
	var[10]=eventCounter_;
	var[11]=layer1HitInEta1_;
	var[12]=mismatchedInput[i].getSId();
	var[13]=mismatchedInput[i].getType();
	ntMisMatched->Fill(var);

	if(useDeltaPhi_)
	   cout<<"MisMatched Entries :: "<<"Eta 1 : "<<var[0]<<" ; Eta 2 :  "<<var[1]<<" DeltaR : "<<sqrt(var[4]*var[4]+var[5]*var[5])<<endl;
	else
	   cout<<"MisMatched Entries :: "<<"Eta 1 : "<<var[0]<<" ; Eta 2 :  "<<var[1]<<" DeltaEta : "<<fabs(var[4])<<endl;

     }

     cout<<"MisMatched Entries Size :: "<<mismatchedInput.size()<<"Coorect Entries ::  "<<input.size()<<endl;

}

int PixelTrackletAnalyzer::associateSimhitToTrackingparticle(unsigned int trid )
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

bool PixelTrackletAnalyzer::checkprimaryparticle(TrackingParticleRef tp)
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

 //define this as a plug-in
 DEFINE_FWK_MODULE(PixelTrackletAnalyzer);

void PixelTrackletAnalyzer::fillGeneratorInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Generator Level Information
  //Important information:
  //
  //The ntuple is filled in event-by-event basis
  //therefore the total number of particles in a certain
  //eta bin is recorded by the hardcoded variable names
  //had1,had2,had3 etc...
  //The way to interpret the etabins is this:
  //had1 : -2, -1.5
  //had2 : -1.5, 1
  //...
  //had7 : 1, 1.5
  //had8 : 1.5, 2 
  // had: number of hadrons
  // lep: number of leptons
  // trt: number of tracklets
  // signalTrt: number of signal tracklets

  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<HepMCProduct> mc;
  iEvent.getByLabel("source",mc);


  const HepMC::GenEvent* hepevt = mc->GetEvent();
  for( HepMC::GenEvent::particle_const_iterator pi = hepevt->particles_begin();
       pi != hepevt->particles_end(); pi++ )
  {
     HepMC::GenParticle* p = *pi;
     int partid = p->pdg_id();  
     double pt = p->momentum().perp();
     double eta = p->momentum().eta();
     double phi = p->momentum().phi();
     if (verbose_) cout <<"Event "<<eventCounter_<<" : "<<partid<<" pt: "<<pt<<" eta: "<<eta<<" phi: "<<phi<<" stat:"<<p->status()<<endl;
     if (p->status() != 1) continue;
  
     for(int ieta = 0 ; ieta < etaBins_; ++ieta)
     {
	double etaBin = ieta * (2.*etaMax_/etaBins_) - etaMax_;
	if(eta<etaBin || eta>=etaBin+(2*etaMax_/etaBins_)) continue;

        if(abs( partid ) == 211 ||
	   abs( partid ) == 321 ||
	   abs( partid ) == 2212 ||
	   abs( partid ) == 3122  )
        {
	   particles[ieta]++;
	}
	if(abs(partid)==11 || abs(partid)==13)
        {
	   particles[ieta+8]++;
	}
     }
  }
  cout <<"Hadron/Lepton: "<<particles<<endl;
  ntgen->Fill(particles);
}
