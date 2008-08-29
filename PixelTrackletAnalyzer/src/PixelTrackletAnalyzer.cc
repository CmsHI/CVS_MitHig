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
// Original Author:  Arya Tafvizi
//         Created:  Tue Jul 22 07:59:06 EDT 2008
// $Id: PixelTrackletAnalyzer.cc,v 1.3 2008/08/29 16:52:16 yilmaz Exp $
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

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
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

using namespace std;


//
// class decleration
//

class PixelTrackletAnalyzer : public edm::EDAnalyzer {
   public:
     explicit PixelTrackletAnalyzer(const edm::ParameterSet&);
      ~PixelTrackletAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
   void fillTrackDetails(const edm::Event& iEvent, const edm::EventSetup& iSetup, double z = 0, int layer1hits = 0);

      // ----------member data ---------------------------

   double beta_;
   double alpha_;
   double etaMax_;
   double deltaCut_;
   int etaBins_;

   bool zoomSim_;
   TNtuple* ntevent;
   TNtuple* ntgen;
  TNtuple * ntrechits;
  TNtuple * ntparticle;
  TNtuple * ntmatched;
  TNtuple * ntsim;
  edm::Service<TFileService> fs;
  const TrackerGeometry* trGeo;
  const CaloGeometry *caloGeo;

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
   zoomSim_ = iConfig.getUntrackedParameter<bool>("investigateSimTracks",false);
   beta_ = iConfig.getUntrackedParameter<double>("inputBeta",1);
   etaMax_ = 2.;
   etaBins_ = 8;

   deltaCut_ = iConfig.getUntrackedParameter<double>("deltaCut",0.2);

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
  int eventnumber=iEvent.id().event();
  
  double simvrtxX = 0;
  double simvrtxY = 0;
  double simvrtxZ = 0;
  int greatestvtx = 0;
  vector<const SiPixelRecHit*> layer1;
  vector<const SiPixelRecHit*> layer2;
  
  //Make sure to have bigger arrays than a few times your eta bins
  float particles[48] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float tracklets[48] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  double ptav = 0;
  
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

  Handle<HepMCProduct> mc;
  iEvent.getByLabel("source",mc);
  const HepMC::GenEvent* hepevt = mc->GetEvent();
  for( HepMC::GenEvent::particle_const_iterator pi = hepevt->particles_begin();
       pi != hepevt->particles_end(); pi++ ){
     HepMC::GenParticle* p = *pi;
     if(p->status() != 1) continue;
     double eta = p->momentum().eta();
     for(int ieta = 0 ; ieta < etaBins_; ++ieta){
	double etaBin = ((double)ieta-2*etaMax_)/(double)etaBins_/2*etaMax_;
	if(eta<etaBin || eta>=etaBin+0.5) continue;
	int partid = p->pdg_id();
	double pt = p->momentum().perp();
	if(abs( partid ) == 211 ||
	   abs( partid ) == 321 ||
	   abs( partid ) == 2212 ||
	   abs( partid ) == 3122  ){
	   ptav = ptav + pt;
	   particles[ieta]++;
	}
	if(abs(partid)==11 || abs(partid)==13){
	   particles[ieta+8]++;
	}
     }
  }

  ntgen->Fill(particles);

  //finds the primary vertex coordinates
  Handle<TrackingVertexCollection> vertices;
  iEvent.getByLabel("mergedtruth","MergedTrackTruth", vertices);
  int daughter = 0;
  for (int i = 0 ; i< vertices->size(); ++i){
    daughter = (*vertices)[i].nDaughterTracks();
    if( daughter > (*vertices)[greatestvtx].nDaughterTracks()){
      greatestvtx = i;
    }}
  if(vertices->size()>0){
    simvrtxZ = (*vertices)[greatestvtx].position().z();
    simvrtxY = (*vertices)[greatestvtx].position().y();
    simvrtxX = (*vertices)[greatestvtx].position().x();
  }
  
  TrackerHitAssociator theHitAssociator(iEvent);
  const SiPixelRecHitCollection* rechits;
  Handle<SiPixelRecHitCollection> rchts;
  iEvent.getByLabel("siPixelRecHits",rchts);
  rechits = rchts.product();
  const PixelGeomDetUnit* pixelLayer;
  for(SiPixelRecHitCollection::id_iterator id = rechits->id_begin(); id!= rechits->id_end(); id++)
    {
      if((*id).subdetId() == int(PixelSubdetector::PixelBarrel)){
	PXBDetId pid(*id);
	SiPixelRecHitCollection::range range;
	int layer = pid.layer();
	if(layer == 1 || layer == 2){
	  range = rechits->get(*id);
	  pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(*id));
	}
	
	
	for(SiPixelRecHitCollection::const_iterator recHit = range.first; recHit!= range.second; recHit++)
	  {
	    if(layer == 1)layer1.push_back(&(*recHit));
	    if(layer == 2)layer2.push_back(&(*recHit));
	    
	  }
      }
      
    }

  int layer1hits=layer1.size();    //total number of rechits on the first layer

  for(int i1 = 0; i1 < layer1.size(); ++i1){          //loops over and gets spatial information and associated simhits for each rechit

    bool signalcheck=0;
    double matchedeta=50;
    double matchedphi=50;
    double deltaeta=100;
    double matchedinveta2=50;
    double matchedinvphi2=50;
    double deltainveta=100;
    
    const SiPixelRecHit* recHit1 = layer1[i1];
    
    pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(recHit1->geographicalId()));
    
    GlobalPoint gpos1 = pixelLayer->toGlobal(recHit1->localPosition());
    
    double z1 = gpos1.z() - simvrtxZ;
    double x1 = gpos1.x();
    double y1 = gpos1.y();
    double phi1 = gpos1.phi();
    
    double rh1 = sqrt(x1*x1 + y1*y1);
    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double cos1 = z1/r1;
    double sin1 = rh1/r1;
    double eta1 = -1*log((1-cos1)/sin1);
    vector<PSimHit> simHits1 = theHitAssociator.associateHit(*recHit1);
    const PSimHit * bestSimHit1 = 0;
    
    for(vector<PSimHit>::const_iterator simHit1 = simHits1.begin(); simHit1!= simHits1.end(); simHit1++){   //gets the primary simhit and its specifications for the rechit 
      if(simHit1->processType() == 2)
	{
	  bestSimHit1 = &(*simHit1);
	  break;
	}}
    
    int trid = -9999; 
    if(bestSimHit1!=0){
      trid = bestSimHit1->trackId();  
    }
    for(int i2 = 0; i2 < layer2.size(); ++i2){
      const SiPixelRecHit* recHit2 = layer2[i2];
      pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(recHit2->geographicalId()));
      
      GlobalPoint gpos2 = pixelLayer->toGlobal(recHit2->localPosition());
      double z2 = gpos2.z() - simvrtxZ;
      double x2 = gpos2.x();
      double y2 = gpos2.y();
      double phi2 = gpos2.phi();
      double rh2 = sqrt(x2*x2 + y2*y2);
      double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
      double cos2 = z2/r2;
      double sin2 = rh2/r2;
      double eta2 = -1*log((1-cos2)/sin2);
      double inveta2 = -1*eta2;
      ntrechits->Fill(eta1,eta2,phi1,phi2);
      
      vector<PSimHit> simHits2 = theHitAssociator.associateHit(*recHit2);
      const PSimHit * bestSimHit2 = 0;
      for(vector<PSimHit>::const_iterator simHit2 = simHits2.begin(); simHit2!= simHits2.end(); simHit2++){
	if(simHit2->processType() == 2)
	  {
	    bestSimHit2 = &(*simHit2);
	    break;
	  }
      }
      if(bestSimHit2!=0){	//for each simhit on the first layer, finds the simhits on the second layer with the same trackid. 
	if(bestSimHit2->trackId() == trid){
	  double pabs = bestSimHit2->momentumAtEntry().mag();
	  int particletype = bestSimHit2->particleType();
	  double energyloss = bestSimHit2->energyLoss();
	  double pt = bestSimHit2->momentumAtEntry().perp();
	  ntsim->Fill(eta1,eta2,phi1,phi2,pabs,pt,particletype,energyloss);                         
	}}          
      
      if(fabs(eta1-eta2)<deltaeta){     //for each rechit on the first layer, finds the rechit on the second layer that has the closest eta value. if they are both from the same track, then signalcheck==1.
	signalcheck=0;
	deltaeta=fabs(eta1-eta2);
	matchedeta=eta2;
	matchedphi=phi2;
	if(trid!=-9999){
	  if(bestSimHit2){
	    if(bestSimHit2->trackId() == trid){
	      signalcheck=1;
	    }}}
      }
      if(fabs(eta1-inveta2)<deltainveta){
	deltainveta=fabs(eta1-inveta2);
	matchedinveta2=inveta2;
	matchedinvphi2=phi2;
      }
    }
    ntmatched->Fill(eta1,matchedeta,phi1,matchedphi,matchedinveta2,matchedinvphi2,signalcheck,layer1hits);
    
    if(deltaeta < deltaCut_){
       for(int ietat = 0 ; ietat < etaBins_; ++ietat){
	  double etaBin = ((double)ietat-2*etaMax_)/(double)etaBins_/2*etaMax_;
	  if(eta1<etaBin || eta1>=etaBin+(2*etaMax_/etaBins_)) continue;
	  ++tracklets[ietat];
       }
    }
  }

  ntevent->Fill(tracklets);

  if(zoomSim_) fillTrackDetails(iEvent, iSetup, simvrtxZ, layer1hits);

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

   ntevent =  fs->make<TNtuple>("ntevent","","trt1:trt2:trt3:trt4:trt5:trt6:trt7:trt8");

   ntmatched = fs->make<TNtuple>("ntmatched","","eta1:matchedeta:phi1:matchedphi:matchedinveta2:matchedinvphi2:signalcheck:layer1hits");
   ntrechits =  fs->make<TNtuple>("ntrechits","","eta1:eta2:phi1:phi2");
   ntparticle =  fs->make<TNtuple>("ntparticle","","eta:theta:phi:charge:energy:p:pt:px:py:pz:pdgid:status:useless:vrtxX:vrtxY:vrtxZ:simvrtxX:simvrtxY:simvrtxZ:eta1:eta2:eta3:phi1:phi2:phi3:trackid1:trackid2:trackid3:processtype1:processtype2:processtype3:energyloss1:energyloss2:energyloss3:particletype1:particletype2:particletype3:particlecounter:layer1hits");
   ntsim = fs->make<TNtuple>("ntsim","","eta1:eta2:phi1:phi2:pabs:pt:particletype:energyloss");
   ntgen = fs->make<TNtuple>("ntgen","","had1:had2:had3:had4:had5:had6:had7:had8:lep1:lep2:lep3:lep4:lep5:lep6:lep7:lep8");

}

// ------------ method called once each job just after ending the event loop  ------------  

void
PixelTrackletAnalyzer::endJob() {

}


void
PixelTrackletAnalyzer::fillTrackDetails(const edm::Event& iEvent, const edm::EventSetup& iSetup, double simvrtxZ, int layer1hits){


   /// For debugging / Getting detailed information about the propagation of tracks

   double simvrtxX = 0;
   double simvrtxY = 0;

  //TrackingParticles information  
   edm::Handle<TrackingParticleCollection> TPColl ;
  iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPColl);
  //gets the information from trackingparticle  
  const TrackingParticleCollection TPCProd = *(TPColl.product());
  for (TrackingParticleCollection::size_type i=0; i<TPCProd.size(); i++){
    TrackingParticleRef tp(TPColl, i);
    
    double eta = tp->eta();
    double theta = tp->theta();
    double charge = tp->charge();
    double phi = tp->phi();
    double p = tp->p();
    double pt = tp->pt();
    double px = tp->px();
    double py = tp->py();
    double pz = tp->pz();
    double energy = tp->energy();
    double vrtxX = tp->vx();
    double vrtxY = tp->vy();
    double vrtxZ = tp->vz();
    int pdgid = tp->pdgId();
    int status = tp->status();
    double eta1=20;
    double eta2=20;
    double eta3=20;
    double phi1=20;
    double phi2=20;
    double phi3=20;
    int processtype1;
    int particletype1;
    double energyloss1;
    int trackid1;
    int processtype2;
    int particletype2;
    double energyloss2;
    int trackid2;
    int processtype3;
    int particletype3;
    double energyloss3;
    int trackid3;
    
    
    vector <PSimHit> particlesimhits = tp->trackPSimHit();
    for(vector<PSimHit>::const_iterator simhit = particlesimhits.begin(); simhit != particlesimhits.end(); ++simhit){
      double hiteta=-99;
      double hitphi=-99;
      double hitphiatentry=-99;
      double hitenergyloss=-99;
      int hitparticletype=0;
      int hittrackid=-99;
      int hitprocesstype=0;
      
      int detid = simhit->detUnitId();
      const PixelGeomDetUnit* pxlayer;
      PXBDetId pid(detid);
      pxlayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(pid));
      if(!pxlayer) continue;
      GlobalPoint gpos=pxlayer->toGlobal(simhit->localPosition());
      double x = gpos.x();
      double y = gpos.y();
      double z = gpos.z() - simvrtxZ;
      double rh = sqrt(x*x+y*y);
      double r = sqrt(x*x+y*y+z*z);
      double cos = z/r;
      double sin = rh/r;
      hiteta = -1*log((1-cos)/sin);
      hitphi=gpos.phi();
      hitphiatentry=simhit->phiAtEntry();
      hitparticletype=simhit->particleType();
      hittrackid=simhit->trackId();
      hitenergyloss=simhit->energyLoss();
      hitprocesstype=simhit->processType();
      int layer = 0;
      // for each particle, finds out if it hits any of the pixel barrel layers and if so, records the information  
      if(fabs(sqrt(x*x+y*y)-4.4)<1 && fabs(hiteta)<1){
	layer=1;
	if(hiteta<eta1 && hitprocesstype==2){
	  eta1=hiteta;
	  phi1=hitphi;
	  processtype1=hitprocesstype;
	  particletype1=hitparticletype;
	  energyloss1=hitenergyloss;
	  trackid1=hittrackid;
	}
      } 
      if(fabs(sqrt(x*x+y*y)-7.3)<1 && fabs(hiteta)<1){
	layer=2;
	if(hiteta<eta2 && hitprocesstype==2){
	  eta2=hiteta;
	  phi2=hitphi;
	  processtype2=hitprocesstype;
	  particletype2=hitparticletype;
	  energyloss2=hitenergyloss;
	  trackid2=hittrackid;
	}
	
      }
      if(fabs(sqrt(x*x+y*y)-10.2)<1 && fabs(hiteta)<1){
	layer=3;
	if(hiteta<eta3 && hitprocesstype==2){
	  eta3=hiteta;
	  phi3=hitphi;
	  processtype3=hitprocesstype;
	  particletype3=hitparticletype;
	  energyloss3=hitenergyloss;
	  trackid3=hittrackid;
	}
      }
    }
    Float_t fillarray[39];
    fillarray[0]=eta;                //eta of the tracking particle
    fillarray[1]=theta;              //theta of the tracking particle  
    fillarray[2]=phi;                //phi of the tracking particle
    fillarray[3]=charge;             //charge of the tracking particle
    fillarray[4]=energy;             //energy of the tracking particle
    fillarray[5]=p;                  //absolute momentum of the tracking particle 
    fillarray[6]=pt;                 //pt of the tracking particle
    fillarray[7]=px;                 //x, y, and z components of the momentum of the tracking particle
    fillarray[8]=py;
    fillarray[9]=pz;
    fillarray[10]=pdgid;             //pdgid for the tracking particle
    fillarray[11]=status;            //status of the tracking particle
    fillarray[12]=-222;
    fillarray[13]=vrtxX;             //Vertex position of the tracking particle
    fillarray[14]=vrtxY;
    fillarray[15]=vrtxZ;
    fillarray[16]=simvrtxX;            //SimVertex position -- obtained by finding the vertex with the most number of daughter tracks
    fillarray[17]=simvrtxY;
    fillarray[18]=simvrtxZ;
    fillarray[19]=eta1;                  //information of the hits on the 3 layers
    fillarray[20]=eta2;
    fillarray[21]=eta3;
    fillarray[22]=phi1;
    fillarray[23]=phi2;
    fillarray[24]=phi3;
    fillarray[25]=trackid1;
    fillarray[26]=trackid2;
    fillarray[27]=trackid3;
    fillarray[28]=processtype1;
    fillarray[29]=processtype2;
    fillarray[30]=processtype3;
    fillarray[31]=energyloss1;
    fillarray[32]=energyloss2;	       
    fillarray[33]=energyloss3;
    fillarray[34]=particletype1;
    fillarray[35]=particletype2;
    fillarray[36]=particletype3;
    fillarray[37]=i;                   //Particle Counter
    fillarray[38]=layer1hits;
    ntparticle->Fill(fillarray);
    
  }
}

 //define this as a plug-in
 DEFINE_FWK_MODULE(PixelTrackletAnalyzer);
