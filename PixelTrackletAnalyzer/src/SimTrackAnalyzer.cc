// -*- C++ -*-
//
// Package:    SimTrackAnalyzer
// Class:      SimTrackAnalyzer
// 
/**\class SimTrackAnalyzer SimTrackAnalyzer.cc SimTrackAnalyzer/SimTrackAnalyzer/src/SimTrackAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Arya Tafvizi, Yen-Jie Lee
//         Created:  Tue Jul 22 07:59:06 EDT 2008
// $Id: SimTrackAnalyzer.cc,v 1.1 2008/09/18 10:36:13 yilmaz Exp $
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
#include "TH2D.h"
#include "TFile.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

// ROOT includes
#include <Math/VectorUtil.h>

using namespace std;
using namespace reco;
using namespace edm;

//
// class decleration
//

class SimTrackAnalyzer : public edm::EDAnalyzer {
   public:
     explicit SimTrackAnalyzer(const edm::ParameterSet&);
      ~SimTrackAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
 
      bool analyzeHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, double z = 0, int layer1hits = 0);
   void analyzeTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup, double z = 0, int layer1hits = 0);
   void fillGeneratorInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      int associateSimhitToTrackingparticle(unsigned int trid );
      bool checkprimaryparticle(TrackingParticleRef tp);

   // ----------member data -------------------------------------------------------------

   double beta_;
   double alpha_;
   double etaMax_;
   double deltaCut_;

   int etaBins_;
   int eventCounter_;

   bool zoomSim_;
   bool doMC_;
   bool checkSecondLayer_;
   bool verbose_;
   bool useRecoVertex_;
  bool doOneEvent_;
  bool switch_;


   TNtuple* ntevent;
   TNtuple* ntgen;
   TNtuple* ntrechits;
   TNtuple* ntparticle;
   TNtuple* ntvertex;

   TNtuple* ntpixelsim1;
   TNtuple* ntpixelsim2;
   TNtuple* ntpixelrec;

   edm::Service<TFileService> fs;           
   const TrackerGeometry* trGeo;
   const CaloGeometry *caloGeo;
   const PixelGeomDetUnit* pixelLayer;
   edm::Handle<TrackingParticleCollection> trackingParticles ;
   
   float particles[48];
   float tracklets[48];
   float layer1Hits[48];
   float signalTracklets[48];
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
SimTrackAnalyzer::SimTrackAnalyzer(const edm::ParameterSet& iConfig)

{
   zoomSim_          = iConfig.getUntrackedParameter<bool>  ("investigateSimTracks",false);
   beta_             = iConfig.getUntrackedParameter<double>("inputBeta",1);
   doMC_             = iConfig.getUntrackedParameter<bool>  ("doMC",true);
   useRecoVertex_             = iConfig.getUntrackedParameter<bool>  ("useRecoVertex",true);
   checkSecondLayer_ = iConfig.getUntrackedParameter<bool>  ("checkSecondLayer", true);
   verbose_          = iConfig.getUntrackedParameter<bool>  ("verbose",true);
   doOneEvent_ = iConfig.getUntrackedParameter<bool>  ("doOneEvent",true);

   etaMax_ = 2.;
   etaBins_ = 8;

   deltaCut_ = iConfig.getUntrackedParameter<double>("deltaCut",0.2);
   eventCounter_ = 0;

   switch_ = true;

}


SimTrackAnalyzer::~SimTrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SimTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(switch_){

  using namespace edm;
  using namespace std;
  using namespace reco;
  
  math::XYZVector vertex(0,0,2);
  int greatestvtx = 0;

  vector<const SiPixelRecHit*> layer1;
  vector<const SiPixelRecHit*> layer2;

  eventCounter_++;

  // Get reconstructed vertices
  const reco::VertexCollection * recoVertices;
  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByLabel("pixelVertices",vertexCollection);
  recoVertices = vertexCollection.product();

  //Get MonteCarlo information
  Handle<TrackingVertexCollection> vertices;

  if (doMC_) 
  {
     iEvent.getByLabel("mergedtruth","MergedTrackTruth", vertices);
     iEvent.getByLabel("mergedtruth","MergedTrackTruth",trackingParticles);
  }

  // reset counters
  for (int i = 0; i < 48; i++) 
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
 
  if (doMC_ && !useRecoVertex_) {
     nVertex = vertices->size();

     for (unsigned int i = 0 ; i< vertices->size(); ++i)
     {
        daughter = (*vertices)[i].nDaughterTracks();
        if( daughter > (*vertices)[greatestvtx].nDaughterTracks())
        {
           greatestvtx = i;
        }
     }

     if(vertices->size()>0)
     {
        vertex = math::XYZVector((*vertices)[greatestvtx].position().x(),
                                 (*vertices)[greatestvtx].position().y(),
                                 (*vertices)[greatestvtx].position().z());
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

     if(recoVertices->size()>0)
     {
        vertex = math::XYZVector((*recoVertices)[greatestvtx].position().x(),
                                 (*recoVertices)[greatestvtx].position().y(),
                                 (*recoVertices)[greatestvtx].position().z());
     }

  }

  if (verbose_) cout <<"vertex: "<<vertex<<endl;

  ntvertex->Fill(vertex.x(),vertex.y(),vertex.z(),daughter,nVertex);
  
  // Prepare the reconstructed hits
  for(SiPixelRecHitCollection::id_iterator id = rechits->id_begin(); id!= rechits->id_end(); id++)
  {
     if((*id).subdetId() == int(PixelSubdetector::PixelBarrel))
     {
	PXBDetId pid(*id);
	SiPixelRecHitCollection::range range;
	int layer = pid.layer();
	if(layer == 1 || layer == 2)
        {
	   range = rechits->get(*id);
	   pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(*id));
	}
	
	for(SiPixelRecHitCollection::const_iterator recHit = range.first; recHit!= range.second; recHit++)
	{
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
	double etaBin = ietat * 0.5 - 2;
	if(eta1<etaBin || eta1>=etaBin+(2*etaMax_/etaBins_)) continue;
        ++layer1Hits[ietat];
     }
     if (fabs(eta1)<1) layer1HitInEta1_++;
      
  }

  Float_t tmpvar[32];
  for (int it = 0 ; it <8; it++) 
  {
     tmpvar[it] = tracklets[it];
     tmpvar[it+8] = signalTracklets[it];
     tmpvar[it+16] = layer1Hits[it];
  }
  ntevent->Fill(tmpvar);

  switch_ = (!analyzeHits(iEvent, iSetup, vertex.z(), layer1.size()));

  if(zoomSim_) analyzeTracks(iEvent, iSetup, vertex.z(), layer1.size());
  }
}

// ------------ method called once each job just before starting event loop  ------------   

void
SimTrackAnalyzer::beginJob(const edm::EventSetup& iSetup){

   edm::ESHandle<CaloGeometry> pGeo;
   iSetup.get<CaloGeometryRecord>().get(pGeo);
   caloGeo = pGeo.product();

   edm::ESHandle<TrackerGeometry> tGeo;
   iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
   trGeo = tGeo.product();

   ntevent =  fs->make<TNtuple>("ntevent","","trt1:trt2:trt3:trt4:trt5:trt6:trt7:trt8:strt1:strt2:strt3:strt4:strt5:strt6:strt7:strt8:hit1:hit2:hit3:hit4:hit5:hit6:hit7:hit8");
   ntparticle =  fs->make<TNtuple>("ntparticle","","eta:theta:phi:charge:energy:p:pt:px:py:pz:pdgid:status:useless:vrtxX:vrtxY:vrtxZ:simvrtxX:simvrtxY:simvrtxZ:eta1:eta2:eta3:phi1:phi2:phi3:trackid1:trackid2:trackid3:processtype1:processtype2:processtype3:energyloss1:energyloss2:energyloss3:particletype1:particletype2:particletype3:particlecounter:layer1hits");
   ntgen = fs->make<TNtuple>("ntgen","","had1:had2:had3:had4:had5:had6:had7:had8:lep1:lep2:lep3:lep4:lep5:lep6:lep7:lep8");
   ntvertex = fs->make<TNtuple>("ntvertex","","x:y:z:n:nvtx");

   ntpixelsim1 = fs->make<TNtuple>("ntl","","x:y:z:r:type:part:id");
   ntpixelsim2 = fs->make<TNtuple>("nth","","x:y:z:r:type:part:id");
   ntpixelrec = fs->make<TNtuple>("ntr","","x:y:z:r:layer");



}

// ------------ method called once each job just after ending the event loop  ------------  

void
SimTrackAnalyzer::endJob() {

}


bool
SimTrackAnalyzer::analyzeHits(const edm::Event& iEvent, const edm::EventSetup& iSetup, double simvrtxZ, int layer1hits){

   double simvrtxX = 0;
   double simvrtxY = 0;

   TH2D h("focus","",100,-0.5,0.5,100,4.,5.0);

   //Get reconstructed hits and geometry                                                                                                                 
   const SiPixelRecHitCollection* rechits;
   Handle<SiPixelRecHitCollection> rchts;
   iEvent.getByLabel("siPixelRecHits",rchts);
   rechits = rchts.product();

   Handle<PSimHitContainer> simhits1;
   iEvent.getByLabel(InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),simhits1);

   Handle<PSimHitContainer> simhits2;
   iEvent.getByLabel(InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),simhits2);

   for(int i1 = 0; i1 < simhits1->size(); ++i1){
      const PSimHit & hit = (*simhits1)[i1];
      int detid = hit.detUnitId();
      const PixelGeomDetUnit* pxlayer;
      PXBDetId pid(detid);
      pxlayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(pid));
      if(!pxlayer) continue;
      GlobalPoint gpos=pxlayer->toGlobal(hit.localPosition());
      double x = gpos.x();
      double y = gpos.y();
      double z = gpos.z() - simvrtxZ;
      double r = sqrt(x*x+y*y);
      int part = hit.particleType();
      int type = hit.processType();
      int id = hit.trackId();
      h.Fill(x,y);
      ntpixelsim1->Fill(x,y,z,r,type,part,id);
   }

   for(int i2 = 0; i2 < simhits2->size(); ++i2){
      const PSimHit & hit = (*simhits2)[i2];
      int detid = hit.detUnitId();
      const PixelGeomDetUnit* pxlayer;
      PXBDetId pid(detid);
      pxlayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(pid));
      if(!pxlayer) continue;
      GlobalPoint gpos=pxlayer->toGlobal(hit.localPosition());
      double x = gpos.x();
      double y = gpos.y();
      double z = gpos.z() - simvrtxZ;
      double r = sqrt(x*x+y*y);
      int part = hit.particleType();
      int type = hit.processType();
      int id = hit.trackId();
      h.Fill(x,y);
      ntpixelsim2->Fill(x,y,z,r,type,part,id);
   }
   for(SiPixelRecHitCollection::id_iterator id = rechits->id_begin(); id!= rechits->id_end(); id++)
      {
         if((*id).subdetId() == int(PixelSubdetector::PixelBarrel))
            {
               PXBDetId pid(*id);
	       SiPixelRecHitCollection::range range;
               int layer = pid.layer();
               if(layer == 1 || layer == 2)
                  {
                     range = rechits->get(*id);
                     pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(*id));

                     for(SiPixelRecHitCollection::const_iterator recHit = range.first; recHit!= range.second; recHit++)
                        {
                           const SiPixelRecHit & hit = *recHit;

                           pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (trGeo->idToDet(hit.geographicalId()));

                           GlobalPoint gpos1 = pixelLayer->toGlobal(hit.localPosition());

                           double x = gpos1.x();
                           double y = gpos1.y();
                           double z = gpos1.z() - simvrtxZ;
                           double r = sqrt(x*x+y*y);
                           ntpixelrec->Fill(x,y,z,r,layer);
			   h.Fill(x,y);
                        }
                  }
            }
      }
   
   if(h.GetMean(2) < 4){
      ntpixelsim1->Reset();
      ntpixelsim2->Reset();
      ntpixelrec->Reset();
      return false;
      
   }else return true;
   
}

void
SimTrackAnalyzer::analyzeTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup, double simvrtxZ, int layer1hits){

   /// For debugging / Getting detailed information about the propagation of tracks

   double simvrtxX = 0;
   double simvrtxY = 0;


  //gets the information from trackingparticle  
  const TrackingParticleCollection TPCProd = *(trackingParticles.product());
  for (TrackingParticleCollection::size_type i=0; i<TPCProd.size(); i++)
  {
     TrackingParticleRef tp(trackingParticles, i);
    
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
     int processtype1 = 0;
     int particletype1 = 0;
     double energyloss1 = 0;
     int trackid1 = 0;
     int processtype2 = 0;
     int particletype2 = 0;
     double energyloss2 = 0;
     int trackid2 = 0;
     int processtype3 = 0;
     int particletype3 = 0;
     double energyloss3 = 0;
     int trackid3 = 0;
     
    
     vector <PSimHit> particlesimhits = tp->trackPSimHit();
     for(vector<PSimHit>::const_iterator simhit = particlesimhits.begin(); simhit != particlesimhits.end(); ++simhit)
     {
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
        if(fabs(sqrt(x*x+y*y)-4.4)<1 && fabs(hiteta)<1)
        {
   	   layer=1;
	   if(hiteta<eta1 && hitprocesstype==2)
           {
	      eta1=hiteta;
	      phi1=hitphi;
	      processtype1=hitprocesstype;
	      particletype1=hitparticletype;
	      energyloss1=hitenergyloss;
	      trackid1=hittrackid;
	   }
        }  
        if(fabs(sqrt(x*x+y*y)-7.3)<1 && fabs(hiteta)<1)
        {  
	   layer=2;
	   if(hiteta<eta2 && hitprocesstype==2)
           {
	      eta2=hiteta;
	      phi2=hitphi;
	      processtype2=hitprocesstype;
	      particletype2=hitparticletype;
   	      energyloss2=hitenergyloss;
	      trackid2=hittrackid;
	   }
	
        }
        if(fabs(sqrt(x*x+y*y)-10.2)<1 && fabs(hiteta)<1)
        {  
	   layer=3;
	   if(hiteta<eta3 && hitprocesstype==2)
           {
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
     fillarray[8]=py;                 //
     fillarray[9]=pz;                 //
     fillarray[10]=pdgid;             //pdgid for the tracking particle
     fillarray[11]=status;            //status of the tracking particle
     fillarray[12]=-222;              //
     fillarray[13]=vrtxX;             //Vertex position of the tracking particle
     fillarray[14]=vrtxY;             //
     fillarray[15]=vrtxZ;             //
     fillarray[16]=simvrtxX;          //SimVertex position -- obtained by finding the vertex with the most number of daughter tracks
     fillarray[17]=simvrtxY;          //
     fillarray[18]=simvrtxZ;          //
     fillarray[19]=eta1;              //information of the hits on the 3 layers
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
     fillarray[37]=i;                 //Particle Counter
     fillarray[38]=layer1hits;
     ntparticle->Fill(fillarray);
     
  }
}

int SimTrackAnalyzer::associateSimhitToTrackingparticle(unsigned int trid )
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

bool SimTrackAnalyzer::checkprimaryparticle(TrackingParticleRef tp)
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
 DEFINE_FWK_MODULE(SimTrackAnalyzer);

void SimTrackAnalyzer::fillGeneratorInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
        double etaBin = ieta * 0.5 - 2;
        if(eta<etaBin || eta>=etaBin+0.5) continue;
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
