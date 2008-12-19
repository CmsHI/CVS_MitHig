// -*- C++ -*-
//
// Package:    HeavyIonJetAnalyzer
// Class:      HeavyIonJetAnalyzer
// 
/**\class HeavyIonJetAnalyzer HeavyIonJetAnalyzer.cc yetkin/HeavyIonJetAnalyzer/src/HeavyIonJetAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yetkin Yilmaz
//         Created:  Tue Dec 18 09:44:41 EST 2007
// $Id: HeavyIonJetAnalyzer.cc,v 1.2 2008/08/29 17:24:35 yilmaz Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"


#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// root include file
#include "TFile.h"
#include "TNtuple.h"

using namespace std;


#define PI 3.14159265358979

#define MAXPARTICLES 5000000
#define MAXJETS 500000

#define MAXHITS 50000
#define MAXVTX 1000
#define ETABINS 3 // Fix also in branch string

//
// class decleration
//

struct HydjetEvent{

   int event;
   float b;
   float npart;
   float ncoll;
   float nhard;

   int n[ETABINS];
   float ptav[ETABINS];

   int np;
   float par_pt[MAXPARTICLES];
   float par_eta[MAXPARTICLES];
   float par_phi[MAXPARTICLES];
   int pdg[MAXPARTICLES];
   int chg[MAXPARTICLES];

   int algos;
   int njet;

   float et[MAXJETS];
   float eta[MAXJETS];
   float phi[MAXJETS];
   float area[MAXJETS];
   
   float vx;
   float vy;
   float vz;
   float vr;

};

class HeavyIonJetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HeavyIonJetAnalyzer(const edm::ParameterSet&);
      ~HeavyIonJetAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

   std::ofstream out_b;
   std::string fBFileName;

   std::ofstream out_n;
   std::string fNFileName;

   std::ofstream out_m;
   std::string fMFileName;

  
   TTree* hydjetTree_;
   HydjetEvent hev_;

   TNtuple *ntpart;

   std::string output;           // Output filename
 
   bool doAnalysis_;
   bool doParticles_;
   bool printLists_;
   bool doCF_;
   double etaMax_;
   double ptMin_;

   vector<string> jetSrc_;

   edm::ESHandle < ParticleDataTable > pdt;
   edm::Service<TFileService> f;


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
HeavyIonJetAnalyzer::HeavyIonJetAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   fBFileName = iConfig.getUntrackedParameter<std::string>("output_b", "b_values.txt");
   fNFileName = iConfig.getUntrackedParameter<std::string>("output_n", "n_values.txt");
   fMFileName = iConfig.getUntrackedParameter<std::string>("output_m", "m_values.txt");
   doAnalysis_ = iConfig.getUntrackedParameter<bool>("doAnalysis", true);
   doParticles_ = iConfig.getUntrackedParameter<bool>("doParticles", true);
   printLists_ = iConfig.getUntrackedParameter<bool>("printLists", false);
   doCF_ = iConfig.getUntrackedParameter<bool>("doMixed", false);
   jetSrc_ = iConfig.getParameter<vector<string> >("jetSrc");

   etaMax_ = iConfig.getUntrackedParameter<double>("etaMax", 2);
   ptMin_ = iConfig.getUntrackedParameter<double>("ptMin", 0);

   // Output

}


HeavyIonJetAnalyzer::~HeavyIonJetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HeavyIonJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace HepMC;
  
   hev_.event = iEvent.id().event();
   for(int ieta = 0; ieta < ETABINS; ++ieta){
      hev_.n[ieta] = 0;
      hev_.ptav[ieta] = 0;
   }

   hev_.njet = 0;
   hev_.np = 0;

   double b = -1;
   int npart = -1;
   int ncoll = -1;
   int nhard = -1;
   double vx = -99;
   double vy = -99;
   double vz = -99;
   double vr = -99;
   const GenEvent* evt;
   
   if(doCF_){

     Handle<CrossingFrame<HepMCProduct> > cf;
     iEvent.getByLabel(InputTag("mix","source"),cf);

     /*

     MixCollection<HepMCProduct> mix(cf.product());

     int mixsize = mix.size();

     cout<<"Mix Collection Size: "<<mixsize<<endl;
     evt = mix.getObject(mixsize-1).GetEvent();

     MixCollection<HepMCProduct>::iterator begin = mix.begin();
     MixCollection<HepMCProduct>::iterator end = mix.end();
     
     for(MixCollection<HepMCProduct>::iterator mixit = begin; mixit != end; ++mixit){

       const GenEvent* subevt = (*mixit).GetEvent();
       int all = subevt->particles_size();
       HepMC::GenEvent::particle_const_iterator begin = subevt->particles_begin();
       HepMC::GenEvent::particle_const_iterator end = subevt->particles_end();
       for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
	 if((*it)->status() == 1){
           float pdg_id = (*it)->pdg_id();
           float eta = (*it)->momentum().eta();
           float pt = (*it)->momentum().perp();
           const ParticleData * part = pdt->particle(pdg_id );
           float charge = part->charge();
	   }
	 }
       }
       
     }     

     */

   }else{
      
      Handle<HepMCProduct> mc;
      iEvent.getByLabel("source",mc);
      evt = mc->GetEvent();

      if(doParticles_){
      int all = evt->particles_size();
      HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
      HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
      for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
	if((*it)->status() == 1){
	   int pdg_id = (*it)->pdg_id();
	   float eta = (*it)->momentum().eta();
           float phi = (*it)->momentum().phi();
	   float pt = (*it)->momentum().perp();
	  const ParticleData * part = pdt->particle(pdg_id );
	  int charge = part->charge();

	  hev_.par_pt[hev_.np] = pt;
          hev_.par_eta[hev_.np] = eta;
          hev_.par_phi[hev_.np] = phi;
          hev_.pdg[hev_.np] = pdg_id;
          hev_.chg[hev_.np] = charge;

	  eta = fabs(eta);
	  int etabin = 0;
	  if(eta > 0.5) etabin = 1; 
	  if(eta > 1.) etabin = 2;
	  if(eta < 2.){
	     hev_.ptav[etabin] += pt;
	     ++(hev_.n[etabin]);
	  }
	  ++(hev_.np);
	}
      }
      }

   }

   const HeavyIon* hi = evt->heavy_ion();
   if(hi){
      b = hi->impact_parameter();
      npart = hi->Npart_proj()+hi->Npart_targ();

      if(printLists_){
	 out_b<<b<<endl;
	 out_n<<npart<<endl;
      }

   }
   
   //   edm::Handle<reco::GenJetCollection> genjets;
   edm::Handle<reco::JetView> genjets;
   //   edm::Handle<vector<reco::Jet> > genjets;
   iEvent.getByLabel(jetSrc_[0],genjets);
   for(int ijet = 0; ijet < genjets->size(); ++ijet){

      const reco::Jet* jet = &((*genjets)[ijet]);
      //	 const reco::GenJet* jet = dynamic_cast<const reco::GenJet*>(&((*genjets)[ijet]));

      cout<<"Jet Quantities : "<<jet->pt()<<" "<<jet->eta()<<" "<<jet->phi()<<" "<<jet->jetArea()<<endl;
      
      hev_.et[hev_.njet] = jet->pt();
      hev_.eta[hev_.njet] = jet->eta();
      hev_.phi[hev_.njet] = jet->phi();
      //      hev_.area[hev_.njet] = jet->jetArea();

      ++(hev_.njet);
      
   }
   
   /*
     edm::Handle<edm::SimVertexContainer> simVertices;
   iEvent.getByType<edm::SimVertexContainer>(simVertices);

   if (! simVertices.isValid() ) throw cms::Exception("FatalError") << "No vertices found\n";
   int inum = 0;

   edm::SimVertexContainer::const_iterator it=simVertices->begin();
   SimVertex vertex = (*it);
   cout<<" Vertex position "<< inum <<" " << vertex.position().rho()<<" "<<vertex.position().z()<<endl;
   vx = vertex.position().x();
   vy = vertex.position().y();
   vz = vertex.position().z();
   vr = vertex.position().rho();

   for(int i = 0; i<3; ++i){
      hev_.ptav[i] = hev_.ptav[i]/hev_.n[i];
   }
   */

   hev_.b = b;
   hev_.npart = npart;
   hev_.ncoll = ncoll;
   hev_.nhard = nhard;
   hev_.vx = vx;
   hev_.vy = vy;
   hev_.vz = vz;
   hev_.vr = vr;

   hydjetTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HeavyIonJetAnalyzer::beginJob(const edm::EventSetup& iSetup)
{
   iSetup.getData(pdt);

   if(printLists_){
      out_b.open(fBFileName.c_str());
      if(out_b.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fBFileName;
      out_n.open(fNFileName.c_str());
      if(out_n.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fNFileName;
      out_m.open(fMFileName.c_str());
      if(out_m.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fMFileName;
   }   
   
   if(doAnalysis_){
      hydjetTree_ = f->make<TTree>("hi","Tree of Hydjet Events");

      hydjetTree_->Branch("event",&hev_.event,"event/I");
      hydjetTree_->Branch("b",&hev_.b,"b/F");
      hydjetTree_->Branch("npart",&hev_.npart,"npart/F");
      hydjetTree_->Branch("ncoll",&hev_.ncoll,"ncoll/F");
      hydjetTree_->Branch("nhard",&hev_.nhard,"nhard/F");

      if(doParticles_){
	 hydjetTree_->Branch("n",hev_.n,"n[3]/I");
	 hydjetTree_->Branch("ptav",hev_.ptav,"ptav[3]/F");
	 hydjetTree_->Branch("np",&hev_.np,"np/I");
	 hydjetTree_->Branch("par_pt",hev_.par_pt,"par_pt[np]/F");
	 hydjetTree_->Branch("par_eta",hev_.par_eta,"par_eta[np]/F");
	 hydjetTree_->Branch("par_phi",hev_.par_phi,"par_phi[np]/F");
	 hydjetTree_->Branch("pdg",hev_.pdg,"pdg[np]/I");
	 hydjetTree_->Branch("chg",hev_.chg,"chg[np]/I");
      }

      hydjetTree_->Branch("vx",&hev_.vx,"vx/F");
      hydjetTree_->Branch("vy",&hev_.vy,"vy/F");
      hydjetTree_->Branch("vz",&hev_.vz,"vz/F");
      hydjetTree_->Branch("vr",&hev_.vr,"vr/F");

      hydjetTree_->Branch("njet",&hev_.njet,"njet/I");
      hydjetTree_->Branch("et",hev_.et,"et[njet]/F");
      hydjetTree_->Branch("eta",hev_.eta,"eta[njet]/F");
      hydjetTree_->Branch("phi",hev_.phi,"phi[njet]/F");
      hydjetTree_->Branch("area",hev_.area,"area[njet]/F");

   }
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyIonJetAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyIonJetAnalyzer);
