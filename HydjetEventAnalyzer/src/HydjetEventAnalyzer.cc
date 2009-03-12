// -*- C++ -*-
//
// Package:    HydjetEventAnalyzer
// Class:      HydjetEventAnalyzer
// 
/**\class HydjetEventAnalyzer HydjetEventAnalyzer.cc MitHig/HydjetEventAnalyzer/src/HydjetEventAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yen-Jie Lee and Yetkin Yilmaz
//         Created:  Mon Mar  9 06:19:40 EDT 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"   

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"   
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// root include file
#include "TFile.h"  
#include "TNtuple.h"
#include "TH1F.h"

using namespace std;
using namespace edm;

struct HydjetEvent{

   int event;
   int mult;
   float b;
   float npart;
   float ncoll;
   float nhard;
   float phi0;

};


//
// class decleration
//

class HydjetEventAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HydjetEventAnalyzer(const edm::ParameterSet&);
      ~HydjetEventAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

   TTree* hydjetTree_;
   HydjetEvent hev_;

   TNtuple *nt;
   TH1F *hDNDeta;

   std::string output;           // Output filename
 
//   edm::ESHandle < ParticleDataTable > pdt;
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
HydjetEventAnalyzer::HydjetEventAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


HydjetEventAnalyzer::~HydjetEventAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HydjetEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace HepMC;
   using namespace reco;

   Handle<GenParticleCollection> pInput1;
   iEvent.getByLabel(InputTag("genParticles"), pInput1);
 
   const GenParticleCollection *collection1 = pInput1.product();
   if(collection1 == 0)
      return;


   hev_.event = iEvent.id().event();
   hev_.mult = 0;
      
   double phi0 = 0;
   double b = -1;  
   int npart = -1; 
   int ncoll = -1; 
   int nhard = -1; 

   const GenEvent* evt;
  
   Handle<HepMCProduct> mc;
   iEvent.getByLabel("source",mc);
   evt = mc->GetEvent();
      

   const HeavyIon* hi = evt->heavy_ion();
   if(hi){
      b = hi->impact_parameter();
      npart = hi->Npart_proj()+hi->Npart_targ();
      ncoll = hi->Ncoll();
      nhard = hi->Ncoll_hard();
      phi0 = hi->event_plane_angle();
   }   

   hev_.b = b;
   hev_.npart = npart;
   hev_.ncoll = ncoll;

   int maxindex = (int)collection1->size();
   int multiplicity = 0;
   int chargedMultiplicity = 0;
   for(int i = 0; i < maxindex; i++)
   {
      const Candidate &c1 = (*collection1)[i];

      if (c1.status()==1 && fabs(c1.eta())<1 && c1.charge()!=0) hev_.mult++;
      if (c1.status()==1) multiplicity++;
      if (c1.status()==1 && c1.charge()!=0) {
         chargedMultiplicity++;
         if (b>16.2215) hDNDeta->Fill(c1.eta()); // 3%
      }
   }



   nt->Fill(b,npart,ncoll,hev_.mult,nhard,phi0,multiplicity,chargedMultiplicity);
   
}


// ------------ method called once each job just before starting event loop  ------------
void  HydjetEventAnalyzer::beginJob(const edm::EventSetup& iSetup)
{

   nt = f->make<TNtuple>("nt","Mixing Analysis","b:npart:ncoll:meta1:nhard:phi0:m:cm");
   hDNDeta = f->make<TH1F>("hDNDeta","dN/d#eta",200,-10,10);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HydjetEventAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HydjetEventAnalyzer);
