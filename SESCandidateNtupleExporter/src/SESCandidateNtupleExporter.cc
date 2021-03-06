// -*- C++ -*-
//
// Package:    SESCandidateNtupleExporter
// Class:      SESCandidateNtupleExporter
// 
/**\class SESCandidateNtupleExporter SESCandidateNtupleExporter.cc Aha/SESCandidateNtupleExporter/src/SESCandidateNtupleExporter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yi Chen, Yen-Jie Lee
//         Created:  Tue Jan  8 16:02:03 EST 2008
// $Id: SESCandidateNtupleExporter.cc,v 1.1 2008/08/18 12:26:03 yjlee Exp $
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

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include "TNtuple.h"
#include "TFile.h"

//
// class decleration
//

class SESCandidateNtupleExporter : public edm::EDAnalyzer
{
   public:
      explicit SESCandidateNtupleExporter(const edm::ParameterSet&);
      ~SESCandidateNtupleExporter();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      std::string input;
      std::string output;
      TNtuple *datatemp;
      TFile *f;
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
SESCandidateNtupleExporter::SESCandidateNtupleExporter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

   input = iConfig.getUntrackedParameter<std::string>("input", "genParticleCandidates");
   output = iConfig.getUntrackedParameter<std::string>("output", "ughuu.root");
}


SESCandidateNtupleExporter::~SESCandidateNtupleExporter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SESCandidateNtupleExporter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

   Handle<GenParticleCollection> pInput1;
   iEvent.getByLabel(InputTag(input), pInput1);
   const GenParticleCollection *collection1 = pInput1.product();
   if(collection1 == 0)
      return;

   int maxindex = (int)collection1->size();
   for(int i = 0; i < maxindex; i++)
   {
      const Candidate &c1 = (*collection1)[i];

      double mid=0,gmid=0,ggmid=0;
      if ( c1.mother(0)!=0 ) {
         mid =c1.mother(0)->pdgId();
         if ( c1.mother(0)->mother(0)!=0 ){ 
               gmid =c1.mother(0)->mother(0)->pdgId();
               if ( c1.mother(0)->mother(0)->mother(0)!=0 ) ggmid =c1.mother(0)->mother(0)->mother(0)->pdgId();
         }
      }

      if (c1.status()==1) datatemp->Fill(c1.phi(), c1.eta(), c1.p(), c1.energy(), c1.charge(), c1.status(), c1.pdgId(), mid,gmid,ggmid);
      if (c1.pdgId()==111) datatemp->Fill(c1.phi(), c1.eta(), c1.p(), c1.energy(), c1.charge(), c1.status(), c1.pdgId(), mid,gmid,ggmid);
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
SESCandidateNtupleExporter::beginJob(const edm::EventSetup&)
{
   TFile::TContext context(0);
   f = TFile::Open(output.c_str(), "recreate");
   datatemp = new TNtuple("name", "title", "phi1:eta1:p1:E1:q1:s1:id1:mid:gmid:ggmid");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SESCandidateNtupleExporter::endJob()
{
   TFile::TContext context(f);

   datatemp->Write();

   f->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(SESCandidateNtupleExporter);
