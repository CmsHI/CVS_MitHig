// -*- C++ -*-
//
// Package:    SESDoubleCandidateNtupleExporter
// Class:      SESDoubleCandidateNtupleExporter
// 
/**\class SESDoubleCandidateNtupleExporter SESDoubleCandidateNtupleExporter.cc Aha/SESDoubleCandidateNtupleExporter/src/SESDoubleCandidateNtupleExporter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yi Chen
//         Created:  Tue Nov  6 14:04:33 EST 2007
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

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "TNtuple.h"
#include "TFile.h"

//
// class decleration
//

class SESDoubleCandidateNtupleExporter : public edm::EDAnalyzer
{
   public:
      explicit SESDoubleCandidateNtupleExporter(const edm::ParameterSet&);
      ~SESDoubleCandidateNtupleExporter();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      std::string input1;
      std::string input2;
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
SESDoubleCandidateNtupleExporter::SESDoubleCandidateNtupleExporter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

   input1 = iConfig.getUntrackedParameter<std::string>("input1", "genParticleCandidates");
   input2 = iConfig.getUntrackedParameter<std::string>("input2", "superclusters");
   output = iConfig.getUntrackedParameter<std::string>("output", "ughuu.root");
}


SESDoubleCandidateNtupleExporter::~SESDoubleCandidateNtupleExporter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SESDoubleCandidateNtupleExporter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

   Handle<CandidateCollection> pInput1;
   iEvent.getByLabel(InputTag(input1), pInput1);
   const CandidateCollection *collection1 = pInput1.product();
   if(collection1 == 0)
      return;

   Handle<CandidateCollection> pInput2;
   iEvent.getByLabel(InputTag(input2), pInput2);
   const CandidateCollection *collection2 = pInput2.product();
   if(collection2 == 0)
      return;

   int maxindex = (int)collection1->size();
   if(maxindex > (int)collection2->size()) {
      std::cout <<"Different Size!!*(@&!(@&!(@!^@( Gu! "<<std::endl;
      maxindex = (int)collection2->size();
   }
   for(int i = 0; i < maxindex; i++)
   {
      const Candidate &c1 = (*collection1)[i];
      const Candidate &c2 = (*collection2)[i];

      datatemp->Fill(c1.phi(), c1.eta(), c1.p(), c1.energy(), c1.charge(),c1.status(),c1.pdgId(),
         c2.phi(), c2.eta(), c2.p(), c2.energy(), c2.charge(),c2.status(),c2.pdgId());
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
SESDoubleCandidateNtupleExporter::beginJob(const edm::EventSetup&)
{
   TFile::TContext context(0);
   f = TFile::Open(output.c_str(), "recreate");

   datatemp = new TNtuple("name", "title", "phi1:eta1:p1:E1:q1:s1:id1:phi2:eta2:p2:E2:q2:s2:id2");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SESDoubleCandidateNtupleExporter::endJob()
{
   TFile::TContext context(f);

   datatemp->Write();

   f->Close();

   // delete datatemp;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SESDoubleCandidateNtupleExporter);

