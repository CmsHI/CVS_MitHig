// -*- C++ -*-
//
// Package:    VtxAnalyzer
// Class:      VtxAnalyzer
// 
/**\class VtxAnalyzer VtxAnalyzer.cc RecoHI/VtxAnalyzer/src/VtxAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Edward Wenger
//         Created:  Fri May 22 08:11:00 EDT 2009
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/HeavyIon.h"
#include "HepMC/GenEvent.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// root include file
#include "TFile.h"  
#include "TNtuple.h"
#include "TH1F.h"


//
// class declaration
//

class VtxAnalyzer : public edm::EDAnalyzer {
   public:
      explicit VtxAnalyzer(const edm::ParameterSet&);
      ~VtxAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
	

	TNtuple *nt;
	edm::Service<TFileService> f;
	
	int evtno;
	float b;
	float vzr_peak;
	float vzr_div;
	float vz_true;
	int vs_peak; // vertex collection size
	int vs_div;
	
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
VtxAnalyzer::VtxAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


VtxAnalyzer::~VtxAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VtxAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace HepMC;
	using namespace reco;

	evtno=iEvent.id().event();
	
	// Get generated info
	edm::Handle<edm::HepMCProduct> hepEv;
	iEvent.getByLabel("source",hepEv);
	const HepMC::GenEvent * inev = hepEv->GetEvent();
	HepMC::HeavyIon* hi = inev->heavy_ion();
	b=hi->impact_parameter();
	
	// Get reconstructed vertices
	edm::Handle<reco::VertexCollection> vertexCollection;
	iEvent.getByLabel("pixelVertices",vertexCollection);
	const reco::VertexCollection * vertices = vertexCollection.product();
	vs_div=vertices->size();
	if(vs_div>0) vzr_div=vertices->begin()->position().z();
	else vzr_div=-999.9;
	
	edm::Handle<reco::VertexCollection> vertexCollection2;
	iEvent.getByLabel("pixelVertexPeak",vertexCollection2);
	const reco::VertexCollection * vertices2 = vertexCollection2.product();
	vs_peak=vertices2->size();
	if(vs_peak>0) vzr_peak=vertices2->begin()->position().z();
	else vzr_peak=-999.9;
	
	// Get signal process vertex
	HepMC::GenVertex* genvtx = inev->signal_process_vertex();
	HepMC::FourVector* vtx_;
	
	if(!genvtx){
		HepMC::GenEvent::particle_const_iterator pt=inev->particles_begin();
		HepMC::GenEvent::particle_const_iterator ptend=inev->particles_end();
		while(!genvtx || ( genvtx->particles_in_size() == 1 && pt != ptend ) ){
			++pt;
			genvtx = (*pt)->production_vertex();
		}
	}
	
	vtx_ = &(genvtx->position());
	vz_true = 0.1 * vtx_->z(); // hepMC gen vtx is in mm.  everything else is cm so we divide by 10 ;)
	
	nt->Fill(evtno,b,vzr_peak,vzr_div,vz_true,vs_peak,vs_div);

   
}


// ------------ method called once each job just before starting event loop  ------------
void 
VtxAnalyzer::beginJob(const edm::EventSetup&)
{
	
	nt = f->make<TNtuple>("nt","Vertex Testing","evtno:b:vzr_peak:vzr_div:vz_true:vs_peak:vs_div");
	
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VtxAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(VtxAnalyzer);
