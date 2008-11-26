#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "MitHig/MatchedTrack/interface/MatchedTrack.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

/*
 class TransientTrackFromFTSFactory {
 public:
 
 reco::TransientTrack build (const FreeTrajectoryState & fts) const;
 reco::TransientTrack build (const FreeTrajectoryState & fts,
 const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
 };
 */

#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"

#include <fstream>
using namespace std;
using namespace reco;

/*****************************************************************************/
class HighPtTrackAnalyzer : public edm::EDAnalyzer
{
public:
	explicit HighPtTrackAnalyzer(const edm::ParameterSet& pset);
	~HighPtTrackAnalyzer();
	virtual void beginJob(const edm::EventSetup& es);
	virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
	virtual void endJob();
	
private:
	int layerFromDetid(const DetId& detId);
	int getNumberOfSimHits(const TrackingParticle& simTrack);
	int getDetLayerId(const PSimHit& simHit);
	int getNumberOfPixelHits(const TrackingParticle& simTrack,float *);
	int getNumberOfRecPixelHits(const reco::Track & recTrack, float *);
	void checkSimTracks (edm::Handle<TrackingParticleCollection>& simCollection,reco::SimToRecoCollection& q);
	
	pair<float,float> refitWithVertex(const reco::Track & recTrack,const reco::VertexCollection* vertices);
	int getParticleId(edm::RefToBase<reco::Track>& recTrack, int & ptype);
	void checkRecTracks(edm::Handle<edm::View<reco::Track> >& recCollection,const reco::VertexCollection* vertices,reco::RecoToSimCollection& p);
	
	const TrackerGeometry * theTracker;
	const TrackAssociatorByHits * theAssociatorByHits;
	const TransientTrackBuilder * theTTBuilder;
	TrackerHitAssociator * theHitAssociator;
	
	vector<string> trackCollectionLabels;
	string resultFileLabel;
	bool plotEvent, zipFiles;
	int proc;
	Int_t iTrkSim,iTrkReco,iVtx;
	Int_t iEvent,iRun;
	Float_t fSimPxlLayerHit,dMinSimPt,fRecPxlLayerHit;
	
	TFile * resultFile; 
	TTree *recInfoTree;
	TClonesArray *CAReco,*CASim;
	
	
};

/*****************************************************************************/
HighPtTrackAnalyzer::HighPtTrackAnalyzer(const edm::ParameterSet& pset)
{
	trackCollectionLabels = pset.getParameter<vector<string> >("trackCollection");
	resultFileLabel       = pset.getParameter<string>("resultFile");
	
	//  plotEvent = pset.getParameter<bool>("plotEvent");
	// zipFiles  = pset.getParameter<bool>("zipFiles");
}

/*****************************************************************************/
HighPtTrackAnalyzer::~HighPtTrackAnalyzer()
{
}

/*****************************************************************************/
void HighPtTrackAnalyzer::beginJob(const edm::EventSetup& es)
{
	// Get tracker geometry
	edm::ESHandle<TrackerGeometry> tracker;
	es.get<TrackerDigiGeometryRecord>().get(tracker);
	theTracker = tracker.product();
	
	// Get associator
	edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
	es.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",
										theHitsAssociator);
	theAssociatorByHits =
	(const TrackAssociatorByHits*)theHitsAssociator.product();
	
	// Get transient track builder
	edm::ESHandle<TransientTrackBuilder> builder;
	es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
	theTTBuilder = builder.product();
	
	// Root
	resultFile = new TFile(resultFileLabel.c_str(),"recreate");
	resultFile->cd();
	
	dMinSimPt=1.0;
	
	CAReco=new TClonesArray("MatchedTrack",10000);
	CASim=new TClonesArray("MatchedTrack",10000);
	
	recInfoTree=new TTree("RecoStudyTree","RecoStudyTree");
	recInfoTree->Branch("RunNo",&iRun,"RunNumber/I");
	recInfoTree->Branch("EventNo",&iEvent,"EventNumber/I");
	recInfoTree->Branch("RecTracks",&CAReco);
	recInfoTree->Branch("SimTracks",&CASim);
	recInfoTree->Branch("TotalRecoTracks",&iTrkReco,"TotalRecTracks/I");
	recInfoTree->Branch("TotalSimTracks",&iTrkSim,"TotalSimTracks/I");
	recInfoTree->Branch("TotalVtx",&iVtx,"TotalVtx/I");
}

/*****************************************************************************/
void HighPtTrackAnalyzer::endJob()
{
	resultFile->cd();
	recInfoTree->Write();
	resultFile->Close();
}

/*****************************************************************************/
int HighPtTrackAnalyzer::layerFromDetid(const DetId& detId)
{
	int layerNumber=0;
	unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());
	if ( subdetId == StripSubdetector::TIB)
    {
		TIBDetId tibid(detId.rawId());
		layerNumber = tibid.layer();
    }
	else if ( subdetId ==  StripSubdetector::TOB )
    {
		TOBDetId tobid(detId.rawId());
		layerNumber = tobid.layer();
    }
	else if ( subdetId ==  StripSubdetector::TID)
    {
		TIDDetId tidid(detId.rawId());
		layerNumber = tidid.wheel();
    }
	else if ( subdetId ==  StripSubdetector::TEC )
    {
		TECDetId tecid(detId.rawId());
		layerNumber = tecid.wheel();
    }
	else if ( subdetId ==  PixelSubdetector::PixelBarrel )
    {
		PXBDetId pxbid(detId.rawId());
		layerNumber = pxbid.layer();
    }
	else if ( subdetId ==  PixelSubdetector::PixelEndcap )
    {
		PXFDetId pxfid(detId.rawId());
		layerNumber = pxfid.disk();
    }
	else
		edm::LogVerbatim("TrackValidator") << "Unknown subdetid: " << subdetId;
	
	return layerNumber;
}

/*****************************************************************************/
int HighPtTrackAnalyzer::getNumberOfSimHits(const TrackingParticle& simTrack)
{
	int oldlay = 0; int newlay = 0;
	int olddet = 0; int newdet = 0;
	
	int nhit = 0;
	
	for(std::vector<PSimHit>::const_iterator
		simHit = simTrack.pSimHit_begin();
		simHit!= simTrack.pSimHit_end(); simHit++)
	{
		const DetId detId = DetId(simHit->detUnitId());
		oldlay = newlay; newlay = layerFromDetid(detId);
		olddet = newdet; newdet = detId.subdetId();
		if(oldlay != newlay || (oldlay == newlay && olddet != newdet) ) nhit++;
	}
	
	return nhit;
}

/*****************************************************************************/
int HighPtTrackAnalyzer::getDetLayerId(const PSimHit& simHit)
{
	int layerId;
	
	DetId id = DetId(simHit.detUnitId());
	LocalPoint lpos = simHit.localPosition();
	GlobalPoint gpos = theTracker->idToDetUnit(id)->toGlobal(lpos);
	
	if(theTracker->idToDetUnit(id)->subDetector() ==
	   GeomDetEnumerators::PixelBarrel)
	{ // barrel
		if(gpos.perp2() < 6 * 6) layerId = 0;
		else
		{
			if(gpos.perp2() < 9 * 9) layerId = 1;
			else layerId = 2;
		}
	}
	else
	{ // endcap
		if(fabsf(gpos.z()) < 40) layerId = 3;
		else layerId = 4;
	}
	
	return layerId;
}

/*****************************************************************************/
int HighPtTrackAnalyzer::getNumberOfPixelHits(const TrackingParticle& simTrack,float *fSimPxlLayerHit)
{
	// How many pixel hits?
	const int nLayers = 5;
	vector<bool> filled(nLayers,false);
	
	int numberOfPixelHits = 0;
	
	for(std::vector<PSimHit>::const_iterator simHit = simTrack.pSimHit_begin();simHit!= simTrack.pSimHit_end();simHit++){
		//DetId id = DetId(simHit->detUnitId());
		unsigned int id = simHit->detUnitId();

		if (id > 490000000){  //exclude simHits from outside the tracker
			cout << "-----------------id " << id << endl;
		} else {
			
			DetId detId(id);
			if(detId.subdetId() ==PixelSubdetector::PixelBarrel||detId.subdetId() ==PixelSubdetector::PixelEndcap){
				filled[getDetLayerId(*simHit)] = true;
				numberOfPixelHits++;
			}
		}
	}
	
	// Count the number of filled pixel layers
	int fLayers = 0;
	for(int i=0; i<nLayers; i++)
		if(filled[i] == true) fLayers++;
	
	*fSimPxlLayerHit=(float)fLayers;
	
	return numberOfPixelHits;
}



int HighPtTrackAnalyzer::getNumberOfRecPixelHits(const reco::Track& recTrack,float *fRecPxlLayerHit){
	
	const int nLayers = 5;
	vector<bool> filled(nLayers,false);
	
	int numberOfPixelHits = 0;
	
	for(trackingRecHit_iterator recHit = recTrack.recHitsBegin();recHit!= recTrack.recHitsEnd(); recHit++){
		if((*recHit)->isValid()){
			
			DetId id = (*recHit)->geographicalId();
			if(!theTracker->idToDet(id))
				continue;
			
			if(theTracker->idToDet(id)->subDetector() ==GeomDetEnumerators::PixelBarrel ||theTracker->idToDet(id)->subDetector() ==GeomDetEnumerators::PixelEndcap){
				int layerId;
				
				GlobalPoint gpos = theTracker->idToDet((*recHit)->geographicalId())->toGlobal((*recHit)->localPosition());
				
				if(theTracker->idToDet(id)->subDetector()==GeomDetEnumerators::PixelBarrel)
				{ // barrel
					if(gpos.perp2() < 6 * 6) layerId = 0;
					else
					{
						if(gpos.perp2() < 9 * 9) layerId = 1;
						else layerId = 2;
					}
				}
				else
				{ // endcap
					if(fabsf(gpos.z()) < 40) layerId = 3;
					else layerId = 4;
				}
				filled[layerId] = true;
				numberOfPixelHits++;
			}
		}
		
	}
	int fLayers = 0;
	for(int i=0; i<nLayers; i++)
		if(filled[i] == true) fLayers++;
	
	*fRecPxlLayerHit=(float)fLayers;
	
	return numberOfPixelHits;
}



/*****************************************************************************/
void HighPtTrackAnalyzer::checkSimTracks(edm::Handle<TrackingParticleCollection>& simCollection,reco::SimToRecoCollection& q){
	Int_t iSimCount=-1;
	TClonesArray &CASimTemp = *((TClonesArray*)CASim);
	
	//std::cout<<"entering"<<std::endl;
	
	for(TrackingParticleCollection::size_type i=0;i < simCollection.product()->size(); ++i){
		const TrackingParticleRef simTrack(simCollection, i);
		
		if(simTrack->charge() != 0){
			if(simTrack->pt()<2.0)
				continue;
			
			
			iSimCount++;
			MatchedTrack *TrackTemp=new(CASimTemp[iSimCount]) MatchedTrack();
			
			//std::cout<<"in the loop : "<<simTrack->pt()<<std::endl;
			
			// sim
			TrackTemp->iPID=simTrack->pdgId();
			//result.push_back(simTrack->parentVertex()->position().T()); // ?
			TrackTemp->fPt=simTrack->pt();
			TrackTemp->fPx=simTrack->px();
			TrackTemp->fPy=simTrack->py();
			TrackTemp->fPz=simTrack->pz();
			TrackTemp->fPhi=simTrack->phi();
			TrackTemp->fEta=simTrack->eta();
			TrackTemp->fPxlHit=getNumberOfPixelHits(*simTrack,(float *)&fSimPxlLayerHit);
			TrackTemp->fHit=getNumberOfSimHits(*simTrack);
			TrackTemp->fD0=simTrack->vertex().Rho();
			TrackTemp->iq=simTrack->charge();
			TrackTemp->fPxlLayerHit=fSimPxlLayerHit;
			
			
			// reco
			edm::RefToBase<reco::Track> matchedRecTrack;
			int nRec=0;
			int nSharedT=0;
			
			if(q.find(simTrack)!=q.end()){
				try{
					
					vector<pair<edm::RefToBase<reco::Track>, double> > recTracks = q[simTrack];
					for(vector<pair<edm::RefToBase<reco::Track>,double> >::const_iterator it = recTracks.begin(); it != recTracks.end(); ++it){
						edm::RefToBase<reco::Track> recTrack = it->first;
						int nShared=(int)(it->second * TrackTemp->fHit + 0.5);
						
						if(nSharedT<nShared){
							nSharedT=nShared;
							matchedRecTrack=recTrack; 
						}
						nRec++;
					}
				}
				catch (cms::Exception& event){
				}
				
				TrackTemp->iMatches=nRec;
				
				if(nSharedT > 0){
					TrackTemp->fMPt=matchedRecTrack->pt();
					TrackTemp->fMPx=matchedRecTrack->px();
					TrackTemp->fMPy=matchedRecTrack->py();
					TrackTemp->fMPz=matchedRecTrack->pz();
					TrackTemp->fMPhi=matchedRecTrack->phi();
					TrackTemp->fMEta=matchedRecTrack->eta();
					TrackTemp->fMD0=matchedRecTrack->d0();
					TrackTemp->fMD0Err=matchedRecTrack->d0Error();
					TrackTemp->iMq=matchedRecTrack->charge();
					TrackTemp->fMHit=matchedRecTrack->recHitsSize();
					TrackTemp->fMPxlHit=getNumberOfRecPixelHits(*matchedRecTrack,(float *)&fRecPxlLayerHit);
					TrackTemp->fMChi2=matchedRecTrack->chi2();
					TrackTemp->fMChi2Norm=matchedRecTrack->normalizedChi2();
					TrackTemp->fHitMatched=(Float_t)nSharedT;
					TrackTemp->fMZ=matchedRecTrack->dz();
					TrackTemp->fMZErr=matchedRecTrack->dzError();
					/* TrackTemp->fMZ=matchedRecTrack->dx();
					 TrackTemp->fMZErr=matchedRecTrack->dxError();
					 TrackTemp->fMZ=matchedRecTrack->dy();
					 TrackTemp->fMZErr=matchedRecTrack->dyError();*/
					TrackTemp->fMPxlLayerHit=fRecPxlLayerHit;
					TrackTemp->fMValidHits=matchedRecTrack->numberOfValidHits();
				}else{
					
				}
			}else{
				TrackTemp->iMatches=0;
			}
		}//end if charge!=0;
	}
}

/*****************************************************************************/
pair<float,float> HighPtTrackAnalyzer::refitWithVertex(const reco::Track & recTrack,const reco::VertexCollection* vertices){
	TransientTrack theTransientTrack = theTTBuilder->build(recTrack);
	
	// If there are vertices found
	if(vertices->size() > 0)
	{
		float dzmin = -1.;
		const reco::Vertex * closestVertex = 0;
		
		// Look for the closest vertex in z
		for(reco::VertexCollection::const_iterator
			vertex = vertices->begin(); vertex!= vertices->end(); vertex++)
		{
			float dz = fabs(recTrack.vertex().z() - vertex->position().z());
			if(vertex == vertices->begin() || dz < dzmin)
			{ dzmin = dz ; closestVertex = &(*vertex); }
		}
		
		// Get vertex position and error matrix
		GlobalPoint vertexPosition(closestVertex->position().x(),
								   closestVertex->position().y(),
								   closestVertex->position().z());
		
		float beamSize = 15e-4; // 15 um
		GlobalError vertexError(beamSize*beamSize, 0,
								beamSize*beamSize, 0,
								0,closestVertex->covariance(2,2));
		
		// Refit track with vertex constraint
		SingleTrackVertexConstraint stvc;
		pair<TransientTrack, float> result =
		stvc.constrain(theTransientTrack, vertexPosition, vertexError);
		
		return pair<float,float>(result.first.impactPointTSCP().pt(),
								 result.second);
	}
	else
		return pair<float,float>(recTrack.pt(), -9999);
	
}

/*****************************************************************************/
int HighPtTrackAnalyzer::getParticleId(edm::RefToBase<reco::Track>& recTrack, int & ptype)
{
	int pid = 0;
	ptype = 0;
	double tmin = 0.;
	
	for(trackingRecHit_iterator recHit = recTrack->recHitsBegin();
		recHit!= recTrack->recHitsEnd(); recHit++)
	{
		vector<PSimHit> simHits = theHitAssociator->associateHit(**recHit);
		
		for(vector<PSimHit>::const_iterator simHit = simHits.begin(); 
			simHit!= simHits.end(); simHit++)
			if(simHit == simHits.begin() || simHit->tof() < tmin )
			{
				pid   = simHit->particleType();
				ptype = simHit->processType();
				tmin  = simHit->tof();
			}
	}  
	
	return pid;
}

/*****************************************************************************/
void HighPtTrackAnalyzer::checkRecTracks(edm::Handle<edm::View<reco::Track> >& recCollection,const reco::VertexCollection* vertices,reco::RecoToSimCollection& p){
	
	Int_t iRecCount=-1;
	TClonesArray &CARecTemp = *((TClonesArray*)CAReco);
	
	for(edm::View<reco::Track> ::size_type i=0;i < recCollection.product()->size(); ++i){
		edm::RefToBase<reco::Track> recTrack(recCollection, i);
		iRecCount++;
		MatchedTrack *TrackTemp=new(CARecTemp[iRecCount]) MatchedTrack();
		
		
		TrackTemp->fPt=recTrack->pt();
		TrackTemp->fPx=recTrack->px();
		TrackTemp->fPy=recTrack->py();
		TrackTemp->fPz=recTrack->pz();
		TrackTemp->fPhi=recTrack->phi();
		TrackTemp->fEta=recTrack->eta();
		TrackTemp->fD0=recTrack->d0();
		TrackTemp->fD0Err=recTrack->d0Error();
		TrackTemp->iq=recTrack->charge();
		TrackTemp->fHit=recTrack->recHitsSize();
		TrackTemp->fPxlHit=getNumberOfRecPixelHits(*recTrack,(float *)&fRecPxlLayerHit);
		TrackTemp->fChi2=recTrack->chi2();
		TrackTemp->fChi2Norm=recTrack->normalizedChi2();
		TrackTemp->fZ=recTrack->dz();
		TrackTemp->fZErr=recTrack->dzError();
		/* TrackTemp->fMZ=matchedRecTrack->dx();
		 TrackTemp->fMZErr=matchedRecTrack->dxError();
		 TrackTemp->fMZ=matchedRecTrack->dy();
		 TrackTemp->fMZErr=matchedRecTrack->dyError();*/
		TrackTemp->fPxlLayerHit=fRecPxlLayerHit;
		TrackTemp->fRefitPt=refitWithVertex(*recTrack,vertices).first;
		TrackTemp->fRefitChi2=refitWithVertex(*recTrack,vertices).second;
		TrackTemp->fValidHits=recTrack->numberOfValidHits();
		
		
		// sim 
		TrackingParticleRef matchedSimTrack;
		int nSim = 0;
		
		Float_t fHitMatch=0.0;
		
		try{
			vector<pair<TrackingParticleRef, double> > simTracks = p[recTrack];
			
			for(vector<pair<TrackingParticleRef, double> >::const_iterator it = simTracks.begin(); it != simTracks.end(); ++it){
				TrackingParticleRef simTrack = it->first;
				float fraction = it->second;
				
				// If more than half is shared
				if(fraction>0.0){
					if(fraction>fHitMatch){ 
						matchedSimTrack = simTrack;
						fHitMatch=fraction;
					}
					nSim++;
				}
				
			}
		}catch (cms::Exception& event){ 
		}
		
		
		if(nSim > 0){
			int parentId;
			float T;
			
			int ptype;
			//int ids = getParticleId(recTrack, ptype); // EDIT: gave compiler warning
			getParticleId(recTrack, ptype);
			
			if(matchedSimTrack->parentVertex()->nSourceTracks() == 0){
				// track is primary, has no parent
				// recTrack can be a true primary, or an untracked daughter
				if(ptype == 2) parentId = 0;                        // primary
				else parentId = matchedSimTrack->pdgId(); // hadronic, decay
			}else{
				// track is not primary, has a parent
				TrackingVertex::tp_iterator iv =matchedSimTrack->parentVertex()->sourceTracks_begin();
				parentId = (*iv)->pdgId();
			}
			
			T = matchedSimTrack->parentVertex()->position().T(); // ?
			
			TrackTemp->iPID=matchedSimTrack->pdgId();
			TrackTemp->iMParentPID=parentId;
			//result.push_back(simTrack->parentVertex()->position().T()); // ?
			TrackTemp->fMPt=matchedSimTrack->pt();
			TrackTemp->fMPx=matchedSimTrack->px();
			TrackTemp->fMPy=matchedSimTrack->py();
			TrackTemp->fMPz=matchedSimTrack->pz();
			TrackTemp->fMPhi=matchedSimTrack->phi();
			TrackTemp->fMEta=matchedSimTrack->eta();
			TrackTemp->fMPxlHit=getNumberOfPixelHits(*matchedSimTrack,(float *)&fSimPxlLayerHit);
			TrackTemp->fMHit=getNumberOfSimHits(*matchedSimTrack);
			TrackTemp->fMD0=matchedSimTrack->vertex().Rho();
			TrackTemp->iMq=matchedSimTrack->charge();
			TrackTemp->fMPxlLayerHit=fSimPxlLayerHit;
			TrackTemp->fHitMatched=fHitMatch;
			TrackTemp->iMatches=nSim;
		}else {
		}  
	}
}

/*****************************************************************************/
void HighPtTrackAnalyzer::analyze(const edm::Event& ev, const edm::EventSetup& es){
	
	std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
	
	// Get associator
	theHitAssociator = new TrackerHitAssociator::TrackerHitAssociator(ev);
	
	// Get generated
	edm::Handle<edm::HepMCProduct> hepEv;
	ev.getByType(hepEv);
	proc = hepEv->GetEvent()->signal_process_id();
	cout<<"Event Number : "<<ev.id().event()<<endl;
	cout<<"[HighPtTrackAnalyzer] process = "<<proc<<endl;
	
	// Get simulated
	edm::Handle<TrackingParticleCollection> simCollection;
	//ev.getByLabel("trackingtruthprod",simCollection); //name changed in trackingParticles_cfi between 2_0_5 and 2_1_7
	ev.getByLabel("mergedtruth",simCollection);
	//  ev.getByType(simCollection);
	
	cout<<"[HighPtTrackAnalyzer] simTracks = "<<simCollection.product()->size()<<endl;
	
	// Get reconstructed
	edm::Handle<edm::View<reco::Track> >  recCollection;
	ev.getByLabel(trackCollectionLabels[0], recCollection); // !!
	
	cout<<"[HighPtTrackAnalyzer] recTracks = "<<recCollection.product()->size()<<endl;
	
	
	// Get vertices
	edm::Handle<reco::VertexCollection> vertexCollection;
	ev.getByLabel("pixelVertices",vertexCollection);
	const reco::VertexCollection * vertices = vertexCollection.product();
	
	iTrkReco = recCollection.product()->size();
	iTrkSim=simCollection.product()->size();
	iVtx = vertexCollection.product()->size();
	iEvent=ev.id().event();
	iRun=ev.id().run();
	
	CAReco->Clear();
	CASim->Clear();
	
	
	// Associators
	reco::SimToRecoCollection simToReco=theAssociatorByHits->associateSimToReco(recCollection, simCollection,&ev);
	
	
	reco::RecoToSimCollection recoToSim=theAssociatorByHits->associateRecoToSim(recCollection, simCollection,&ev);
	
	// Analyze
	checkSimTracks(simCollection,simToReco);
	
	checkRecTracks(recCollection, vertices, recoToSim);
	
	recInfoTree->Fill();
	
	
	cout<<"[HighPtTrackAnalyzer] done, "<<ev.id()<< endl;
	cerr<<"----------------------------------------------"<<endl;
	
	delete theHitAssociator;
}


DEFINE_FWK_MODULE(HighPtTrackAnalyzer);
