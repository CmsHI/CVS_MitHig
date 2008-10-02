#include "MitHig/PixelTrackletAnalyzer/interface/TrackletFinder.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "TTree.h"

using namespace edm;
using namespace std;

namespace {
bool compareTracklet(Tracklet a,Tracklet b) {    return fabs(a.dR2())<fabs(b.dR2()); }
}

void TrackletFinder::setEvent(const edm::Event & iEvent){

  event_ = &iEvent;
  layer1_.clear();
  layer2_.clear();

}


void TrackletFinder::sortLayers(){
  
  const SiPixelRecHitCollection* rechits;
  Handle<SiPixelRecHitCollection> rchts;
  event_->getByLabel("siPixelRecHits",rchts);
  rechits = rchts.product();

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
	    }
	  
	  for(SiPixelRecHitCollection::const_iterator recHit = range.first; recHit!= range.second; recHit++)
	    {
	      if(layer == 1) layer1_.push_back(&(*recHit));
	      if(layer == 2) layer2_.push_back(&(*recHit));
	    }
	}
    }

}

std::vector<Tracklet> TrackletFinder::makeTracklets(bool invert){
  
  vector<Tracklet> recoTracklets;
  TrackerHitAssociator theHitAssociator(*event_);
    
  for(unsigned int i1 = 0; i1 < layer1_.size(); ++i1){   //loops over and gets spatial information and associated simhits for each rechit
    
    // Ids
    int rechit1Type = 0;
    int rechit2Type = 0;
    int trackletType = -1;
    int signalExistCheck = 0;
    
    const SiPixelRecHit* recHit1 = layer1_[i1];
    
    const PixelGeomDetUnit* pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit1->geographicalId()));
    
    GlobalPoint gpos1 = pixelLayer->toGlobal(recHit1->localPosition());
    
    // Calculate the rechit position with respect to the vertex
    math::XYZVector rechitPos(gpos1.x(),gpos1.y(),gpos1.z()-vertex_.z());
    double phi1 = rechitPos.phi();
    double eta1 = rechitPos.eta();
    
    // Match with second layer reconstructed hits
    for(unsigned int i2 = 0; i2 < layer2_.size(); ++i2){
      const SiPixelRecHit* recHit2 = layer2_[i2];
      pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit2->geographicalId()));
      GlobalPoint gpos2 = pixelLayer->toGlobal(recHit2->localPosition());
	  
      // Calculate the rechit position with respect to the vertex
      math::XYZVector rechit2Pos(gpos2.x(),gpos2.y(),gpos2.z()-vertex_.z());
      double phi2 = rechit2Pos.phi();
      double eta2 = rechit2Pos.eta();
      if (invert==1) eta2*=-1;
      
      Tracklet mytracklet(eta1,eta2,phi1,phi2);
      
      mytracklet.setIt1(i1);
      mytracklet.setIt2(i2);
      mytracklet.setId1(rechit1Type);
      mytracklet.setId2(rechit2Type);
      mytracklet.setId(trackletType);
      mytracklet.setSId(signalExistCheck);

      if(invert) mytracklet.setType(-1);
      else mytracklet.setType(1);

      recoTracklets.push_back(mytracklet);
      
    }
    
  }
  
  return recoTracklets;
  
}



std::vector<Tracklet> TrackletFinder::cleanTracklets(vector<Tracklet> input){

  vector<Tracklet> output;
  sort( input.begin() , input.end() , compareTracklet);

  int used1[1000];
  int used2[1000];

  for (int i=0;i<1000;i++) { 
    used1[i]=0;
    used2[i]=0;
  } 

  for (unsigned int i = 0; i < input.size(); i++){
      int i1=input[i].getIt1();
      int i2=input[i].getIt2();
      if (used1[i1]==0&&used2[i2]==0) {
	Tracklet tmp = input[i];
	output.push_back(tmp);
	used1[i1]=1;
	if (checkSecondLayer_) used2[i2]=1;
      }
  }

  return output;

}

void TrackletFinder::fillPixelEvent(TNtuple* nt1, TNtuple* nt2){

  TrackerHitAssociator theHitAssociator(*event_);

  for(unsigned int i1 = 0; i1 < layer1_.size(); ++i1){   //loops over and gets spatial information and associated simhits for each rechit

    // Ids
    int rechit1Type = 0;
    int rechit2Type = 0;
    int trackletType = -1;
    int signalExistCheck = 0;

    const SiPixelRecHit* recHit1 = layer1_[i1];

    const PixelGeomDetUnit* pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit1->geographicalId()));

    GlobalPoint gpos1 = pixelLayer->toGlobal(recHit1->localPosition());

    // Calculate the rechit position with respect to the vertex
    //    math::XYZVector rechitPos(gpos1.x(),gpos1.y(),gpos1.z()-vertex_.z());

    nt1->Fill(gpos1.x(),gpos1.y(),gpos1.z());

    math::XYZVector rechitPos(gpos1.x(),gpos1.y(),gpos1.z());
    double phi1 = rechitPos.phi();
    double eta1 = rechitPos.eta();

  }

    // Match with second layer reconstructed hits
    for(unsigned int i2 = 0; i2 < layer2_.size(); ++i2){
      const SiPixelRecHit* recHit2 = layer2_[i2];

      const PixelGeomDetUnit* pixelLayer = dynamic_cast<const PixelGeomDetUnit*> (geo_->idToDet(recHit2->geographicalId()));

      GlobalPoint gpos2 = pixelLayer->toGlobal(recHit2->localPosition());

      // Calculate the rechit position with respect to the vertex
      //      math::XYZVector rechit2Pos(gpos2.x(),gpos2.y(),gpos2.z()-vertex_.z());

      nt1->Fill(gpos2.x(),gpos2.y(),gpos2.z());

      math::XYZVector rechit2Pos(gpos2.x(),gpos2.y(),gpos2.z()-vertex_.z());
      double phi2 = rechit2Pos.phi();
      double eta2 = rechit2Pos.eta();
    }


}











