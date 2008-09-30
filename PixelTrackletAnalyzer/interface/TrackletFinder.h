#ifndef PIXELTRACKLETANA_FINDER
#define PIXELTRACKLETANA_FINDER

#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MitHig/PixelTracklet/interface/TrackletCorrections.h"
#include "MitHig/PixelTracklet/interface/Tracklet.h"

class TrackletFinder{
public :

  TrackletFinder(
		 TrackletCorrections* corr, const TrackerGeometry* geo, bool count = false, bool checkSecondLayer = false
		 ): 
  corr_(corr), 
  geo_(geo), 
  countMode_(count), 
  checkSecondLayer_(checkSecondLayer) 
  {;}
  void sortLayers();
  std::vector<Tracklet> makeTracklets(bool invert = false);
  std::vector<Tracklet> cleanTracklets(std::vector<Tracklet> input);
  std::vector<Tracklet> getTracklets() {return cleanTracklets(makeTracklets());}
  void setEvent(const edm::Event & iEvent);
  //  bool compareTracklet(Tracklet a,Tracklet b);
  
private:

  const edm::Event* event_;
  std::vector<const SiPixelRecHit*> layer1_;
  std::vector<const SiPixelRecHit*> layer2_;
  math::XYZVector vertex_;
  TrackletCorrections* corr_;  
  const TrackerGeometry* geo_;
  bool countMode_;
  bool checkSecondLayer_;

};

//bool compareTracklet(Tracklet a,Tracklet b) {    return fabs(a.dR2())<fabs(b.dR2()); }

#endif
