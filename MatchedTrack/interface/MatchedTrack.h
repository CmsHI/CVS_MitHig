#ifndef ROOT_MATCHEDTRACK_H
#define ROOT_MATCHEDTRACK_H

#include "TObject.h"
#include <vector>

using namespace std;

class MatchedTrack : public TObject {
	
public:
	MatchedTrack();
	virtual ~MatchedTrack();
	
	float fPt,fEta,fPhi;
	float fMPt,fMEta,fMPhi;
	float fPx,fPy,fPz;
	float fMPx,fMPy,fMPz;
	
	int iPID;
	int iParentPID;
	int iq,iMq;
	
	vector<int> PXBLayerHits;
	vector<int> PXFLayerHits;
	vector<int> TIBLayerHits;
	vector<int> TOBLayerHits;
	vector<int> TIDLayerHits;
	vector<int> TECLayerHits;
	
	vector<int> PXBLayerMHits;
	vector<int> PXFLayerMHits;
	vector<int> TIBLayerMHits;
	vector<int> TOBLayerMHits;
	vector<int> TIDLayerMHits;
	vector<int> TECLayerMHits;
	
	float fHitLayers;
	float fMHitLayers;

	float fD0,fMD0;
	float fD0Err,fMD0Err;
	
	float fHit;
	float fPxlHit;
	float fPxlLayerHit; 
	float fHitMatched;
	float fMPxlHit;
	float fMPxlLayerHit; 
	float fMHit;
	
	int iMatches;  
	
	float fChi2,fChi2Norm;
	float fMChi2,fMChi2Norm;
	float fZ,fX,fY;
	float fMZ,fMX,fMY;
	
	int iMPID;
	int iMParentPID;
	
	float fZErr,fXErr,fYErr;
	float fMZErr,fMXErr,fMYErr;
	float fRefitPt,fRefitChi2;
	float fMRefitPt,fMRefitChi2;
	float fValidHits,fMValidHits;
	
	ClassDef(MatchedTrack,4)
	
};

#endif
