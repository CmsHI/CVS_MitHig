#include "edwenger/MatchedTrack/interface/MatchedTrack.h"

MatchedTrack::MatchedTrack()
{
	fPt=0.0;
	fPx=0.0;
	fPy=0.0;
	fPz=0.0;
	fPhi=0.0;
	fEta=0.0;
	fD0=0.0;
	iq=0;
	fHit=0;
	fChi2=0;
	fChi2Norm=0;
	fMPt=0.0;
	fMPx=0.0;
	fMPy=0.0;
	fMPz=0.0;
	fMPhi=0.0;
	fMEta=0.0;
	fMD0=0.0;
	iMq=0;
	fMHit=0;
	fMChi2=0;
	fMChi2Norm=0;
	fZ=0;fY=0;fX=0;
	fMZ=0;fMY=0;fMX=0;
	iMatches=0;
	iPID=0;
	iParentPID=0;
	fPxlHit=0.0;
	fPxlLayerHit=0.0; 
	fHitMatched=0.0;
	fD0Err=0.0;fMD0Err=0.0;
	
}

MatchedTrack::~MatchedTrack()
{
	Clear();	
}

ClassImp(MatchedTrack)
