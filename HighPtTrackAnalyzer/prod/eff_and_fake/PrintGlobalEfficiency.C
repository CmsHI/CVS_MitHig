


Float_t fHitsCutOff=8;
Float_t fPxlLayersCutOff=3;
Float_t fPtLowCutOff=2.0;
Float_t fEtaCutOff=0.7;
Float_t fHitMatchCut=0.0;
Float_t fb=0;


Float_t fPadX=800,fPadY=500;

PrintGlobalEfficiency(Float_t fPtLowCutOffT=1.0,Float_t fEtaCutOffT=0.5,Float_t fHitMatchCutT=0.80,Float_t fbT=0,Float_t fHitsCutOffT=8,Float_t fPxlLayersCutOffT=3){

  fHitsCutOff=fHitsCutOffT;
  fPxlLayersCutOff=fPxlLayersCutOffT;
  fb=fbT;
  fPtLowCutOff=fPtLowCutOffT;
  fEtaCutOff=fEtaCutOffT;
  fHitMatchCut=fHitMatchCutT;


  Char_t cTemp[200];
  
  sprintf(cTemp,"GE_b%2.1f_PtMin%2.2f_EtaCutOff%2.2f_HitMch%2.2f_Hit%d_PxlLayers%1.0f.root",fb,fPtLowCutOff,fEtaCutOff,fHitMatchCut,fHitsCutOff,fPxlLayersCutOff);
  TFile *fSummary=new TFile(cTemp);
  fSummary->ls();


  // hSimPt->Rebin(2);
  //hSimPtAcc->Rebin(2);

  //hSimPxlLyPt->Rebin(2);
  //hSimPxlLyPtAcc->Rebin(2);
 

  TGraphErrors *gSimPtAll = new TGraphErrors();
  TGraphErrors *gSimPtPxlLy = new TGraphErrors();
  TGraphErrors *gSimPtAllReco = new TGraphErrors();
 
  Float_t fRangeMin=1.0;
  Float_t fRangeMax=2.0;
  Float_t fRangeIncr=1.0;

  Double_t dAllTracks,dAllTracksAcc;
  
  Double_t dPxlLyTracks,dPxlLyTracksAcc;

  Double_t dAllR,dPxlLyR,dAllRecoR;
  
  Double_t dAllRErr,dPxlLyRErr,dAllRecoRErr;
  for(Int_t i=0;i<16;i++){

    sprintf(cTemp,"   %03d - %03d   ",fRangeMin,fRangeMax);
    cout<<cTemp;
    
    dAllTracks=hSimPt->Integral(fRangeMin,fRangeMax);
    dAllTracksAcc=hSimPtAcc->Integral(fRangeMin,fRangeMax);
    
    dPxlLyTracks=hSimHitsPt->Integral(fRangeMin,fRangeMax);
    dPxlLyTracksAcc=hSimHitsPtAcc->Integral(fRangeMin,fRangeMax);

    dAllR=0.0;
    dPxlLyR=0.0;
    dAllRecoR=0.0;
   
    dAllRErr=0.0;
    dPxlLyRErr=0.0;
    dAllRecoRErr=0.0;                                                     ;

    if(dAllTracks>0.0)
      dAllR=dAllTracksAcc/dAllTracks;
    if(dPxlLyTracks>0.0)
      dPxlLyR=dPxlLyTracksAcc/dPxlLyTracks;
    if(dAllTracks>0.0)
      dAllRecoR=dPxlLyTracks/dAllTracks;
      
    if(dAllTracks>0.0&&dAllTracksAcc>0.0)
      dAllRErr=TMath::Sqrt((1./dAllTracksAcc)**2.0+(1.0/dAllTracks)**2.0);
    if(dPxlLyTracks>0.0&&dPxlLyTracksAcc>0.0)
      dPxlLyRErr=TMath::Sqrt((1./dPxlLyTracksAcc)**2.0+(1.0/dPxlLyTracks)**2.0);
    if(dAllTracks>0.0&&dPxlLyTracks>0.0)
      dAllRecoRErr=TMath::Sqrt((1./dPxlLyTracks)**2.0+(1.0/dAllTracks)**2.0);
      

    gSimPtAll->SetPoint(i,(fRangeMax-fRangeMin)/2.0,dAllR);
    gSimPtAll->SetPointError(i,0.0,dAllRErr);

    gSimPtPxlLy->SetPoint(i,(fRangeMax-fRangeMin)/2.0,dPxlLyR);
    gSimPtPxlLy->SetPointError(i,0.0,dPxlLyRErr);

    gSimPtAllReco->SetPoint(i,(fRangeMax-fRangeMin)/2.0,dAllRecoR);
    gSimPtAllReco->SetPointError(i,0.0,dAllRecoRErr);

    sprintf(cTemp,"   %1.6f +/- %1.6f   ",dAllRecoR,dAllRecoRErr);
    cout<<cTemp;

    sprintf(cTemp,"   %1.6f +/- %1.6f   ",dAllR,dAllRErr);
    cout<<cTemp;
    
    sprintf(cTemp,"   %1.6f +/- %1.6f   ",dPxlLyR,dPxlLyRErr);
    cout<<cTemp;
    
    cout<<endl;

    fRangeMin=fRangeMin+fRangeIncr;
    if(fRangeMax>=5.){
      fRangeIncr=5.0;
      if(fRangeMax>=30.){
	fRangeIncr=10.0;
      }
      if(fRangeMax>=70.){
	fRangeIncr=20.0;
      }
      if(fRangeMax>=110.){
	fRangeIncr=90.0;
      }
    }
    fRangeMax=fRangeMax+fRangeIncr;
  }

  

  

  fSummary->Close();


}
