TH2F *hSimEtaPt,*hSimEtaPtEff,*hSimEtaPtAcc;
TH2F *hRecEtaPt,*hRecEtaPtFR;

TH3F *hMomResEtaPt;

Float_t fPtMin=0.0,fPtMax=200;
Int_t iPtBins=200;

Float_t fPtMinRes=-50,fPtMaxRes=50;
Int_t iPtBinsRes=400;


Float_t fEtaMin=-13,fEtaMax=13;
Int_t iEtaBins=261;


Float_t fRadMin=0.0,fRadMax=15;
Int_t iRadBins=150;

Float_t fRadTemp;

Float_t fb=0;
Float_t fProbCutOff;
Float_t fDCACutOff;
Float_t fHitsCutOff;
Float_t fHitMatchCutOff;

Float_t fSimHitsCutOff=0;
Float_t fSimPxlLayersCutOff=3;


Float_t fPtMaxT,fPhiMaxT,fEtaMaxT;

MakeTrackingEfficiency(Float_t fbT=0,Float_t fProbCutOffT=0.01,Float_t fDCACutOffT=3.0,Float_t fHitsCutOffT=12,Float_t fHitMatchCutOffT=0.75){

  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMitHigMatchedTrack.so");

  fb=fbT;
  fProbCutOff=fProbCutOffT;
  fDCACutOff=fDCACutOffT;
  fHitsCutOff=fHitsCutOffT;
  fHitMatchCutOff=fHitMatchCutOffT;


  hSimEtaPt=new TH2F("hSimEtaPt","",iEtaBins,fEtaMin,fEtaMax,iPtBins,fPtMin,fPtMax); 
  hSimEtaPtEff=new TH2F("hSimEtaPtEff","",iEtaBins,fEtaMin,fEtaMax,iPtBins,fPtMin,fPtMax);
  hSimEtaPtAcc=new TH2F("hSimEtaPtAcc","",iEtaBins,fEtaMin,fEtaMax,iPtBins,fPtMin,fPtMax);
  hMomResEtaPt=new TH3F("hMomResEtaPt","",iEtaBins,fEtaMin,fEtaMax,iPtBins,fPtMin,fPtMax,iPtBinsRes,fPtMinRes,fPtMaxRes);

  
  hRecEtaPt=new TH2F("hRecEtaPt","",iEtaBins,fEtaMin,fEtaMax,iPtBins,fPtMin,fPtMax); 
  hRecEtaPtFR=new TH2F("hRecEtaPtFR","",iEtaBins,fEtaMin,fEtaMax,iPtBins,fPtMin,fPtMax);
  

  Int_t iEvent,iRun;

  Int_t iTracks;
  Double_t dphi;
  
  Char_t cFileName[256];
  ifstream fNames;

  TClonesArray *CAReco,*CASim;
  CAReco=new TClonesArray("MatchedTrack");
  CASim=new TClonesArray("MatchedTrack");

  fNames.open("FileNames.txt");
  while(!fNames.eof()){
    fNames>>cFileName;
    cout<<cFileName<<endl;
    TChain *chain=new TChain("RecoStudyTree");
    chain->Add(cFileName);
   
    chain->SetBranchAddress("RunNo",&iRun);
    chain->SetBranchAddress("EventNo",&iEvent);
    chain->SetBranchAddress("RecTracks",&CAReco);
    chain->SetBranchAddress("SimTracks",&CASim);
    
    Float_t fPtMin=0.0;
    fPtMaxT=0.0;
    
    cout<<"Entries "<<chain->GetEntries()<<endl;    
    
    for(Long_t i=0;i<chain->GetEntries();i++){
      chain->GetEntry(i);

      Int_t iRecPart=CAReco->GetEntriesFast();
      for(Int_t j=0;j<iRecPart;j++){
	MatchedTrack *recTempTrack=(MatchedTrack*)CAReco->At(j);
	Float_t fProb=0.00;
	if(recTempTrack->fChi2Norm>0.0)
	  fProb=TMath::Prob(recTempTrack->fChi2,recTempTrack->fChi2/recTempTrack->fChi2Norm);
	if(fProb>fProbCutOffT&&TMath::Abs(recTempTrack->fD0/recTempTrack->fD0Err)<fDCACutOffT&&recTempTrack->fHit>fHitsCutOffT){
	  hRecEtaPt->Fill(recTempTrack->fEta,recTempTrack->fPt);
	  if(recTempTrack->fHitMatched<fHitMatchCutOff){
	    hRecEtaPtFR->Fill(recTempTrack->fEta,recTempTrack->fPt);
	  }else{
	    hMomResEtaPt->Fill(recTempTrack->fEta,recTempTrack->fPt,recTempTrack->fPt-recTempTrack->fMPt);
	  }
	}
      }

      Int_t iSimPart=CASim->GetEntriesFast();
      for(Int_t j=0;j<iSimPart;j++){
	MatchedTrack *simTempTrack=(MatchedTrack*)CASim->At(j);
       	if(simTempTrack->fD0>0.01)
	  continue;
	if(simTempTrack->iq==0)
	  continue;
	/*	if(simTempTrack->fHit==0)
		continue;*/
	/*	if(TMath::Abs(simTempTrack->iPID)!=211)
		continue;*/
	

	hSimEtaPt->Fill(simTempTrack->fEta,simTempTrack->fPt);
	if(simTempTrack->fHit>=fSimHitsCutOff&&simTempTrack->fPxlLayerHit>=fSimPxlLayersCutOff){	
	  hSimEtaPtAcc->Fill(simTempTrack->fEta,simTempTrack->fPt);
	}

	if(simTempTrack->fHit>=fSimHitsCutOff&&simTempTrack->fPxlLayerHit>=fSimPxlLayersCutOff){
	  if(simTempTrack->iMatches>0){
	    Float_t fProb=0.00;
	    /*if(TMath::Abs(simTempTrack->fMPxlHit-simTempTrack->fPxlHit)<=2)
	      hSimEtaPtEff->Fill(simTempTrack->fEta,simTempTrack->fPt);

	    if((Float_t)simTempTrack->fMPxlHit/(Float_t)simTempTrack->fPxlHit>0.9)
	       hSimEtaPtEff->Fill(simTempTrack->fEta,simTempTrack->fPt);
	    */
	      
	    if(simTempTrack->fMChi2Norm>0.0)
	      fProb=TMath::Prob(simTempTrack->fMChi2,simTempTrack->fMChi2/simTempTrack->fMChi2Norm);
	    
	    if(fProb>fProbCutOffT&&simTempTrack->fMHit>fHitsCutOffT&&TMath::Abs(simTempTrack->fMD0/simTempTrack->fMD0Err)<fDCACutOffT){
	       hSimEtaPtEff->Fill(simTempTrack->fEta,simTempTrack->fPt);
	    
	    }  
	    /* if(simTempTrack->fMHit>fHitsCutOffT&&TMath::Abs(simTempTrack->fMD0/simTempTrack->fMD0Err)<fDCACutOffT){
	       hSimEtaPtEff->Fill(simTempTrack->fEta,simTempTrack->fPt);}*/
	  }
	}
      }

    }
  
    delete chain;
  }  
  Char_t cTemp[200];
  
  sprintf(cTemp,"TrkEff_b%2.1f_Prob%0.5f_DCA%1.2f_Hits%f_Match%0.2f.root",fb,fProbCutOff,fDCACutOff,fHitsCutOff,fHitMatchCutOff);
  TFile *fSummary=new TFile(cTemp,"recreate");
  
  hSimEtaPt->Write();
  hSimEtaPtEff->Write();
  hSimEtaPtAcc->Write();   
  
  
  hRecEtaPt->Write();
  hRecEtaPtFR->Write();  


  hMomResEtaPt->Write();
  
  fSummary->Close();
  
  
  
  
}

