

TH1F *hSimPt,*hSimEta;
TH1F *hSimPtAcc,*hSimEtaAcc;

TH1F *hSimPxlLyPt,SimPxlLyEta;
TH1F *hSimPxlLyPtAcc,*hSimPxlLyEtaAcc;


TH1F *hSimHitsPt,*hSimHitsEta;
TH1F *hSimHitsPtAcc,*hSimHitsEtaAcc;

TH1F *hSimPxlLyHitsPt,*hSimPxlLyHitsEta;
TH1F *hSimPxlLyHitsPtAcc,*hSimPxlLyHitsEtaAcc;



Float_t fPtMin=0.0,fPtMax=200;
Int_t iPtBins=200;

Float_t fEtaMin=-13,fEtaMax=13;
Int_t iEtaBins=130;

Float_t fHitsCutOff=8;
Float_t fPxlLayersCutOff=3;
Float_t fPtLowCutOff=2.0;
Float_t fEtaCutOff=0.7;
Float_t fHitMatchCut=0.0;
Float_t fb=0;


Bool_t dBug=kFALSE;

MakeGlobalEfficiecny(Float_t fPtLowCutOffT=1.0,Float_t fEtaCutOffT=0.5,Float_t fHitMatchCutT=0.80,Float_t fbT=0,Float_t fHitsCutOffT=8,Float_t fPxlLayersCutOffT=3){


  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMitHigMatchedTrack.so");

  fHitsCutOff=fHitsCutOffT;
  fPxlLayersCutOff=fPxlLayersCutOffT;
  fb=fbT;
  fPtLowCutOff=fPtLowCutOffT;
  fEtaCutOff=fEtaCutOffT;
  fHitMatchCut=fHitMatchCutT;

  hSimPt=new TH1F("hSimPt","",iPtBins,fPtMin,fPtMax); 
  hSimEta=new TH1F("hSimEta","",iEtaBins,fEtaMin,fEtaMax);
  hSimPtAcc=new TH1F("hSimPtAcc","",iPtBins,fPtMin,fPtMax); 
  hSimEtaAcc=new TH1F("hSimEtaAcc","",iEtaBins,fEtaMin,fEtaMax);

  hSimPxlLyPt=new TH1F("hSimPxlLyPt","",iPtBins,fPtMin,fPtMax); 
  hSimPxlLyEta=new TH1F("hSimPxlLyEta","",iEtaBins,fEtaMin,fEtaMax);
  hSimPxlLyPtAcc=new TH1F("hSimPxlLyPtAcc","",iPtBins,fPtMin,fPtMax); 
  hSimPxlLyEtaAcc=new TH1F("hSimPxlLyEtaAcc","",iEtaBins,fEtaMin,fEtaMax);
  // hSimPxlLyPtReco=new TH1F("hSimPxlLyPtReco","",iPtBins,fPtMin,fPtMax); 
  //  hSimPxlLyEtaReco=new TH1F("hSimPxlLyEtaReco","",iEtaBins,fEtaMin,fEtaMax);


  hSimHitsPt=new TH1F("hSimHitsPt","",iPtBins,fPtMin,fPtMax); 
  hSimHitsEta=new TH1F("hSimHitsEta","",iEtaBins,fEtaMin,fEtaMax);
  hSimHitsPtAcc=new TH1F("hSimHitsPtAcc","",iPtBins,fPtMin,fPtMax); 
  hSimHitsEtaAcc=new TH1F("hSimHitsEtaAcc","",iEtaBins,fEtaMin,fEtaMax);



  hSimPxlLyHitsPt=new TH1F("hSimPxlLyHitsPt","",iPtBins,fPtMin,fPtMax); 
  hSimPxlLyHitsEta=new TH1F("hSimPxlLyHitsEta","",iEtaBins,fEtaMin,fEtaMax);
  hSimPxlLyHitsPtAcc=new TH1F("hSimPxlLyHitsPtAcc","",iPtBins,fPtMin,fPtMax); 
  hSimPxlLyHitsEtaAcc=new TH1F("hSimEtaPxlLyHitsEtaAcc","",iEtaBins,fEtaMin,fEtaMax);


  if(dBug)
    cout<<"1"<<endl;


  

  TChain *chain=new TChain("RecoStudyTree");  
  chain->Add("RecoStudyOutput.root");


  TClonesArray *CAReco,*CASim;
  CAReco=new TClonesArray("MatchedTrack");
  CASim=new TClonesArray("MatchedTrack");

  Int_t iTrkSim,iTrkReco;
  Int_t iEvent,iRun;
  
  
  chain->SetBranchAddress("RunNo",&iRun);
  chain->SetBranchAddress("EventNo",&iEvent);
  chain->SetBranchAddress("RecTracks",&CAReco);
  chain->SetBranchAddress("SimTracks",&CASim);
  
  

  for(Long_t i=0;i<chain->GetEntries();i++){
    chain->GetEntry(i);
    if(i%20==0)
      cout<<endl;
    cout<<i<<" : ";
    Int_t iSimPart=CASim->GetEntriesFast();
    for(Int_t j=0;j<iSimPart;j++){
      MatchedTrack *simTempTrack=(MatchedTrack*)CASim->At(j);
      
      if(simTempTrack->fPt>fPtLowCutOff){
	hSimEta->Fill(simTempTrack->fEta);
	
	if(simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	  hSimPxlLyEta->Fill(simTempTrack->fEta);
	}
	
	if(simTempTrack->fHit>=fHitsCutOff){
	  hSimHitsEta->Fill(simTempTrack->fEta);
	}
	
	if(simTempTrack->fHit>=fHitsCutOff&&simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	  hSimPxlLyHitsEta->Fill(simTempTrack->fEta);
	}
	
	if(simTempTrack->iMatches>0){
	  
	  hSimEtaAcc->Fill(simTempTrack->fEta);
	 
	  if(simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	    hSimPxlLyEtaAcc->Fill(simTempTrack->fEta);
	  }
	 
	  if(simTempTrack->fHit>=fHitsCutOff){
	    hSimHitsEtaAcc->Fill(simTempTrack->fEta);
	  }
	 
	  if(simTempTrack->fHit>=fHitsCutOff&&simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	    hSimPxlLyHitsEtaAcc->Fill(simTempTrack->fEta);
	  }
	 
	}
      }
      

      if(TMath::Abs(simTempTrack->fEta)<fEtaCutOff&&TMath::Abs(simTempTrack->fEta)>fEtaCutOff-0.5){
	hSimPt->Fill(simTempTrack->fPt);
	
	if(simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	  hSimPxlLyPt->Fill(simTempTrack->fPt);
	}
	
	if(simTempTrack->fHit>=fHitsCutOff){
	  hSimHitsPt->Fill(simTempTrack->fPt);
	}
	
	if(simTempTrack->fHit>=fHitsCutOff&&simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	  hSimPxlLyHitsPt->Fill(simTempTrack->fPt);
	}

	if(simTempTrack->iMatches>0){
	  hSimPtAcc->Fill(simTempTrack->fPt);

	  if(simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	    hSimPxlLyPtAcc->Fill(simTempTrack->fPt);
	  }
	  if(simTempTrack->fHit>=fHitsCutOff){
	    hSimHitsPtAcc->Fill(simTempTrack->fPt);
	  }
	  
	  if(simTempTrack->fHit>=fHitsCutOff&&simTempTrack->fPxlLayerHit>=fPxlLayersCutOff){
	    hSimPxlLyHitsPtAcc->Fill(simTempTrack->fPt);
	  }
    
	}
      }
      //      delete simTempTrack;
    }
    //  cout<<endl;
  }
  
  Char_t cTemp[200];

  sprintf(cTemp,"GE_b%2.1f_PtMin%2.2f_EtaCutOff%2.2f_HitMch%2.2f_Hit%d_PxlLayers%1.0f.root",fb,fPtLowCutOff,fEtaCutOff,fHitMatchCut,fHitsCutOff,fPxlLayersCutOff);
  TFile *fSummary=new TFile(cTemp,"recreate");
   
  hSimPt->Write();
  hSimEta->Write();
  hSimPtAcc->Write();
  hSimEtaAcc->Write();

  hSimPxlLyPt->Write();
  hSimPxlLyEta->Write();
  hSimPxlLyPtAcc->Write();
  hSimPxlLyEtaAcc->Write();

  hSimHitsPt->Write();
  hSimHitsEta->Write();
  hSimHitsPtAcc->Write();
  hSimHitsEtaAcc->Write();
  
  hSimPxlLyHitsPt->Write();
  hSimPxlLyHitsEta->Write();
  hSimPxlLyHitsPtAcc->Write();
  hSimPxlLyHitsEtaAcc->Write();
  

  fSummary->Close();
}
