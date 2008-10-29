

Float_t fb=0;
Float_t fProbCutOff;
Float_t fDCACutOff;
Float_t fHitsCutOff;
Float_t fHitMatchCutOff;

Float_t fSimHitsCutOff=0;
Float_t fSimPxlLayersCutOff=3;


Float_t fPtMaxT,fPhiMaxT,fEtaMaxT;

Float_t fPtBins[]={3.,7.0,10.,15.,20.,25.,35.,45.,75.,200.};
Int_t iPtBins=10;

//Float_t fPtBins[]={7.0,13.0,18.0,28.0,38.0,60.0,100.0,150.,200};
//Int_t iPtBins=9;

TCanvas *EventProfile;

Float_t fPadX=350,fPadY=350;


TH1F *hDenominator,*hNumerator;


PrintTrackingEfficiency(Float_t fbT=0,Float_t fProbCutOffT=0.01,Float_t fDCACutOffT=3.0,Float_t fHitsCutOffT=12,Float_t fHitMatchCutOffT=0.75){
  //  gROOT->ProcessLine(".x ~/rootlogon.C");
  gStyle->SetErrorX(0.);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();

  fb=fbT;
  fProbCutOff=fProbCutOffT;
  fDCACutOff=fDCACutOffT;
  fHitsCutOff=fHitsCutOffT;
  fHitMatchCutOff=fHitMatchCutOffT;

  Char_t cTemp[200];   
  sprintf(cTemp,"TrkEff_b%2.1f_Prob%0.5f_DCA%1.2f_Hits%f_Match%0.2f.root",fb,fProbCutOff,fDCACutOff,fHitsCutOff,fHitMatchCutOff);
  /*  cout<<GetPtTrackingEff(cTemp,1.0,2.0,0.0,0.5)<<" : "<<GetPtTrackingEffErr(cTemp,1.0,2.0,0.0,0.5)<<endl;
  cout<<GetPtAlgorithmicEff(cTemp,1.0,2.0,0.0,0.5)<<" : "<<GetPtAlgorithmicEffErr(cTemp,1.0,2.0,0.0,0.5)<<endl;
  cout<<GetPtGeometricAcc(cTemp,1.0,2.0,0.0,0.5)<<" : "<<GetPtGeometricAccErr(cTemp,1.0,2.0,0.0,0.5)<<endl;
  cout<<GetPtGlobalTrackingEff(cTemp,1.0,2.0,0.0,0.5)<<" : "<<GetPtGlobalTrackingEffErr(cTemp,1.0,2.0,0.0,0.5)<<endl;
  cout<<GetPtGlobalTrackingEffWithoutFakeRate(cTemp,1.0,2.0,0.0,0.5)<<" : "<<GetPtGlobalTrackingEffErrWithoutFakeRate(cTemp,1.0,2.0,0.0,0.5)<<endl;
  */
  EventProfile = new TCanvas("EventProfile", "",0,0,(Int_t)fPadX,(Int_t)fPadY);
  EventProfile->SetHighLightColor(0);
  EventProfile->SetFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  EventProfile->SetBorderMode(0);
  EventProfile->SetBorderSize(0);
  EventProfile->SetFrameBorderMode(0);
  EventProfile->SetFrameBorderSize(0);

  Int_t iIndex=1;

  for(Float_t fEta=0.0;fEta<2.5;fEta+=1.0){
   Float_t fPtMin=2.0;
   cout<<"Eta bin "<<fEta<<" : "<<fEta+0.5<<endl;
   EventProfile->cd();
   TPad *padMainPlot=new TPad("MainPlot","",0.0,0,1.0,1.0);
   padMainPlot->SetRightMargin(0.03);
   padMainPlot->SetBottomMargin(0.15);
   padMainPlot->SetLeftMargin(0.15);
   padMainPlot->SetTopMargin(0.03);
   padMainPlot->SetBorderMode(0);
   padMainPlot->SetFillColor(0);
   padMainPlot->SetBorderSize(2);
   
  
   padMainPlot->SetTicks();
   padMainPlot->Draw();
   padMainPlot->cd();
   
   TH1F *hFrameMain=padMainPlot->DrawFrame(0.0,-0.05,55,1.00);
   hFrameMain->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   // hFrameMain->GetYaxis()->SetTitle("% ");
   hFrameMain->GetYaxis()->SetLabelSize(0.0);
   
   hFrameMain->GetXaxis()->CenterTitle();
   hFrameMain->GetYaxis()->CenterTitle();
   
   EventProfile->cd();
   

   TPad *padYaxis=new TPad("Yaxis","",0.0,0.15,0.15,1.0);
   padYaxis->SetRightMargin(0.0);
   padYaxis->SetBottomMargin(0.0);
   padYaxis->SetLeftMargin(0.0);
   padYaxis->SetTopMargin(0.0);
   padYaxis->SetBorderMode(0);
   padYaxis->SetFillColor(0);
   padYaxis->SetBorderSize(2);
   
   
   padYaxis->SetTicks();
   padYaxis->Draw();
   padYaxis->cd();
   

   TGaxis *gYaxis=new TGaxis(1.0,0.00,1.0,0.96,-5,100,510,"-");
   
   gYaxis->SetTitle("Percentage [%]");
   gYaxis->CenterTitle();
   gYaxis->SetTitleOffset(1.1 );
   gYaxis->SetLabelOffset(0.02);
   gYaxis->SetLabelSize(0.3);
   gYaxis->SetTitleSize(0.3);
   gYaxis->SetLabelFont(42);
   gYaxis->SetTitleFont(42);
   
   
   gYaxis->Draw();

   gStyle->SetErrorX(0.);

   for(Int_t i=0;i<iPtBins;i++){

     padMainPlot->cd();
     cout<<fPtMin<<" - "<<fPtBins[i]<<"  ";
     hDenominator=new TH1F("hDenominator","",1,fPtMin,fPtBins[i]);
     hNumerator=new TH1F("hNumerator","",1,fPtMin,fPtBins[i]);
     
     TGraphAsymmErrors *gEfficiency = new TGraphAsymmErrors();
     
     cout<<GetPtAlgorithmicEff(cTemp,fPtMin,fPtBins[i],fEta,fEta+0.5,hDenominator,hNumerator)<<endl;
     
     //cout<<GetPtGlobalTrackingEffWithoutFakeRate(cTemp,fPtMin,fPtBins[i],fEta,fEta+0.5,hDenominator,hNumerator)<<endl;
     //cout<<GetPtGeometricAcc(cTemp,fPtMin,fPtBins[i],fEta,fEta+0.5,hDenominator,hNumerator)<<endl;
     
     hNumerator->Sumw2();
     hDenominator->Sumw2();
     
     gEfficiency->BayesDivide(hNumerator,hDenominator);
     Double_t dYhigh=gEfficiency->GetErrorYhigh(0);
     Double_t dYlow=gEfficiency->GetErrorYlow(0);
     
     gEfficiency->SetPointError(0,0.0,0.0,dYlow,dYhigh);
     gEfficiency->SetMarkerStyle(20);
     gEfficiency->SetMarkerColor(1);
     gEfficiency->SetMarkerSize(1.0);
     
     gEfficiency->Draw("P");  
     delete hDenominator;
     delete hNumerator;

     hDenominator=new TH1F("hDenominator","",1,fPtMin,fPtBins[i]);
     hNumerator=new TH1F("hNumerator","",1,fPtMin,fPtBins[i]);
     
     TGraphAsymmErrors *gEfficiency = new TGraphAsymmErrors();
     
     
     
     cout<<FakeRate(cTemp,fPtMin,fPtBins[i],fEta,fEta+0.5,hDenominator,hNumerator)<<endl;
     
     hNumerator->Sumw2();
     hDenominator->Sumw2();
     
     gEfficiency->BayesDivide(hNumerator,hDenominator);
     dYhigh=gEfficiency->GetErrorYhigh(0);
     dYlow=gEfficiency->GetErrorYlow(0);
     
     gEfficiency->SetPointError(0,0.0,0.0,dYlow,dYhigh);
     
     
     gEfficiency->SetMarkerStyle(24);
     gEfficiency->SetMarkerColor(1);
     gEfficiency->SetMarkerSize(1.);
     
     gEfficiency->Draw("P");
          

     delete hDenominator;
     delete hNumerator;
       
     fPtMin=fPtBins[i]+0.5;   
   }


   
   
   TLegend *legend1=new TLegend();
   legend1->SetX1NDC(0.20);
   legend1->SetX2NDC(0.45);
   legend1->SetY1NDC(0.95);
   legend1->SetY2NDC(0.80);
   legend1->SetFillColor(0);
   legend1->SetFillStyle(0);
   legend1->SetBorderSize(0);
   legend1->SetTextFont(62);
   
  
   TLegendEntry *entry=legend1->AddEntry("Temp","Fake rate","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(2.0);
   entry->SetTextFont(62);

   TLegendEntry *entry1=legend1->AddEntry("Temp","Algorithmic efficiency","P");
   entry1->SetLineColor(1);
   entry1->SetLineStyle(1);
   entry1->SetLineWidth(2);
   entry1->SetMarkerStyle(20);
   entry1->SetMarkerSize(2.0);
   entry1->SetTextFont(62);

   //   legend1->Draw();


   Char_t cTemp1[100];
   
   sprintf(cTemp1,"%1.1f < |#\eta| < %1.1f",fEta,fEta+0.5);
   
   TLatex *lTemp=new TLatex();
   lTemp->SetNDC(kTRUE);
   lTemp->SetTextSize(0.06);
   lTemp->DrawLatex(0.20,0.5,cTemp1);


   sprintf(cTemp1,"CMS Preliminary");
   //   lTemp->DrawLatex(0.20,0.30,cTemp1);
   cout<<endl;
   /*   sprintf(cTemp1,"Geometric Acce");
   lTemp->DrawLatex(0.40,0.45,cTemp1);
   cout<<endl;
   */
   Float_t EtaMin=fEta;
   Float_t EtaMax=fEta+0.5;

  sprintf(cTemp1,"EffFR_b%2.1f_Prob%0.3f_DCA%1.2f_Hits%0.0f_Match%0.2f_Eta%0.2fto%0.2f_Pt.gif",fb,fProbCutOff,fDCACutOff,fHitsCutOff,fHitMatchCutOff,EtaMin,EtaMax);
  // sprintf(cTemp1,"GlobalFR_b%2.1f_Prob%0.3f_DCA%1.2f_Hits%0.0f_Match%0.2f_Eta%0.2fto%0.2f_Pt.gif",fb,fProbCutOff,fDCACutOff,fHitsCutOff,fHitMatchCutOff,EtaMin,EtaMax);
   //sprintf(cTemp1,"AlgoEff_b%2.1f_Eta%0.2fto%0.2f_Pt.gif",fb,EtaMin,EtaMax);
   //sprintf(cTemp1,"LowFakeRate%d.gif",iIndex);
   EventProfile->Print(cTemp1,"gif");
   
    // sprintf(cTemp1,"EffFR_b%2.1f_Prob%0.3f_DCA%1.2f_Hits%0.0f_Match%0.2f_Eta%0.2fto%0.2f_Pt.eps",fb,fProbCutOff,fDCACutOff,fHitsCutOff,fHitMatchCutOff,EtaMin,EtaMax);
    // sprintf(cTemp1,"GeomAcc_b%2.1f_Eta%0.2fto%0.2f_Pt.eps",fb,EtaMin,EtaMax);
   //sprintf(cTemp1,"LowFakeRate%d.eps",iIndex);
   //  EventProfile->Print(cTemp1,"eps");
   
   //sprintf(cTemp1,"EffFR_b%2.1f_Prob%0.3f_DCA%1.2f_Hits%0.0f_Match%0.2f_Eta%0.2fto%0.2f_Pt.C",fb,fProbCutOff,fDCACutOff,fHitsCutOff,fHitMatchCutOff,EtaMin,EtaMax);
   //sprintf(cTemp1,"GeomAcc_b%2.1f_Eta%0.2fto%0.2f_Pt.C",fb,EtaMin,EtaMax);
   //sprintf(cTemp1,"LowFakeRate%d.C",iIndex);
   // EventProfile->Print(cTemp1);
   iIndex++;

  }
}




//Tracking Efficiecny

Float_t GetPtTrackingEff(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin,Float_t fEtaMax){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   fSummary->Close();
   if(dAllTracks!=0)
  	return dAllTracksEff/dAllTracks;
   else
	 return 0;
  
}



Float_t GetPtTrackingEffErr(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   fSummary->Close();
   if(dAllTracks==0||dAllTracksEff==0){
  	return 0;
   }
	
  if(dAllTracksEff!=0){
  	Float_t fAllErr=TMath::Sqrt((1./dAllTracksEff)+(1.0/dAllTracks));	
	return fAllErr*dAllTracksEff/dAllTracks;
  }


}


//Geometric Acceptance

Float_t GetPtGeometricAcc(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5,TH1F *hDenominator,TH1F *hNumerator){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   hDenominator->Fill((fPTMin+fPTMax)/2.0,dAllTracks);
   hNumerator->Fill((fPTMin+fPTMax)/2.0,dAllTracksEff);

   fSummary->Close();
   if(dAllTracks!=0)
     return dAllTracksEff/dAllTracks;
   else
     return 0;
   
}

Float_t GetPtGeometricAccErr(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   //   hDenominator->Fill((fPtMin+fPtMax)/2.0,dAllTracks);
   //hNumerator->Fill((fPtMin+fPtMax)/2.0,dAllTracksEff);


   fSummary->Close();
   if(dAllTracks==0||dAllTracksEff==0){
  	return 0;
   }
	
  if(dAllTracksEff!=0){
  	Float_t fAllErr=TMath::Sqrt((1./dAllTracksEff)+(1.0/dAllTracks));	
	return fAllErr*dAllTracksEff/dAllTracks;
  }


}



//Tracking Efficiecny

Float_t GetPtAlgorithmicEff(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5,TH1F *hDenominator,TH1F *hNumerator){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

   cout<<iPTMin<<" : "<<iPTMax<<" : "<<iEtaMin<<" : "<<iEtaMax<<endl;

	
   Float_t dAllTracks=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   cout<<iPTMin<<" : "<<iPTMax<<" : "<<iEtaMin<<" : "<<iEtaMax<<endl;

   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   hDenominator->Fill((fPTMin+fPTMax)/2.0,dAllTracks);
   hNumerator->Fill((fPTMin+fPTMax)/2.0,dAllTracksEff);

   fSummary->Close();
   if(dAllTracks!=0)
  	return dAllTracksEff/dAllTracks;
   else
	 return 0;
  
}



Float_t GetPtAlgorithmicEffErr(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPtAcc->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hSimEtaPtEff->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   fSummary->Close();
   if(dAllTracks==0||dAllTracksEff==0){
  	return 0;
   }
	
  if(dAllTracksEff!=0){
  	Float_t fAllErr=TMath::Sqrt((1./dAllTracksEff)+(1.0/dAllTracks));	
	return fAllErr*dAllTracksEff/dAllTracks;
  }

}




//Global Tracking Efficiecny

Float_t GetPtGlobalTrackingEff(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5,TH1F *hDenominator,TH1F *hNumerator){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksFR=hRecEtaPtFR->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksFR+=hRecEtaPtFR->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   dAllTracksEff=dAllTracksEff-dAllTracksFR;
   if(dAllTracksEff<0)
      dAllTracksEff=0;
   
   hDenominator->Fill((fPTMin+fPTMax)/2.0,dAllTracks);
   hNumerator->Fill((fPTMin+fPTMax)/2.0,dAllTracksEff);
   
   fSummary->Close();
   if(dAllTracks!=0)
  	return dAllTracksEff/dAllTracks;
   else
	 return 0;
  
}



Float_t GetPtGlobalTrackingEffErr(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   fSummary->Close();
   if(dAllTracks==0||dAllTracksEff==0){
  	return 0;
   }
	
  if(dAllTracksEff!=0){
  	Float_t fAllErr=TMath::Sqrt((1./dAllTracksEff)+(1.0/dAllTracks));	
	return fAllErr*dAllTracksEff/dAllTracks;
  }
  
}


//Global Tracking Efficiecny with out fake rate

Float_t GetPtGlobalTrackingEffWithoutFakeRate(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5,TH1F *hDenominator,TH1F *hNumerator){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   hDenominator->Fill((fPTMin+fPTMax)/2.0,dAllTracks);
   hNumerator->Fill((fPTMin+fPTMax)/2.0,dAllTracksEff);


   fSummary->Close();
   if(dAllTracks!=0)
  	return dAllTracksEff/dAllTracks;
   else
	 return 0;
  
}



Float_t GetPtGlobalTrackingEffErrWithoutFakeRate(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hSimEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   fSummary->Close();
    if(dAllTracks==0||dAllTracksEff==0){
  	return 0;
   }
	
  if(dAllTracksEff!=0){
  	Float_t fAllErr=TMath::Sqrt((1./dAllTracksEff)+(1.0/dAllTracks));	
	return fAllErr*dAllTracksEff/dAllTracks;
  }
 
}



Float_t FakeRate(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5,TH1F *hDenominator,TH1F *hNumerator){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hRecEtaPtFR->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hRecEtaPtFR->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

    hDenominator->Fill((fPTMin+fPTMax)/2.0,dAllTracks);
   hNumerator->Fill((fPTMin+fPTMax)/2.0,dAllTracksEff);

   fSummary->Close();
   if(dAllTracks!=0)
  	return dAllTracksEff/dAllTracks;
   else
	 return 0;
  
}



Float_t FakeRateErr(Char_t *cFileName,Float_t fPTMin,Float_t fPTMax,Float_t fEtaMin=0.0,Float_t fEtaMax=0.5){

   if(fPTMin>fPTMax){
      cout<<"Incorrect pt bin size"<<endl;
	return -1;
   }  
   if(fEtaMin>fEtaMax){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }
   if(fEtaMin<0||fEtaMax<0){
      cout<<"Incorrect eta bin size"<<endl;
	return -1;
   }

	
   TFile *fSummary=new TFile(cFileName);
   Int_t iPTMin=hSimEtaPt->GetYaxis()->FindFixBin(fPTMin);
   Int_t iPTMax=hSimEtaPt->GetYaxis()->FindFixBin(fPTMax);

   Int_t iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   Int_t iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);

	
   Float_t dAllTracks=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   Float_t dAllTracksEff=hRecEtaPtFR->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   fEtaMin=-1.0*fEtaMin;
   fEtaMax=-1.0*fEtaMax;
   iEtaMin=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMax);
   iEtaMax=hSimEtaPt->GetXaxis()->FindFixBin(fEtaMin);
   if(fEtaMin==0)
      iEtaMax--;
   dAllTracks+=hRecEtaPt->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);
   dAllTracksEff+=hRecEtaPtFR->Integral(iEtaMin,iEtaMax,iPTMin,iPTMax);

   fSummary->Close();
    if(dAllTracks==0||dAllTracksEff==0){
  	return 0;
   }
	
  if(dAllTracksEff!=0){
  	Float_t fAllErr=TMath::Sqrt((1./dAllTracksEff)+(1.0/dAllTracks));	
	return fAllErr*dAllTracksEff/dAllTracks;
  }
 
}


