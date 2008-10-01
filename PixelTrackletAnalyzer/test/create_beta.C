#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include "TrackletData.h"
#include "../../PixelTracklet/interface/TrackletCorrections.h"

using namespace std;

void formatHist(TH1* h, int col = 1, double norm = 1);
void saveCanvas(TCanvas* c, int date = 20080829);

void create_beta(const char* infile = "/server/03a/yetkin/data/pythia900GeV_d20080920.root", double etaMax = 1., const char* outfile = "histograms.root"){

  gROOT->Reset();
  gROOT->ProcessLine(".x rootlogon.C");
  gStyle->SetErrorX(0.);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  gStyle->SetPadLeftMargin(0.132);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.02);

  TH1::SetDefaultSumw2();

  double hitmax = 40;
  int hitbins = 1;
  int etabins = 1;

  TrackletCorrections* corr = new TrackletCorrections(hitbins, etabins,1);

  corr->setMidEtaCut(0);
  corr->setDeltaRCut(0.2);
  corr->setCPhi(43);
  corr->setNormDRMin(0.6);
  corr->setNormDRMin(5);

  corr->setHistBins(20);
  corr->setHistMax(5);


  corr->setHitBins(hitbins);
  corr->setHitMax(hitmax);
  corr->setEtaMax(etaMax);
  corr->setZMax(20);


  corr->start();

  int di=hitmax/hitbins;

  int size  =  corr->size();
 
  cout<<"Corrections size : "<<size<<endl;

  TFile *f = new TFile(infile);

  TrackletData data(f,corr);

  TFile *of = new TFile(outfile,"recreate");

  for(int ibin = 0; ibin < size; ++ibin){
    corr->setBeta(ibin,data.getBeta(ibin,true));
  }

  of->Write();
  of->Close();

  corr->save("corrections.root");

}

double TrackletData::getBeta(int bin,bool saveplots)
{

  /// Parameters

  corr->getBin(bin);

  double etaMin = corr->binEtaMin();
  double etaMax = corr->binEtaMax();
  double MinHit = corr->binHitMin();
  double MaxHit = corr->binHitMax();
  double MinZ = corr->binZMin();
  double MaxZ = corr->binZMax();
  
  cout<<"Bin : "<<bin<<endl;
  cout<<"      "<<endl;


  cout<<"Hit Min = "<<MinHit<<endl;
  cout<<"Hit Max = "<<MaxHit<<endl;
  cout<<"      "<<endl;

  cout<<"Eta Min = "<<etaMin<<endl;
  cout<<"Eta Max = "<<etaMax<<endl;
  cout<<"      "<<endl;

  cout<<"Z Min   = "<<MinZ<<endl;
  cout<<"Z Max   = "<<MaxZ<<endl;
  cout<<"      "<<endl;


  cout<<"-------------"<<endl;  
  cout<<"      "<<endl;
  cout<<"      "<<endl;

  cout<<"      "<<endl;

  double deltaCut = corr->getDeltaRCut();
  int etaBins = corr->getEtaBins();
  int nBins = corr->getHistBins(); //20;
  int max = corr->getHistMax(); //20;

  double normMin = corr->getNormDRMin();
  double normMax = corr->getNormDRMax();



  TH1F * h1 = new TH1F(Form("h1_%d",counter),Form("Everything;D;#_{pixel pairs}/event/%.2f",1/(double)nBins),nBins,0,max);
  TH1F * h2 = new TH1F(Form("h2_%d",counter),"Signal",nBins,0,max);
  TH1F * h3 = new TH1F(Form("h3_%d",counter),"Background",nBins,0,max);
  TH1F * h4 = new TH1F(Form("h4_%d",counter),"Normalized Reproduced Background",nBins,0,max);
  TH1F * h5 = new TH1F(Form("h5_%d",counter),"",nBins,0,max);
  TH1F * h6 = new TH1F(Form("h6_%d",counter),"Reproduced Background Subtracted",nBins,0,max);
  TH1F * h7 = new TH1F(Form("h7_%d",counter),"",nBins,0,max);
  TH1F * h8 = new TH1F(Form("h8_%d",counter),"",nBins,0,max);
  TH1F * h9 = new TH1F(Form("h9_%d",counter),"",nBins,0,max);
  TH1F * h10 = new TH1F(Form("h10_%d",counter),"",nBins,0,max);

  TProfile * dNdEtaHadron = new TProfile(Form("dNdEtaHadron_%d",counter),"",32,-2.1,2.1);
  TProfile * dNdEtaLepton = new TProfile(Form("dNdEtaLepton_%d",counter),"",32,-2.1,2.1);
  TProfile * dNdEtaTracklet = new TProfile(Form("dNdEtaTracklet_%d",counter),"",32,-2.1,2.1);

  // TH2D * dNdEtaHadron = new TH2D(Form("dNdEtaHadron_%d",counter),"",32,-2.1,2.1,50,0,5);
  // TH2D * dNdEtaLepton = new TH2D(Form("dNdEtaLepton_%d",counter),"",32,-2.1,2.1,50,0,5);

  // TH2D * correlation = new TH2D(Form("correlation_%d",counter),"; #_{hadrons}; #_{tracklets}",100,0,50,100,0,50);

  counter++;

  int nevents = ntgen->GetEntries();
  int partentries = ntparticle->GetEntries();
  int matchedentries = ntmatched->GetEntries();

  for(int i = 0; i<matchedentries;i++){
    ntmatched->GetEntry(i);
    if(fabs(matchedeta1)>=etaMax) continue;
    if(fabs(matchedeta2)>etaMax) continue;
    if(layer1hits>MaxHit) continue;
    if(layer1hits<MinHit) continue;    

    //    float dR= sqrt(deta*deta+dphi*dphi/43./43.);
    // This is not the standard delta R - It is with the phi normalized with corr->getCPhi();
    double dR = deltaR(deta,dphi);

    h1->Fill(fabs(dR));
    h6->Fill(fabs(dR));
    if(signalcheck==1){
      h2->Fill(fabs(dR));
    }
    if(signalcheck==0){
      h3->Fill(fabs(dR));
    } 
  }

  int invmatchedentries = ntInvMatched->GetEntries();

  for(int i = 0; i<invmatchedentries;i++){
    ntInvMatched->GetEntry(i);
    if(fabs(inveta1)>etaMax) continue;
    //    if(fabs(inveta2)>etaMax) continue;
    if(fabs(inveta1)<0.1) continue;
    
    //    float dR= sqrt(invdeta*invdeta+invdphi*invdphi);
    // This is not the standard delta R - It is with the phi normalized with corr->getCPhi();
    double dR = deltaR(invdeta,invdphi);
    h4->Fill(dR);


    //    h4->Fill(fabs(invdeta));

  }
  
  for(int i = 0; i<partentries;i++){
    ntparticle->GetEntry(i);
    if(fabs(eta1)>etaMax || fabs(eta2)>etaMax) continue;
    if(charge==0) continue;
    h5->Fill(fabs(eta1-eta2));
  }


  if(saveplots){
  formatHist(h1,1,nevents);
  formatHist(h2,2,nevents);
  formatHist(h3,3,nevents);
  formatHist(h4,4,nevents);
  formatHist(h5,5,nevents);
  formatHist(h6,6,nevents);
  }

  //// Normalization of background

  Float_t sc = ((h1->Integral((int)(normMin*nBins/max),(int)(normMax*nBins/max),"width"))/(h4->Integral((int)(normMin*nBins/max),(int)(normMax*nBins/max),"width")));

  cout<<"background normalization: "<<sc<<endl;
  h4->Scale(sc);
  h6->Add(h4,-1);

  //// Determination of correction factor beta
  double beta = (h6->Integral(0,(int)(deltaCut*nBins/max),"width")/(h1->Integral(0,(int)(deltaCut*nBins/max),"width")));
  cout<<"beta: "<<beta<<endl;

  return beta;

}

double TrackletData::getAlpha(int bin, bool saveplots){
  return 0;
}





