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

void formatHist(TH1* h, int col = 1, double norm = 1);
void saveCanvas(TCanvas* c, int date = 20080829);

double getBeta(const char* infile, double etaMax, int MinHit, int MaxHit);

void analyze_beta(const char* infile = "p0829.root", double etaMax = 1.){

  int di=4;
  TH1F *beta = new TH1F("beta","",40/di,0,40);

  for (int i=0;i<40;i+=di) {
    beta->SetBinContent(i/di+1,getBeta(infile,etaMax,i,i+di));
  }

  TCanvas *c = new TCanvas ("c","",400,400);
  beta->Fit("pol2","m");
  beta->SetXTitle("N_{Hits}");
  beta->SetYTitle("#beta");
  beta->SetAxisRange(0,0.3,"Y");
  beta->Draw("p");
}

double getBeta(const char* infile, double etaMax, int MinHit, int MaxHit)
{
  /// Parameters
  double normRange = 0.6;
  double deltaCut = 0.1;
  int etaBins = 8;
  int nBins = 20;

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

  TFile *f = new TFile(infile);
  f->cd("ana;1");
  TNtuple * ntparticle= dynamic_cast<TNtuple *>(f->Get("anasim/ntparticle"));
  TNtuple * ntmatched= dynamic_cast<TNtuple *>(f->Get("ana/ntmatched"));
  TNtuple * ntInvMatched= dynamic_cast<TNtuple *>(f->Get("ana/ntInvMatched"));
  TNtuple * ntgen = dynamic_cast<TNtuple *>(f->Get("ana/ntgen"));
  TNtuple * ntevent = dynamic_cast<TNtuple *>(f->Get("ana/ntevent"));
  TNtuple * ntvertex = dynamic_cast<TNtuple *>(f->Get("ana/ntvertex"));
 
  TNtuple * ntcorr = new TNtuple("ntcorr","","nhad:ntrt:nstrt:vz:nvtx:vntrk");
  TH1F * h1 = new TH1F("h1",Form("Everything;D;#_{pixel pairs}/event/%.2f",1/(double)nBins),nBins,0,2); 
  TH1F * h2 = new TH1F("h2","Signal",nBins,0,2);
  TH1F * h3 = new TH1F("h3","Background",nBins,0,2);
  TH1F * h4 = new TH1F("h4","Normalized Reproduced Background",nBins,0,2);
  TH1F * h5 = new TH1F("h5","",nBins,0,2);
  TH1F * h6 = new TH1F("h6","Reproduced Background Subtracted",nBins,0,2);
  TH1F * h7 = new TH1F("h7","",nBins,0,2);
  TH1F * h8 = new TH1F("h8","",nBins,0,2);
  TH1F * h9 = new TH1F("h9","",nBins,0,2);
  TH1F * h10 = new TH1F("h10","",nBins,0,2);

  TProfile * dNdEtaHadron = new TProfile("dNdEtaHadron","",32,-2.1,2.1);
  TProfile * dNdEtaLepton = new TProfile("dNdEtaLepton","",32,-2.1,2.1);
  TProfile * dNdEtaTracklet = new TProfile("dNdEtaTracklet","",32,-2.1,2.1);

  // TH2D * dNdEtaHadron = new TH2D("dNdEtaHadron","",32,-2.1,2.1,50,0,5);
  // TH2D * dNdEtaLepton = new TH2D("dNdEtaLepton","",32,-2.1,2.1,50,0,5);

  TH2D * corr = new TH2D("correlation","; #_{hadrons}; #_{tracklets}",100,0,50,100,0,50);

  Float_t matchedeta1;
  Float_t matchedeta2;
  Float_t matchedinveta2;
  Float_t matchedinvphi2;
  Float_t signalcheck;
  Float_t layer1hits;
  Float_t charge;
  Float_t ntpartlayer1hits;
  Float_t matchedphi1;
  Float_t matchedphi2;
 
  Float_t eta1;
  Float_t eta2;
  Float_t deta;
  Float_t dphi;
  Float_t invdeta;
  Float_t invdphi;
  Float_t inveta1;
  Float_t inveta2;
  Float_t invphi1;
  Float_t invphi2;
  Float_t vz;
  Float_t vntrk;
  Float_t nvtx;


  ntvertex->SetBranchAddress("z",&vz);
  ntvertex->SetBranchAddress("n",&vntrk);
  ntvertex->SetBranchAddress("nvtx",&nvtx);

  ntmatched->SetBranchAddress("eta1",&matchedeta1);
  ntmatched->SetBranchAddress("phi1",&matchedphi1);
  ntmatched->SetBranchAddress("matchedeta",&matchedeta2);
  ntmatched->SetBranchAddress("matchedphi",&matchedphi2);
  //  ntmatched->SetBranchAddress("signalCheck",&signalcheck);
  ntmatched->SetBranchAddress("sid",&signalcheck);
  ntmatched->SetBranchAddress("deta",&deta);
  ntmatched->SetBranchAddress("dphi",&dphi);
  ntmatched->SetBranchAddress("nhit1",&layer1hits);

  ntInvMatched->SetBranchAddress("eta1",&inveta1);
  ntInvMatched->SetBranchAddress("phi1",&invphi1);
  ntInvMatched->SetBranchAddress("matchedeta",&inveta2);
  ntInvMatched->SetBranchAddress("matchedphi",&invphi2);
  ntInvMatched->SetBranchAddress("deta",&invdeta);
  ntInvMatched->SetBranchAddress("dphi",&invdphi);

  ntparticle->SetBranchAddress("eta1",&eta1);
  ntparticle->SetBranchAddress("eta2",&eta2);
  ntparticle->SetBranchAddress("charge",&charge);
  ntparticle->SetBranchAddress("layer1hits",&ntpartlayer1hits);


  int nevents = ntgen->GetEntries();
  int partentries = ntparticle->GetEntries();
  int matchedentries = ntmatched->GetEntries();


  for(int i = 0; i<matchedentries;i++){
    ntmatched->GetEntry(i);
    if(fabs(matchedeta1)>=etaMax) continue;
    if(fabs(matchedeta2)>etaMax) continue;
    if(layer1hits>MaxHit) continue;
    if(layer1hits<MinHit) continue;    
    float dR= sqrt(deta*deta+dphi*dphi/43./43.);


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
    

    float dR= sqrt(invdeta*invdeta+invdphi*invdphi);
    h4->Fill(dR);

    //    h4->Fill(fabs(invdeta));

  }
  
  for(int i = 0; i<partentries;i++){
    ntparticle->GetEntry(i);
    if(fabs(eta1)>etaMax || fabs(eta2)>etaMax) continue;
    if(charge==0) continue;
    h5->Fill(fabs(eta1-eta2));
  }

  formatHist(h1,1,nevents);
  formatHist(h2,2,nevents);
  formatHist(h3,3,nevents);
  formatHist(h4,4,nevents);
  formatHist(h5,5,nevents);
  formatHist(h6,6,nevents);

  //// Normalization of background

  Float_t sc = ((h1->Integral((int)(normRange*nBins),nBins,"width"))/(h4->Integral((int)(normRange*nBins),nBins,"width")));

  cout<<"background normalization: "<<sc<<endl;
  h4->Scale(sc);
  h6->Add(h4,-1);

  //// Determination of correction factor beta
  double beta = 1-((h2->Integral(0,(int)(deltaCut*nBins),"width"))/(h6->Integral(0,(int)(deltaCut*nBins),"width")));
  cout<<"beta: "<<beta<<endl;

  /*  
  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->SetLogy();

  TLegend * leg1 = new TLegend(0.25,0.66,0.56,0.84);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(h1,"Everything","l");
//  leg1->AddEntry(h2,"Signal","l");
//  leg1->AddEntry(h3,"Background","l");
  leg1->AddEntry(h4,"Normalized Reproduced Background","l");
  leg1->AddEntry(h6,"Reproduced Background Subtracted","l");

  h1->Draw("e");
//  h2->Draw("e same");
//  h3->Draw("e same");
  h4->Draw("e same");
  h6->Draw("e same");
  leg1->Draw();
  */
  return beta;
}

void formatHist(TH1* h, int col, double norm){

  h->Scale(1/norm);
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

}

void saveCanvas(TCanvas* c, int date){
  c->Write();
  c->Draw();
  c->Print(Form("./figures/%s_d%d.gif",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.eps",c->GetName(),date));
  c->Print(Form("./figures/%s_d%d.C",c->GetName(),date));
}



