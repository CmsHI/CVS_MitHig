#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"

using namespace std;

void analyze_tree(){
  TFile* inf = new TFile("pythia_mb_900GeV_vtxFlat_d20081015_1.root");
  TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));

  int zbins = 1;
  int hitbins = 100;
  int nbins = zbins*hitbins;
  
  TFile* outf = new TFile("pdf.root","recreate");
  vector<TH1D*> layer1HitEta;
  layer1HitEta.reserve(nbins);
  vector<TH1D*> layer1HitPhi;
  layer1HitPhi.reserve(nbins);

  vector<TH1D*> layer2HitEta;
  layer2HitEta.reserve(nbins);
  vector<TH1D*> layer2HitPhi;
  layer2HitPhi.reserve(nbins);
  
  for(int i = 0; i< nbins; ++i){
    layer1HitEta[i] = new TH1D(Form("dNdEtaHits1_%02d",i),"dNdEta Hits Layer 1",500,-3,3);
    layer2HitEta[i] = new TH1D(Form("dNdEtaHits2_%02d",i),"dNdEta Hits Layer 2",500,-3,3);
    layer1HitPhi[i] = new TH1D(Form("dNdPhiHits1_%02d",i),"dNdPhi Hits Layer 1",500,-3,3);
    layer2HitPhi[i] = new TH1D(Form("dNdPhiHits2_%02d",i),"dNdPhi Hits Layer 2",500,-3,3);
  }
  
  TH1D* h1 = new TH1D("h1","dNdEta Hits Layer 1",500,-3,3);
  TH1D* h2 = new TH1D("h2","dNdEta Hits Layer 2",500,-3,3);
  TH1D* hm1 = new TH1D("hm1","Number of Hits Layer 1",50,0,50);

  TH3F* nhits = new TH3F("nhits","",100,0,100,100,0,100,100,0,100);

  float eta1[1000],phi1[1000],r1[1000],eta2[1000],phi2[1000],r2[1000],vz[100];

  int nhits1,nhits2,mult,nv;

  t->SetBranchAddress("eta1",eta1);
  t->SetBranchAddress("phi1",phi1);
  t->SetBranchAddress("r1",r1);
  t->SetBranchAddress("eta2",eta2);
  t->SetBranchAddress("phi2",phi2);
  t->SetBranchAddress("r2",r2);
  t->SetBranchAddress("nhits1",&nhits1);
  t->SetBranchAddress("nhits2",&nhits2);
  t->SetBranchAddress("mult",&mult);
  t->SetBranchAddress("vz",vz);
  t->SetBranchAddress("nv",&nv);

  cout <<t->GetEntries()<<endl;
  for(int i = 0; i< t->GetEntries(); ++i){
    t->GetEntry(i);
    if (fabs(vz[0])>10) continue;
    nhits->Fill(mult,nhits1,nhits2);
    //std::cout<<"Number of Hits : "<<nhits1<<endl;
    //std::cout<<"Hit Multiplicity : "<<mult<<endl;

    hm1->Fill(mult);
    for(int ihit = 0; ihit < nhits1; ++ihit){
      int hitbin = (int)nhits1;
      if (hitbin>99) hitbin=99;
      layer1HitEta[hitbin]->Fill(eta1[ihit]);
      layer1HitPhi[hitbin]->Fill(phi1[ihit]);
    }
    for(int ihit = 0; ihit < nhits2; ++ihit){
      //std::cout<<"Hit Eta : "<<eta2[ihit]<<endl;
      int hitbin = (int)nhits2;
      if (hitbin>99) hitbin=99;
      layer2HitEta[hitbin]->Fill(eta2[ihit]);
      layer2HitPhi[hitbin]->Fill(phi2[ihit]);
    }
  }
  h1->Draw();
  hm1->Draw();

  t->Draw("nhits1");
  outf->Write();
  outf->Close(); 
}
