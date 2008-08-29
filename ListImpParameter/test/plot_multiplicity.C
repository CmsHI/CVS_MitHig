
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>


void plot_multiplicity(){

  bool normalize = true;

  TFile* file1 = new TFile("pythia_10TeV_eta1.nt.root");
  TFile* file2 = new TFile("pythia_10TeV_eta2.nt.root");
  
  TH1F* h1 = new TH1F("h1","All Charged particles with |#eta|<1;n;events",200,0,200);
  TH1F* h2 = new TH1F("h2","All Charged particles with |#eta|<2;n;events",200,0,200);

  TNtuple * nt1 = dynamic_cast<TNtuple *>(file1->Get("ana/event"));
  TNtuple * nt2 = dynamic_cast<TNtuple *>(file2->Get("ana/event"));

  float ch1,ch2,all1,all2;

  nt1->SetBranchAddress("nc2",&ch1); //charged
  nt1->SetBranchAddress("m1",&all1); //all
  nt2->SetBranchAddress("nc2",&ch2); //charged
  nt2->SetBranchAddress("m1",&all2); //all

  int nevent1 = nt1->GetEntries();
  int nevent2 = nt2->GetEntries();
  
  for(int i1 = 0; i1<nevent1;++i1){
    nt1->GetEntry(i1);
    h1->Fill(ch1);
  }
  
  for(int i2 = 0; i2<nevent2;++i2){
    nt2->GetEntry(i2);
    h2->Fill(ch2);
  }

  if(normalize){
    h1->Scale(1/(double)nevent1);
    h2->Scale(1/(double)nevent2);
  }
  
  h2->SetLineColor(4);
  h1->Draw();  
  h2->Draw("same");

}
