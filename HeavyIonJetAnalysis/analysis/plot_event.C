#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "MitHigFunctions.h"
using namespace std;


void plot_event(char * infile = "genjets_test1.root", char * outfile = "event.root"){

  TFile* inf = new TFile(infile);
  TTree* tsub = dynamic_cast<TTree*>(inf->Get("subevent/hi"));
  TTree* tall = dynamic_cast<TTree*>(inf->Get("allevent/hi"));

  TFile* outf = new TFile(outfile,"recreate");

  TCanvas* c1 = new TCanvas();
  tsub->Draw("eta:phi","et> 5","colz");
  c1->Print("subeventjets_et5.gif");
  
  TCanvas* c2 = new TCanvas();
  tsub->Draw("eta:phi","et> 15","colz");
  c2->Print("subeventjets_et15.gif");
  
  TCanvas* c3 = new TCanvas();
  tall->Draw("eta:phi","et> 5","colz");
  c3->Print("globaleventjets_et5.gif");
  
  TCanvas* c4 = new TCanvas();
  tall->Draw("eta:phi","et> 35","colz");
  c4->Print("globaleventjets_et35.gif");


}
