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

#define MAXPARTICLES 50000
#define MAXJETS 5000
#define MAXHITS 5000
#define MAXVTX 100
#define ETABINS 3 

using namespace std;

struct HydjetEvent{

   int event;
   float b;
   float npart;
   float ncoll;
   float nhard;

   int n[ETABINS];
   float ptav[ETABINS];

   int np;
   float par_pt[MAXPARTICLES];
   float par_eta[MAXPARTICLES];
   float par_phi[MAXPARTICLES];
   int pdg[MAXPARTICLES];
   int chg[MAXPARTICLES];

   int algos;
   int njet;

   float et[MAXJETS];
   float eta[MAXJETS];
   float phi[MAXJETS];
   float area[MAXJETS];

   float vx;
   float vy;
   float vz;
   float vr;

};

double analyze_with_cut(TTree* tsub, TTree* tall, double jetEtCut = 20);

void analyze_genjets(char * infile = "genjets.root", char * outfile = "output.root"){

  TFile* inf = new TFile(infile);
  TTree* tsub = dynamic_cast<TTree*>(inf->Get("subevent/hi"));
  TTree* tall = dynamic_cast<TTree*>(inf->Get("recoevent/hi"));

  TFile* outf = new TFile(outfile,"recreate");

  double x[6];
  double y[6];

  for(int i=0; i < 3; ++i){
    x[i] = 5*i;
    cout<<"Et Cut  : "<<x[i]<<endl;

    y[i] = analyze_with_cut(tsub,tall,5*i);
    cout<<"Overlap : "<<y[i]<<endl;

  }  

  TGraph* g1 = new TGraph(6,x,y);

  TCanvas* c1 = new TCanvas();
  g1->Draw("a*");

  g1->Write();
  c1->Write();

  TCanvas* c5 = new TCanvas();
  tsub->Draw("njet");
  c5->Print("nsubjets.gif");

  TCanvas* c6 = new TCanvas();
  tall->Draw("njet");
  c6->Print("nglobjets.gif");

  TCanvas* c7 = new TCanvas();
  tsub->Draw("et");
  c7->Print("subjet_et.gif");

  TCanvas* c8 = new TCanvas();
  tall->Draw("et");
  c8->Print("globjet_et.gif");




}

double analyze_with_cut(TTree* tsub, TTree* tall, double jetEtCut){
  cout<<"Begin"<<endl;

   TH1F* h1 = new TH1F(Form("h1_et%02d",(int)jetEtCut),"Self Correlation;#Delta R;jets",200,0,6);
   TH1F* h2 = new TH1F(Form("h2_et%02d",(int)jetEtCut),"Relation between Globally reconstructed and Sub-Event based GenJets;#Delta R;jets",200,0,6);
   TH2F* het = new TH2F(Form("het_et%02d",(int)jetEtCut),";E_{T}^{sub-event genjet};E_{T}^{global genjet}",50,0,200,50,0,200);

   double cone = 0.5;
   double match = cone/2.;

   int maxEvents = 100000;

   cout<<"A"<<endl;

   HydjetEvent jet;
   HydjetEvent jet2;

   tsub->SetBranchAddress("b",&jet.b);
 
   tsub->SetBranchAddress("npart",&jet.npart);
   tsub->SetBranchAddress("ncoll",&jet.ncoll);
   tsub->SetBranchAddress("nhard",&jet.nhard);

   tsub->SetBranchAddress("n",jet.n);
   tsub->SetBranchAddress("ptav",jet.ptav);
   tsub->SetBranchAddress("np",&jet.np);
   tsub->SetBranchAddress("par_pt",jet.par_pt);
   tsub->SetBranchAddress("par_eta",jet.par_eta);
   tsub->SetBranchAddress("par_phi",jet.par_phi);
   tsub->SetBranchAddress("pdg",jet.pdg);
   tsub->SetBranchAddress("chg",jet.chg);
   tsub->SetBranchAddress("njet",&jet.njet);
   tsub->SetBranchAddress("et",jet.et);
   tsub->SetBranchAddress("eta",jet.eta);
   tsub->SetBranchAddress("phi",jet.phi);
   tsub->SetBranchAddress("area",jet.area);
   tsub->SetBranchAddress("vx",&jet.vx);
   tsub->SetBranchAddress("vy",&jet.vy);
   tsub->SetBranchAddress("vz",&jet.vz);
   tsub->SetBranchAddress("vr",&jet.vr);

   tall->SetBranchAddress("njet",&jet2.njet);
   tall->SetBranchAddress("et",jet2.et);
   tall->SetBranchAddress("eta",jet2.eta);
   tall->SetBranchAddress("phi",jet2.phi);
   tall->SetBranchAddress("area",jet2.area);


   cout<<"B"<<endl;


   // Event Loop
   for(int i = 0; i< tsub->GetEntries() && i < maxEvents; ++i){    
     tsub->GetEntry(i);
     tall->GetEntry(i);
     if (i % 1000 == 0) cout <<"Event "<<i<<endl;    
     // Selection on Events

     // Loop over GenJets

     for(int j = 0; j < jet.njet; ++j){
       if(jet.et[j] < jetEtCut) continue;

       for(int j1 = 0; j1< jet.njet; ++j1){
         if(jet.et[j1] < jetEtCut) continue;
         double dR = deltaR(jet.eta[j],jet.phi[j],jet.eta[j1],jet.phi[j1]);

	 if(dR != 0) 
	   h1->Fill(dR);
         if(dR < cone){

         }
       }

       for(int j2 = 0; j2< jet2.njet; ++j2){
	 if(jet2.et[j2] < jetEtCut) continue;
	 double dR = deltaR(jet.eta[j],jet.phi[j],jet2.eta[j2],jet2.phi[j2]);

         h2->Fill(dR);
	 if(dR < match){
	   het->Fill(jet.et[j],jet2.et[j2]);
	 }
       }
      
     }





   }


   cout<<"End."<<endl;

   TCanvas* c1 = new TCanvas();
   h1->Draw();

   TCanvas* c2 = new TCanvas();
   h2->Draw();

   TCanvas* c3 = new TCanvas();
   het->Draw("colz");


   c1->Write();
   c2->Write();
   c3->Write();


   h1->Write();
   h2->Write();
   het->Write();

   double overlap = h1->Integral(0,cone/(h1->GetBinWidth(1)));
   overlap /= h1->Integral();

   return overlap;
}
