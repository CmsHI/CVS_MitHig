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

double analyze_with_cut(TTree* tsub, TTree* treco, double jetEtCut = 20);

void analyze_jets(char * infile = "jets_d20081210.root", char * outfile = "output.root"){

  TH1::SetDefaultSumw2();

  TFile* inf = new TFile(infile);
  TTree* tsub = dynamic_cast<TTree*>(inf->Get("subevent/hi"));
  TTree* treco = dynamic_cast<TTree*>(inf->Get("recoevent/hi"));

  TFile* outf = new TFile(outfile,"recreate");

  double x[6];
  double y[6];

  for(int i=0; i < 3; ++i){
    x[i] = 5*i;
    cout<<"Et Cut  : "<<x[i]<<endl;

    y[i] = analyze_with_cut(tsub,treco,5*i);
    cout<<"Overlap : "<<y[i]<<endl;

  }  

  TGraph* g1 = new TGraph(6,x,y);

  TCanvas* c1 = new TCanvas();
  g1->Draw("a*");

  g1->Write();
  c1->Write();

  TCanvas* c5 = new TCanvas("c5","c5",400,400);
  tsub->Draw("njet");
  c5->Print("nsubjets.gif");

  TCanvas* c6 = new TCanvas("c6","c6",400,400);
  treco->Draw("njet");
  c6->Print("nrecojets.gif");

  TCanvas* c7 = new TCanvas("c7","c7",400,400);
  tsub->Draw("et");
  c7->Print("subjet_et.gif");

  TCanvas* c8 = new TCanvas("c8","c8",400,400);
  treco->Draw("et");
  c8->Print("recojet_et.gif");




}

double analyze_with_cut(TTree* tsub, TTree* treco, double jetEtCut){
  cout<<"Begin"<<endl;

  int nsize[3] = {2000,200,50};

   TH1F* h1 = new TH1F(Form("h1_et%02d",(int)jetEtCut),"Self Correlation;#Delta R;jets",200,0,6);
   TH1F* h2 = new TH1F(Form("h2_et%02d",(int)jetEtCut),"Relation between Globally reconstructed and Sub-Event based Jets;#Delta R;jets",200,0,6);

   TH1F* hres = new TH1F(Form("hres_et%02d",(int)jetEtCut),"Jet Energy Resolution;E_{T}^{CaloJet}-E_{T}^{GenJet} [GeV];jets",200,-50,50);

   //   TH2F* het = new TH2F(Form("het_et%02d",(int)jetEtCut),";E_{T}^{GenJet};E_{T}^{CaloJet} [GeV]",50,0,200,50,0,200);
   TH2F* het = new TH2F(Form("het_et%02d",(int)jetEtCut),";E_{T}^{GenJet};E_{T}^{CaloJet} [GeV]",100,0,50,100,0,50);
   TH2F* hnjet = new TH2F(Form("hnjet_et%02d",(int)jetEtCut),";#_{GenJet};#_{CaloJet} [GeV]",100,0,nsize[(int)(jetEtCut/5.)],100,0,50);

   TH1F* heff1 = new TH1F(Form("heff1_et%02d",(int)jetEtCut),"Jet Finding Efficiency;E_{T} [GeV];efficiency",100,0,50);
   TH1F* heff2 = new TH1F(Form("heff2_et%02d",(int)jetEtCut),"Jet Finding Efficiency;E_{T} [GeV];efficiency 2",100,0,50);
   TH1F* heff3 = new TH1F(Form("heff2_et%02d",(int)jetEtCut),";E_{T} [GeV];efficiency",100,0,50);

   TH1F* hfake1 = new TH1F(Form("hfake1_et%02d",(int)jetEtCut),"Fake Jets;E_{T} [GeV];fake",100,0,50);

   TH1F* hcaloet = new TH1F(Form("hcaloet_et%02d",(int)jetEtCut),";E_{T} [GeV];jets",100,0,50);
   TH1F* hgenet = new TH1F(Form("hgenet_et%02d",(int)jetEtCut),";E_{T} [GeV];jets",100,0,50);



   double cone = 0.5;
   double match = cone/4.;

   int maxEvents = 100000;

   cout<<"A"<<endl;

   HydjetEvent jet; //RecoJets
   HydjetEvent jet2; //GenJets

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

   treco->SetBranchAddress("njet",&jet.njet);
   treco->SetBranchAddress("et",jet.et);
   treco->SetBranchAddress("eta",jet.eta);
   treco->SetBranchAddress("phi",jet.phi);
   treco->SetBranchAddress("area",jet.area);

   tsub->SetBranchAddress("vx",&jet.vx);
   tsub->SetBranchAddress("vy",&jet.vy);
   tsub->SetBranchAddress("vz",&jet.vz);
   tsub->SetBranchAddress("vr",&jet.vr);

   tsub->SetBranchAddress("njet",&jet2.njet);
   tsub->SetBranchAddress("et",jet2.et);
   tsub->SetBranchAddress("eta",jet2.eta);
   tsub->SetBranchAddress("phi",jet2.phi);
   tsub->SetBranchAddress("area",jet2.area);


   cout<<"B"<<endl;

   // Event Loop
   for(int i = 0; i< tsub->GetEntries() && i < maxEvents; ++i){    
     tsub->GetEntry(i);
     treco->GetEntry(i);
     int ngenjet = 0;
     int nrecojet = 0;
     bool counted = false;
     int jetmatch[MAXJETS];

     if (i % 1000 == 0) cout <<"Event "<<i<<endl;    
     // Selection on Events

     // Loop over RecoJets
     for(int j = 0; j < jet.njet; ++j){
       if(jet.et[j] < jetEtCut) continue;
       hcaloet->Fill(jet.et[j]);
       nrecojet++;
  
       for(int j1 = 0; j1< jet.njet; ++j1){
         if(jet.et[j1] < jetEtCut) continue;
         double dR = deltaR(jet.eta[j],jet.phi[j],jet.eta[j1],jet.phi[j1]);

	 if(dR != 0) 
	   h1->Fill(dR);
         if(dR < cone){

         }
       }
       int j2match = -99;
       double etmatch = 0;
       for(int j2 = 0; j2< jet2.njet; ++j2){
	 if(jet2.et[j2] < jetEtCut) continue;
	 double dR = deltaR(jet.eta[j],jet.phi[j],jet2.eta[j2],jet2.phi[j2]);
         h2->Fill(dR);
	 if(dR < match && jet2.et[j2] > etmatch){ 
	   j2match = j2;	 
	   etmatch = jet2.et[j2match];
	 }
       }
       jetmatch[j] = j2match;
       if(j2match > -99){
	 het->Fill(jet2.et[j2match],jet.et[j]);
	 hres->Fill(jet.et[j]-jet2.et[j2match]);
       }else{
	 //Jet is a fake
	 hfake1->Fill(jet.et[j]);
       }      
     }
 
     //Loop over GenJets
     for(int j2 = 0; j2< jet2.njet; ++j2){
       if(jet2.et[j2] < jetEtCut) continue;
       hgenet->Fill(jet2.et[j2]);
       ngenjet++;
       int counter  = 0;       
       int j1match = -99;
       double etmatch = 0;
       for(int j1 = 0; j1< jet.njet; ++j1){
	 if(jet.et[j1] < jetEtCut) continue;       
	 double dR = deltaR(jet.eta[j1],jet.phi[j1],jet2.eta[j2],jet2.phi[j2]);
	 if(dR < match && jet2.et[j2] > etmatch){
	   j1match = j1;
	   etmatch = jet.et[j1match];
	 }
	 if(jetmatch[j1] == j2){
	   heff2->Fill(jet2.et[j2]);
	   counter++;
	   if(counter > 1){
	   cout<<"j2 is "<<j2<<endl;
	   cout<<" Jet Energy is "<<jet2.et[j2]<<endl;
	   }
	   break;
	 }
       }     
       if(j1match > -99){
	 //	 cout<<"j2 matched is "<<j2<<endl;
         heff1->Fill(jet2.et[j2]);
       }else{

       }
     }

     hnjet->Fill(ngenjet,nrecojet);
   }

   hfake1->Divide(hfake1,hcaloet,1,1,"B");
   heff1->Divide(heff1,hgenet,1,1,"B");
   heff2->Divide(heff2,hgenet,1,1,"B");

   heff3->Divide(heff2,heff1,1,1,"B");

   
   cout<<"End."<<endl;

   TCanvas* c1 = new TCanvas(Form("c1_et%02d",(int)jetEtCut),Form("c1_et%02d",(int)jetEtCut),400,400);
   h1->Draw();

   TCanvas* c2 = new TCanvas(Form("c2_et%02d",(int)jetEtCut),Form("c2_et%02d",(int)jetEtCut),400,400);
   h2->Draw();

   TCanvas* c3 = new TCanvas(Form("c3_et%02d",(int)jetEtCut),Form("c3_et%02d",(int)jetEtCut),400,400);
   het->Draw("colz");
   c3->Print(Form("EnergyScatter_et%02d.gif",(int)jetEtCut));

   TCanvas* c4 = new TCanvas(Form("c4_et%02d",(int)jetEtCut),Form("c4_et%02d",(int)jetEtCut),400,400);
   hres->Draw("");
   c4->Print(Form("EnergyResolution_et%02d.gif",(int)jetEtCut));

   TCanvas* c5 = new TCanvas(Form("c5_et%02d",(int)jetEtCut),Form("c5_et%02d",(int)jetEtCut),400,400);
   hnjet->Draw("colz");
   c5->Print(Form("NumberOfJets_et%02d.gif",(int)jetEtCut));

   TCanvas* c6 = new TCanvas(Form("c6_et%02d",(int)jetEtCut),Form("c6_et%02d",(int)jetEtCut),400,400);
   heff1->Draw("");
   c6->Print(Form("JetEfficiency_et%02d.gif",(int)jetEtCut));

   TCanvas* c7 = new TCanvas(Form("c7_et%02d",(int)jetEtCut),Form("c7_et%02d",(int)jetEtCut),400,400);
   hfake1->Draw("");
   c7->Print(Form("JetFakes_et%02d.gif",(int)jetEtCut));

   TCanvas* c8 = new TCanvas(Form("c8_et%02d",(int)jetEtCut),Form("c8_et%02d",(int)jetEtCut),400,400);
   heff2->Draw("");
   c8->Print(Form("JetEfficiency2_et%02d.gif",(int)jetEtCut));

   TCanvas* c9 = new TCanvas(Form("c9_et%02d",(int)jetEtCut),Form("c9_et%02d",(int)jetEtCut),400,400);
   heff3->Draw("");
   c9->Print(Form("JetEfficiencyComparison_et%02d.gif",(int)jetEtCut));



   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   c5->Write();


   h1->Write();
   h2->Write();
   het->Write();
   hres->Write();
   hnjet->Write();

   double overlap = h1->Integral(0,cone/(h1->GetBinWidth(1)));
   overlap /= h1->Integral();

   return overlap;
}
