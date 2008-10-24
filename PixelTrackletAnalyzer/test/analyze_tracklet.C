#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "RecoHit.h"
#include "tracklet.h"
#include "Math/Vector3D.h"

using namespace std;

class SelectionCriteria {
   public:
       
   double drCut   ;       // to remove double hit
   double dPhiCut ;       // to remove double hit
   double vzCut   ;       // vertex cut

   bool verbose_ ;
   bool useDeltaPhi_;
   bool checkSecondLayer_;
};

class Parameters {
   public:

   float eta1[1000],phi1[1000],r1[1000],eta2[1000],phi2[1000],r2[1000],vz[100];
   int nhits1,nhits2,mult,nv;
};

bool compareEta(RecoHit a,RecoHit b) { return a.eta<b.eta;}
bool comparePhi(RecoHit a,RecoHit b) { return a.phi<b.phi;}
bool compareDeltaR(Tracklet a,Tracklet b) { return fabs(a.dR2())<fabs(b.dR2());}
bool compareDeltaEta(Tracklet a,Tracklet b) {return fabs(a.deta())<fabs(b.deta());}

vector<RecoHit> removeDoubleHits(Parameters par, SelectionCriteria cuts,Int_t layer);
vector<Tracklet> cleanTracklets(vector<Tracklet> input, int matchNumber, SelectionCriteria cuts);
vector<Tracklet> recoProtoTracklets(vector<RecoHit> firstLayerHits, vector<RecoHit> secondLayerHits);

void analyze_tracklet(char * infile){
  TFile* inf = new  TFile(infile);
  TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));

  int zbins = 1;
  int hitbins = 100;
  int nbins = zbins*hitbins;
  
  // Selection on Hits and events
  SelectionCriteria cuts;
  cuts.drCut   = 0.4;      // to remove double hit
  cuts.dPhiCut = 0.04;     // to remove double hit
  cuts.vzCut   = 10;       // vertex cut

  cuts.verbose_ = false;
  cuts.useDeltaPhi_ = false;
  cuts.checkSecondLayer_ = true;
  
  // Output PDF
  TFile* outf = new TFile("output.root","recreate");
  TNtuple *ntmatched = new TNtuple("ntmatched","","eta1:matchedeta:phi1:matchedphi:deta:dphi:signalCheck:tid:r1id:r2id:evtid:nhit1:sid:ptype:vz");
  TNtuple *ntmult = new TNtuple("ntmult","","mult:nhit1:nhit2");

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
    layer1HitPhi[i] = new TH1D(Form("dNdPhiHits1_%02d",i),"dNdPhi Hits Layer 1",500,-3.2,3.2);
    layer2HitPhi[i] = new TH1D(Form("dNdPhiHits2_%02d",i),"dNdPhi Hits Layer 2",500,-3.2,3.2);
  }
  
  TH1D* hm1 = new TH1D("hm1","Number of Hits Layer 1",50,0,50);

  TH3F* nhits = new TH3F("nhits","",100,0,100,100,0,100,100,0,100);
  
  // Parameters for the tree:
  Parameters par;  

  t->SetBranchAddress("eta1",par.eta1);
  t->SetBranchAddress("phi1",par.phi1);
  t->SetBranchAddress("r1",par.r1);
  t->SetBranchAddress("eta2",par.eta2);
  t->SetBranchAddress("phi2",par.phi2);
  t->SetBranchAddress("r2",par.r2);
  t->SetBranchAddress("nhits1",&par.nhits1);
  t->SetBranchAddress("nhits2",&par.nhits2);
  t->SetBranchAddress("mult",&par.mult);
  t->SetBranchAddress("vz",par.vz);
  t->SetBranchAddress("nv",&par.nv);

  cout <<"Number of Events: "<<t->GetEntries()<<endl;

  // Main loop
  for(int i = 0; i< t->GetEntries(); ++i){    
    t->GetEntry(i);
    if (i % 1000 == 0) cout <<"Event "<<i<<endl;    
    // Selection on Events
    if (fabs(par.vz[1])>cuts.vzCut) continue;

    hm1->Fill(par.mult);

    // Process the first layer
    vector<RecoHit> layer1 = removeDoubleHits(par, cuts,1);
    double mult = 0;
    for(int ihit = 0; ihit< (int)layer1.size(); ++ihit) {
      int hitbin1 = (int)layer1.size();
      if (hitbin1 > 99) hitbin1 = 99;
      layer1HitEta[hitbin1]->Fill(layer1[ihit].eta);
      layer1HitPhi[hitbin1]->Fill(layer1[ihit].phi);
      if(fabs(layer1[ihit].eta)<1) mult++;
    }

    // Process the second layer
    vector<RecoHit> layer2 = removeDoubleHits(par, cuts,2);

    for(int ihit = 0; ihit< (int)layer2.size(); ++ihit) {
      int hitbin2 = (int)layer2.size();
      if (hitbin2 > 99) hitbin2 = 99;
      layer2HitEta[hitbin2]->Fill(layer2[ihit].eta);
      layer2HitPhi[hitbin2]->Fill(layer2[ihit].phi);
    }

    // Form Tracklets        
    vector<Tracklet> protoTracklets = recoProtoTracklets(layer1,layer2);
    vector<Tracklet> recoTracklets = cleanTracklets(protoTracklets,0,cuts);

    // Fill Ntuple
    for (int j=0;j<(int)recoTracklets.size();j++)
    {
        float var[100];
	var[0] = recoTracklets[j].eta1();
	var[1] = recoTracklets[j].eta2();
	var[2] = recoTracklets[j].phi1();
	var[3] = recoTracklets[j].phi2();
	var[4] = recoTracklets[j].deta();
	var[5] = recoTracklets[j].dphi();
	var[6] = 0;
	var[7] = recoTracklets[j].getId();
	var[8] = recoTracklets[j].getId1();
	var[9] = recoTracklets[j].getId2();
	var[10] = i;
	var[11] = (int)par.mult;
	var[12] = recoTracklets[j].getSId();
	var[13] = recoTracklets[j].getType();
	var[14] = par.vz[1];
        ntmatched->Fill(var);
    }

    nhits->Fill(mult,layer1.size(),layer2.size());
    ntmult->Fill(mult,layer1.size(),layer2.size());
  }

  outf->Write();
  outf->Close(); 
}

vector<RecoHit> removeDoubleHits(Parameters par, SelectionCriteria cuts,Int_t layer)
{
    vector<RecoHit> hits;
    vector<RecoHit> cleanedHits;

    if (layer == 1) {
       for(int ihit = 0; ihit < par.nhits1; ++ihit){
         RecoHit tmp(par.eta1[ihit],par.phi1[ihit],par.r1[ihit]);
         hits.push_back(tmp);
       }
    } else {
       for(int ihit = 0; ihit < par.nhits2; ++ihit){
         RecoHit tmp(par.eta2[ihit],par.phi2[ihit],par.r2[ihit]);
         hits.push_back(tmp);
       }
    }
    sort (hits.begin(),hits.end(),comparePhi);
    
    for(int ihit = 0; ihit < (int)hits.size(); ++ihit) {
      double dr=0;
      double dphi=10;
      if (ihit !=0) {
         dphi = fabs(hits[ihit-1].phi - hits[ihit].phi);
	 dr   = fabs(hits[ihit-1].r - hits[ihit].r);
      }
      
      if (dr>cuts.drCut && dphi < cuts.dPhiCut) continue;
      
      // recalculate eta and phi
      double x = hits[ihit].r*cos(hits[ihit].phi);
      double y = hits[ihit].r*sin(hits[ihit].phi);
      double z = hits[ihit].r/tan(atan(exp(-hits[ihit].eta))*2);
      ROOT::Math::XYZVector tmpVector(x,y,z-par.vz[1]);
      RecoHit tmpHit(tmpVector.eta(),tmpVector.phi(),tmpVector.rho());
      cleanedHits.push_back(tmpHit);      
    }
    return cleanedHits;
}

vector<Tracklet> recoProtoTracklets(vector<RecoHit> hits1, vector<RecoHit> hits2)
{
   vector<Tracklet> protoTracklets;

   for (int i = 0; i < (int) hits1.size(); i++)
   {
      for (int j = 0; j < (int) hits2.size(); j++)
      {
         Tracklet mytracklet(hits1[i].eta,hits2[j].eta,hits1[i].phi,hits2[j].phi);
         mytracklet.setIt1(i);
	 mytracklet.setIt2(j);
   	 protoTracklets.push_back(mytracklet);
      }
   }
   
   return protoTracklets;
}


vector<Tracklet> cleanTracklets(vector<Tracklet> input, int matchNumber,SelectionCriteria cuts)
{
   vector<Tracklet> output;

   if(cuts.useDeltaPhi_)
      sort( input.begin() , input.end() , compareDeltaR);
   else
      sort( input.begin() , input.end() , compareDeltaEta);

   if (cuts.verbose_) {
      for (unsigned int i = 0; i < input.size(); i++)
      {
         cout <<input[i].deta()<<" "<<input[i].getIt1()<<" "<<input[i].getIt2()<<endl;
      }
   }

   int used1[1000];
   int used2[1000];

   for (int i=0;i<1000;i++) { 
      used1[i]=0;
      used2[i]=0;
   } 

   if (cuts.verbose_) cout<<"Printing Hits"<<endl;
   
   for (unsigned int i = 0; i < input.size(); i++){
      

      if(cuts.useDeltaPhi_)
	 if (cuts.verbose_) cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta R : "<<input[i].dR()<<endl;
      else
	 if (cuts.verbose_) cout<<"Eta 1 : "<<input[i].eta1()<<"  ; Eta 2 : "<<input[i].eta2()<<" ;  Delta Eta : "<<input[i].deta()<<endl; 
      
      int i1=input[i].getIt1();
      int i2=input[i].getIt2();

      if (used1[i1]==0&&used2[i2]==matchNumber) {
	 Tracklet tmp = input[i];
	 output.push_back(tmp);
	 used1[i1]++;
	 if (cuts.checkSecondLayer_) used2[i2]++;
      }
      if (used1[i1]==0&&used2[i2]<matchNumber) {
	 if (cuts.checkSecondLayer_) used2[i2]++;
      }
   }
   if (cuts.verbose_) {
      cout <<"Output:"<<endl;
      for (unsigned int i = 0; i < output.size(); i++)
      {
         cout <<output[i].deta()<<" "<<output[i].getIt1()<<" "<<output[i].getIt2()<<endl;
      }
   }
   
   return output;
}

