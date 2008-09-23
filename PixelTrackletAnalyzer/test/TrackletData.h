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
#include "../../PixelTracklet/interface/TrackletCorrections.h"

struct TrackletData
{
  TrackletData(
               TNtuple *nt1=0,
               TNtuple *nt2=0,
               TNtuple *nt3=0,
               TNtuple *nt4=0,
               TNtuple *nt5=0,
               TNtuple *nt6=0
               ) :
    ntparticle(nt1),
    ntmatched(nt2),
    ntInvMatched(nt3),
    ntgen(nt4),
    ntevent(nt5),
    ntvertex(nt6)
    {;}

  TrackletData(TFile* f){

    f->cd("ana;1");
    ntparticle= dynamic_cast<TNtuple *>(f->Get("anasim/ntparticle"));
    ntmatched= dynamic_cast<TNtuple *>(f->Get("ana/ntmatched"));
    ntInvMatched= dynamic_cast<TNtuple *>(f->Get("ana/ntInvMatched"));
    ntgen = dynamic_cast<TNtuple *>(f->Get("ana/ntgen"));
    ntevent = dynamic_cast<TNtuple *>(f->Get("ana/ntevent"));
    ntvertex = dynamic_cast<TNtuple *>(f->Get("ana/ntvertex"));
    ntcorr = new TNtuple("ntcorr","","nhad:ntrt:nstrt:vz:nvtx:vntrk");

    counter = 0;

    ntvertex->SetBranchAddress("z",&vz);
    ntvertex->SetBranchAddress("n",&vntrk);
    ntvertex->SetBranchAddress("nvtx",&nvtx);

    ntmatched->SetBranchAddress("eta1",&matchedeta1);
    ntmatched->SetBranchAddress("phi1",&matchedphi1);
    ntmatched->SetBranchAddress("matchedeta",&matchedeta2);
    ntmatched->SetBranchAddress("matchedphi",&matchedphi2);
    ntmatched->SetBranchAddress("signalCheck",&signalcheck);
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



  }

  ~TrackletData(){;}
  TNtuple * ntparticle;
  TNtuple * ntmatched;
  TNtuple * ntInvMatched;
  TNtuple * ntgen;
  TNtuple * ntevent;
  TNtuple * ntvertex;
  TNtuple * ntcorr;

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

  double getBeta(TrackletCorrections* corr, int bin, bool saveplots = false);
  double getAlpha(TrackletCorrections* corr, int bin, bool saveplots = false);

  int counter;
};



