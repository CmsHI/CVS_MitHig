#ifndef PIXELTRACKLETCORR_H
#define PIXELTRACKLETCORR_H

#include <iostream>
#include <string>
#include "Rtypes.h"
#include "TH3D.h"
#include "TFile.h"
#include "TNamed.h"


class TrackletCorrections : public TNamed
{
   public :
     TrackletCorrections();
   TrackletCorrections(TFile* f);
   TrackletCorrections(int hitbins, int etabins, int zbins);

      ~TrackletCorrections(){};

   

      double alpha(double hits, double eta, double z, bool find = true);
      double beta(double hits, double eta, double z, bool find = true);

      double getDeltaR() {return deltaRCut_;}

      //      void setAlpha(double hits, double eta, double z, double value, bool find = true);
      void setBeta(double hits, double eta, double z, double value, bool find = true);

      void setDeltaR(double value) {deltaRCut_ = value;}


      void save(const char* output);

   

   private:

      TH3D* alphas_;
      TH3D* betas_;

      std::string dataTag_;
      std::string mcTag_;

      int hitBins_;
      int etaBins_;

      int zBins_;

      double deltaRCut_;
      double midEtaCut;
      double etaMax_;
      double zMax_;
      double hitMin_;
      double hitMax_;


      // Optional data

      std::string comments_;
      double comEnergy_;

      ClassDef(TrackletCorrections,1)


};



TrackletCorrections::TrackletCorrections() : alphas_(0), betas_(0) {
  std::cout<<"Please make sure the alpha and beta histograms are set correctly. TrackletCorrections::TrackletCorrections(TFile* f) is suggested."<<std::endl;
}

TrackletCorrections::TrackletCorrections(TFile* f) {
  // Get corrections from a file
  
  TrackletCorrections* corr = dynamic_cast<TrackletCorrections*>(f->Get("corrections"));
  alphas_ = dynamic_cast<TH3D*>(f->Get("alpha"));
  betas_ = dynamic_cast<TH3D*>(f->Get("beta"));
  
  *this = *corr;

}

TrackletCorrections::TrackletCorrections(int hitbins, int etabins, int zbins){
   // Create new correction object

  this->SetName("corrections");
  etaMax_ = 2;
  zMax_ = 10;
  hitMin_ = 0;
  hitMax_ = 50;

  alphas_ = new TH3D("alpha","",hitbins,hitMin_,hitMax_,etabins,-etaMax_,etaMax_,zbins,-zMax_,zMax_);
  betas_ = new TH3D("beta","",hitbins,hitMin_,hitMax_,etabins,-etaMax_,etaMax_,zbins,-zMax_,zMax_);

}

double TrackletCorrections::alpha(double hits, double eta, double z, bool find) 
{


  std::cout<<"function working"<<std::endl;

   double alpha = -1;
   if(find){
   int bin = alphas_->FindBin(hits,eta,z);
   alpha = alphas_->GetBinContent(bin);
   }else{
      alpha = alphas_->GetBinContent(hits,eta,z);
   }

   return alpha; 
}

double TrackletCorrections::beta(double hits, double eta, double z, bool find)
{

   double beta = -1;
   if(find){
      int bin = betas_->FindBin(hits,eta,z);
      beta = betas_->GetBinContent(bin);
   }else{
      beta = betas_->GetBinContent(hits,eta,z);
   }

   return beta;

}

void TrackletCorrections::setBeta(double hits, double eta, double z, double value, bool find)
{

  double beta = -1;
  if(find){
    int bin = betas_->FindBin(hits,eta,z);
    beta = betas_->GetBinContent(bin);
  }else{
    beta = betas_->GetBinContent(hits,eta,z);
  }


}


void TrackletCorrections::save(const char* output){

  TFile* out = new TFile(output,"recreate");
  out->cd();
  alphas_->Write();
  betas_->Write();
  this->Write();
  out->Close();

}


ClassImp(TrackletCorrections)


#endif

