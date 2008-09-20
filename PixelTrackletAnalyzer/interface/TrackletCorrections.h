#ifndef PIXELTRACKLETCORR_H
#define PI 3.14159265358979
#define PIXELTRACKLETCORR_H


#include "TH3D.h"
#include "TFile.h"
#include "TNamed.h"

class TrackletCorrections : public TNamed
{
   public :
      TrackletCorrections(TFile* f);
      TrackletCorrections(int hitbins, int etabins, int zbins);

      ~TrackletCorrections(){};

      double alpha(double hits, double eta, double z, bool val = true);
      double beta(double hits, double eta, double z, bool val = true);

      void save();

   private:

      TH3D* alphas_;
      TH3D* betas_;

      string dataTag_;
      string mcTag_;

      int hitBins_;
      int etaBins_;

      int zBins_;

      double deltaRCut_;
      double midEtaCut;
      double etaMin_;
      double etaMax_;
      double zMax_;
      double hitMin_;
      double hitMax_;


      // Optional data

      string comments_;
      double comEnergy_;

};

TrackletCorrections::TrackletCorrections(TFile* f) {
   // Get corrections from a file

}

TrackletCorrections::TrackletCorrections(int hitbins, int etabins, int zbins){
   // Create new correction object

}

double TrackletCorrections::alpha(double hits, double eta, double z, bool val) 
{
   double alpha = -1;
   if(val){
   int bin = alphas_->FindBin(hits,eta,z);
   alpha = alphas_->GetBinContent(bin);
   }else{
      alpha = alphas_->GetBinContent(hits,eta,z);
   }

   return alpha; 
}


double TrackletCorrections::beta(double hits, double eta, double z, bool val)
{
   double beta = -1;
   if(val){
      int bin = betas_->FindBin(hits,eta,z);
      beta = betas_->GetBinContent(bin);
   }else{
      beta = betas_->GetBinContent(hits,eta,z);
   }

   return beta;

}

void TrackletCorrections::save(){

}

#endif

