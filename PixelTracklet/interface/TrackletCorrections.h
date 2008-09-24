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
   TrackletCorrections(int hitbins, int etabins = 1, int zbins = 1);

      ~TrackletCorrections(){};
   
      double alpha(double hits, double eta, double z, bool find = true);
      double beta(double hits, double eta, double z, bool find = true);

      double alpha(int bin);
      double beta(int bin);

      double getDeltaRCut() {return deltaRCut_;}

      double getNormDRMin() {return normMax_;}
      double getNormDRMax() {return normMin_;}

      double getMidEtaCut() { return  midEtaCut_;}

      double getCPhi() { return cPhi_;}

      int getHitBins() { return  hitBins_;}
      int getEtaBins() { return  etaBins_;}
      int getZBins() { return  zBins_;}

      int size();  // Total number of divisions in 3D

      void getBin(int value);

      double binEtaMax();
      double binEtaMin();
      double binHitMax();
      double binHitMin();
      double binZMax();
      double binZMin();

      double getHitMax() { return  hitMax_;}
      double getEtaMax() { return  etaMax_;}
      double getZMax() { return  zMax_;}

      std::string getDataTag() { return  dataTag_;}
      std::string getMCTag() { return  mcTag_;}

      void start();

      void setAlpha(double hits, double eta, double z, double value, bool find = true);
      void setBeta(double hits, double eta, double z, double value, bool find = true);

      void setAlpha(int bin, double value);
      void setBeta(int bin, double value);

      void setDeltaRCut(double value) {deltaRCut_ = value;}
      void setMidEtaCut(double value) { midEtaCut_ = value;}

      void setCPhi(double value) {cPhi_ = value;}

      void setHitBins(int value) {hitBins_ = value;}
      void setEtaBins(int value) {etaBins_ = value;}
      void setZBins(int value) {zBins_ = value;}

      void setHitMax(double value) {hitMax_ = value;}
      void setEtaMax(double value) {etaMax_ = value;}
      void setZMax(double value) {zMax_ = value;}

      void setDataTag(std::string value) { dataTag_= value;}
      void setMCTag(std::string value) { mcTag_= value;}

      void save(const char* output);

   

   private:

      TH3D* alphas_;
      TH3D* betas_;

      std::string dataTag_;
      std::string mcTag_;

      int hitBins_;
      int etaBins_;

      int zBins_;

      int bin_;
      int binx_;
      int biny_;
      int binz_;



      // Analysis Parameters

      double deltaRCut_;
      double midEtaCut_;
      double normMax_;
      double normMin_;
      double cPhi_;

      // Analysis Regions

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
  bin_ = -1;
  binx_ = -1;
  biny_ = -1;
  binz_ = -1;
  std::cout<<"Please make sure the alpha and beta histograms are set correctly. TrackletCorrections::TrackletCorrections(TFile* f) is suggested."<<std::endl;
}

TrackletCorrections::TrackletCorrections(TFile* f) {
  // Get corrections from a file
  
  TrackletCorrections* corr = dynamic_cast<TrackletCorrections*>(f->Get("corrections"));
  alphas_ = dynamic_cast<TH3D*>(f->Get("alpha"));
  betas_ = dynamic_cast<TH3D*>(f->Get("beta"));

  bin_ = -1;
  binx_ = -1;
  biny_ = -1;
  binz_ = -1;

  *this = *corr;

}

TrackletCorrections::TrackletCorrections(int hitbins, int etabins, int zbins) : alphas_(0), betas_(0) {
   // Create new correction object

  this->SetName("corrections");

  // Deafults - should be removed
  etaMax_ = 2;
  zMax_ = 10;
  hitMin_ = 0;
  hitMax_ = 50;

  hitBins_ = hitbins;
  etaBins_ = etabins;
  zBins_ = zbins;

  bin_ = -1;
  binx_ = -1;
  biny_ = -1;
  binz_ = -1;

  alphas_ = new TH3D("alpha","",hitbins,hitMin_,hitMax_,etabins,-etaMax_,etaMax_,zbins,-zMax_,zMax_);
  betas_ = new TH3D("beta","",hitbins,hitMin_,hitMax_,etabins,-etaMax_,etaMax_,zbins,-zMax_,zMax_);

}

double TrackletCorrections::alpha(int bin){
  this->getBin(bin);
  return this->alpha(binx_,biny_,binz_,false);
}

double TrackletCorrections::beta(int bin){
  this->getBin(bin);
  return this->beta(binx_,biny_,binz_,false);
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

int TrackletCorrections::size(){
  return (alphas_->GetNbinsX())*(alphas_->GetNbinsY())*(alphas_->GetNbinsZ());
}


void TrackletCorrections::start()
{
  delete alphas_;
  delete betas_;
alphas_ = new TH3D("alpha","",hitBins_,hitMin_,hitMax_,etaBins_,-etaMax_,etaMax_,zBins_,-zMax_,zMax_);
betas_ = new TH3D("beta","",hitBins_,hitMin_,hitMax_,etaBins_,-etaMax_,etaMax_,zBins_,-zMax_,zMax_);
}



void TrackletCorrections::setAlpha(double hits, double eta, double z, double value, bool find)
{
  if(find){
    int bin = alphas_->FindBin(hits,eta,z);
    alphas_->SetBinContent(bin,value);
  }else{
    alphas_->SetBinContent(hits,eta,z,value);
  }
}

void TrackletCorrections::setAlpha(int bin, double value){
  this->getBin(bin);
  this->setAlpha(binx_,biny_,binz_,value,false);
}

void TrackletCorrections::setBeta(int bin, double value){
  this->getBin(bin);
  this->setBeta(binx_,biny_,binz_,value,false);
}

void TrackletCorrections::setBeta(double hits, double eta, double z, double value, bool find)
{
  if(find){
    int bin = betas_->FindBin(hits,eta,z);
    betas_->SetBinContent(bin,value);
  }else{
    betas_->SetBinContent(hits,eta,z,value);
  }
}

void TrackletCorrections::getBin(int value){

  // Modified version of TH1::GetBinXYZ() with overflows skipped

  bin_ = value;

  int nx = hitBins_;
  int ny = etaBins_;

  binx_ = (bin_%nx);
  biny_ = (((bin_-binx_)/nx)%ny);
  binz_ = (((bin_-binx_)/nx -biny_)/ny);

  binx_++;
  biny_++;
  binz_++;

}

double TrackletCorrections::binHitMax(){
  double result = -99999;
  if(bin_ < 0){
    std::cout<<"Bin not set!"<<std::endl;
  }else{
    result = alphas_->GetXaxis()->GetBinUpEdge(binx_);
  }
  return result;
}

double TrackletCorrections::binEtaMax(){
  double result = -99999;
  if(bin_ < 0){
    std::cout<<"Bin not set!"<<std::endl;
  }else{
    result = alphas_->GetYaxis()->GetBinUpEdge(biny_);
  }
  return result;
}

double TrackletCorrections::binZMax(){
  double result = -99999;
  if(bin_ < 0){
    std::cout<<"Bin not set!"<<std::endl;
  }else{
    result = alphas_->GetZaxis()->GetBinUpEdge(binz_);
  }
  return result;
}


double TrackletCorrections::binHitMin(){
  double result = -99999;
  if(bin_ < 0){
    std::cout<<"Bin not set!"<<std::endl;
  }else{
    result = alphas_->GetXaxis()->GetBinLowEdge(binx_);
  }
  return result;
}

double TrackletCorrections::binEtaMin(){
  double result = -99999;
  if(bin_ < 0){
    std::cout<<"Bin not set!"<<std::endl;
  }else{
    result = alphas_->GetYaxis()->GetBinLowEdge(biny_);
  }
  return result;
}

double TrackletCorrections::binZMin(){
  double result = -99999;
  if(bin_ < 0){
    std::cout<<"Bin not set!"<<std::endl;
  }else{
    result = alphas_->GetZaxis()->GetBinLowEdge(binz_);
  }
  return result;
}

void TrackletCorrections::save(const char* output){

  TFile* out = new TFile(output,"recreate");
  out->cd();
  alphas_->Write();
  betas_->Write();
  this->Write();

  std::cout<<"Tracklet Corrections saved with values: "<<std::endl;
  std::cout<<"Delta R Cut                 : "<<deltaRCut_<<std::endl;
  std::cout<<"Number of Layer 1 Hit Bins  : "<<hitBins_<<std::endl;
  std::cout<<"Maximum Hits in Layer 1     : "<<hitMax_<<std::endl;
  std::cout<<"Number of Eta Bins          : "<<etaBins_<<std::endl;
  std::cout<<"Maximum Eta                 : "<<etaMax_<<std::endl;
  std::cout<<"   : "<<std::endl;

  out->Close();

}


ClassImp(TrackletCorrections)


#endif

