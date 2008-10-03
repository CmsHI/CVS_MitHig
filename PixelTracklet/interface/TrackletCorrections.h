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

      int getHistBins() { return histBins_;}
      int getHistMax() { return histMax_;}

      int getHitBins() { return  hitBins_;}
      int getEtaBins() { return  etaBins_;}
      int getZBins() { return  zBins_;}

      int size();  // Total number of divisions in 3D

      int findBin(double hits, double eta, double z);
      void getBin(int value);
      void printBin();

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

      void setNormDRMin(double value) {normMin_ = value;}
      void setNormDRMax(double value) {normMax_ = value;}

      void setMidEtaCut(double value) { midEtaCut_ = value;}

      void setCPhi(double value) {cPhi_ = value;}

      void setHistBins(int value) { histBins_ = value;}
      void setHistMax(int value) { histMax_ = value;}

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

      int histBins_;
      int histMax_;

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

      ClassDef(TrackletCorrections,2)

};

#endif

