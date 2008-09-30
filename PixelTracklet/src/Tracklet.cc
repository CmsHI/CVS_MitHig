#include <math.h>
#include "MitHig/PixelTracklet/interface/Tracklet.h"

Tracklet::Tracklet(double eta1, double eta2, double phi1, double phi2) {

   eta1_ = eta1;
   eta2_ = eta2;

   phi1_ = phi1;
   while (phi1_>2*PI) phi1_-=2*PI;

   phi2_ = phi2;
   while (phi2_>2*PI) phi2_-=2*PI;
}

double Tracklet::dphi() 
{
   double dphi=phi1_-phi2_;

   if (dphi>0){
      while (dphi>2*PI) dphi-=2*PI;
      if (dphi>PI) dphi=2*PI-dphi;
   } else {
      while (dphi<-2*PI) dphi+=2*PI;
      if (dphi<-PI) dphi=-2*PI-dphi;
   }


   return dphi; 
}

double Tracklet::dR() 
{
   return sqrt(dR2());
}

double Tracklet::dR2() 
{
   double dPhi=dphi();
   double dEta=deta();

   return dPhi*dPhi+dEta*dEta;
}





