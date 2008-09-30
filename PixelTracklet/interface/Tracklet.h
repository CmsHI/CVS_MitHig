#ifndef PIXELTRACKLET_H
#define PI 3.14159265358979
#define PIXELTRACKLET_H

class Tracklet
{
   public:
      Tracklet(double eta1,double eta2,double phi1,double phi2);
      ~Tracklet(){};

      double eta1() { return eta1_; }
      double eta2() { return eta2_; }
      double phi1() { return phi1_; }
      double phi2() { return phi2_; }

      double deta() { return eta1_-eta2_; }
      double dphi();
      double dR();
      double dR2();
      
 
      void   setIt1(int i) { it1_ = i; }
      void   setIt2(int i) { it2_ = i; }
      int    getIt1() { return it1_; }
      int    getIt2() { return it2_; }

      void   setId(int i) { id_ = i; }
      void   setSId(int i) { sid_ = i; }
      void   setType(int i) { ptype_ = i; }
      void   setId1(int i) { id1_ = i; }
      void   setId2(int i) { id2_ = i; }
      int    getId() { return id_; }
      int    getSId() { return sid_; }
      int    getType() { return ptype_; }
      int    getId1() { return id1_; }
      int    getId2() { return id2_; }


   private:

      double eta1_;
      double eta2_;

      double phi1_;
      double phi2_;

      int it1_;   // first iterator
      int it2_;   // second iterator
      
      int sid_;   // signal exist?
      int ptype_;  // process type
      int id_;
      int id1_;
      int id2_;
      
};

#endif

