class RecoHit {
   public:

   RecoHit(double _eta,double _phi,double _r) 
   { 
      eta = _eta;
      phi = _phi;
      r = _r;
   }; 
   ~RecoHit(){};
   
      double eta;
      double phi;
      double r;
};
