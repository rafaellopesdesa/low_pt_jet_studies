#include "TMath.h"

#include "RooDoubleCBFast.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;

 ClassImp(RooDoubleCBFast) 

 RooDoubleCBFast::RooDoubleCBFast(){}

 RooDoubleCBFast::RooDoubleCBFast(const char *name, const char *title, 
                    RooAbsReal& _x,
                    RooAbsReal& _mean,
                    RooAbsReal& _width,
                    RooAbsReal& _alpha1,
                    RooAbsReal& _n1,
                    RooAbsReal& _alpha2,
                    RooAbsReal& _n2
                    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   alpha1("alpha1","alpha1",this,_alpha1),
   n1("n1","n1",this,_n1),
   alpha2("alpha2","alpha2",this,_alpha2),
   n2("n2","n2",this,_n2)
 { 
 } 


 RooDoubleCBFast::RooDoubleCBFast(const RooDoubleCBFast& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   alpha1("alpha1",this,other.alpha1),
   n1("n1",this,other.n1),
   alpha2("alpha2",this,other.alpha2),
   n2("n2",this,other.n2)

 { 
 } 
 
 double RooDoubleCBFast::evaluate() const 
 { 
   double t = (x-mean)/width;
   double val = -99.;
   if(t>-alpha1 && t<alpha2){
     val = TMath::Exp(-0.5*t*t);
   }else if(t<=-alpha1){
     double alpha1invn1 = alpha1/n1;
     val = TMath::Exp(-0.5*alpha1*alpha1)*TMath::Power(1. - alpha1invn1*(alpha1+t), -n1);     
   }else if(t>=alpha2){
     double alpha2invn2 = alpha2/n2;
     val = TMath::Exp(-0.5*alpha2*alpha2)*TMath::Power(1. - alpha2invn2*(alpha2-t), -n2);     
   } 
   return val;
    
 } 

 Int_t RooDoubleCBFast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const 
 {
   if (matchArgs(allVars,analVars,x)) return 1;
   return 0;
 }

 Double_t RooDoubleCBFast::analyticalIntegral(Int_t code, const char* rangeName) const 
 {
   assert(code==1) ;
 
   double central=0;
   double left=0;
   double right=0;
   
   double xmin = x.min(rangeName);
   double xmax = x.max(rangeName);
   
   //   printf("xmin = %5e, xmax = %5e\n",xmin,xmax);
   
   //bool isfullrange = xmin<=-RooNumber::infinity() && xmax>=RooNumber::infinity();
 
   static const double rootPiBy2 = TMath::Sqrt(atan2(0.0,-1.0)/2.0);
   static const double invRoot2 = 1.0/TMath::Sqrt(2);   
   
   double invwidth = 1.0/width;
   
   double tmin = (xmin-mean)*invwidth;
   double tmax = (xmax-mean)*invwidth;
   
   bool isfullrange = (tmin<-1000. && tmax>1000.);
   
   //compute gaussian contribution
   double central_low =std::max(xmin,mean - alpha1*width );
   double central_high=std::min(xmax,mean + alpha2*width );
   
   double tcentral_low = (central_low-mean)*invwidth;
   double tcentral_high = (central_high-mean)*invwidth;
   if(central_low < central_high)  {// is the gaussian part in range?
     central = rootPiBy2*width*(TMath::Erf(tcentral_high*invRoot2)-TMath::Erf(tcentral_low*invRoot2));
   }
   //compute left tail;
   if (isfullrange  && (n1-1.0)>1.e-5) {
     left = width*TMath::Exp(-0.5*alpha1*alpha1)*n1/(alpha1*(n1-1.)); 
   } else {
  
    double left_low=xmin;
    double left_high=std::min(xmax,mean - alpha1*width);
    double thigh = (left_high-mean)*invwidth;
    
    if(left_low < left_high){ //is the left tail in range?
     double n1invalpha1 = n1/(TMath::Abs(alpha1));
      if(TMath::Abs(n1-1.0)>1.e-5) {
	double invn1m1 = 1./(n1-1.);
	double leftpow = TMath::Power(n1invalpha1,-n1*invn1m1);
	double left0 = width*TMath::Exp(-0.5*alpha1*alpha1)*invn1m1;
	double left1, left2;
	
	if (xmax>(mean-alpha1*width)) left1 = n1invalpha1;
	else left1 = TMath::Power( leftpow*(n1invalpha1 - alpha1 - thigh), 1.-n1);
	
	if (tmin<-1000.) left2 = 0.;
	else left2 = TMath::Power( leftpow*(n1invalpha1 - alpha1 - tmin ), 1.-n1);
	
	left = left0*(left1-left2);
	
      }
      else {
	double A1 = TMath::Power(n1invalpha1,n1)*TMath::Exp(-0.5*alpha1*alpha1);
	double B1 = n1invalpha1-TMath::Abs(alpha1);	
	left = A1*width*(TMath::Log(B1-(left_low-mean)*invwidth) - TMath::Log(B1-(left_high-mean)*invwidth) );
      }
    }
   }
 
   //compute right tail;
   if (isfullrange && (n2-1.0)>1.e-5) {
     right = width*TMath::Exp(-0.5*alpha2*alpha2)*n2/(alpha2*(n2-1.));
   }
   else {    
    double right_low=std::max(xmin,mean + alpha2*width);
    double right_high=xmax;
    double tlow = (right_low - mean)*invwidth;
    
    if(right_low < right_high){ //is the right tail in range?
      double n2invalpha2 = n2/(TMath::Abs(alpha2)); 
      if(TMath::Abs(n2-1.0)>1.e-5) {
	double invn2m2 = 1./(n2-1.);
	double rightpow = TMath::Power(n2invalpha2,-n2*invn2m2);
	double right0 = width*TMath::Exp(-0.5*alpha2*alpha2)*invn2m2;
	double right1, right2;
	
	if (xmin<(mean+alpha2*width)) right1 = n2invalpha2;
	else right1 = TMath::Power( rightpow*(n2invalpha2 - alpha2 + tlow), 1.-n2);
	
	if (tmax>1000.) right2 = 0.;
	else right2 = TMath::Power( rightpow*(n2invalpha2 - alpha2 + tmax), 1.-n2);
	
	right = right0*(right1-right2);	
	
      }
      else {
	double A2 = TMath::Power(n2invalpha2,n2)*TMath::Exp(-0.5*alpha2*alpha2);
	double B2 = n2invalpha2-TMath::Abs(alpha2);
	right = A2*width*(TMath::Log(B2+(right_high-mean)*invwidth) - TMath::Log(B2+(right_low-mean)*invwidth) );
      }
    }
   }
   
   double sum = left + central + right;
   
   //if (!std::isnormal(left) || !std::isnormal(central) || !std::isnormal(right)) {
   if (!std::isnormal(sum)) {
     // printf("bad int: mean = %5f, sigma = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f\n",double(mean),double(width),double(alpha1),double(n1),double(alpha2), double(n2));
     //// printf("left = %5f, central = %5f, right = %5f, A1 = %5f, B1 = %5f, A2 = %5f, B2 = %5f, integral = %5f\n",left,central,right,A1,B1,A2,B2,left+central+right);
     // printf("left = %5f, central = %5f, right = %5f, integral = %5f\n",left,central,right,sum);
   }
     
   return sum;
 
 }
 
