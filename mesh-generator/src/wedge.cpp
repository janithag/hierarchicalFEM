#include "basis.hpp"
//#include "main.hpp"

double wedge::eval_phi(const int* I,const double* X, const int n) const{
  return evalTriPhi(I[0], I[1], X[0], X[1], n)*evalPhi(I[2], X[2], n); 
}

double wedge::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalTriDPhidx(I[0], I[1], X[0], X[1], n)*evalPhi(I[2], X[2], n); 
}

double wedge::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalTriDPhidy(I[0], I[1], X[0], X[1], n)*evalPhi(I[2], X[2], n); 
}

double wedge::eval_dphidz(const int* I,const double* X, const int n) const{
  return evalTriPhi(I[0], I[1], X[0], X[1], n)*evalDphi(I[2], X[2], n); 
}

int wedge::ElementDOF(const int p) const{
  return pow(p+1,2)*(p+2)/2;;
}

void wedge::RefinedNodes(const int p, double** x) const{
return;
}

////************************************************************
//double  wedgeth::eval_phi(const int *I,const double* x) const{ 
//  return wed_th(x[0],x[1],x[2],I[0],I[1],I[2]);
//}
//
//double  wedgeth::eval_dphidx(const int *I,const double* x) const{ 
//  return dwed_thdx(x[0],x[1],x[2],I[0],I[1],I[2]);
//}
//
//double  wedgeth::eval_dphidy(const int *I,const double* x) const{ 
//  return dwed_thdy(x[0],x[1],x[2],I[0],I[1],I[2]);
//}
//
//double  wedgeth::eval_dphidz(const int *I,const double* x) const{ 
//  return dwed_thdz(x[0],x[1],x[2],I[0],I[1],I[2]);
//}
//
//double wedgeth::wed_th(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const{
//  double f=0.;
//  double t=x+y;
//  switch(i){
//  case 0:
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=(1-t)*(-2*t-z)*(1-z)*0.5;
//	break;
//      case 1:
//	f=(1-t)*(1-z*z);
//	break;
//      case 2:
//	f=(1-t)*(-2*t+z)*(1+z)*0.5;
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=2*y*(1-z)*(1-t);
//	break;
//      case 2:
//	f=2*y*(1+z)*(1-t);
//	break;
//      }
//      break;
//    case 2:
//      switch(k){
//      case 0:
//	f=y*(-2+2*y-z)*(1-z)*0.5;
//	break;
//      case 1:
//	f=y*(1-z*z);
//	break;
//      case 2:
//	f=y*(-2+2*y+z)*(1+z)*0.5;
//	break;
//      }
//      break;
//    }
//    break;
//  case 1: 
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=2*x*(1-z)*(1-t);
//	break;
//      case 2:
//	f=2*x*(1+z)*(1-t);
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=2*x*y*(1-z);
//	break;
//      case 2:
//	f=2*x*y*(1+z);
//	break;
//      }
//      break;
//    }
//    break;
//  case 2:
//    switch(k){
//    case 0:
//      f=x*(-2+2*x-z)*(1-z)*0.5;
//      break; 
//    case 1:
//      f=x*(1-z*z);
//      break;
//    case 2:
//      f=x*(-2+2*x+z)*(1+z)*0.5;
//      break;
//    }
//    break;
//  }
//  return f;
//}

//
//double wedgeth::dwed_thdx(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const{
//  double f=0.;
//  double t=x+y;
//  switch(i){
//  case 0:
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=(-1.+2.*t+0.5*z)*(1.-z);
//	break;
//      case 1:
//	f=(z*z-1.);
//	break;
//      case 2:
//	f=(-1.+2.*t-0.5*z)*(1.+z);
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=-2.*y*(1-z);
//	break;
//      case 2:
//	f=-2.*y*(1+z);
//	break;
//      }
//      break;
//    case 2:
//      switch(k){
//      case 0:
//	f=0.;
//	break;
//      case 1:
//	f=0.;
//	break;
//      case 2:
//	f=0.;
//	break;
//      }
//      break;
//    }
//    break;
//  case 1: 
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=2*(1.-z)*(1.-2.*x-y);
//	break;
//      case 2:
//	f=2*(1.+z)*(1.-2.*x-y);
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=2.*y*(1.-z);
//	break;
//      case 2:
//	f=2.*y*(1.+z);
//	break;
//      }
//      break;
//    }
//    break;
//  case 2:
//    switch(k){
//    case 0:
//      f=(-2.+4.*x-z)*(1.-z)*0.5;
//      break; 
//    case 1:
//      f=(1.-z*z);
//      break;
//    case 2:
//      f=(-2.+4.*x+z)*(1+z)*0.5;
//      break;
//    }
//    break;
//  }
//  return f;
//}
//
//
//double wedgeth::dwed_thdy(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const{
//  double f=0.;
//  double t=x+y;
//  switch(i){
//  case 0:
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=(-1.+2.*t+0.5*z)*(1.-z);
//	break;
//      case 1:
//	f=(z*z-1.);
//	break;
//      case 2:
//	f=(-1.+2.*t-0.5*z)*(1.+z);
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=2.*(1-z)*(1.-x-2.*y);
//	break;
//      case 2:
//	f=2.*(1+z)*(1.-x-2.*y);
//	break;
//      }
//      break;
//    case 2:
//      switch(k){
//      case 0:
//	f=(-2.+4.*y-z)*(1.-z)*0.5;
//	break;
//      case 1:
//	f=(1.-z*z);
//	break;
//      case 2:
//	f=(-2.+4.*y+z)*(1.+z)*0.5;
//	break;
//      }
//      break;
//    }
//    break;
//  case 1: 
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=-2.*x*(1.-z);
//	break;
//      case 2:
//	f=-2.*x*(1.+z);
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=2.*x*(1.-z);
//	break;
//      case 2:
//	f=2.*x*(1.+z);
//	break;
//      }
//      break;
//    }
//    break;
//  case 2:
//    switch(k){
//    case 0:
//      f=0.;
//      break; 
//    case 1:
//      f=0.;
//      break;
//    case 2:
//      f=0.;
//      break;
//    }
//    break;
//  }
//  return f;
//}
//
//
//double wedgeth::dwed_thdz(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const{
//  double f=0.;
//  double t=1.-(x+y);
//  switch(i){
//  case 0:
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=t*(0.5-t+z);
//	break;
//      case 1:
//	f=-2.*t*z;
//	break;
//      case 2:
//	f=t*(-0.5+t+z);
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=-2.*y*t;
//	break;
//      case 2:
//	f=2.*y*t;
//	break;
//      }
//      break;
//    case 2:
//      switch(k){
//      case 0:
//	f=y*(0.5-y+z);
//	break;
//      case 1:
//	f=-2.*y*z;
//	break;
//      case 2:
//	f=y*(-0.5+y+z);
//	break;
//      }
//      break;
//    }
//    break;
//  case 1: 
//    switch(j){
//    case 0:
//      switch(k){
//      case 0:
//	f=-2.*x*t;
//	break;
//      case 2:
//	f=2.*x*t;
//	break;
//      }
//      break;
//    case 1:
//      switch(k){
//      case 0:
//	f=-2.*x*y;
//	break;
//      case 2:
//	f=2.*x*y;
//	break;
//      }
//      break;
//    }
//    break;
//  case 2:
//    switch(k){
//    case 0:
//      f=x*(0.5-x+z);
//      break; 
//    case 1:
//      f=-2.*x*z;
//      break;
//    case 2:
//      f=x*(-0.5+x+z);
//      break;
//    }
//    break;
//  }
//  return f;
//}


