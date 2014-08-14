#include "basish.hpp"

double  lineh::eval_phi(const int *I,const double* X, const int n) const{
  return evalLinehPhi(I[0], I[1], I[2], X[0], n);
}

double lineh::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalLinehDPhi(I[0], I[1], I[2], X[0], n);
} 
     
double lineh::eval_dphidy(const int *I,const double* X, const int n) const{
  return 0.;
}

double lineh::eval_dphidz(const int *I,const double* X, const int n) const{
  return 0.;
}

void lineh::eval_sign(const int* o, vector <int>& orientationFlag) const{ 
  for(int i=0;i<dimOfBasis();i++)
    orientationFlag[i]=1;
}

void lineh::eval_type(const int* referenceElementType, vector <int>& type) const{ 
  for(int i=0;i<dimOfBasis();i++)
    type[i] = 1;
}