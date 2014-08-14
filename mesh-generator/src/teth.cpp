#include "basish.hpp"

double teth::eval_phi(const int* I,const double* X, const int n) const{
  return evalTethPhi(I[0],I[1],I[2], X[0], X[1], X[2], n, type);
}

double teth::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalTethDPhidx(I[0],I[1],I[2], X[0], X[1], X[2], n, type);
}

double teth::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalTethDPhidy(I[0],I[1],I[2], X[0], X[1], X[2], n, type);
}

double teth::eval_dphidz(const int* I,const double* X, const int n) const{
  return evalTethDPhidz(I[0],I[1],I[2], X[0], X[1], X[2], n, type);
}

void teth::eval_sign(const int* o, vector <int>& orientationFlag) const{ 
 for(int i=0;i<dimOfBasis();i++)
    orientationFlag[i]=1;
}

void teth::eval_type(const int* referenceElementType, vector <int>& type) const{ 
  
  for(int i=0;i<dimOfBasis();i++)
    type[i] = referenceElementType[0];
  
}
