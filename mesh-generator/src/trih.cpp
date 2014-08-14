#include "basish.hpp"
#include <cmath>
#include <iostream>
double trih::eval_phi(const int* I,const double* X, const int n) const{
  return evalTrihPhi(I[0], I[1], I[2], X[0], X[1], n);
}

double trih::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalTrihDPhidx(I[0], I[1], I[2], X[0], X[1], n);
}

double trih::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalTrihDPhidy(I[0], I[1], I[2], X[0], X[1], n);
}

double trih::eval_dphidz(const int* I,const double* X, const int n) const{
  return 0; 
}

void trih::eval_sign(const int* o, vector <int>& orientationFlag) const{ 

   unsigned count=0;

  // vertices
  for(int i=0;i<3;i++){
    orientationFlag[count] = 1;
    count++;
  }
   
  for(int k=0; k<p-1;k++){
  // edges 
  for(int i=0;i<3;i++){
	orientationFlag[count] = pow(o[i],k);
	count++;
  }
  //interior modes
  if(k>0)
    for(int j=0; j<k; j++){
    orientationFlag[count] = 1;
    count++;
    } 
}

}

void trih::eval_type(const int* referenceElementType, vector <int>& type) const{ 
  for(int i=0;i<dimOfBasis();i++)
    type[i] = 0;
}
