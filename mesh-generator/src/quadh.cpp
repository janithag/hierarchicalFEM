#include "basish.hpp"
#include <cmath>

double quadh::eval_phi(const int* I,const double* X, const int n) const{
  return evalQuadhPhi(I[0], I[1], I[2], X[0], X[1], n);
}

double quadh::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalQuadhDPhiDx(I[0], I[1], I[2], X[0], X[1], n);
}

double quadh::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalQuadhDPhiDy(I[0], I[1], I[2], X[0], X[1], n);
}

double quadh::eval_dphidz(const int* I,const double* X, const int n) const{
  return 0.;
}

/* Evaluate the sign as -1 or +1 based on the orientationFlag of the spatial element. */
void quadh::eval_sign(const int* o, vector <int>& orientationFlag) const{ 

int count=0;
  
//vertices
  for(int i=0;i<4;i++){
    orientationFlag[count] = 1;
  count++;
  }
  
  for(int k=0;k<p-1;k++){
//edges
    for(int i=0;i<4;i++){
	orientationFlag[count] = pow(o[i],k);
	count++;
      }
//interior modes
  if(k>1){
    for(int i=0;i<k-1;i++){    
      orientationFlag[count] = 1;
      count++;
    }
  
  }
 }
}

void quadh::eval_type(const int* referenceElementType, vector <int>& type) const{ 
  for(int i=0;i<dimOfBasis();i++)
    type[i] = 0;
}


