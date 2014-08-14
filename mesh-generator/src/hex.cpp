#include "basis.hpp"
//#include "main.hpp"

double hex::eval_phi(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalPhi(I[1], X[1], n)*evalPhi(I[2], X[2], n);
}

double hex::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalDphi(I[0], X[0], n)*evalPhi(I[1], X[1], n)*evalPhi(I[2], X[2], n);
}

double hex::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalDphi(I[1], X[1], n)*evalPhi(I[2], X[2], n);
}

double hex::eval_dphidz(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalPhi(I[1], X[1], n)*evalDphi(I[2], X[2], n);
}

int hex::ElementDOF(const int p) const{
  return pow(p+1,3);
}

void hex::RefinedNodes(const int p, double** x) const{
    double h = 2./p;

  for(int i=0;i<p+1;i++)
    for(int j=0;j<p+1;j++)
      for(int k=0;k<p+1;k++){
	x[i*(p+1)*(p+1)+j*(p+1)+k][0]=-1+k*h;
	x[i*(p+1)*(p+1)+j*(p+1)+k][1]=-1+j*h;
	x[i*(p+1)*(p+1)+j*(p+1)+k][2]=-1+i*h;
      }
return;
}



