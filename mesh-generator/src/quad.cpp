#include "basis.hpp"
//#include <iostream>
//#include "main.hpp"

double quad::eval_phi(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalPhi(I[1], X[1], n);
}

double quad::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalDphi(I[0], X[0], n)*evalPhi(I[1], X[1], n);
}

double quad::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalDphi(I[1], X[1], n);
}

double quad::eval_dphidz(const int* I,const double* X, const int n) const{
  return 0.;
}

int quad::ElementDOF(const int p) const{ // p is 1-based
  return (p+1)*(p+1);
}

void quad::RefinedNodes(const int p, double** x) const{ // x first then y, bottom to top
  
  double h = 2./p;

  for(int i=0;i<p+1;i++)
    for(int j=0;j<p+1;j++){
      x[i*(p+1)+j][0]=-1+j*h;
      x[i*(p+1)+j][1]=-1+i*h;
      x[i*(p+1)+j][2]=0.;
    }
  return;
}


