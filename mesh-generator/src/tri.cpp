#include "basis.hpp"
//#include "main.hpp"

double tri::eval_phi(const int* I,const double* X, const int n) const{
  return evalTriPhi(I[0], I[1], X[0], X[1], n);
}

double tri::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalTriDPhidx(I[0], I[1], X[0], X[1], n);
}

double tri::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalTriDPhidy(I[0], I[1], X[0], X[1], n);
}

double tri::eval_dphidz(const int* I,const double* X, const int n) const{
  return 0; 
}

int tri::ElementDOF(const int p) const{
  return (p+1)*(p+2)/2;
}

void tri::RefinedNodes(const int p, double** x) const{ // This is for equilateral triangle used in hierarchical basis.
  double hx = 2./p;
  double hy = sqrt(3)/p;
  int count=0;
  
  for(int i=0;i<p+1;i++)
    for(int j=0;j<p+1-i;j++){
      x[count][0] = (i*hy-sqrt(3))/sqrt(3) + j*hx;
      x[count][1] = i*hy;
      x[count][2] = 0.;
      count++;
    }
  return;
}

