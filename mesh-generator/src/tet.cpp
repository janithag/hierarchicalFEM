#include "basis.hpp"
#include "main.hpp"

double tet::eval_phi(const int* I,const double* X, const int n) const{
  return evalTetraPhi(I[0],I[1],I[2], X[0], X[1], X[2], n);
}

double tet::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalTetraDPhidx(I[0],I[1],I[2], X[0], X[1], X[2], n);
}

double tet::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalTetraDPhidy(I[0],I[1],I[2], X[0], X[1], X[2], n);
}

double tet::eval_dphidz(const int* I,const double* X, const int n) const{
  return evalTetraDPhidz(I[0],I[1],I[2], X[0], X[1], X[2], n);
}

int tet::ElementDOF(const int p) const{
  return (p+1)*(p+2)*(p+3)/6;
}

void tet::RefinedNodes(const int p, double** x) const{ // obtained by mapping unit tet to the reference tet in hierarchic basis

  double h = 1./p;
  double xi, eta, zeta; // coordinates of unit tet
  int count=0;
  
  for(int i=0;i<p+1;i++)
    for(int j=0;j<p+1-i;j++)
      for(int k=0;k<p+1-i-j;k++){
	xi = k*h;
	eta = j*h;
	zeta= i*h;
	x[count][0] = 2*xi+eta+zeta-1;
	x[count][1] = sqrt(3)*eta+zeta/sqrt(3);
	x[count][2] = sqrt(8./3)*zeta;
	count++;
    }
  
return;
}