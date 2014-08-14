#include "basis.hpp"
#include <iostream>
using namespace std;

double  line::eval_phi(const int *I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n);
}
double line::eval_dphidx(const int *I,const double* X, const int n) const{
  return evalDphi(I[0], X[0], n);
}      
double line::eval_dphidy(const int *I,const double* X, const int n) const{
  return 0.;
}

double line::eval_dphidz(const int *I,const double* X, const int n) const{
  return 0.;
}

int line::ElementDOF(const int p) const{
  return p+1;
}

void line::RefinedNodes(const int p, double** x) const{
  
int localDOF=p+1;
double h=2./p;

for(int i=0;i<localDOF;i++){
  x[i][0]=-1+i*h;
  x[i][1]=0.;
  x[i][2]=0.;
}

return;
  
}
/////////////////////////////

void line::eval_phi(const double x, double *phi, const int n) const{
  
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else{
    double h;
    
    int np=n+1;
    double* node = new double [np];
    
    h=2./n;
    node[0]=-1.;
    for (int i=1; i<np; i++){
      node[i]=node[i-1]+h;
    }
    
    for(int i=0;i<np;i++){ 
      phi[i]=1;
      for (int j=0;j<np;j++)
	if (j!=i)
	  phi[i]*=(x-node[j])/(node[i]-node[j]); 
    }
      
    delete [] node;
  }
  return;
}


void line::eval_dphidx(const double x, double* dphidx, const int n) const{//n: order of polynomial     
  
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else{
    double h;
    int np = n+1;
    double* node = new double[np];
    
    h=2./n;
    node[0]=-1.;
    for (int i=1; i<np; i++){
      node[i]=node[i-1]+h;
    }
    
    for(int i=0;i<np;i++){ 
      dphidx[i]=0.;
      double prod;
      for (int k=0;k<np;k++){
	if (k!=i){
	  prod=1./(node[i]-node[k]);
	  for(int j=0;j<np;j++)
	    if( (j!=k) && (j!=i) )
	      prod*=(x-node[j])/(node[i]-node[j]);
	  dphidx[i]+=prod; 
	}
      } 
    }
    
    delete [] node;
    
  }
  
  return;
}
