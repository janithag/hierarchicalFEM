#include "main.hpp"
#include <iostream>
using namespace std;

double  evalPhi(const int i, const double x, const int n) {

  double phi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i>n)
    cout<<"Error! test function index i > n \n";
  else{
    double h;
    
    int np=n+1;
    double* node = new double [np];
    
    h=2./n;
    node[0]=-1.;
    
    for (int j=1; j<np; j++){
      node[j]=node[j-1]+h;
    }
    
    for (int j=0;j<np;j++)
      if (j!=i)
	phi*=(x-node[j])/(node[i]-node[j]); 
    delete [] node;
  }
  return phi;
}


double  evalDphi(const int i, const double x, const int n) {

  double dphidx=0.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i>n)
    cout<<"Error! test function index i  > n \n";
  else{
    double h;
    int np = n+1;
    double* node = new double[np];
  
    h=2./n;
    node[0]=-1.;
    
    for (int j=1; j<np; j++){
      node[j]=node[j-1]+h;
    }
    
    double prod;
    for (int k=0;k<np;k++){
      if (k!=i){
        prod=1./(node[i]-node[k]);
	for(int j=0;j<np;j++)
	  if( (j!=k) && (j!=i) )
	    prod*=(x-node[j])/(node[i]-node[j]);
	dphidx+=prod; 
      }
    } 
    
    delete [] node;
    
  }
  
  return dphidx;
}
