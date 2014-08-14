//#include "basis.hpp"
#include "main.hpp"
#include <iostream>
using namespace std;

double evalTetraPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n) {
 double phi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i+j+k>n)
    cout<<"Error! test function indices i+j+k > n \n";
  else{
    double h;
    int np=n+1;
    int r=n-i-j-k;
    
/********** Create Nodes *************/
    double* nodex = new double[np];
    double* nodey = new double[np];
    double* nodez = new double[np];
    h=1./n;
    nodex[0]=0.;
    nodey[0]=0.;
    nodez[0]=0.;
    for (int l=1; l<np; l++){
      nodex[l]=nodex[l-1]+h;
      nodey[l]=nodey[l-1]+h;
      nodez[l]=nodez[l-1]+h;
    }
    
/********* Evaluate phi ********************/
    double prodx=1.; 
    double prody=1.; 
    double prodz=1.;
    double prodxyz=1.;
    
    for (int l=0;l<i;l++)
	  prodx*=(x-nodex[l])/(nodex[i]-nodex[l]); 
    for (int l=0;l<j;l++)
      prody*=(y-nodey[l])/(nodey[j]-nodey[l]); 
    for (int l=0;l<k;l++)
      prodz*=(z-nodez[l])/(nodez[k]-nodez[l]);
    for (int l=0;l<r;l++)
      prodxyz*=(x+y+z-nodex[n-l])/(nodex[i]+nodey[j]+nodez[k]-nodex[n-l]); 
    phi = prodx*prody*prodz*prodxyz;
    
    delete [] nodex;
    delete [] nodey;
    delete [] nodez;
  }
  return phi;
}

double evalTetraDPhidx(const int i, const int j, const int k, const double x, const double y, const double z, const int n){

  double dphi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i+j+k>n)
    cout<<"Error! test function indices i+j+k > n \n";
  else{
    double h;
    int np=n+1;
    double prodx=1.; 
    double prody=1.; 
    double prodz=1.;
    double prodxyz=1.;
    int r=n-i-j-k;
    
/********** Create Nodes *************/
    double* nodex = new double[np];
    double* nodey = new double[np];
    double* nodez = new double[np];
    h=1./n;
    nodex[0]=0.;
    nodey[0]=0.;
    nodez[0]=0.;
    for (int l=1; l<np; l++){
      nodex[l]=nodex[l-1]+h;
      nodey[l]=nodey[l-1]+h;
      nodez[l]=nodez[l-1]+h;
    }
    
/********** Evaluate dphidx *************/
    for (int l=0;l<i;l++)
	  prodx*=(x-nodex[l])/(nodex[i]-nodex[l]); 
    for (int l=0;l<j;l++)
      prody*=(y-nodey[l])/(nodey[j]-nodey[l]); 
    for (int l=0;l<k;l++)
      prodz*=(z-nodez[l])/(nodez[k]-nodez[l]);
    for (int l=0;l<r;l++)
      prodxyz*=(x+y+z-nodex[n-l])/(nodex[i]+nodey[j]+nodez[k]-nodex[n-l]); 
      
    double p1;
    double dprodxdx=0;
    for (int m=0;m<i;m++){
        p1=1./(nodex[i]-nodex[m]);
	for(int l=0;l<i;l++)
	  if(l!=m)
	    p1*=(x-nodex[l])/(nodex[i]-nodex[l]);
	dprodxdx+=p1; 
    }
    
    double p3;
    double dprodxyzdx=0;
    for (int m=0;m<r;m++){
        p3=1./(nodex[i]+nodey[j]+nodez[k]-nodex[n-m]);
	for(int l=0;l<r;l++)
	  if(l!=m)
	    p3*=(x+y+z-nodex[n-l])/(nodex[i]+nodey[j]+nodez[k]-nodex[n-l]);
	dprodxyzdx+=p3; 
    }
    
    dphi=prody*prodz*(prodxyz*dprodxdx+prodx*dprodxyzdx);
    
    delete [] nodex;
    delete [] nodey;
    delete [] nodez;
  }

  return dphi;
}

double evalTetraDPhidy(const int i, const int j, const int k, const double x, const double y, const double z, const int n){

  double dphi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i+j+k>n)
    cout<<"Error! test function indices i+j+k > n \n";
  else{
    double h;
    int np=n+1;
    double prodx=1.; 
    double prody=1.; 
    double prodz=1.;
    double prodxyz=1.;
    int r=n-i-j-k;
    
/********** Create Nodes *************/
    double* nodex = new double[np];
    double* nodey = new double[np];
    double* nodez = new double[np];
    h=1./n;
    nodex[0]=0.;
    nodey[0]=0.;
    nodez[0]=0.;
    for (int l=1; l<np; l++){
      nodex[l]=nodex[l-1]+h;
      nodey[l]=nodey[l-1]+h;
      nodez[l]=nodez[l-1]+h;
    }
    
/********** Evaluate dphidx *************/
    for (int l=0;l<i;l++)
	  prodx*=(x-nodex[l])/(nodex[i]-nodex[l]); 
    for (int l=0;l<j;l++)
      prody*=(y-nodey[l])/(nodey[j]-nodey[l]); 
    for (int l=0;l<k;l++)
      prodz*=(z-nodez[l])/(nodez[k]-nodez[l]);
    for (int l=0;l<r;l++)
      prodxyz*=(x+y+z-nodex[n-l])/(nodex[i]+nodey[j]+nodez[k]-nodex[n-l]); 
      
    double p2;
    double dprodydy=0;
    for (int m=0;m<j;m++){
        p2=1./(nodey[j]-nodey[m]);
	for(int l=0;l<j;l++)
	  if(l!=m)
	    p2*=(y-nodey[l])/(nodey[j]-nodey[l]);
	dprodydy+=p2; 
    }
    
    double p3;
    double dprodxyzdx=0;
    for (int m=0;m<r;m++){
        p3=1./(nodex[i]+nodey[j]+nodez[k]-nodex[n-m]);
	for(int l=0;l<r;l++)
	  if(l!=m)
	    p3*=(x+y+z-nodex[n-l])/(nodex[i]+nodey[j]+nodez[k]-nodex[n-l]);
	dprodxyzdx+=p3; 
    }
    
    dphi=prodx*prodz*(prodxyz*dprodydy+prody*dprodxyzdx);
    
    delete [] nodex;
    delete [] nodey;
    delete [] nodez;
  }

  return dphi;
}

double evalTetraDPhidz(const int i, const int j, const int k, const double x, const double y, const double z, const int n){

  double dphi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i+j+k>n)
    cout<<"Error! test function indices i+j+k > n \n";
  else{
    double h;
    int np=n+1;
    double prodx=1.; 
    double prody=1.; 
    double prodz=1.;
    double prodxyz=1.;
    int r=n-i-j-k;
    
/********** Create Nodes *************/
    double* nodex = new double[np];
    double* nodey = new double[np];
    double* nodez = new double[np];
    h=1./n;
    nodex[0]=0.;
    nodey[0]=0.;
    nodez[0]=0.;
    for (int l=1; l<np; l++){
      nodex[l]=nodex[l-1]+h;
      nodey[l]=nodey[l-1]+h;
      nodez[l]=nodez[l-1]+h;
    }
    
/********** Evaluate dphidx *************/
    for (int l=0;l<i;l++)
	  prodx*=(x-nodex[l])/(nodex[i]-nodex[l]); 
    for (int l=0;l<j;l++)
      prody*=(y-nodey[l])/(nodey[j]-nodey[l]); 
    for (int l=0;l<k;l++)
      prodz*=(z-nodez[l])/(nodez[k]-nodez[l]);
    for (int l=0;l<r;l++)
      prodxyz*=(x+y+z-nodex[n-l])/(nodex[i]+nodey[j]+nodez[k]-nodex[n-l]); 

    double p3;
    double dprodxyzdx=0;
    for (int m=0;m<r;m++){
        p3=1./(nodex[i]+nodey[j]+nodez[k]-nodex[n-m]);
	for(int l=0;l<r;l++)
	  if(l!=m)
	    p3*=(x+y+z-nodex[n-l])/(nodex[i]+nodey[j]+nodez[k]-nodex[n-l]);
	dprodxyzdx+=p3; 
    }
    
    double p4;
    double dprodzdz=0;
    for (int m=0;m<k;m++){
        p4=1./(nodez[k]-nodez[m]);
	for(int l=0;l<k;l++)
	  if(l!=m)
	    p4*=(z-nodez[l])/(nodez[k]-nodez[l]);
	dprodzdz+=p4; 
    }
    
    dphi=prodx*prody*(prodxyz*dprodzdz+prodz*dprodxyzdx);
    
    delete [] nodex;
    delete [] nodey;
    delete [] nodez;
  }

  return dphi;
}
