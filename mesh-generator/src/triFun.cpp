#include <iostream>
using namespace std;


double evalTriPhi(const int i, const int j, const double x, const double y, const int n) {
 double phi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i+j>n)
    cout<<"Error! test function indices i+j > n \n";
  else{
    double h;
    int np=n+1;
    double prod1=1.; 
    double prod2=1.; 
    double prod3=1.;
    int k=n-i-j;
/********** Create Node Matrix *************/
     
    double* nodex = new double[np];
    double* nodey = new double[np];
    h=1./n;
    nodex[0]=0.;
    nodey[0]=0.;
    for (int l=1; l<np; l++){
      nodex[l]=nodex[l-1]+h;
      nodey[l]=nodey[l-1]+h;
    }
    
/********* Evaluate phi ********************/
    for (int l=0;l<i;l++)
	  prod1*=(x-nodex[l])/(nodex[i]-nodex[l]); 
    for (int l=0;l<j;l++)
      prod2*=(y-nodey[l])/(nodey[j]-nodey[l]); 
    for (int l=0;l<k;l++)
      prod3*=(x+y-nodex[n-l])/(nodex[i]+nodey[j]-nodex[n-l]); 
   phi = prod1*prod2*prod3;
    delete [] nodex;
    delete [] nodey;
  }
  return phi;
}

double evalTriDPhidx(const int i, const int j, const double x, const double y, const int n) {

  double dphi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i+j>n)
    cout<<"Error! test function indices i+j > n \n";
  else{
    double h;
    int np=n+1;
    double prod1=1.; 
    double prod2=1.; 
    double prod3=1.;
    int k=n-i-j;
    
/********** Create Node Matrix *************/
    double* nodex = new double[np];
    double* nodey = new double[np];
    h=1./n;
    nodex[0]=0.;
    nodey[0]=0.;
    for (int l=1; l<np; l++){
      nodex[l]=nodex[l-1]+h;
      nodey[l]=nodey[l-1]+h;
    }
    
/********** Evaluate dphidx *************/
    for (int l=0;l<i;l++)
	  prod1*=(x-nodex[l])/(nodex[i]-nodex[l]); 
    for (int l=0;l<j;l++)
      prod2*=(y-nodey[l])/(nodey[j]-nodey[l]); 
    for (int l=0;l<k;l++)
      prod3*=(x+y-nodex[n-l])/(nodex[i]+nodey[j]-nodex[n-l]); 

    double p1;
    double dprod1dx=0;
    for (int m=0;m<i;m++){
        p1=1./(nodex[i]-nodex[m]);
	for(int l=0;l<i;l++)
	  if(l!=m)
	    p1*=(x-nodex[l])/(nodex[i]-nodex[l]);
	dprod1dx+=p1; 
    }
    
    double p3;
    double dprod3dx=0;
    for (int m=0;m<k;m++){
        p3=1./(nodex[i]+nodey[j]-nodex[n-m]);
	for(int l=0;l<k;l++)
	  if(l!=m)
	    p3*=(x+y-nodex[n-l])/(nodex[i]+nodey[j]-nodex[n-l]);
	dprod3dx+=p3; 
    }
    
    dphi=prod2*(prod3*dprod1dx+prod1*dprod3dx);
    
    delete [] nodex;
    delete [] nodey;
  }

  return dphi;
}

double evalTriDPhidy(const int i, const int j, const double x, const double y, const int n) {

  double dphi=1.;
  if(n<1)
    cout<<"Error! number of element segments < 1 \n";
  else if(i+j>n)
    cout<<"Error! test function indices i+j > n \n";
  else{
    double h;
    int np=n+1;
    double prod1=1.; 
    double prod2=1.; 
    double prod3=1.;
    int k=n-i-j;
    
/********** Create Node Matrix *************/
    double* nodex = new double[np];
    double* nodey = new double[np];
    h=1./n;
    nodex[0]=0.;
    nodey[0]=0.;
    for (int l=1; l<np; l++){
      nodex[l]=nodex[l-1]+h;
      nodey[l]=nodey[l-1]+h;
    }
    
/********** Evaluate dphidy *************/
    for (int l=0;l<i;l++)
	  prod1*=(x-nodex[l])/(nodex[i]-nodex[l]); 
    for (int l=0;l<j;l++)
      prod2*=(y-nodey[l])/(nodey[j]-nodey[l]); 
    for (int l=0;l<k;l++)
      prod3*=(x+y-nodex[n-l])/(nodex[i]+nodey[j]-nodex[n-l]); 
    
    double p2;
    double dprod2dy=0;
    for (int m=0;m<j;m++){
        p2=1./(nodey[j]-nodey[m]);
	for(int l=0;l<j;l++)
	  if(l!=m)
	    p2*=(y-nodey[l])/(nodey[j]-nodey[l]);
	dprod2dy+=p2; 
    }
    
    double p3;
    double dprod3dy=0;
    for (int m=0;m<k;m++){
        p3=1./(nodex[i]+nodey[j]-nodex[n-m]);
	for(int l=0;l<k;l++)
	  if(l!=m)
	    p3*=(x+y-nodex[n-l])/(nodex[i]+nodey[j]-nodex[n-l]);
	dprod3dy+=p3; 
    }
    
    dphi=prod1*(prod2*dprod3dy+prod3*dprod2dy);
    
    delete [] nodex;
    delete [] nodey;
  }

  return dphi;
}


