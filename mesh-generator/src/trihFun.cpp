// This files generates the hierarchic shape functions for triangle, tetrahedron and pentahedron.
// i: type (vertex:0, edge:1, face:2, internal:3), j: location, k: basis function number (we use 0-based index)
// p is 1-based.

#include <stdexcept>
#include <iostream>
#include <cmath>
#include "hFun.hpp"
#include <algorithm>

using std::runtime_error;
using std::cout;
using std::swap;
using std::endl;

// This function generates heirarchic basis functions for triangle.
double evalTrihPhi(const int i, const int j, const int k, const double x, const double y, const int p) {

  double phi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>2) throw runtime_error("Invalid type of shape functions for triangle");

  if ((i==0)&&(k==0))// vertex
    phi = L(j, x, y);
   
  else if(i==1){// edge  
    if ((k>p-2)||(p<2)||(j>2))throw runtime_error("Invalid for edge shape functions in triangle");
    int j2=(j+1)%3; 
    phi = L(j,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j,x,y));
  }       
   
  else if (i==2){// face 

    if((p<3) || (k>(p-1)*(p-2)*0.5-1))throw runtime_error("Invalid for face shape functions in triangle");
    int count = 0;
    for(int r=3;r<=p;r++)
      for(int r1=r-3;r1>=0;r1--)
        for(int r2=0;r2<=r-3;r2++)
          if (r1+r2==r-3){
            if (count == k)
              phi = L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1);
            count++;
          }
  }
  else 
    throw runtime_error("Invalid indices for triangle");
  return phi;
}

// This function generates the partial derivatives w.r.t. x for the heirarchic basis functions for triangle.
double evalTrihDPhidx(const int i, const int j, const int k,  const double x, const double y, const int p) {

  double dphi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>2) throw runtime_error("Invalid type of shape functions for triangle");
  
  if ((i==0)&&(k==0))// vertex
    dphi = dLdx(j, x, y);
         
  else if(i==1){ // edge
    if ((k>p-2)||(p<2)||(j>2))throw runtime_error("Invalid for edge shape functions in triangle");
    int j2=(j+1)%3; 
    dphi = dLdx(j,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j,x,y))+L(j,x,y)*dLdx(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j,x,y))+L(j,x,y)*L(j2,x,y)*dkernel(k+1,L(j2,x,y)-L(j,x,y))*(dLdx(j2,x,y)-dLdx(j,x,y));
  }     

  else if (i==2){// face 

    if((p<3) || (k>(p-1)*(p-2)*0.5-1))throw runtime_error("Invalid for face shape functions in triangle");
    int count = 0;
    for(int r=3;r<=p;r++)
      for(int r1=r-3;r1>=0;r1--)
        for(int r2=0;r2<=r-3;r2++)
	  if (r1+r2==r-3){
	    if (count == k)
	      dphi = dLdx(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)+L(0,x,y)*dLdx(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)
		+L(0,x,y)*L(1,x,y)*dLdx(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)+L(0,x,y)*L(1,x,y)*L(2,x,y)*dlegendre(r1,L(1,x,y)-L(0,x,y))*(dLdx(1,x,y)-dLdx(0,x,y))*legendre(r2,2*L(2,x,y)-1)
		+L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*dlegendre(r2,2*L(2,x,y)-1)*(2*dLdx(2,x,y));
	    count++;
	  }
  }

  else 
    throw runtime_error("Invalid indices for triangle");  
  
  return dphi;
}

// This function generates heirarchic basis functions for triangle.
double evalTrihDPhidy(const int i, const int j,  const int k, const double x, const double y, const int p) {

  double dphi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>2) throw runtime_error("Invalid type of shape functions for triangle");
  
  if ((i==0)&&(k==0))// vertex
    dphi = dLdy(j, x, y);

  else if(i==1){// edge  
    if ((k>p-2)||(p<2)||(j>2))throw runtime_error("Invalid for edge shape functions in triangle");
    int j2=(j+1)%3; 
    dphi = dLdy(j,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j,x,y))+L(j,x,y)*dLdy(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j,x,y))+L(j,x,y)*L(j2,x,y)*dkernel(k+1,L(j2,x,y)-L(j,x,y))*(dLdy(j2,x,y)-dLdy(j,x,y));
  }   
         
  else if (i==2){// face 

    if((p<3) || (k>(p-1)*(p-2)*0.5-1))throw runtime_error("Invalid for face shape functions in triangle");
    int count = 0;
    for(int r=3;r<=p;r++)
      for(int r1=r-3;r1>=0;r1--)
        for(int r2=0;r2<=r-3;r2++)
          if (r1+r2==r-3){
            if (count == k)
	      dphi = dLdy(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)+L(0,x,y)*dLdy(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)
		+L(0,x,y)*L(1,x,y)*dLdy(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)+L(0,x,y)*L(1,x,y)*L(2,x,y)*dlegendre(r1,L(1,x,y)-L(0,x,y))*(dLdy(1,x,y)-dLdy(0,x,y))*legendre(r2,2*L(2,x,y)-1)
		+L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*dlegendre(r2,2*L(2,x,y)-1)*(2*dLdy(2,x,y));
	    count++;
	  }
  }
  else 
    throw runtime_error("Invalid indices for triangle");
     
  return dphi;
}

// This function generates heirarchic basis functions for Tetrahedron.
double evalTethPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {

  double phi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for tetrahedron");
  
  if ((i==0)&&(k==0))// vertex
        phi = tetL(j, x, y, z);      
         
  else if(i==1){// edge  
    if ((k>p-2)||(p<2)||(j>5)) throw runtime_error("Invalid for edge shape functions in tetrahedron");
    int j1, j2;
    if(j<=2){
      j1=j;
      j2=(j+1)%3;
    }
    else{
      j1=(j+1)%4;
      j2=3;
    }
    
    if((j==2)||((j==1)&&(type==1))) swap(j1,j2); // Fix orientation. See Ainsworth and Coyle
    
    phi = tetL(j1,x,y,z)*tetL(j2,x,y,z)*kernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z));
    
  }
  else if (i==2){// face 
    if ((k>=(p-1)*(p-2)/2)||(p<3)||(j>3)) throw runtime_error("Invalid for face shape functions for tetrahedron");
    int J[3];

    int indices[2][4][3]={{{0,1,2},{0,1,3},{1,2,3},{0,2,3}},{{0,2,1},{0,1,3},{2,1,3},{0,2,3}}};
    for(int c=0;c<3;c++)
      J[c]=indices[type][j][c];
    
    int count = 0;
    for(int K=3;K<=p;K++)
      for(int r1=K-3;r1>=0;r1--)
	for(int r2=0;r2<=K-3;r2++)
	  if(r1+r2==K-3){
	    if (k==count) 
              phi = tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z)); 
            count++;
      }
    }
           
  else if((i==3)&&(j==0)){ // internal
    if((p<4)||(k>=(p-1)*(p-2)*(p-3)/6.)) throw runtime_error("Invalid internal modes for tetrahedron");
    int count = 0;
    for(int N=4;N<=p;N++)
      for(int r1=0;r1<=N-4;r1++)
        for(int r2=0;r2<=N-4;r2++)
          for(int r3=0;r3<=N-4;r3++){
            if (k==count)
              phi = tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1);
	    if (r1+r2+r3 == N-4)
              count++;
          }
  }
  else 
    throw runtime_error("Invalid indices for tetrahedron: eval_phi");
  return phi;
}

double evalTethDPhidx(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {
  double dphi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for tetrahedron");
  
  if ((i==0)&&(k==0))// vertex
    dphi = dtetLdx(j, x, y, z);      
  
  else if(i==1) {// edge
    if ((k>p-2)||(p<2)||(j>5)) throw runtime_error("Invalid indices for edge shape functions in tetrahedron");
    int j1, j2;
   if(j<=2){
      j1=j;
      j2=(j+1)%3;
    }
    else{
      j1=(j+1)%4;
      j2=3;
    }
   if((j==2)||((j==1)&&(type==1))) swap(j1,j2); // Fix orientation. See Ainsworth and Coyle

  dphi = dtetLdx(j1,x,y,z)*tetL(j2,x,y,z)*kernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z)) + tetL(j1,x,y,z)*dtetLdx(j2,x,y,z)*kernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z))
     + tetL(j1,x,y,z)*tetL(j2,x,y,z)*dkernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z))*(dtetLdx(j2,x,y,z)-dtetLdx(j1,x,y,z));
 }
 
  else if(i==2){// Face 
    if ((k>=(p-1)*(p-2)/2)||(p<3)||(j>3)) throw runtime_error("Invalid for face shape functions for tetrahedron");
    int J[3];
    
    int indices[2][4][3]={{{0,1,2},{0,1,3},{1,2,3},{0,2,3}},{{0,2,1},{0,1,3},{2,1,3},{0,2,3}}};
    for(int c=0;c<3;c++)
      J[c]=indices[type][j][c];

    int count = 0;
    for(int K=3;K<=p;K++)
        for(int r1=K-3;r1>=0;r1--)
          for(int r2=0;r2<=K-3;r2++)
            if(r1+r2==K-3){
              if (k==count) 
               dphi = dtetLdx(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z)) 
	        + tetL(J[0],x,y,z)*dtetLdx(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))
                + tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*dtetLdx(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))
		+ tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*dlegendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))*(dtetLdx(J[1],x,y,z)-dtetLdx(J[0],x,y,z))
                + tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*dlegendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))*(dtetLdx(J[2],x,y,z)-dtetLdx(J[0],x,y,z)); 
              count++;
            }
  }         

  else if ((i==3)&&(j==0)){// internal
    if((p<4)||(k>=(p-1)*(p-2)*(p-3)/6.)) throw runtime_error("Invalid internal modes for tetrahedron");
    int count = 0;
    for(int N=4;N<=p;N++)
      for(int r1=0;r1<=N-4;r1++)
        for(int r2=0;r2<=N-4;r2++)
          for(int r3=0;r3<=N-4;r3++){
            if (k==count)
            dphi = dtetLdx(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*dtetLdx(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*dtetLdx(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*dtetLdx(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*dlegendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*(dtetLdx(1,x,y,z)-dtetLdx(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*dlegendre(r2,2*tetL(2,x,y,z)-1)*2*dtetLdx(2,x,y,z)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*dlegendre(r1,2*tetL(3,x,y,z)-1)*2*dtetLdx(3,x,y,z);
          if (r1+r2+r3 == N-4)
              count++; 
	  }
    }
  
else 
    throw runtime_error("Invalid indices for tetrahedron: eval_dphidx");
  return dphi;
}

double evalTethDPhidy(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {
  double dphi=0.;

  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for tetrahedron");
  
  if ((i==0)&&(k==0))// vertex
      dphi = dtetLdy(j, x, y, z);      
    
  else if((i==1)) {// Edge 
    if ((k>p-2)||(p<2)||(j>5)) throw runtime_error("Invalid indices for edge shape functions in tetrahedron");
    int j1, j2;
    if(j<=2){
      j1=j;
      j2=(j+1)%3;
    }
    else{
      j1=(j+1)%4;
      j2=3;
    }
    
    if((j==2)||((j==1)&&(type==1))) swap(j1,j2); // Fix orientation. See Ainsworth and Coyle

    dphi = dtetLdy(j1,x,y,z)*tetL(j2,x,y,z)*kernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z)) + tetL(j1,x,y,z)*dtetLdy(j2,x,y,z)*kernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z))
      + tetL(j1,x,y,z)*tetL(j2,x,y,z)*dkernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z))*(dtetLdy(j2,x,y,z)-dtetLdy(j1,x,y,z));
  }
  
   else if(i==2){// Face 
    if ((k>=(p-1)*(p-2)/2)||(p<3)||(j>3)) throw runtime_error("Invalid for face shape functions for tetrahedron");
    int J[3];
    
    int indices[2][4][3]={{{0,1,2},{0,1,3},{1,2,3},{0,2,3}},{{0,2,1},{0,1,3},{2,1,3},{0,2,3}}};
    for(int c=0;c<3;c++)
      J[c]=indices[type][j][c];
    
    int count = 0;
    for(int K=3;K<=p;K++)
        for(int r1=K-3;r1>=0;r1--)
          for(int r2=0;r2<=K-3;r2++)
            if(r1+r2==K-3){
              if (k==count) 
                dphi = dtetLdy(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z)) 
	        + tetL(J[0],x,y,z)*dtetLdy(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))
                + tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*dtetLdy(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))
		+ tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*dlegendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))*(dtetLdy(J[1],x,y,z)-dtetLdy(J[0],x,y,z))
                + tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*dlegendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))*(dtetLdy(J[2],x,y,z)-dtetLdy(J[0],x,y,z));  
              count++;
            }
  }         
      
  else if ((i==3)&&(j==0)){// internal 
    if((p<4)||(k>=(p-1)*(p-2)*(p-3)/6.)) throw runtime_error("Invalid internal modes for tetrahedron");
    int count = 0;
    for(int N=4;N<=p;N++)
      for(int r1=0;r1<=N-4;r1++)
        for(int r2=0;r2<=N-4;r2++)
          for(int r3=0;r3<=N-4;r3++){
            if (k==count)
              dphi = dtetLdy(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*dtetLdy(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*dtetLdy(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*dtetLdy(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*dlegendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*(dtetLdy(1,x,y,z)-dtetLdy(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*dlegendre(r2,2*tetL(2,x,y,z)-1)*2*dtetLdy(2,x,y,z)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*dlegendre(r1,2*tetL(3,x,y,z)-1)*2*dtetLdy(3,x,y,z);
           if (r1+r2+r3 == N-4)
              count++; 
	    
	  }
    }
  
else 
    throw runtime_error("Invalid indices for tetrahedron");

  return dphi;
}

double evalTethDPhidz(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {
  double dphi=0.;
 
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for tetrahedron");
  
  if ((i==0)&&(k==0))// vertex
    dphi = dtetLdz(j, x, y, z);      
   
  else if(i==1) {// edge 
    
      if ((k>p-2)||(p<2)||(j>5)) throw runtime_error("Invalid indices for edge shape functions in tetrahedron");
      int j1, j2;
      if(j<=2){
      j1=j;
      j2=(j+1)%3;
      }
      else{
	j1=(j+1)%4;
	j2=3;
      }
      
      if((j==2)||((j==1)&&(type==1))) swap(j1,j2); // Fix orientation. See Ainsworth and Coyle

      dphi = dtetLdz(j1,x,y,z)*tetL(j2,x,y,z)*kernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z)) + tetL(j1,x,y,z)*dtetLdz(j2,x,y,z)*kernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z))
      + tetL(j1,x,y,z)*tetL(j2,x,y,z)*dkernel(k+1,tetL(j2,x,y,z)-tetL(j1,x,y,z))*(dtetLdz(j2,x,y,z)-dtetLdz(j1,x,y,z));
    }
  
  else if(i==2){// Face
    if ((k>=(p-1)*(p-2)/2)||(p<3)||(j>3)) throw runtime_error("Invalid for face shape functions for tetrahedron");
    int J[3];

    int indices[2][4][3]={{{0,1,2},{0,1,3},{1,2,3},{0,2,3}},{{0,2,1},{0,1,3},{2,1,3},{0,2,3}}};
    for(int c=0;c<3;c++)
      J[c]=indices[type][j][c];
    
    int count = 0;
    for(int K=3;K<=p;K++)
        for(int r1=K-3;r1>=0;r1--)
          for(int r2=0;r2<=K-3;r2++)
            if(r1+r2==K-3){
              if (k==count) 
               dphi = dtetLdz(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z)) 
	        + tetL(J[0],x,y,z)*dtetLdz(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))
                + tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*dtetLdz(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))
		+ tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*dlegendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*legendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))*(dtetLdz(J[1],x,y,z)-dtetLdz(J[0],x,y,z))
                + tetL(J[0],x,y,z)*tetL(J[1],x,y,z)*tetL(J[2],x,y,z)*legendre(r1,tetL(J[1],x,y,z)-tetL(J[0],x,y,z))*dlegendre(r2,tetL(J[2],x,y,z)-tetL(J[0],x,y,z))*(dtetLdz(J[2],x,y,z)-dtetLdz(J[0],x,y,z));  
              count++;
            }
  }   
  
  else if((i==3)&&(j==0)){// internal 
  if((p<4)||(k>=(p-1)*(p-2)*(p-3)/6.)) throw runtime_error("Invalid internal modes for tetrahedron");
    int count = 0;
    for(int N=4;N<=p;N++)
      for(int r1=0;r1<=N-4;r1++)
        for(int r2=0;r2<=N-4;r2++)
          for(int r3=0;r3<=N-4;r3++){
            if (k==count)
              dphi = dtetLdz(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*dtetLdz(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*dtetLdz(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*dtetLdz(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*dlegendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*(dtetLdz(1,x,y,z)-dtetLdz(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*dlegendre(r2,2*tetL(2,x,y,z)-1)*2*dtetLdz(2,x,y,z)*legendre(r1,2*tetL(3,x,y,z)-1)
              + tetL(0,x,y,z)*tetL(1,x,y,z)*tetL(2,x,y,z)*tetL(3,x,y,z)*legendre(r3,tetL(1,x,y,z)-tetL(0,x,y,z))*legendre(r2,2*tetL(2,x,y,z)-1)*dlegendre(r1,2*tetL(3,x,y,z)-1)*2*dtetLdz(3,x,y,z);
          if (r1+r2+r3 == N-4)
              count++; 
	    
	  }
    }
   
else 
    throw runtime_error("Invalid indices for tetrahedron");     

  return dphi;
}

double evalWedgehPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {

  double phi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for wedge");
  
  if((i==0)&&(k==0))// vertex 
    phi =(j<3)?L(j, x, y)*(1-z)/2.:L(j-3, x, y)*(1+z)/2.;   
 
  else if(i==1) {// Edge 

    if ((j>8)||(k>p-2)||(p<2))
      throw runtime_error("Invalid indices for edge modes in wedge");

      int j1, j2;
      if(j<3){
        j1=j;
        j2=(1+j)%3;
        phi = L(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1-z)/2.;
      }
      else if(j<6){
        j1=j-3;
        j2=(1+j1)%3;
        phi = L(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1+z)/2.;
      }
      else
        phi = L(j-6,x,y)*lobatto(k+2,z);
 }

  else if((i==2)&&(j>2)){// triangular face 
  
    if((j>4)||(k>=(p-1)*(p-2)/2)||(p<3))
      throw runtime_error("Invalid indices for triangular face modes in wedge");
    
    int J[3];

    int indices[6][3]={{0,1,2},{1,2,0},{2,0,1},{0,2,1},{1,0,2},{2,1,0}};
    for(int c=0;c<3;c++)
      J[c]=indices[type][c];
    
      int count = 0;
      for(int K=3;K<=p;K++)
        for(int r1=K-3;r1>=0;r1--)
          for(int r2=0;r2<=K-3;r2++)
            if (r1+r2==K-3){
              if (count == k){
		phi=(j==3)?(1-z)/2:(1+z)/2;
                phi *= L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y));
	    }
              count++;
            }
   }
   
   else if((i==2)&&(j<=2)){// quadrilateral face

     if ((k>=(p-2)*(p-3)/2)||(p<4))
       throw runtime_error("Invalid indices for quadrilateral face modes in wedge");
     
    int count = 0;
    if(type==0){
      for(int K=4;K<=p;K++)
        for(int r1=K-2;r1>=2;r1--)
          for(int r2=2;r2<=K-2;r2++)
            if (r1+r2==K){
              if (count == k)
                phi = lobatto(r2,z)*L(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y));
              count++;
            }
    }
    else if(type==1){
      for(int K=4;K<=p;K++)
        for(int r2=K-2;r2>=2;r2--)
          for(int r1=2;r1<=K-2;r1++)
            if (r1+r2==K){
              if (count == k)
                phi = lobatto(r2,z)*L(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y));
              count++;
            }
    }
    
  }
    
                   
  else if((i==3)&&(j==0)){// Internal 

    if((p<5)||(k>=(p-2)*(p-3)*(p-4)/6.))
      throw runtime_error("Invalid indices for internal modes in wedge");
    
    int count = 0;
    for(int K=5;K<=p;K++)
      for(int m=K-3;m>=2;m--)
        for(int r1=K-5;r1>=0;r1--)
          for(int r2=0;r2<=K-5;r2++)
            if (r1+r2+m==K-3){
              if (count==k)
                phi = lobatto(m,z)*L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1);
              count++;
            }
  }

  else 
    throw runtime_error("Invalid indices for wedge!");  

  return phi;
}

double evalWedgehDPhiDx(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {

  double dphi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for wedge");
  
  if((i==0)&&(k==0))// vertex 
    dphi=(j<3)?dLdx(j, x, y)*(1-z)/2.:dLdx(j-3, x, y)*(1+z)/2.;   
  
  else if(i==1){// Edge 

    if ((j>8)||(k>p-2)||(p<2))
      throw runtime_error("Invalid indices for edge modes in wedge");
    
    int j1, j2;
    if(j<3){
        j1=j;
        j2=(1+j)%3;
        dphi = dLdx(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1-z)/2. 
        + L(j1,x,y)*dLdx(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1-z)/2.
        + L(j1,x,y)*L(j2,x,y)*dkernel(k+1,L(j2,x,y)-L(j1,x,y))*(dLdx(j2,x,y)-dLdx(j1,x,y))*(1-z)/2.;
      }
      else if(j<6){
        j1=j-3;
        j2=(1+j1)%3;
        dphi = dLdx(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1+z)/2. 
        + L(j1,x,y)*dLdx(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1+z)/2.
        + L(j1,x,y)*L(j2,x,y)*dkernel(k+1,L(j2,x,y)-L(j1,x,y))*(dLdx(j2,x,y)-dLdx(j1,x,y))*(1+z)/2.;
      }
      else
        dphi = dLdx(j-6,x,y)*lobatto(k+2,z);
  }
  
  else if((i==2)&&(j>2)){//triangular face
  
    if((j>4)||(k>=(p-1)*(p-2)/2)||(p<3))
      throw runtime_error("Invalid indices for triangular face modes in wedge");
        
    int J[3];

    int indices[6][3]={{0,1,2},{1,2,0},{2,0,1},{0,2,1},{1,0,2},{2,1,0}};
    for(int c=0;c<3;c++)
      J[c]=indices[type][c];
    
    int count = 0;
    for(int K=3;K<=p;K++)
      for(int r1=K-3;r1>=0;r1--)
        for(int r2=0;r2<=K-3;r2++)
          if (r1+r2==K-3){
            if (count == k){
	      dphi=(j==3)?(1-z)/2:(1+z)/2;
              dphi *= dLdx(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y)) 
	      + L(0,x,y)*dLdx(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y))
              + L(0,x,y)*L(1,x,y)*dLdx(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y)) 
	      + L(0,x,y)*L(1,x,y)*L(2,x,y)*dlegendre(r1,L(J[1],x,y)-L(J[0],x,y))*(dLdx(J[1],x,y)-dLdx(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y))
              + L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*dlegendre(r2,L(J[2],x,y)-L(J[0],x,y))*(dLdx(J[2],x,y)-dLdx(J[0],x,y));
	    }
	    count++;
          }
    }
   
  else if((i==2)&&(j<=2)){ //quadrilateral face

    if ((k>=(p-2)*(p-3)/2)||(p<4))
       throw runtime_error("Invalid indices for quadrilateral face modes in wedge");

    int count = 0;
    if(type==0){
    for(int K=4;K<=p;K++)
      for(int r1=K-2;r1>=2;r1--)
        for(int r2=2;r2<=K-2;r2++)
          if (r1+r2==K){
            if (count == k)
              dphi = lobatto(r2,z)*(dLdx(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y)) 
		+ L(j,x,y)*dLdx((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))
		+ L(j,x,y)*L((1+j)%3,x,y)*dkernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))*(dLdx((1+j)%3,x,y)-dLdx(j,x,y)));
            count++;
          }
    }
    else if (type==1){
      for(int K=4;K<=p;K++)
      for(int r2=K-2;r2>=2;r2--)
        for(int r1=2;r1<=K-2;r1++)
          if (r1+r2==K){
            if (count == k)
              dphi = lobatto(r2,z)*(dLdx(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y)) 
		+ L(j,x,y)*dLdx((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))
		+ L(j,x,y)*L((1+j)%3,x,y)*dkernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))*(dLdx((1+j)%3,x,y)-dLdx(j,x,y)));
            count++;
          }
    }
    }
  
  
  else if((i==3)&&(j==0)){// Internal 

    if((p<5)||(k>=(p-2)*(p-3)*(p-4)/6.))
      throw runtime_error("Invalid indices for internal modes in wedge");
    
    int count = 0;
    for(int K=5;K<=p;K++)
      for(int m=K-3;m>=2;m--)
        for(int r1=K-5;r1>=0;r1--)
          for(int r2=0;r2<=K-5;r2++)
            if (r1+r2+m==K-3){
              if (count == k)
                dphi = lobatto(m,z)*(dLdx(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1) 
		+ L(0,x,y)*dLdx(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)
                + L(0,x,y)*L(1,x,y)*dLdx(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1) 
		+ L(0,x,y)*L(1,x,y)*L(2,x,y)*dlegendre(r1,L(1,x,y)-L(0,x,y))*(dLdx(1,x,y)-dLdx(0,x,y))*legendre(r2,2*L(2,x,y)-1)
                + L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*dlegendre(r2,2*L(2,x,y)-1)*2*dLdx(2,x,y));
              count++;
            }
  }

  else 
    throw runtime_error("Invalid indices for wedge!");  
    
  return dphi;
}

double evalWedgehDPhiDy(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {
  
  double dphi=0.;
  
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for wedge");
  
  if((i==0)&&(k==0))// vertex 
     dphi=(j<3)?dLdy(j, x, y)*(1-z)/2.:dLdy(j-3, x, y)*(1+z)/2.; 

  else if(i==1){// Edge 

    if ((j>8)||(k>p-2)||(p<2))
      throw runtime_error("Invalid indices for edge modes in wedge");
    
      int j1, j2;
      if(j<3){
        j1=j;
        j2=(1+j)%3;
        dphi = dLdy(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1-z)/2. 
        + L(j1,x,y)*dLdy(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1-z)/2.
        + L(j1,x,y)*L(j2,x,y)*dkernel(k+1,L(j2,x,y)-L(j1,x,y))*(dLdy(j2,x,y)-dLdy(j1,x,y))*(1-z)/2.;
      }
      else if(j<6){
        j1=j-3;
        j2=(1+j1)%3;
        dphi = dLdy(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1+z)/2. 
        + L(j1,x,y)*dLdy(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))*(1+z)/2.
        + L(j1,x,y)*L(j2,x,y)*dkernel(k+1,L(j2,x,y)-L(j1,x,y))*(dLdy(j2,x,y)-dLdy(j1,x,y))*(1+z)/2.;
      }
      else
        dphi = dLdy(j-6,x,y)*lobatto(k+2,z);
  }
  
  else if((i==2)&&(j>2)){//triangular face
  
    if((j>4)||(k>=(p-1)*(p-2)/2)||(p<3))
      throw runtime_error("Invalid indices for triangular face modes in wedge");
    
    int indices[6][3]={{0,1,2},{1,2,0},{2,0,1},{0,2,1},{1,0,2},{2,1,0}};
    int J[3];
    for(int c=0;c<3;c++)
      J[c]=indices[type][c];
    
    int count = 0;
    for(int K=3;K<=p;K++)
      for(int r1=K-3;r1>=0;r1--)
        for(int r2=0;r2<=K-3;r2++)
          if (r1+r2==K-3){
            if(count == k){
	      dphi=(j==3)?(1-z)/2:(1+z)/2;
              dphi *= dLdy(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y)) 
	      + L(0,x,y)*dLdy(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y))
              + L(0,x,y)*L(1,x,y)*dLdy(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y)) 
	      + L(0,x,y)*L(1,x,y)*L(2,x,y)*dlegendre(r1,L(J[1],x,y)-L(J[0],x,y))*(dLdy(J[1],x,y)-dLdy(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y))
              + L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*dlegendre(r2,L(J[2],x,y)-L(J[0],x,y))*(dLdy(J[2],x,y)-dLdy(J[0],x,y));
	    }
            count++;
          }
    }
   
  else if((i==2)&&(j<=2)){ //quadrilateral face

    if ((k>=(p-2)*(p-3)/2)||(p<4))
       throw runtime_error("Invalid indices for quadrilateral face modes in wedge");

    int count = 0;
    if(type==0){
    for(int K=4;K<=p;K++)
      for(int r1=K-2;r1>=2;r1--)
        for(int r2=2;r2<=K-2;r2++)
          if (r1+r2==K){
            if (count == k)
              dphi = lobatto(r2,z)*(dLdy(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))
		+ L(j,x,y)*dLdy((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))
		+ L(j,x,y)*L((1+j)%3,x,y)*dkernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))*(dLdy((1+j)%3,x,y)-dLdy(j,x,y)));
            count++;
          }
    }
    else if(type==1){
      for(int K=4;K<=p;K++)
      for(int r2=K-2;r2>=2;r2--)
        for(int r1=2;r1<=K-2;r1++)
          if (r1+r2==K){
            if (count == k)
              dphi = lobatto(r2,z)*(dLdy(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))
		+ L(j,x,y)*dLdy((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))
		+ L(j,x,y)*L((1+j)%3,x,y)*dkernel(r1-1,L((1+j)%3,x,y)-L(j,x,y))*(dLdy((1+j)%3,x,y)-dLdy(j,x,y)));
            count++;
          }
    }
    }

  else if((i==3)&&(j==0)){// Internal 
   
    if((p<5)||(k>=(p-2)*(p-3)*(p-4)/6.))
      throw runtime_error("Invalid indices for internal modes in wedge");

    int count = 0;
    for(int K=5;K<=p;K++)
      for(int m=K-3;m>=2;m--)
        for(int r1=K-5;r1>=0;r1--)
          for(int r2=0;r2<=K-5;r2++)
            if (r1+r2+m==K-3){
              if (count == k)
                dphi = lobatto(m,z)*(dLdy(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1) 
		+ L(0,x,y)*dLdy(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1)
                + L(0,x,y)*L(1,x,y)*dLdy(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1) 
		+ L(0,x,y)*L(1,x,y)*L(2,x,y)*dlegendre(r1,L(1,x,y)-L(0,x,y))*(dLdy(1,x,y)-dLdy(0,x,y))*legendre(r2,2*L(2,x,y)-1)
                + L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*dlegendre(r2,2*L(2,x,y)-1)*2*dLdy(2,x,y));
              count++;
            }
  }

  else 
    throw runtime_error("Invalid indices for wedge!"); 
  
  return dphi;
}

double evalWedgehDPhiDz(const int i, const int j, const int k, const double x, const double y, const double z, const int p, const int type) {
  
  double dphi=0.;

  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>3) throw runtime_error("Invalid type of shape functions for wedge");
  
  if((i==0)&&(k==0))// vertex 
    dphi =(j<3)?-L(j, x, y)/2.:L(j-3, x, y)/2.;   

  else if(i==1){// Edge 

    if ((j>8)||(k>p-2)||(p<2))
      throw runtime_error("Invalid indices for edge modes in wedge");
    
      int j1, j2;
      if(j<3){
        j1=j;
        j2=(1+j)%3;
        dphi = -L(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))/2.;
      }
      else if(j<6){
        j1=j-3;
        j2=(1+j1)%3;
        dphi = L(j1,x,y)*L(j2,x,y)*kernel(k+1,L(j2,x,y)-L(j1,x,y))/2.;
      }
      else
        dphi = L(j-6,x,y)*dlobatto(k+2,z);
  }


  else if((i==2)&&(j>2)){//triangular face
  
    if((j>4)||(k>=(p-1)*(p-2)/2)||(p<3))
      throw runtime_error("Invalid indices for triangular face modes in wedge");
    
    int indices[6][3]={{0,1,2},{1,2,0},{2,0,1},{0,2,1},{1,0,2},{2,1,0}};
    int J[3];
    for(int c=0;c<3;c++)
      J[c]=indices[type][c];
    
    int count = 0;
    for(int K=3;K<=p;K++)
      for(int r1=K-3;r1>=0;r1--)
        for(int r2=0;r2<=K-3;r2++)
          if (r1+r2==K-3){
            if (count == k){
	      dphi=(j==3)?-0.5:0.5;
              dphi *= L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(J[1],x,y)-L(J[0],x,y))*legendre(r2,L(J[2],x,y)-L(J[0],x,y));
	    }
            count++;
          }
    }
   
    else if((i==2)&&(j<=2)){ //quadrilateral face

    if ((k>=(p-2)*(p-3)/2)||(p<4))
       throw runtime_error("Invalid indices for quadrilateral face modes in wedge");
    
    int count = 0;
    if(type==0){
    for(int K=4;K<=p;K++)
      for(int r1=K-2;r1>=2;r1--)
        for(int r2=2;r2<=K-2;r2++)
          if (r1+r2==K){
            if (count == k)
              dphi = dlobatto(r2,z)*L(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y));
            count++;
          }
    }
    else if(type==1)
      {
    for(int K=4;K<=p;K++)
      for(int r2=K-2;r2>=2;r2--)
        for(int r1=2;r1<=K-2;r1++)
          if (r1+r2==K){
            if (count == k)
              dphi = dlobatto(r2,z)*L(j,x,y)*L((1+j)%3,x,y)*kernel(r1-1,L((1+j)%3,x,y)-L(j,x,y));
            count++;
          }
    }
    }
   
  else if((i==3)&&(j==0)){// Internal 
   
    if((p<5)||(k>=(p-2)*(p-3)*(p-4)/6.))
      throw runtime_error("Invalid indices for internal modes in wedge");
    
    int count = 0;
     for(int K=5;K<=p;K++)
      for(int m=K-3;m>=2;m--)
        for(int r1=K-5;r1>=0;r1--)
          for(int r2=0;r2<=K-5;r2++)
            if (r1+r2+m==K-3){
              if (count == k)
                dphi = dlobatto(m,z)*L(0,x,y)*L(1,x,y)*L(2,x,y)*legendre(r1,L(1,x,y)-L(0,x,y))*legendre(r2,2*L(2,x,y)-1);
              count++;
            }
  }    

  else 
    throw runtime_error("Invalid indices for wedge!"); 

  return dphi;
}
