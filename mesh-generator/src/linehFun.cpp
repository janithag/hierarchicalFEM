 // This file contains the functions needed for heirarchic basis functions for line, Quadrilateral, and Hexahedron.
// i: type (vertex:0, edge:1, face:2, internal:3), j: location, k: shape function number (we use 0-based index)

#include <stdexcept>
#include <iostream>
#include <cmath>
#include "hFun.hpp"

using std::runtime_error;
using std::cout;
using std::swap;
// This function generates heirarchic basis for line.
double evalLinehPhi(const int i, const int j, const int k, const double x, const int p) {

  double phi=0.;
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>1) throw runtime_error("Invalid type of shape functions for line");

  if ((i==0)&&(k==0)){ // vertex
    if (j==0)
      phi=(1-x)/2;
    else if (j==1)
      phi=(1+x)/2;
    else throw runtime_error("Invalid location for vertex shape functions in line");
  }
  else if((i==1)&&(j==0)){ // edge
    if (k<p-1) phi=lobatto(k+2,x);
    else throw runtime_error("Invalid mode for edge shape functions in line");
  }
  else 
    throw runtime_error("Invalid indices for line");
 
  return phi;
}

// This function computes the derivative of the basis functions for line.
double evalLinehDPhi(const int i,  const int j, const int k, const double x, const int p){
  double dphi=1.;
  if(p<1) throw runtime_error("Number of element segments < 1");
  if(i>1) throw runtime_error("Invalid type of shape functions for line");

  if ((i==0)&&(k==0)){ //vertex
    if (j==0)
      dphi=-0.5;
    else if (j==1)
      dphi=0.5;
    else throw runtime_error("Invalid location for vertex shape functions in line");
  }
  else if((i==1)&&(j==0)){ //edge
    if(k<p-1)   dphi=dlobatto(k+2,x);
    else throw runtime_error("Invalid mode for edge shape functions in line");
  }
  else 
    throw runtime_error("Invalid indices for line");
  return dphi;
}

// This function generates heirarchic basis for quadrilateral.
double evalQuadhPhi(const int i, const int j, const int k, const double x, const double y, const int n) {

  double phi=0.;
  
  if(n<1) throw runtime_error("number of element segments < 1 ");
  
  if ((i==0)&&(k==0))// Vertex 
       switch(j){
         case 0: phi = (1-x)*(1-y)/4.; break;
         case 1: phi = (1+x)*(1-y)/4.; break;
         case 2: phi = (1+x)*(1+y)/4.; break;
         case 3: phi = (1-x)*(1+y)/4.; break;
         default: throw runtime_error("Invalid vertex number for quadrilateral.");
       }
         
  else if ((i==1)&&(k<n-1)&&(n>1)){ // Edge 
    switch(j){ 
    case 0: phi = (1-y)*lobatto(k+2,x)/2.; break;
    case 1: phi = (1+x)*lobatto(k+2,y)/2.; break;
    case 2: phi = (1+y)*lobatto(k+2,x)/2.; break;
    case 3: phi = (1-x)*lobatto(k+2,y)/2.; break;
    default: throw runtime_error("Invalid edge number for quadrilateral. ");
    }
  }     
       
  else if ((i==2)&&(j==0)&&(n>=4)&&(k<(n-2)*(n-3)*0.5)){ // Face modes
    int count = 0;
    for(int N=4;N<=n;N++)
      for(int l=N-2;l>=2;l--)
	for(int m=2;m<=N-2;m++)
	  if (l+m==N){
	    if (count == k)
	      phi = lobatto(m,x)*lobatto(l,y);
	    count++;
	  }
  }

  else throw runtime_error("Invalid indices for quadrilateral in \"linehFun\"");
  return phi;
}
 
double evalQuadhDPhiDx(const int i, const int j, const int k, const double x, const double y, const int n) {

  double dphidx=0.;
  
  if(n<1) throw runtime_error("number of element segments < 1 ");
   
  if ((i==0)&&(k==0))// vertex 
    switch(j){
    case 0: dphidx = -(1-y)/4.; break;
    case 1: dphidx = (1-y)/4.; break;
    case 2: dphidx = (1+y)/4.; break;
    case 3: dphidx = -(1+y)/4.; break;
    default: throw runtime_error("Invalid vertex number for quadrilateral.");
    }

  else if ((i==1)&&(k<n-1)&&(n>1)) // edge
    switch(j){ 
    case 0: dphidx = (1-y)*dlobatto(k+2,x)/2.;break;
    case 1: dphidx = lobatto(k+2,y)/2.;break;
    case 2: dphidx = (1+y)*dlobatto(k+2,x)/2.;break;
    case 3: dphidx = -lobatto(k+2,y)/2.;break;
    default: throw runtime_error("Invalid edge number for quadrilateral. ");
    }
              
  else if ((i==2)&&(j==0)&&(n>=4)&&(k<(n-2)*(n-3)*0.5)){ // face
    int count = 0;
    for(int N=4;N<=n;N++)
      for(int l=N-2;l>=2;l--)
        for(int m=2;m<=N-2;m++)
          if (l+m==N){
            if (count == k)
              dphidx = dlobatto(m,x)*lobatto(l,y);
            count++;
          }
    }  

  else throw runtime_error("Invalid indices for quadrilateral");
  return dphidx;
}

double evalQuadhDPhiDy(const int i, const int j, const int k, const double x, const double y, const int n) {

  double dphidy = 0.;
  
  if(n<1) throw runtime_error("number of element segments < 1 ");
  
  if ((i==0)&&(k==0)) // vertex
    switch(j){
    case 0: dphidy = -(1-x)/4.;break;
    case 1: dphidy = -(1+x)/4.;break;
    case 2: dphidy = (1+x)/4.;break;
    case 3: dphidy = (1-x)/4.;break;
    default: throw runtime_error("Invalid vertex number for quadrilateral.");
    }

  else if ((i==1)&&(k<n-1)&&(n>1)) // edge
    switch(j){ 
    case 0: dphidy = -lobatto(k+2,x)/2.;break;
    case 1: dphidy = (1+x)*dlobatto(k+2,y)/2.;break;
    case 2: dphidy = lobatto(k+2,x)/2.;break;
    case 3: dphidy = (1-x)*dlobatto(k+2,y)/2.;break;
    default: throw runtime_error("Invalid edge number for quadrilateral. ");
    }

  else if ((i==2)&&(j==0)&&(n>=4)&&(k<(n-2)*(n-3)*0.5)){ // face  
    int count = 0;
    for(int N=4;N<=n;N++)
      for(int l=N-2;l>=2;l--)
        for(int m=2;m<=N-2;m++)
          if (l+m==N){
            if (count == k)
              dphidy = lobatto(m,x)*dlobatto(l,y);
            count++;
          }
  }         
  
  else throw runtime_error("Invalid indices for quadrilateral");
  return dphidy;
}

// This function generates heirarchic basis for hexahedron.
double evalHexhPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type) {

  double phi=0.;
  
  if(n<1) throw runtime_error("number of element segments < 1 ");
   
  if((i==0)&&(k==0))// Vertex
    switch(j){
    case 0: phi = (1-x)*(1-y)*(1-z)/8.; break;
    case 1: phi = (1+x)*(1-y)*(1-z)/8.; break;
    case 2: phi = (1+x)*(1+y)*(1-z)/8.; break;
    case 3: phi = (1-x)*(1+y)*(1-z)/8.; break;
    case 4: phi = (1-x)*(1-y)*(1+z)/8.; break;
    case 5: phi = (1+x)*(1-y)*(1+z)/8.; break;
    case 6: phi = (1+x)*(1+y)*(1+z)/8.; break;
    case 7: phi = (1-x)*(1+y)*(1+z)/8.; break;
    default: throw runtime_error("Invalid vertex number for hexahedron.");
    }
       
  else if((i==1)&&(k<n-1)&&(n>1))// Edge
    switch(j){
    case 0: phi = (1-y)*(1-z)*lobatto(k+2,x)/4.;break;
    case 1: phi = (1+x)*(1-z)*lobatto(k+2,y)/4.;break;
    case 2: phi = (1+y)*(1-z)*lobatto(k+2,x)/4.;break;
    case 3: phi = (1-x)*(1-z)*lobatto(k+2,y)/4.;break;
    case 4: phi = (1-y)*(1+z)*lobatto(k+2,x)/4.;break;
    case 5: phi = (1+x)*(1+z)*lobatto(k+2,y)/4.;break;
    case 6: phi = (1+y)*(1+z)*lobatto(k+2,x)/4.;break;
    case 7: phi = (1-x)*(1+z)*lobatto(k+2,y)/4.;break;
    case 8: phi = (1-x)*(1-y)*lobatto(k+2,z)/4.;break;
    case 9: phi = (1+x)*(1-y)*lobatto(k+2,z)/4.;break;
    case 10: phi = (1+x)*(1+y)*lobatto(k+2,z)/4.;break;
    case 11: phi = (1-x)*(1+y)*lobatto(k+2,z)/4.;break;
    default: throw runtime_error("Invalid edge number for hexahedron. ");
    }   
    
  else if ((i==2)&&(n>=4)&&(k<(n-2)*(n-3)*0.5)){ // Face 
    int count = 0;
    if(type==1){
    for(int N=4;N<=n;N++)
      for(int r2=N-2;r2>=2;r2--)
        for(int r1=2;r1<=N-2;r1++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: phi = (1-y)*lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 1: phi = (1+x)*lobatto(r1,y)*lobatto(r2,z)/2.; break;
	      case 2: phi = (1+y)*lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 3: phi = (1-x)*lobatto(r1,y)*lobatto(r2,z)/2.; break;
	      case 4: phi = (1-z)*lobatto(r1,x)*lobatto(r2,y)/2.; break;
	      case 5: phi = (1+z)*lobatto(r1,x)*lobatto(r2,y)/2.; break;
	      default: throw runtime_error("Invalid face number for hexahedron. ");
              }
            count++;
	  }
	} 
    }
   else if(type==0){
    for(int N=4;N<=n;N++)
      for(int r1=N-2;r1>=2;r1--)
        for(int r2=2;r2<=N-2;r2++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: phi = (1-y)*lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 1: phi = (1+x)*lobatto(r1,y)*lobatto(r2,z)/2.; break;
	      case 2: phi = (1+y)*lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 3: phi = (1-x)*lobatto(r1,y)*lobatto(r2,z)/2.; break;
	      case 4: phi = (1-z)*lobatto(r1,x)*lobatto(r2,y)/2.; break;
	      case 5: phi = (1+z)*lobatto(r1,x)*lobatto(r2,y)/2.; break;
	      default: throw runtime_error("Invalid face number for hexahedron. ");
              }
            count++;
	  }
	} 
   }
  }

  else if((i==3)&&(j==0)&&(n>=6)&&(k<(n-3)*(n-4)*(n-5)/6)){ // Internal 
    int count = 0;
    for(int N=6;N<=n;N++)
      for(int r=2;r<=n-4;r++)
	for(int q=2;q<=n-4;q++)
	  for(int p=2;p<=n-4;p++){
	    if (p+q+r == N){
	      if (k==count) phi = lobatto(p,x)*lobatto(q,y)*lobatto(r,z);
	      count++;
	    }
	  }
  }

  else throw runtime_error("Invalid indices for hexahedron");
  return phi;
}

double evalHexhDPhiDx(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type) {

  double dphi=0.;
  
  if(n<1) throw runtime_error("number of element segments < 1 ");
    
  if((i==0)&&(k==0)) // vertex
    switch(j){
    case 0: dphi = -(1-y)*(1-z)/8.;break;
    case 1: dphi = (1-y)*(1-z)/8.;break;
    case 2: dphi = (1+y)*(1-z)/8.;break;
    case 3: dphi = -(1+y)*(1-z)/8.;break;
    case 4: dphi = -(1-y)*(1+z)/8.;break;
    case 5: dphi = (1-y)*(1+z)/8.;break;
    case 6: dphi = (1+y)*(1+z)/8.;break;
    case 7: dphi = -(1+y)*(1+z)/8.;break;
    default: throw runtime_error("Invalid vertex number for hexahedron.");
    }

  else if((i==1)&&(k<n-1)&&(n>1))// Edge
    switch(j){
    case 0: dphi = (1-y)*(1-z)*dlobatto(k+2,x)/4.;break;
    case 1: dphi = (1-z)*lobatto(k+2,y)/4.;break;
    case 2: dphi = (1+y)*(1-z)*dlobatto(k+2,x)/4.;break;
    case 3: dphi = -(1-z)*lobatto(k+2,y)/4.;break;
    case 4: dphi = (1-y)*(1+z)*dlobatto(k+2,x)/4.;break;
    case 5: dphi = (1+z)*lobatto(k+2,y)/4.;break;
    case 6: dphi = (1+y)*(1+z)*dlobatto(k+2,x)/4.;break;
    case 7: dphi = -(1+z)*lobatto(k+2,y)/4.;break;
    case 8: dphi = -(1-y)*lobatto(k+2,z)/4.;break;  
    case 9: dphi = (1-y)*lobatto(k+2,z)/4.;break;
    case 10: dphi = (1+y)*lobatto(k+2,z)/4.;break;
    case 11: dphi = -(1+y)*lobatto(k+2,z)/4.;break;
    default: throw runtime_error("Invalid edge number for hexahedron. ");
    }

  else if((i==2)&&(n>=4)&&(k<(n-2)*(n-3)*0.5)){// Face 
    int count = 0;
    if(type==1){
     for(int N=4;N<=n;N++)
      for(int r2=N-2;r2>=2;r2--)
        for(int r1=2;r1<=N-2;r1++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: dphi = (1-y)*dlobatto(r1,x)*lobatto(r2,z)/2.;break;
	      case 1: dphi = lobatto(r1,y)*lobatto(r2,z)/2.;break; 
	      case 2: dphi = (1+y)*dlobatto(r1,x)*lobatto(r2,z)/2.;break;
	      case 3: dphi = -lobatto(r1,y)*lobatto(r2,z)/2.;break;
	      case 4: dphi = (1-z)*dlobatto(r1,x)*lobatto(r2,y)/2.;break;
	      case 5: dphi = (1+z)*dlobatto(r1,x)*lobatto(r2,y)/2.;break;
	      default: throw runtime_error("Invalid face number for hexahedron. ");
              }
	    count++;
	  }
	}      
    }
  else if(type==0){
    for(int N=4;N<=n;N++)
      for(int r1=N-2;r1>=2;r1--)
        for(int r2=2;r2<=N-2;r2++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: dphi = (1-y)*dlobatto(r1,x)*lobatto(r2,z)/2.;break;
	      case 1: dphi = lobatto(r1,y)*lobatto(r2,z)/2.;break; 
	      case 2: dphi = (1+y)*dlobatto(r1,x)*lobatto(r2,z)/2.;break;
	      case 3: dphi = -lobatto(r1,y)*lobatto(r2,z)/2.;break;
	      case 4: dphi = (1-z)*dlobatto(r1,x)*lobatto(r2,y)/2.;break;
	      case 5: dphi = (1+z)*dlobatto(r1,x)*lobatto(r2,y)/2.;break;
	      default: throw runtime_error("Invalid face number for hexahedron. ");
              }
	    count++;
	  }
	}   
  }
  }

  else if((i==3)&&(j==0)&&(n>=6)&&(k<(n-3)*(n-4)*(n-5)/6)){ // Internal
    int count = 0;
    for(int N=6;N<=n;N++)
      for(int r=2;r<=n-4;r++)
        for(int q=2;q<=n-4;q++)
          for(int p=2;p<=n-4;p++){
            if (p+q+r == N){
	      if (k==count) dphi = dlobatto(p,x)*lobatto(q,y)*lobatto(r,z);
	      count++;
	    }
	  }
  }

  else throw runtime_error("Invalid indices for hexahedron");
  return dphi;
}

double evalHexhDPhiDy(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type) {

  double dphi=0.;
  
  if(n<1) throw runtime_error("number of element segments < 1 ");
  
  if((i==0)&&(k==0))// vertex
    switch(j){
    case 0: dphi = -(1-x)*(1-z)/8.;break;
    case 1: dphi = -(1+x)*(1-z)/8.;break;
    case 2: dphi = (1+x)*(1-z)/8.;break;
    case 3: dphi = (1-x)*(1-z)/8.;break;
    case 4: dphi = -(1-x)*(1+z)/8.;break;
    case 5: dphi = -(1+x)*(1+z)/8.;break;
    case 6: dphi = (1+x)*(1+z)/8.;break;
    case 7: dphi = (1-x)*(1+z)/8.;break;
    default: throw runtime_error("Invalid vertex number for hexahedron.");
    }

  else if ((i==1)&&(k<n-1)&&(n>1)) //edge
    switch(j){
    case 0: dphi = -(1-z)*lobatto(k+2,x)/4.;break;
    case 1: dphi = (1+x)*(1-z)*dlobatto(k+2,y)/4.;break;
    case 2: dphi = (1-z)*lobatto(k+2,x)/4.;break;
    case 3: dphi = (1-x)*(1-z)*dlobatto(k+2,y)/4.; break;
    case 4: dphi = -(1+z)*lobatto(k+2,x)/4.;break;
    case 5: dphi = (1+x)*(1+z)*dlobatto(k+2,y)/4.;break;
    case 6: dphi = (1+z)*lobatto(k+2,x)/4.;break;
    case 7: dphi = (1-x)*(1+z)*dlobatto(k+2,y)/4.;break;
    case 8: dphi = -(1-x)*lobatto(k+2,z)/4.; break;
    case 9: dphi = -(1+x)*lobatto(k+2,z)/4.;break;
    case 10: dphi = (1+x)*lobatto(k+2,z)/4.;break;
    case 11: dphi = (1-x)*lobatto(k+2,z)/4.;break;
    default: throw runtime_error("Invalid edge number for hexahedron. ");
    }  
   
  else if((i==2)&&(n>=4)&&(k<(n-2)*(n-3)*0.5)){ // Face
    int count = 0;
    if(type==1){
    for(int N=4;N<=n;N++)
      for(int r2=N-2;r2>=2;r2--)
        for(int r1=2;r1<=N-2;r1++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: dphi = -lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 1: dphi = (1+x)*dlobatto(r1,y)*lobatto(r2,z)/2.;break;
	      case 2: dphi =  lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 3: dphi = (1-x)*dlobatto(r1,y)*lobatto(r2,z)/2.;break;
	      case 4: dphi = (1-z)*lobatto(r1,x)*dlobatto(r2,y)/2.;break;
	      case 5: dphi = (1+z)*lobatto(r1,x)*dlobatto(r2,y)/2.;break;
	      default: throw runtime_error("Invalid face index for hexahedron. ");
              }
	    count++;
	  }
	}
    }
    else if(type==0){
    for(int N=4;N<=n;N++)
      for(int r1=N-2;r1>=2;r1--)
        for(int r2=2;r2<=N-2;r2++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: dphi = -lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 1: dphi = (1+x)*dlobatto(r1,y)*lobatto(r2,z)/2.;break;
	      case 2: dphi =  lobatto(r1,x)*lobatto(r2,z)/2.; break;
	      case 3: dphi = (1-x)*dlobatto(r1,y)*lobatto(r2,z)/2.;break;
	      case 4: dphi = (1-z)*lobatto(r1,x)*dlobatto(r2,y)/2.;break;
	      case 5: dphi = (1+z)*lobatto(r1,x)*dlobatto(r2,y)/2.;break;
	      default: throw runtime_error("Invalid face index for hexahedron. ");
              }
	    count++;
	  }
	}
    }
  }         

  else if((i==3)&&(j==0)&&(n>=6)&&(k<(n-3)*(n-4)*(n-5)/6)){// Internal
    int count = 0;
    for(int N=6;N<=n;N++)
      for(int r=2;r<=n-4;r++)
        for(int q=2;q<=n-4;q++)
          for(int p=2;p<=n-4;p++){
            if (p+q+r == N){
	      if (k==count) dphi = lobatto(p,x)*dlobatto(q,y)*lobatto(r,z);
	      count++;
	    }
          }
  }

  else throw runtime_error("Invalid indices for hexahedron");
  return dphi;
}


double evalHexhDPhiDz(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type) {

  double dphi=0.;
  
  if(n<1) throw runtime_error("number of element segments < 1 ");
  
  if((i==0)&&(k==0))// vertex
       switch(j){
         case 0: dphi = -(1-x)*(1-y)/8.;break;
         case 1: dphi = -(1+x)*(1-y)/8.;break;
         case 2: dphi = -(1+x)*(1+y)/8.;break;
         case 3: dphi = -(1-x)*(1+y)/8.;break;
         case 4: dphi = (1-x)*(1-y)/8.; break;
         case 5: dphi = (1+x)*(1-y)/8.; break;
         case 6: dphi = (1+x)*(1+y)/8.; break;
         case 7: dphi = (1-x)*(1+y)/8.; break;
         default: throw runtime_error("Invalid vertex number for hexahedron.");
       }

  else if((i==1)&&(k<n-1)&&(n>1))// Edge
    switch(j){
    case 0: dphi = -(1-y)*lobatto(k+2,x)/4.;break;
    case 1: dphi = -(1+x)*lobatto(k+2,y)/4.;break;
    case 2: dphi = -(1+y)*lobatto(k+2,x)/4.;break;
    case 3: dphi = -(1-x)*lobatto(k+2,y)/4.;break;
    case 4: dphi = (1-y)*lobatto(k+2,x)/4.;break;
    case 5: dphi = (1+x)*lobatto(k+2,y)/4.;break;
    case 6: dphi = (1+y)*lobatto(k+2,x)/4.;break;
    case 7: dphi = (1-x)*lobatto(k+2,y)/4.;break;
    case 8: dphi = (1-x)*(1-y)*dlobatto(k+2,z)/4.;break;  
    case 9: dphi = (1+x)*(1-y)*dlobatto(k+2,z)/4.;break;
    case 10: dphi = (1+x)*(1+y)*dlobatto(k+2,z)/4.;break;
    case 11: dphi = (1-x)*(1+y)*dlobatto(k+2,z)/4.;break;
    default: throw runtime_error("Invalid edge index for hexahedron. ");
    }

  else if((i==2)&&(n>=4)&&(k<(n-2)*(n-3)*0.5)){// Face
    int count = 0;
    if(type==1){
    for(int N=4;N<=n;N++)
      for(int r2=N-2;r2>=2;r2--)
        for(int r1=2;r1<=N-2;r1++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: dphi = (1-y)*lobatto(r1,x)*dlobatto(r2,z)/2.;break;
	      case 1: dphi = (1+x)*lobatto(r1,y)*dlobatto(r2,z)/2.;break;
	      case 2: dphi = (1+y)*lobatto(r1,x)*dlobatto(r2,z)/2.;break;
	      case 3: dphi = (1-x)*lobatto(r1,y)*dlobatto(r2,z)/2.;break; 
	      case 4: dphi = -lobatto(r1,x)*lobatto(r2,y)/2.;break;
	      case 5: dphi = lobatto(r1,x)*lobatto(r2,y)/2.;break;
	      default: throw runtime_error("Invalid face index for hexahedron. ");
              }
	    count++;
	  }
	}
    }
    else if(type==0){
      for(int N=4;N<=n;N++)
      for(int r1=N-2;r1>=2;r1--)
        for(int r2=2;r2<=N-2;r2++){
	  if (r1+r2 == N){
            if (k==count) 
              switch(j){ 
	      case 0: dphi = (1-y)*lobatto(r1,x)*dlobatto(r2,z)/2.;break;
	      case 1: dphi = (1+x)*lobatto(r1,y)*dlobatto(r2,z)/2.;break;
	      case 2: dphi = (1+y)*lobatto(r1,x)*dlobatto(r2,z)/2.;break;
	      case 3: dphi = (1-x)*lobatto(r1,y)*dlobatto(r2,z)/2.;break; 
	      case 4: dphi = -lobatto(r1,x)*lobatto(r2,y)/2.;break;
	      case 5: dphi = lobatto(r1,x)*lobatto(r2,y)/2.;break;
	      default: throw runtime_error("Invalid face index for hexahedron. ");
              }
	    count++;
	  }
	}
    }
  }         

  else if((i==3)&&(j==0)&&(n>=6)&&(k<(n-3)*(n-4)*(n-5)/6)){// Internal
    int count = 0;
    for(int N=6;N<=n;N++)
      for(int r=2;r<=n-4;r++)
        for(int q=2;q<=n-4;q++)
          for(int p=2;p<=n-4;p++){
            if (p+q+r == N){
	      if (k==count) dphi = lobatto(p,x)*lobatto(q,y)*dlobatto(r,z);
	      count++;
	    }
          }
  }

  else throw runtime_error("Invalid indices for hexahedron");
  return dphi;
}
