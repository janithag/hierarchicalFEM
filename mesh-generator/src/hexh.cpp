#include "basish.hpp"
#include "main.hpp"

double hexh::eval_phi(const int* I,const double* X, const int n) const{
  return evalHexhPhi(I[0], I[1], I[2], X[0], X[1], X[2], n, type);
}

double hexh::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalHexhDPhiDx(I[0], I[1], I[2], X[0], X[1], X[2], n, type);
}

double hexh::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalHexhDPhiDy(I[0], I[1], I[2], X[0], X[1], X[2], n, type);
}

double hexh::eval_dphidz(const int* I,const double* X, const int n) const{
    return evalHexhDPhiDz(I[0], I[1], I[2], X[0], X[1], X[2], n, type);
}

void hexh::eval_sign(const  int* o, vector <int>& orientationFlag) const{ 
  for(int i=0;i<dimOfBasis();i++)
    orientationFlag[i]=1;

 /* assigning the orienation flags */ 
  
  int count=8; 
  for(unsigned k=0;k<p-1;k++){
    //edges
    for(int i=0;i<12;i++){
	orientationFlag[count] = pow(o[i],k);
	count++;
    }
    //faces

    if(k>1){ 
     
	for(int r1=0;r1<=k-2;r1++){
	  for(int r2=k-2;r2>=0;r2--){
	    if(r1+r2==k-2)
	     for(int i=0;i<6;i++){
	      orientationFlag[count] = pow(o[12+i*2],r1)* pow(o[12+i*2+1],r2);
	      count++;
	    }
	}
      }   
    }
    //interior modes
    if(k>3)   
      count += (k-3)*(k-2)/2;  
  }

}

void hexh::eval_type(const int* referenceElementType, vector <int>& type) const{ 

  for(int i=0;i<dimOfBasis();i++){
    type[i] = 1;
    }

  unsigned count = 8; // 8 vertices
  for(int i=0;i<p-1;i++){ 
   count += 12; //edges
   if(i>1){ // face modes
    for(int j=0;j<i-1;j++)
      for(int iFace=0; iFace<6; iFace++){
	type[count]=referenceElementType[iFace];count++;
      }
   }
   if(i>3)
      count += (i-2)*(i-3)/2; //internal modes 
    }
  }