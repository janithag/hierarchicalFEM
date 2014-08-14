#include "basish.hpp"
#include "main.hpp"

double wedgeh::eval_phi(const int* I,const double* X, const int n) const{
  return evalWedgehPhi(I[0],I[1],I[2], X[0], X[1], X[2], n, type); 
}

double wedgeh::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalWedgehDPhiDx(I[0],I[1],I[2], X[0], X[1], X[2], n, type); 
}

double wedgeh::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalWedgehDPhiDy(I[0],I[1],I[2], X[0], X[1], X[2], n, type); 
}

double wedgeh::eval_dphidz(const int* I,const double* X, const int n) const{
  return evalWedgehDPhiDz(I[0],I[1],I[2], X[0], X[1], X[2], n, type); 
}

void wedgeh::eval_sign(const int* o, vector <int>& orientationFlag) const{ 
 for(int i=0;i<dimOfBasis();i++)
    orientationFlag[i]=1;

  int count=6; 
  for(unsigned k=0;k<p-1;k++){
    //edges
    for(int i=0;i<9;i++){
	orientationFlag[count] = pow(o[i],k);
	count++;
    }
   
    //faces 
    if(k>1) 
      	for(int r1=0;r1<=k-2;r1++){
	  for(int r2=k-2;r2>=0;r2--){
	    if(r1+r2==k-2){
	     for(int i=0;i<3;i++){
	      orientationFlag[count] = pow(o[9+i*2],r1)*pow(o[9+i*2+1],r2);
	      count++;
	    }
	    count+=2;  
	    }
	}
      }  
      
    //additional tri face modes
    if(k>0) 
     count += 2;
    
    //interior modes
    if(k>2)   
      count += (k-2)*(k-1)/2;  
  }
 
}

void wedgeh::eval_type(const int* referenceElementType, vector <int>& type) const{ 
  
  for(int i=0;i<dimOfBasis();i++)
    type[i] = 1;

  unsigned count = 6; // 6 vertices
  for(int i=0;i<p-1;i++){ 
   count += 9; //edges
   if(i>1){ // face modes
    for(int j=0;j<i-1;j++)
      for(int iFace=0; iFace<5; iFace++){
	type[count]=referenceElementType[iFace];count++;
      }
   }
   if(i>0)//tri face modes
     for(int iFace=3; iFace<5; iFace++){
	type[count]=referenceElementType[iFace];count++;
      } 
   
   if(i>2)
      count += (i-2)*(i-1)/2; //internal modes 
    }
  }
  
