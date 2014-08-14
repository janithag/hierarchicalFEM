#ifndef BASISH_HPP
#define BASISH_HPP

#include <vector>
#include "friendFuns.hpp"
#include "main.hpp"

using std::vector;

class basish{

protected:  
  int* numOfModes;
  int p;

  
public:
  basish(){
    p=1;
    numOfModes = new int[8];
  }

  basish(const int P){
    numOfModes = new int[8];
    p=P;
  }

  ~basish(){
    delete[] numOfModes;
  }
  
  virtual int spatialDim() const=0;
  virtual int dimOfBasis() const=0; 
  virtual int* getModes() const=0;
  virtual double eval_phi(const int *I,const double* x, const int n=0 )const = 0;
  virtual double eval_dphidx(const int *I,const double* x, const int n=0) const = 0;
  virtual double eval_dphidy(const int *I,const double* x, const int n=0) const = 0;
  virtual double eval_dphidz(const int *I,const double* x, const int n=0) const = 0;
  
  virtual void eval_sign(const int* o, vector <int>& orientationFlag) const = 0;
  virtual void eval_type(const int* referenceElementType, vector <int>& type_) const=0;
};

class lineh: public basish{
  friend double evalLinehPhi(const int i, const int j, const int k, const double x, const int n);
  friend double evalLinehDPhi(const int i, const int j, const int k, const double x, const int n);

public:
  lineh():basish(){};

  lineh(const int P):basish(P){
  }
  
  int spatialDim() const{
    return 1;
  }
  
  int*  getModes()const{
    numOfModes[0]=2; //number of vertices
    numOfModes[1]=1; //number of edges
    numOfModes[2]=0; //number of faces
    numOfModes[3]=0; //number of regions
    numOfModes[4]=1; //number of shape functions in each vertex
    numOfModes[5]=p-1; //number of shape functions in each edge
    numOfModes[6]=0; //number of shape functions in each face
    numOfModes[7]=0; //number of internal modes.

    return numOfModes;
  }

  int dimOfBasis() const{
    return p+1;
  }

  double eval_phi(const int *I,const double* x, const int n=0) const;
  double eval_dphidx(const int *I,const double* x, const int n=0) const;
  double eval_dphidy(const int *I,const double* x, const int n=0) const;
  double eval_dphidz(const int *I,const double* x, const int n=0) const;
  
  void eval_sign(const int* o, vector <int>& orientationFlag)const;
  void eval_type(const int* referenceElementType, vector <int>& type_) const;
};

class quadh: public basish{
  friend double evalQuadhPhi(const int i, const int j, const int k, const double x, const double y, const int n);
  friend double evalQuadhDPhiDx(const int i, const int j, const int k, const double x, const double y, const int n);
  friend double evalQuadhDPhiDy(const int i, const int j, const int k, const double x, const double y, const int n);

public:
  quadh():basish(){};

  quadh(const int P):basish(P){
  }

  int spatialDim() const{
    return 2;
  }
  
  int*  getModes()const{
    numOfModes[0]=4; //number of vertices
    numOfModes[1]=4; //number of edges
    numOfModes[2]=1; //number of faces
    numOfModes[3]=0; //number of regions
    numOfModes[4]=1; //number of shape functions in each vertex
    numOfModes[5]=p-1; //number of shape functions in each edge
    if (p>=4)    //number of shape functions in each face
      numOfModes[6]=(p-2)*(p-3)/2;
    else 
      numOfModes[6]=0;
    numOfModes[7]=0; //number of internal modes.
    return numOfModes;
  }

  int dimOfBasis() const{
    return (p<4) ? 4*p : 4*p+(p-2)*(p-3)/2;
  }

  
   double eval_phi(const int *I,const double* x, const int n=0) const;
   double eval_dphidx(const int *I,const double* x, const int n=0) const;
   double eval_dphidy(const int *I,const double* x, const int n=0) const;
   double eval_dphidz(const int *I,const double* x, const int n=0) const;
   
   void  eval_sign(const int* o, vector <int>& orientationFlag) const;
   void eval_type(const int* referenceElementType, vector <int>& type_) const;
 };

class hexh: public basish{
   friend double evalHexhPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
   friend double evalHexhDPhiDx(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
   friend double evalHexhDPhiDy(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
   friend double evalHexhDPhiDz(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);

 public:
  hexh(const int t):basish(){
    if(t>1||t<0)
      throw runtime_error("invalid hexh type in basish.hpp");
   type=t;
  }

  hexh(const int P,const int t):basish(P){
    if(t>1||t<0)
      throw runtime_error("invalid hexh type in basish.hpp");
    type=t;
  }
  
  int spatialDim() const{
    return 3;
  }
  
  int*  getModes()const{
     numOfModes[0]=8; //number of vertices
     numOfModes[1]=12; //number of edges
     numOfModes[2]=6; //number of quad faces
     numOfModes[3]=0; //number of tri faces
     numOfModes[4]=p-1; //number of shape functions in each edge
     numOfModes[5]=(p>=4) ? (p-2)*(p-3)/2 : 0;  //number of shape functions in each face
     numOfModes[6]=0; // number of basis functions in each tri face
     numOfModes[7]=(p>=6) ? (p-3)*(p-4)*(p-5)/6 : 0; ; //number of internal modes.

    return numOfModes;
  }

  int dimOfBasis() const{
    return  (p<4) ? 12*p-4 : 12*p-4 + 3*(p-2)*(p-3) + (p-3)*(p-4)*(p-5)/6;
  }

   double eval_phi(const int *I,const double* x, const int n=0) const;
   double eval_dphidx(const int *I,const double* x, const int n=0) const;
   double eval_dphidy(const int *I,const double* x, const int n=0) const;
   double eval_dphidz(const int *I,const double* x, const int n=0) const;
   
   void eval_sign(const int* o, vector <int>& orientationFlag)const;
   void eval_type(const int* referenceElementType, vector <int>& type_) const;
   
private:
  int type;
 };
 
class trih: public basish{
   friend double evalTrihPhi(const int i, const int j, const int k, const double x, const double y, const int n);
   friend double evalTrihDPhidx(const int i, const int j,  const int k, const double x, const double y, const int n);
   friend double evalTrihDPhidy(const int i, const int j,  const int k, const double x, const double y, const int n);

 public:
 trih():basish(){};

  trih(const int P):basish(P){
  }
  
  int spatialDim() const{
    return 2;
  }
  
   int*  getModes()const{
     numOfModes[0]=3; //number of vertices
     numOfModes[1]=3; //number of edges
     numOfModes[2]=1; //number of faces
     numOfModes[3]=0; //number of regions
     numOfModes[4]=1; //number of shape functions in each vertex
     numOfModes[5]=p-1; //number of shape functions in each edge
     numOfModes[6]=(p>=3) ? (p-2)*(p-1)/2 : 0;  //number of shape functions in each face
     numOfModes[7]=0; //number of internal modes.
     return numOfModes;
   }

  int dimOfBasis() const{
    return 3*p+(p-2)*(p-1)/2;
  }

   double eval_phi(const int* I,const double* x, const int n=0) const;
   double eval_dphidx(const int* I,const double* x, const int n=0) const;
   double eval_dphidy(const int* I,const double* x, const int n=0) const;
   double eval_dphidz(const int* I,const double* x, const int n=0) const;
   
   void eval_sign(const int* o, vector <int>& orientationFlag)const;
   void eval_type(const int* referenceElementType, vector <int>& type_) const;
 };

class wedgeh: public basish{
  //friend int evalSignW(const int global[], int i, int j, int k);     
   friend double evalWedgehPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
   friend double evalWedgehDPhiDx(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
   friend double evalWedgehDPhiDy(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
   friend double evalWedgehDPhiDz(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);

 public:
 wedgeh(const int t):basish(){
   if(t>5||t<0)
      throw runtime_error("invalid wedgeh type in basish.hpp");
   type=t;
  }

  wedgeh(const int P, const int t):basish(P){
    if(t>5||t<0)
      throw runtime_error("invalid wedgeh type in basish.hpp");
    type=t;
  }

  int spatialDim() const{
    return 3;
  }
  int*  getModes()const{

     numOfModes[0] = 6; //number of vertices
     numOfModes[1] = 9; //number of edges
     numOfModes[2] = 3; //number of quad faces
     numOfModes[3] = 2; //number of tri faces
     numOfModes[4] = p-1; //number of shape functions in each edge
     numOfModes[5] = (p>=4) ? (p-2)*(p-3)/2 : 0; //number of shape functions in each quad face
     numOfModes[6] = (p>=3) ? (p-1)*(p-2)/2 : 0; //number of shape functions in each tri face
     numOfModes[7] = (p>=5) ? (p-4)*(p-3)*(p-2)/6 : 0; //number of internal modes.

    return numOfModes;
  }

  int dimOfBasis() const{
    return 6 + 9*(p-1) + (p-1)*(p-2) + (p>3)*3*(p-2)*(p-3)/2 + (p>4)*(p-2)*(p-3)*(p-4)/6;
  }

   double eval_phi(const int* I,const double* x, const int n=0) const;
   double eval_dphidx(const int* I,const double* x, const int n=0) const;
   double eval_dphidy(const int* I,const double* x, const int n=0) const;
   double eval_dphidz(const int* I,const double* x, const int n=0) const;
   
   void eval_sign(const int* o, vector <int>& orientationFlag) const;
   void eval_type(const int* referenceElementType, vector <int>& type_) const;
   
  private:
    int type; 
 };

class teth: public basish{
  friend double evalTethPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
  friend double evalTethDPhidx(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
  friend double evalTethDPhidy(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
  friend double evalTethDPhidz(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
 
 public:
  teth(const int t):basish(){
    if(t>1||t<0)
      throw runtime_error("invalid teth type in basish.hpp");
    type=t;
  }

  teth(const int P, const int t):basish(P){
    if(t>1||t<0)
      throw runtime_error("invalid teth type in basish.hpp");
    type=t;
  }

  int spatialDim() const{
    return 3;
  }
  
  int*  getModes()const{

     numOfModes[0] = 4; //number of vertices
     numOfModes[1] = 6; //number of edges
     numOfModes[2] = 0; //number of quad faces
     numOfModes[3] = 4; //number of tri faces
     numOfModes[4] = p-1; //number of shape functions in each edge
     numOfModes[5] = 0; //number of quad shape functions in each face
     numOfModes[6] = (p>=3)? (p-1)*(p-2)/2 : 0;  //number of tri shape functions in each face
     numOfModes[7] = (p>=4)? (p-1)*(p-2)*(p-3)/6 : 0; //number of internal modes.

    return numOfModes;
  }

  int dimOfBasis() const{
    return 4 + 6*(p-1) + 2*(p-1)*(p-2) + (p-1)*(p-2)*(p-3)/6;
  }

   double eval_phi(const int* I,const double* x, const int n=0) const;
   double eval_dphidx(const int* I,const double* x, const int n=0) const;
   double eval_dphidy(const int* I,const double* x, const int n=0) const;
   double eval_dphidz(const int* I,const double* x, const int n=0) const;
   
   void eval_sign(const int* o, vector <int>& orientationFlag) const;
   void eval_type(const int* referenceElementType, vector <int>& type_) const;

private:
  int type; 
 };

#endif
