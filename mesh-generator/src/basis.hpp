#ifndef BASIS_HPP
#define BASIS_HPP

#include "friendFuns.hpp"
#include <vector>
#include <cmath>

using std::vector;

class basis {
public:
  virtual double eval_phi(const int *I,const double* x, const int n=0) const = 0;
  virtual double eval_dphidx(const int *I,const double* x, const int n=0) const = 0;
  virtual double eval_dphidy(const int *I,const double* x, const int n=0) const = 0;
  virtual double eval_dphidz(const int *I,const double* x, const int n=0) const = 0;
  virtual void RefinedNodes(const int p, double** x) const = 0;
  virtual int ElementDOF(const int p) const = 0;
};

class line: public basis{
  friend double evalPhi(const int i, const double x, const int n);
  friend double evalDphi(const int i, const double x, const int n);

private:
  void eval_phi(const double x, double *phi, const int n) const;
  void eval_dphidx(const double x, double* dphi_dx, const int n) const;

public:
  double eval_phi(const int *I,const double* x, const int n=0) const;
  double eval_dphidx(const int *I,const double* x, const int n=0) const;
  double eval_dphidy(const int *I,const double* x, const int n=0) const;
  double eval_dphidz(const int *I,const double* x, const int n=0) const;
  void RefinedNodes(const int p, double** x) const;
  int ElementDOF(const int p) const;
};

class quad: public basis{
   friend double evalPhi(const int i, const double x, const int n);
   friend double evalDphi(const int i, const double x, const int n);
 public:
   double eval_phi(const int *I,const double* x, const int n=0) const;
   double eval_dphidx(const int *I,const double* x, const int n=0) const;
   double eval_dphidy(const int *I,const double* x, const int n=0) const;
   double eval_dphidz(const int *I,const double* x, const int n=0) const;
   void RefinedNodes(const int p, double** x) const;
   int ElementDOF(const int p) const;
 };

class hex: public basis{
   friend double evalPhi(const int i, const double x, const int n);
   friend double evalDphi(const int i, const double x, const int n);
 public:
   double eval_phi(const int *I,const double* x, const int n=0) const;
   double eval_dphidx(const int *I,const double* x, const int n=0) const;
   double eval_dphidy(const int *I,const double* x, const int n=0) const;
   double eval_dphidz(const int *I,const double* x, const int n=0) const;
   void RefinedNodes(const int p, double** x) const;
   int ElementDOF(const int p) const;
 };
 

 class tri: public basis{
   friend double evalTriPhi(const int i, const int j, const double x, const double y, const int n);
   friend double evalTriDPhidx(const int i, const int j, const double x, const double y, const int n);
   friend double evalTriDPhidy(const int i, const int j, const double x, const double y, const int n);
 public:
   double eval_phi(const int* I,const double* x, const int n=0) const;
   double eval_dphidx(const int* I,const double* x, const int n=0) const;
   double eval_dphidy(const int* I,const double* x, const int n=0) const;
   double eval_dphidz(const int* I,const double* x, const int n=0) const;
   void RefinedNodes(const int p, double** x) const;
   int ElementDOF(const int p) const;
 };

 class wedge: public basis{
   friend double evalPhi(const int i, const double x, const int n);
   friend double evalDphi(const int i, const double x, const int n);
   friend double evalTriPhi(const int i, const int j, const double x, const double y, const int n);
   friend double evalTriDPhidx(const int i, const int j, const double x, const double y, const int n);
   friend double evalTriDPhidy(const int i, const int j, const double x, const double y, const int n);
 public:
   double eval_phi(const int* I,const double* x, const int n=0) const;
   double eval_dphidx(const int* I,const double* x, const int n=0) const;
   double eval_dphidy(const int* I,const double* x, const int n=0) const;
   double eval_dphidz(const int* I,const double* x, const int n=0) const;
   void RefinedNodes(const int p, double** x) const;
   int ElementDOF(const int p) const;
 };
  
class tet: public basis{
   friend double evalTetraPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n);
   friend double evalTetraDPhidx(const int i, const int j, const int k, const double x, const double y, const double z, const int n);
   friend double evalTetraDPhidy(const int i, const int j, const int k, const double x, const double y, const double z, const int n);
   friend double evalTetraDPhidz(const int i, const int j, const int k, const double x, const double y, const double z, const int n);
 public:
   double eval_phi(const int* I,const double* x, const int n=0) const;
   double eval_dphidx(const int* I,const double* x, const int n=0) const;
   double eval_dphidy(const int* I,const double* x, const int n=0) const;
   double eval_dphidz(const int* I,const double* x, const int n=0) const;
   void RefinedNodes(const int p, double** x) const;
   int ElementDOF(const int p) const;
 };

#endif
