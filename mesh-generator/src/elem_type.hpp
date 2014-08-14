#ifndef ELEM_TYPE_HPP
#define ELEM_TYPE_HPP
#include "basis.hpp"
#include "basish.hpp"
#include <string>
#include <algorithm>
#include <stdexcept>
#include "main.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::max;

class elem;

class elem_type{

public:
  /* Note on gauss order n: if the polynomial to integrate is x^n1 y^n2 z^n3, then following rule apply:
     line, quad, hex: n = max(n1, n2, n3).  tri, tet, wedge: n = n1 + n2 + n3 */ 
  elem_type(const int solid, const int p, const int n, const int orderOfMapping);
  ~elem_type();

  void ProjectionMatrix2D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const;
  void ProjectionMatrix3D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const;
  void (elem_type::*projectionMatrix)(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const;

  int** getIndices() const{
    return IND;
  }

  unsigned ElementDOF() const{
    return numOfShapeFuns;  
  }

  double* GaussPoint(int gaussIndex)const{
    return X[gaussIndex];
  }
  
  unsigned NumOfGaussPts() const{
    return numOfGaussPts;
  }

  double* GaussWeight() const{
    return gaussWeights;
  }

  double** Phi() const{
    return phi[0];  
  }

  double** DphiDxi() const{
    return dphidxi[0];  
  }

  double** DphiDeta() const{
    return dphideta[0];  
  }

  double** DphiDzeta() const{
    return dphidzeta[0];  
  }

  void mapSpatialToRefElem1D(vector<vector<double> > node, const unsigned& gaussIndex, const int* sign, double& weight, double* phi_, double** gradPhi, double* X, const int* type) const;
  void mapSpatialToRefElem2D(vector<vector<double> > node, const unsigned& gaussIndex, const int* sign, double& weight, double* phi_, double** gradPhi, double* X, const int* type) const;
  void mapSpatialToRefElem3D(vector<vector<double> > node, const unsigned& gaussIndex, const int* sign, double& weight, double* phi_, double** gradPhi, double* X, const int* type) const;
  void (elem_type::*spatialToRefPtr)(vector<vector<double> > node, const unsigned& gaussIndex, const int* sign, double& weight, double* phi_, double** gradPhi, double* X, const int* type) const;

private:
  int numOfShapeFuns;
  int numOfGaussPts;
  int element;
  double** X;
  int** IND;
  int** I; //lagrange basis indices.
  int numOfNodes;  //modify: store in the class.
  int dim;
  basish* pt_h[6];
  basis* pt_lag;
  double* x;
  unsigned hierarchicP;
  
  
  double* gaussWeights;
  double*** phi;
  double*** dphidxi;
  double*** dphideta;
  double*** dphidzeta;
  double*** Dphi; //hierarchical basis derivatives
  double*** dphi; //lagrange basis derivatives
  double** lphi; //lagrange basis functions
 
  void error(const string& msg) const;

};

#endif
