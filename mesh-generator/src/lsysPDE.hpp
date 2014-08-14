#ifndef LSYSPDE_HPP
#define LSYSPDE_HPP

#include "mesh.hpp"

class elem_type;

class lsysPDE: public mesh{
  friend class MultiGrid;
 
private: 
  vector <Vec> Res;  // residual at each level
  vector <Vec> Eps;  // solution of the residual equation at each level (error)
 
  vector <int> SolType;
  vector <char*> SolName;
  
  Mat K,CC; // K: stiffness matrix
  Mat R; // Restriction matrix
  Mat P; // Prolongation matrix
  Vec RES, RESC; // Residual at each level: RES
  Vec EPS, EPSC; // Error at each level: EPS (K.EPS = RES)
  vector <int> KIndex;  
  PetscErrorCode ierr;
  vector < vector <int> > sign; // orientationFlag for each element
  vector < vector <int> > type; // type of the reference element for each element 
  
public:
  unsigned p; //0-based index : order of each level
  lsysPDE(const char infile[], vector < vector < double> > &vt, unsigned order);
  lsysPDE(const unsigned& igrid, elem* elc);

  ~lsysPDE(){ 
    for(unsigned i=0; i<SolName.size(); i++)
      delete[] SolName[i];
  };
 
  int printsol_vtk_ASCII(const  vector <vector <double> > &vt, const unsigned& time, const char type[]="linear");
  void AddSolutionVector(const char name[], const unsigned order);
  int FreeSolutionVectors();
  int ResizeSolutionVector(const char name[]);
  unsigned GetIndex(const char name[]);
  int SetResZero(const vector <unsigned> &MGIndex); 
  int SetEpsZero(const vector <unsigned> &MGIndex);
  int SumEpsCToEps(const vector <unsigned> &MGIndex);
  int SumEpsToSol(vector <vector <double> > &Sol,const vector <unsigned> &MGIndex);
  int FindResMax(const vector < vector <unsigned short> >  &BC, const vector <unsigned> &MGIndex);
  int UpdateResidual();
  vector <double> ResInf;// residual infinity norm at each level for each variable //TODO make private and provide function
  void NodeGlobalDofOfModes(const unsigned iCell, vector <unsigned>& iMode) const;
  void SetSign();
    
  int* sign_(unsigned iel){
    return &sign[iel][0];
  }
  
  int* type_(unsigned iel){
    return &type[iel][0];
  }
  
  int InitMultigrid(const vector <unsigned> &MGIndex);
  int DeallocateMatrix(const vector <unsigned> &MGIndex);
  
  int AssembleMatrix(vector < vector <const elem_type*> >  type_elem, vector <vector <double> > &vt,  const unsigned& gridn, const char type[]="All");
  int AssembleResidual(vector < vector <const elem_type*> > type_elem, vector <vector <double> > &vt, const vector <vector <double> > &Sol,  const unsigned& gridn, const char type[]="All", const unsigned& pOrder=0);
  int VankaPETSCSmoother(const vector < vector <unsigned short> > &BC, const vector <unsigned> &MGIndex, const int& DofstoSolve);
  int SpaceDecompositionSmoother(const vector  <vector <unsigned short> >& BC, const vector <unsigned>& MGIndex, const int& DofstoSolve);

private:
  
  static PetscInt* node;
  static double A[27*27];
  //static double*  B[3][3];//double B[4][4][27*27];
  vector < vector <double> > vertex;//,vy[27],vz[27]; vx and vertex are the same. introduced two to resolve conflicts.
  static double vx[3][27];
  
  static double Weight;
  static double phi1[27],gradphi1[27][3],Weight1;  
  static double _beta,_gamma,_alpha,_Td;
};

#endif
