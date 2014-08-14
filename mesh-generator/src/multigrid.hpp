#ifndef MULTIGRID_HPP
#define MULTIGRID_HPP

#include "lsysPDE.hpp"

class elem_type;

class MultiGrid{ 
  
private:
  vector <lsysPDE*> D; // Ku=f at each level
  vector <vector <double> > vt; //grid vertex, vy, vz
  vector <vector <double> > LNode; // Lagrange nodes in the projected mesh.
  vector <vector <double> > Sol; //solution vector (u)
  vector <vector <double> > LSol; //projected solution to Lagrange basis
  vector <vector <unsigned short> > BC; // boundary conditions
  vector <vector <unsigned> > Cell; // h-refined cells
  vector < unsigned> CellType; // h-refined cell types
  vector < map <unsigned,bool> > index;
  vector <int> SolType;
  vector <char*> SolName;
  vector <unsigned> MGIndex;
  
  vector <lsysPDE*> D_old;
  vector <vector <double> > vt_old;
  vector <vector <double> > Sol_old;
  vector <unsigned> elr_old;
  unsigned maxOrder;
  unsigned p;
  unsigned lagrangeP;
  bool time_test;
  unsigned short gridn, gridr;
  vector < vector <const elem_type*> > type_elem; //const elem_type* type_elem[3][3]
  vector <unsigned> DOF; // Number of degrees of freedom at each grid level.
  vector <unsigned> DofstoSolve; //number of dofs to solve at each level

public:
  MultiGrid(const char mesh_file[], const unsigned short& gridP);
  
  ~MultiGrid();

  void setDOF(); // sets the number of degrees of freedom at each grid level.
  void solve(const unsigned short &time); 
  void GenerateBC(const char name[]); //const unsigned short& fGrid, const unsigned short& cGrid=0
  void generate_vt(const unsigned short &grid_end, const unsigned short& grid_start=1);
  
  void Print_vtk_ASCII(const unsigned& time, const int type=0);
  void L2Error() const;
  void GenerateLagrangeSolution2D(int LagrangeOrder);
  void GenerateLagrangeSolution3D(int LagrangeOrder);
  int LagrangeHRefinement2D();
  int LagrangeHRefinement3D();
  
  void AddSolutionVector(const char name[], const unsigned order=0); 
  void ResizeSolutionVector( const char name[]);
  void CheckVectorSize(const unsigned &i);
  void Initialize(const char name[]);
  unsigned GetIndex(const char name[]);
  
  void ClearMGIndex();
  void AddToMGIndex(const char name[]);
  void InitMultigrid();

  int Restrictor(unsigned gridf);
  int Prolongator(unsigned gridf);
  void ProlongatorSol(unsigned gridf);

  int BuildProlongatorMatrix(unsigned gridf);
};

#endif
