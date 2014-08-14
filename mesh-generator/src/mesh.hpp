#ifndef MESH_HPP
#define MESH_HPP

#include "elem.hpp"

class mesh {

private:
  unsigned nvt,nel,grid; //nvt: total number of nodes in the mesh as extracted from the Gambit file
  
public:
  elem* el;  // elements

  mesh(const char [], vector < vector < double> > &vt);
  mesh(const unsigned& igrid, elem* elc);

  ~mesh(){ 
    delete el;
  }
  
  void Read1D(const char infile [], vector < vector < double> > &vt);
  void ReadGambit(const char infile [], vector < vector < double> > &vt);
  void AddAdditionalNodes();   
  void SetMeshNodeInformation();
  void BuildAdjVtx();
  void Buildkmid();
  void Buildkel();
  void set_elr(const unsigned& test=100);
  void set_elr(const  vector <vector <double> > &vt);
  void copy_elr(vector <unsigned>& other_vec) const;
  unsigned MeshNodeNum (const unsigned type) const; 
  unsigned MeshElementNum() const;

  //      <<D[0]->MeshNodeNum (0)<<"\t" // Returns the number of vertices in the mesh.
  //      <<D[0]->MeshNodeNum (1)<<"\t" //Returns the #vertices + #edges in the mesh.
  //      <<D[0]->MeshNodeNum (2)<<"\t" //Returns the #vertices + #edges + #faces in the mesh.
  //      <<D[0]->MeshElementNum()//Returns the total number of elements in the entire mesh.
  unsigned GetGridNumber() const;

};

#endif
