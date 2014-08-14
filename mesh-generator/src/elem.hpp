// mesh is stored here
#ifndef __elem_hpp__
#define __elem_hpp__

#include "main.hpp"
#include "basish.hpp"

// change to vertices, vertices+edges, vertices+edges+quadfaces, vertices+edges+faces, vertices+edges+faces+interior
const unsigned NVE[6][5]={{8,20,26,26,27},//hex
			  {4,10,10,14,15}, //tet
			  {6,15,18,20,21}, //wedge
			  {4,8,8,8,9},//quad
			  {3,6,6,6,7},//tri 
			  {2,3,3,3,3}};//line


class elem{
  
private:
  
  int** kel; // adjacent element to a facet
  int* kel_memory;
  
  unsigned** kvtel;
  unsigned* kvtel_memory;
  
  unsigned** kvert; // global indices of each node of each element
  unsigned* kvert_memory;
  
  unsigned* elr;
  unsigned* nve;
  unsigned nvt,nv0,nv1,nv2,nv3; // nv0: number of vertices of the mesh, nv1: number of edges, nv2: number of quad faces, nv3: number of tri faces.
  unsigned nel,nelt[6];
  unsigned nelr,nelrt[6];

  short unsigned* elt;
  vector <vector <int> > globalDOF; // global degree of freedom 
  vector < vector <int> > referenceElementType;
  vector <vector <int> > o;
public:

  elem(const unsigned & other_nel);
  elem(const elem *elc);
  ~elem(){ 
    delete [] kvert_memory;
    delete [] kvert;
    delete [] kel_memory;
    delete [] kel;
    delete [] elt;
    delete [] elr;
    delete [] kvtel_memory; 
    delete [] kvtel;
    delete [] nve;
  }
  
 
  unsigned ElementNodeNum(const unsigned& iel,const unsigned& type=3) const;// returns the number of 0-vertices, 1-vertices+edges, 2-vertices+edges+faces, 3-vertices+edges+faces+interiors for a type of element. See function definition, see array NVE[][] above. 
  unsigned GlobalNodeIndex(const unsigned &iel,const unsigned &inode)const; 
  const unsigned* NodeAddress(const unsigned &iel,const unsigned &inode)const;
  void SetElementNodeGlobalIndex(const unsigned& iel,const unsigned& inode, const unsigned& value);
 
  unsigned FacetNodeGlobalIndex(const unsigned &iel,const unsigned &iface, const unsigned &inode) const;
  short unsigned ElementType(const unsigned &iel) const;
  void SetElementType(const unsigned &iel, const short unsigned &value);
 
  int FacetAdjacentElemIndex(const unsigned &iel,const unsigned &iface) const;
  void SetFacetAdjacentElemIndex(const unsigned &iel,const unsigned &iface, const int &value); 
  unsigned GetIndex(const char name[]) const;

  unsigned MeshElementNum(const char* name="All") const;// Gives the total number of elements in the entire mesh.
  void AddToElementNumber(const unsigned &value, const char name[]);  
  void AddToElementNumber(const unsigned &value, short unsigned ielt);

  unsigned GetRefinedElementNumber(const char name[]="All") const;
  unsigned GetRefinedElementNumber(short unsigned ielt) const;
  void AddToRefinedElementNumber(const unsigned &value, const char name[]="All");
  void AddToRefinedElementNumber(const unsigned &value, short unsigned ielt);
  void InitRefinedToZero();
  unsigned GetRefinedElementIndex(const unsigned &iel) const;
  void SetRefinedElementIndex(const unsigned &iel, const unsigned &value);
 
  unsigned MeshNodeNum()const;
  unsigned MeshVertexNum()const;//Gives the total number of vertices in the entire mesh.
  unsigned MeshEdgeNum()const;// Gives the total number of edges in the entire mesh.
  unsigned MeshQuadFaceNum()const;//Gives the total number of quadfaces in the entire mesh. 
  unsigned MeshTriFaceNum()const; //Gives the total number of trifaces in the entire mesh.
  const unsigned NumOfDof(const int p) const; // returns the number of degrees of freedom for a given p.
  
  int TriIndex(const int iel) const;
  int TetIndex(const int iel) const;
  int WedgeIndex(const int iel) const;
  int HexIndex(const int iel) const;
  void SetGlobalDOF(const int P);

  void SetNodeNumber(const unsigned &value);
  void SetMeshNumOfVertices(const unsigned &value);
  void SetMeshNumOfEdges(const unsigned &value);
  void SetMeshQuadFaceNum(const unsigned &value);
  void SetMeshTriFaceNum(const unsigned& value);
  void SetOrientationFlagsAndType();
  vector <int> RefElementType(const unsigned& iel)const;
  
  unsigned ElementFacetNum(const unsigned &iel,const unsigned &type=1)const;
  unsigned GetElementSquareFaceNumber(const unsigned &iel)const;
  unsigned GetElementTriangleFaceNumber(const unsigned &iel)const;
  
  void AllocateVertexElementMemory();
  unsigned NodeNumOfAdjacentElements(const unsigned &inode)const;
  unsigned NodeAdjacentElementIndex(const unsigned &inode,const unsigned &jnode)const;
  const unsigned* NodeAdjacentElementAddress(const unsigned &inode,const unsigned &jnode)const;
  void SetVertexAdjacentElementIndex(const unsigned &inode,const unsigned &jnode, const unsigned &value);
  
  int* o_(unsigned iel){
    return &o[iel][0];
  }
  
  int* type_(unsigned iel){
    return &referenceElementType[iel][0];
  }
  
  int* globalDOF_(const unsigned iel){
    return &globalDOF[iel][0];
  }
};

#endif
