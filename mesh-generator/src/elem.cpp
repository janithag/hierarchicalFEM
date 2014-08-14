#include "elem.hpp" 
#include "basis.hpp"
/**
 * Number of elements obtained with one refinement
**/
const unsigned NRE[6]={8,8,8,4,4,2};
 
/** 
 * Number of facets: FACES(3D), edges(2D) or point-etrema(1D) for each considered element
 **/
// first entry: number of quadrilateral facets, 2nd entry: number of facets
const unsigned NFC[6][2]={{6,6}, //hex
 			  {0,4}, //tet
			  {3,5}, //wedge
			  {0,4}, //quad
			  {0,3}, //tri
			  {0,2}}; //line
			  
/**
 * Node ordering for each element facet position for each considered element
 * (vertices first)
 **/
const unsigned ig[6][6][9]={{{0,1,5,4,8,17,12,16,20}, //front
			     {1,2,6,5,9,18,13,17,21}, //right
			     {2,3,7,6,10,19,14,18,22}, //back
			     {3,0,4,7,11,16,15,19,23}, //left
			     {0,3,2,1,11,10,9,8,24}, //bottom
			     {4,5,6,7,12,13,14,15,25}}, //top
			     
			    {{0,2,1,6,5,4,10}, //bottom 
			     {0,1,3,4,8,7,11},
			     {1,2,3,5,9,8,12},
			     {2,0,3,6,7,9,13}},
			     
			    {{0,1,4,3,6,13,9,12,15},
			     {1,2,5,4,7,14,10,13,16},
			     {2,0,3,5,8,12,11,14,17},
			     {0,2,1,8,7,6,18},
			     {3,4,5,9,10,11,19}},
			     
			    {{0,1,4},
			     {1,2,5},
			     {2,3,6},
			     {3,0,7}},
			     
			    {{0,1,3},
			     {1,2,4},
			     {2,0,5}},
			     
			    {{0},
			     {1}}};
/**
 * This constructor allocates the memory for the \textit{coarsest elem}
 **/
elem::elem(const unsigned& other_nel){

  nelt[0]=nelt[1]=nelt[2]=nelt[3]=nelt[4]=nelt[5]=0;
  nel=other_nel;
  
  elt = new unsigned short [nel]; //element type
  elr = new unsigned [nel]; //if the element is refined
  
  kvert=new unsigned* [nel]; //topology
  kel=new int* [nel]; //guess a neighbour element list for a facet
  
  kvert_memory = new unsigned [nel*NVE[0][4]];  
  kel_memory = new int [nel*NFC[0][1]];
  for(unsigned i=0; i<nel*NFC[0][1]; i++)
    kel_memory[i]=-1;
  
  unsigned* pt_u = kvert_memory;
  int* pt_i = kel_memory; 
  
  //local indices
  for (unsigned i=0;i<nel;i++){
    kvert[i] = pt_u;
    pt_u += NVE[0][4]; //NVE[4]: number of vertices+edges+faces+interiors of hex 
    kel[i] = pt_i;
    pt_i += NFC[0][1]; //NFC: number of facets of hex
  }
}

/**
 * This constructor allocates the memory for the \textit{finer elem} 
 * starting from the paramenters of the \textit{coarser elem}
 **/
elem::elem(const elem* elc){
  nelt[0]=nelt[1]=nelt[2]=nelt[3]=nelt[4]=nelt[5]=0;
  nel=elc->GetRefinedElementNumber()*REF_INDEX;
    
  elt=new unsigned short [nel];
  elr=new unsigned [nel];
  
  kvert=new unsigned*[nel];
  kel=new int* [nel];
  
  unsigned kvert_size=0;
  unsigned kel_size=0;
  for(unsigned i=0;i<6;i++){
    kvert_size+=elc->GetRefinedElementNumber(i)*NVE[i][4];
    kel_size+=elc->GetRefinedElementNumber(i)*NFC[i][1];
  }
  kvert_size*=REF_INDEX;
  kel_size*=REF_INDEX;
  
  kvert_memory=new unsigned [kvert_size];    
  kel_memory=new int [kel_size];
  for(unsigned i=0;i<kel_size; i++)
    kel_memory[i]=0;
  
  int *pt_i=kel_memory;
  unsigned *pt_u=kvert_memory;
  unsigned jel=0;
  for (unsigned iel=0;iel<elc->MeshElementNum();iel++){
    if( elc->GetRefinedElementIndex(iel) ){
      short unsigned elemt=elc->ElementType(iel);
      for(unsigned j=0;j<NRE[elemt];j++){
	kvert[jel+j]=pt_u;
	pt_u+=elc->ElementNodeNum(iel);
	
	kel[jel+j]=pt_i;
	pt_i+=elc->ElementFacetNum(iel);
      }
      jel+=NRE[elemt];
    } 
  }
}

/**
 * Return the number of vertices(type=0) + vertices+edges(type=1) + vertices+edges+quadfaces(type=2) 
 * vertices+edges+faces(type=3) + vertices+edges+faces+interior(type=4) 
 **/
unsigned elem::ElementNodeNum(const unsigned& iel, const unsigned& type) const{
  return NVE[elt[iel]][type];
}

/**
 * Return the local->global node number
 **/
unsigned elem::GlobalNodeIndex(const unsigned& iel, const unsigned& inode)const{
  return kvert[iel][inode];
}

vector <int> elem::RefElementType(const unsigned& iel)const{
return referenceElementType[iel];
}
/**
 * Return the local->global node address
 **/
const unsigned* elem::NodeAddress(const unsigned& iel, const unsigned& inode)const{
  return &kvert[iel][inode];
}

/**
 * Set the local->global node number in kvert
 **/
void elem::SetElementNodeGlobalIndex(const unsigned& iel, const unsigned& inode, const unsigned& value){
  kvert[iel][inode]=value;
}

/**
 * Return the global node index of the facet ifacet for local node inode
 **/
unsigned elem::FacetNodeGlobalIndex(const unsigned& iel, const unsigned& ifacet, const unsigned& inode)const{
  return kvert[iel][ig[elt[iel]][ifacet][inode]];
}

/* Based on the orientation of the spatial element, this function stores the orientation flags and referenceElementType for each element */
void elem::SetOrientationFlagsAndType(){ // tested 
  
  int minIndex, maxIndex;
  unsigned nNodes;
  referenceElementType.resize(MeshElementNum());
  o.resize(MeshElementNum());
  
  // Some information needed to set face orientation
  const unsigned iFace[6][4] = {{0,1,4,5}, {1,2,5,6}, {3,2,7,6}, {0,3,4,7}, {0,1,3,2}, {4,5,7,6}}; //hex faces
  const unsigned wFace[5][4] = {{0,1,3,4},{1,2,4,5},{2,0,5,3},{0,1,2},{3,4,5}}; //wedge faces
  const unsigned possibleCombinations[8][3] = {{0,1,2},{1,0,3},{2,3,0},{3,2,1},{0,2,1},{1,3,0},{2,0,3},{3,1,2}};
  const int signs[8][3] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1}};
  
  for(unsigned iel=0;iel<MeshElementNum();iel++){
    //cout<<"element: "<<iel<<endl;
    if(ElementType(iel)==0)
      referenceElementType[iel].resize(6);
    else if (ElementType(iel)==2)
      referenceElementType[iel].resize(5);
    else {
      referenceElementType[iel].resize(1);
      referenceElementType[iel][0]=0; 
    }
    //cout<<"\n\nElement index: "<<iel<<endl;    
    switch(ElementType(iel)) {
      case 0: {// Hexahedron 
	  nNodes=8;
	  o[iel].resize(24); // edge+face orientation flags 
	  int globalIndex[nNodes];
	  
	  for(int i=0;i<24;i++)
	    o[iel][i]=1;
	  
	  for(unsigned i=0;i<nNodes;i++)  
	    globalIndex[i] = GlobalNodeIndex(iel,i)-1u;
	  
	  /* Edges */
	  if(globalIndex[0] > globalIndex[1])
	    o[iel][0]=-1;
  
	  if(globalIndex[1] > globalIndex[2])
	    o[iel][1]=-1;
  
	  if(globalIndex[3] > globalIndex[2])
	    o[iel][2]=-1;
  
	  if(globalIndex[0] > globalIndex[3])
	    o[iel][3]=-1;
  
	  if(globalIndex[4] > globalIndex[5])
	    o[iel][4]=-1;
  
	  if(globalIndex[5] > globalIndex[6])
	    o[iel][5]=-1;
  
	  if(globalIndex[7] > globalIndex[6])
	    o[iel][6]=-1;
  
	  if(globalIndex[4] > globalIndex[7])
	    o[iel][7]=-1;
  
	  if(globalIndex[0] > globalIndex[4])
	    o[iel][8]=-1;
  
	  if(globalIndex[1] > globalIndex[5])
	    o[iel][9]=-1;
  
	  if(globalIndex[2] > globalIndex[6])
	    o[iel][10]=-1;
  
	  if(globalIndex[3] > globalIndex[7])
	    o[iel][11]=-1;
  
  /* Faces */ 
    
  int minIndex[6][3];
  int min;
  
  // Find A
  for(int i=0; i<6; i++){
    min=globalIndex[iFace[i][0]];
    minIndex[i][0]=0;
    for(int j=1;j<4;j++)
      if(globalIndex[iFace[i][j]]<min){
	minIndex[i][0]=j;
        min = globalIndex[iFace[i][j]];
      }
  }

  int k;
  
  for(int i=0;i<6;i++){
	k=minIndex[i][0];
	minIndex[i][1]=possibleCombinations[k][1];
	minIndex[i][2]=possibleCombinations[k][2];
	
	if(globalIndex[iFace[i][minIndex[i][1]]] > globalIndex[iFace[i][minIndex[i][2]]])
	  k+=4;
	o[iel][12+i*2]=signs[k][0];
	o[iel][12+i*2+1]=signs[k][1];
	referenceElementType[iel][i]=(signs[k][2]==1); // type 0: 03 = -1, type 1: 03 = +1//
      
  }
      break;
    }
  case 1:{
    nNodes = 10;
    int index[nNodes], temp[nNodes];
    int fIndex[4], adjElem[4], fTemp[4];
    int globalIndex[nNodes];
    
    /* vertex global indices of the tetrahedron */
    for(unsigned i=0;i<nNodes;i++){   
      globalIndex[i] = GlobalNodeIndex(iel,i)-1u;
      index[i]=i;
    }
    
    for(unsigned i=0;i<4;i++)
      fIndex[i]=i;
  
    minIndex = std::min_element(globalIndex, globalIndex + 4)-globalIndex;
  
    /* align local index 0 with minimum global index */
    if (minIndex==3){
      index[0]=3; //vertices
      index[1]=0;
      index[3]=1;
      index[4]=7; //edges
      index[5]=6;
      index[6]=9;
      index[7]=8;
      index[8]=4;
      index[9]=5;
      fIndex[0]=3;
      fIndex[2]=0;
      fIndex[3]=2;
    }
    
    else
      while(index[0]!=minIndex){ 
	for(int i=0;i<3;i++){
  	index[i]=(index[i]+1)%3; // vertices
	index[i+4]=(index[4+i])%3+4; //edges
	index[i+7]=(index[7+i])%3+7;	
	fIndex[i+1]=(fIndex[i+1])%3+1;//faces
	}
      }
    
    maxIndex = std::max_element(globalIndex,globalIndex + 4)-globalIndex;
    int Index=0;
    while(index[Index]!=maxIndex) Index++;
    maxIndex=Index;
    
    /* align local index 3 with maximum global index by rotating the face opposite to vertex 0 */ 
    if (3 != maxIndex)
    {
      for(unsigned i=0;i<nNodes;i++)
	temp[i]=index[i];
      
      for(unsigned i=0;i<4;i++)
	fTemp[i]=fIndex[i];
    
      index[1]=temp[3];
      index[2]=temp[1];
      index[3]=temp[2];
      index[4]=temp[7];
      index[5]=temp[8];
      index[6]=temp[4];
      index[7]=temp[6];
      index[8]=temp[9];
      index[9]=temp[5];
      fIndex[0]=fTemp[1];
      fIndex[1]=fTemp[3];
      fIndex[3]=fTemp[0];
  }
   
  maxIndex = std::max_element(globalIndex,globalIndex + 4)-globalIndex;
  Index=0;
  while(index[Index]!=maxIndex) Index++;
    maxIndex=Index;
  
  if (3 != maxIndex)
  {
    for(unsigned i=0;i<nNodes;i++)
      temp[i]=index[i];
    
    for(unsigned i=0;i<4;i++)
	fTemp[i]=fIndex[i];
    
      index[1]=temp[3];
      index[2]=temp[1];
      index[3]=temp[2];
      index[4]=temp[7];
      index[5]=temp[8];
      index[6]=temp[4];
      index[7]=temp[6];
      index[8]=temp[9];
      index[9]=temp[5];
      fIndex[0]=fTemp[1];
      fIndex[1]=fTemp[3];
      fIndex[3]=fTemp[0];
  }
  
  if(globalIndex[index[1]]>globalIndex[index[2]])
    referenceElementType[iel][0]=1;
  
  for(unsigned inode=0;inode<nNodes;inode++)
    SetElementNodeGlobalIndex(iel, inode, globalIndex[index[inode]]+1u);
  
  for(unsigned iface=0;iface<4;iface++)
    adjElem[iface]=FacetAdjacentElemIndex(iel,iface);
    
  for(unsigned iface=0;iface<4;iface++)
    SetFacetAdjacentElemIndex(iel,iface,adjElem[fIndex[iface]]);
  break;
  }
  case 2:{//wedge
    nNodes=6;
    o[iel].resize(15);  //9 + 3x2 
    
    for(int i=0;i<15;i++)
	o[iel][i]=1;
    
    int globalIndex[nNodes];
    //cout<<"global indices"<<endl;
    for(unsigned i=0;i<nNodes;i++) {
      globalIndex[i] = GlobalNodeIndex(iel,i)-1u;
    //cout<<globalIndex[i]<<"\t";
  }
      if(globalIndex[0] > globalIndex[1])
	o[iel][0]=-1;
  
      if(globalIndex[1] > globalIndex[2])
	o[iel][1]=-1;
  
      if(globalIndex[2] > globalIndex[0])
	o[iel][2]=-1;
  
      if(globalIndex[3] > globalIndex[4])
	o[iel][3]=-1;
      
      if(globalIndex[4] > globalIndex[5])
	o[iel][4]=-1;
  
      if(globalIndex[5] > globalIndex[3])
	o[iel][5]=-1;
  
      if(globalIndex[0] > globalIndex[3])
	o[iel][6]=-1;
  
      if(globalIndex[1] > globalIndex[4])
	o[iel][7]=-1;
      
      if(globalIndex[2] > globalIndex[5])
	o[iel][8]=-1;
      
  /* Faces */ 
    
  int minIndex[5][3];
  int min;
  
  // quadrilateral faces
  for(int i=0; i<3; i++){
    min=globalIndex[wFace[i][0]];
    minIndex[i][0]=0;
    for(int j=1;j<4;j++)
      if(globalIndex[wFace[i][j]]<min){
	minIndex[i][0]=j;
        min = globalIndex[wFace[i][j]];
      }
  }
  
  int k;
  
  for(int i=0;i<3;i++){
	k=minIndex[i][0];
	minIndex[i][1]=possibleCombinations[k][1];
	minIndex[i][2]=possibleCombinations[k][2];
	
	if(globalIndex[wFace[i][minIndex[i][1]]] > globalIndex[wFace[i][minIndex[i][2]]])
	  k+=4;
	o[iel][9+i*2]=signs[k][0];
	o[iel][9+i*2+1]=signs[k][1];
	referenceElementType[iel][i]=(signs[k][2]==1); // type 1: 03 = +1, type 0: 03 = -1//  
  }
  //cout<<"ref element type"<<endl;
  // triangular faces
  for(int i=3; i<5; i++){
    min=globalIndex[wFace[i][0]];
    minIndex[i][0]=0;
    for(int j=1;j<3;j++)
      if(globalIndex[wFace[i][j]]<min){
	minIndex[i][0]=j;
        min = globalIndex[wFace[i][j]];
      }
      
    referenceElementType[iel][i]=minIndex[i][0];
    if(globalIndex[wFace[i][(minIndex[i][0]+1)%3]] > globalIndex[wFace[i][(minIndex[i][0]+2)%3]])
      referenceElementType[iel][i]+=3;
    //cout<<referenceElementType[iel][i]<<"\t";
  }
  
      
    break;
  }
  case 3:{ //quad
      nNodes=4;
      o[iel].resize(nNodes);  
      for(int i=0;i<nNodes;i++)
	o[iel][i]=1;
      
      int globalIndex[nNodes];
	  
      for(unsigned i=0;i<nNodes;i++)  
	globalIndex[i] = GlobalNodeIndex(iel,i)-1u;
	  
      if(globalIndex[0] > globalIndex[1])
	o[iel][0]=-1;
  
      if(globalIndex[1] > globalIndex[2])
	o[iel][1]=-1;
  
      if(globalIndex[3] > globalIndex[2])
	o[iel][2]=-1;
  
      if(globalIndex[0] > globalIndex[3])
	o[iel][3]=-1;
      break;
    }
  case 4:{ //triangle
      nNodes=3;
      o[iel].resize(nNodes); 
      for(int i=0;i<nNodes;i++)
	o[iel][i]=1;
      
      int globalIndex[nNodes];
      
      for(unsigned i=0;i<nNodes;i++)  
	globalIndex[i] = GlobalNodeIndex(iel,i)-1u;
      
      if(globalIndex[0] > globalIndex[1])
	o[iel][0]=-1;
      if(globalIndex[1] > globalIndex[2])
	o[iel][1]=-1;
      if(globalIndex[2] > globalIndex[0])
	o[iel][2]=-1;
      break;
  }
  case 5: break;
  }
  }
}

/**
 * Return the total node number
 **/
unsigned elem::MeshNodeNum()const{
  return nvt;
}

/**
 * Returns the total number of vertices in the mesh
 **/
unsigned elem::MeshVertexNum()const{
  return nv0;
}

/**
 * Returns the total number of edges in the mesh
 **/
unsigned elem::MeshEdgeNum()const{
  return nv1;
}

/**
 * Returns the total number of quadrilateral faces
 **/
unsigned elem::MeshQuadFaceNum()const{
  return nv2;
}

/**
 * Returns the total number of triangle faces in the mesh 
 **/
unsigned elem::MeshTriFaceNum()const{
  return nv3;
}


/**
 * Set the total number of Lagrange nodes
 **/
void elem::SetNodeNumber(const unsigned& value){
  nvt=value;
}

/**
 * Set the number of total vertices
 **/
void elem::SetMeshNumOfVertices(const unsigned& value){
  nv0=value;
}

/**
 * Set the number of total edges
 **/
void elem::SetMeshNumOfEdges(const unsigned& value){
  nv1=value;
}

/**
 * Set the total number of quad faces
 **/
void elem::SetMeshQuadFaceNum(const unsigned& value){
  nv2=value;
}

/**
 * Set the total number of tri faces
 **/
void elem::SetMeshTriFaceNum(const unsigned& value){
  nv3=value;
}

/**
 * Return the total number of the element to refine 
 **/
unsigned elem::GetRefinedElementNumber(const char* name) const{
  if(!strcmp(name,"All")){
    return nelr;
  }
  unsigned i;
  i=GetIndex(name);
  return nelrt[i];
}
unsigned  elem::GetRefinedElementNumber(short unsigned ielt)const{
  return nelrt[ielt];
}

/**
 * Add value to the total number of the refined element
 **/
void elem::AddToRefinedElementNumber(const unsigned &value, const char name[]){
  if(!strcmp(name,"All")){
    nelr+=value;
    return;
  }
  unsigned i;
  i=this->GetIndex(name);
  nelrt[i]+=value;
}
void elem::AddToRefinedElementNumber(const unsigned &value, short unsigned ielt){
  nelrt[ielt]+=value;
}


unsigned elem::GetRefinedElementIndex(const unsigned &iel) const{
  return elr[iel];
}
void elem::SetRefinedElementIndex(const unsigned &iel, const unsigned &value){
  elr[iel]=value;
}
void elem::InitRefinedToZero(){
  nelr=nelrt[0]=nelrt[1]=nelrt[2]=nelrt[3]=nelrt[4]=nelrt[5]=0;
  for(unsigned i=0;i<nel;i++) elr[i]=0;
}

/**
 * Return the total number of the element 
 **/
unsigned elem::MeshElementNum(const char* name) const{
  if(!strcmp(name,"All")){
	return nel;
  }
  unsigned i;
  i=this->GetIndex(name);
  return nelt[i];
}

/**
 * Add value to the total number of the element
 **/
void elem::AddToElementNumber(const unsigned &value, const char name[]){
  unsigned i;
  i=GetIndex(name);
  nelt[i]+=value;
}

void elem::AddToElementNumber(const unsigned &value, short unsigned ielt){
  nelt[ielt]+=value;
}

unsigned elem::ElementFacetNum(const unsigned &iel, const unsigned &type)const{
  return NFC[ elt[iel] ][type];
}

unsigned elem::GetElementSquareFaceNumber(const unsigned &iel)const{
  return NFC[ elt[iel] ][0];
}
unsigned elem::GetElementTriangleFaceNumber(const unsigned &iel)const{
  return NFC[ elt[iel] ][1]-NFC[ elt[iel] ][0];
}

/**
 * Return the global adjacent-to-facet element index
 **/
int elem::FacetAdjacentElemIndex(const unsigned &iel,const unsigned &ifacet) const{
  return kel[iel][ifacet];
}

/**
 * Set the global adjacent-to-facet element index
 **/
void elem::SetFacetAdjacentElemIndex(const unsigned &iel,const unsigned &ifacet, const int &value) {
  kel[iel][ifacet]=value;
}


/**
 * Return element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
short unsigned elem::ElementType(const unsigned &iel) const{
  return elt[iel];
}

/**
 * Set element type: 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
void elem::SetElementType(const unsigned &iel, const short unsigned &value){
  elt[iel]=value;
}
  
/**
 * Set the memory storage and initialize nve and kvtel (node->element vectors)
 **/   
void elem::AllocateVertexElementMemory(){
  unsigned counter=(nelt[0]*NVE[0][0]+nelt[1]*NVE[1][0]+nelt[2]*NVE[2][0]+
		    nelt[3]*NVE[3][0]+nelt[4]*NVE[4][0]+nelt[5]*NVE[5][0]);
    
  kvtel=new unsigned* [nv0];
  nve = new unsigned[nv0];  
  
  for(unsigned inode=0;inode<nv0;inode++) // loop over all vertices
    nve[inode]=0;
  
  for(unsigned iel=0;iel<nel;iel++)
    for(unsigned inode=0;inode<ElementNodeNum(iel,0);inode++)
      nve[GlobalNodeIndex(iel,inode)-1u]++;

  kvtel_memory=new unsigned[counter];
  unsigned *pt= kvtel_memory;
  for(unsigned inode=0; inode<nv0; inode++){
    kvtel[inode]=pt;
    pt+=nve[inode];
  }
  memset(kvtel_memory, 0, counter*sizeof(unsigned));
}

/**
 * Return the number of elements which share the given node
 **/  
unsigned elem::NodeNumOfAdjacentElements(const unsigned& inode)const{
  return nve[inode];
}

/**
 * Return the element index for the given i-node in the j-position with 0<=j<nve(i)
 **/  
unsigned elem::NodeAdjacentElementIndex(const unsigned &inode,const unsigned &jnode)const{
  return kvtel[inode][jnode];
}

/**
 * Return the element address for the given i-node in the j-position with 0<=j<nve(i)
 **/  
const unsigned* elem::NodeAdjacentElementAddress(const unsigned& inode, const unsigned& jnode)const{
  return &kvtel[inode][jnode];
}

/**
 * Set the element index for the given i-node in the j-position with 0<=j<nve(i)
 **/  
void elem::SetVertexAdjacentElementIndex(const unsigned &inode,const unsigned &jnode, const unsigned &value){
  kvtel[inode][jnode]=value;
}

/**
 * return the index 0=hex, 1=Tet, 2=Wedge, 3=Quad, 4=Triangle and 5=Line
 **/
unsigned elem::GetIndex(const char name[]) const{
  unsigned index=0;
  if(!strcmp(name,"Hex")){
    index=0;
  }
  else if(!strcmp(name,"Tet")){
    index=1;
  }
  else if(!strcmp(name,"Wedge")){
    index=2;
  }
  else if(!strcmp(name,"Quad")){
    index=3;
  }
  else if(!strcmp(name,"Triangle")){
    index=4;
  }
  else if(!strcmp(name,"Line")){
    index=5;
  }
  else{
    cout<<"error! invalid Element Shape in elem::GetIndex(...)"<<endl;
    exit(0);
  }
  return index;
}

/* returns the number of degrees of freedom of the mesh for a given p. 
 The number of degree of freedom of the mesh = the total number of vertices in the mesh + 
 the total number of edges in the mesh x number of edge functions per edge + the total number of faces in the mesh x number of face functions per face 
 + number of elements x number of internal modes per element.
 */
const unsigned elem::NumOfDof(const int p) const{ //  linear: p=0.
  unsigned dof;
 
  dof = MeshVertexNum() + p*MeshEdgeNum() + (MeshQuadFaceNum())*(p>2)*(p-1)*(p-2)/2 + (MeshTriFaceNum())*(p>1)*(p-1)*p/2 
    + DIM[1]*MeshElementNum("Triangle")*(p>1)*p*(p-1)/2 + DIM[1]*MeshElementNum("Quad")*(p>2)*(p-1)*(p-2)/2 
    + DIM[2]*MeshElementNum("Hex")*(p>4)*(p-2)*(p-3)*(p-4)/6 + DIM[2]*MeshElementNum("Tet")*(p>2)*p*(p-1)*(p-2)/6 
    + DIM[2]*MeshElementNum("Wedge")*(p>3)*(p-1)*(p-2)*(p-3)/6;
  return dof; //1-based
}


/* returns the index scheme of triangles only for a 2D mesh */
int elem::TriIndex(const int iel) const
{
  int index=-1;
  for(int i=0;i<=iel;i++){
    if (ElementType(i)==4)
      index++;
  }
  
  if(index==-1)
    throw runtime_error("Element not a triangle in lsysPDE::TriIndex()"); 
  
  return index;
}

/* returns the index scheme of tetrahedra for a 3D mesh */
int elem::TetIndex(const int iel) const
{
  int index=-1;
  for(int i=0;i<=iel;i++){
    if (ElementType(i)==1)
      index++;
  }
  
  if(index==-1)
    throw runtime_error("Element not a tetrahedron in lsysPDE::TetIndex()"); 
  
  return index;
}

/* returns the index scheme of wedges for a 3D mesh */
int elem::WedgeIndex(const int iel) const
{
  int index=-1;
  for(int i=0;i<=iel;i++){
    if (ElementType(i)==2)
      index++;
  }
  
  if(index==-1)
    throw runtime_error("Element not a wedge in lsysPDE::WedgeIndex()"); 
  
  return index;
}

/* returns the index scheme of hexahedra for a 3D mesh */
int elem::HexIndex(const int iel) const
{
  int index=-1;
  for(int i=0;i<=iel;i++){
    if (ElementType(i)==0)
      index++;
  }
  
  if(index==-1)
    throw runtime_error("Element not a hex in lsysPDE::HexIndex()"); 
  
  return index;
}

/* The local to global dof mapping function an element 
 
 
The global dofs are numbered hierarchically. That means, all the vertices are added first. Then modes are added in the order p is increased.
Having global dofs for $p-1^{\text{st}}$ level assigned, the modes for pth level are added as follows.
First all the edge modes are added. Then all the face modes are added. If there are more than one face mode, first mode
is added to allfaces before adding the second one and so on. Also quadrialeral faces are added before triangle faces. Finally all the interior modes are added.


 */ 
void elem::SetGlobalDOF(const int P){ //0-based index P=maxOrder
  basish* pt[6];
    
  pt[0] = new hexh(P,0);
  pt[1] = new teth(P,0);
  pt[2] = new wedgeh(P,0);
  pt[3] = new quadh(P);
  pt[4] = new trih(P);
  pt[5] = new lineh(P);
  
  globalDOF.resize(MeshElementNum());
  
for(int iel=0;iel<MeshElementNum();iel++){
  unsigned ielt = ElementType(iel); 
  unsigned count=0;

  globalDOF[iel].resize(pt[ielt]->dimOfBasis());

  for(unsigned i=0; i < ElementNodeNum(iel,0); i++){ //vertex modes
    globalDOF[iel][count]=GlobalNodeIndex(iel,i)-1u;
    count++;
  }

  for(int i=1; i<P;i++){ 
    
    for(unsigned j = ElementNodeNum(iel,0); j < ElementNodeNum(iel,1); j++){ //edge modes
      globalDOF[iel][count] = NumOfDof(i-1) + GlobalNodeIndex(iel,j)-1u - MeshVertexNum();
      count++;
    }

    if (DIMENSION == 0) continue;
    
    if((i>=3)&&(ielt>=3)) //interior modes for tri and quad
      for(unsigned j=0; j<i-2; j++){
	globalDOF[iel][count] = NumOfDof(i-1) + MeshEdgeNum() + MeshElementNum()*j + iel;
  	count++;
      }
    
    if((i>=2)&&(ielt==4)){ //additional tri interior mode (if element=triangle)
     globalDOF[iel][count] = NumOfDof(i-1) + MeshEdgeNum() + MeshElementNum()*(i-2) + TriIndex(iel);
      count++;
    }
    
    if (DIMENSION == 1) continue;
    
     if(i>=3) //face modes for tri and quad
      for(unsigned j=0; j<i-2; j++){
	for(unsigned k =  ElementNodeNum(iel,1); k < ElementNodeNum(iel,3); k++){
	globalDOF[iel][count] = NumOfDof(i-1) + GlobalNodeIndex(iel,k)-1u + j*(MeshQuadFaceNum() + MeshTriFaceNum()) - MeshVertexNum();
  	count++;
        }
      }
      
    if(i>=2) //additional tri face mode 
      for(unsigned j =  ElementNodeNum(iel,2); j < ElementNodeNum(iel,3); j++){
	globalDOF[iel][count] = NumOfDof(i-1) + GlobalNodeIndex(iel,j)-1u - MeshVertexNum() + (i-3)*MeshQuadFaceNum()+(i-2)*(MeshTriFaceNum());
	count++;
    }
    
    // Add interior modes 
    if((i>=3)&&(ielt==1)) // tet 
      for(int j=0;j<(i-2)*(i-1)/2;j++){
	globalDOF[iel][count] = NumOfDof(i-1) + MeshEdgeNum() + (i-2)*MeshQuadFaceNum() + (i-1)*MeshTriFaceNum() + j*MeshElementNum("Tet")+TetIndex(iel);
	count++;
    }
    if((i>=4)&&(ielt==2)) // wedge do wedge index or wedge + tets index for hybrid meshes and replace iel
      for(unsigned j=0;j<(i-3)*(i-2)/2;j++){
	globalDOF[iel][count] = NumOfDof(i-1) + MeshEdgeNum() + (i-2)*MeshQuadFaceNum() + (i-1)*MeshTriFaceNum() + (i-2)*(i-1)/2*MeshElementNum("Tet")+j*MeshElementNum("Wedge") + WedgeIndex(iel);
	count++;
    }
    if((i>=5)&&(ielt==0)) // hex
      for(unsigned j=0;j<(i-3)*(i-4)/2;j++){
	globalDOF[iel][count] = NumOfDof(i-1) + MeshEdgeNum() + (i-2)*MeshQuadFaceNum() + (i-1)*MeshTriFaceNum() + (i-2)*(i-1)/2*MeshElementNum("Tet") + (i-3)*(i-2)/2*MeshElementNum("Wedge") + MeshElementNum("Hex")*j + HexIndex(iel);
	count++;
    }
  }
}
//  // view results
//  cout<<"Global dofs"<<endl;
//  for(unsigned i=0;i<globalDOF.size();i++){
//    for(unsigned j=0;j<globalDOF[i].size();j++)
// 	cout<<globalDOF[i][j]<<" ";
//     cout<<endl;
//    }
 
  for(int p=0;p<6;p++)
    delete pt[p];

}
