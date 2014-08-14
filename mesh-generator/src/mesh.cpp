#include "main.hpp"
#include "mesh.hpp"

/**  
 *  This constructor generates the coarsest mesh level, $l_0$, from the gambit data file  
 **/ 
mesh::mesh(const char infile[], vector < vector < double> > &vt){
  grid=0;
 
  if(DIM[0]==1)
    Read1D(infile,vt);
  else
    ReadGambit(infile,vt);
};

/**  
 *  This constructor generates a fine mesh level, $l_i$, from a coarse mesh level $l_{i-1}$, $i>0$ 
 **/ 
mesh::mesh(const unsigned& igrid, elem* elc){
  const unsigned vertex_index[6][8][8]={{{1,9,25,12,17,21,27,24},{9,2,10,25,21,18,22,27},
					 {25,10,3,11,27,22,19,23},{12,25,11,4,24,27,23,20},
					 {17,21,27,24,5,13,26,16},{21,18,22,27,13,6,14,26},
					 {27,22,19,23,26,14,7,15},{24,27,23,20,16,26,15,8}},
					{{1,5,7,8},{2,6,5,9},{3,7,6,10},{8,9,10,4},
					 {5,6,8,9},{6,10,8,9},{5,8,6,7},{6,8,10,7}},
					{{1,7,9,13,16,18},{2,8,7,14,17,16},
					 {3,9,8,15,18,17},{7,8,9,16,17,18},
					 {13,16,18,4,10,12},{14,17,16,5,11,10},
					 {15,18,17,6,12,11},{16,17,18,10,11,12}},
					{{1,5,9,8},{5,2,6,9},{9,6,3,7},{8,9,7,4}},
					{{1,4,6},{2,5,4},{3,6,5},{4,5,6}},
					{{1,3},{3,2}}};
  
  const unsigned face_index[6][6][4][2]={{{{0,0},{1,0},{4,0},{5,0}},
					  {{1,1},{2,1},{5,1},{6,1}},
					  {{2,2},{3,2},{6,2},{7,2}},
					  {{3,3},{0,3},{7,3},{4,3}},
					  {{0,4},{1,4},{2,4},{3,4}},
					  {{4,5},{5,5},{6,5},{7,5}}},
					 {{{0,0},{1,0},{2,0},{6,3}},
					  {{0,1},{1,3},{3,1},{4,3}},
					  {{1,1},{2,3},{3,2},{5,1}},
					  {{2,1},{0,3},{3,3},{7,2}}},
					 {{{0,0},{1,2},{4,0},{5,2}},
					  {{1,0},{2,2},{5,0},{6,2}},
					  {{2,0},{0,2},{6,0},{4,2}},
					  {{0,3},{1,3},{2,3},{3,3}},
					  {{4,4},{5,4},{6,4},{7,4}}},
					 {{{0,0},{1,0}},
					  {{1,1},{2,1}},
					  {{2,2},{3,2}},
					  {{3,3},{0,3}}},
					 {{{0,0},{1,2}},
					  {{1,0},{2,2}},
					  {{2,0},{0,2}}},
					 {{{0,0}},
					  {{1,1}}} };
  
  const unsigned midpoint_index[6][12][2]={{{0,1},{1,2},{2,3},{3,0},
					    {4,5},{5,6},{6,7},{7,4},
					    {0,4},{1,5},{2,6},{3,7}},
					   {{0,1},{1,2},{2,0},
					    {0,3},{1,3},{2,3}},
					   {{0,1},{1,2},{2,0},
					    {3,4},{4,5},{5,3},
					    {0,3},{1,4},{2,5}},
					   {{0,1},{1,2},{2,3},{3,0}},
					   {{0,1},{1,2},{2,0}},
					   {{0,1}}};

  const unsigned midpoint_index2[6][7][8]={{{0,8,0,11,16},
					    {0,0,9,0,0,17},
					    {0,0,0,10,0,0,18},
					    {0,0,0,0,0,0,0,19},
					    {0,0,0,0,0,12,0,15},
					    {0,0,0,0,0,0,13},
					    {0,0,0,0,0,0,0,14}},
					   {{0,4,6,7},
					    {0,0,5,8},
					    {0,0,0,9}},
					   {{0,6,8,12},
					    {0,0,7,0,13},
					    {0,0,0,0,0,14},
					    {0,0,0,0,9,11},
					    {0,0,0,0,0,10}},
					   {{0,4,0,7},
					    {0,0,5},
					    {0,0,0,6}},
					   {{0,3,5},
					    {0,0,4}},
					   {{0,2}}};
  grid=igrid;
 
  nel=elc->GetRefinedElementNumber()*REF_INDEX;
  el=new elem(elc);
  
  
  unsigned jel=0;
  
  //divide each coarse element in 8(3D), 4(2D) or 2(1D) fine elemets and find all the vertices
  for(unsigned iel=0; iel<elc->MeshElementNum(); iel++){
    if( elc->GetRefinedElementIndex(iel) ){
      elc->SetRefinedElementIndex(iel,jel+1u);
      unsigned elt=elc->ElementType(iel);
      // project element type
      for(unsigned j=0;j<REF_INDEX;j++){
	el->SetElementType(jel+j,elt);
      }
      // project vertex indices
      for(unsigned j=0;j<REF_INDEX;j++)
	for(unsigned inode=0;inode<elc->ElementNodeNum(iel,0);inode++)
	  el->SetElementNodeGlobalIndex(jel+j,inode,elc->GlobalNodeIndex(iel,vertex_index[elt][j][inode]-1u));
      // project face indeces
      for(unsigned iface=0;iface<elc->ElementFacetNum(iel);iface++){
	int value=elc->FacetAdjacentElemIndex(iel,iface);
	if(0>value)
	  for(unsigned jface=0;jface<FACE_INDEX;jface++)
	    el->SetFacetAdjacentElemIndex(jel+face_index[elt][iface][jface][0],face_index[elt][iface][jface][1], value);
      }
      // update element numbers
      jel+=REF_INDEX; 
      el->AddToElementNumber(REF_INDEX,elt);
    }
  }
  
  nvt=elc->MeshNodeNum ();
  el->SetMeshNumOfVertices(nvt);
  
  //find all the elements near each vertex 
  BuildAdjVtx();
  
  //initialize to zero all the middle edge points
  for(unsigned iel=0;iel<nel;iel++) 
    for(unsigned inode=el->ElementNodeNum(iel,0);inode<el->ElementNodeNum(iel,1);inode++)
      el->SetElementNodeGlobalIndex(iel,inode,0);
    
  //find all the middle edge points
  for(unsigned iel=0;iel<nel;iel++){ 
    unsigned ielt=el->ElementType(iel);
    unsigned istart=el->ElementNodeNum(iel,0);
    unsigned iend=el->ElementNodeNum(iel,1);
    for(unsigned inode=istart; inode<iend; inode++)
      if(0==el->GlobalNodeIndex(iel,inode)){
	nvt++;
	el->SetElementNodeGlobalIndex(iel,inode,nvt);
	unsigned im=el->GlobalNodeIndex(iel,midpoint_index[ielt][inode-istart][0]);
	unsigned ip=el->GlobalNodeIndex(iel,midpoint_index[ielt][inode-istart][1]);
	
	//find all the near elements which share the same middle edge point
	for(unsigned j=0;j<el->NodeNumOfAdjacentElements(im-1u);j++){
	  unsigned jel=el->NodeAdjacentElementIndex(im-1u,j)-1u;
	  if(jel>iel){
	    unsigned jm=0,jp=0;
	    unsigned jelt=el->ElementType(jel);
	    for(unsigned jnode=0;jnode<el->ElementNodeNum(jel,0);jnode++){
	      if(el->GlobalNodeIndex(jel,jnode)==im) {
		jm=jnode+1u;
		break;
	      }
	    }
	    if(jm!=0){ 
	      for(unsigned jnode=0;jnode<el->ElementNodeNum(jel,0);jnode++){
		if(el->GlobalNodeIndex(jel,jnode)==ip){
		  jp=jnode+1u;
		  break;
		}
	      }
	      if(jp!=0){
		if(jp<jm){ 
		  unsigned tp=jp; jp=jm; jm=tp; 
		}
		el->SetElementNodeGlobalIndex(jel,midpoint_index2[jelt][--jm][--jp],nvt);	
	      }
	    }
	  }
	}
      }
  }
  el->SetMeshNumOfEdges(nvt - el->MeshVertexNum());
  // Now build kmid 
  Buildkmid();
  Buildkel();

}

/**
 * This function read the data from the custom 1D file.
 * It is used in the constructor of the coarse mesh.  
 **/
void mesh::Read1D(const char infile [], vector < vector < double> > &vt){
  // set the file
  std::ifstream inf;
  std::string str2;
 
  grid=0;
  /* read NUMBER OF ELEMENTS */ 
  inf.open(infile); 
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " can not read parameters\n";
    exit(0);
  }
  str2="0"; 
  while (str2.compare("NEL") != 0) inf >> str2;
  inf >> nel;
  inf.close();
  
  nvt=2*nel+1;
  el = new elem(nel);
  el->AddToElementNumber(nel,"Line");
  for(unsigned iel=0;iel<nel;iel++){    
    el->SetElementType(iel,5);
    el->SetElementNodeGlobalIndex(iel,0,iel+1); 
    el->SetElementNodeGlobalIndex(iel,1,iel+2); 
    el->SetElementNodeGlobalIndex(iel,2,(nel+1)+(iel+1)); 
  }
  // end read NUMBER OF ELEMENTS **************** A

  // read NODAL COORDINATES **************** B
  inf.open(infile); if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read nodes\n";
    exit(0);
  }
  double x0,x1,h;
  while (str2.compare("COORDINATES") != 0) inf >> str2;   
  inf >> x0>>x1; 
  inf.close(); 
  
  if(x1<x0){
    double temp=x1; x1=x0; x0=temp;
  }
  vt[0].resize(nvt);
  vt[1].resize(nvt);
  vt[2].resize(nvt);
  
  h=(x1-x0)/nel;
  for(unsigned iel=0;iel<nel;iel++){   
    vt[0][iel]=x0+iel*h;  
    vt[1][iel]=0.;  
    vt[2][iel]=0.;
    vt[0][nel+(iel+1)]=vt[0][iel]+h/2.;
    vt[1][nel+(iel+1)]=0.;  
    vt[2][nel+(iel+1)]=0.;
  }
  vt[0][nel]=x1;
  vt[1][nel]=0.;  
  vt[2][nel]=0.;
  // end read NODAL COORDINATES ************* B
 
  // set boundary **************** C
  el->SetFacetAdjacentElemIndex(0,0,-2);
  el->SetFacetAdjacentElemIndex(nel-1,1,-3);
  // end set boundary **************** 
  el->SetNodeNumber(nvt);
  el->SetMeshNumOfVertices(nel+1);
  el->SetMeshNumOfEdges(nel);
  el->SetMeshQuadFaceNum(0);
}



/**
 * This function read the data form the gambit file.
 * It is used in the constructor of the coarse mesh.  
 **/
void mesh::ReadGambit(const char infile [], vector < vector < double> > &vt){
 
 /* Map gambit element nodes into FEM geometry */ 
  const unsigned GambitNodeIndex[6][27] = {{4,16,0,15,23,11,7,19,3,12,20,8,25,26,24,14,22,10,5,17,1,13,21,9,6,18,2}, //hex
					   {0,4,1,6,5,2,7,8,9,3,10,11,12,13,14}, //tet
					   {3,11,5,9,10,4,12,17,14,15,16,13,0,8,2,6,7,1,18,19,20}, //wedge
					   {0,4,1,5,2,6,3,7,8}, //quad
					   {0,3,1,4,2,5,6}, //tri
					   {0,2,1}}; //line
 /* facets */
  const unsigned GambitFacetIndex[6][6]={{0,4,2,5,3,1}, // hex 
					{0,1,2,3}, // tet
					{2,1,0,4,3}, // wedge
					{0,1,2,3}, // quad
					{0,1,2}, // tri
					{0,1}}; // line

  std::ifstream inf;
  std::string str2;
  unsigned ngroup;
  unsigned nbcd;

  grid=0;
 
  /* read control data */
  inf.open(infile); 
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " can not read parameters\n";
    exit(0);
  }
  str2="0"; 
  while (str2.compare("NDFVL") != 0) inf >> str2; 
  inf >> nvt >> nel >>  ngroup >> nbcd >>str2 >> str2 ; 
  inf >> str2;
  if(str2.compare("ENDOFSECTION") != 0){
    cout<<"error control data mesh"<<endl;
    exit(0);
  }
  inf.close();
  /* end read control data */
  
  /* read ELEMENTS/CELLS */
  inf.open(infile); if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read elements\n";
    exit(0);
  }
  
  el = new elem(nel);
  
  while (str2.compare("ELEMENTS/CELLS") != 0)
    inf >> str2; 
  inf >> str2; 
  
  for(unsigned iel=0;iel<nel;iel++){
    unsigned nve;
    inf >> str2 >> str2 >> nve;
    if(nve==27){
      el->AddToElementNumber(1,"Hex");
      el->SetElementType(iel,0);
    }  
    else if(nve==10){
      el->AddToElementNumber(1,"Tet");   
      el->SetElementType(iel,1);
    } 
    else if(nve==18){
      el->AddToElementNumber(1,"Wedge");
      el->SetElementType(iel,2);
    }
    else if(nve==9){
      el->AddToElementNumber(1,"Quad");
      el->SetElementType(iel,3);
    }
    else if(nve==6 && DIM[1]==1){
      el->AddToElementNumber(1,"Triangle");
      el->SetElementType(iel,4);
    }
    else if(nve==3 && DIM[0]==1){
      el->AddToElementNumber(1,"Line");
      el->SetElementType(iel,5);
    }
    else{
      cout<<"Error! Invalid element type in reading Gambit File!"<<endl;
      cout<<"Error! Use a second order discretization"<<endl;
      exit(0);
    }
    
    /* Set global node indices in kvert */
    for(unsigned i=0;i<nve;i++){ //nve: number of nodes in the gambit element
      unsigned inode=GambitNodeIndex[el->ElementType(iel)][i];
      double value;
      inf>>value; // node in the gambit file} 
      el->SetElementNodeGlobalIndex(iel,inode,value);
    }
  }
  
  inf >> str2;
  if(str2.compare("ENDOFSECTION") != 0){
    cout<<"error element data mesh"<<endl;
    exit(0);
  }
  inf.close();
  
  /* end read  ELEMENTS/CELLS */
  
  unsigned gambitNvt = nvt;
  
  // read NODAL COORDINATES **************** C
  inf.open(infile); if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read nodes\n";
    exit(0);
  }
  while (str2.compare("COORDINATES") != 0)
    inf >> str2;   
  
  inf >> str2;  // 2.0.4
  vt[0].resize(nvt);
  vt[1].resize(nvt);
  vt[2].resize(nvt);
  
  if(DIM[2]==1){
    for(unsigned j=0;j<gambitNvt;j++) 
      inf >> str2 >> vt[0][j]>>vt[1][j]>>vt[2][j];
  }
  else if(DIM[1]==1){
    for(unsigned j=0;j<gambitNvt;j++){ 
      inf >> str2 >> vt[0][j]>>vt[1][j];
      vt[2][j]=0.;
    }
  }
  else if(DIM[0]==1){
    for(unsigned j=0;j<gambitNvt;j++){ 
      inf >> str2 >> vt[0][j];
      vt[1][j]=0.;
      vt[2][j]=0.;
    }
  }
  inf >> str2; // "ENDOFSECTION"
  if(str2.compare("ENDOFSECTION") != 0){
    cout<<"error node data mesh 1"<<endl;
    exit(0);
  }
  inf.close(); 
  /* end read NODAL COORDINATES */
  
  /* read boundary */  
  inf.open(infile); 
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read boundary\n";
    exit(0);
  }
  for(unsigned k=0; k<nbcd; k++){
    while (str2.compare("CONDITIONS") != 0) inf >> str2; 
    inf >> str2;
    int value;
    unsigned nface;
    inf >>value>>str2>>nface>>str2>>str2;
    value=-value-1;
    for(unsigned i=0;i<nface;i++){
      unsigned iel,iface;
      inf>>iel>>str2>>iface;
      iel--;
      iface=GambitFacetIndex[el->ElementType(iel)][iface-1u];
      el->SetFacetAdjacentElemIndex(iel,iface,value);
    }
    inf >> str2;
    if(str2.compare("ENDOFSECTION") != 0){
      cout<<"error boundary data mesh"<<endl;
      exit(0);
    }
  }
  inf.close();
  /* end read boundary */

  /* Reorder existing nodes to comply with FEM topology */
  unsigned short gambitNodeNum[6][4]={{8,20,26,27},//hex
				      {4,10,10,10}, //tet
				      {6,15,18,18}, //wedge
				      {4,8,8,9}, //quad
				      {3,6,6,6}, //tri
				      {2,3,3,3}}; //line
  
  vector <unsigned> dof_index;
  dof_index.resize(nvt);
  for(unsigned i=0;i<nvt;i++){
    dof_index[i]=i+1;
  }
  //reorder vertices and mid-points vs central points
  for (unsigned iel=0; iel<nel; iel++) {
    unsigned ielt = el->ElementType(iel);
    for (unsigned inode=0; inode<gambitNodeNum[ielt][2]; inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
	unsigned jelt = el->ElementType(jel);
    for (unsigned jnode=gambitNodeNum[jelt][2];jnode<gambitNodeNum[jelt][3]; jnode++) {
      unsigned ii=el->GlobalNodeIndex(iel,inode)-1;
      unsigned jj=el->GlobalNodeIndex(jel,jnode)-1;
      unsigned i0=dof_index[ii];
      unsigned i1=dof_index[jj];
      if(i0>i1){
        dof_index[ii]=i1;
        dof_index[jj]=i0;
// 	double x0 = vt[0][i0-1u]; vt[0][i0-1u]=vt[0][i1-1u]; vt[0][i1-1u]=x0; 
// 	double y0 = vt[1][i0-1u]; vt[1][i0-1u]=vt[1][i1-1u]; vt[1][i1-1u]=y0;
// 	double z0 = vt[2][i0-1u]; vt[2][i0-1u]=vt[2][i1-1u]; vt[2][i1-1u]=z0;
      }
    }
      }
    }
  }
  
/*   //( vertices + mid points + face points ) before interiors
  for(unsigned iel=0;iel<nel;iel++){
    unsigned ielt = el->ElementType(iel);
    for(unsigned inode=0; inode<gambitNodeNum[ielt][2]; inode++){
      for(unsigned jel=0;jel<nel;jel++){
	unsigned jelt = el->ElementType(jel);
	for(unsigned jnode=gambitNodeNum[jelt][2];jnode<gambitNodeNum[jelt][3];jnode++){
	  unsigned i0=el->GlobalNodeIndex(iel,inode);
	  unsigned i1=el->GlobalNodeIndex(jel,jnode);
	  if(i0>i1){
	    double x0 = vt[0][i0-1u]; vt[0][i0-1u]=vt[0][i1-1u]; vt[0][i1-1u]=x0; 
	    double y0 = vt[1][i0-1u]; vt[1][i0-1u]=vt[1][i1-1u]; vt[1][i1-1u]=y0;
	    double z0 = vt[2][i0-1u]; vt[2][i0-1u]=vt[2][i1-1u]; vt[2][i1-1u]=z0;
	    for(unsigned kel=0;kel<nel;kel++){
	      unsigned kelt= el->ElementType(kel);
	      for(unsigned knode=0;knode<gambitNodeNum[kelt][3];knode++){
		unsigned i2=el->GlobalNodeIndex(kel,knode);
		if(i2==i0) el->SetElementNodeGlobalIndex(kel,knode,i1);
		else if(i2==i1) el->SetElementNodeGlobalIndex(kel,knode,i0);
	      }
	    }
	  }
	}
      }
    }
  }
	*/	

  
  for (unsigned iel=0; iel<nel; iel++) {
    unsigned ielt = el->ElementType(iel);
    for (unsigned inode=0; inode<gambitNodeNum[ielt][1]; inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
	unsigned jelt = el->ElementType(jel);
    for (unsigned jnode=gambitNodeNum[jelt][1];jnode<gambitNodeNum[jelt][2]; jnode++) {
      unsigned ii=el->GlobalNodeIndex(iel,inode)-1;
      unsigned jj=el->GlobalNodeIndex(jel,jnode)-1;
      unsigned i0=dof_index[ii];
      unsigned i1=dof_index[jj];
      if(i0>i1){
        dof_index[ii]=i1;
        dof_index[jj]=i0;
// 	double x0 = vt[0][i0-1u]; vt[0][i0-1u]=vt[0][i1-1u]; vt[0][i1-1u]=x0; 
// 	double y0 = vt[1][i0-1u]; vt[1][i0-1u]=vt[1][i1-1u]; vt[1][i1-1u]=y0;
// 	double z0 = vt[2][i0-1u]; vt[2][i0-1u]=vt[2][i1-1u]; vt[2][i1-1u]=z0;
      }
    }
      }
    }
  }
/*  
 //( vertices + mid points ) before face points
  for(unsigned iel=0;iel<nel;iel++){
    unsigned ielt = el->ElementType(iel);
    for(unsigned inode=0; inode<gambitNodeNum[ielt][1]; inode++){
      for(unsigned jel=0;jel<nel;jel++){
	unsigned jelt = el->ElementType(jel);
	for(unsigned jnode=gambitNodeNum[jelt][1];jnode<gambitNodeNum[jelt][2];jnode++){
	  unsigned i0=el->GlobalNodeIndex(iel,inode);
	  unsigned i1=el->GlobalNodeIndex(jel,jnode);
	  if(i0>i1){
	    double x0 = vt[0][i0-1u]; vt[0][i0-1u]=vt[0][i1-1u]; vt[0][i1-1u]=x0; 
	    double y0 = vt[1][i0-1u]; vt[1][i0-1u]=vt[1][i1-1u]; vt[1][i1-1u]=y0;
	    double z0 = vt[2][i0-1u]; vt[2][i0-1u]=vt[2][i1-1u]; vt[2][i1-1u]=z0;
	    for(unsigned kel=0;kel<nel;kel++){
	      unsigned kelt= el->ElementType(kel);
	      for(unsigned knode=0;knode<gambitNodeNum[kelt][2];knode++){
		unsigned i2=el->GlobalNodeIndex(kel,knode);
		if(i2==i0) el->SetElementNodeGlobalIndex(kel,knode,i1);
		else if(i2==i1) el->SetElementNodeGlobalIndex(kel,knode,i0);
	      }
	    }
	  }
	}
      }
    }
  }
  */

  
  for (unsigned iel=0; iel<nel; iel++) {
    
    for (unsigned inode=0; inode<el->ElementNodeNum(iel,0); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
	
    for (unsigned jnode=el->ElementNodeNum(jel,0);jnode<el->ElementNodeNum(jel,1); jnode++) {
      unsigned ii=el->GlobalNodeIndex(iel,inode)-1;
      unsigned jj=el->GlobalNodeIndex(jel,jnode)-1;
      unsigned i0=dof_index[ii];
      unsigned i1=dof_index[jj];
      if(i0>i1){
        dof_index[ii]=i1;
        dof_index[jj]=i0;
// 	double x0 = vt[0][i0-1u]; vt[0][i0-1u]=vt[0][i1-1u]; vt[0][i1-1u]=x0; 
// 	double y0 = vt[1][i0-1u]; vt[1][i0-1u]=vt[1][i1-1u]; vt[1][i1-1u]=y0;
// 	double z0 = vt[2][i0-1u]; vt[2][i0-1u]=vt[2][i1-1u]; vt[2][i1-1u]=z0;
      }
    }
      }
    }
  }
  
//   // vertices before mid points
//   for(unsigned iel=0;iel<nel;iel++){
//     //cout<<iel<<endl;
//     for(unsigned inode=0;inode<el->ElementNodeNum(iel,0);inode++){
//       for(unsigned jel=0;jel<nel;jel++){
// 	for(unsigned jnode=el->ElementNodeNum(jel,0);jnode<el->ElementNodeNum(jel,1);jnode++){
// 	  unsigned i0=el->GlobalNodeIndex(iel,inode);
// 	  unsigned i1=el->GlobalNodeIndex(jel,jnode);
// 	  if(i0>i1){
// 	    double x0=vt[0][i0-1u];  vt[0][i0-1u]=vt[0][i1-1u]; vt[0][i1-1u]=x0; 
// 	    double y0=vt[1][i0-1u];  vt[1][i0-1u]=vt[1][i1-1u]; vt[1][i1-1u]=y0;
// 	    double z0=vt[2][i0-1u];  vt[2][i0-1u]=vt[2][i1-1u]; vt[2][i1-1u]=z0;
// 	    for(unsigned kel=0;kel<nel;kel++){
// 	      for(unsigned knode=0;knode<el->ElementNodeNum(kel,1);knode++){
// 		unsigned i2=el->GlobalNodeIndex(kel,knode);
// 		if(i2==i0) el->SetElementNodeGlobalIndex(kel,knode,i1);
// 		else if(i2==i1) el->SetElementNodeGlobalIndex(kel,knode,i0);
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
  // update all
  for (unsigned iel=0; iel<nel; iel++) {
    unsigned ielt= el->ElementType(iel);
    for (unsigned inode=0; inode<gambitNodeNum[ielt][3]; inode++) {
      unsigned ii=el->GlobalNodeIndex(iel,inode)-1;
      el->SetElementNodeGlobalIndex(iel,inode,dof_index[ii]);
    }
  }
  vector <double> vt_temp;
  for(int i=0;i<3;i++){
    vt_temp=vt[i];
    for(unsigned j=0;j<nvt;j++){
      vt[i][dof_index[j]-1]=vt_temp[j];
    }
  }
  
  el->SetNodeNumber(nvt);
  
  unsigned nv0=0;
  for(unsigned iel=0;iel<nel;iel++)
    for(unsigned inode=0;inode<el->ElementNodeNum(iel,0);inode++){
      unsigned i0=el->GlobalNodeIndex(iel,inode);
      if(nv0<i0) nv0=i0;
    }
    
  el->SetMeshNumOfVertices(nv0);
  el->SetOrientationFlagsAndType();
  BuildAdjVtx();
  Buildkel();
  AddAdditionalNodes();
 
  // Add and initialize additional nodes
  vt[0].resize(nvt);
  vt[1].resize(nvt);
  vt[2].resize(nvt);
  
   for(unsigned j=gambitNvt;j<nvt;j++){ 
      vt[0][j]=0.;
      vt[1][j]=0.;
      vt[2][j]=0.;
    }
    
  // Triangular faces before hex interior
  for(unsigned iel=0;iel<nel;iel++){
    for(unsigned inode=inode<el->ElementNodeNum(iel,2); inode<el->ElementNodeNum(iel,3); inode++){
      for(unsigned jel=0;jel<nel;jel++){
	if(el->ElementType(jel)!=0)continue;
	for(unsigned jnode=el->ElementNodeNum(jel,3);jnode<el->ElementNodeNum(jel,4);jnode++){
	  unsigned i0=el->GlobalNodeIndex(iel,inode);
	  unsigned i1=el->GlobalNodeIndex(jel,jnode);
	  if(i0>i1){
	    double x0 = vt[0][i0-1u]; vt[0][i0-1u]=vt[0][i1-1u]; vt[0][i1-1u]=x0; 
	    double y0 = vt[1][i0-1u]; vt[1][i0-1u]=vt[1][i1-1u]; vt[1][i1-1u]=y0;
	    double z0 = vt[2][i0-1u]; vt[2][i0-1u]=vt[2][i1-1u]; vt[2][i1-1u]=z0;
	    for(unsigned kel=0;kel<nel;kel++){
	      for(unsigned knode=el->ElementNodeNum(kel,2);knode<el->ElementNodeNum(kel,4);knode++){
		unsigned i2=el->GlobalNodeIndex(kel,knode);
		if(i2==i0) el->SetElementNodeGlobalIndex(kel,knode,i1);
		else if(i2==i1) el->SetElementNodeGlobalIndex(kel,knode,i0);
	      }
	    }
	  }
	}
      }
    }
  }
  
  SetMeshNodeInformation();  
 
};

void mesh::AddAdditionalNodes() //tested!
{
  /* face node for triangular faces in tetrahedron and wedge */
  for(unsigned iel=0;iel<el->MeshElementNum();iel++)
    for(unsigned inode=el->ElementNodeNum(iel,2);inode<el->ElementNodeNum(iel,3);inode++) // triangular faces
      el->SetElementNodeGlobalIndex(iel,inode,0);  
    
  if (DIMENSION == 2){
  
  for(unsigned iel=0;iel<el->MeshElementNum();iel++){
    for(unsigned ifacet=0; ifacet<el->ElementFacetNum(iel,1)-el->ElementFacetNum(iel,0); ifacet++){ // loop over triangular facets
      unsigned inode=el->ElementNodeNum(iel,2)+ifacet;
      unsigned numOfQuadFacets = el->ElementNodeNum(iel,2)-el->ElementNodeNum(iel,1);
      if( 0==el->GlobalNodeIndex(iel,inode) ) {
	el->SetElementNodeGlobalIndex(iel,inode,++nvt);
	
	unsigned i1=el->FacetNodeGlobalIndex(iel,numOfQuadFacets+ifacet,0);
	unsigned i2=el->FacetNodeGlobalIndex(iel,numOfQuadFacets+ifacet,1);
	unsigned i3=el->FacetNodeGlobalIndex(iel,numOfQuadFacets+ifacet,2);
	
	for(unsigned j=0; j < el->NodeNumOfAdjacentElements(i1-1u); j++){ 
	  
	  unsigned jel= el->NodeAdjacentElementIndex(i1-1u,j)-1u;
	  unsigned JnumOfQuadFacets = el->ElementNodeNum(jel,2)-el->ElementNodeNum(jel,1);
	  if(jel>iel){
	    for(unsigned jfacet=0;jfacet<el->ElementFacetNum(jel,1)-el->ElementFacetNum(jel,0);jfacet++){
	      unsigned jnode=el->ElementNodeNum(jel,2)+jfacet;
	      if( 0==el->GlobalNodeIndex(jel,jnode) ) {
	      	unsigned j1=el->FacetNodeGlobalIndex(jel,JnumOfQuadFacets+jfacet,0); 
		unsigned j2=el->FacetNodeGlobalIndex(jel,JnumOfQuadFacets+jfacet,1);
		unsigned j3=el->FacetNodeGlobalIndex(jel,JnumOfQuadFacets+jfacet,2);
		
		if((i1==j1 || i1==j2 || i1==j3  )&&
		   (i2==j1 || i2==j2 || i2==j3  )&&
		   (i3==j1 || i3==j2 || i3==j3  )){
		  el->SetElementNodeGlobalIndex(jel,jnode,nvt);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  }
  
  /* Add interior node for triangle, wedge and tetrahedron. */ 
  for(unsigned iel=0;iel<el->MeshElementNum();iel++)
    if( (el->ElementType(iel)==1) || (el->ElementType(iel)==2) || (el->ElementType(iel)==4))
      el->SetElementNodeGlobalIndex(iel,el->ElementNodeNum(iel,4)-1u,++nvt);
  
}

/* Set the number of total nodes, vertices, edges, and faces in the mesh */ 
void mesh::SetMeshNodeInformation() //tested!
{
  unsigned nv0=el->MeshVertexNum();
  
  int nv1=0;
  for(unsigned iel=0;iel<nel;iel++)
    for(unsigned inode=el->ElementNodeNum(iel,0);inode<el->ElementNodeNum(iel,1);inode++){
      int i1=el->GlobalNodeIndex(iel,inode);
      if(nv1<i1) nv1=i1;
    }
  el->SetMeshNumOfEdges(nv1-nv0);
  
  int nv2=nv1;
  for(unsigned iel=0;iel<nel;iel++)
    for(unsigned inode=el->ElementNodeNum(iel,1);inode<el->ElementNodeNum(iel,2);inode++){
      int i2=el->GlobalNodeIndex(iel,inode);
      if(nv2<i2) nv2=i2;
    }
    
  el->SetMeshQuadFaceNum(nv2-nv1);
  el->SetMeshTriFaceNum(nvt-nel-nv2);   
  el->SetNodeNumber(nvt); 
  
//   cout<<"Mesh information"<<endl;
//   cout<<"Elements: "<<el->MeshElementNum()<<"\t Vertices: "<<el->MeshVertexNum()<<"\t Edges: "<<el->MeshEdgeNum()<<endl;
//   cout<<"Quad faces: "<<el->MeshQuadFaceNum()<<"\t Tri faces: "<<el->MeshTriFaceNum()<<"\t Total nodes: "<<el->MeshNodeNum()<<endl;
  
 }

/**
 * This function searches all the elements around all the vertices 
 **/
void mesh::BuildAdjVtx(){
  el->AllocateVertexElementMemory();
  for(unsigned iel=0;iel<nel;iel++)
    for(unsigned inode=0;inode < el->ElementNodeNum(iel,0);inode++){
      unsigned ii=el->GlobalNodeIndex(iel,inode)-1u;
      unsigned jj=0;
      while( 0 != el->NodeAdjacentElementIndex(ii,jj) ) jj++;
      el->SetVertexAdjacentElementIndex(ii,jj,iel+1u);
    }
}


/**
 * This function generates kmid for hex and wedge elements 
 **/
void mesh::Buildkmid(){ 
  
  for(unsigned iel=0;iel<el->MeshElementNum();iel++)
    for(unsigned inode=el->ElementNodeNum(iel,1);inode<el->ElementNodeNum(iel,2);inode++)
      el->SetElementNodeGlobalIndex(iel,inode,0);  
  
  for(unsigned iel=0;iel<el->MeshElementNum();iel++){
    for(unsigned ifacet=0; ifacet<el->ElementFacetNum(iel,0); ifacet++){ // number of hex facets
      unsigned inode=el->ElementNodeNum(iel,1)+ifacet;
      if( 0==el->GlobalNodeIndex(iel,inode) ) {
	el->SetElementNodeGlobalIndex(iel,inode,++nvt);
	unsigned i1=el->FacetNodeGlobalIndex(iel,ifacet,0);
	unsigned i2=el->FacetNodeGlobalIndex(iel,ifacet,1);
	unsigned i3=el->FacetNodeGlobalIndex(iel,ifacet,2);
	
	for(unsigned j=0;j< el->NodeNumOfAdjacentElements(i1-1u);j++){
	  unsigned jel= el->NodeAdjacentElementIndex(i1-1u,j)-1u;
	  if(jel>iel){
	    for(unsigned jface=0;jface<el->ElementFacetNum(jel,0);jface++){
	      unsigned jnode=el->ElementNodeNum(jel,1)+jface;
	      if( 0==el->GlobalNodeIndex(jel,jnode) ) {
	      	unsigned j1=el->FacetNodeGlobalIndex(jel,jface,0); 
		unsigned j2=el->FacetNodeGlobalIndex(jel,jface,1);
		unsigned j3=el->FacetNodeGlobalIndex(jel,jface,2);
		unsigned j4=el->FacetNodeGlobalIndex(jel,jface,3);
		if((i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
		   (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
		   (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 )){
		  el->SetElementNodeGlobalIndex(jel,jnode,nvt);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  for(unsigned iel=0;iel<el->MeshElementNum();iel++){
    if(0==el->ElementType(iel)){
      el->SetElementNodeGlobalIndex(iel,26,++nvt);
    }
    if(3==el->ElementType(iel)){
      el->SetElementNodeGlobalIndex(iel,8,++nvt);
    }
  }
  el->SetNodeNumber(nvt);
  
  unsigned nv0= el->MeshVertexNum();
  unsigned nv1= el->MeshEdgeNum();
  el->SetMeshQuadFaceNum(nvt-nv0-nv1);
 
}

/**
 * This function stores the element adjacent to the element facet (iel,ifacet) 
 * and stores it in kel[iel][ifacet]   
 **/
void mesh::Buildkel(){ 
  for(unsigned iel=0;iel<el->MeshElementNum();iel++){
    for(unsigned iface=0;iface<el->ElementFacetNum(iel);iface++){
      if( el->FacetAdjacentElemIndex(iel,iface) <= 0){
	unsigned i1=el->FacetNodeGlobalIndex(iel,iface,0);
	unsigned i2=el->FacetNodeGlobalIndex(iel,iface,1);
	unsigned i3=el->FacetNodeGlobalIndex(iel,iface,2);
	for(unsigned j=0;j< el->NodeNumOfAdjacentElements(i1-1u);j++){
	  unsigned jel= el->NodeAdjacentElementIndex(i1-1u,j)-1u;
	  if(jel>iel){
	    for(unsigned jface=0;jface<el->ElementFacetNum(jel);jface++){
	      if( el->FacetAdjacentElemIndex(jel,jface) <= 0){
		unsigned j1=el->FacetNodeGlobalIndex(jel,jface,0); 
		unsigned j2=el->FacetNodeGlobalIndex(jel,jface,1);
		unsigned j3=el->FacetNodeGlobalIndex(jel,jface,2);
		unsigned j4=el->FacetNodeGlobalIndex(jel,jface,3);
		if((DIM[2]==1 &&
		    (i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
		    (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
		    (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 ))||
		   (DIM[1]==1 &&
		    (i1==j1 || i1==j2 )&&
		    (i2==j1 || i2==j2 ))||
		   (DIM[0]==1 &&
		    (i1==j1))
		  ){
		  el->SetFacetAdjacentElemIndex(iel,iface,jel+1u);
		  el->SetFacetAdjacentElemIndex(jel,jface,iel+1u);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

/**
 * Returns the number of different type of mesh nodes.
 **/
unsigned mesh::MeshNodeNum(const unsigned type) const{
  switch(type){
  case 0:
    return el->MeshVertexNum(); break;
  case 1:
    return el->MeshVertexNum()+el->MeshEdgeNum(); break;
  case 2:
    return el->MeshVertexNum()+el->MeshEdgeNum()+el->MeshQuadFaceNum(); break;
  case 3:
    return el->MeshVertexNum()+el->MeshEdgeNum()+el->MeshQuadFaceNum()+el->MeshTriFaceNum(); break;
  case 4:
    return el->MeshVertexNum()+el->MeshEdgeNum()+el->MeshQuadFaceNum()+el->MeshTriFaceNum()+el->MeshElementNum(); break;
  }
  return 0;
}

unsigned mesh::MeshElementNum() const{
  return nel;
}

unsigned mesh::GetGridNumber() const{
  return grid;
}

void mesh::set_elr(const  vector <vector <double> > &vt){
  el->InitRefinedToZero();   
  //refine all next grid even elements
  for(unsigned iel=0;iel<nel;iel+=1){
    unsigned nve=el->ElementNodeNum(iel,0);
    double vtx=0.,vty=0.,vtz=0.;
    
    for( unsigned i=0;i<nve;i++){
      unsigned inode=el->GlobalNodeIndex(iel,i)-1u;	
      vtx+=vt[0][inode];
      vty+=vt[1][inode];
      vtz+=vt[2][inode];
    }
    
    vtx/=nve;  vty/=nve;  vtz/=nve;
  
    if( vty<1. && vtx<=0.25){
      el->SetRefinedElementIndex(iel,1);
      el->AddToRefinedElementNumber(1);
      short unsigned elt=el->ElementType(iel);
      el->AddToRefinedElementNumber(1,elt);
    }
  }
}

/**
 * This function flags the elements that will be refined   
 **/
void mesh::set_elr(const unsigned &test){
  el->InitRefinedToZero();

  if(test==0){
    //refine all next grid elements
    for(unsigned iel=0;iel<nel;iel++){
      el->SetRefinedElementIndex(iel,1);
      el->AddToRefinedElementNumber(1);
      short unsigned elt=el->ElementType(iel);
      el->AddToRefinedElementNumber(1,elt);
    }
  }
  else if(test==1){
    //refine all next grid even elements
    for(unsigned iel=0;iel<nel;iel+=2){
      el->SetRefinedElementIndex(iel,1);
      el->AddToRefinedElementNumber(1);
      short unsigned elt=el->ElementType(iel);
      el->AddToRefinedElementNumber(1,elt);
    }
  }
  else if(test==2){ 
    //refine all next grid odd elements
    for(unsigned iel=1;iel<nel;iel+=2){
      el->SetRefinedElementIndex(iel,1);
      el->AddToRefinedElementNumber(1);
      short unsigned elt=el->ElementType(iel);
      el->AddToRefinedElementNumber(1,elt);
    } 
  }
}


/**
 * This function copies the refined element index vector in other_vector   
 **/
void mesh::copy_elr(vector <unsigned> &other_vec) const{
  for(unsigned i=0;i<nel;i++)
    other_vec[i]=el->GetRefinedElementIndex(i);
}
//TODO Eugnenio's code for reordering 
// //*************** start reorder mesh dofs **************
//   //(1)linear (2)quadratic (3)biquaratic
//  
//   vector <unsigned> dof_index;
//   dof_index.resize(nvt);
//   for(unsigned i=0;i<nvt;i++){
//     dof_index[i]=i+1;
//   }
//   //reorder vertices and mid-points vs central points
//   for (unsigned iel=0; iel<nel; iel++) {
//     for (unsigned inode=0; inode<el->GetElementDofNumber(iel,1); inode++) {
//       for (unsigned jel=0; jel<nel; jel++) {
//     for (unsigned jnode=el->GetElementDofNumber(jel,1); jnode<el->GetElementDofNumber(jel,3); jnode++) {
//       unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
//       unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
//       unsigned i0=dof_index[ii];
//       unsigned i1=dof_index[jj];
//       if(i0>i1){
//         dof_index[ii]=i1;
//         dof_index[jj]=i0;
//       }
//     }
//       }
//     }
//   }
//   //reorder vertices vs mid-points
//   for (unsigned iel=0; iel<nel; iel++) {
//     for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
//       for (unsigned jel=0; jel<nel; jel++) {
//         for (unsigned jnode=el->GetElementDofNumber(jel,0); jnode<el->GetElementDofNumber(jel,1); jnode++) {
//           unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
//       unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
//       unsigned i0=dof_index[ii];
//           unsigned i1=dof_index[jj];
//       if(i0>i1){
//         dof_index[ii]=i1;
//         dof_index[jj]=i0;
//       }
//     }
//       }
//     }
//   }
//  
//   // update all
//   for (unsigned iel=0; iel<nel; iel++) {
//     for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++) {
//       unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
//       el->SetElementVertexIndex(iel,inode,dof_index[ii]);
//     }
//   }
//   vector <double> vt_temp;
//   for(int i=0;i<3;i++){
//     vt_temp=vt[i];
//     for(unsigned j=0;j<nvt;j++){
//       vt[i][dof_index[j]-1]=vt_temp[j];
//     }
//   }
//   // **************  end reoreder mesh dofs **************