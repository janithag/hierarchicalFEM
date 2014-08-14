#include "main.hpp"
#include "multigrid.hpp"
#include "elem_type.hpp"
#include <ctime>

MultiGrid::~MultiGrid(){ 
  for(unsigned i=0; i<p; i++) 
    delete D[i];

  // if(time_test)
  //   for(unsigned i=gridr; i<gridn; i++) delete D_old[i];
  
  for(unsigned i=0; i<6; i++)
    for(unsigned j=0; j<maxOrder; j++)      
      delete type_elem[i][j];
  
  for(unsigned i=0; i<SolName.size(); i++)  
    delete[] SolName[i];
};

MultiGrid::MultiGrid(const char mesh_file[], const unsigned short& gridP): // 0-based index used throughout the class.
  gridr(gridP){
  time_test = 0;
  int gaussOrder = 18;
  p = gridP;
  gridn=p;
  
  try{
    DOF.resize(p);
    D.resize(p);
    vt.resize(3); 
    
  clock_t start_time, end_time, presmoothing, assembleRes, prolongate, resT0, resT1;
  cout<<"\nLevel: "<<p<<endl;
  
    for(unsigned i=0;i<p;i++){
      //start_time = clock();
      D[i] = new lsysPDE(mesh_file, vt, i);  
      //end_time = clock();
      //cout<<"Time elapsed for mesh at level "<<i<<": "<<(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;
      if(i<p-1) 
	D[i]->set_elr(0);
      else 
	D[i]->set_elr();
    }
    
    // generate solution vector
    AddSolutionVector("U", 0);
    //AddSolutionVector("V", 0);
    //AddSolutionVector("W", 0);
  
    if(SolType.size()>0)
      maxOrder = *max_element( SolType.begin(), SolType.end() );
    
    for(unsigned i=0;i<p;i++){
      D[i]->el->SetGlobalDOF(maxOrder);
    }  
    
    setDOF();
    //cout<<"# of dofs: "<<DOF[D[p-1]->SolType[0]]<<endl;
    //construct 1D, 2D and 3D elements.
    type_elem.resize(6);  
    for(unsigned i=0;i<6;i++){
      type_elem[i].resize(maxOrder);
      for(unsigned j=0;j<maxOrder;j++)
	type_elem[i][j] = new const elem_type(i, j+1, gaussOrder, 1);
    }
    
   // Initialize 
    Initialize("U");  
//     Initialize("V");
//     Initialize("W");
   
   // Set Boundary conditions
    GenerateBC("U");
//     GenerateBC("V");
//     GenerateBC("W");

    //start Multigrid for UVW
    ClearMGIndex();
    AddToMGIndex("U"); 
    //AddToMGIndex("V");
      //AddToMGIndex("W");

    InitMultigrid();
    for(unsigned g=1;g<p;g++){
      BuildProlongatorMatrix(g);
    }
  //start_time = clock();
   D[gridn-1]->AssembleMatrix(type_elem,vt,gridn-1u,"All");
    for(int ig=gridn-1;ig>0;ig--){
      MatPtAP(D[ig]->K, D[ig]->P, MAT_INITIAL_MATRIX,1.0,&D[ig-1]->CC); 
      D[ig-1u]->AssembleMatrix(type_elem,vt,gridn-1u,"All");
      MatAXPY(D[ig-1u]->K,1,D[ig-1u]->CC,DIFFERENT_NONZERO_PATTERN);
      MatDestroy(&(D[ig-1u]->CC)); 
    }
//     D[gridn-1u]->SetResZero(MGIndex);  
//     D[gridn-1u]->SetEpsZero(MGIndex);
//     D[gridn-1u]->AssembleResidual(type_elem,vt,Sol,gridn-1u,"All",p-1u);
//     D[gridn-1]->VankaPETSCSmoother(BC,MGIndex);
    //end_time = clock();
   //cout<<"Time elapsed for assembly: "<<(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;
    
     /* Full Multigrid */   
    start_time = clock();
    presmoothing=0.;assembleRes=0.;//prolongate=0.;
    for(unsigned igridn=1; igridn<=p; igridn++){
       for(unsigned icycle=0; icycle<2; icycle++){//<5
         
    	D[igridn-1u]->SetResZero(MGIndex);  
    	D[igridn-1u]->SetEpsZero(MGIndex);
	//resT0 = clock();
    	D[igridn-1u]->AssembleResidual(type_elem,vt,Sol,igridn-1u,"All",p-1u);
    	//resT1 = clock();
	//assembleRes += resT1 -resT0;
        //presmoothing
	for(unsigned ig=igridn-1u;ig>0;ig--){
	  //start_time = clock();
    	  for(unsigned k=0;k<2;k++)	
    	    //D[ig]->SpaceDecompositionSmoother(BC,MGIndex, DofstoSolve[ig]);	
	    D[ig]->VankaPETSCSmoother(BC,MGIndex,DofstoSolve[ig]);
	  //end_time = clock(); 
	  //presmoothing += end_time-start_time;
	  //cout<<"Time elapsed presmoothing step: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;
    	  D[ig-1u]->SetResZero(MGIndex);
    	  D[ig-1u]->SetEpsZero(MGIndex);
    	  Restrictor(ig);
	  
	 //assemble residual only on the part of the coarse grid that is not refined
    	 //D[ig-1u]->AssembleResidual(type_elem,vt,Sol,igridn-1u,"All",p-1u);
    	  
    	  //cout<<"RES TIME "<<ig<<"-->"<<ig-1<<" ="<<(end_time-start_time)/CLOCKS_PER_SEC<<endl;
    	}
    	
        //direct solver on the coarse level
    	for(unsigned k=0;k<1;k++)
    	  //D[0]->SpaceDecompositionSmoother(BC,MGIndex,DofstoSolve[0]);
	  D[0]->VankaPETSCSmoother(BC,MGIndex,DofstoSolve[0]);
    	
         //post smoothing
    	for(unsigned ig=1;ig<igridn;ig++){
    	  Prolongator(ig);
    	  D[ig]->UpdateResidual();
    	  D[ig]->SumEpsCToEps(MGIndex);
	 
    	  for(unsigned k=0;k<2;k++)
    	    //D[ig]->SpaceDecompositionSmoother(BC,MGIndex,DofstoSolve[ig]);
	    D[ig]->VankaPETSCSmoother(BC,MGIndex,DofstoSolve[ig]);
    	}
          
    	for(unsigned ig=0;ig<igridn;ig++)
    	  D[ig]->FindResMax(BC,MGIndex);
    	
    	D[igridn-1u]->SumEpsToSol(Sol,MGIndex);
    	if(igridn==1) icycle += 100;
	
	//stopping criteria
// 	if(D[igridn-1u]->ResInf[0]<pow(10,-6)){
// 	  //cout<<"cycles at level "<<igridn<<" : "<<icycle+1<<endl; 
// 	  break;
// 	}
        }
          
    //prolongate the solution into the next finer level: automatically done by the way Sol is defined and interpolation operator works.
       
        }
	end_time = clock(); 
	cout<<"FMG: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;
	//cout<<"assembleRes: "<<static_cast<double>(assembleRes)/CLOCKS_PER_SEC<<"s"<<endl;
	//cout<<"presmoothing: "<<static_cast<double>(presmoothing)/CLOCKS_PER_SEC<<"s"<<endl;
	//cout<<"prolongate, restrict: "<<static_cast<double>(prolongate)/CLOCKS_PER_SEC<<"s"<<endl;
	
        L2Error();
	
	for(int i=0;i<p;i++)
	  Print_vtk_ASCII(0,i);
	for(unsigned ig=0;ig<p;ig++)
	  D[ig]->FreeSolutionVectors();
        for(unsigned ig=0;ig<gridn;ig++)
          D[ig]->DeallocateMatrix(MGIndex);
   
//      PetscViewer viewer;
//        for(int ig= p-1;ig>0;ig--){
//        PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);
//        MatView(D[ig]->K,viewer);
//        double f;
//        cin>>f;
//      }
 }
 
  catch(exception& ex){
    cerr<<"Detected exception: "<<ex.what()<<endl;
  }
}

/* For a given p, this function sets the number of DOF at each grid level. */
void MultiGrid::setDOF(){
  DOF.resize(p+maxOrder);
  for (unsigned j=0;j<maxOrder;j++){
    DOF[j] = D[p-1]->el->NumOfDof(j); 
    //cout<<DOF[j]<<"\t";
  }
} 

void MultiGrid::AddSolutionVector(const char name[], const unsigned order){
  unsigned n = Sol.size(); // number of solution vectors
  Sol.resize(n+1u);
  SolType.resize(n+1u);
  SolName.resize(n+1u);
  BC.resize(n+1u);

  SolType[n]=order+p;
 
  SolName[n]= new char[8];
  strcpy(SolName[n],name);
  
  for(unsigned i=0; i<p; i++)
    D[i]->AddSolutionVector(name,order+i);
}

void MultiGrid::Initialize(const char name[]){
  unsigned i_start;
  unsigned i_end;
  if(!strcmp(name,"All")){  
    i_start=0;
    i_end=Sol.size();
  }
  else{
    i_start=GetIndex(name);
    i_end=i_start+1u;
  }
  
  for(unsigned i=i_start; i<i_end; i++){
    CheckVectorSize(i); 
    for(unsigned j=0; j< D[p-1]->el->NumOfDof(SolType[i]); j++) 
      Sol[i][j]=0.;
  }
}

void MultiGrid::CheckVectorSize(const unsigned& i){
  if(Sol[i].size() < D[p-1]->el->NumOfDof(SolType[i])) 
    ResizeSolutionVector(SolName[i]);
}

void MultiGrid::ResizeSolutionVector(const char name[]){
  unsigned i_start;
  unsigned i_end;
  if(!strcmp(name,"All")){  
    i_start=0;
    i_end=Sol.size();
  }
  else{
    i_start=GetIndex(name);
    i_end=i_start+1u;
  }
  
  for(unsigned i=i_start;i<i_end;i++){
    Sol[i].resize(D[p-1]->el->NumOfDof(SolType[i])); 
    BC[i].resize(D[p-1]->el->NumOfDof(SolType[i]));	
    for(unsigned j=0; j<p; j++)
      D[j]->ResizeSolutionVector(SolName[i]);
  }
}

unsigned MultiGrid::GetIndex(const char name[]){
  unsigned index = 0;
  while(strcmp(SolName[index],name)){
    index++;
    if(index==Sol.size()) 
      throw runtime_error("Invalid name entry GetIndex(...)" );
  }
  return index;
}

void MultiGrid::ClearMGIndex(){
  MGIndex.clear();
};

void MultiGrid::InitMultigrid(){

  for(unsigned i=0;i<gridn; i++)
    D[i]->InitMultigrid(MGIndex);  
}
    
void MultiGrid::AddToMGIndex(const char name[]){
  unsigned n=MGIndex.size();
  MGIndex.resize(n+1u);
  MGIndex[n]=GetIndex(name);
};

int MultiGrid::Restrictor(unsigned gridf){
  PetscErrorCode ierr;
  ierr = MatMult(D[gridf]->R,D[gridf]->RES,D[gridf-1]->RES); CHKERRQ(ierr);
  return 1;
}

int MultiGrid::Prolongator(unsigned gridf){
  PetscErrorCode ierr;
  ierr = MatMult(D[gridf]->P,D[gridf-1]->EPS,D[gridf]->EPSC); CHKERRQ(ierr);
  return 1;
}

int MultiGrid::BuildProlongatorMatrix(unsigned gridf){ 

  if(gridf<1)
    throw runtime_error("Invalid input argument in function \"BuildProlongatorMatrix\" ");
 
  PetscErrorCode ierr;
  unsigned gridc = gridf-1;
  PetscInt nf = D[gridf]->KIndex[D[gridf]->KIndex.size()-1u];
  PetscInt nc = D[gridc]->KIndex[D[gridc]->KIndex.size()-1u];
  const PetscScalar value = 1;
  PetscInt i, j;
  
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nf,nc,1,PETSC_NULL,&D[gridf]->P); CHKERRQ(ierr);

  for(unsigned k=0; k<MGIndex.size(); k++){
    unsigned SolIndex = MGIndex[k];
    for(unsigned m=0; m<DOF[D[gridc]->SolType[SolIndex]]; m++){
      i = D[gridf]->KIndex[k]+m;
      j = D[gridc]->KIndex[k]+m;
      ierr = MatSetValues(D[gridf]->P,1,&i,1,&j,&value,INSERT_VALUES); CHKERRQ(ierr);
    } 
  }

  ierr = MatAssemblyBegin(D[gridf]->P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(D[gridf]->P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatTranspose(D[gridf]->P,MAT_INITIAL_MATRIX,&D[gridf]->R);

  /* Visualize the output */
  // PetscViewer viewer;
  // PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);
  // MatView(D[gridf]->P,viewer);
  // double f;
  // cin>>f;
  // MatView(D[gridf]->P,PETSC_VIEWER_STDOUT_WORLD);

  return 1;
 
}

void MultiGrid::solve(const unsigned short &time){
  //copy the old meshes from 0 to gridn-1
  time_test=1;
  for(unsigned ig=gridr;ig<gridn;ig++){
    if(time-1u) delete D_old[ig];
    D_old[ig]=D[ig];
  }
  D[gridr-1u]->copy_elr(elr_old);
  
  //crete the new meshes from gridr to gridn-1
  for(unsigned i=gridr;i<gridn;i++){
    D[i-1u]->set_elr(1+time%2);
    D[i] = new lsysPDE(i,D[i-1u]->el);
    generate_vt(i+1u,i);
    GenerateBC("All"); //GenerateBC("All",i+1u,i);
  }
}

void MultiGrid::generate_vt(const unsigned short &grid_end, const unsigned short &grid_start){
  if( (unsigned) vt[0].size()< D[grid_end-1u]->MeshNodeNum(4)){
    vt[0].resize( D[grid_end-1u]->MeshNodeNum(4));
    vt[1].resize( D[grid_end-1u]->MeshNodeNum(4));
    vt[2].resize( D[grid_end-1u]->MeshNodeNum(4));
  }
}

void MultiGrid::GenerateBC(const char name[]){ // grid numbers in 0-based index

  unsigned fGrid=p-1;
  unsigned cGrid=0;
  
  /* hex-0, tet-1, wedge-2, quad-3, tri-4, line-5 */
  const short unsigned NV1[6][2]={{9,9},{7,7},{9,7},{3,3},{3,3},{1,1}}; // number of nodes on a facet.
  
  unsigned startDOF = (0==cGrid) ? 0 : DOF[cGrid-1];
  unsigned i_start;
  unsigned i_end;

  if(!strcmp(name,"All")){  
    i_start=0;
    i_end=Sol.size();
  }

  else{
    i_start=GetIndex(name);
    i_end=i_start+1u;
  }
  
  for(unsigned i=i_start;i<i_end;i++){ //loop on solution vectors

    BC[i].resize(DOF[D[fGrid]->SolType[i]]);
    
    /* Set Neumann boudary + interior dofs */
    for(unsigned j=startDOF; j<DOF[D[fGrid]->SolType[i]]; j++)
      BC[i][j]=20;

     unsigned g=fGrid; // we build BC for fine grid
     double x,y,z;
    
    // flag the dofs to update only
    for(unsigned iel=0; iel < D[g]->MeshElementNum(); iel++){
      
      // find centroid
      unsigned numOfVertices=D[g]->el->ElementNodeNum(iel,0);
      int vIndex[numOfVertices];
      x=0.;y=0.;z=0.;
      
      for(unsigned k=0;k<numOfVertices;k++){
	vIndex[k] = D[g]->el->GlobalNodeIndex(iel,k)-1u;	
	x += vt[0][vIndex[k]];
	y += vt[1][vIndex[k]];
	z += vt[2][vIndex[k]];
      }
      
      int localP = localP_(x/numOfVertices,y/numOfVertices,z/numOfVertices,p); 
      
      for(unsigned k=type_elem[D[g]->el->ElementType(iel)][localP-1]->ElementDOF(); k<type_elem[D[g]->el->ElementType(iel)][p-1]->ElementDOF();k++){
	int j=D[g]->el->globalDOF_(iel)[k];
	if (localP<BC[i][j])
	  BC[i][j]=localP;
      }
    }
    
    /* Set Dirichlet BC nodes*/
   
    //for(unsigned g=cGrid; g <= fGrid; g++)// loop on each grid for Domain decomposition
    for(unsigned iel=0; iel < D[g]->MeshElementNum(); iel++)// loop on each element of the grid
      for(unsigned facet=0; facet < D[g]->el->ElementFacetNum(iel); facet++){ // loop on each facet of each element
	 
	if( D[g]->el->FacetAdjacentElemIndex(iel,facet) < 0 ){ // if the facet is on the boundary

	  short unsigned ielt = D[g]->el->ElementType(iel);
	  unsigned nNodes = NV1[ielt][facet>=D[g]->el->ElementFacetNum(iel,0)]; 
	   
	  for(unsigned k=0; k<nNodes; k++){
 
	    unsigned iNode = D[g]->el->FacetNodeGlobalIndex(iel,facet,k)-1u; 
	    double value;
	    bool test;
	    if( iNode<D[g]->MeshNodeNum(0) ){ //if the the node is a vertex 
	      test = Boundary[DIMENSION](vt[0][iNode],vt[1][iNode],vt[2][iNode],SolName[i],value,-(D[g]->el->FacetAdjacentElemIndex(iel,facet)+1));
	
	      if(test){
		BC[i][iNode] = 0; // Dirichlet boundary
		Sol[i][iNode] = value; // Dirichlet BC (We consider only linear case here.)
	      }
	    }
	     
	    else { // if the node is an edge or a face
	      
	      test = Boundary[DIMENSION](vt[0][iNode],vt[1][iNode],vt[2][iNode],SolName[i],value,-(D[g]->el->FacetAdjacentElemIndex(iel,facet)+1));

	      if(!test) break;
	      vector<unsigned> iMode;
	      iMode.resize(0);
	      D[g]->NodeGlobalDofOfModes(iNode, iMode);
	      for(unsigned c=0;c<iMode.size();c++){
		BC[i][iMode[c]] = 0; //Dirichlet boundary
		Sol[i][iMode[c]] = 0; //Dirichlet BC
	      }
	    }
	  }
	}
	 
      }
    // count dofs to solve
    DofstoSolve.resize(p);
    for(int level=0;level<p;level++){
      int count=0;
      for(unsigned j=startDOF; j<DOF[D[level]->SolType[i]]; j++)
      {
	//cout<<BC[i][j]<<" ";
	if(BC[i][j]>level) count++;
      }
      DofstoSolve[level]=count;
    }cout<<"#dofs: "<<DofstoSolve[p-1]<<endl;
  }
}

/* L2 error */
void MultiGrid::L2Error() const{ 

  double L2Err;
  double ElementL2Err=0.;
  double uExact=0.;
  double  u=0.;
  int order = SolType[0];
  short unsigned kelt;
  double* phi;
  double** gradphi;
  double* X;
  double w=0.;
  vector < vector <double> > vertex;
  int elemDOF;
 
  unsigned maxDOF = type_elem[0][order-1]->ElementDOF();
  for(unsigned i=1;i<6;i++)
    if (type_elem[i][order-1]->ElementDOF() > maxDOF)
      maxDOF = type_elem[i][order-1]->ElementDOF();

  phi = new double[maxDOF];
  X = new double[3];
  
  gradphi = new double*[maxDOF];
  for(unsigned c=0; c<maxDOF; c++)
    gradphi[c]= new double[3];
  
  int* s;
  int* t;
  //for(int side=0;side<2;side++){
    L2Err=0.;
  for (unsigned kel=0; kel < D[p-1]->MeshElementNum(); kel++){ 
    
    /* vertices of the element */
    unsigned numOfVertices = D[p-1]->el->ElementNodeNum(kel,0);
    vertex.resize(numOfVertices);
    int vIndex[numOfVertices];
    double x=0.;
    for(unsigned i=0; i< numOfVertices; i++){
      vertex[i].resize(3);
      vIndex[i] = D[p-1]->el->GlobalNodeIndex(kel,i)-1u;	
      vertex[i][0]=vt[0][vIndex[i]];
      vertex[i][1]=vt[1][vIndex[i]];
      vertex[i][2]=vt[2][vIndex[i]];
      x += vt[0][vIndex[i]];
      //cout<<"vertices: "<<vertex[i][0]<<" "<<vertex[i][1]<<" "<<vertex[i][2]<<endl;
    }
    //if(((x/numOfVertices>0.5)&&(side==0))||((x/numOfVertices<0.5) && (side==1))) continue;
    
    kelt = D[p-1]->el->ElementType(kel);
    elemDOF = type_elem[kelt][order-1]->ElementDOF();
    ElementL2Err=0.;

    s=D[p-1]->sign_(kel);t=D[p-1]->type_(kel);
    
    for(unsigned ig=0; ig < type_elem[kelt][order-1]->NumOfGaussPts(); ig++){
      (type_elem[kelt][order-1]->*type_elem[kelt][order-1]->spatialToRefPtr)(vertex,ig,s,w,phi,gradphi,X,t);
  
     //uExact = sin(2*PI*X[0])*sin(2*PI*X[1]); // unit square 2D - trig
     //uExact = X[0]*X[1]*(1-X[0])*(1-X[1]);// unit square 2D - polynomial
      //uExact=X[0]*sin(4*PI*pow(X[0],2))*sin(4*PI*X[1]);//2D
      
     uExact = sin(2*PI*(X[0]))*sin(2*PI*(X[1]))*sin(2*PI*X[2]); // 3D unit cube trig
     // uExact = pow(X[0]*X[1]*X[2]*(1-X[0]-X[1]-X[2]),2);//3D unit tet order 8 polynomial
      //uExact = X[0]*X[1]*X[2]*(1-X[0])*(1-X[1])*(1-X[2]);// 3D unit cube polynomial
       //uExact=X[0]*sin(4*PI*pow(X[0],2))*sin(2*PI*X[1])*sin(2*PI*X[2]);
     
       // uExact = X[1]*X[2]*(X[0]-1)*(X[1]-1)*(X[2]-X[0]);//unit wedge all dirichlet
      //uExact = X[2]*(X[0]-1)*(X[2]-X[0]);//unit wedge Neumann in tri faces 4,5
      //uExact =  X[1]*(1-X[1]); //unit wedge Neumann on all quad faces

      u=0.;
      
      for(int i=0;i<elemDOF;i++)
	u += Sol[0][D[p-1]->el->globalDOF_(kel)[i]]*phi[i];
	
      ElementL2Err += pow(uExact-u,2)*w;
    }
    //cout<<"element l2 error: "<<ElementL2Err<<endl;
    L2Err += ElementL2Err;
  }
  cout << "L2 Error: "<<sqrt(L2Err)<< endl;
//}
//     ofstream fout;
//     fout.open("L2 error.txt");
//     fout.precision(14);
//     fout<<p<<"\t"<<DOF[p-1]<<"\t"<<sqrt(L2Err)<<endl;
//     fout.close();
  delete[] phi;
  delete[] X;
  for(unsigned c=0;c<maxDOF;c++)
    delete[] gradphi[c];
  delete[] gradphi;
}

void MultiGrid::GenerateLagrangeSolution2D(int LagrangeOrder){ // extend this for multiple variables. LagrangeOrder is 1-based
  
  lagrangeP=LagrangeOrder;
  
  LSol.resize(Sol.size());
  unsigned nel = D[p-1]->MeshElementNum();
  unsigned LElemDof;
  unsigned HElemDof; 
  unsigned LSolSize=0; 
  unsigned kelt;
  vector<unsigned> order;
  order.resize(Sol.size());
  vector < vector <double> > vertex; // coordinates of the vertices of each element
  vector < vector <double> > x; // global lagrange nodes
  
  basis* pt[6];
  line line_ = line();
  quad quad_ = quad();
  tri tri_ = tri();
  hex hex_ = hex();
  tet tet_ = tet();
  wedge wedge_ = wedge();
  
  pt[0]=&hex_;
  pt[1]=&tet_;
  pt[2]=&wedge_;
  pt[3]=&quad_;
  pt[4]=&tri_;
  pt[5]=&line_;
    
  for(unsigned kel=0; kel < nel; kel++){
    kelt = D[p-1]->el->ElementType(kel);
    LElemDof = pt[kelt]->ElementDOF(LagrangeOrder); // LagrangeOrder should be 1-based here.
    LSolSize += LElemDof;
  }
 
  for(unsigned i=0;i<LSol.size();i++){
    LSol[i].resize(LSolSize);
    order[i] = SolType[i];
  }
  
  LNode.resize(LSolSize);
  for(unsigned i=0;i<LSolSize; i++)
    LNode[i].resize(3);
  
  int offset=0;
  
  unsigned HexElemDof=pt[0]->ElementDOF(LagrangeOrder);
  double** phi = new double*[HexElemDof]; // create projection matrix to the size of hex
  
  for(unsigned i=0; i < HexElemDof; i++)
    phi[i] = new double[type_elem[3][order[0]-1]->ElementDOF()]; //size of quad. change this to hex when you do 3-D
  
  for(unsigned kel=0; kel < nel; kel++){
    
    kelt = D[p-1]->el->ElementType(kel);
    LElemDof = pt[kelt]->ElementDOF(LagrangeOrder);
    HElemDof = type_elem[kelt][order[0]-1]->ElementDOF(); 
    x.resize(LElemDof);
     
    for(unsigned i=0;i<LElemDof;i++)
      x[i].resize(3);
       
    vertex.resize(D[p-1]->el->ElementNodeNum(kel,0));
    
    int vIndex[D[p-1]->el->ElementNodeNum(kel,0)]; // global indices of the vertices of the element
    for(unsigned i=0;i<D[p-1]->el->ElementNodeNum(kel,0);i++){
      vertex[i].resize(3);
      vIndex[i] = D[p-1]->el->GlobalNodeIndex(kel,i)-1u;	
      vertex[i][0]=vt[0][vIndex[i]];
      vertex[i][1]=vt[1][vIndex[i]];
      vertex[i][2]=vt[2][vIndex[i]];
    }
  
    (type_elem[kelt][order[0]-1]->*type_elem[kelt][order[0]-1]->projectionMatrix)(LagrangeOrder, D[p-1]->sign_(kel), vertex, phi, x, D[p-1]->type_(kel)); 
        
    for(unsigned i=0;i<LElemDof;i++){
      LSol[0][i+offset]=0.; 
      for(unsigned j=0;j<HElemDof;j++)
	LSol[0][i+offset] += Sol[0][D[p-1]->el->globalDOF_(kel)[j]]*phi[i][j];
    }
     
    for(unsigned i=0;i<LElemDof;i++){ 
      LNode[i+offset][0]=x[i][0];
      LNode[i+offset][1]=x[i][1];
      LNode[i+offset][2]=x[i][2];
    }
     
    offset += LElemDof;
     
  }
  
  for(unsigned i=0;i<HexElemDof;i++)
    delete[] phi[i];
  delete[] phi;
}

void MultiGrid::GenerateLagrangeSolution3D(int LagrangeOrder){ // extend this for multiple variables. LagrangeOrder is 1-based
  
  lagrangeP=LagrangeOrder;
  
  LSol.resize(Sol.size());
  unsigned nel = D[p-1]->MeshElementNum();
  unsigned LElemDof;
  unsigned HElemDof; 
  unsigned LSolSize=0; 
  unsigned kelt;
  vector<unsigned> order;
  order.resize(Sol.size());
  vector < vector <double> > vertex; // coordinates of the vertices of each element
  vector < vector <double> > x; // global lagrange nodes
  
  basis* pt[3];

  hex hex_ = hex();
  tet tet_ = tet();
  wedge wedge_ = wedge();
  
  pt[0]=&hex_;
  pt[1]=&tet_;
  pt[2]=&wedge_;

    LElemDof = pt[0]->ElementDOF(LagrangeOrder); // LagrangeOrder should be 1-based here.
    LSolSize = LElemDof*nel;
 
  for(unsigned i=0;i<LSol.size();i++){
    LSol[i].resize(LSolSize);
    order[i] = SolType[i];
  }
  
  LNode.resize(LSolSize);
  for(unsigned i=0;i<LSolSize; i++)
    LNode[i].resize(3);
  
  x.resize(LElemDof);
  for(unsigned i=0;i<LElemDof;i++)
    x[i].resize(3);
    
  int offset=0;
  
  double** phi = new double*[LElemDof]; // create projection matrix to the size of hex
  
  for(unsigned i=0; i < LElemDof; i++)
    phi[i] = new double[type_elem[0][order[0]-1]->ElementDOF()]; //size of quad. change this to hex when you do 3-D
  
  for(unsigned kel=0; kel < nel; kel++){
    
    kelt = D[p-1]->el->ElementType(kel);
    HElemDof = type_elem[kelt][order[0]-1]->ElementDOF(); 
       
    vertex.resize(D[p-1]->el->ElementNodeNum(kel,0));
    
    int vIndex[D[p-1]->el->ElementNodeNum(kel,0)]; // global indices of the vertices of the element
    for(unsigned i=0;i<D[p-1]->el->ElementNodeNum(kel,0);i++){
      vertex[i].resize(3);
      vIndex[i] = D[p-1]->el->GlobalNodeIndex(kel,i)-1u;	
      vertex[i][0]=vt[0][vIndex[i]];
      vertex[i][1]=vt[1][vIndex[i]];
      vertex[i][2]=vt[2][vIndex[i]];
    }
  
    (type_elem[kelt][order[0]-1]->*type_elem[kelt][order[0]-1]->projectionMatrix)(LagrangeOrder, D[p-1]->sign_(kel), vertex, phi, x, D[p-1]->type_(kel)); 
    
    for(unsigned i=0;i<LElemDof;i++){
      LSol[0][i+offset]=0.; 
      for(unsigned j=0;j<HElemDof;j++)
	LSol[0][i+offset] += Sol[0][D[p-1]->el->globalDOF_(kel)[j]]*phi[i][j];
    }
     
    for(unsigned i=0;i<LElemDof;i++){ 
      LNode[i+offset][0]=x[i][0];
      LNode[i+offset][1]=x[i][1];
      LNode[i+offset][2]=x[i][2];
    }
     
    offset += LElemDof;
     
  }
  
  for(unsigned i=0;i<LElemDof;i++)
    delete[] phi[i];
  delete[] phi;
}
int MultiGrid::LagrangeHRefinement2D(){ 
  // for line: Number of splitted elements per single element = lagrangeP
  // for quad, tri: number of splitted elements per single element = lagrangeP^2
  
  unsigned nel = D[p-1]->MeshElementNum();
  unsigned kelt;
  unsigned short nVertices;
  int row, col;
  int drow;
  unsigned blockElementNum=pow(lagrangeP,2); //lagrangeP is 1-based
  const int cellType[6]={12,10,13,9,5,3};
  unsigned elemOffset = 0;
  unsigned nodeOffset = 0;
  int rowChanged;
   
  unsigned CellSize = nel*blockElementNum;
  Cell.resize(CellSize);
  CellType.resize(CellSize);
   
  basis* pt[2];

  quad quad_ = quad();
  tri tri_ = tri();
 
  pt[0]=&quad_;
  pt[1]=&tri_;
   
  for(unsigned kel=0; kel < nel; kel++){
     
    nVertices=D[p-1]->el->ElementNodeNum(kel,0);
    kelt=D[p-1]->el->ElementType(kel);
     
    row = 0; drow=0;
    for(unsigned i=0;i<blockElementNum;i++){
      Cell[i+elemOffset].resize(nVertices);
      
      if(kelt==3){ //quadrilateral
	row = i/lagrangeP; col=i%lagrangeP;
	Cell[i+elemOffset][0]=row*(lagrangeP+1)+col+nodeOffset;
	Cell[i+elemOffset][1]=row*(lagrangeP+1)+col+1+nodeOffset;
	Cell[i+elemOffset][2]=(row+1)*(lagrangeP+1)+col+1+nodeOffset;
	Cell[i+elemOffset][3]=(row+1)*(lagrangeP+1)+col+nodeOffset;
      }
      
      else{ //triangle
	rowChanged=false;
	
	if (i==0){
	  Cell[i+elemOffset][0]= 0+nodeOffset;
	  Cell[i+elemOffset][1]= 1+nodeOffset;
	  Cell[i+elemOffset][2]= 1+lagrangeP+nodeOffset;
	}
	
	else if(i < lagrangeP*(lagrangeP+1)/2){// upright triangles
	  unsigned sum = 0;
	  for(int k=0;k<row+1;k++)
	    sum += lagrangeP-k;
	
	  if(i==sum){rowChanged=1; row++;}
	 
	  Cell[i+elemOffset][0]= Cell[i-1+elemOffset][1]+rowChanged;
	  Cell[i+elemOffset][1]= Cell[i+elemOffset][0]+1;
	  Cell[i+elemOffset][2]= Cell[i+elemOffset][1]+lagrangeP-row;
	}
	else {// upside down triangles
	  
	  unsigned k = i-(lagrangeP)*(lagrangeP+1)/2;
	  if (k==0){
	    Cell[i+elemOffset][0]= 1+nodeOffset;
	    Cell[i+elemOffset][1]= 2+lagrangeP+nodeOffset;
	    Cell[i+elemOffset][2]= Cell[i+elemOffset][1]-1;
	  } 
	  else{
	    unsigned sum = 0;
	    for(int kk=0;kk<drow+1;kk++)
	      sum += lagrangeP-1-kk;
	
	    if(k==sum){rowChanged=2; drow++;}
	 
	    Cell[i+elemOffset][0]=Cell[i-1+elemOffset][0]+1+rowChanged;
	    Cell[i+elemOffset][1]=Cell[i+elemOffset][0]+lagrangeP-drow+1;
	    Cell[i+elemOffset][2]=Cell[i+elemOffset][1]-1;
	  }
	  
	}
      }
      CellType[i+elemOffset]=cellType[kelt];
      
    }
    elemOffset += blockElementNum;
    nodeOffset += pt[kelt-3]->ElementDOF(lagrangeP);
  }
   
  return D[p-1]->el->MeshElementNum("Quad")*blockElementNum*5 + D[p-1]->el->MeshElementNum("Triangle")*blockElementNum*4;
}

int MultiGrid::LagrangeHRefinement3D(){ 
  
  unsigned nel = D[p-1]->MeshElementNum();
  unsigned short nVertices=8;
  int row, col, height;
  unsigned blockElementNum=pow(lagrangeP,3); //lagrangeP is 1-based
  
  unsigned elemOffset = 0;
  unsigned nodeOffset = 0;
   
  unsigned CellSize = nel*blockElementNum;
  Cell.resize(CellSize);
  CellType.resize(CellSize);
   
  basis* pt[1];
  hex hex_ = hex();
  pt[0]=&hex_;
  unsigned hexElemDof=pt[0]->ElementDOF(lagrangeP);
  
  for(unsigned kel=0; kel < nel; kel++){

    for(unsigned i=0;i<blockElementNum;i++){
      Cell[i+elemOffset].resize(nVertices);
 
	height=i/pow(lagrangeP,2);  row = (i-pow(lagrangeP,2)*height)/lagrangeP; col=i%lagrangeP;
	
	Cell[i+elemOffset][0]=height*pow(lagrangeP+1,2)+row*(lagrangeP+1)+col+nodeOffset;
	Cell[i+elemOffset][1]=height*pow(lagrangeP+1,2)+row*(lagrangeP+1)+col+1+nodeOffset;
	Cell[i+elemOffset][2]=height*pow(lagrangeP+1,2)+(row+1)*(lagrangeP+1)+col+1+nodeOffset;
	Cell[i+elemOffset][3]=height*pow(lagrangeP+1,2)+(row+1)*(lagrangeP+1)+col+nodeOffset;
	Cell[i+elemOffset][4]=(height+1)*pow(lagrangeP+1,2)+row*(lagrangeP+1)+col+nodeOffset;
	Cell[i+elemOffset][5]=(height+1)*pow(lagrangeP+1,2)+row*(lagrangeP+1)+col+1+nodeOffset;
	Cell[i+elemOffset][6]=(height+1)*pow(lagrangeP+1,2)+(row+1)*(lagrangeP+1)+col+1+nodeOffset;
	Cell[i+elemOffset][7]=(height+1)*pow(lagrangeP+1,2)+(row+1)*(lagrangeP+1)+col+nodeOffset;

      CellType[i+elemOffset]=12;
      
    }
    elemOffset += blockElementNum;
    nodeOffset += hexElemDof;
  }
   
  return nel*blockElementNum*9; // 9 = number of vertices of a linear element +1
}


void  MultiGrid::Print_vtk_ASCII(const unsigned& time, const int type){ //type is 0-based

  int count;
  if(DIMENSION==1){
    GenerateLagrangeSolution2D(type+1);
    count = LagrangeHRefinement2D();
  }
  else{
    GenerateLagrangeSolution3D(type+1);
    count = LagrangeHRefinement3D();
  }
  
  char* filename= new char[60];
  sprintf(filename,"../mesh-generator/output/mesh.%d.%d.vtk",time,type);
  std::ofstream fout;
  fout.open(filename);
  if (!fout) {
    cout << "Output mesh file "<<filename<<" cannot be opened.\n";
    exit(0);
  }

  /* head */
  fout<<"# vtk DataFile Version 3.0\nAMR mesh\nASCII\n";
  fout<<"DATASET UNSTRUCTURED_GRID\n";
  
  /* nodes */
  fout<<"POINTS "<<LNode.size()<<" double\n";
  for(unsigned i=0;i<LNode.size();i++){
    fout<<LNode[i][0]<<" "<<LNode[i][1]<<" "<<LNode[i][2]<<endl; 
  }

  /* cells */
  fout<<"\nCELLS "<<Cell.size()<<" "<<count<<endl; // counter: total number of integers required to represent the list. 
  
  for(unsigned i=0;i<Cell.size();i++){
    fout<<Cell[i].size()<<" ";
    for(unsigned j=0;j<Cell[i].size();j++)
      fout<<Cell[i][j]<<" ";
    fout<<endl;
  }
  
  fout<<"CELL_TYPES "<<CellType.size()<<endl;

  for(unsigned i=0;i<CellType.size();i++)
    fout<<CellType[i]<<endl;
  
  /* solution */
  fout<<"\nPOINT_DATA "<<LNode.size()<<endl;
  for(unsigned i=0;i<LSol.size();i++){
    fout<<"SCALARS Sol"<<SolName[i]<<" double 1"<<endl<<"LOOKUP_TABLE default\n";
    for(unsigned j=0;j<LSol[i].size();j++){
      fout<<LSol[i][j]<<endl;
    } 
    fout<<endl; 
  }

  //   /* boundary condition */
  //   for(unsigned i=0;i<Sol.size();i++){
  //     if(BC[i].size()>=nvt){
  //       fout<<"SCALARS BC"<<SolName[i]<<" int 1"<<endl<<"LOOKUP_TABLE default\n";
  //       for(unsigned j=0;j<nvt;j++){
  // 	fout<<BC[i][j]<<endl;
  //       } 
  //       fout<<endl;
  //     }
  //   }
 
  fout.close();
  delete [] filename;
}
//WHAT Eugenio told about stopping criteria
// Generally I do not have a tolerance on a single level,
// 
// I sweep only once in the post-smoothing and 1 or 0 in the pre-smoothing.
// 
// What I check are the eps at the upper level, at the end of each cycle,
// if in l_infity norm of eps<10^4 I stop, why eps?
// because this is the difference between U_old and U,
// if this does not increase of a reasonable quantity it means the solution as converged.
// 
// If you want to be more precise you can check the relative
// ||eps||_{l_infty}/||sol||_{l_infty}<0.0001 or whatever small number you choose.
// 
// I also check that the residual norm of res=Ax-b is small,
// but this is less indicative than before norm, especially for non-linear problem.