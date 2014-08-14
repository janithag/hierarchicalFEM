#include "main.hpp"
#include "lsysPDE.hpp"
#include "elem_type.hpp"

PetscInt* lsysPDE::node;
double lsysPDE::A[27*27];
double lsysPDE::vx[3][27];

double lsysPDE::Weight;

double lsysPDE::phi1[27];
double lsysPDE::gradphi1[27][3];
double lsysPDE::Weight1;
 
double lsysPDE::_beta=0.;
double lsysPDE::_gamma=1.;
double lsysPDE::_alpha=10000;
double lsysPDE::_Td=0.2;
 
const unsigned END_IND[3]={0,1,3};

lsysPDE::lsysPDE(const char infile[], vector < vector<double> > &vt, unsigned order): 
  mesh(infile, vt){
  p=order;
  SetSign();
};

lsysPDE::lsysPDE(const unsigned& igrid, elem* elc):
  mesh(igrid,elc){
  p=igrid;
  SetSign();
}

void lsysPDE::AddSolutionVector(const char name[], const unsigned order){
  unsigned n = Res.size();
  Res.resize(n+1u);
  Eps.resize(n+1u);
  SolType.resize(n+1u);
  SolName.resize(n+1u);

  SolType[n]=order;
  SolName[n]=new char[8];
  strcpy(SolName[n],name);
}

int lsysPDE::ResizeSolutionVector(const char name[]){
 
  unsigned i = GetIndex(name);
  
  ierr=VecCreate(PETSC_COMM_SELF,&Res[i]); CHKERRQ(ierr);
  ierr=VecSetSizes(Res[i],PETSC_DECIDE,el->NumOfDof(SolType[i])); CHKERRQ(ierr);  
  ierr=VecSetFromOptions(Res[i]); CHKERRQ(ierr);
  ierr=VecDuplicate(Res[i],&Eps[i]); CHKERRQ(ierr);
  
  return 1;
}

/* Returns all the modes for a specified node of a facet */
void lsysPDE::NodeGlobalDofOfModes(const unsigned iNode, vector <unsigned>& iMode) const{ // p: 0-based
   unsigned count=0;

  if ((iNode<MeshNodeNum(0))||(p==0))return; //vertex;
 
  else if (iNode<MeshNodeNum(1)){ //edge
    iMode.resize(p);
    iMode[0]=iNode;
    for(unsigned i=1;i<p;i++)
      iMode[i]=el->NumOfDof(i) - el->MeshVertexNum() + iNode;
  }
 
  else if (iNode<MeshNodeNum(2)){ //quadrilateral face
    if(p<3) return;
    iMode.resize((p-1)*(p-2)/2);
    for(unsigned i=0;i<p-2;i++)
      for(unsigned j=0;j<i+1;j++){
	iMode[count]=el->NumOfDof(i+2) - el->MeshVertexNum() + j*(el->MeshQuadFaceNum()+el->MeshTriFaceNum()) + iNode;
	count++;
      }
  }
  
  else if (iNode<MeshNodeNum(3)){ //triangle face
    if(p<2)return;
    iMode.resize(p*(p-1)/2);
    for(unsigned i=0;i<p-1;i++)
      for(unsigned j=0;j<i+1;j++){
        iMode[count]=el->NumOfDof(i+1) - el->MeshVertexNum() + j*(el->MeshQuadFaceNum() + el->MeshTriFaceNum()) - (j==i)*el->MeshQuadFaceNum() + iNode;
	count++;
      }
  }
  
  else
    throw runtime_error("Invalid node index for a facet in lsysPDE::NodeGlobalDofOfModes()"); 
  
}

void lsysPDE::SetSign(){
  unsigned ielt;
  basish* pt[6];
  pt[0] = new hexh(p+1,0);
  pt[1] = new teth(p+1,0);
  pt[2] = new wedgeh(p+1,0);
  pt[3] = new quadh(p+1);
  pt[4] = new trih(p+1);
  pt[5] = new lineh(p+1);

  sign.resize(MeshElementNum());
  type.resize(MeshElementNum());
  
  for(unsigned iel=0;iel<el->MeshElementNum();iel++){
    ielt=el->ElementType(iel);
    
    sign[iel].resize(pt[ielt]->dimOfBasis());
    pt[ielt]->eval_sign(el->o_(iel), sign[iel]);
    
    type[iel].resize(pt[ielt]->dimOfBasis());
    pt[ielt]->eval_type(el->type_(iel), type[iel]);
  }

//cout<<pt[0]->dimOfBasis()<<endl;
// if(p==3)
//   for(unsigned iel=0;iel<1;iel++)
//     for(int i=0;i<pt[2]->dimOfBasis();i++)
//       cout<<i<<" "<<sign[iel][i]<<"\n";
    
  for(int i=0;i<6;i++)
    delete pt[i];
 
}
unsigned lsysPDE::GetIndex(const char name[]){
  unsigned index=0;
  while(strcmp(SolName[index],name)){
    index++;
    if(index==Res.size()) {
      cout<<"error! invalid name entry GetIndex(...)"<<endl;
      exit(0);
    }
  }
  return index;
}

int lsysPDE::InitMultigrid(const vector <unsigned> &MGIndex){
   
  KIndex.resize(MGIndex.size()+1u);
  KIndex[0]=0;

  for(unsigned i=1;i<KIndex.size();i++)
    KIndex[i] = KIndex[i-1] + el->NumOfDof(SolType[MGIndex[i-1]]);
  
  PetscInt RESsize = KIndex[ KIndex.size()-1];
  
  ierr=VecCreate(PETSC_COMM_SELF,&RES); CHKERRQ(ierr);
  ierr=VecSetSizes(RES,PETSC_DECIDE,RESsize); CHKERRQ(ierr);  
  ierr=VecSetFromOptions(RES); CHKERRQ(ierr);

  ierr=VecDuplicate(RES,&RESC); CHKERRQ(ierr);
  ierr=VecDuplicate(RES,&EPS); CHKERRQ(ierr);
  ierr=VecDuplicate(RES,&EPSC); CHKERRQ(ierr);
  
  return 1;
}

int lsysPDE::SetResZero(const vector <unsigned> &MGIndex){
  ierr=VecSet(RES,0.);  CHKERRQ(ierr);
  return 1;
}

int lsysPDE::SetEpsZero(const vector <unsigned> &MGIndex){
  ierr=VecSet(EPS,0.); CHKERRQ(ierr);
  ierr=VecSet(EPSC,0.); CHKERRQ(ierr);

  return 1;
}

int lsysPDE::SumEpsCToEps(const vector <unsigned> &MGIndex){
  ierr=VecAXPBY(EPS,1,1,EPSC); CHKERRQ(ierr);
  return 1;
}

int lsysPDE::UpdateResidual(){
    
  ierr = MatMult(K,EPSC,RESC); CHKERRQ(ierr);
  ierr = VecAXPBY(RES,-1.,1.,RESC);CHKERRQ(ierr);
   
  return 1;
}


int lsysPDE::FindResMax(const vector < vector <unsigned short> > &BC, const vector <unsigned> &MGIndex){

  PetscScalar *R[1], *E[1]; 
  ResInf.resize(MGIndex.size());
  
  ierr = VecGetArray(RES,R); CHKERRQ(ierr);
  ierr = VecGetArray(EPS,E); CHKERRQ(ierr);
  
  for(unsigned k=0;k<MGIndex.size();k++){ 
    unsigned indexSol=MGIndex[k];
    
    for(PetscInt i=KIndex[k];i<KIndex[k+1];i++){
      PetscInt inode=i-KIndex[k];
      if(BC[indexSol][inode]>p){
	ierr = VecSetValues(Res[indexSol],1,&inode,&R[0][i],INSERT_VALUES);CHKERRQ(ierr);
      }
      else{
	PetscScalar zero=0.;
        ierr = VecSetValues(Res[indexSol],1,&inode,&zero,INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = VecSetValues(Eps[indexSol],1,&inode,&E[0][i],INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(Res[indexSol]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res[indexSol]);CHKERRQ(ierr);
  }
  
  ierr = VecRestoreArray(RES,R); CHKERRQ(ierr);
  ierr = VecRestoreArray(EPS,E); CHKERRQ(ierr);
  
  
  for(unsigned k=0;k<MGIndex.size();k++){
    unsigned indexSol=MGIndex[k];
    double ResMax=0.;
    ierr=VecNorm(Res[indexSol],NORM_INFINITY,&ResMax);   CHKERRQ(ierr);
    ResInf[k]=ResMax;
    //cout<<"level="<<p<<"\tResMax"<<SolName[indexSol]<<"="<<ResMax<<endl;
  }
  
  return 0;
}

int lsysPDE::SumEpsToSol(vector <vector <double> >&Sol, const vector <unsigned> &MGIndex){
  
  PetscScalar *E[1]; 
  for(unsigned k=0;k<MGIndex.size();k++){
    unsigned kk=MGIndex[k];
    ierr = VecGetArray(Eps[kk],E); CHKERRQ(ierr);
    for(unsigned i=0;i<el->NumOfDof(SolType[kk]);i++){
      Sol[kk][i] += E[0][i];
    } 
    ierr = VecRestoreArray(Eps[kk],E); CHKERRQ(ierr); 
  }

  return 1;
}

int  lsysPDE::DeallocateMatrix(const vector <unsigned> &MGIndex){
  ierr=MatDestroy(&K);CHKERRQ(ierr);

   if(p>0){
     ierr=MatDestroy(&P);CHKERRQ(ierr);
     ierr=MatDestroy(&R);CHKERRQ(ierr);
  } 

   ierr=VecDestroy(&RES);CHKERRQ(ierr);
   ierr=VecDestroy(&RESC);CHKERRQ(ierr);
   
   ierr=VecDestroy(&EPSC);CHKERRQ(ierr);
   ierr=VecDestroy(&EPS);CHKERRQ(ierr);
   
   return 1;
}

int lsysPDE::FreeSolutionVectors(){
   for(unsigned i=0;i<Res.size();i++){
      ierr=VecDestroy(&Res[i]); CHKERRQ(ierr);
      ierr=VecDestroy(&Eps[i]); CHKERRQ(ierr);
   }
   return 1;
}

/* We assemble the stiffness matrix for finest level in this function */
int lsysPDE::AssembleMatrix(vector < vector <const elem_type*> >  type_elem, vector <vector <double> > &vt, const unsigned& gridn, const char type[]){

  clock_t AssemblyTime = 0; 
  
  unsigned N=1; //number of variables
  unsigned order[N];
  unsigned elemDOF[N];
  double* phi;
  double** gradphi;
  
  order[0] = SolType[GetIndex("U")];
 // order[1] = SolType[GetIndex("V")];
  //order[2] = SolType[GetIndex("W")];
  unsigned maxOrder = *max_element(order,order+N-1);

  vertex.resize(3);
  
  /* size of global stiffness matrix */
  PetscInt Ksize = KIndex[KIndex.size()-1u]; 
  PetscErrorCode ierr;

  unsigned maxDOF = type_elem[0][maxOrder]->ElementDOF();
  for(unsigned i=1;i<6;i++)
    if (type_elem[i][maxOrder]->ElementDOF() > maxDOF)
      maxDOF = type_elem[i][maxOrder]->ElementDOF();
  
  /* Global stiffness matrix initialization */
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,Ksize,Ksize,20*maxDOF,PETSC_NULL,&K); CHKERRQ(ierr); //TODO number of nz per row 
  ierr = MatSetFromOptions(K); CHKERRQ(ierr);
  ierr = MatZeroEntries(K); CHKERRQ(ierr); 
 
  double B[N][N][maxDOF*maxDOF];
  node = new PetscInt[maxDOF];
  phi = new double[maxDOF];
  double* X=new double[3];
  gradphi = new double*[maxDOF];
  for(unsigned c=0; c<maxDOF; c++)
    gradphi[c]= new double[3];  
  
  int* s;
  int* t;
  //cout<<"signs"<<endl;
  for(unsigned kel=0;kel<MeshElementNum();kel++){
//   for(int i=0;i<type_elem[2][p]->ElementDOF();i++)
//       cout<<i<<": "<<sign_(kel)[i]<<"\t";
//   cout<<endl;
    short unsigned kelt = el->ElementType(kel);
  
    for(unsigned k=0;k<N;k++){
      elemDOF[k] = type_elem[kelt][order[k]]->ElementDOF(); 
      memset(B[k][k],0.,maxDOF*maxDOF*sizeof(double));
    }

    vertex.resize(el->ElementNodeNum(kel,0));

    /* vertices of the element */
    int vIndex[el->ElementNodeNum(kel,0)];
    for(unsigned i=0;i<el->ElementNodeNum(kel,0);i++){
      vertex[i].resize(3);
      vIndex[i] = el->GlobalNodeIndex(kel,i)-1u;	
      vertex[i][0]=vt[0][vIndex[i]];
      vertex[i][1]=vt[1][vIndex[i]];
      vertex[i][2]=vt[2][vIndex[i]];
    }
    s=sign_(kel);t=type_(kel);
    /* Build element stiffness matrix */
    if(p==gridn || !el->GetRefinedElementIndex(kel) ){
    for(unsigned ig=0; ig < type_elem[kelt][maxOrder]->NumOfGaussPts(); ig++){
      (type_elem[kelt][maxOrder]->*type_elem[kelt][maxOrder]->spatialToRefPtr)(vertex,ig,s,Weight,phi,gradphi,X,t);
    
      for(unsigned k=0;k<N;k++)

	/* Laplace operator */ //tested!
	for(unsigned i=0; i<elemDOF[k]; i++)
	  for(unsigned j=0; j<elemDOF[k]; j++) 
	    B[k][k][i*elemDOF[k]+j] += (gradphi[i][0]*gradphi[j][0]+gradphi[i][1]*gradphi[j][1]+gradphi[i][2]*gradphi[j][2])*Weight;
    }
    }
    /* some additional code to see the sparse structure */
    // for(unsigned k=0;k<N;k++)
    //   for(unsigned i=0; i<elemDOF[k]; i++)
    // 	for(unsigned j=0; j<elemDOF[k]; j++) 
    // 	  if(B[k][k][i*elemDOF[k]+j]<pow(10,-14)) B[k][k][i*elemDOF[k]+j]=0;
    
    /* assemble global stiffness matrix */
    for(unsigned k=0;k<N;k++){
      for(unsigned i=0;i<elemDOF[k];i++) 
      	node[i] = el->globalDOF_(kel)[i]+KIndex[k];
        
      ierr = MatSetValuesBlocked(K,elemDOF[k],node,elemDOF[k],node,B[k][k],ADD_VALUES);CHKERRQ(ierr);
    }
  }

  delete[] node; 
  delete[] phi;
  delete[] X;
  
  for(unsigned c=0;c<maxDOF;c++)
    delete[] gradphi[c];
  delete[] gradphi;

  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Visualize the matrix */
 
//    PetscViewer viewer;
//    PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);
//    MatView(K,viewer);
//    double f;
//    cin>>f;
 
//MatView(K,PETSC_VIEWER_STDOUT_WORLD);

  /* Computational info */
//   cout<<"Level = "<<p<<endl;
//   cout<<"ASSEMBLY time="<<static_cast<double>(AssemblyTime)/CLOCKS_PER_SEC<<endl;  
  
  return 1;
  
}

/* We assemble the residual for the finest level in this function */
int lsysPDE::AssembleResidual(vector < vector <const elem_type*> > type_elem, vector <vector <double> > &vt, const vector <vector <double> > &Sol, const unsigned& gridn, const char type[], const unsigned& pOrder){
  
  clock_t AssemblyTime = 0; 
  clock_t start_time, end_time;
  start_time = clock(); 
  unsigned N=1; //number of variables
  unsigned order[N];
  unsigned elemDOF[N];
  unsigned solIndex[N];
  double sum[3];
  double* phi;
  double** gradphi;
  
  order[0] = SolType[GetIndex("U")];
  //order[1] = SolType[GetIndex("V")];
  //order[2] = SolType[GetIndex("W")]; 
  unsigned maxOrder = *max_element(order,order+N-1);
  //PetscReal value;
  double f; // Poisson's equation
  double* X;

  unsigned maxDOF = type_elem[0][maxOrder]->ElementDOF();

  PetscScalar b[N][maxDOF];
  solIndex[0] = GetIndex("U");
  //solIndex[1] = GetIndex("V");
  //solIndex[2] = GetIndex("W");

  X = new double[3];
  node = new PetscInt[maxDOF];
  phi = new double[maxDOF];

  gradphi = new double*[maxDOF];
  for(unsigned c=0; c<maxDOF; c++)
    gradphi[c]= new double[3];
  int* s;
  int* t;
  
  for (unsigned kel=0;kel<MeshElementNum();kel++){
//     for(int i=0;i<10;i++)
// 	cout<<kel<<" "<<sign_(kel)[i]<<endl;
    f=1;
   if(p==gridn || !el->GetRefinedElementIndex(kel) ){
    short unsigned kelt = el->ElementType(kel);
    
    unsigned numOfVertices = el->ElementNodeNum(kel,0);
    vertex.resize(numOfVertices);

    /* vertices of the element */
    int vIndex[numOfVertices];
    double x=0.;
    for(unsigned i=0;i<numOfVertices;i++){
      vertex[i].resize(3);
      vIndex[i] = el->GlobalNodeIndex(kel,i)-1u;	
      vertex[i][0]=vt[0][vIndex[i]];
      vertex[i][1]=vt[1][vIndex[i]];
      vertex[i][2]=vt[2][vIndex[i]];
      x += vt[0][vIndex[i]];
    }
    //bool flag = (p==pOrder && x/numOfVertices<0.5);
    for(unsigned k=0;k<N;k++) {
      elemDOF[k] = type_elem[kelt][order[k]]->ElementDOF();
      memset(b[k],0.,maxDOF*sizeof(double)); 
    }
    //if(p==pOrder && x/numOfVertices<0.5) continue;
    s=sign_(kel);t=type_(kel);
  
    /* Build element residual vector for Poisson's equation */ //TODO Changed the order of loops. recheck this. 
    for(unsigned k=0;k<N;k++){
 
    for(unsigned ig=0; ig < type_elem[kelt][maxOrder]->NumOfGaussPts(); ig++){
      (type_elem[kelt][maxOrder]->*type_elem[kelt][maxOrder]->spatialToRefPtr)(vertex,ig,s,Weight,phi,gradphi,X,t);    
      

	for(unsigned n=0;n<3;n++){
	  sum[n]=0.;
	  for (unsigned j=0;j<elemDOF[k];j++)
	    sum[n] += gradphi[j][n]*Sol[solIndex[k]][el->globalDOF_(kel)[j]];
	}
        //f=8*PI*X[0]*sin(4*PI*X[1])*(2*PI*(1+4*pow(X[0],2))*sin(4*PI*pow(X[0],2))-3*cos(4*PI*pow(X[0],2)));
	//f=-2*pow(X[0],2)+2*X[0]-2*pow(X[1],2)+2*X[1]; // uEx in polynomial form for unit square 2D
	//f=8*pow(PI,2)*sin(2*PI*X[0])*sin(2*PI*X[1]); // uEx=sin(2*pi*x)sin(2*pi*y) for unit square 2D
	
	f=12*pow(PI,2)*sin(2*PI*X[0])*sin(2*PI*X[1])*sin(2*PI*X[2]);  //3D unit cube trig
	//f=2*X[1]*X[2]*(1-X[1])*(1-X[2])+2*X[2]*X[0]*(1-X[2])*(1-X[0])+2*X[0]*X[1]*(1-X[0])*(1-X[1]);//polynomial 3D unit cube
	//f=8*PI*X[0]*sin(2*PI*X[1])*sin(2*PI*X[2])*(PI*(1+8*pow(X[0],2))*sin(4*PI*pow(X[0],2))-3*cos(4*PI*pow(X[0],2)));
	
	//f=2*X[2]*(pow(X[1],2)-X[1])-2*(pow(X[2],2)*X[0]-pow(X[0],2)*X[2]-pow(X[2],2)+X[2]*X[0])-2*(pow(X[1],2)-X[1])*(X[0]-1);//unit wedge
	 /* Mixed boundary conditions for unit wedge */
	
	 //f=2*X[2]+2-2*X[0]; // Neumann on all tri faces 4,5
	 //f=2;//Neumann on all quad faces
	 
	 for(unsigned i=0; i<elemDOF[k]; i++)	 
	  b[k][i] += (f*phi[i]-gradphi[i][0]*sum[0]-gradphi[i][1]*sum[1]-gradphi[i][2]*sum[2])*Weight;
	
      }
    }
    
    /* Assemble global residual vector */
    for(unsigned k=0;k<N;k++){
      
      for(unsigned i=0;i<elemDOF[k];i++)
	node[i] = el->globalDOF_(kel)[i]+KIndex[k]; 
	
      ierr = VecSetValues(RES,elemDOF[k],node,b[k],ADD_VALUES);CHKERRQ(ierr);
    }
  
  }
}
  delete[] X;
  delete[] node; 
  delete[] phi;

  for(unsigned c=0;c<maxDOF;c++)
    delete[] gradphi[c];
  delete[] gradphi;

  ierr = VecAssemblyBegin(RES);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(RES);CHKERRQ(ierr);

  end_time = clock();
  AssemblyTime += (end_time-start_time);
  //cout<<"RHS ASSEMBLY time="<<static_cast<double>(AssemblyTime)/CLOCKS_PER_SEC<<endl;

 /* Visualize the vector */
 
//  cout<<"RES in assembly function \n";
//  VecView(RES,PETSC_VIEWER_STDOUT_WORLD);
 
  return 1;
}

int lsysPDE::SpaceDecompositionSmoother(const vector  <vector <unsigned short> >& BC, const vector <unsigned>& MGIndex, const int& DofstoSolve){

  PetscErrorCode ierr;
  //clock_t SearchTime=0, AssemblyTime=0, SolveTime0=0, UpdateTime=0; 
  //clock_t start_time, end_time;
  
  static PetscInt A_max=0;

  unsigned nvt = KIndex[KIndex.size()-1u];
  unsigned nel = MeshElementNum();
    
  int p_max=p;
 
  double tolerance=1.0e-6;
    
  unsigned IndexaSize=KIndex[KIndex.size()-1];
  unsigned IndexbSize=KIndex[KIndex.size()-1];
  std::vector < PetscInt > indexai; 
  indexai.resize(IndexaSize);
  std::vector < PetscInt > indexbi;
  indexbi.resize(IndexbSize);
  

  /* Search for elements and modes to be solved/updated for the block */
  for (int level=0; level<p_max+1; level++){ // loop on the spaces V^l
    
    //start_time = clock();
    PetscInt counta=0;
    PetscInt countb=0;
    
    for(unsigned indexSol=0; indexSol<KIndex.size()-1u; indexSol++){
      for(unsigned idof = KIndex[indexSol]; idof<KIndex[indexSol+1u]; idof++){
	unsigned jnode=idof-KIndex[indexSol];
	
	if( (level==0 || jnode>=el->NumOfDof(level-1))  && jnode<el->NumOfDof(level) && p<BC[MGIndex[indexSol]][jnode] ){// if in the fine mesh Dirichlet
	    indexai[counta] = idof;
	    counta++;
	}
	
	else if (0<BC[MGIndex[indexSol]][jnode]){// if not Dirichlet
	  indexbi[countb] = idof;
	  countb++;
	}
      }
    }
  
    if(0==counta) continue;
    
    PetscInt PBsize=counta;     
    PetscInt PAmBsize=countb;
    
// /* test */    
// cout<<"\np: "<<p;
// cout<<"\nnodes to solve"<<endl;
// cout<<"size: "<<PBsize<<endl;
// // for(int i=0;i<PBsize;i++){
// //   cout<<indexai[i]<<" ";
// // }
// cout<<"\nnodes to solve or update"<<endl;
// cout<<"size: "<<counterb<<endl;
// for(int i=0;i<counterb;i++)
//   cout<<indexbi[i]<<" ";
 
    if(counta>A_max) A_max=counta;
    
    //end_time=clock();
    //SearchTime+=(end_time-start_time); 
    
    // ***************** ASSEMBLY *******************
    //start_time=clock();
    Vec Pw, Pf, Pr;     /* exact solution, RHS*/
    Mat PA    ;         /* linear system matrix */
    KSP ksp;            /* linear solver context */
    PC  pc;             /* preconditioner context */
    
    // generate IS
    IS isA;
    ierr=ISCreateGeneral(PETSC_COMM_SELF,PBsize,&indexai[0],PETSC_USE_POINTER,&isA);CHKERRQ(ierr);
    ierr=ISSort(isA); 
    
    // Solution Vector Pw
    ierr = VecCreate(PETSC_COMM_SELF,&Pw); CHKERRQ(ierr);
    ierr = VecSetSizes(Pw,PETSC_DECIDE,PBsize); CHKERRQ(ierr);
    ierr = VecSetFromOptions(Pw); CHKERRQ(ierr);

    // RHS Pf
    ierr = VecDuplicate(Pw,&Pf); CHKERRQ(ierr);

    PetscInt *indl=new PetscInt[PBsize];
    for(int i=0;i<PBsize;i++) indl[i]=i;
    PetscScalar* y=new PetscScalar [PBsize];
    const PetscInt *ind[1];
    
    ierr = ISGetIndices(isA,ind); CHKERRQ(ierr);
    ierr = VecGetValues(RES,PBsize,ind[0],y);
    ierr = VecSetValues(Pf,PBsize,indl,y,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(Pf);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Pf);CHKERRQ(ierr);
    
    ierr = ISRestoreIndices(isA,ind); CHKERRQ(ierr);
    
    delete [] y;
    delete [] indl;
    
    
    // RHS Pr
    ierr = VecDuplicate(Pw,&Pr); CHKERRQ(ierr);
    
    //Matrix PA
    ierr = MatGetSubMatrix(K,isA,isA,MAT_INITIAL_MATRIX,&PA); CHKERRQ(ierr);
    
    //end_time=clock();
    //AssemblyTime+=(end_time-start_time);
    
   
    // ***************** SOLVE *******************
    //start_time=clock();
    //set KSP and solve    
      
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,PA,PA,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
     
    if(p==0){
      ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    }
    else{
      ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
    }
    ierr = KSPSetTolerances(ksp,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    
    ierr = KSPSolve(ksp,Pf,Pw);CHKERRQ(ierr); 
    //end_time=clock();
    //SolveTime0 += (end_time-start_time);
    
    // ***************** UPDATE ****************** 
    //start_time=clock();
    
    //update Residual and solution for the nodes that are solved for
    ierr = MatMult(PA,Pw,Pr); CHKERRQ(ierr);
    
    ierr = VecScale (Pr, -1.); CHKERRQ(ierr);
    PetscScalar *R[1]; 
    ierr = VecGetArray(Pr,R); CHKERRQ(ierr);
    
    PetscScalar *W[1];
    ierr = VecGetArray(Pw,W); CHKERRQ(ierr);
    
    ierr = ISGetIndices(isA,ind); CHKERRQ(ierr);
    
    ierr = VecSetValues(RES,PBsize,ind[0],R[0],ADD_VALUES);  CHKERRQ(ierr);
    ierr = VecAssemblyBegin(RES);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(RES);CHKERRQ(ierr);
    
    ierr = VecSetValues(EPS,PBsize,ind[0],W[0],ADD_VALUES);  CHKERRQ(ierr);
    ierr = VecAssemblyBegin(EPS);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(EPS);CHKERRQ(ierr);

    ierr = ISRestoreIndices(isA,ind); CHKERRQ(ierr);
    ierr = VecRestoreArray(Pr,R); CHKERRQ(ierr);
      
    ierr = VecRestoreArray(Pw,W); CHKERRQ(ierr);
    
    if(PAmBsize){ //update residual for the nodes that are not solved for
      IS isB;
      ierr=ISCreateGeneral(PETSC_COMM_SELF,PAmBsize,&indexbi[0],PETSC_USE_POINTER,&isB);CHKERRQ(ierr);
      ierr=ISSort(isB);
      
      Vec Ps;
      ierr = VecCreate(PETSC_COMM_SELF,&Ps); CHKERRQ(ierr);
      ierr = VecSetSizes(Ps,PETSC_DECIDE,PAmBsize); CHKERRQ(ierr);
      ierr = VecSetFromOptions(Ps); CHKERRQ(ierr);
      
      Mat PB;
      ierr = MatGetSubMatrix(K,isB,isA,MAT_INITIAL_MATRIX,&PB); CHKERRQ(ierr);
      ierr = MatMult(PB,Pw,Ps); CHKERRQ(ierr);

      ierr = VecScale (Ps, -1.); CHKERRQ(ierr);
      PetscScalar *S[1]; 
      ierr = VecGetArray(Ps,S); CHKERRQ(ierr);
      ierr = ISGetIndices(isB,ind); CHKERRQ(ierr);

      ierr = VecSetValues(RES,PAmBsize,ind[0],S[0],ADD_VALUES);  CHKERRQ(ierr);
      ierr = VecAssemblyBegin(RES);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(RES);CHKERRQ(ierr);

      ierr = ISRestoreIndices(isB,ind); CHKERRQ(ierr);
      ierr = VecRestoreArray(Ps,S); CHKERRQ(ierr);
      ierr = VecDestroy(&Ps);CHKERRQ(ierr);
      ierr = MatDestroy(&PB);CHKERRQ(ierr); 
      ierr = ISDestroy(&isB);CHKERRQ(ierr); 
    }
        
    ierr = VecDestroy(&Pw);CHKERRQ(ierr);   
    ierr = VecDestroy(&Pf);CHKERRQ(ierr);  
    ierr = VecDestroy(&Pr);CHKERRQ(ierr);
    ierr = MatDestroy(&PA);CHKERRQ(ierr); 
    ierr = ISDestroy(&isA);CHKERRQ(ierr);
    
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr); 
    
    //end_time=clock();
    //UpdateTime+=(end_time-start_time);
      
  }
  //cout<<"EPS inside Vanka"<<endl;
  //VecView(EPS,PETSC_VIEWER_STDOUT_WORLD);
  // *** Computational info ***
//   cout<<"Grid="<<p<<endl;
//   cout<<"SEARCH   time="<<static_cast<double>(SearchTime)/CLOCKS_PER_SEC<<endl;
//   cout<<"ASSEMBLY time="<<static_cast<double>(AssemblyTime)/CLOCKS_PER_SEC<<endl;
//   cout<<"SOLVE 0  time="<<static_cast<double>(SolveTime0)/CLOCKS_PER_SEC<<endl;
//   cout<<"UPDATE   time="<<static_cast<double>(UpdateTime)/CLOCKS_PER_SEC<<endl;
//   cout<<"A_MAX="<<A_max<<" C_MAX="<<C_max<<endl;
  //cout<<"number of blocks: "<<blocks<<endl;
  return 1;
}

int lsysPDE::VankaPETSCSmoother(const vector  <vector <unsigned short> >& BC, const vector <unsigned>& MGIndex, const int& DofstoSolve){

  PetscErrorCode ierr;
  //clock_t SearchTime=0, AssemblyTime=0, SolveTime0=0, UpdateTime=0; 
  //clock_t start_time, end_time;
  vector <short unsigned> indexa;
  vector <short unsigned> indexb;
  vector <short unsigned> indexc; 
  vector <short unsigned> indexd; 
  
  vector <PetscInt> indexai;
  vector <PetscInt> indexbi;
  vector <PetscInt> indexci;
  vector <PetscInt> indexdi;
  
  PetscInt A_max=0;
  PetscInt C_max=0;

  unsigned nvt = KIndex[KIndex.size()-1u];
  unsigned nel = MeshElementNum();
  
  if(indexa.size()<nvt || indexc.size()<nel){ 
    indexa.resize(nvt); // modes to solve
    indexb.resize(nvt); // modes to solve or update
    indexc.resize(nel); // elements to solve
    indexd.resize(nel); // elements to solve or update
  }
  
  unsigned IndexaSize = 50000;
  unsigned IndexbSize = 50000;
  unsigned IndexcSize = 10000;
  unsigned IndexdSize = 10000;
  
  if(indexai.size() < IndexaSize){
    indexai.resize(IndexaSize);
    indexbi.resize(IndexbSize);
    indexci.resize(IndexcSize); 
    indexdi.resize(IndexdSize);
  }
  
  for(unsigned i=0;i<nvt;i++){
    indexa[i] = IndexaSize;
    indexb[i] = IndexbSize;
  }

  for(unsigned i=0;i<nel;i++){
    indexc[i] = IndexcSize;
    indexd[i] = IndexdSize;
  }
   
  /* Start Vanka Block  */
  int Dgel; // number of elements in the block in the level
  int grid = p;
 
  int blocks=DofstoSolve/9000+1;//10000
  Dgel = el->MeshElementNum()/blocks;  
  
  double tolerance=1.0e-6;
  //cout<<"grid="<<grid<<" tolerance="<<tolerance<<endl;

  /* Search for elements and modes to be solved/updated for the block */
  for (int gel=0; gel>=0 && gel<static_cast <int> (nel); gel += Dgel){ // loop on all blocks 
    //start_time = clock();
    PetscInt Asize=0;
    PetscInt counterb=0;
    PetscInt Csize=0; 
    PetscInt Dsize=0;

    for(int iel=gel; iel<gel+Dgel && iel< static_cast <int> (nel); iel++){ // loop on all elements in each block

      for(unsigned i=0; i<el->ElementNodeNum(iel,0); i++){ // loop on all the vertices of each element
	unsigned inode = el->GlobalNodeIndex(iel,i)-1u; // global dof of the vertex
	unsigned nvei = el->NodeNumOfAdjacentElements(inode); // number of adjacent elements to the vertex
	const unsigned* pt_jel = el->NodeAdjacentElementAddress(inode,0); // address of 0-direction adjacent element to the vertex

	for(unsigned j=0;j<nvei;j++){ // loop on all adjacent elements
	  unsigned jel = *(pt_jel++)-1u; // index of the adjacent element
	  unsigned jelt = el->ElementType(jel);
	  
	  /* add elements to be solved */
	  if(indexc[jel] == IndexcSize){
	     //cout<<"Element: "<<indexci[Csize]<<endl;
	    indexci[Csize] = jel;
	    indexc[jel] = Csize++;

	    /* add nodes to be solved */
	    for(unsigned indexSol=0; indexSol<KIndex.size()-1u; indexSol++){ // loop on all the variables
	      elem_type E(jelt,SolType[MGIndex[indexSol]]+1,1,1);
	      for(unsigned J=0; J<E.ElementDOF(); J++){ // loop on all the modes of the adjacent elements
		unsigned jnode = el->globalDOF_(jel)[J];
		
		if (indexa[jnode+KIndex[indexSol]]==IndexaSize && p<BC[MGIndex[indexSol]][jnode] ){ //TEST p// if in the fine level and if not Dirichlet
		  indexai[Asize] = jnode + KIndex[indexSol]; // add the global dof in indexai
		  indexa[jnode + KIndex[indexSol]] = Asize++; // add the search information for that dof
		}
	      }
	    }

	    for(unsigned J=0; J<el->ElementNodeNum(jel,0); J++){
	      unsigned jnode = el->GlobalNodeIndex(jel,J)-1u;
	      unsigned nvej = el->NodeNumOfAdjacentElements(jnode);
	      const unsigned* pt_kel = el->NodeAdjacentElementAddress(jnode,0);

	      for(unsigned k=0;k<nvej;k++){
		unsigned kel = *(pt_kel++)-1u;
		unsigned kelt = el->ElementType(kel);
		
		/* add elements for all variables to be updated to indexdi*/
		if(indexd[kel] == IndexdSize){
		  indexdi[Dsize] = kel;
		  indexd[kel] = Dsize++;

		  /* add modes for all variables to be updated to indexbi */
		  for(unsigned indexSol=0;indexSol<KIndex.size()-1u;indexSol++){
		    elem_type F(kelt,SolType[MGIndex[indexSol]]+1,1,1);
		    
		    for(unsigned K=0; K<F.ElementDOF(); K++){
			unsigned knode = el->globalDOF_(kel)[K];
			if (indexb[knode+KIndex[indexSol]] == IndexbSize && 0<BC[MGIndex[indexSol]][knode]){// if not Dirichlet
			  indexbi[counterb] = knode + KIndex[indexSol];
			  indexb[knode + KIndex[indexSol]] = counterb++;
			}
		      }
		  }
		}
	      }
	    }  
	  }	  
	}
      }
    } 
  
    if(0==Asize) continue;
    
    PetscInt PBsize=Asize;
    for(PetscInt i=0;i<counterb;i++){
      unsigned jnode=indexbi[i];
      if(indexa[jnode]==IndexaSize){
	indexai[Asize]=jnode;
	indexa[jnode]=Asize++;
      }
      // reinitialize indexb
      indexb[jnode]=IndexbSize;
    }
    
    // re-initialize indices
    for(PetscInt i=0;i<Asize;i++){
      indexa[indexai[i]]=IndexaSize;
    }
    
    for(PetscInt i=0;i<Csize;i++){
      indexc[indexci[i]]=IndexcSize;
    }
    
    for(PetscInt i=0;i<Dsize;i++){
      indexd[indexdi[i]]=IndexdSize;
    }
    
    PetscInt PAmBsize=Asize-PBsize;
    
// /* test */    
// cout<<"\np: "<<p;
// cout<<"\nnodes to solve"<<endl;
// cout<<"size: "<<PBsize<<endl;
// for(int i=0;i<PBsize;i++){
//   cout<<indexai[i]<<" ";
// }
// cout<<"\nnodes to solve or update"<<endl;
// cout<<"size: "<<counterb<<endl;
// for(int i=0;i<counterb;i++)
//   cout<<indexbi[i]<<" ";
 
    if(Asize>A_max) A_max=Asize;
    if(Csize>C_max) C_max=Csize;
    
    //end_time=clock();
    //SearchTime+=(end_time-start_time); 
    // ***************** END NODE/ELEMENT SEARCH *******************
    
    // ***************** ASSEMBLY *******************
    //start_time=clock();
    Vec Pw, Pf, Pr;     /* exact solution, RHS*/
    Mat PA    ;         /* linear system matrix */
    KSP ksp;            /* linear solver context */
    PC  pc;             /* preconditioner context */
    
    // generate IS
    IS isA;
    ierr=ISCreateGeneral(PETSC_COMM_SELF,PBsize,&indexai[0],PETSC_USE_POINTER,&isA);CHKERRQ(ierr);
    ierr=ISSort(isA); 
    
    // Solution Vector Pw
    ierr = VecCreate(PETSC_COMM_SELF,&Pw); CHKERRQ(ierr);
    ierr = VecSetSizes(Pw,PETSC_DECIDE,PBsize); CHKERRQ(ierr);
    ierr = VecSetFromOptions(Pw); CHKERRQ(ierr);

    // RHS Pf
    ierr = VecDuplicate(Pw,&Pf); CHKERRQ(ierr);

    PetscInt *indl=new PetscInt[PBsize];
    for(int i=0;i<PBsize;i++) indl[i]=i;
    PetscScalar* y=new PetscScalar [PBsize];
    const PetscInt *ind[1];
    
    ierr = ISGetIndices(isA,ind); CHKERRQ(ierr);
    ierr = VecGetValues(RES,PBsize,ind[0],y);
    ierr = VecSetValues(Pf,PBsize,indl,y,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(Pf);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Pf);CHKERRQ(ierr);
    
    ierr = ISRestoreIndices(isA,ind); CHKERRQ(ierr);
    
    delete [] y;
    delete [] indl;
    
    
    // RHS Pr
    ierr = VecDuplicate(Pw,&Pr); CHKERRQ(ierr);
    
    //Matrix PA
    ierr = MatGetSubMatrix(K,isA,isA,MAT_INITIAL_MATRIX,&PA); CHKERRQ(ierr);
    
    //end_time=clock();
    //AssemblyTime+=(end_time-start_time);
    
    // ***************** END ASSEMBLY ******************   
    
    // ***************** SOLVE *******************
    //start_time=clock();
    //set KSP and solve    
      
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,PA,PA,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
     
    if(p==0){
      ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    }
    else{
      ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
    }
    ierr = KSPSetTolerances(ksp,tolerance,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    
    ierr = KSPSolve(ksp,Pf,Pw);CHKERRQ(ierr); 
    //end_time=clock();
    //SolveTime0 += (end_time-start_time);
    // ***************** END SOLVE ****************** 
    
    // ***************** START UPDATE ****************** 
    //start_time=clock();
    
    //update Residual and solution for the nodes that are solved for
    ierr = MatMult(PA,Pw,Pr); CHKERRQ(ierr);
    
    ierr = VecScale (Pr, -1.); CHKERRQ(ierr);
    PetscScalar *R[1]; 
    ierr = VecGetArray(Pr,R); CHKERRQ(ierr);
    
    PetscScalar *W[1];
    ierr = VecGetArray(Pw,W); CHKERRQ(ierr);
    
    ierr = ISGetIndices(isA,ind); CHKERRQ(ierr);
    
    ierr = VecSetValues(RES,PBsize,ind[0],R[0],ADD_VALUES);  CHKERRQ(ierr);
    ierr = VecAssemblyBegin(RES);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(RES);CHKERRQ(ierr);
    
    ierr = VecSetValues(EPS,PBsize,ind[0],W[0],ADD_VALUES);  CHKERRQ(ierr);
    ierr = VecAssemblyBegin(EPS);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(EPS);CHKERRQ(ierr);

    ierr = ISRestoreIndices(isA,ind); CHKERRQ(ierr);
    ierr = VecRestoreArray(Pr,R); CHKERRQ(ierr);
      
    ierr = VecRestoreArray(Pw,W); CHKERRQ(ierr);
    
    if(PAmBsize){ //update residual for the nodes that are not solved for
      IS isB;
      ierr=ISCreateGeneral(PETSC_COMM_SELF,PAmBsize,&indexai[PBsize],PETSC_USE_POINTER,&isB);CHKERRQ(ierr);
      ierr=ISSort(isB);
      
      Vec Ps;
      ierr = VecCreate(PETSC_COMM_SELF,&Ps); CHKERRQ(ierr);
      ierr = VecSetSizes(Ps,PETSC_DECIDE,PAmBsize); CHKERRQ(ierr);
      ierr = VecSetFromOptions(Ps); CHKERRQ(ierr);
      
      Mat PB;
      ierr = MatGetSubMatrix(K,isB,isA,MAT_INITIAL_MATRIX,&PB); CHKERRQ(ierr);
      ierr = MatMult(PB,Pw,Ps); CHKERRQ(ierr);

      ierr = VecScale (Ps, -1.); CHKERRQ(ierr);
      PetscScalar *S[1]; 
      ierr = VecGetArray(Ps,S); CHKERRQ(ierr);
      ierr = ISGetIndices(isB,ind); CHKERRQ(ierr);

      ierr = VecSetValues(RES,PAmBsize,ind[0],S[0],ADD_VALUES);  CHKERRQ(ierr);
      ierr = VecAssemblyBegin(RES);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(RES);CHKERRQ(ierr);

      ierr = ISRestoreIndices(isB,ind); CHKERRQ(ierr);
      ierr = VecRestoreArray(Ps,S); CHKERRQ(ierr);
      ierr = VecDestroy(&Ps);CHKERRQ(ierr);
      ierr = MatDestroy(&PB);CHKERRQ(ierr); 
      ierr = ISDestroy(&isB);CHKERRQ(ierr); 
    }
        
    ierr = VecDestroy(&Pw);CHKERRQ(ierr);   
    ierr = VecDestroy(&Pf);CHKERRQ(ierr);  
    ierr = VecDestroy(&Pr);CHKERRQ(ierr);
    ierr = MatDestroy(&PA);CHKERRQ(ierr); 
    ierr = ISDestroy(&isA);CHKERRQ(ierr);
    
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr); 
    
    //end_time=clock();
    //UpdateTime+=(end_time-start_time);
      
    // ***************** END SOLVE AND UPDATE *******************
  }
  //cout<<"EPS inside Vanka"<<endl;
  //VecView(EPS,PETSC_VIEWER_STDOUT_WORLD);
  // *** Computational info ***
//   cout<<"Grid="<<p<<endl;
//   cout<<"SEARCH   time="<<static_cast<double>(SearchTime)/CLOCKS_PER_SEC<<endl;
//   cout<<"ASSEMBLY time="<<static_cast<double>(AssemblyTime)/CLOCKS_PER_SEC<<endl;
//   cout<<"SOLVE 0  time="<<static_cast<double>(SolveTime0)/CLOCKS_PER_SEC<<endl;
//   cout<<"UPDATE   time="<<static_cast<double>(UpdateTime)/CLOCKS_PER_SEC<<endl;
//   cout<<"A_MAX="<<A_max<<" C_MAX="<<C_max<<endl;
  //cout<<"number of blocks: "<<blocks<<endl;
  return 1;
}

int  lsysPDE::printsol_vtk_ASCII(const  vector <vector <double> > &vt, const unsigned &time, const char type[]) {
  unsigned index=0;
  if(strcmp(type,"linear"))
    index=1; 
  const int eltp[2][6]={{12,10,13,9,5,3},{25,24,26,23,22,21}};
    
  char *filename= new char[60];
  sprintf(filename,"../mesh-generator/output/mesh.%d.%d.%s.vtk",time,p,type);
  std::ofstream fout;
  fout.open(filename);
  if (!fout) {
    cout << "Output mesh file "<<filename<<" cannot be opened.\n";
    exit(0);
  }
  // head ************************************************
  fout<<"# vtk DataFile Version 3.0\nAMR mesh\nASCII\n";
  fout<<"DATASET UNSTRUCTURED_GRID\n";

  // nodes ***********************************************
  //unsigned nvt_index=(index==0)?1:0;
  PetscInt nvt=MeshNodeNum (index);
  fout<<"POINTS "<<nvt<<" double\n";
  for(PetscInt i=0;i<nvt;i++){
    fout<<vt[0][i]<<" "<<vt[1][i]<<" "<<vt[2][i]<<endl;
  }

  // cells ***********************************************
  unsigned nel, counter;
  nel=MeshElementNum(); 
  counter=el->MeshElementNum();
  counter+=el->MeshElementNum("Hex")*NVE[0][index];
  counter+=el->MeshElementNum("Tet")*NVE[1][index];
  counter+=el->MeshElementNum("Wedge")*NVE[2][index];
  counter+=el->MeshElementNum("Quad")*NVE[3][index];
  counter+=el->MeshElementNum("Triangle")*NVE[4][index];
  counter+=el->MeshElementNum("Line")*NVE[5][index];
    
  
  fout<<"\nCELLS "<<nel<<" "<<counter<<endl;
  for(unsigned iel=0;iel<nel;iel++){
    fout<<el->ElementNodeNum(iel,index)<<" ";
    for(unsigned j=0;j<el->ElementNodeNum(iel,index);j++)
      fout<<el->GlobalNodeIndex(iel,j)-1u<<" ";
    fout<<endl;
  }
  fout<<"CELL_TYPES "<<nel<<endl;
  for(unsigned ii=0;ii<nel;ii++){
    short unsigned ielt=el->ElementType(ii);
    fout<<eltp[index][ielt]<<endl;
  }
  
  // data ************************************************
  fout<<"\nPOINT_DATA "<< nvt<<endl;
  for(unsigned i=0;i<Res.size();i++){
    PetscInt size;
    ierr=VecGetSize(Res[i],&size); CHKERRQ(ierr);
    if(size>=nvt){
      PetscScalar *R[0];
      ierr = VecGetArray(Res[i],R); CHKERRQ(ierr);
      fout<<"SCALARS Res"<<SolName[i]<<" double 1"<<endl<<"LOOKUP_TABLE default\n";
      for(PetscInt j=0;j<nvt;j++){
	fout<<R[0][j]<<endl;
      } 
      fout<<endl;
      ierr = VecRestoreArray(Res[i],R); CHKERRQ(ierr);
    }
  }
  for(unsigned i=0;i<Eps.size();i++){
    PetscInt size;
    ierr=VecGetSize(Eps[i],&size); CHKERRQ(ierr);
    PetscScalar *E[0];
    ierr = VecGetArray(Eps[i],E); CHKERRQ(ierr);
    if(size>=nvt){
      fout<<"SCALARS Eps"<<SolName[i]<<" double 1"<<endl<<"LOOKUP_TABLE default\n";
      for(PetscInt j=0;j<nvt;j++){
	fout<<E[0][j]<<endl;
      } 
      fout<<endl;
      ierr = VecRestoreArray(Eps[i],E); CHKERRQ(ierr);
    }
  }
  
  fout.close();
  delete [] filename;
  
  return 1;
}


