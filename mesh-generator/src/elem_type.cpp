
#include "hGaussPoints.hpp"
#include "elem_type.hpp"
#include <cmath>


const double* hGauss[6][12]={ {hGaussCube0[0], hGaussCube1[0], hGaussCube2[0], hGaussCube3[0], hGaussCube4[0], hGaussCube5[0],
			        hGaussCube6[0], hGaussCube7[0], hGaussCube8[0], hGaussCube9[0], hGaussCube10[0], hGaussCube11[0]},

			      {hGaussTetra0[0], hGaussTetra1[0], hGaussTetra2[0], hGaussTetra3[0], hGaussTetra4[0], hGaussTetra5[0],
			        hGaussTetra6[0], hGaussTetra7[0], hGaussTetra8[0], hGaussTetra9[0], hGaussTetra10[0], hGaussTetra11[0]},

			      {hGaussWedge0[0], hGaussWedge1[0], hGaussWedge2[0], hGaussWedge3[0], hGaussWedge4[0], hGaussWedge5[0],
			        hGaussWedge6[0], hGaussWedge7[0], hGaussWedge8[0], hGaussWedge9[0], hGaussWedge10[0], hGaussWedge11[0]},

			      {hGaussSquare0[0], hGaussSquare1[0], hGaussSquare2[0], hGaussSquare3[0], hGaussSquare4[0], hGaussSquare5[0], 
			        hGaussSquare6[0], hGaussSquare7[0], hGaussSquare8[0], hGaussSquare9[0], hGaussSquare10[0], hGaussSquare11[0]},

			      {hGaussTriangle0[0], hGaussTriangle1[0], hGaussTriangle2[0], hGaussTriangle3[0], hGaussTriangle4[0], hGaussTriangle5[0],
			        hGaussTriangle6[0], hGaussTriangle7[0], hGaussTriangle8[0], hGaussTriangle9[0], hGaussTriangle10[0], hGaussTriangle11[0]},

			      {hGaussLine0[0], hGaussLine1[0], hGaussLine2[0], hGaussLine3[0], hGaussLine4[0], hGaussLine5[0],
			        hGaussLine6[0], hGaussLine7[0], hGaussLine8[0], hGaussLine9[0], hGaussLine10[0], hGaussLine11[0]},			    
};
const short unsigned typeSize[6] = {2,2,6,1,1,1};
elem_type::~elem_type(){

  for(int i=0;i<numOfNodes;i++)
    delete[] I[i];
  delete [] I;

  for(int i=0;i<3;i++){
    for(int j=0; j<numOfGaussPts;j++)
      delete[] dphi[i][j];
    delete[] dphi[i];
  }
  delete[] dphi;
  
  for(int i=0;i<numOfGaussPts;i++)
    delete[] lphi[i];
  delete[] lphi;

  for(int i=0;i<numOfGaussPts;i++)
    delete[] X[i];
  delete[] X;
  
  
  delete[] x;

  delete[] gaussWeights;

  for(int i=0;i<numOfShapeFuns;i++)
    delete[] IND[i];
  delete [] IND;

  for(int i=0;i<typeSize[element];i++){
    for(int j=0;j<numOfGaussPts;j++)
      delete[] phi[i][j];
    delete [] phi[i];
  }
  delete [] phi;

  for(int i=0;i<typeSize[element];i++){
    for(int j=0;j<numOfGaussPts;j++)
      delete[] dphidxi[i][j];
    delete [] dphidxi[i];
  }
  delete[] dphidxi; 

  for(int i=0;i<typeSize[element];i++){
    for(int j=0;j<numOfGaussPts;j++)
      delete[] dphideta[i][j];
    delete [] dphideta[i]; 
  }
  delete [] dphideta; 

  for(int i=0;i<typeSize[element];i++){
    for(int j=0;j<numOfGaussPts;j++)
      delete[] dphidzeta[i][j];
    delete [] dphidzeta[i]; 
  }
  delete [] dphidzeta; 

  for(int i=0;i<typeSize[element];i++)
    delete pt_h[i];
  
  delete pt_lag;
  delete [] Dphi;
};

elem_type::elem_type(const int solid, const int p, const int n, const int orderOfMapping){ //1-based index

  hierarchicP = p;
  if(p<1) error("basis function order p < 1.");

  element=solid;
  basish* pt_[6];
  pt_[0] = new hexh(p,0);
  pt_[1] = new teth(p,0);
  pt_[2] = new wedgeh(p,0);
  pt_[3] = new quadh(p);
  pt_[4] = new trih(p);
  pt_[5] = new lineh(p);

  pt_h[0] = pt_[element];
 
  if(element==0) pt_h[1] = new hexh(p,1);
  if(element==1) pt_h[1] = new teth(p,1);
  if(element==2) {
    for(int i=1;i<6;i++)
    pt_h[i] = new wedgeh(p,i);
  }
    
  for(int i=0;i<6;i++)
    if(i != element) delete pt_[i];

  int gaussOrder=0;
  numOfShapeFuns = pt_h[0]->dimOfBasis();
  dim = pt_h[0]->spatialDim();
  
  switch(element){

  case 0: //hex
    projectionMatrix = &elem_type::ProjectionMatrix3D;
    spatialToRefPtr = &elem_type::mapSpatialToRefElem3D;
    gaussOrder = ceil( (n-1)/2. )+1;
    if(gaussOrder>11)error("Integration not exact! ");
    numOfGaussPts = hGaussPointsCube[gaussOrder];
    break;

  case 1: //teth
    projectionMatrix = &elem_type::ProjectionMatrix3D;
    spatialToRefPtr = &elem_type::mapSpatialToRefElem3D;
    gaussOrder = (n<=5) ? ceil( (n-1)/2. ) : ceil( (n+1)/2. );
    if(gaussOrder>11)error("Integration not exact! ");
    numOfGaussPts = hGaussPointsTetra[gaussOrder];
    break;

  case 2: //wedgeh
    projectionMatrix = &elem_type::ProjectionMatrix3D;
    spatialToRefPtr = &elem_type::mapSpatialToRefElem3D;
    gaussOrder = (n<8) ? ceil( (n-1)/2. ) : ceil( n/2. );
    if(gaussOrder>11)error("Integration not exact! ");
    numOfGaussPts = hGaussPointsWedge[gaussOrder];
    break;

  case 3: //quadh
    projectionMatrix = &elem_type::ProjectionMatrix2D;
    spatialToRefPtr = &elem_type::mapSpatialToRefElem2D;
    gaussOrder = ceil( (n-1)/2. );
    if(gaussOrder>11)error("Integration not exact! ");
    numOfGaussPts = hGaussPointsSquare[gaussOrder]; 
    break;

  case 4: //trih
    projectionMatrix = &elem_type::ProjectionMatrix2D;
    spatialToRefPtr = &elem_type::mapSpatialToRefElem2D;
    gaussOrder = (n<8) ? ceil( (n-1)/2. ) : ceil( n/2. );
    if(gaussOrder>11)error("Integration not exact! ");
    numOfGaussPts = hGaussPointsTriangle[gaussOrder];
    break;

  case 5: //lineh
    projectionMatrix = &elem_type::ProjectionMatrix2D;
    spatialToRefPtr = &elem_type::mapSpatialToRefElem1D;
    gaussOrder = ceil( (n-1)/2. );
    if(gaussOrder>11)error("Integration not exact! ");
    numOfGaussPts = hGaussPointsLine[gaussOrder];
    break;
 
  default: error("Invalid input for solid!");
  }

  IND = new int* [numOfShapeFuns];
  for(int i=0;i<numOfShapeFuns;i++)
    IND[i]=new int[3];

  int* pt = pt_h[0]->getModes();
  int count=0;

  /* Element DOF map: I[type][location][index] listing to I[count]  */ //tested FOR ALL DIMENSIONS!
  // order: 1-vertices, 2-edges, 3-quad faces (if any), 4-tri faces (if any), 5- interior
  int interior=0; 
  int quadFace=0;
  int triFace=0;

  for(int i=0; i < pt[0]; i++){ //vertex modes
    IND[count][0]=0;
    IND[count][1]=i;
    IND[count][2]=0;
    count++;
  }

  if (dim != 3){ // for 1D, 2D elements.
  for(int i=0; i<p-1;i++){ // p: 1-based
    for(int location=0;location<pt[1];location++){
      IND[count][0]=1;
      IND[count][1]=location;
      IND[count][2]=i;
      count++;
    }
    if (element == 5) continue;
    
    if(i>1) //interior modes for tri/quad   
      for(int j=0; j<i-1; j++){
	IND[count][0]=2;
	IND[count][1]=0;
	IND[count][2]=interior;
	count++;interior++;
      }
   
    if((i>0)&&(element==4)){ //additional tri interior mode (if element=triangle)
      IND[count][0]=2;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++; interior++;
    }
  }
  }
  
  else { // 3-D tested!!!
  for(int i=0; i<p-1;i++){ 
    for(int location=0;location<pt[1];location++){//edge modes
      IND[count][0]=1;
      IND[count][1]=location;
      IND[count][2]=i;
      count++;
    }
    
    if(i>1) // quad + tri face modes 
    for(int j=0;j<i-1;j++){
    for(int location=0; location<pt[2]; location++){ 
      IND[count][0]=2;
      IND[count][1]=location;
      IND[count][2]=quadFace;
      count++;
    }
    quadFace++;
    
    for(int location=pt[2]; location<pt[2]+pt[3]; location++){ 
      IND[count][0]=2;
      IND[count][1]=location;
      IND[count][2]=triFace;
      count++;
    }
    triFace++;

    }
    
  if(i>0){//additional tri face mode
      for(int location=pt[2]; location<pt[2]+pt[3]; location++){
	IND[count][0]=2;
	IND[count][1]=location;
	IND[count][2]=triFace;
	count++;
      }triFace++;
  }
    
    if((i>1)&&(element==1)) // internal modes for tet
    for(int j=0;j<(i-1)*i/2;j++){
      IND[count][0]=3;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++;interior++;
    }
    
    else if((i>2)&&(element==2)) // internal modes for wedge
    for(int j=0;j<(i-1)*(i-2)/2;j++){
      IND[count][0]=3;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++;interior++;
    }
    
    else if((i>3)&&(element==0)) // internal modes for hex
    for(int j=0;j<(i-2)*(i-3)/2;j++){
      IND[count][0]=3;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++;interior++;
    }
    
  }
  }
 
  /* Gauss points and Gauss weights */
  X = new double* [numOfGaussPts];
  for(int i=0;i<numOfGaussPts;i++)
    X[i]=new double[3];

  x = new double[3];
  
  gaussWeights = new double [numOfGaussPts];

  for(int i=0;i<numOfGaussPts;i++){
    gaussWeights[i]=hGauss[solid][gaussOrder][i];
    X[i][0] = hGauss[solid][gaussOrder][i+numOfGaussPts];
    X[i][1] = (solid==5) ? 0 : hGauss[solid][gaussOrder][i+numOfGaussPts*2];
    X[i][2] = (solid<=2) ? hGauss[solid][gaussOrder][i+numOfGaussPts*3]:0;
  }

  phi=new double**[typeSize[element]];
  dphidxi=new double**[typeSize[element]];
  dphideta=new double**[typeSize[element]];
  dphidzeta=new double**[typeSize[element]];

  for(int i=0;i<typeSize[element];i++){
    phi[i]=new double*[numOfGaussPts];
    dphidxi[i]=new double*[numOfGaussPts];
    dphideta[i]=new double*[numOfGaussPts];
    dphidzeta[i]=new double*[numOfGaussPts];

    for(int j=0;j<numOfGaussPts;j++){
      phi[i][j]=new double[numOfShapeFuns];
      dphidxi[i][j]=new double[numOfShapeFuns];
      dphideta[i][j]=new double[numOfShapeFuns];
      dphidzeta[i][j]=new double[numOfShapeFuns];
    }
  }

  for(int k=0;k<typeSize[element];k++)
    for(int i=0;i<numOfGaussPts;i++)
      for (int j=0;j<numOfShapeFuns;j++){
	phi[k][i][j] = pt_h[k]->eval_phi(IND[j],X[i],p);
	dphidxi[k][i][j] = pt_h[k]->eval_dphidx(IND[j],X[i],p);
	dphideta[k][i][j] = pt_h[k]->eval_dphidy(IND[j],X[i],p);
	dphidzeta[k][i][j] = pt_h[k]->eval_dphidz(IND[j],X[i],p);
      }
     
  Dphi = new double**[3];
  
  Dphi[0] = dphidxi[0];
  Dphi[1] = dphideta[0];
  Dphi[2] = dphidzeta[0];
  
  /* Lagrage basis: Only affine and quadratic mappings are allowed */
  int p1 = orderOfMapping;
  if(p1>2)error("Order of mapping should be linear or quadratic only!");
  basis* ptLag[6];
  ptLag[0] = new hex();
  ptLag[1] = new tet();
  ptLag[2] = new wedge();
  ptLag[3] = new quad();
  ptLag[4] = new tri();
  ptLag[5] = new line();
  
  pt_lag = ptLag[element];
  
  for(int i=0;i<6;i++)
    if(i!=element)
      delete ptLag[i];

  numOfNodes=pt_lag->ElementDOF(p1);
  
  /* Create indices for Lagrange basis upto order 2 */
  I = new int* [numOfNodes];
  for(int i=0;i<numOfNodes;i++)
    I[i]=new int[3];

  count = 0;

  /* Initialize indices to zero */
 
  for(int k=0;k<numOfNodes;k++){
    I[count][0]=0;
    I[count][1]=0;
    I[count][2]=0;
    count++;
  }

  count=0;
  switch(element){
  case 0://hex
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	I[count][0]=p1*((2+j+i)%2);
	I[count][1]=i*p1;
	count++;
      }
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	I[count][0]=p1*((2+j+i)%2);
	I[count][1]=i*p1;
	I[count][2]=1;
	count++;
     }
     break;
  case 1:break;//tet
  case 2:break;//wedge

  case 3: // tested!
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	I[count][0]=p1*((2+j+i)%2);
	I[count][1]=i*p1;
	count++;
      }
    if(p1==1) break;
    else 
      for(int i=0;i<2;i++)
	for(int j=0;j<2;j++){
	  I[count][0]=(4+j-i)%(3-i);
	  I[count][1]=(2+i+j)%2+i;
	  count++;
	}
    I[count][0]=1;I[count][1]=1; break;

  case 4: //tested!
    for(int i=0;i<2;i++)
      for(int j=0;j<2-i;j++){
	I[count][0]=j*p1;
	I[count][1]=i*p1;
	count++;
      }
    if(p1==1) break;
    else
      for(int i=0;i<2;i++)
	for(int j=0;j<=i;j++){
	  I[count][0]=(j+1)%2;
	  I[count][1]=i;
	  count++;
	}
    break;

  case 5: //tested!
    for(int i=0; i<=p1; i++)
      I[(p1-i+1)%(p1+1)][0]=i; break; 
  }

  /* Lagrange basis functions */
  lphi = new double*[numOfGaussPts];
  for(int i=0;i<numOfGaussPts;i++)
    lphi[i]= new double[numOfNodes];

  /* Partial derivatives of Lagrange basis functions */
  dphi = new double** [3];
  for(int i=0;i<3;i++){
    dphi[i]=new double*[numOfGaussPts];
    for(int j=0;j<numOfGaussPts;j++)
      dphi[i][j]=new double[numOfNodes];
  }

  for(int i=0;i<numOfGaussPts;i++)
    for (int j=0;j<numOfNodes;j++){
      lphi[i][j] = pt_lag->eval_phi(I[j],X[i],p1);
      dphi[0][i][j] = pt_lag->eval_dphidx(I[j],X[i],p1);
      dphi[1][i][j] = pt_lag->eval_dphidy(I[j],X[i],p1);
      dphi[2][i][j] = pt_lag->eval_dphidz(I[j],X[i],p1);
    } 
  
}

void elem_type::ProjectionMatrix2D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const{
  
  int numOfVertices = pt_lag->ElementDOF(1);
  
  double** x_; // refined Lagrange nodes of the element
 
  /* Build the local projection matrix */
  int lagrangeNodes = pt_lag->ElementDOF(p_);
  
  x_ = new double* [lagrangeNodes];
  for(int i=0;i<lagrangeNodes;i++)
    x_[i] = new double[3];
  pt_lag->RefinedNodes(p_,x_);
  
  /* Build the global projection matrix */  
  for(int i=0;i<lagrangeNodes;i++)
    for(int j=0;j<numOfShapeFuns;j++)
      P[i][j]=sign[j]*pt_h[0]->eval_phi(IND[j],x_[i],hierarchicP);
     
  /* Global nodes of the refined mesh */ // Only affine map is considered. But do for a general one.
  for(int i=0;i<lagrangeNodes;i++)
    for(int coordinate=0;coordinate<3;coordinate++){
      node[i][coordinate]=0.;
      for(int k=0;k<numOfVertices;k++) 
	node[i][coordinate] += vertex[k][coordinate]*pt_h[0]->eval_phi(IND[k],x_[i],hierarchicP);//TODO check if this works for hierarchicP=1
    }
    
  for(int i=0;i<lagrangeNodes;i++)
      delete[] x_[i];
    
  delete[] x_;
}

/* 3D visualization is done by projecting to Lagrange hexahedral mesh */
void elem_type::ProjectionMatrix3D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const{
  
  hex hex1;
  int numOfVertices = hex1.ElementDOF(1); // change to hex pointer
  double xi,eta,zeta;
 
  double** x_; // refined Lagrange nodes of hex
  double** X_; // refined Lagrange nodes mapped to hierarchic tet
  vector < vector < double > > v; // vertices of the spatial hexahedral element.
  
  /* Build the local projection matrix */
  int lagrangeNodes = hex1.ElementDOF(p_);
  
  x_ = new double* [lagrangeNodes]; X_ = new double* [lagrangeNodes]; 
  for(int i=0;i<lagrangeNodes;i++){
    x_[i] = new double[3];
    X_[i] = new double[3];
  }
  hex1.RefinedNodes(p_,x_);
  
  for(int i=0;i<lagrangeNodes;i++)
    for(int j=0;j<3;j++)
      X_[i][j]=x_[i][j];
  
  /* Map hexahedral Lagrange nodes to hierarchic reference tetrahedron */
  if(element==1) 
  for(int i=0;i<lagrangeNodes;i++){
    
    // hex to unit tet
    xi = (1+x_[i][0])*(1-x_[i][1])*(1-x_[i][2])/8.;
    eta = (1+x_[i][1])*(1-x_[i][2])/4.;
    zeta = (1+x_[i][2])/2.;
    
    // unit tet to hierarchic tet
    X_[i][0] = 2*xi+eta+zeta-1;
    X_[i][1] = sqrt(3)*eta+zeta/sqrt(3);
    X_[i][2] = sqrt(8./3)*zeta;
  
  }

  if(element==2) 
  for(int i=0;i<lagrangeNodes;i++){
    
    // hex to unit wedge
    xi = (1+x_[i][0])*(1-x_[i][1])/4.;
    eta = (1+x_[i][1])/2.;
    zeta = x_[i][2];
    
    // unit wedge to hierarchic wedge 
    X_[i][0] = 2*xi+eta-1;
    X_[i][1] = sqrt(3)*eta;
    X_[i][2] = zeta;
  
  }
  
  /* Build the global projection matrix */
  for(int i=0;i<lagrangeNodes;i++)
    for(int j=0;j<numOfShapeFuns;j++)
      P[i][j]=sign[j]*pt_h[type[j]]->eval_phi(IND[j],X_[i],hierarchicP);
  
  /* Global nodes of the refined mesh */ // Only affine map is considered. But do for a general one.
  v.resize(numOfVertices);
  for(int i=0;i<numOfVertices;i++)
    v[i].resize(3);

  for(int i=0;i<vertex.size();i++)
    for(int j=0;j<3;j++)
      v[i][j]=vertex[i][j];
    
  if(element==1)
  for(int j=0;j<3;j++){
    v[0][j]=vertex[0][j];
    v[1][j]=vertex[1][j];
    v[2][j]=vertex[2][j];
    v[3][j]=vertex[2][j];
    v[4][j]=vertex[3][j];
    v[5][j]=vertex[3][j];
    v[6][j]=vertex[3][j];
    v[7][j]=vertex[3][j];
  }
  
  if(element==2) 
  for(int j=0;j<3;j++){
    v[0][j]=vertex[0][j];
    v[1][j]=vertex[1][j];
    v[2][j]=vertex[2][j];
    v[3][j]=vertex[2][j];
    v[4][j]=vertex[3][j];
    v[5][j]=vertex[4][j];
    v[6][j]=vertex[5][j];
    v[7][j]=vertex[5][j];
  }
 
 int ind[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
		  
  for(int i=0;i<lagrangeNodes;i++)
    for(int coordinate=0;coordinate<3;coordinate++){
      node[i][coordinate]=0.;
      for(int k=0;k<numOfVertices;k++) 
	node[i][coordinate] += v[k][coordinate]*hex1.eval_phi(ind[k],x_[i],1);
    }
    
  for(int i=0;i<lagrangeNodes;i++){
      delete[] x_[i];
      delete[] X_[i];
  }
  delete[] X_;
  delete[] x_;
}

void elem_type::mapSpatialToRefElem1D(vector<vector<double> >node, const unsigned& gaussIndex, const int* sign, double& weight, double* phi_, double** gradPhi, double* X, const int* type) const{
  double det;
  double J = 0;
     
  for(int k=0;k<numOfNodes;k++)
    J += node[k][0]*dphi[0][gaussIndex][k];
 
  det = J;
  if(det==0)
    throw runtime_error("No unique solution! Zero determinant at elem_type::mapSpatialToRefElem1D");
  weight = fabs(det)*gaussWeights[gaussIndex];

  for(int i=0;i<numOfShapeFuns;i++){
    phi_[i] = phi[0][gaussIndex][i];
    gradPhi[i][0] = dphidxi[0][gaussIndex][i]/det;
    gradPhi[i][1] = 0.;
    gradPhi[i][2] = 0.;
  }
}

void elem_type::mapSpatialToRefElem2D(vector<vector<double> > node, const unsigned& gaussIndex, const int* sign, double& weight, double* phi_, double** gradPhi, double* x, const int* type) const{
 
  double det;
  double J[2][2]={{0.,0.},{0.,0.}};

  /* Determinant for arbitrary order mapping. */
  for(int j=0;j<2;j++)
    for(int i=0;i<2;i++){
      J[i][j]=0.; 
      for(int k=0;k<numOfNodes;k++)
	J[i][j] += node[k][i]*dphi[j][gaussIndex][k]; // node[index][x/y/z]
    }

  /* Lagrange basis reference triangle to hierarchic basis reference triangle mapping */
  if(element==4) {
    J[0][1] = -J[0][0]/(2*sqrt(3)) + J[0][1]/sqrt(3);
    J[0][0] *= 0.5;
    J[1][1] = -J[1][0]/(2*sqrt(3)) + J[1][1]/sqrt(3);
    J[1][0] *= 0.5;
  }

  det = J[0][0]*J[1][1] - J[0][1]*J[1][0]; 
  if(det==0)
    throw runtime_error("No unique solution! Zero determinant at elem_type::mapSpatialToRefElem2D");
  
  weight = fabs(det)*gaussWeights[gaussIndex];

  for(int i=0;i<numOfShapeFuns;i++){
    phi_[i] = sign[i]*phi[0][gaussIndex][i];
    gradPhi[i][0] = sign[i]*(dphidxi[0][gaussIndex][i]*J[1][1] - dphideta[0][gaussIndex][i]*J[1][0])/det;
    gradPhi[i][1] = sign[i]*(-dphidxi[0][gaussIndex][i]*J[0][1] + dphideta[0][gaussIndex][i]*J[0][0])/det;
    gradPhi[i][2] = 0.;
  }
  
   if(element==4){ // correct for linear map only. Do for general case.
    x[0]=0.;x[1]=0;x[2]=0.;
    for(int coordinate=0;coordinate<3;coordinate++)
      for(int k=0;k<3;k++)
	x[coordinate] += node[k][coordinate]*Phi()[gaussIndex][k];
   }
    else{
    x[0]=0.;x[1]=0;x[2]=0.;
    for(int coordinate=0;coordinate<3;coordinate++)
      for(int k=0;k<numOfNodes;k++)
	x[coordinate] += node[k][coordinate]*lphi[gaussIndex][k];
      
    }
}

void elem_type::mapSpatialToRefElem3D(vector<vector<double> > node, const unsigned& gaussIndex, const int* sign, double& weight, double* phi_, double** gradPhi, double* x, const int* type) const{
  
  double det;
  double J[3][3];
  double Jinv[3][3];
  
  for(int j=0;j<3;j++){
    for(int i=0;i<3;i++){
      J[i][j]=0.; 
      for(int k=0;k<numOfNodes;k++)
	J[i][j] += node[k][i]*Dphi[j][gaussIndex][k];
    }
  }

  det = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0]) +  J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
  
  if(det==0)
    throw runtime_error("No unique solution! Zero determinant at elem_type::mapSpatialToRefElem3D");
  
  Jinv[0][0]= (-J[1][2]*J[2][1] + J[1][1]*J[2][2])/det; 
  Jinv[0][1]= ( J[0][2]*J[2][1] - J[0][1]*J[2][2])/det;
  Jinv[0][2]= (-J[0][2]*J[1][1] + J[0][1]*J[1][2])/det;
  Jinv[1][0]= ( J[1][2]*J[2][0] - J[1][0]*J[2][2])/det;
  Jinv[1][1]= (-J[0][2]*J[2][0] + J[0][0]*J[2][2])/det; 
  Jinv[1][2]= ( J[0][2]*J[1][0] - J[0][0]*J[1][2])/det;
  Jinv[2][0]= (-J[1][1]*J[2][0] + J[1][0]*J[2][1])/det; 
  Jinv[2][1]= ( J[0][1]*J[2][0] - J[0][0]*J[2][1])/det;
  Jinv[2][2]= (-J[0][1]*J[1][0] + J[0][0]*J[1][1])/det;
  
  weight = fabs(det)*gaussWeights[gaussIndex];

  for(int i=0;i<numOfShapeFuns;i++){
    phi_[i] = sign[i]*phi[type[i]][gaussIndex][i];
    
    gradPhi[i][0] = sign[i]*(Jinv[0][0]*dphidxi[type[i]][gaussIndex][i]+Jinv[1][0]*dphideta[type[i]][gaussIndex][i]+Jinv[2][0]*dphidzeta[type[i]][gaussIndex][i]);
    gradPhi[i][1] = sign[i]*(Jinv[0][1]*dphidxi[type[i]][gaussIndex][i]+Jinv[1][1]*dphideta[type[i]][gaussIndex][i]+Jinv[2][1]*dphidzeta[type[i]][gaussIndex][i]);
    gradPhi[i][2] = sign[i]*(Jinv[0][2]*dphidxi[type[i]][gaussIndex][i]+Jinv[1][2]*dphideta[type[i]][gaussIndex][i]+Jinv[2][2]*dphidzeta[type[i]][gaussIndex][i]);
  }

  // gauss point in spatial domain
    x[0]=0.;x[1]=0;x[2]=0.;
    for(int coordinate=0;coordinate<3;coordinate++)
      for(int k=0;k<numOfNodes;k++)
	x[coordinate] += node[k][coordinate]*Phi()[gaussIndex][k];
}

void elem_type::error(const string& msg) const{
  throw std::runtime_error(msg);
}


