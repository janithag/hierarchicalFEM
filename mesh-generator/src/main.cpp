#include "main.hpp"
#include "elem_type.hpp"
#include "multigrid.hpp"
#include <sstream>

// Program usage in parallel:  mpiexec -np 4 mesh.out [-help] [all PETSc options] 
static char help[] = "Multigrid solver with PETSC\n\n";

int main(int argc, char **argv){

  PetscErrorCode ierr;
  PetscMPIInt    size;
  PetscInitialize(&argc, &argv, (char *)0, help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
  try{
    unsigned short numOfMGLevels, nr,p;
    numOfMGLevels = 1;  
    p=POrder; //1-based
    nr=0; // number of maximum refinement levels
    int tmp=numOfMGLevels;  numOfMGLevels += nr;  nr = tmp;

    if( DIM[0]*DIM[1] || DIM[1]*DIM[2] || DIM[0]*DIM[2])
      throw runtime_error("Invalid spatial dimension");

    char* infile = new char [50];
  
    std::stringstream s;
    s << "../mesh-generator/input/mesh.mesh_" << DIMENSION + 1  <<  "D_wedge2";//;_square10by10;_tri10by10_hybrid_finalHybrid
    //										_tetCube_hexCube_wedge128_hybridY,_hybridZ,hex16by16
    char* a = new char[s.str().size()+1];
    
    a[s.str().size()]=0;
    memcpy(a,s.str().c_str(),s.str().size());

    sprintf(infile,a);
    delete[] a;
//     ofstream fout;
//     fout.open("L2 error.txt");
//     fout.precision(14);
//     fout<<"p \t DOF \t L2 error";
    MultiGrid* mg;
    clock_t start_time, end_time;
    
    for(int i=p-1;i<p;i++){
      start_time = clock();
      mg = new MultiGrid(infile, i+1);
      delete mg;
      end_time = clock();
      //cout<<"Time elapsed: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;
    }
    
    //fout.close();
    ierr = PetscFinalize(); CHKERRQ(ierr);
  
    delete[] infile;
  }

  catch(exception& ex){
    cerr<<"Detected exception: "<<ex.what()<<endl;
  }

  return(0);
}

unsigned localP_(const double& x, const double& y, const double& z, unsigned& P){// centroid of the element
  int localp=P;
  if(x<0.5) localp=P;
    return localp;
}

bool BoundaryND(const double& x, const double& y, const double& z, const char name[], double& value, const int facet){
  bool test=0; //test=1 :essential BC, test=0 : natural BC 
  value=0.;

  if(!strcmp(name,"U")){
    if(facet==3){  //z
     test=1.;
     value=0.;
    }
    
    else if (facet==1){ //z
      test=1.;
     value=0.;
    }
      
    else if(facet==4){//x
      test=1.;
      value=0.;
    }
      
    else if(facet==5){//y
      test=1.;
      value=0.;	
    }
    
    else if(facet==2){//x
      test=1.;
      value=0.;
    }  
    
    else if(facet==6){//y
      test=1.;
      value=0.;
    }  
  }
  
  else if(!strcmp(name,"V")){
    if(facet==1){
      test=1;
      value=0.;
    }  
    else if(facet==2){  
      test=1;
     value=0.;
    }
    else if(facet==3){  
      test=1;
      value=0.;	
    }
    else if(facet==4){  
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"W")){
    if(1==facet){
      test=1;
      value=0.;
    }  
    else if(2==facet ){  
      test=1;
      value=0.;
    }
    else if(3==facet ){  
      test=1;
      value=0.;
    }
    else if(4==facet ){  
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(1==facet){
      test=0;
      value=0.;
    }  
    else if(2==facet ){  
      test=0;
      value=0.;
    }
    else if(3==facet ){  
      test=0;
      value=0.;
    }
    else if(4==facet ){  
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"T")){
    if(1==facet){
      test=0;
      value=0.;
    }  
    else if(2==facet ){  
      test=1;
      value=0.;
    }
    else if(3==facet ){  
      test=1;
      value=0.;
    }
    else if(4==facet ){  
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"L")){
    if(1==facet){
      test=0;
      value=0.;
    }  
    else if(2==facet ){  
      test=1;
      value=0.;
    }
    else if(3==facet ){  
      test=1;
      value=0.;
    }
    else if(4==facet ){  
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"T0")){
    if(1==facet){
      test=0;
      value=0.;
    }  
    else if(2==facet ){  
      test=1;
      value=1.;
      if(y>1.9999) value=0.;
    }
    else if(3==facet ){  
      test=0;
      value=0.;
    }
    else if(4==facet ){  
      test=1;
      value=0.;
    }
  }
  return test;
}

bool Boundary1D(const double &x, const double &y, const double &z,const char name[], double &value, const int FaceName){
  bool test=1; //Dirichlet
  value=0.;
  if(!strcmp(name,"T")){
    if(1==FaceName){
      test=1;
      value=1.;
    }  
    else if(2==FaceName ){  
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"L")){
    if(1==FaceName){
      test=0;
      value=0.;
    }  
    else if(2==FaceName ){  
      test=1;
      value=3.;
    }
  }
  else if(!strcmp(name,"T0")){
    if(1==FaceName){
      test=0;
      value=0.;
    }  
    else if(2==FaceName ){  
      test=1;
      value=-3.;
    }
  }
  return test;
}


double Init(const double &x, const double &y, const double &z,const char name[]){
  double value=0.;
  if(!strcmp(name,"U")){
    double r=sqrt(y*y+z*z);
    value=0*(r-0.1)*(r+0.1);
  }
  else if(!strcmp(name,"V")){
    value=0;
  } 
  else if(!strcmp(name,"W")){
    value=0;
  }
  else if(!strcmp(name,"P")){
    value=0;
  }
  return value;
}




