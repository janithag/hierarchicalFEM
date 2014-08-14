#include <iostream>
#include <cmath>
#include <stdexcept>

using std::runtime_error;
using std::cout;
// This is a recursive function to generate Legendre polynomial of degree n.
double legendre(const int n, const double x) {
  if (n<0) throw runtime_error("negative order for legendre polynomial!");
  if(n==0)
    return 1;
  else if(n==1) 
    return x;
  else
    return (2*n-1)*x*legendre(n-1,x)/n-(n-1)*legendre(n-2,x)/n;      
}

// Computes the derivative of the legendre polynomial.
double dlegendre(const int m, const double x){
  if (m<0) throw runtime_error("negative order for the derivative of legendre polynomial!");
  if(m==0)
    return 0;
  else if(m==1)
    return 1;
  else
    return (2*m-1)*legendre(m-1,x)+dlegendre(m-2,x);
}

// Computes the 2nd derivative of the legendre polynomial.
double d2legendre(const int m, const double x){
  if (m<0) throw runtime_error("negative order for the second derivative of legendre polynomial!");
  if((m==0)||(m==1))
    return 0;
  else
    return (2*m-1)*dlegendre(m-1,x)+d2legendre(m-2,x);
}

// This function evaluates the lobatto function declared in 3.12c in Szabo/Babuska.
double lobatto(const int i, const double x) {
    return (legendre(i,x)-legendre(i-2,x))/sqrt(2*(2*i-1));
}

// This function evaluates the derivative of the lobatto function declared in 3.12c in Szabo/Babuska.
double dlobatto(const int i, const double x) {
    return (dlegendre(i,x)-dlegendre(i-2,x))/sqrt(2*(2*i-1));  
}

// This function evaluates the barycentric functions L1, L2, L3 for triangle. (defined 6.6a-c in Szabo/Babuska)
double L(const int i, const double x, const double y){
  switch(i){
    case 0:return (1-x-y/sqrt(3))/2.;break;
    case 1:return (1+x-y/sqrt(3))/2.;break;
    case 2:return y/sqrt(3);break;
    default:throw runtime_error("Invalid index for triangle barycentric functions.");
  }
}

// This function evaluates the partial derivatives of the barycentric functions L1, L2, L3 w.r.t x. (defined 6.6a-c in the textbook)
double dLdx(const int i, const double x, const double y){
  switch(i){
    case 0:return -1/2.;break;
    case 1:return 1/2.;break;
    case 2:return 0;break;
    default:throw runtime_error("Invalid index for triangle barycentric functions.");
  }
}

// This function evaluates the partial derivatives of the barycentric functions L1, L2, L3 w.r.t y. (defined 6.6a-c in the textbook)
double dLdy(const int i, const double x, const double y){
  switch(i){
    case 0:return -1/(sqrt(3)*2);break;
    case 1:return -1/(sqrt(3)*2);break;
    case 2:return 1/sqrt(3);break;
    default:throw runtime_error("Invalid index for triangle barycentric functions.");
  }
}

// kernal function of order i-1 (This function implements the function defined by 2.3c in Adjerid, Aiffa, Flaherty). 
double kernel(const int i, const double x){
  return -sqrt(16*i+8)*dlegendre(i,x)/(i*(i+1)); 
}

// This function evaluates the derivative of the function defined by 2.3c in Adjerid, Aiffa, Flaherty.
double dkernel(const int i, const double x){
  return -sqrt(16*i+8)*d2legendre(i,x)/(i*(i+1));
}

// This function evaluates the tetrahedral coordinates as defined in 13.2a-d in the textbook.
double tetL(const int i, const double x, const double y,  const double z){
  switch(i){
    case 0:return (1-x-y/sqrt(3)-z/sqrt(6))/2.;break;
    case 1:return (1+x-y/sqrt(3)-z/sqrt(6))/2.;break;
    case 2:return (y-z/sqrt(8))/sqrt(3);break;
    case 3:return z*sqrt(3/8.);break;
    default:throw runtime_error("Invalid index for tetrahedron barycentric functions.");
  }
}

// This function evaluates the derivatives of the tetrahedral coordinates as defined in 13.2a-d in the textbook.
double dtetLdx(const int i, const double x, const double y,  const double z){
  switch(i){
    case 0:return -1/2.;break;
    case 1:return 1/2.;break;
    case 2:return 0.;break;
    case 3:return 0.;break;
    default:throw runtime_error("Invalid index for tetrahedron barycentric functions.");
  }
}

// This function evaluates the derivatives of the tetrahedral coordinates as defined in 13.2a-d in the textbook.
double dtetLdy(const int i, const double x, const double y,  const double z){
  switch(i){
    case 0:return -1/(sqrt(3)*2.);break;
    case 1:return -1/(sqrt(3)*2.);break;
    case 2:return 1/sqrt(3);break;
    case 3:return 0.;break;
    default:throw runtime_error("Invalid index for tetrahedron barycentric functions.");
  }
}

// This function evaluates the derivatives of the tetrahedral coordinates as defined in 13.2a-d in the textbook.
double dtetLdz(const int i, const double x, const double y,  const double z){
  switch(i){
    case 0:return -1/(sqrt(6)*2.);break;
    case 1:return -1/(sqrt(6)*2.);break;
    case 2:return -1/(sqrt(8)*sqrt(3));break;
    case 3:return sqrt(3/8.);break;
    default:throw runtime_error("Invalid index for tetrahedron barycentric functions.");
  }
}

// This function evaluates the orthogonal basis functions using Gram-Schmidt process as defined in fig 3 of the paper Adjerid, Aiffa, Flaherty.

// This function evaluates the bilinear form associated with the Laplacian as defined in 3.1 of the paper Adjerid, Aiffa, Flaherty.
