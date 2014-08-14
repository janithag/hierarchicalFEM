#ifndef HFUN_HPP
#define HFUN_HPP

double legendre(const int n, const double x);
double dlegendre(const int m, const double x);
double lobatto(const int i, const double x);
double dlobatto(const int i, const double x);
double L(const int i, const double x, const double y);
double kernel(const int i, const double x);
double dLdx(const int i, const double x, const double y);
double dLdy(const int i, const double x, const double y);
double dkernel(const int i, const double x);
double tetL(const int i, const double x, const double y,  const double z);
double dtetLdx(const int i, const double x, const double y,  const double z);
double dtetLdy(const int i, const double x, const double y,  const double z);
double dtetLdz(const int i, const double x, const double y,  const double z);

#endif

// This function evaluates the orthogonal basis functions using Gram-Schmidt process as defined in fig 3 of the paper Adjerid, Aiffa, Flaherty.

// This function evaluates the bilinear form associated with the Laplacian as defined in 3.1 of the paper Adjerid, Aiffa, Flaherty.
