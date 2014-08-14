#ifndef FRIENDFUNS_HPP
#define FRIENDFUNS_HPP

double evalWedgehPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalWedgehDPhiDx(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalWedgehDPhiDy(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalWedgehDPhiDz(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);

double evalTethPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalTethDPhidx(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalTethDPhidy(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalTethDPhidz(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);

double evalTrihPhi(const int i, const int j, const int k, const double x, const double y, const int n);
double evalTrihDPhidx(const int i, const int j,  const int k, const double x, const double y, const int n);
double evalTrihDPhidy(const int i, const int j,  const int k, const double x, const double y, const int n);

double evalHexhPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalHexhDPhiDx(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalHexhDPhiDy(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);
double evalHexhDPhiDz(const int i, const int j, const int k, const double x, const double y, const double z, const int n, const int type);

double evalQuadhDPhiDy(const int i, const int j, const int k, const double x, const double y, const int n);
double evalQuadhDPhiDx(const int i, const int j, const int k, const double x, const double y, const int n);
double evalQuadhPhi(const int i, const int j, const int k, const double x, const double y, const int n);

double evalLinehPhi(const int i,const int j,const int k, const double x, const int n);
double evalLinehDPhi(const int i, const int j, const int k, const double x, const int n);

double evalPhi(const int i, const double x, const int n);
double evalDphi(const int i, const double x, const int n);

double evalTriPhi(const int i, const int j, const double x, const double y, const int n);
double evalTriDPhidx(const int i, const int j, const double x, const double y, const int n);
double evalTriDPhidy(const int i, const int j, const double x, const double y, const int n);

double evalTetraPhi(const int i, const int j, const int k, const double x, const double y, const double z, const int n);
double evalTetraDPhidx(const int i, const int j, const int k, const double x, const double y, const double z, const int n);
double evalTetraDPhidy(const int i, const int j, const int k, const double x, const double y, const double z, const int n);
double evalTetraDPhidz(const int i, const int j, const int k, const double x, const double y, const double z, const int n);

#endif
