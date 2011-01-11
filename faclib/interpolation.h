#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

void PolyBasis(int n, double *c, double x, double logx);
void PolyFit(int n, double *c, int nd, double *x, double *y);
void SVDFit(int np, double *coeff, double *chisq, double tol,
	    int nd, double *x, double *logx, double *y, double *sig,
	    void Basis(int, double *, double, double));
int NLSQFit(int np, double *p, double tol, int *ipvt,
	    double *fvec, double *fjac, int ldfjac, double *wa, int lwa,
	    int n, double *x, double *logx, double *y, double *sig,
	    void func(int, double *, int , double *, double *, 
		      double *, double *, int, void *), 
	    void *extra);
double Simpson(double *y, int ia, int ib);
int NewtonCotes(double r[], const double x[], int ilast,
                int last_only, int id);

#define UVIP3P uvip3p
void uvip3p(int nd, const double *xd, const double *yd,
		 int ni, const double *xi, double *yi);

#endif
