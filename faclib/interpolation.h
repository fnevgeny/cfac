#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

void SVDFit(int np, double *coeff, double *chisq, double tol,
	    int nd, double *x, double *logx, double *y, double *sig,
	    void Basis(int, double *, double, double));
int NLSQFit(int np, double *p, double tol,
	    double *fvec, double *fjac,
	    int n, double *x, double *logx, double *y, double *sigma,
	    void Func(int np, double *p, int n, double *x, double *logx, 
		     double *y, double *dy, int ndy, void *extra), 
	    void *extra);
double Simpson(double *y, int ia, int ib);
int NewtonCotes(double r[], const double x[], int ilast,
                int last_only, int id);

#define UVIP3P uvip3p
void uvip3p(int nd, const double *xd, const double *yd,
		 int ni, const double *xi, double *yi);

#define UVIP3C uvip3c
void uvip3c(int nd, const double xd[], const double yd[],
		 double c1[], double c2[], double c3[]);

#endif
