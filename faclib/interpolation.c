#include <string.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlin.h>

#include "interpolation.h"

/* closed Newton-Cotes formulae coeff. */
static double _CNC[5][5] = {
  {0, 0, 0, 0, 0},
  {0.5, 0.5, 0, 0, 0},
  {1.0/3.0, 4.0/3.0, 1.0/3.0, 0, 0},
  {3.0/8, 9.0/8, 9.0/8, 3.0/8, 0},
  {14.0/45, 64.0/45, 24.0/45, 64.0/45, 14.0/45,}
};


void SVDFit(int np, double *coeff, double *chisq, double tol,
	    int nd, double *x, double *logx, double *y, double *sig,
	    void Basis(int, double *, double, double)) {
  int i, j;
  double *afunc;
  double thresh, wmax;
  int infor;
  gsl_matrix *Am, *Vm;
  gsl_vector *bv, *Sv, *wv;
  gsl_vector_view xv;
 
  Am = gsl_matrix_alloc(nd, np);
  Vm = gsl_matrix_alloc(np, np);
  bv = gsl_vector_alloc(nd);
  Sv = gsl_vector_alloc(np);
  
  wv = gsl_vector_alloc(np);

  xv = gsl_vector_view_array(coeff, np);

  afunc = malloc(sizeof(double)*np);

  for (i = 0; i < nd; i++) {
    double weight;
    
    if (logx) {
      Basis(np, afunc, x[i], logx[i]);
    } else {
      Basis(np, afunc, x[i], 0.0);
    }
    
    if (sig) {
      weight = 1.0/sig[i];
      printf("%g %g\n", sig[i], weight);
    } else {
      weight = 1.0;
    }
    
    for (j = 0; j < np; j++) {
      gsl_matrix_set(Am, i, j, afunc[j]*weight);
    }
    gsl_vector_set(bv, i, y[i]*weight);
  }
  
  free(afunc);

  infor = gsl_linalg_SV_decomp(Am, Vm, Sv, wv);
  
  gsl_vector_free(wv);
    
  if (infor != 0) {
    fprintf(stderr, "gsl_linalg_SV_decomp() failed with %d\n", infor);
    abort();
  }
    
  wmax = gsl_vector_get(Sv, 0);
  thresh = tol*wmax;
  for (j = 0; j < np; j++) {
    if (gsl_vector_get(Sv, j) < thresh) gsl_vector_set(Sv, j, 0.0);
  }

  gsl_linalg_SV_solve(Am, Vm, Sv, bv, &xv.vector);
  
  gsl_vector_free(bv);
  gsl_vector_free(Sv);
  gsl_matrix_free(Am);
  gsl_matrix_free(Vm);
}

typedef struct {
    double *x;
    double *logx;
    double *y;
    double *sigma;
    double *fvec;
    double *fjac;
    void *extra;
    void (*function)(int np, double *p, int n, double *x, double *logx, 
		     double *y, double *dy, int ndy, void *extra);
} nlfit_data;

static int nlsqfit_f(const gsl_vector *p, void *data, gsl_vector *f)
{
    nlfit_data *nldata = (nlfit_data *) data;
    
    unsigned int i, n, np, ndy;
    
    double *fvec = nldata->fvec;
    double *x    = p->data;
    
    np = p->size;
    n  = f->size;

    ndy = 0;
    nldata->function(np, x, n, nldata->x,
			    nldata->logx, fvec, 
			    NULL, ndy, nldata->extra);

    for (i = 0; i < n; i++) {
        fvec[i] = fvec[i] - nldata->y[i];
        if (nldata->sigma) {
            fvec[i] /= nldata->sigma[i];
        }
        gsl_vector_set(f, i, fvec[i]);
    }

    return GSL_SUCCESS;
}

static int nlsqfit_df(const gsl_vector *p, void *data, gsl_matrix *J)
{
    nlfit_data *nldata = (nlfit_data *) data;

    unsigned int i, j, k, n, np, ndy;
    
    double *x    = p->data;
    double *fjac = nldata->fjac;
    
    np = p->size;
    n  = J->size1;

    ndy = n;
    nldata->function(np, x, n, nldata->x,
			    nldata->logx, NULL, 
			    fjac, ndy, nldata->extra);

    /* Jacobian matrix J(i,j) = dfi / dxj */
    for (i = 0; i < n; i++) {
	k = i;
	for (j = 0; j < np; j++) {
            if (nldata->sigma) {
                fjac[k] /= nldata->sigma[i];
            }
            gsl_matrix_set(J, i, j, fjac[k]); 
	    k += ndy;
	}
    }

    return GSL_SUCCESS;
}

static int nlsqfit_fdf(const gsl_vector *x,
    void *data, gsl_vector *f, gsl_matrix *J)
{
    nlsqfit_f(x, data, f);
    nlsqfit_df(x, data, J);

    return GSL_SUCCESS;
}

int NLSQFit(int np, double *p, double tol,
	    double *fvec, double *fjac,
	    int n, double *x, double *logx, double *y, double *sigma,
	    void Func(int np, double *p, int n, double *x, double *logx, 
		     double *y, double *dy, int ndy, void *extra), 
	    void *extra) {
  
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status;
    unsigned int iter = 0;

    gsl_multifit_function_fdf f;
    gsl_vector_view xv = gsl_vector_view_array(p, np);

    unsigned int i, maxfev = 5000*np;
    
    nlfit_data nldata;
    nldata.x = x;
    nldata.logx = logx;
    nldata.y = y;
    nldata.sigma = sigma;
    nldata.extra = extra;
    nldata.function = Func;
    nldata.fvec = fvec;
    nldata.fjac = fjac;
    
    f.f = &nlsqfit_f;
    f.df = &nlsqfit_df;
    f.fdf = &nlsqfit_fdf;
    f.n = n;
    f.p = np;
    f.params = &nldata;

    T = gsl_multifit_fdfsolver_lmder;
    s = gsl_multifit_fdfsolver_alloc(T, n, np);
    gsl_multifit_fdfsolver_set(s, &f, &xv.vector);

    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);
        
        if (status) {
            break;
        }

        status = gsl_multifit_test_delta(s->dx, s->x, 0.0, tol);
    } while (status == GSL_CONTINUE && iter < maxfev);

    for (i = 0; i < np; i++) {
        p[i] = gsl_vector_get(s->x, i);
    }
    
    gsl_multifit_fdfsolver_free(s);

    return status;
}

double Simpson(double *x, int i0, int i1) {
  int i, k;
  double a, b;

  b = x[i0];
  a = 0.0;
  for (i = i0+1; i < i1; i += 2) {
    a += x[i];
  }
  b += 4.0*a;
  a = 0.0;
  k = i1-1;
  for (i = i0+2; i < k; i += 2) {
    a += x[i];
  }
  b += 2.0*a;
  if (i == i1) {
    b += x[i1];
    b /= 3.0;
  } else {
    b += x[k];
    b /= 3.0;
    b += 0.5*(x[k] + x[i1]);
  }

  return b;
}

/*
 * Integration by Newton-Cotes formula
 * input: x[]
 * limits (indices): 0 .. ilast
 * output: r[]
 * last_only: only the finite integral is needed
 * id - direction (<0 => inward)
 */
int NewtonCotes(double r[], const double x[], int ilast,
                int last_only, int id) {
  int i, k;
  double a;

  if (id >= 0) {
    if (last_only) {
      r[ilast] = x[0];
      a = 0.0;
      for (i = 1; i < ilast; i += 2) {
	a += x[i];
      }
      r[ilast] += 4.0*a;
      a = 0.0;
      k = ilast-1;
      for (i = 2; i < k; i += 2) {
	a += x[i];
      }
      r[ilast] += 2.0*a;
      if (i == ilast) {
	r[ilast] += x[ilast];
	r[ilast] /= 3.0;
      } else {
	r[ilast] += x[k];
	r[ilast] /= 3.0;
	r[ilast] += 0.5*(x[k] + x[ilast]);
      }
      r[ilast] += r[0];
    } else {
      for (i = 2; i <= ilast; i += 2) {
	r[i] = r[i-2] + _CNC[2][0]*(x[i-2]+x[i]) + _CNC[2][1]*x[i-1];
	r[i-1] = r[i-2] + 0.5*(x[i-2] + x[i-1]);
      }
      if (i == ilast+1) {
	k = ilast-1;
	r[ilast] = r[k] + 0.5*(x[k]+x[ilast]);
      }
    }
  } else {
    if (last_only) {
      r[ilast] = x[ilast];
      a = 0.0;
      for (i = ilast-1; i > 0; i -= 2) {
	a += x[i];
      }
      r[ilast] += 4.0*a;
      a = 0.0;
      k = 1;
      for (i = ilast-2; i > k; i -= 2) {
	a += x[i];
      }
      r[ilast] += 2.0*a;
      if (i == 0) {
	r[ilast] += x[0];
	r[ilast] /= 3.0;
      } else {
	r[ilast] += x[k];
	r[ilast] /= 3.0;
	r[ilast] += 0.5*(x[k] + x[0]);
      }
      r[ilast] += r[0];
    } else {
      for (i = ilast-2; i >= 0; i -= 2) {
	r[i] = r[i+2] + _CNC[2][0]*(x[i+2]+x[i]) + _CNC[2][1]*x[i+1];
	r[i+1] = r[i+2] + 0.5*(x[i+2] + x[i+1]);
      }
      if (i == -1) {
	k = 1;
	r[0] = r[k] + 0.5*(x[k]+x[0]);
      }
    }
  }

  return 0;
}

/* Imitate F77 UVIP3P by Hiroshi Akima (NP=3) */
void uvip3p(int nd, const double *xd, const double *yd,
		 int ni, const double *xi, double *yi)
{
    int i;
    const gsl_interp_type *type;
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double xmin, xmax;
    
    switch (nd) {
    case 2:
        type = gsl_interp_linear;
        break;
    case 3:
    case 4:
        type = gsl_interp_cspline;
        break;
    default:
        type = gsl_interp_akima;
        break;
    }
    
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(type, nd);

    gsl_spline_init(spline, xd, yd, nd);

    xmin = xd[0];
    xmax = xd[nd - 1];
    
    for (i = 0; i < ni; i++) {
        if (xi[i] < xmin && nd >= 2) {
            yi[i] = yd[0];
        } else
        if (xi[i] > xmax && nd >= 2) {
            yi[i] = yd[nd - 1];
        } else {
            int gsl_status = gsl_spline_eval_e(spline, xi[i], acc, &yi[i]);
            if (gsl_status != GSL_SUCCESS) {
                fprintf(stderr,
                    "gsl_spline_eval_e failed with status %d\n", gsl_status);
                abort();
            }
        }
    }
    
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return;
}

/* private GSL structures... */
typedef struct
{
  double * c;
  double * g;
  double * diag;
  double * offdiag;
} cspline_state_t;

typedef struct
{
  double * b;
  double * c;
  double * d;
  double * _m;
} akima_state_t;

static inline void
coeff_calc (const double c_array[], double dy, double dx, size_t index,  
            double * b, double * c, double * d)
{
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  *c = c_i;
  *d = (c_ip1 - c_i) / (3.0 * dx);
}

/* Imitate F77 UVIP3C (NP=3) */
void uvip3c(int nd, const double xd[], const double yd[],
		 double c1[], double c2[], double c3[])
{
    int i;
    const gsl_interp_type *type;
    gsl_spline *spline;
    gsl_interp *interp;
    
    
    if (nd > 4) {
        type = gsl_interp_akima;
    } else
    if (nd > 2) {
        type = gsl_interp_cspline;
    } else
    if (nd == 2) {
        type = gsl_interp_linear;
    } else {
        fprintf(stderr, "uvip3c() called with nd = %d\n", nd);
        abort();
    }
    
    spline = gsl_spline_alloc(type, nd);

    gsl_spline_init(spline, xd, yd, nd);
    
    interp = spline->interp;
    
    if (nd > 4) {
        const akima_state_t* as = (akima_state_t *) interp->state;
        memcpy(c1, as->b, sizeof(double)*(nd - 1));
        memcpy(c2, as->c, sizeof(double)*(nd - 1));
        memcpy(c3, as->d, sizeof(double)*(nd - 1));
    } else
    if (nd > 2) {
        const cspline_state_t* cs = (cspline_state_t *) interp->state;
        for (i = 0; i < nd - 1; i++) {
            double dx = xd[i + 1] - xd[i];
            double dy = yd[i + 1] - yd[i];
            
            coeff_calc(cs->c, dy, dx, i, &c1[i], &c2[i], &c3[i]);
        }
    } else {
        c1[0] = (yd[1] - yd[0])/(xd[1] - xd[0]);
        c2[0] = 0.0;
        c3[0] = 0.0;
    }

    gsl_spline_free(spline);

    return;
}
