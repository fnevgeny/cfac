#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_const.h>

#include "coulrad.h"

/* #define BIEDENHARN_LOUCK_SIGNS */

static double A(unsigned n, unsigned l)
{
    if (n*l == 0) {
        abort();
    }
    
    return sqrt(n*n - l*l)/(n*l);
}

static double R1_last(unsigned n, unsigned np)
{
    unsigned nd, ns;
    double a, b, c, R;

    nd = n - np;
    ns = n + np;

    a = (2*np + 2.5)*M_LN2;
    b = gsl_sf_lnchoose(ns, nd) + log(nd*np);
    c = (np + 2.0)*log(n*np) + (nd - 2.0)*log(nd) - (ns + 2.0)*log(ns);

    R = exp(a + 0.5*b + c);

#ifdef BIEDENHARN_LOUCK_SIGNS
    if (nd%2 == 0) {
        R = -R;
    }
#endif
    
    return R;
}

gsl_coulomb_me *gsl_coulomb_me_alloc(unsigned int n, unsigned int np)
{
    gsl_coulomb_me *r;
    unsigned int l;
    
    if (np > n) {
        GSL_ERROR_VAL("n' > n", GSL_EINVAL, NULL);
    }

    r = malloc(sizeof(gsl_coulomb_me));
    if (!r) {
        GSL_ERROR_VAL("failed to allocate struct", GSL_ENOMEM, NULL);
    }
    
    r->np = np;
    r->n  = n;
    
    if (n == np) {
        /* no need to allocate table */
        r->table = NULL;
        return r;
    }
    
    r->table = malloc((2*np)*sizeof(double));
    if (!r->table) {
        free(r);
        GSL_ERROR_VAL("failed to allocate table", GSL_ENOMEM, NULL);
    }


    r->table[np - 1] = R1_last(n, np);
    r->table[2*np - 1] = 0.0;
    
    for (l = np - 1; l > 0; l--) {
        r->table[l - 1]  = ((2*l + 1)*A(n, l + 1)*r->table[l] +
                            A(np, l + 1)*r->table[np + l])/(2*l*A(np,l));
        r->table[np + l - 1] = ((2*l + 1)*A(np, l + 1)*r->table[np + l] +
                            A(n, l + 1)*r->table[l])/(2*l*A(n,l));
    } 
    
    return r;
}

void gsl_coulomb_me_free(gsl_coulomb_me *r)
{
    if (r) {
        if (r->table) {
            free(r->table);
        }
        free(r);
    }
}


double
gsl_coulomb_me_get(const gsl_coulomb_me *r, unsigned int l, unsigned int lp)
{
    unsigned int n = r->n, np = r->np;
    
    if (l >= n) {
        GSL_ERROR_VAL("l >= n", GSL_EINVAL, 0.0);
    }
    if (lp >= np) {
        GSL_ERROR_VAL("l' >= n'", GSL_EINVAL, 0.0);
    }
    
    if (l != lp + 1 && lp != l + 1) {
        return 0.0;
    }
    
    if (n == np) {
        unsigned int lmax;
        double R;
        if (l > lp) {
            lmax = l;
        } else {
            lmax = lp;
        }
        R = -1.5*n*sqrt(n*n - lmax*lmax);
#ifdef BIEDENHARN_LOUCK_SIGNS
        R = -R;
#endif
        return R;
    } else
    if (l > lp) {
        return r->table[lp];
    } else {
        return r->table[np + l];
    }
}

double gsl_coulomb_me_scale(double Z, double amass)
{
    static const double au =
        GSL_CONST_MKS_UNIFIED_ATOMIC_MASS/GSL_CONST_MKS_MASS_ELECTRON;

    double scale;
    
    if (Z <= 0.0) {
        GSL_ERROR_VAL("Z is not positive", GSL_EDOM, 0.0);
    }
    if (amass < 0.0) {
        GSL_ERROR_VAL("mass is negative", GSL_EDOM, 0.0);
    }
    
    scale = 1/Z;
    
    if (amass) {
        scale *= 1.0 + 1.0/(au*amass);
    }
    
    return scale;
}
