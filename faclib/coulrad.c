/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 *   Portions Copyright (C) 2010-2015 Evgeny Stambulchik
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

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


static double lnR1_fb_last0(unsigned int n)
{
    double lR;

    lR = 0.5*(log(M_PI/32) - gsl_sf_lnfact(2*n - 1)) + (n + 2)*log(4*n) - 2*n;

    return lR;
}

static double R1_fb_last(unsigned int n, double kappa)
{
    unsigned int s;
    double lR0, R, a, b, c;

    lR0 = lnR1_fb_last0(n);

    if (kappa <= 0.0) {
        return exp(lR0);
    }

    a = 0;
    for (s = 1; s <= n; s++) {
        a += log(1 + s*s*kappa*kappa);
    }
    a -= log(1 - exp(-2*M_PI/kappa));

    b = 2*n - 2/kappa*atan(n*kappa);
    c = -log(1 + n*n*kappa*kappa)*(n + 2);

    R = exp(0.5*a + b + c + lR0);

    return R;
}

static double C(double kappa, unsigned int l)
{
    return sqrt(1 + l*l*kappa*kappa)/l;
}

gsl_coulomb_fb *
gsl_coulomb_fb_alloc(unsigned int nb, unsigned int ne, const  double *e)
{
    gsl_coulomb_fb *r;
    unsigned int ie, l;

    r = malloc(sizeof(gsl_coulomb_fb));
    if (!r) {
        GSL_ERROR_VAL("failed to allocate struct", GSL_ENOMEM, NULL);
    }

    r->nb = nb;
    r->ne = ne;
    r->e  = e;

    r->table = malloc((2*nb*ne)*sizeof(double));
    if (!r->table) {
        free(r);
        GSL_ERROR_VAL("failed to allocate table", GSL_ENOMEM, NULL);
    }

    for (ie = 0; ie < ne; ie++) {
        double *t = &r->table[ie*2*nb];
        double kappa = sqrt(2*e[ie]);

        /* lp = nb - 1, lk = nb */
        t[nb - 1]   = R1_fb_last(nb, kappa);
        /* lp = nb, lk = nb - 1 */
        t[2*nb - 1] = 0.0;

        for (l = nb - 1; l > 0; l--) {
            /* lp = l - 1, lk = l */
            t[l - 1]  = ((2*l + 1)*C(kappa, l + 1)*t[l] +
                                A(nb, l + 1)*t[nb + l])/(2*l*A(nb,l));
            /* lp = l, lk = l - 1 */
            t[nb + l - 1] = ((2*l + 1)*A(nb, l + 1)*t[nb + l] +
                                C(kappa, l + 1)*t[l])/(2*l*C(kappa,l));
        }
    }

    return r;
}

void gsl_coulomb_fb_free(gsl_coulomb_fb *r)
{
    if (r) {
        if (r->table) {
            free(r->table);
        }
        free(r);
    }
}

double
gsl_coulomb_fb_get(const gsl_coulomb_fb *rfb, unsigned int lp, unsigned int l,
    unsigned int ie)
{
    unsigned int nb = rfb->nb;
    double *t;

    if (lp > nb) {
        GSL_ERROR_VAL("lp > nb", GSL_EINVAL, 0.0);
    }

    if (l != lp + 1 && lp != l + 1) {
        return 0.0;
    }

    t = &rfb->table[ie*2*nb];

    if (l > lp) {
        return t[lp];
    } else {
        return t[nb + l];
    }
}


double gsl_coulomb_fb_get_dfdE(const gsl_coulomb_fb *rfb,
    unsigned int lb, int lk, unsigned int ie)
{
    double r, S, xs, eph;
    unsigned int nb = rfb->nb;

    if (lk >= 0) {
        r = gsl_coulomb_fb_get(rfb, lb, lk, ie);
    } else {
        r = gsl_coulomb_fb_get(rfb, lb, lb + 1, ie);
    }
    S = (lb + 1)*r*r;

    if (lk < 0 && lb > 0) {
        r = gsl_coulomb_fb_get(rfb, lb, lb - 1, ie);

        S += lb*r*r;
    }

    /* photon energy */
    eph = rfb->e[ie] + 1.0/(2*nb*nb);

    xs = 4.0/3*eph/(2*lb + 1)*S;

    return xs;
}

double gsl_coulomb_fb_get_xs(const gsl_coulomb_fb *rfb,
    unsigned int lb, int lk, unsigned int ie)
{
    double dfdE = gsl_coulomb_fb_get_dfdE(rfb, lb, lk, ie);

    return 2*M_PI*GSL_CONST_NUM_FINE_STRUCTURE*dfdE;
}
