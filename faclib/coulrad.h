#ifndef COULRAD_H
#define COULRAD_H

typedef struct {
    unsigned int n;
    unsigned int np;
    double *table;
} gsl_coulomb_me;

typedef struct {
    unsigned int nb;
    unsigned int ne;
    const double *e;
    double *table;
} gsl_coulomb_fb;


gsl_coulomb_me *
gsl_coulomb_me_alloc(unsigned int n, unsigned int np);

void
gsl_coulomb_me_free(gsl_coulomb_me *r);

double
gsl_coulomb_me_get(const gsl_coulomb_me *r, unsigned int l, unsigned int lp);

double
gsl_coulomb_me_scale(double Z, double amass);

gsl_coulomb_fb *
gsl_coulomb_fb_alloc(unsigned int nb, unsigned int ne, const  double *e);

void
gsl_coulomb_fb_free(gsl_coulomb_fb *r);

double
gsl_coulomb_fb_get(const gsl_coulomb_fb *rfb, unsigned int lp, unsigned int l,
    unsigned int ie);

double
gsl_coulomb_fb_get_dfdE(const gsl_coulomb_fb *rfb,
    unsigned int lb, int lk, unsigned int ie);

double
gsl_coulomb_fb_get_xs(const gsl_coulomb_fb *rfb,
    unsigned int lb, int lk, unsigned int ie);

#endif /* COULRAD_H */
