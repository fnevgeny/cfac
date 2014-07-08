#ifndef COULRAD_H
#define COULRAD_H

typedef struct {
    unsigned int n;
    unsigned int np;
    double *table;
} gsl_coulomb_me;

gsl_coulomb_me *
gsl_coulomb_me_alloc(unsigned int n, unsigned int np);

void
gsl_coulomb_me_free(gsl_coulomb_me *r);

double
gsl_coulomb_me_get(const gsl_coulomb_me *r, unsigned int l, unsigned int lp);

double
gsl_coulomb_me_scale(double Z, double amass);

#endif /* COULRAD_H */
