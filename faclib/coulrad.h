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
