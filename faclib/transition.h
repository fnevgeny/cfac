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

#ifndef _TRANSITION_H_
#define _TRANSITION_H_

typedef struct _TRANSITION_ {
    int nup;
    int nlo;
    LEVEL *lup;
    LEVEL *llo;
    double e;
} TRANSITION;

int GetTransition(const cfac_t *cfac,
    int nlo, int nup, TRANSITION *tr, int *swapped);
void SetTransitionMode(cfac_t *cfac, int m);
void SetTransitionGauge(cfac_t *cfac, int m);
void SetTransitionOptions(cfac_t *cfac, int gauge, int mode);
int GetTransitionGauge(const cfac_t *cfac);
int GetTransitionMode(const cfac_t *cfac);
int TRMultipole(cfac_t *cfac, double *rme, double *energy,
                int m, int low, int up);
int OverlapLowUp(int nlow, int *low, int nup, int *up);
double OscillatorStrength(int m, double e, double s, double *ga);
int SaveTransitionEB(cfac_t *cfac, int nlow, int *low, int nup, int *up,
                     char *fn, int multipole);
int GetLowUpEB(const cfac_t *cfac, int *nlow, int **low, int *nup, int **up,
               int nlow0, const int *low0, int nup0, const int *up0);

int TRMultipoleEB(cfac_t *cfac, cfac_w3j_cache_t *w3j_cache,
    double *strength, double *energy, int m, int lower, int upper);

#endif
