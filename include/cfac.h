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

#ifndef __CFAC_H_
#define __CFAC_H_

/* Versioning */
#define CFAC_VERSION        1
#define CFAC_SUBVERSION     5
#define CFAC_SUBSUBVERSION  0

#define CFAC_SUCCESS    0
#define CFAC_FAILURE    1

typedef struct _cfac_t cfac_t;

/* cfac.c */
cfac_t *
cfac_new(void);

void 
cfac_free(cfac_t *cfac);

/* nucleus.c */
int
cfac_set_atom(cfac_t *cfac, const char *s, int z, double mass, double rn);
unsigned int
cfac_get_atomic_number(const cfac_t *cfac);
double
cfac_get_atomic_mass(const cfac_t *cfac);
double
cfac_get_atomic_rn(const cfac_t *cfac);
const char *
cfac_get_atomic_symbol(const cfac_t *cfac);
double
cfac_get_nucleus_potential(const cfac_t *cfac, double r);

/* config.c */
int
cfac_add_config(cfac_t *cfac, const char *gname, const char *cfg_str, int uta);
int
cfac_get_config_gid(const cfac_t *cfac, const char *cname);
int
cfac_set_uta(cfac_t *cfac, int uta);

/* structure.c */
int
cfac_get_num_levels(const cfac_t *cfac);
int
cfac_calculate_structure(cfac_t *cfac,
    int ng, const int *gids, int npg, const int *pgids, int no_ci);

/* transition.c */
typedef struct {
    unsigned int ii, fi; /* initial (upper) and final (lower) level indices */
    double rme;          /* reduced matrix element                          */
    
    double uta_de;       /* UTA energy shift                                */
    double uta_sd;       /* UTA std. dev.                                   */
} cfac_rtrans_data_t;

typedef int
(*cfac_tr_sink_t)(const cfac_t *cfac,
    const cfac_rtrans_data_t *rtdata, void *udata);

int
crac_calculate_rtrans(cfac_t *cfac,
    unsigned nlow, unsigned *low, unsigned nup, unsigned *up,
    int mpole, int mode,
    cfac_tr_sink_t sink, void *udata);

#endif /* __CFAC_H_ */
