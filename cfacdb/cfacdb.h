/* 
 * Copyright (C) 2013-2015 Evgeny Stambulchik
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef _CFACDB_H
#define _CFACDB_H

#define CFACDB_FALSE 0
#define CFACDB_TRUE  1

#define CFACDB_SUCCESS 0
#define CFACDB_FAILURE 1

/* TODO: include a public cfac header when it exists... */
#define CFACDB_CS_CE    1
#define CFACDB_CS_CI    2
#define CFACDB_CS_PI    3

typedef struct _cfacdb_t cfacdb_t;

typedef struct {
    unsigned long ndim;     /* number of levels                */
    unsigned long rtdim;    /* number of radiative transitions */
    unsigned long aidim;    /* number of AI transitions        */
    unsigned long cedim;    /* number of CE transitions        */
    unsigned long cidim;    /* number of CI transitions        */
    unsigned long pidim;    /* number of PI transitions        */
} cfacdb_stats_t;

typedef struct {
    unsigned long int sid;
    const char *sym;
    unsigned int anum;
    double mass;
    unsigned int nele_min;
    unsigned int nele_max;
    int uta;
} cfacdb_sessions_data_t;

typedef struct {
    unsigned int nele;
    double e_gs;
    unsigned long nlevels;
} cfacdb_cstates_data_t;

typedef struct {
    unsigned int i;
    
    unsigned int ifac;
    
    double energy;
    unsigned int nele;
    unsigned int g;
    unsigned int vn, vl, p;
    const char *name, *ncmplx, *sname;    
} cfacdb_levels_data_t;

typedef struct {
    unsigned int ii, fi;
    
    int mpole;
    
    double de;
    
    double gf;
    
    double uta_de;
    double uta_sd;
} cfacdb_rtrans_data_t;

typedef struct {
    unsigned int ii, fi;
    
    double rate;
} cfacdb_aitrans_data_t;

typedef struct {
    unsigned int ii, fi;
    
    unsigned int type;
    
    double de;
    unsigned int kl;
    double ap0, ap1, ap2, ap3;
    
    unsigned int nd;
    double *e;
    double *d;
} cfacdb_ctrans_data_t;

typedef struct {
    unsigned int ii, fi;
    
    unsigned int type;
    
    double de;
    
    double ratec;
} cfacdb_crates_data_t;

typedef struct {
    unsigned int ndata; /* number of data points */
    double      *e;     /* energy grid           */
    double      *d;     /* data                  */

    double       ap[5]; /* asymptote parameters  */

    double       d0;    /* threshold limit       */
    double       let;   /* low-e tangent         */

    int          cube;  /* cubic interpolation   */
    double     (*f_asymptote)(double x, const double *ap);
} cfacdb_intext_t;

typedef int (*cfacdb_sessions_sink_t)(const cfacdb_t *cdb,
    cfacdb_sessions_data_t *cbdata, void *udata);
typedef int (*cfacdb_cstates_sink_t)(const cfacdb_t *cdb,
    cfacdb_cstates_data_t *cbdata, void *udata);
typedef int (*cfacdb_levels_sink_t)(const cfacdb_t *cdb,
    cfacdb_levels_data_t *cbdata, void *udata);
typedef int (*cfacdb_rtrans_sink_t)(const cfacdb_t *cdb,
    cfacdb_rtrans_data_t *cbdata, void *udata);
typedef int (*cfacdb_aitrans_sink_t)(const cfacdb_t *cdb,
    cfacdb_aitrans_data_t *cbdata, void *udata);
typedef int (*cfacdb_ctrans_sink_t)(const cfacdb_t *cdb,
    cfacdb_ctrans_data_t *cbdata, void *udata);
typedef int (*cfacdb_crates_sink_t)(const cfacdb_t *cdb,
    cfacdb_crates_data_t *cbdata, void *udata);


cfacdb_t *cfacdb_open(const char *fname);
void cfacdb_close(cfacdb_t *cdb);
int cfacdb_init(cfacdb_t *cdb, unsigned long sid, int nele_min, int nele_max);

int cfacdb_set_udata(cfacdb_t *cdb, void *udata);
void *cfacdb_get_udata(cfacdb_t *cdb);

unsigned int cfacdb_get_nsessions(const cfacdb_t *cdb);

int cfacdb_sessions(const cfacdb_t *cdb,
    cfacdb_sessions_sink_t sink, void *udata);

int cfacdb_get_species(const cfacdb_t *cdb, unsigned int *anum, double *mass);
int cfacdb_get_stats(const cfacdb_t *cdb, cfacdb_stats_t *stats);

int cfacdb_cstates(cfacdb_t *cdb, cfacdb_cstates_sink_t sink, void *udata);
int cfacdb_levels(cfacdb_t *cdb, cfacdb_levels_sink_t sink, void *udata);
int cfacdb_rtrans(cfacdb_t *cdb, cfacdb_rtrans_sink_t sink, void *udata);
int cfacdb_aitrans(cfacdb_t *cdb, cfacdb_aitrans_sink_t sink, void *udata);
int cfacdb_ctrans(cfacdb_t *cdb, cfacdb_ctrans_sink_t sink, void *udata);

int cfacdb_crates(cfacdb_t *cdb,
    double T, cfacdb_crates_sink_t sink, void *udata);

int cfacdb_prepare_intext(const cfacdb_t *cdb,
    const cfacdb_ctrans_data_t *cbdata, cfacdb_intext_t *intext);
double cfacdb_intext(const cfacdb_intext_t *intext, double x);

#endif /* _CFACDB_H */
