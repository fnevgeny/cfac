/* 
 * Copyright (C) 2013-2014 Evgeny Stambulchik
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

#include <sqlite3.h>

/* TODO: include a public cfac header when it exists... */
#define DB_SQL_CS_CE    1
#define DB_SQL_CS_CI    2
#define DB_SQL_CS_PI    3

typedef struct {
    sqlite3 *db;
    int db_format;
    
    int nele_min;
    int nele_max;
    
    unsigned long int sid;
    
    unsigned long ndim;
    unsigned long rtdim;
    unsigned long aidim;
    unsigned long cedim;
    unsigned long cidim;
    unsigned long pidim;
    
    unsigned int anum;
    double mass;
    
    unsigned int id_min;
    unsigned int id_max;

    unsigned int *lmap;
} cfacdb_t;

typedef struct {
    unsigned int nele;
    double e_gs;
    unsigned long nlevels;
} cstates_cb_data_t;

typedef struct {
    unsigned int i;
    
    unsigned int ifac;
    
    double energy;
    unsigned int nele;
    unsigned int g;
    unsigned int vn, vl, p;
    const char *name, *ncmplx, *sname;    
} levels_cb_data_t;

typedef struct {
    unsigned int ii, fi;
    
    int mpole;
    
    double gf;
} rtrans_cb_data_t;

typedef struct {
    unsigned int ii, fi;
    
    double rate;
} aitrans_cb_data_t;

typedef struct {
    unsigned int ii, fi;
    
    unsigned int type;
    
    double de;
    unsigned int kl;
    double ap0, ap1, ap2, ap3;
    
    unsigned int nd;
    double *e;
    double *d;
} ctrans_cb_data_t;

typedef struct {
    unsigned int ii, fi;
    
    unsigned int type;
    
    double de;
    
    double ratec;
} crates_cb_data_t;

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

typedef void (*cfacdb_crates_sink_t)(const cfacdb_t *cdb,
    crates_cb_data_t *cbdata, void *udata);


cfacdb_t *cdb_init(const char *fname, int nele_min, int nele_max);
void cfacdb_close(cfacdb_t *cdb);

int cfacdb_cstates(cfacdb_t *cdb,
    void (*sink)(const cfacdb_t *cdb, cstates_cb_data_t *cbdata, void *udata),
    void *udata);

int cfacdb_levels(cfacdb_t *cdb,
    void (*sink)(const cfacdb_t *cdb, levels_cb_data_t *cbdata, void *udata),
    void *udata);

int cfacdb_rtrans(cfacdb_t *cdb,
    void (*sink)(const cfacdb_t *cdb, rtrans_cb_data_t *cbdata, void *udata),
    void *udata);

int cfacdb_aitrans(cfacdb_t *cdb,
    void (*sink)(const cfacdb_t *cdb, aitrans_cb_data_t *cbdata, void *udata),
    void *udata);

int cfacdb_ctrans(cfacdb_t *cdb,
    void (*sink)(const cfacdb_t *cdb, ctrans_cb_data_t *cbdata, void *udata),
    void *udata);

int cfacdb_prepare_intext(const cfacdb_t *cdb,
    const ctrans_cb_data_t *cbdata, cfacdb_intext_t *intext);
double cfacdb_intext(const cfacdb_intext_t *intext, double x);

int cfacdb_crates(cfacdb_t *cdb, double T,
    cfacdb_crates_sink_t sink,
    void *udata);

#endif /* _CFACDB_H */
