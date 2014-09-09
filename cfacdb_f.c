/* F77 wrapper functions for CFACDB */

/* 
 * Copyright (C) 2013 Evgeny Stambulchik
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

#include <stdlib.h>
#include <string.h>

#include "cfacdb.h"

typedef void (*cfacdb_levels_fsink_t)(int *i,
    double *e, int *nele, int *g, int *vn, int *vl, int *p,
    char *name, char *ncmplx, char *sname,
    int _namelen, int _ncmplxlen, int _snamelen);  

typedef void (*cfacdb_rtrans_fsink_t)(int *i, int *j,
    int *mpole, double *gf);  

typedef void (*cfacdb_aitrans_fsink_t)(int *i, int *j,
    double *rate);  

typedef void (*cfacdb_ctrans_fsink_t)(int *i, int *j,
    int *type, int *kl, double *ap0, double *ap1, double *ap2, double *ap3,
    int *nd, double *e, double *d);  

typedef void (*cfacdb_crates_fsink_t)(int *i, int *j,
    int *type, double *ratec);  


static char *f77char2str(const char *fstr, unsigned int fstrlen)
{
    char *s = malloc(fstrlen + 1);
    if (s) {
        unsigned int i;
        
        strncpy(s, fstr, fstrlen);
        s[fstrlen] = '\0';
        
        for (i = fstrlen - 1; i >= 0; i--) {
            if (s[i] == ' ') {
                s[i] = '\0';
            } else {
                break;
            }
        }
    }
    
    return s;
}


/* Fortran API functions below */
static cfacdb_t *cdb = NULL;


void cfacdb_init_(const char *fname, int *nele_min, int *nele_max,
    int *ndim, int *rtdim, int *aidim, int *cedim, int *cidim, int *pidim,
    int *ierr,
    int _fnamelen)
{
    char *s = f77char2str(fname, _fnamelen);
    
    *ierr = 0;
    
    if (!s) {
        *ierr = 1;
        return;
    }
    
    cdb = cdb_init(s, *nele_min, *nele_max);
    free(s);
    
    if (!cdb) {
        *ierr = 1;
        return;
    }
    
    *ndim = cdb->ndim;
    *rtdim = cdb->rtdim;
    *aidim = cdb->aidim;
    *cedim = cdb->cedim;
    *cidim = cdb->cidim;
    *pidim = cdb->pidim;
}


void cfacdb_close_(void)
{
    cfacdb_close(cdb);
    cdb = NULL;
}


void cfacdb_species_(int *anum, double *mass, int *ierr)
{
    if (!cdb) {
        *ierr = 1;
        return;
    } else {
        *ierr = 0;
    }

    *anum = cdb->anum;
    *mass = cdb->mass;
}


static void levels_fsink(const cfacdb_t *cdb,
    levels_cb_data_t *cbdata, void *udata)
{
    struct {
        cfacdb_levels_fsink_t sink;
    } *fdata = udata;
    
    int fi = cbdata->i + 1;
    
    fdata->sink(&fi, &cbdata->energy, (int *) &cbdata->nele, (int *) &cbdata->g,
                (int *) &cbdata->vn, (int *) &cbdata->vl, (int *) &cbdata->p,
                (char *) cbdata->name, (char *) cbdata->ncmplx, (char *) cbdata->sname,
                strlen(cbdata->name), strlen(cbdata->ncmplx),
                strlen(cbdata->sname));
}

void cfacdb_levels_(cfacdb_levels_fsink_t sink, int *ierr)
{

    struct {
        cfacdb_levels_fsink_t sink;
    } fdata;
    
    fdata.sink = sink;
    
    *ierr = cfacdb_levels(cdb, levels_fsink, &fdata);
}


static void rtrans_fsink(const cfacdb_t *cdb,
    rtrans_cb_data_t *cbdata, void *udata)
{
    struct {
        cfacdb_rtrans_fsink_t sink;
    } *fdata = udata;
    
    int fi = cbdata->ii + 1;
    int fj = cbdata->fi + 1;
    
    fdata->sink(&fi, &fj, &cbdata->mpole, &cbdata->gf);
}

void cfacdb_rtrans_(cfacdb_rtrans_fsink_t sink, int *ierr)
{

    struct {
        cfacdb_rtrans_fsink_t sink;
    } fdata;
    
    fdata.sink = sink;
    
    *ierr = cfacdb_rtrans(cdb, rtrans_fsink, &fdata);
}


static void aitrans_fsink(const cfacdb_t *cdb,
    aitrans_cb_data_t *cbdata, void *udata)
{
    struct {
        cfacdb_aitrans_fsink_t sink;
    } *fdata = udata;
    
    int fi = cbdata->ii + 1;
    int fj = cbdata->fi + 1;
    
    fdata->sink(&fi, &fj, &cbdata->rate);
}

void cfacdb_aitrans_(cfacdb_aitrans_fsink_t sink, int *ierr)
{

    struct {
        cfacdb_aitrans_fsink_t sink;
    } fdata;
    
    fdata.sink = sink;
    
    *ierr = cfacdb_aitrans(cdb, aitrans_fsink, &fdata);
}


static void ctrans_fsink(const cfacdb_t *cdb,
    ctrans_cb_data_t *cbdata, void *udata)
{
    struct {
        cfacdb_ctrans_fsink_t sink;
    } *fdata = udata;
    
    int fi, fj;
    
    fi = cbdata->ii + 1;
    fj = cbdata->fi + 1;
    
    fdata->sink(&fi, &fj, (int *) &cbdata->type, (int *) &cbdata->kl,
        &cbdata->ap0, &cbdata->ap1, &cbdata->ap2, &cbdata->ap3,
        (int *) &cbdata->nd, cbdata->e, cbdata->d);
}

void cfacdb_ctrans_(cfacdb_ctrans_fsink_t sink, int *ierr)
{
    struct {
        cfacdb_ctrans_fsink_t sink;
    } fdata;
    
    fdata.sink = sink;
    
    *ierr = cfacdb_ctrans(cdb, ctrans_fsink, &fdata);
}


static void crates_fsink(const cfacdb_t *cdb,
    crates_cb_data_t *cbdata, void *udata)
{
    struct {
        cfacdb_crates_fsink_t sink;
    } *fdata = udata;
    
    int fi, fj;
    
    fi = cbdata->ii + 1;
    fj = cbdata->fi + 1;
    
    fdata->sink(&fi, &fj, (int *) &cbdata->type, &cbdata->ratec);
}


void cfacdb_crates_(double *T, cfacdb_crates_fsink_t sink, int *ierr)
{
    struct {
        cfacdb_crates_fsink_t sink;
    } fdata;
    
    fdata.sink = sink;
    
    *ierr = cfacdb_crates(cdb, *T, crates_fsink, &fdata);
}
