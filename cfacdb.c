/* 
 * Callback-based access to CFAC SQLite DB, based on corresponding
 * code from CRaC. Most of it will eventially become part, and distributed
 * with, CFAC
 */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_const_num.h>
#define ALPHA GSL_CONST_NUM_FINE_STRUCTURE

#include "cfacdb.h"

#define SQR(a) ((a)*(a))

#include "cfac_schema.i"


static int sid_cb(void *udata,
    int argc, char **argv, char **colNames)
{
    unsigned long int *sid = (unsigned long int *) udata;

    if (argc != 1 || !argv[0]) {
        return -1;
    }
    
    *sid = atol(argv[0]);

    return 0;
}

cfac_db_t *cdb_init(const char *fname, int nele_min, int nele_max)
{
    cfac_db_t *cdb = NULL;
    
    sqlite3_stmt *stmt;
    const char *sql;
    char *errmsg;
    int rc;
    
    unsigned int i;
    unsigned long ntot;

    cdb = malloc(sizeof(cfac_db_t));
    if (!cdb) {
        return NULL;
    }
    memset(cdb, 0, sizeof(cfac_db_t));   

    rc = sqlite3_open_v2(fname, &cdb->db, SQLITE_OPEN_READONLY, NULL);
    if (rc) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_close(cdb->db);
        return NULL;
    }

    /* create temporary views etc */
    i = 0;
    while ((sql = schema_str[i])) {
        rc = sqlite3_exec(cdb->db, sql, NULL, NULL, &errmsg);
        if (rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", errmsg);
            sqlite3_free(errmsg);
            sqlite3_close(cdb->db);
            return NULL;
        }
        i++;
    }

    /* select latest cFAC session */
    sql = "SELECT MAX(sid) FROM sessions";
    
    rc = sqlite3_exec(cdb->db, sql, sid_cb, &cdb->sid, &errmsg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", errmsg);
        sqlite3_free(errmsg);
        sqlite3_close(cdb->db);
        return NULL;
    }

    /* get dimension of the database subset */
    sql = "SELECT SUM(nlevels) AS ndim" \
          " FROM _cstates_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?";
    
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }

    cdb->ndim = sqlite3_column_int64(stmt, 0);
    
    sqlite3_reset(stmt);


    sql = "SELECT COUNT(sid) AS rtdim" \
          " FROM _rtransitions_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?";
    
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }

    cdb->rtdim = sqlite3_column_int64(stmt, 0);
    
    sqlite3_reset(stmt);

    sql = "SELECT COUNT(sid) AS aidim" \
          " FROM _aitransitions_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?";
    
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }

    cdb->aidim = sqlite3_column_int64(stmt, 0);
    
    sqlite3_reset(stmt);

    sql = "SELECT COUNT(cid) AS cedim" \
          " FROM _ctransitions_v" \
          " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ? AND type = 1";
    
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }

    cdb->cedim = sqlite3_column_int64(stmt, 0);
    
    sqlite3_reset(stmt);

    sql = "SELECT COUNT(cid) AS cidim" \
          " FROM _ctransitions_v" \
          " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ? AND type = 2";
    
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }

    cdb->cidim = sqlite3_column_int64(stmt, 0);
    
    sqlite3_reset(stmt);

    sql = "SELECT COUNT(cid) AS pidim" \
          " FROM _ctransitions_v" \
          " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ? AND type = 3";
    
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }

    cdb->pidim = sqlite3_column_int64(stmt, 0);
    
    sqlite3_reset(stmt);


    /* get element properties */
    sql = "SELECT symbol, anum, mass, id_min, id_max" \
          " FROM _species_v" \
          " WHERE sid = ?";
    
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }
    
    cdb->anum   = sqlite3_column_int   (stmt, 1);
    cdb->mass   = sqlite3_column_double(stmt, 2);
    cdb->id_min = sqlite3_column_int64 (stmt, 3);
    cdb->id_max = sqlite3_column_int64 (stmt, 4);
    
    ntot = cdb->id_max - cdb->id_min + 1;
    
    cdb->nele_min = nele_min;
    cdb->nele_max = nele_max;
    
    cdb->lmap = malloc(ntot*sizeof(unsigned int));
    if (!cdb->lmap) {
        fprintf(stderr, "Failed allocating memory for ndim=%lu\n", cdb->ndim);
        sqlite3_finalize(stmt);
        sqlite3_close(cdb->db);
        return NULL;
    }
    memset(cdb->lmap, 0, ntot*sizeof(unsigned int));

    sqlite3_finalize(stmt);
    
    return cdb;
}

void cfac_db_close(cfac_db_t *cdb)
{
    if (cdb) {
        if (cdb->lmap) {
            free(cdb->lmap);
        }
        
        sqlite3_close(cdb->db);
        
        free(cdb);
    }
    
    cdb = NULL;
}


int cfac_db_levels(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, levels_cb_data_t *cbdata, void *udata),
    void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;
    
    unsigned int i;

    if (!cdb) {
        return 1;
    }
    
    sql = "SELECT id, name, nele, e, g, vn, vl, p, ncomplex, sname" \
          " FROM _levels_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?" \
          " ORDER BY zsp, e";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    i = 0;
    do {
        levels_cb_data_t cbdata;
        
        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            cbdata.i      = i;
            cbdata.ifac   = sqlite3_column_int64 (stmt, 0);
            cbdata.name   = (char *) sqlite3_column_text  (stmt, 1);
            cbdata.nele   = sqlite3_column_int   (stmt, 2);
            cbdata.energy = sqlite3_column_double(stmt, 3);
            cbdata.g      = sqlite3_column_int   (stmt, 4);
            cbdata.vn     = sqlite3_column_int   (stmt, 5);
            cbdata.vl     = sqlite3_column_int   (stmt, 6);
            cbdata.p      = sqlite3_column_int   (stmt, 7);
            cbdata.ncmplx = (char *) sqlite3_column_text  (stmt, 8);
            cbdata.sname  = (char *) sqlite3_column_text  (stmt, 9);
            
            sink(cdb, &cbdata, udata);
            
            cdb->lmap[cbdata.ifac - cdb->id_min] = i; i++;
            
            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return 2;
            break;
        }
    } while (rc == SQLITE_ROW);
    
    sqlite3_finalize(stmt);
    
    return 0;
}


int cfac_db_rtrans(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, rtrans_cb_data_t *cbdata, void *udata),
    void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;
    
    if (!cdb) {
        return 1;
    }
    
    sql = "SELECT ini_id, fin_id, mpole, rme, de" \
          " FROM _rtransitions_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ? AND de > 0" \
          " ORDER BY ini_id, fin_id";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    do {
        double de, rme;
        unsigned int ilfac, iufac, m2;
        int mpole;
        
        rtrans_cb_data_t cbdata;
        
        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            iufac = sqlite3_column_int   (stmt, 0);
            ilfac = sqlite3_column_int   (stmt, 1);
            mpole = sqlite3_column_int   (stmt, 2);
            rme   = sqlite3_column_double(stmt, 3);
            de    = sqlite3_column_double(stmt, 4);
            
    
            m2 = 2*abs(mpole);
            cbdata.gf = SQR(rme)*de*pow(ALPHA*de, m2 - 2)/(m2 + 1);
            cbdata.mpole = mpole;
            
            cbdata.ii = cdb->lmap[ilfac - cdb->id_min];
            cbdata.fi = cdb->lmap[iufac - cdb->id_min];
            
            sink(cdb, &cbdata, udata);
            
            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return 2;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);
    
    return 0;
}


int cfac_db_aitrans(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, aitrans_cb_data_t *cbdata, void *udata),
    void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;
    
    if (!cdb) {
        return 1;
    }
    
    sql = "SELECT ini_id, fin_id, rate" \
          " FROM _aitransitions_v" \
          " WHERE sid = ? AND nele <= ? AND nele > ?" \
          " ORDER BY ini_id, fin_id";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    do {
        double rate;
        unsigned int ilfac, iufac;

        aitrans_cb_data_t cbdata;
        
        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            iufac = sqlite3_column_int   (stmt, 0);
            ilfac = sqlite3_column_int   (stmt, 1);
            rate  = sqlite3_column_double(stmt, 2);
            
            cbdata.ii = cdb->lmap[ilfac - cdb->id_min];
            cbdata.fi = cdb->lmap[iufac - cdb->id_min];
            
            cbdata.rate = rate;
            
            sink(cdb, &cbdata, udata);
            
            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return 2;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);
    
    return 0;
}


int cfac_db_ctrans(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, ctrans_cb_data_t *cbdata, void *udata),
    void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int nd, rc;
    unsigned int ilfac_prev, iufac_prev, type_prev;
    
    double es[256], ds[256]; /* TODO: make truly allocatable ? */
    ctrans_cb_data_t cbdata;
    
    cbdata.e = es;
    cbdata.d = ds;

    if (!cdb) {
        return 1;
    }
    
    sql = "SELECT ini_id, fin_id, type, e, strength, de, ap0, ap1" \
          " FROM _cstrengths_v" \
          " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ?" \
          " ORDER BY ini_id, fin_id, type, e";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);
    
    ilfac_prev = 0;
    iufac_prev = 0;
    type_prev = 0;
    nd = 0;
    do {
        unsigned int ilfac, iufac, type;
        double de, ap0, ap1, e, strength;
        
        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
            if (nd) {
                cbdata.nd   = nd;

                sink(cdb, &cbdata, udata);
            }
            break;
        case SQLITE_ROW:
            ilfac    = sqlite3_column_int(stmt, 0);
            iufac    = sqlite3_column_int(stmt, 1);
            
            type     = sqlite3_column_int(stmt, 2);
            e        = sqlite3_column_double(stmt, 3);
            strength = sqlite3_column_double(stmt, 4);
            de       = sqlite3_column_double(stmt, 5);
            ap0      = sqlite3_column_double(stmt, 6);
            ap1      = sqlite3_column_double(stmt, 7);
            
            if (ilfac != ilfac_prev ||
                iufac != iufac_prev ||
                type  != type_prev) {

                if (nd) {
                    cbdata.nd   = nd;
                    
                    sink(cdb, &cbdata, udata);
                }

                nd = 0;

                ilfac_prev = ilfac;
                iufac_prev = iufac;
                type_prev  = type;
            }
            
            cbdata.ii   = cdb->lmap[ilfac - cdb->id_min];
            cbdata.fi   = cdb->lmap[iufac - cdb->id_min];

            cbdata.type = type;

            cbdata.de   = de;

            cbdata.ap0  = ap0;
            cbdata.ap1  = ap1;

            es[nd] = 1 + e/de;
            ds[nd] = strength;
            
            nd++;

            if (strength <= 0.0) {
                fprintf(stderr,
                    "ignoring non-positive cstrength %g at e=%g\n",
                        strength, e);
                continue;
            }
            
            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return 2;
            break;
        }

    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);
    
    return 0;
}
