/*
 * Callback-based access to CFAC SQLite DB, based on corresponding
 * code from CRaC.
 */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_const_num.h>
#define ALPHA GSL_CONST_NUM_FINE_STRUCTURE

#include "cfacdbP.h"

#define SQR(a) ((a)*(a))

#include "cfac_schema.i"
#include "cfac_schema_v1.i"
#include "cfac_schema_v2.i"
#include "cfac_schema_v3.i"
#include "cfac_schema_v4.i"

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

static int format_cb(void *udata,
    int argc, char **argv, char **colNames)
{
    int *db_format = udata;

    if (argc != 1 || !argv[0]) {
        return -1;
    }

    *db_format = atol(argv[0]);

    return 0;
}

static int session_cb(void *udata,
    int argc, char **argv, char **colNames)
{
    cfacdb_t *cdb = udata;

    if (argc != 1 || !argv[0]) {
        return -1;
    }

    cdb->nsessions = atol(argv[0]);

    return 0;
}

static int field_cb(void *udata,
    int argc, char **argv, char **colNames)
{
    cfacdb_t *cdb = udata;

    if (argc != 1 || !argv[0]) {
        return -1;
    }

    cdb->nfields = atol(argv[0]);

    return 0;
}

cfacdb_t *cfacdb_open(const char *fname, cfacdb_temp_t temp_store)
{
    cfacdb_t *cdb = NULL;

    const char *sql;
    char *errmsg;
    int rc;

    unsigned int i, ns;

    const char **schemas[3];

    cdb = malloc(sizeof(cfacdb_t));
    if (!cdb) {
        return NULL;
    }
    memset(cdb, 0, sizeof(cfacdb_t));

    rc = sqlite3_open_v2(fname, &cdb->db, SQLITE_OPEN_READONLY, NULL);
    if (rc) {
        fprintf(stderr, "Cannot open database \"%s\": %s\n",
            fname, sqlite3_errmsg(cdb->db));
        cfacdb_close(cdb);
        return NULL;
    }

    /* Set temporal storage */
    switch (temp_store) {
    case CFACDB_TEMP_FILE:
        sql = "PRAGMA temp_store = FILE";
        break;
    case CFACDB_TEMP_MEMORY:
        sql = "PRAGMA temp_store = MEMORY";
        break;
    default:
        sql = "PRAGMA temp_store = DEFAULT";
        break;
    }

    rc = sqlite3_exec(cdb->db, sql, NULL, NULL, &errmsg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", errmsg);
        sqlite3_free(errmsg);
        cfacdb_close(cdb);
        return NULL;
    }

    /* verify the version/format is compatible */
    sql = "SELECT value FROM cfacdb WHERE property = 'format'";

    rc = sqlite3_exec(cdb->db, sql, format_cb, &cdb->db_format, &errmsg);
    if (rc != SQLITE_OK) {
        /* assume the first version, without the cfacdb table */
        cdb->db_format = 1;
    }

    if (cdb->db_format < 1 || cdb->db_format > 5) {
        fprintf(stderr, "Unsupported database format %d\n", cdb->db_format);
        cfacdb_close(cdb);
        return NULL;
    }

    schemas[0] = cfac_schema;
    if (cdb->db_format < 2) {
        schemas[1] = cfac_schema_v1;
    } else {
        schemas[1] = cfac_schema_v2;
    }
    if (cdb->db_format < 4) {
        schemas[2] = cfac_schema_v3;
    } else {
        schemas[2] = cfac_schema_v4;
    }

    /* create temporary views etc */
    for (ns = 0; ns < 3; ns++) {
        const char **schema = schemas[ns];
        i = 0;
        while ((sql = schema[i])) {
            rc = sqlite3_exec(cdb->db, sql, NULL, NULL, &errmsg);
            if (rc != SQLITE_OK) {
                fprintf(stderr, "SQL error: %s\n", errmsg);
                sqlite3_free(errmsg);
                cfacdb_close(cdb);
                return NULL;
            }
            i++;
        }
    }

    sql = "SELECT COUNT(sid) FROM sessions";
    rc = sqlite3_exec(cdb->db, sql, session_cb, cdb, &errmsg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", errmsg);
        sqlite3_free(errmsg);
        cfacdb_close(cdb);
        return NULL;
    }

    if (!cdb->nsessions) {
        fprintf(stderr, "%s contains no valid data\n", fname);
        cfacdb_close(cdb);
        return NULL;
    }

    if (cdb->db_format > 4) {
        sql = "SELECT COUNT(id) FROM fields";
        rc = sqlite3_exec(cdb->db, sql, field_cb, cdb, &errmsg);
        if (rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", errmsg);
            sqlite3_free(errmsg);
            cfacdb_close(cdb);
            return NULL;
        }
    }

    return cdb;
}

int cfacdb_map_allocate(cfacdb_map_t *map,
    unsigned int id_min, unsigned int id_max)
{
    unsigned int ntot;

    if (id_min >= id_max) {
        return CFACDB_FAILURE;
    }

    ntot = id_max - id_min + 1;

    map->map = malloc(ntot*sizeof(unsigned int));
    if (!map->map) {
        fprintf(stderr, "Failed allocating memory for ntot=%u\n",
            ntot);
        return CFACDB_FAILURE;
    }
    memset(map->map, 0, ntot*sizeof(unsigned int));

    map->id_min = id_min;
    map->id_max = id_max;

    return CFACDB_SUCCESS;
}

void cfacdb_map_free(cfacdb_map_t *map)
{
    if (map && map->map) {
        free(map->map);
        map->map = NULL;
    }
}

unsigned int cfacdb_map_get(const cfacdb_map_t *map, unsigned int id)
{
    if (id >= map->id_min && id <= map->id_max) {
        return map->map[id - map->id_min];
    } else {
        return 0;
    }
}

int cfacdb_map_set(cfacdb_map_t *map, unsigned int id, unsigned int mapped_id)
{
    if (id >= map->id_min && id <= map->id_max) {
        map->map[id - map->id_min] = mapped_id;
        return CFACDB_SUCCESS;
    } else {
        return CFACDB_FAILURE;
    }
}

void cfacdb_close(cfacdb_t *cdb)
{
    if (cdb) {
        cfacdb_map_free(&cdb->lmap);
        cfacdb_map_free(&cdb->smap);

        sqlite3_close(cdb->db);

        if (cdb->cached) {
            sqlite3_close(cdb->cache_db);
        }

        free(cdb);
    }

    cdb = NULL;
}

int cfacdb_set_udata(cfacdb_t *cdb, void *udata)
{
    if (!cdb) {
        return CFACDB_FAILURE;
    }

    cdb->udata = udata;

    return CFACDB_SUCCESS;
}

void *cfacdb_get_udata(cfacdb_t *cdb)
{
    if (!cdb) {
        return NULL;
    } else {
        return cdb->udata;
    }
}

unsigned int cfacdb_get_nsessions(const cfacdb_t *cdb)
{
    if (!cdb) {
        return 0;
    } else {
        return cdb->nsessions;
    }
}

int cfacdb_sessions(const cfacdb_t *cdb,
    cfacdb_sessions_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    if (!cdb) {
        return CFACDB_FAILURE;
    }

    if (cdb->db_format <= 3) {
        sql = "SELECT sid, symbol, anum, mass, nele_min, nele_max" \
              " FROM _sessions_v" \
              " ORDER BY sid";
    } else {
        sql = "SELECT sid, symbol, anum, mass, nele_min, nele_max, tstamp" \
              " FROM _sessions_v" \
              " ORDER BY sid";
    }
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);

    do {
        cfacdb_sessions_data_t cbdata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            cbdata.sid      = sqlite3_column_int   (stmt, 0);
            cbdata.sym      = (char *) sqlite3_column_text(stmt, 1);
            cbdata.anum     = sqlite3_column_int   (stmt, 2);
            cbdata.mass     = sqlite3_column_double(stmt, 3);
            cbdata.nele_min = sqlite3_column_int   (stmt, 4);
            cbdata.nele_max = sqlite3_column_int   (stmt, 5);
            if (cdb->db_format >= 4) {
                cbdata.tstamp = sqlite3_column_int (stmt, 6);
            } else {
                cbdata.tstamp = cbdata.sid;
            }

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}

int cfacdb_init(cfacdb_t *cdb, unsigned long sid, int nele_min, int nele_max)
{
    sqlite3_stmt *stmt;
    const char *sql;
    char *errmsg;
    int rc;

    unsigned int id_min, id_max, i;

    if (!cdb) {
        return CFACDB_FAILURE;
    }

    cdb->nele_min = nele_min;
    cdb->nele_max = nele_max;

    /* free from possible previous invocation of cfacdb_init() */
    cfacdb_map_free(&cdb->lmap);
    cfacdb_map_free(&cdb->smap);

    if (sid) {
        cdb->sid = sid;
    } else {
        if (cdb->nsessions > 1) {
            fprintf(stderr,
                "Warning: the DB contains %d sessions, choosing the latest\n",
                cdb->nsessions);
        }
        /* select the latest cFAC session */
        sql = "SELECT MAX(sid) FROM sessions";

        rc = sqlite3_exec(cdb->db, sql, sid_cb, &cdb->sid, &errmsg);
        if (rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", errmsg);
            sqlite3_free(errmsg);
            return CFACDB_FAILURE;
        }
    }

    /* get dimension of the database subset */
    sql = "SELECT COUNT(*) AS ndim" \
          " FROM _levels_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?";

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        return CFACDB_FAILURE;
    }

    cdb->stats.ndim = sqlite3_column_int64(stmt, 0);
    sqlite3_finalize(stmt);
    if (cdb->stats.ndim == 0) {
        fprintf(stderr, "Empty or non-existing session id %lu\n", cdb->sid);
        return CFACDB_FAILURE;
    }

    if (cdb->db_format >= 4) {
        int nuta;

        /* get number of UTA levels */
        sql = "SELECT COUNT(*) AS nuta" \
              " FROM _levels_v" \
              " WHERE sid = ? AND nele <= ? AND nele >= ? AND uta = TRUE";

        sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, cdb->sid);
        sqlite3_bind_int(stmt, 2, nele_max);
        sqlite3_bind_int(stmt, 3, nele_min);

        rc = sqlite3_step(stmt);
        if (rc != SQLITE_ROW) {
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
        }

        nuta = sqlite3_column_int64(stmt, 0);

        sqlite3_finalize(stmt);

        if (nuta == 0) {
            cdb->stats.mode = CFACDB_MODE_DETAILED;
        } else
        if (nuta == cdb->stats.ndim) {
            cdb->stats.mode = CFACDB_MODE_UTA;
        } else {
            cdb->stats.mode = CFACDB_MODE_MIXED;
        }
    } else {
        cdb->stats.mode = CFACDB_MODE_UNKNOWN;
    }

    sql = "SELECT COUNT(sid) AS rtdim" \
          " FROM _rtransitions_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ? AND de <> 0";

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        return CFACDB_FAILURE;
    }

    cdb->stats.rtdim = sqlite3_column_int64(stmt, 0);

    sqlite3_finalize(stmt);

    sql = "SELECT COUNT(sid) AS aidim" \
          " FROM _aitransitions_v" \
          " WHERE sid = ? AND nele <= ? AND nele > ?";

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, nele_max);
    sqlite3_bind_int(stmt, 3, nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        return CFACDB_FAILURE;
    }

    cdb->stats.aidim = sqlite3_column_int64(stmt, 0);

    sqlite3_finalize(stmt);

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
        return CFACDB_FAILURE;
    }

    cdb->stats.cedim = sqlite3_column_int64(stmt, 0);

    sqlite3_finalize(stmt);

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
        return CFACDB_FAILURE;
    }

    cdb->stats.cidim = sqlite3_column_int64(stmt, 0);

    sqlite3_finalize(stmt);

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
        return CFACDB_FAILURE;
    }

    cdb->stats.pidim = sqlite3_column_int64(stmt, 0);

    sqlite3_finalize(stmt);


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
        return CFACDB_FAILURE;
    }

    cdb->anum   = sqlite3_column_int   (stmt, 1);
    cdb->mass   = sqlite3_column_double(stmt, 2);

    id_min = sqlite3_column_int64 (stmt, 3);
    id_max = sqlite3_column_int64 (stmt, 4);

    sqlite3_finalize(stmt);

    if (cfacdb_map_allocate(&cdb->lmap, id_min, id_max) != CFACDB_SUCCESS) {
        return CFACDB_FAILURE;
    }

    sql = "SELECT id" \
          " FROM levels" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?" \
          " ORDER BY nele DESC, e ASC";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    i = 0;
    do {
        int ifac;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            ifac = sqlite3_column_int64(stmt, 0);

            cfacdb_map_set(&cdb->lmap, ifac, i); i++;
            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    cdb->initialized = CFACDB_TRUE;

    return CFACDB_SUCCESS;
}

int cfacdb_get_species(const cfacdb_t *cdb, unsigned int *anum, double *mass)
{
    if (cdb->initialized) {
        *anum = cdb->anum;
        *mass = cdb->mass;
        return CFACDB_SUCCESS;
    } else {
        return CFACDB_FAILURE;
    }
}

int cfacdb_get_stats(const cfacdb_t *cdb, cfacdb_stats_t *stats)
{
    if (cdb->initialized) {
        *stats = cdb->stats;
        return CFACDB_SUCCESS;
    } else {
        return CFACDB_FAILURE;
    }
}

int cfacdb_cstates(cfacdb_t *cdb, cfacdb_cstates_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    if (!cdb || !cdb->initialized) {
        return CFACDB_FAILURE;
    }

    sql = "SELECT nele, e_gs, nlevels" \
          " FROM _cstates_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?" \
          " ORDER BY nele DESC";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    do {
        cfacdb_cstates_data_t cbdata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            cbdata.nele    = sqlite3_column_int   (stmt, 0);
            cbdata.e_gs    = sqlite3_column_double(stmt, 1);
            cbdata.nlevels = sqlite3_column_int   (stmt, 2);

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}

int cfacdb_levels(cfacdb_t *cdb, cfacdb_levels_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    unsigned int i;

    if (!cdb || !cdb->initialized) {
        return CFACDB_FAILURE;
    }

    sql = "SELECT id, name, nele, e, g, vn, vl, p, ncomplex, sname, uta" \
          " FROM _levels_v" \
          " WHERE sid = ? AND nele <= ? AND nele >= ?" \
          " ORDER BY zsp, e";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    i = 0;
    do {
        cfacdb_levels_data_t cbdata;

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
            cbdata.uta    = sqlite3_column_int   (stmt, 10);

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }
            i++;

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}


int cfacdb_rtrans(cfacdb_t *cdb, cfacdb_rtrans_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    if (!cdb || !cdb->initialized) {
        return CFACDB_FAILURE;
    }

    if (cdb->db_format < 3) {
        sql = "SELECT ini_id, fin_id, mpole, rme, de" \
              " FROM _rtransitions_v" \
              " WHERE sid = ? AND nele <= ? AND nele >= ? AND de <> 0" \
              " ORDER BY ini_id, fin_id";
    } else {
        sql = "SELECT ini_id, fin_id, mpole, rme, de, uta_de, uta_sd" \
              " FROM _rtransitions_v" \
              " WHERE sid = ? AND nele <= ? AND nele >= ? AND de <> 0" \
              " ORDER BY ini_id, fin_id";
    }

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    do {
        double de, rme, uta_de, uta_sd;
        unsigned int ilfac, iufac, m2;
        int mpole;

        cfacdb_rtrans_data_t cbdata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            ilfac = sqlite3_column_int   (stmt, 0);
            iufac = sqlite3_column_int   (stmt, 1);
            mpole = sqlite3_column_int   (stmt, 2);
            rme   = sqlite3_column_double(stmt, 3);
            de    = sqlite3_column_double(stmt, 4);
            if (cdb->db_format >= 3) {
                uta_de = sqlite3_column_double(stmt, 5);
                uta_sd = sqlite3_column_double(stmt, 6);
            } else {
                uta_de = 0.0;
                uta_sd = 0.0;
            }

            /* in the original DB format, ini/fin levels were swapped */
            if (cdb->db_format < 2) {
                unsigned int ibuf;
                ibuf = ilfac; ilfac = iufac; iufac = ibuf;
                de = -de;
            }

            m2 = 2*abs(mpole);
            cbdata.gf = SQR(rme)*de*pow(ALPHA*de, m2 - 2)/(m2 + 1);
            cbdata.mpole = mpole;
            cbdata.de = de;

            cbdata.ii = cfacdb_map_get(&cdb->lmap, ilfac);
            cbdata.fi = cfacdb_map_get(&cdb->lmap, iufac);

            cbdata.de = de;

            cbdata.uta_de = uta_de;
            cbdata.uta_sd = uta_sd;

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}


int cfacdb_aitrans(cfacdb_t *cdb, cfacdb_aitrans_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    if (!cdb || !cdb->initialized) {
        return CFACDB_FAILURE;
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

        cfacdb_aitrans_data_t cbdata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            iufac = sqlite3_column_int   (stmt, 0);
            ilfac = sqlite3_column_int   (stmt, 1);
            rate  = sqlite3_column_double(stmt, 2);

            cbdata.ii = cfacdb_map_get(&cdb->lmap, iufac);
            cbdata.fi = cfacdb_map_get(&cdb->lmap, ilfac);

            cbdata.rate = rate;

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}


int cfacdb_ctrans(cfacdb_t *cdb, cfacdb_ctrans_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int nd, rc;
    unsigned int ilfac_prev, iufac_prev, type_prev;

    double es[256], ds[256]; /* TODO: make truly allocatable ? */
    cfacdb_ctrans_data_t cbdata;

    cbdata.e = es;
    cbdata.d = ds;

    if (!cdb || !cdb->initialized) {
        return CFACDB_FAILURE;
    }

    if (cdb->db_format == 1) {
        if (cdb->cached) {
            sql = "SELECT cid, ini_id, fin_id, type, e, strength, de," \
                  "       ap0, ap1" \
                  " FROM _cstrengths_v" \
                  " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ?" \
                  " AND cid NOT IN (SELECT cid FROM _cache_temp)" \
                  " ORDER BY ini_id, fin_id, type, e";
        } else {
            sql = "SELECT cid, ini_id, fin_id, type, e, strength, de," \
                  "       ap0, ap1" \
                  " FROM _cstrengths_v" \
                  " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ?" \
                  " ORDER BY ini_id, fin_id, type, e";
        }
    } else {
        if (cdb->cached) {
            sql = "SELECT cid, ini_id, fin_id, type, e, strength, de," \
                  "       kl, ap0, ap1, ap2, ap3" \
                  " FROM _cstrengths_v" \
                  " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ?" \
                  " AND cid NOT IN (SELECT cid FROM _cache_temp)" \
                  " ORDER BY ini_id, fin_id, type, e";
        } else {
            sql = "SELECT cid, ini_id, fin_id, type, e, strength, de," \
                  "       kl, ap0, ap1, ap2, ap3" \
                  " FROM _cstrengths_v" \
                  " WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ?" \
                  " ORDER BY ini_id, fin_id, type, e";
        }
    }

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    ilfac_prev = 0;
    iufac_prev = 0;
    type_prev = 0;
    nd = 0;
    do {
        unsigned int cid, ilfac, iufac, type, kl;
        double de, ap0, ap1, ap2, ap3, e, strength;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
            if (nd) {
                cbdata.nd   = nd;

                if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                    sqlite3_finalize(stmt);
                    return CFACDB_FAILURE;
                }
            }
            break;
        case SQLITE_ROW:
            cid      = sqlite3_column_int   (stmt, 0);

            ilfac    = sqlite3_column_int   (stmt, 1);
            iufac    = sqlite3_column_int   (stmt, 2);

            type     = sqlite3_column_int   (stmt, 3);
            e        = sqlite3_column_double(stmt, 4);
            strength = sqlite3_column_double(stmt, 5);
            de       = sqlite3_column_double(stmt, 6);
            if (cdb->db_format == 1) {
                kl   = 0;
                ap0  = sqlite3_column_double(stmt, 7);
                ap1  = sqlite3_column_double(stmt, 8);
                ap2  = 0.0;
                ap3  = 0.0;
            } else {
                kl   = sqlite3_column_int   (stmt,  7);
                ap0  = sqlite3_column_double(stmt,  8);
                ap1  = sqlite3_column_double(stmt,  9);
                ap2  = sqlite3_column_double(stmt, 10);
                ap3  = sqlite3_column_double(stmt, 11);
            }

            if (ilfac != ilfac_prev ||
                iufac != iufac_prev ||
                type  != type_prev) {

                if (nd) {
                    cbdata.nd   = nd;

                    if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                        sqlite3_finalize(stmt);
                        return CFACDB_FAILURE;
                    }
                }

                nd = 0;

                ilfac_prev = ilfac;
                iufac_prev = iufac;
                type_prev  = type;
            }

            cbdata.cid  = cid;

            cbdata.ii   = cfacdb_map_get(&cdb->lmap, ilfac);
            cbdata.fi   = cfacdb_map_get(&cdb->lmap, iufac);

            cbdata.type = type;

            cbdata.de   = de;

            cbdata.kl   = kl;
            cbdata.ap0  = ap0;
            cbdata.ap1  = ap1;
            cbdata.ap2  = ap2;
            cbdata.ap3  = ap3;

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
            return CFACDB_FAILURE;
            break;
        }

    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}

unsigned int cfacdb_get_nfields(const cfacdb_t *cdb)
{
    if (!cdb) {
        return 0;
    } else {
        return cdb->nfields;
    }
}

int cfacdb_fields(const cfacdb_t *cdb, cfacdb_fields_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    if (!cdb) {
        return CFACDB_FAILURE;
    }

    if (cdb->db_format < 4) {
        return CFACDB_FAILURE;
    }

    sql = "SELECT id, bf, ef, angle" \
          " FROM fields" \
          " ORDER BY id";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);

    do {
        cfacdb_fields_data_t cbdata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            cbdata.fid   = sqlite3_column_int   (stmt, 0);
            cbdata.bf    = sqlite3_column_double(stmt, 1);
            cbdata.ef    = sqlite3_column_double(stmt, 2);
            cbdata.angle = sqlite3_column_double(stmt, 3);

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}

int cfacdb_init_field(cfacdb_t *cdb, unsigned int fid)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    unsigned int id_min, id_max, i;

    if (!cdb || !cdb->initialized || cdb->db_format < 4) {
        return CFACDB_FAILURE;
    }

    /* free from possible previous invocation of cfacdb_init_field() */
    cfacdb_map_free(&cdb->smap);

    /* get dimension of the database subset */
    sql = "SELECT MIN(id) AS id_min, MAX(id) AS id_max" \
          " FROM _states_v" \
          " WHERE sid = ? AND fid = ? AND nele <= ? AND nele >= ?";

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, fid);
    sqlite3_bind_int(stmt, 3, cdb->nele_max);
    sqlite3_bind_int(stmt, 4, cdb->nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        return CFACDB_FAILURE;
    }

    id_min = sqlite3_column_int64(stmt, 0);
    id_max = sqlite3_column_int64(stmt, 1);
    sqlite3_finalize(stmt);

    if (cfacdb_map_allocate(&cdb->smap, id_min, id_max) != CFACDB_SUCCESS) {
        return CFACDB_FAILURE;
    }

    sql = "SELECT id" \
          " FROM _states_v" \
          " WHERE sid = ? AND fid = ? AND nele <= ? AND nele >= ?" \
          " ORDER BY nele DESC, e ASC";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, fid);
    sqlite3_bind_int(stmt, 3, cdb->nele_max);
    sqlite3_bind_int(stmt, 4, cdb->nele_min);

    i = 0;
    do {
        int ifac;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            ifac = sqlite3_column_int64(stmt, 0);

            cfacdb_map_set(&cdb->smap, ifac, i); i++;
            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    cdb->initialized_fid = fid;

    return CFACDB_SUCCESS;
}

int cfacdb_states(const cfacdb_t *cdb, cfacdb_states_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;
    unsigned int i;

    if (!cdb || !cdb->initialized || !cdb->initialized_fid) {
        return CFACDB_FAILURE;
    }

    sql = "SELECT id, de, mj, name, nele, e, g, vn, vl, p, ncomplex, sname, uta" \
          " FROM _states_v" \
          " WHERE sid = ? AND fid = ? AND nele <= ? AND nele >= ?" \
          " ORDER BY nele DESC, e ASC";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->initialized_fid);
    sqlite3_bind_int(stmt, 3, cdb->nele_max);
    sqlite3_bind_int(stmt, 4, cdb->nele_min);

    i = 0;
    do {
        cfacdb_states_data_t cbdata;
        cfacdb_levels_data_t ldata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            cbdata.i     = i;
            cbdata.ifac  = sqlite3_column_int64 (stmt, 0);
            cbdata.de    = sqlite3_column_double(stmt, 1);
            cbdata.mj    = sqlite3_column_int   (stmt, 2);

            cbdata.level = &ldata;

            ldata.name   = (char *) sqlite3_column_text  (stmt, 3);
            ldata.nele   = sqlite3_column_int   (stmt, 4);
            ldata.energy = sqlite3_column_double(stmt, 5);

            ldata.g      = sqlite3_column_int   (stmt, 6);
            ldata.vn     = sqlite3_column_int   (stmt, 7);
            ldata.vl     = sqlite3_column_int   (stmt, 8);
            ldata.p      = sqlite3_column_int   (stmt, 9);
            ldata.ncmplx = (char *) sqlite3_column_text  (stmt, 10);
            ldata.sname  = (char *) sqlite3_column_text  (stmt, 11);
            ldata.uta    = sqlite3_column_int   (stmt, 12);

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
        i++;
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}

int cfacdb_rtrans_m(cfacdb_t *cdb, cfacdb_rtrans_m_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    if (!cdb || !cdb->initialized || !cdb->initialized_fid) {
        return CFACDB_FAILURE;
    }

    sql = "SELECT ini_id, fin_id, mpole, q, rme, de" \
          " FROM _rtransitions_m_v" \
          " WHERE sid = ? AND fid = ? AND nele <= ? AND nele >= ?" \
          " ORDER BY ini_id, fin_id";

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->initialized_fid);
    sqlite3_bind_int(stmt, 3, cdb->nele_max);
    sqlite3_bind_int(stmt, 4, cdb->nele_min);

    do {
        double de, rme;
        unsigned int ilfac, iufac, m2;
        int q, mpole;

        cfacdb_rtrans_m_data_t cbdata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            ilfac = sqlite3_column_int   (stmt, 0);
            iufac = sqlite3_column_int   (stmt, 1);
            mpole = sqlite3_column_int   (stmt, 2);
            q     = sqlite3_column_int   (stmt, 3);
            rme   = sqlite3_column_double(stmt, 4);
            de    = sqlite3_column_double(stmt, 5);

            m2 = 2*abs(mpole);
            cbdata.gf = SQR(rme)*de*pow(ALPHA*de, m2 - 2)/(m2 + 1);
            cbdata.mpole = mpole;
            cbdata.q  = q;
            cbdata.de = de;

            cbdata.ii = cfacdb_map_get(&cdb->smap, ilfac);
            cbdata.fi = cfacdb_map_get(&cdb->smap, iufac);

            cbdata.de = de;

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}
