/*
 * Callback-based access to CFAC SQLite DB, based on corresponding
 * code from CRaC.
 */

/*
 * Copyright (C) 2015 Evgeny Stambulchik
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

#include "cfacdbP.h"

#include "cache_schema.i"

int cfacdb_attach_cache(cfacdb_t *cdb, const char *fname)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc, i;
    char *errmsg;

    if (!fname) {
        return CFACDB_FAILURE;
    }

    /* open the cache DB */
    rc = sqlite3_open(fname, &cdb->cache_db);
    if (rc) {
        fprintf(stderr, "Cannot open '%s' database: %s\n",
            fname, sqlite3_errmsg(cdb->cache_db));
        return CFACDB_FAILURE;
    }

    /* ... and initialize it if needed */
    i = 0;
    while ((sql = cache_schema[i])) {
        rc = sqlite3_exec(cdb->cache_db, sql, NULL, NULL, &errmsg);
        if (rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", errmsg);
            sqlite3_free(errmsg);
            return CFACDB_FAILURE;
        }
        i++;
    }

    /* now attach it to the primary database */
    sql = "ATTACH DATABASE ? AS cache";
    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_text(stmt, 1, fname, -1, SQLITE_STATIC);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        return CFACDB_FAILURE;
    }
    sqlite3_finalize(stmt);

    /* we're ready for cached mode */
    cdb->cached = CFACDB_TRUE;

    return CFACDB_SUCCESS;
}

int cfacdb_crates_cached(cfacdb_t *cdb, double T,
    cfacdb_crates_sink_t sink, void *udata)
{
    sqlite3_stmt *stmt;
    const char *sql;
    int rc;

    sql = "INSERT INTO _cache_temp" \
          " SELECT DISTINCT (cid) cid, rate" \
          "  FROM crates" \
          "  WHERE t >= ? AND t <= ? AND cid IN (" \
          "   SELECT cid FROM _ctransitions_v" \
          "    WHERE sid = ? AND ini_nele <= ? AND fin_nele >= ?)";


    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_double(stmt, 1, 0.99*T);
    sqlite3_bind_double(stmt, 2, 1.01*T);

    sqlite3_bind_int(stmt, 3, cdb->sid);
    sqlite3_bind_int(stmt, 4, cdb->nele_max);
    sqlite3_bind_int(stmt, 5, cdb->nele_min);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE) {
        fprintf(stderr, "SQL error %d: %s\n", rc, sqlite3_errmsg(cdb->db));
        sqlite3_finalize(stmt);
        return CFACDB_FAILURE;
    }
    sqlite3_finalize(stmt);

    sql = "SELECT ini_id, fin_id, type, de, rate" \
          " FROM _crates_cached_v" \
          " ORDER BY ini_id, fin_id, type";

    sqlite3_prepare_v2(cdb->db, sql, -1, &stmt, NULL);
    sqlite3_bind_int(stmt, 1, cdb->sid);
    sqlite3_bind_int(stmt, 2, cdb->nele_max);
    sqlite3_bind_int(stmt, 3, cdb->nele_min);

    do {
        double rate, de;
        unsigned int ilfac, iufac, type;

        cfacdb_crates_data_t cbdata;

        rc = sqlite3_step(stmt);
        switch (rc) {
        case SQLITE_DONE:
        case SQLITE_OK:
            break;
        case SQLITE_ROW:
            iufac = sqlite3_column_int   (stmt, 0);
            ilfac = sqlite3_column_int   (stmt, 1);

            type  = sqlite3_column_int   (stmt, 2);
            de    = sqlite3_column_double(stmt, 3);

            rate  = sqlite3_column_double(stmt, 4);

            cbdata.ii = cdb->lmap.map[iufac - cdb->lmap.id_min];
            cbdata.fi = cdb->lmap.map[ilfac - cdb->lmap.id_min];

            cbdata.type = type;
            cbdata.de   = de;

            cbdata.ratec = rate;

            if (sink(cdb, &cbdata, udata) != CFACDB_SUCCESS) {
                sqlite3_finalize(stmt);
                return CFACDB_FAILURE;
            }

            break;
        default:
            fprintf(stderr, "SQL error %d: %s\n", rc, sqlite3_errmsg(cdb->db));
            sqlite3_finalize(stmt);
            return CFACDB_FAILURE;
            break;
        }
    } while (rc == SQLITE_ROW);

    sqlite3_finalize(stmt);

    return CFACDB_SUCCESS;
}
