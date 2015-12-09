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

#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <sys/stat.h>

#include "cfacP.h"
#include "transition.h"
#include "dbase.h"
#include "cfacdb.h"

#include "schema.i"

#define CFACDB_FORMAT_VERSION   3

#define SQLITE3_BIND_STR(stmt, id, txt) \
        sqlite3_bind_text(stmt, id, txt, -1, SQLITE_STATIC)

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

int StoreInit(const cfac_t *cfac,
    const char *fn, int reset, sqlite3 **db, unsigned long *sid)
{
    int retval = 0;
    struct stat sb;
    sqlite3_stmt *stmt;
    int rc;
    const char *sql;
    char *errmsg;
    int need_truncate = 0, need_init = 0;
    
    *sid = (unsigned long) time(NULL);

    if (reset) {
        need_truncate = 1;
        need_init     = 1;
    }
    
    if (stat(fn, &sb) == -1) {
        if (errno == ENOENT) {
            need_truncate = 0;
            need_init     = 1;
        } else {
            perror("stat");
            return -1;
        }
    } else {
        if (sb.st_size == 0) {
            need_truncate = 0;
            need_init = 1;
        }
    }

    if (need_truncate) {
        FILE *fp = fopen(fn, "w");
        if (!fp || fclose(fp) != 0) {
            return -1;
        }
    }

    rc = sqlite3_open(fn, db);
    if (rc) {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(*db));
        sqlite3_close(*db);
        return -1;
    }

    rc = sqlite3_exec(*db, "PRAGMA foreign_keys = ON", NULL, NULL, &errmsg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", errmsg);
        sqlite3_free(errmsg);
        sqlite3_close(*db);
        return -1;
    }


    if (need_init) {
        int i = 0;
        while ((sql = schema_str[i])) {
            rc = sqlite3_exec(*db, sql, NULL, NULL, &errmsg);
            if (rc != SQLITE_OK) {
                fprintf(stderr, "SQL error: %s\n", errmsg);
                sqlite3_free(errmsg);
                sqlite3_close(*db);
                retval = -1;
                break;
            }
            i++;
        }
        
        sql = "INSERT INTO cfacdb" \
              " (property, value)" \
              " VALUES ('format', ?)";

        sqlite3_prepare_v2(*db, sql, -1, &stmt, NULL);

        sqlite3_bind_int(stmt, 1, CFACDB_FORMAT_VERSION);

        rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(*db));
            sqlite3_close(*db);
            retval = -1;
        }
        sqlite3_reset(stmt);
    } else {
        int db_format;
        sql = "SELECT value FROM cfacdb WHERE property = 'format'";
        rc = sqlite3_exec(*db, sql, format_cb, &db_format, &errmsg);
        if (rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", errmsg);
            sqlite3_free(errmsg);
            sqlite3_close(*db);
            return -1;
        }
        
        if (db_format != CFACDB_FORMAT_VERSION) {
            fprintf(stderr, "Incompatible DB format %d, expected %d\n",
                db_format, CFACDB_FORMAT_VERSION);
            sqlite3_free(errmsg);
            sqlite3_close(*db);
            return -1;
        }
    }

    sql = "INSERT INTO sessions" \
          " (sid, version, uta, cmdline)" \
          " VALUES (?, ?, ?, '')";

    sqlite3_prepare_v2(*db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt, 1, *sid);
    sqlite3_bind_int(stmt, 2,
        10000*CFAC_VERSION + 100*CFAC_SUBVERSION + CFAC_SUBSUBVERSION);
    sqlite3_bind_int(stmt, 3, cfac->uta ? 1:0);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(*db));
        sqlite3_close(*db);
        retval = -1;
    }
    sqlite3_reset(stmt);
    
    sql = "INSERT INTO species (sid, symbol, anum, mass) VALUES (?, ?, ?, ?)";
    
    sqlite3_prepare_v2(*db, sql, -1, &stmt, NULL);
    
    sqlite3_bind_int   (stmt, 1, *sid);
    SQLITE3_BIND_STR   (stmt, 2, cfac_get_atomic_symbol(cfac));
    sqlite3_bind_int   (stmt, 3, cfac_get_atomic_number(cfac));
    sqlite3_bind_double(stmt, 4, cfac_get_atomic_mass(cfac));

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(*db));
        sqlite3_close(*db);
        retval = -1;
    }
    
    sqlite3_finalize(stmt);

    return retval;
}

int StoreENTable(sqlite3 *db, unsigned long int sid, FILE *fp, int swp)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "INSERT INTO levels" \
          " (sid, id, nele, name, e, g, vn, vl, p, ibase, ncomplex, sname)" \
          " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt,  1, sid);

    sqlite3_exec(db, "BEGIN", 0, 0, 0);
    
    while (retval == 0) {
        EN_HEADER h;
        int i, n;
        
        n = ReadENHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }
        
        sqlite3_bind_int(stmt,  3, h.nele);

        for (i = 0; i < h.nlevels; i++) {
            EN_RECORD r;
            int p, vnl, vn, vl, g, ibase;
            
            n = ReadENRecord(fp, &r, swp);
            if (n == 0) {
                break;
            }

            if (r.p < 0) {
	      p = 1;
	      vnl = -r.p;
            } else {
	      p = 0;
	      vnl = r.p;
            }
            
            g = JFromENRecord(&r) + 1;
            vn = vnl/100;
            vl = vnl - 100*vn;
            
            ibase = IBaseFromENRecord(&r);
    
            sqlite3_bind_int   (stmt,  2, r.ilev);
            SQLITE3_BIND_STR   (stmt,  4, r.name);
            sqlite3_bind_double(stmt,  5, r.energy);
            sqlite3_bind_int   (stmt,  6, g);
            sqlite3_bind_int   (stmt,  7, vn);
            sqlite3_bind_int   (stmt,  8, vl);
            sqlite3_bind_int   (stmt,  9, p);
            if (ibase >= 0) {
                sqlite3_bind_int (stmt, 10, ibase);
            } else {
                sqlite3_bind_null(stmt, 10);
            }
            SQLITE3_BIND_STR   (stmt, 11, r.ncomplex);
            SQLITE3_BIND_STR   (stmt, 12, r.sname);

            rc = sqlite3_step(stmt);
            if (rc != SQLITE_DONE) {
                fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                retval = -1;
                break;
            }
            sqlite3_reset(stmt);
        }
    }

    sqlite3_exec(db, "COMMIT", 0, 0, 0);
    
    sqlite3_finalize(stmt);

    return retval;
}

int StoreTRTable(sqlite3 *db, unsigned long int sid, FILE *fp, int swp)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;
    
    sql = "INSERT INTO rtransitions" \
          " (sid, ini_id, fin_id, mpole, rme, mode)" \
          " VALUES (?, ?, ?, ?, ?, ?)";

    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt,  1, sid);
    
    sqlite3_exec(db, "BEGIN", 0, 0, 0);

    while (retval == 0) {
        TR_HEADER h;
        int n, i;
        
        n = ReadTRHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }

        /* TODO: h.gauge ? */

        for (i = 0; i < h.ntransitions; i++) {
            TR_RECORD r;
            TR_EXTRA rx;

            n = ReadTRRecord(fp, &r, &rx, swp);
            if (n == 0) {
                break;
            }
            
            sqlite3_bind_int   (stmt,  2, r.lower);
            sqlite3_bind_int   (stmt,  3, r.upper);
            sqlite3_bind_int   (stmt,  4, h.multipole);
            sqlite3_bind_double(stmt,  5, r.rme);
            sqlite3_bind_int   (stmt,  6, h.mode);
            sqlite3_bind_double(stmt,  7, rx.de);
            sqlite3_bind_double(stmt,  8, rx.sdev);

            rc = sqlite3_step(stmt);
            if (rc != SQLITE_DONE) {
                fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                retval = -1;
                break;
            }
            sqlite3_reset(stmt);
        }
    }

    sqlite3_exec(db, "COMMIT", 0, 0, 0);

    sqlite3_finalize(stmt);

    return retval;
}

int StoreAITable(sqlite3 *db, unsigned long int sid, FILE *fp, int swp)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "INSERT INTO aitransitions" \
          " (sid, ini_id, fin_id, rate)" \
          " VALUES (?, ?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt,  1, sid);
    
    sqlite3_exec(db, "BEGIN", 0, 0, 0);

    while (retval == 0) {
        AI_HEADER h;
        int n, i;

        n = ReadAIHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }

        for (i = 0; i < h.ntransitions; i++) {
            AI_RECORD r;
            
            n = ReadAIRecord(fp, &r, swp);
            if (n == 0) {
                break;
            }
            
            sqlite3_bind_int   (stmt,  2, r.b);
            sqlite3_bind_int   (stmt,  3, r.f);
            sqlite3_bind_double(stmt,  4, r.rate);

            rc = sqlite3_step(stmt);
            if (rc != SQLITE_DONE) {
                fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                retval = -1;
                break;
            }
            sqlite3_reset(stmt);
        }
    }

    sqlite3_exec(db, "COMMIT", 0, 0, 0);

    sqlite3_finalize(stmt);

    return retval;
}

static int StoreCTransitionCB(void *udata,
    int argc, char **argv, char **colNames)
{
    unsigned long int *cid = (unsigned long int *) udata;

    if (argc != 1) {
        return -1;
    }
    
    *cid = atol(argv[0]);

    return 0;
}

static int StoreCTransition(sqlite3 *db, unsigned long int sid,
    int type, int qk_mode, int ini_id, int fin_id, int kl,
    double ap0, double ap1, double ap2, double ap3,
    unsigned long int *cid)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *errmsg;
    const char *sql;
    
    *cid = 0;

    sql = "INSERT INTO ctransitions" \
          " (sid, type, qk_mode, ini_id, fin_id, kl, ap0, ap1, ap2, ap3)" \
          " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt, 1, sid);
    sqlite3_bind_int(stmt, 2, type);
    sqlite3_bind_int(stmt, 3, qk_mode);
    sqlite3_bind_int(stmt, 4, ini_id);
    sqlite3_bind_int(stmt, 5, fin_id);

    sqlite3_bind_int(stmt, 6, kl);
    
    sqlite3_bind_double(stmt, 7, ap0);
    sqlite3_bind_double(stmt, 8, ap1);
    sqlite3_bind_double(stmt, 9, ap2);
    sqlite3_bind_double(stmt, 10, ap3);
    
    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
        sqlite3_finalize(stmt);
        return -1;
    }
    sqlite3_finalize(stmt);

    sql = "SELECT MAX(cid) FROM ctransitions";

    rc = sqlite3_exec(db, sql, StoreCTransitionCB, cid, &errmsg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", errmsg);
        sqlite3_free(errmsg);
        retval = -1;
    }

    return retval;
}

int StoreCETable(sqlite3 *db, unsigned long int sid, FILE *fp, int swp)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "INSERT INTO cstrengths" \
          " (cid, e, strength)" \
          " VALUES (?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_exec(db, "BEGIN", 0, 0, 0);

    while (retval == 0) {
        CE_HEADER h;
        int n, i;

        n = ReadCEHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }

        for (i = 0; i < h.ntransitions; i++) {
            CE_RECORD r;
            unsigned long int cid;
            int k, p2;
            
            n = ReadCERecord(fp, &r, swp, &h);
            if (n == 0) {
                break;
            }

            retval = StoreCTransition(db, sid,
                CFACDB_CS_CE, QK_EXACT, r.lower, r.upper,
                0, r.bethe, r.born[0], r.born[1], 0.0,
                &cid);
            if (retval != 0) {
                break;
            }
            
            sqlite3_bind_int(stmt, 1, cid);

            p2 = 0;
            for (k = 0; k < r.nsub && retval == 0; k++) {
                int t;
                /* TODO: msub */

	        for (t = 0; t < h.n_usr; t++, p2++) {
                    sqlite3_bind_double(stmt, 2, h.usr_egrid[t]);
                    sqlite3_bind_double(stmt, 3, r.strength[p2]);
                    
                    rc = sqlite3_step(stmt);
                    if (rc != SQLITE_DONE) {
                        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                        retval = -1;
                        break;
                    }
                    sqlite3_reset(stmt);
	        }
            }
                  
            if (h.msub) {
                free(r.params);
            }
            free(r.strength);
        }

        free(h.tegrid);
        free(h.egrid);
        free(h.usr_egrid);
    }

    sqlite3_exec(db, "COMMIT", 0, 0, 0);

    sqlite3_finalize(stmt);

    return retval;
}

int StoreCITable(sqlite3 *db, unsigned long int sid, FILE *fp, int swp)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "INSERT INTO cstrengths" \
          " (cid, e, strength)" \
          " VALUES (?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_exec(db, "BEGIN", 0, 0, 0);

    while (retval == 0) {
        CI_HEADER h;
        int n, i;

        n = ReadCIHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }

        for (i = 0; i < h.ntransitions && retval == 0; i++) {
            CI_RECORD r;
            unsigned long int cid;
            int t;
            
            n = ReadCIRecord(fp, &r, swp, &h);
            if (n == 0) {
                break;
            }

            retval = StoreCTransition(db, sid,
                CFACDB_CS_CI, h.qk_mode, r.b, r.f,
                r.kl, r.params[0], r.params[1], r.params[2], r.params[3],
                &cid);
            if (retval != 0) {
                break;
            }

            sqlite3_bind_int(stmt, 1, cid);

            for (t = 0; t < h.n_usr; t++) {
                sqlite3_bind_double(stmt, 2, h.usr_egrid[t]);
                sqlite3_bind_double(stmt, 3, r.strength[t]);

                rc = sqlite3_step(stmt);
                if (rc != SQLITE_DONE) {
                    fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                    retval = -1;
                    break;
                }
                sqlite3_reset(stmt);
            }

            free(r.params); 
            free(r.strength);
        }

        free(h.tegrid);
        free(h.egrid);
        free(h.usr_egrid);
    }

    sqlite3_exec(db, "COMMIT", 0, 0, 0);

    sqlite3_finalize(stmt);

    return retval;
}

int StoreRRTable(const cfac_t *cfac,
    sqlite3 *db, unsigned long int sid, FILE *fp, int swp)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "INSERT INTO cstrengths" \
          " (cid, e, strength)" \
          " VALUES (?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt, 1, sid);
    sqlite3_bind_int(stmt, 2, CFACDB_CS_PI);
    
    sqlite3_exec(db, "BEGIN", 0, 0, 0);

    while (retval == 0) {
        RR_HEADER h;
        int n, i;

        n = ReadRRHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }

        sqlite3_bind_int(stmt, 3, h.qk_mode);

        for (i = 0; i < h.ntransitions && retval == 0; i++) {
            RR_RECORD r;
            unsigned long int cid;
            int t;
            double ap0, ap1, ap2, ap3;
            
            n = ReadRRRecord(fp, &r, swp, &h);
            if (n == 0) {
                break;
            }
            
            if (h.qk_mode == QK_FIT && r.params[0]) {
                ap0 = r.params[0];
                ap1 = r.params[1];
                ap2 = r.params[2];
                ap3 = r.params[3];
            } else {
                ap0 = r.strength[h.n_usr - 1]*
                    pow(h.usr_egrid[h.n_usr - 1], 3.5 + r.kl);
                ap1 = ap2 = ap3 = 0.0;
            }

            retval = StoreCTransition(db, sid,
                CFACDB_CS_PI, h.qk_mode, r.b, r.f,
                r.kl, ap0, ap1, ap2, ap3,
                &cid);
            if (retval != 0) {
                break;
            }

            sqlite3_bind_int(stmt, 1, cid);

            for (t = 0; t < h.n_usr; t++) {
                sqlite3_bind_double(stmt, 2, h.usr_egrid[t]);
                sqlite3_bind_double(stmt, 3, r.strength[t]);

                rc = sqlite3_step(stmt);
                if (rc != SQLITE_DONE) {
                    fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                    retval = -1;
                    break;
                }
                sqlite3_reset(stmt);
            }

            free(r.params); 
            free(r.strength);
        }

        free(h.tegrid);
        free(h.egrid);
        free(h.usr_egrid);
    }

    sqlite3_exec(db, "COMMIT", 0, 0, 0);

    sqlite3_finalize(stmt);

    return retval;
}

int StoreClose(sqlite3 *db, unsigned long int sid, const char *cmdline)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "UPDATE sessions" \
          " SET cmdline = ?" \
          " WHERE sid = ?";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    SQLITE3_BIND_STR(stmt, 1, cmdline);
    sqlite3_bind_int(stmt, 2, sid);

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE) {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
        retval = -1;
    }
    
    sqlite3_close(db);
    
    return retval;
}
