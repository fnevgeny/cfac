#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "cfacP.h"
#include "transition.h"
#include "dbase.h"

#include "schema.i"

#define SQLITE3_BIND_STR(stmt, id, txt) \
        sqlite3_bind_text(stmt, id, txt, -1, SQLITE_STATIC)

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
        if (truncate(fn, 0)) {
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
    }

    sql = "INSERT INTO sessions" \
          " (sid, version, fname, config)" \
          " VALUES (?, ?, '', '')";

    sqlite3_prepare_v2(*db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt, 1, *sid);
    sqlite3_bind_int(stmt, 2, 10000*VERSION + 100*SUBVERSION + SUBSUBVERSION);

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
            
            sqlite3_bind_int   (stmt,  2, r.upper);
            sqlite3_bind_int   (stmt,  3, r.lower);
            sqlite3_bind_int   (stmt,  4, h.multipole);
            sqlite3_bind_double(stmt,  5, r.strength);
            sqlite3_bind_int   (stmt,  6, h.mode);

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
    int type, int qk_mode, int ini_id, int fin_id,
    double ap0, double ap1, double a_e,
    unsigned long int *cid)
{
    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *errmsg;
    const char *sql;
    
    *cid = 0;

    sql = "INSERT INTO ctransitions" \
          " (sid, type, qk_mode, ini_id, fin_id, ap0, ap1, a_e)" \
          " VALUES (?, ?, ?, ?, ?, ?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    sqlite3_bind_int(stmt, 1, sid);
    sqlite3_bind_int(stmt, 2, type);
    sqlite3_bind_int(stmt, 3, qk_mode);
    sqlite3_bind_int(stmt, 4, ini_id);
    sqlite3_bind_int(stmt, 5, fin_id);
    
    sqlite3_bind_double(stmt, 6, ap0);
    sqlite3_bind_double(stmt, 7, ap1);
    sqlite3_bind_double(stmt, 8, a_e);
    
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
                DB_SQL_CS_CE, QK_EXACT, r.lower, r.upper,
                r.bethe, r.born[0], r.born[1],
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

            /* TODO: r.kl; full fit? */

            retval = StoreCTransition(db, sid,
                DB_SQL_CS_CI, h.qk_mode, r.b, r.f,
                r.params[0], r.params[1], 0.0,
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
    sqlite3_bind_int(stmt, 2, DB_SQL_CS_RR);
    
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
            double ap0, ap1;
            TRANSITION tr;
            int swapped;
            double dE;
            
            n = ReadRRRecord(fp, &r, swp, &h);
            if (n == 0) {
                break;
            }
            
            GetTransition(cfac, r.b, r.f, &tr, &swapped);
            dE = fabs(tr.e);

            /* TODO: full fit ? */
            
            if (h.qk_mode == QK_FIT && r.params[0]) {
                double p0, p1, p2, p3;
                p0 = r.params[0];
                p1 = r.params[1];
                p2 = r.params[2];
                p3 = r.params[3];
                ap0 = p0*pow(1.0 + p2, p1)*pow(p3/dE, 3.5 + r.kl);
            } else {
                ap0 = r.strength[h.n_usr - 1]*
                    pow(h.usr_egrid[h.n_usr - 1], 3.5 + r.kl);
            }
            ap1 = (double) r.kl;

            retval = StoreCTransition(db, sid,
                DB_SQL_CS_RR, h.qk_mode, r.b, r.f,
                ap0, ap1, 0.0,
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

int StoreClose(sqlite3 *db)
{
    sqlite3_close(db);
    
    return 0;
}
