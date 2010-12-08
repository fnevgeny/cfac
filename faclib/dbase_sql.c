#include "dbase.h"
#include "structure.h"

#define SQLITE3_BIND_STR(stmt, id, txt) \
        sqlite3_bind_text(stmt, id, txt, -1, SQLITE_STATIC)

int StoreENTable(sqlite3 *db, FILE *fp, int swp)
{
    EN_HEADER h;
    EN_RECORD r;
    int i, n;
    int p, vnl, vn, vl, g, ibase;
    double e;

    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "INSERT INTO levels" \
          " (id, name, e, g, vn, vl, p, ibase, ncomplex, sname)" \
          " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    while (1) {
        n = ReadENHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }

        for (i = 0; i < h.nlevels; i++) {
            n = ReadENRecord(fp, &r, swp);
            if (n == 0) {
                break;
            }
            e = r.energy;
	    // e -= mem_en_table[iground].energy;
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
    
            sqlite3_bind_int   (stmt,  1, r.ilev);
            SQLITE3_BIND_STR   (stmt,  2, r.name);
            sqlite3_bind_double(stmt,  3, e);
            sqlite3_bind_int   (stmt,  4, g);
            sqlite3_bind_int   (stmt,  5, vn);
            sqlite3_bind_int   (stmt,  6, vl);
            sqlite3_bind_int   (stmt,  7, p);
            if (ibase >= 0) {
                sqlite3_bind_int (stmt, 8, ibase);
            } else {
                sqlite3_bind_null(stmt, 8);
            }
            SQLITE3_BIND_STR   (stmt,  9, r.ncomplex);
            SQLITE3_BIND_STR   (stmt, 10, r.sname);

            rc = sqlite3_step(stmt);
            if (rc != SQLITE_DONE) {
                fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                retval = -1;
                break;
            }
            sqlite3_reset(stmt);
        }
    }
    
    sqlite3_finalize(stmt);

    return retval;
}


int StoreTRTable(sqlite3 *db, FILE *fp, int swp) {
    TR_HEADER h;
    TR_RECORD r;
    TR_EXTRA rx;
    int n, i;

    int retval = 0;
    int rc;
    sqlite3_stmt *stmt;
    
    char *sql;

    sql = "INSERT INTO transitions" \
          " (ini_id, fin_id, mpole, me, mode)" \
          " VALUES (?, ?, ?, ?, ?)";
    
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    while (1) {
        n = ReadTRHeader(fp, &h, swp);
        if (n == 0) {
            break;
        }

        // fprintf(f2, "GAUGE\t= %d\n", (int)h.gauge);

        for (i = 0; i < h.ntransitions; i++) {
            n = ReadTRRecord(fp, &r, &rx, swp);
            if (n == 0) {
                break;
            }
            
            sqlite3_bind_int   (stmt,  1, r.upper);
            sqlite3_bind_int   (stmt,  2, r.lower);
            sqlite3_bind_int   (stmt,  3, h.multipole);
            sqlite3_bind_double(stmt,  4, r.strength);
            sqlite3_bind_int   (stmt,  5, h.mode);

            rc = sqlite3_step(stmt);
            if (rc != SQLITE_DONE) {
                fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                retval = -1;
                break;
            }
            sqlite3_reset(stmt);
        }
    }

    sqlite3_finalize(stmt);

    return retval;
}
