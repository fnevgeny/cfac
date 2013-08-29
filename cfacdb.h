#ifndef _CFACDB_H
#define _CFACDB_H

#include <sqlite3.h>

/* TODO: include a public cfac header when it exists... */
#define DB_SQL_CS_CE    1
#define DB_SQL_CS_CI    2
#define DB_SQL_CS_PI    3

typedef struct {
    sqlite3 *db;
    
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
} cfac_db_t;

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
    double ap0, ap1;
    
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

typedef void (*cfac_db_crates_sink_t)(const cfac_db_t *cdb,
    crates_cb_data_t *cbdata, void *udata);


cfac_db_t *cdb_init(const char *fname, int nele_min, int nele_max);
void cfac_db_close(cfac_db_t *cdb);

int cfac_db_levels(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, levels_cb_data_t *cbdata, void *udata),
    void *udata);

int cfac_db_rtrans(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, rtrans_cb_data_t *cbdata, void *udata),
    void *udata);

int cfac_db_aitrans(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, aitrans_cb_data_t *cbdata, void *udata),
    void *udata);

int cfac_db_ctrans(cfac_db_t *cdb,
    void (*sink)(const cfac_db_t *cdb, ctrans_cb_data_t *cbdata, void *udata),
    void *udata);

int cfac_db_crates(cfac_db_t *cdb, double T,
    cfac_db_crates_sink_t sink,
    void *udata);

#endif /* _CFACDB_H */
