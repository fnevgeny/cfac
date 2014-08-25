#ifndef __CFAC_H_
#define __CFAC_H_

typedef struct _cfac_t cfac_t;

/* cfac.c */
cfac_t *
cfac_new(void);

void 
cfac_free(cfac_t *cfac);

/* nucleus.c */
int
cfac_set_atom(cfac_t *cfac, const char *s, double z, double mass, double rn);
double
cfac_get_atomic_number(const cfac_t *cfac);
double
cfac_get_atomic_mass(const cfac_t *cfac);
double
cfac_get_atomic_rn(const cfac_t *cfac);
const char *
cfac_get_atomic_symbol(const cfac_t *cfac);
double
cfac_get_atomic_effective_z(const cfac_t *cfac, double r);

/* config.c */
int
cfac_add_config(cfac_t *cfac, const char *gname, const char *cfg_str);
int
cfac_get_config_gid(const cfac_t *cfac, const char *cname);

/* structure.c */
int
cfac_get_num_levels(const cfac_t *cfac);
int
cfac_calculate_structure(cfac_t *cfac,
    int ng, const int *gids, int npg, const int *pgids, int no_ci);

/* transition.c */
typedef struct {
    unsigned int ii, fi; /* initial (upper) and final (lower) level indices */
    double rme;          /* reduced matrix element                          */
} cfac_rtrans_data_t;

typedef int
(*cfac_tr_sink_t)(const cfac_t *cfac,
    const cfac_rtrans_data_t *rtdata, void *udata);

int
crac_calculate_rtrans(cfac_t *cfac,
    unsigned nlow, unsigned *low, unsigned nup, unsigned *up,
    int mpole, int mode,
    cfac_tr_sink_t sink, void *udata);

#endif /* __CFAC_H_ */
