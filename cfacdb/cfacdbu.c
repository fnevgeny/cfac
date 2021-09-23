#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "cfac.h"
#include "cfacdb.h"

typedef struct {
    const char *db_fname;
    const char *cache_fname;

    int print_info;

    long sid;

    int nele_min;
    int nele_max;

    unsigned long *sids;
    unsigned long nsid;

    double T;
} cfacdbu_t;

static void verinfo(void)
{
    printf("cfacdbu - CFACDB utility (part of cFAC-%d.%d.%d).\n\n",
        CFAC_VERSION, CFAC_SUBVERSION, CFAC_SUBSUBVERSION);

    puts("License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.");
    puts("This is free software: you are free to change and redistribute it.");
    puts("There is NO WARRANTY, to the extent permitted by law.\n");

    puts("Written by Evgeny Stambulchik.");
}

static void usage(FILE *fp, const char *progname)
{
    fprintf(fp, "Usage: %s [OPTION]... <DB file>\n", progname);
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -i, --info             print basic information about the DB\n");
    fprintf(fp, "  -s, --session ID       select session ID [none]\n");
    fprintf(fp, "      --nele-min N       min number of electrons in charge state [0]\n");
    fprintf(fp, "      --nele-max N       max number of electrons in charge state [100]\n");
    fprintf(fp, "  -c, --cache FILE       create (if needed) and attach a cache DB [none]\n");
    fprintf(fp, "  -T, --temperature T    populate the cache DB with rate coefficients,\n" \
                "                         calculated at temperature T (a.u.)\n");

    fprintf(fp, "  -V, --version          print version info and exit\n");
    fprintf(fp, "  -h, --help             display this help and exit\n");
}

static int parse_args(cfacdbu_t *u, unsigned int argc, char *const *argv)
{
    int optc;

    while (CFACDB_TRUE) {
        static struct option long_options[] = {
            {"session",          required_argument, NULL,  's'},
            {"cache",            required_argument, NULL,  'c'},
            {"temperature",      required_argument, NULL,  'T'},
            {"nele-min",         required_argument, NULL,  128},
            {"nele-max",         required_argument, NULL,  129},
            {"info",             no_argument,       NULL,  'i'},
            {"version",          no_argument,       NULL,  'V'},
            {"help",             no_argument,       NULL,  'h'},
            {NULL,               0,                 NULL,    0}
        };

        /* `getopt_long' stores the option index here. */
        int option_index = 0;

        optc = getopt_long(argc, argv,
            "is:c:T:Vh",
            long_options, &option_index);

        /* Detect the end of the options. */
        if (optc == -1) {
            break;
        }

        switch (optc) {
        case 'i':
            u->print_info = CFACDB_TRUE;
            break;
        case 'c':
            u->cache_fname = optarg;
            break;
        case 's':
            if (!strcmp(optarg, "all")) {
                u->sid = -1;
            } else {
                u->sid = atoi(optarg);
            }
            break;
        case 'T':
            u->T = atof(optarg);
            if (u->T <= 0.0) {
                fprintf(stderr, " Temperature must be positive!\n");
                return CFACDB_FAILURE;
            }
            break;
        case 128:
            u->nele_min = atoi(optarg);
            if (u->nele_min < 0) {
                fprintf(stderr, " nele-min must be non-negative!\n");
                return CFACDB_FAILURE;
            }
            break;
        case 129:
            u->nele_max = atoi(optarg);
            if (u->nele_max < 0) {
                fprintf(stderr, " nele-max must be non-negative!\n");
                return CFACDB_FAILURE;
            }
            break;
        case 'V':
            verinfo();
            exit(0);
            break;
        case 'h':
            usage(stdout, argv[0]);
            exit(0);
            break;
        case '?':
            /* `getopt_long' already printed an error message. */
            usage(stderr, argv[0]);
            return CFACDB_FAILURE;
        default:
            return CFACDB_FAILURE;
        }
    }

    if (optind == argc - 1) {
        u->db_fname = argv[optind];
        optind++;
    } else {
        usage(stderr, argv[0]);
        return CFACDB_FAILURE;
    }

    /* Print any remaining command line arguments (not options). */
    if (optind < argc - 1) {
        fprintf(stderr, "unrecognized argument(s): ");
        while (optind < argc - 1) {
            fprintf(stderr, "%s ", argv[++optind]);
        }
        fprintf(stderr, "\n");
        usage(stderr, argv[0]);
        return CFACDB_FAILURE;
    }

    if (u->nele_min > u->nele_max) {
        fprintf(stderr, "nele-min > nele-max: %d > %d!\n",
            u->nele_min, u->nele_max);
        return CFACDB_FAILURE;
    }

    return CFACDB_SUCCESS;
}

static int sessions_sink(const cfacdb_t *cdb,
    cfacdb_sessions_data_t *cbdata, void *udata)
{
    cfacdbu_t *cdu = udata;

    cdu->sids[cdu->nsid] = cbdata->sid;
    cdu->nsid++;

    if (cdu->print_info) {
        printf("Session #%lu (sid = %ld):\n",
                  cdu->nsid, cbdata->sid);
        printf("\t%s (Z = %d, mass = %.2f) nele = %d ... %d, UTA = %s \n",
                  cbdata->sym, cbdata->anum, cbdata->mass,
                  cbdata->nele_min, cbdata->nele_max,
                  cbdata->uta ? "true":"false");
    }

    if (cdu->nele_min > cbdata->nele_max) {
        fprintf(stderr,
            "\tWarning: resetting nele-min to %d!\n", cbdata->nele_max);
        cdu->nele_min = cbdata->nele_max;
    }
    if (cdu->nele_max < cbdata->nele_min) {
        fprintf(stderr,
            "\tWarning: resetting nele-max to %d!\n", cbdata->nele_min);
        cdu->nele_max = cbdata->nele_min;
    }

    return CFACDB_SUCCESS;
}

static int crates_sink(const cfacdb_t *cdb,
    cfacdb_crates_data_t *cbdata, void *udata)
{
    return CFACDB_SUCCESS;
}

int main(int argc, char *const *argv)
{
    cfacdb_t *cdb;
    cfacdbu_t cdu;
    cfacdb_stats_t stats;
    unsigned int i, nsessions;

    memset(&cdu, 0, sizeof(cdu));
    cdu.nele_min = 0;
    cdu.nele_max = 100;

    if (parse_args(&cdu, argc, argv) != CFACDB_SUCCESS) {
        exit(1);
    }

    cdb = cfacdb_open(cdu.db_fname, CFACDB_TEMP_DEFAULT);
    if (!cdb) {
        exit(1);
    }

    if (cdu.cache_fname) {
        if (cfacdb_attach_cache(cdb, cdu.cache_fname) != CFACDB_SUCCESS) {
            fprintf(stderr,
                "Failed attaching cache \"%s\" to DB \"%s\"\n",
                cdu.cache_fname, cdu.db_fname);
        }
    }

    nsessions = cfacdb_get_nsessions(cdb);
    if (cdu.print_info) {
        printf("%s: %d session%s\n",
            cdu.db_fname, nsessions, nsessions == 1 ? "":"s");
    }

    if (nsessions == 0) {
        cfacdb_close(cdb);
        exit(0);
    }

    cdu.sids = malloc(nsessions*sizeof(unsigned long));

    cfacdb_sessions(cdb, sessions_sink, &cdu);

    /* choose the latest session by default */
    if (cdu.sid == 0) {
        cdu.sid = cdu.sids[nsessions - 1];
    }

    /* override sid array if user wants a specific session */
    if (cdu.sid >= 0) {
        nsessions = 1;
        cdu.sids[0] = cdu.sid;
    }

    for (i = 0; i < nsessions; i++) {
        unsigned long sid = cdu.sids[i];

        if (cfacdb_init(cdb, sid, cdu.nele_min, cdu.nele_max)
            != CFACDB_SUCCESS) {
            fprintf(stderr, "Initialization of session ID %lu failed\n", sid);
            exit(1);
        }

        if (cdu.print_info) {
            if (cfacdb_get_stats(cdb, &stats) != CFACDB_SUCCESS) {
                fprintf(stderr,
                    "Failed getting statistics of DB \"%s\"\n", cdu.db_fname);
                cfacdb_close(cdb);
                exit(1);
            }

            printf("Stats of session ID %ld with nele = %d ... %d:\n",
                sid, cdu.nele_min, cdu.nele_max);
            printf("\tLevels: %lu, RT: %lu, AI: %lu, CE: %lu, CI: %lu, RR: %lu\n",
                stats.ndim, stats.rtdim, stats.aidim, stats.cedim, stats.cidim,
                stats.pidim);
        }

        if (cdu.cache_fname && cdu.T > 0.0) {
            cfacdb_crates(cdb, cdu.T, crates_sink, NULL);
        }
    }

    free(cdu.sids);

    cfacdb_close(cdb);

    exit(0);
}
