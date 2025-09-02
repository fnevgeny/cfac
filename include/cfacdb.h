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

/*! \file cfacdb.h
 *  \brief The public cFACdb API
 *
 *  \warning All physical data and parameters are expressed in the atomic units.
*/

#if defined(__cplusplus)
extern "C" {
#endif

#include <time.h>

#ifndef _CFACDB_H
#define _CFACDB_H

/* Boolean values */
/*! Boolean false. */
#define CFACDB_FALSE 0
/*! Boolean true. */
#define CFACDB_TRUE  1

/* Return values */
/*! Return value indicating success. */
#define CFACDB_SUCCESS 0
/*! Return value indicating failure. */
#define CFACDB_FAILURE 1

/* Collision process types */
/*! Collision excitation. */
#define CFACDB_CS_CE    1
/*! Collision ionization. */
#define CFACDB_CS_CI    2
/*! Photoionization. */
#define CFACDB_CS_PI    3

/*!
 * \brief Storage type for temporary views, indices, etc .
 */
typedef enum {
    CFACDB_TEMP_DEFAULT = 0,    /*!< SQLite default. */
    CFACDB_TEMP_FILE,           /*!< File.           */
    CFACDB_TEMP_MEMORY          /*!< Memory.         */
} cfacdb_temp_t;

/*!
 * \brief Calculation mode.
 */
typedef enum {
    CFACDB_MODE_UNKNOWN = 0,    /*!< Could not be determined. */
    CFACDB_MODE_DETAILED,       /*!< Detailed calculations.   */
    CFACDB_MODE_UTA,            /*!< UTA calculations.        */
    CFACDB_MODE_MIXED           /*!< Mixed mode.              */
} cfacdb_mode_t;

/*!
 * \brief The main cFACdb structure (used opaquely throughout the API).
 */
typedef struct _cfacdb_t cfacdb_t;

/*!
 * \brief Per-session statistics.
 */
typedef struct {
    unsigned long ndim;     /*!< Number of levels.                */
    unsigned long rtdim;    /*!< Number of radiative transitions. */
    unsigned long aidim;    /*!< Number of AI transitions.        */
    unsigned long cedim;    /*!< Number of CE transitions.        */
    unsigned long cidim;    /*!< Number of CI transitions.        */
    unsigned long pidim;    /*!< Number of PI transitions.        */
    cfacdb_mode_t mode;     /*!< Mode of calculations.            */
} cfacdb_stats_t;

/*!
 * \brief Session basic info.
 */
typedef struct {
    unsigned long int sid;  /*!< Session ID.                  */
    const char *sym;        /*!< Atomic symbol.               */
    unsigned int anum;      /*!< Atomic number.               */
    double mass;            /*!< Atomic mass.                 */
    unsigned int nele_min;  /*!< Minimal number of electrons. */
    unsigned int nele_max;  /*!< Maximal number of electrons. */
    unsigned int version;   /*!< cFAC version.                */
    time_t tstamp;          /*!< Time stamp.                  */
} cfacdb_sessions_data_t;

/*!
 * \brief Charge-state data.
 */
typedef struct {
    unsigned int nele;      /*!< Number of electrons.    */
    double e_gs;            /*!< Ground-state energy.    */
    unsigned long nlevels;  /*!< Total number of levels. */
} cfacdb_cstates_data_t;

/*!
 * \brief Level data.
 */
typedef struct {
    unsigned int i;     /*!< Level index.                            */

    unsigned int ifac;  /*!< cFAC level index.                       */

    double energy;      /*!< Level energy.                           */
    unsigned int nele;  /*!< Number of electrons.                    */
    unsigned int g;     /*!< Degeneracy (2J + 1 if non-UTA).         */
    unsigned int vn;    /*!< PQN of the valence electron.            */
    unsigned int vl;    /*!< Orbital number of the valence electron. */
    unsigned int p;     /*!< Parity.                                 */
    const char *name;   /*!< Name in the jj notation.                */
    const char *ncmplx; /*!< Name of the complex.                    */
    const char *sname;  /*!< Non-relativistic notation.              */
    int uta;            /*!< UTA flag.                               */
} cfacdb_levels_data_t;

/*!
 * \brief Radiative transition data.
 */
typedef struct {
    unsigned int ii;    /*!< Initial level ID.                      */
    unsigned int fi;    /*!< Final level ID (\f$E_f > E_i\f$).      */

    int mpole;          /*!< Multipole type (-1 = E1, 1 = M1, etc). */

    double de;          /*!< Transition energy.                     */

    double gf;          /*!< Symmetrized oscillator strength.       */

    double uta_de;      /*!< UTA shift.                             */
    double uta_sd;      /*!< UTA Gaussian width.                    */
} cfacdb_rtrans_data_t;

/*!
 * \brief Autoionization transition data.
 */
typedef struct {
    unsigned int ii;    /*!< Initial level ID.                */
    unsigned int fi;    /*!< Final level ID (\f$E_f > E_i\f$).*/

    double rate;        /*!< AI rate.                         */
} cfacdb_aitrans_data_t;

/*!
 * \brief Collisional transition data.
 */
typedef struct {
    unsigned int cid;   /*!< Collision transition ID.           */
    unsigned int ii;    /*!< Initial level ID.                  */
    unsigned int fi;    /*!< Final level ID (\f$E_f > E_i\f$).  */

    unsigned int type;  /*!< Process type.                      */

    double de;          /*!< Transition (threshold) energy.     */
    unsigned int kl;    /*!< Dominant L of the ionized shell.   */
    double ap0;         /*!< A fit/extrapolation parameter.     */
    double ap1;         /*!< A fit/extrapolation parameter.     */
    double ap2;         /*!< A fit/extrapolation parameter.     */
    double ap3;         /*!< A fit/extrapolation parameter.     */

    unsigned int nd;    /*!< Number of collision-strength data. */
    double *e;          /*!< Energy mesh (in unis of de).       */
    double *d;          /*!< Collision-strength data.           */
} cfacdb_ctrans_data_t;

/*!
 * \brief Collisional rate data.
 */
typedef struct {
    unsigned int ii;    /*!< Initial level ID.                  */
    unsigned int fi;    /*!< Final level ID (\f$E_f > E_i\f$).  */

    unsigned int type;  /*!< Process type.                      */

    double de;          /*!< Transition energy.                 */

    double ratec;       /*!< Rate coefficient                   */
} cfacdb_crates_data_t;

/*!
 * \brief Interpolationa/extrapolation data structure.
 */
typedef struct {
    unsigned int ndata; /*!< Number of data points.                   */
    double      *e;     /*!< Energy-grid array, length = ndata.       */
    double      *d;     /*!< Data array, length = ndata.              */

    double       ap[5]; /*!< Asymptote parameters.                    */

    double       d0;    /*!< Threshold limit.                         */
    double       let;   /*!< Low-E tangent.                           */

    int          cube;  /*!< Cubic interpolation flag.                */
    double     (*f_asymptote)(double x, const double *ap);
                        /*!< Function for evaluating high-E asymptote.*/
} cfacdb_intext_t;

/*!
 * \brief Field configuration data.
 */
typedef struct {
    unsigned int fid; /*!< Field configuration ID.            */

    double ef;        /*!< Electric field (V/cm).             */
    double bf;        /*!< Magnetic field (gauss).            */
    double angle;     /*!< Angle between EF and MF (degrees). */
} cfacdb_fields_data_t;

/*!
 * \brief State data.
 */
typedef struct {
    unsigned int i;              /*!< State index.                           */

    unsigned int ifac;           /*!< cFAC state index.                      */

    cfacdb_levels_data_t *level; /* Parent level.                            */

    double de;                   /* Energy difference from the parent level. */
    int mj;                      /* Angular momentum projection (leading).   */
} cfacdb_states_data_t;

/*!
 * \brief m-resolved radiative transition data.
 */
typedef struct {
    unsigned int ii;    /*!< Initial state ID.                      */
    unsigned int fi;    /*!< Final state ID (\f$E_f > E_i\f$).      */

    int mpole;          /*!< Multipole type (-1 = E1, 1 = M1, etc). */

    double de;          /*!< Transition energy.                     */

    int q;              /*!< Multipole component.                   */

    double gf;          /*!< Symmetrized oscillator strength.       */
} cfacdb_rtrans_m_data_t;

/*!
 * \brief Session data sink prototype
 */
typedef int (*cfacdb_sessions_sink_t)(const cfacdb_t *cdb,
    cfacdb_sessions_data_t *cbdata, void *udata);
/*!
 * \brief Charge-state data sink prototype
 */
typedef int (*cfacdb_cstates_sink_t)(const cfacdb_t *cdb,
    cfacdb_cstates_data_t *cbdata, void *udata);
/*!
 * \brief Level data sink prototype
 */
typedef int (*cfacdb_levels_sink_t)(const cfacdb_t *cdb,
    cfacdb_levels_data_t *cbdata, void *udata);
/*!
 * \brief Radiative transition data sink prototype
 */
typedef int (*cfacdb_rtrans_sink_t)(const cfacdb_t *cdb,
    cfacdb_rtrans_data_t *cbdata, void *udata);
/*!
 * \brief AI transition data sink prototype
 */
typedef int (*cfacdb_aitrans_sink_t)(const cfacdb_t *cdb,
    cfacdb_aitrans_data_t *cbdata, void *udata);
/*!
 * \brief Collision transition data sink prototype
 */
typedef int (*cfacdb_ctrans_sink_t)(const cfacdb_t *cdb,
    cfacdb_ctrans_data_t *cbdata, void *udata);
/*!
 * \brief Collision rate data sink prototype
 */
typedef int (*cfacdb_crates_sink_t)(const cfacdb_t *cdb,
    cfacdb_crates_data_t *cbdata, void *udata);
/*!
 * \brief Field data sink prototype
 */
typedef int (*cfacdb_fields_sink_t)(const cfacdb_t *cdb,
    cfacdb_fields_data_t *cbdata, void *udata);
/*!
 * \brief State data sink prototype
 */
typedef int (*cfacdb_states_sink_t)(const cfacdb_t *cdb,
    cfacdb_states_data_t *cbdata, void *udata);
/*!
 * \brief m-resolved radiative transition data sink prototype
 */
typedef int (*cfacdb_rtrans_m_sink_t)(const cfacdb_t *cdb,
    cfacdb_rtrans_m_data_t *cbdata, void *udata);

/*!
 * \brief cFACdb constructor.
 * Opens a database and allocates a new cFACdb object.
 * \param fname The SQLite database file to open.
 * \param temp_store Storage type of temporary SQL constructs.
 * \return The object allocated or NULL if failed.
 */
cfacdb_t *cfacdb_open(const char *fname, cfacdb_temp_t temp_store);
/*!
 * \brief cFACdb destructor.
 * Closes the database and deallocates all associated memory.
 * \param cdb The cFACdb object.
 */
void cfacdb_close(cfacdb_t *cdb);
/*!
 * \brief Initialize the database.
 *
 * Select a specific session in the database, optionally choosing a subset
 * of the data according to a range of charge states.
 * \param cdb The cFACdb object.
 * \param sid Session ID. Use 0 to pick the last one.
 * \param nele_min The minimal number of electrons.
 * \param nele_max The maximal number of electrons.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_init(cfacdb_t *cdb, unsigned long sid, int nele_min, int nele_max);

/*!
 * \brief Set arbitrary user-supplied data.
 * \param cdb The cFACdb object.
 * \param udata An opaque pointer to the user-supplied data.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_set_udata(cfacdb_t *cdb, void *udata);
/*!
 * \brief Get the user-supplied data.
 * \param cdb The cFACdb object.
 * \return The user-supplied data.
 */
void *cfacdb_get_udata(cfacdb_t *cdb);

/*!
 * \brief Get number of sessions.
 * \param cdb The cFACdb object.
 * \return The number of sessions in the connected database.
 */
unsigned int cfacdb_get_nsessions(const cfacdb_t *cdb);

/*!
 * \brief Get session info.
 * \param cdb The cFACdb object.
 * \param sink A user-provided function invoked for each session present.
 * \param udata An opaque pointer to arbitrary data, passed to sink.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_sessions(const cfacdb_t *cdb,
    cfacdb_sessions_sink_t sink, void *udata);

/*!
 * \brief Get basic species properties.
 * \param cdb The cFACdb object.
 * \param anum Atomic number.
 * \param mass Atomic mass.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_get_species(const cfacdb_t *cdb, unsigned int *anum, double *mass);
/*!
 * \brief Get some statistics on the data.
 * \param cdb The cFACdb object.
 * \param stats A pointer to the structure to be filled in by the function.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_get_stats(const cfacdb_t *cdb, cfacdb_stats_t *stats);

/*!
 * \brief Get charge state info.
 * \param cdb The cFACdb object.
 * \param sink A user-provided function invoked for each charge state present.
 * \param udata An opaque pointer to arbitrary data, passed to sink.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_cstates(cfacdb_t *cdb, cfacdb_cstates_sink_t sink, void *udata);
/*!
 * \brief Get level data.
 * \param cdb The cFACdb object.
 * \param sink A user-provided function invoked for each level present.
 * \param udata An opaque pointer to arbitrary data, passed to sink.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_levels(cfacdb_t *cdb, cfacdb_levels_sink_t sink, void *udata);
/*!
 * \brief Get radiative transitions.
 * \param cdb The cFACdb object.
 * \param sink A user-provided function invoked for each transition.
 * \param udata An opaque pointer to arbitrary data, passed to sink.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_rtrans(cfacdb_t *cdb, cfacdb_rtrans_sink_t sink, void *udata);
/*!
 * \brief Get autoionization data.
 * \param cdb The cFACdb object.
 * \param sink A user-provided function invoked for each AI transition.
 * \param udata An opaque pointer to arbitrary data, passed to sink.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_aitrans(cfacdb_t *cdb, cfacdb_aitrans_sink_t sink, void *udata);
/*!
 * \brief Get collision processes.
 *
 * Unless interested in raw data, consider using \ref cfacdb_crates instead.
 * \param cdb The cFACdb object.
 * \param sink A user-provided function invoked for each collision process.
 * \param udata An opaque pointer to arbitrary data, passed to sink.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_ctrans(cfacdb_t *cdb, cfacdb_ctrans_sink_t sink, void *udata);

/*!
 * \brief Get Maxwellian-integrated collision rates.
 * \param cdb The cFACdb object.
 * \param T The temperature.
 * \param sink A user-provided function invoked for each collision process.
 * \param udata An opaque pointer to arbitrary data, passed to sink.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_crates(cfacdb_t *cdb,
    double T, cfacdb_crates_sink_t sink, void *udata);

/*!
 * \brief Prepare intext structure.
 * \param cdb The cFACdb object.
 * \param cbdata collision-strength data passed to sink by \ref cfacdb_ctrans.
 * \param intext poniter to the intext structure.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_prepare_intext(const cfacdb_t *cdb,
    const cfacdb_ctrans_data_t *cbdata, cfacdb_intext_t *intext);
/*!
 * \brief Evaluate (using inter/extrapolation) a collision-strength datum.
 * \param intext The intext structure.
 * \param x Projectile energy (in units of the threshold energy).
 * \return Result of the evaluation.
 */
double cfacdb_intext(const cfacdb_intext_t *intext, double x);

/*!
 * \brief Attach a database for caching collision rates.
 * \param cdb The cFACdb object.
 * \param fname Path to the cache DB.
 * \return \ref CFACDB_SUCCESS on success or \ref CFACDB_FAILURE otherwise.
 */
int cfacdb_attach_cache(cfacdb_t *cdb, const char *fname);

/*!
 * \brief Get number of field configurations.
 * \param cdb The cFACdb object.
 * \return The number of field configurations in the connected database.
 */
unsigned int cfacdb_get_nfields(const cfacdb_t *cdb);


int cfacdb_fields(const cfacdb_t *cdb, cfacdb_fields_sink_t sink, void *udata);

int cfacdb_init_field(cfacdb_t *cdb, unsigned int fid);

int cfacdb_states(const cfacdb_t *cdb, cfacdb_states_sink_t sink, void *udata);

int cfacdb_rtrans_m(cfacdb_t *cdb, cfacdb_rtrans_m_sink_t sink, void *udata);

#endif /* _CFACDB_H */

#if defined(__cplusplus)
}
#endif
