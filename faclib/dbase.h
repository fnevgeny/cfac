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

#ifndef _DBASE_H_
#define _DBASE_H_ 1

#include <stdio.h>
#include <sqlite3.h>

#define DB_EN 1
#define DB_TR 2
#define DB_CE 3
#define DB_RR 4
#define DB_AI 5
#define DB_CI 6
#define DB_AIM 7
#define DB_CIM 8
#define DB_ENF 9
#define DB_TRF 10
#define DB_CEF 11
#define DB_CEMF 12
#define NDB   12

#define LNCOMPLEX   32
#define LSNAME      24
#define LNAME       56

typedef struct _F_HEADER_ {
  long int tsession;
  int version;
  int sversion;
  int ssversion;
  int type;
  float atom;
  char symbol[4];
  int nblocks;
} F_HEADER;
#define SIZE_F_HEADER (sizeof(long int)+5*sizeof(int)+sizeof(float)+4)

typedef struct _EN_HEADER_ {
  long int position;
  long int length;
  int nele;
  int nlevels;
} EN_HEADER;

typedef struct _EN_RECORD_ {
  short p;
  short j;
  int ilev;
  int ibase;
  double energy;
  char ncomplex[LNCOMPLEX];
  char sname[LSNAME];
  char name[LNAME];
} EN_RECORD;
#define SIZE_EN_RECORD \
  (sizeof(short)*2+sizeof(int)*2+sizeof(double)+LNCOMPLEX+LSNAME+LNAME)

typedef struct _EN_SRECORD_ {
  int p;
  int j;
  int ibase;
  double energy;
} EN_SRECORD;

typedef struct _ENF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int nlevels;
  double efield;
  double bfield;
  double fangle;
} ENF_HEADER;

typedef struct _ENF_RECORD_ {
  int ilev;
  double energy;
  int pbasis;
} ENF_RECORD;

typedef struct _TR_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int gauge;
  int mode;
  int multipole;
} TR_HEADER;

typedef struct _TRF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int gauge;
  int mode;
  int multipole;
  double efield;
  double bfield;
  double fangle;
} TRF_HEADER;

typedef struct _TR_RECORD_ {
  int lower;
  int upper;
  float rme;
} TR_RECORD;
#define SIZE_TR_RECORD (sizeof(int)+sizeof(int)+sizeof(float))

typedef struct _TR_EXTRA_ {
  float de;
  float sdev;
} TR_EXTRA;

typedef struct _TRF_RECORD_ {
  int lower;
  int upper;
  float *strength;
} TRF_RECORD;

typedef struct _CE_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int n_tegrid;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  int nparams;
  int pw_type;
  int msub;
  float te0;
  double *tegrid;
  double *egrid;
  double *usr_egrid;
} CE_HEADER;

typedef struct _CEF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int n_tegrid;
  int n_egrid;
  float te0;
  double efield;
  double bfield;
  double fangle;
  double *tegrid;
  double *egrid;
} CEF_HEADER;

typedef struct _CE_RECORD_ {
  int lower;
  int upper;
  int nsub;
  float bethe;
  float born[2];
  float *params;
  float *strength;
} CE_RECORD;

typedef struct _CEF_RECORD_ {
  int lower;
  int upper;
  float bethe;
  float born[2];
  float *strength;
} CEF_RECORD;

typedef struct _CEMF_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int n_tegrid;
  int n_egrid;
  int n_thetagrid;
  int n_phigrid;
  float te0;
  double efield;
  double bfield;
  double fangle;
  double *tegrid;
  double *egrid;
  double *thetagrid;
  double *phigrid;
} CEMF_HEADER;

typedef struct _CEMF_RECORD_ {
  int lower;
  int upper;
  float *bethe;
  float *born;
  float *strength;
} CEMF_RECORD;

typedef struct _RR_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int qk_mode;
  int multipole;
  int n_tegrid;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  int nparams;
  double *tegrid;
  double *egrid;
  double *usr_egrid;
} RR_HEADER;

typedef struct _RR_RECORD_ {
  int b;
  int f;
  int kl;
  float *params;
  float *strength;
} RR_RECORD;

typedef struct _AI_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  float emin;
  int n_egrid;
  double *egrid;
} AI_HEADER;

typedef struct _AI_RECORD_ {
  int b;
  int f;
  float rate;
} AI_RECORD;

typedef struct _AIM_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  float emin;
  int n_egrid;
  double *egrid;
} AIM_HEADER;

typedef struct _AIM_RECORD_ {
  int b;
  int f;
  int nsub;
  float *rate;
} AIM_RECORD;

typedef struct _CI_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int qk_mode;
  int n_tegrid;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  int nparams;
  int pw_type;
  double *tegrid;
  double *egrid;
  double *usr_egrid;
} CI_HEADER;

typedef struct _CI_RECORD_ {
  int b;
  int f;
  int kl;
  float *params;
  float *strength;
} CI_RECORD;

typedef struct _CIM_HEADER_ {
  long int position;
  long int length;
  int nele;
  int ntransitions;
  int n_egrid;
  int egrid_type;
  int n_usr;
  int usr_egrid_type;
  double *egrid;
  double *usr_egrid;
} CIM_HEADER;

typedef struct _CIM_RECORD_ {
  int b;
  int f;
  int nsub;
  float *strength;
} CIM_RECORD;

/* these read functions interface with the binary data files.
 * they can be used in custom c/c++ codes to read the binary
 * files directly. to do so, copy consts.h, dbase.h, and dbase.c
 * into a working directory, and compile and link dbase.c against the
 * custom code using these functions.
 */
int ReadFHeader(FILE *f, F_HEADER *fh, int *swp);
int ReadENHeader(FILE *f, EN_HEADER *h, int swp);
int ReadENRecord(FILE *f, EN_RECORD *r, int swp);
int ReadENFHeader(FILE *f, ENF_HEADER *h, int swp);
int ReadENFRecord(FILE *f, ENF_RECORD *r, int swp);
int ReadTRHeader(FILE *f, TR_HEADER *h, int swp);
int ReadTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp);
int ReadTRFHeader(FILE *f, TRF_HEADER *h, int swp);
int ReadTRFRecord(FILE *f, TRF_RECORD *r, int swp, TRF_HEADER *h);
int ReadCEHeader(FILE *f, CE_HEADER *h, int swp);
int ReadCERecord(FILE *f, CE_RECORD *r, int swp, CE_HEADER *h);
int ReadCEFHeader(FILE *f, CEF_HEADER *h, int swp);
int ReadCEFRecord(FILE *f, CEF_RECORD *r, int swp, CEF_HEADER *h);
int ReadCEMFHeader(FILE *f, CEMF_HEADER *h, int swp);
int ReadCEMFRecord(FILE *f, CEMF_RECORD *r, int swp, CEMF_HEADER *h);
int ReadRRHeader(FILE *f, RR_HEADER *h, int swp);
int ReadRRRecord(FILE *f, RR_RECORD *r, int swp, RR_HEADER *h);
int ReadAIHeader(FILE *f, AI_HEADER *h, int swp);
int ReadAIRecord(FILE *f, AI_RECORD *r, int swp);
int ReadAIMHeader(FILE *f, AIM_HEADER *h, int swp);
int ReadAIMRecord(FILE *f, AIM_RECORD *r, int swp);
int ReadCIHeader(FILE *f, CI_HEADER *h, int swp);
int ReadCIMHeader(FILE *f, CIM_HEADER *h, int swp);
int ReadCIRecord(FILE *f, CI_RECORD *r, int swp, CI_HEADER *h);
int ReadCIMRecord(FILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h);

void CEMF2CEFHeader(CEMF_HEADER *mh, CEF_HEADER *h);
void CEMF2CEFRecord(CEMF_RECORD *mr, CEF_RECORD *r, CEMF_HEADER *mh,
                    int ith, int iph);
/* to accommadate for the possible larger statistical weight of UTA levels.
 * the eqivalent 2j value for them are stored in r.ibase of the EN_RECORD for
 * UTA. The two functions here are wrappers to determine the 2j and ibase for
 * the level, depending on whether r.j < 0.
 */
int UTAFromENRecord(const EN_RECORD *r);
int JFromENRecord(const EN_RECORD *r);
int IBaseFromENRecord(const EN_RECORD *r);

int SaveLevels(const cfac_t *cfac, const char *fn, int start, int n);

int SaveTransition(cfac_t *cfac,
    unsigned nlow, unsigned *low, unsigned nup, unsigned *up,
    const char *fn, int mpole);

/* these are the write functions, which shouldn't be of much interest.
 * unless one needs to format the external data into FAC binary format.
 */
int WriteFHeader(FILE *f, F_HEADER *fh);
int WriteENHeader(FILE *f, EN_HEADER *h);
int WriteENFHeader(FILE *f, ENF_HEADER *h);
int WriteTRHeader(FILE *f, TR_HEADER *h);
int WriteTRFHeader(FILE *f, TRF_HEADER *h);
int WriteCEHeader(FILE *f, CE_HEADER *h);
int WriteCEFHeader(FILE *f, CEF_HEADER *h);
int WriteCEMFHeader(FILE *f, CEMF_HEADER *h);
int WriteRRHeader(FILE *f, RR_HEADER *h);
int WriteAIHeader(FILE *f, AI_HEADER *h);
int WriteAIMHeader(FILE *f, AIM_HEADER *h);
int WriteCIHeader(FILE *f, CI_HEADER *h);
int WriteCIMHeader(FILE *f, CIM_HEADER *h);

int CheckEndian(F_HEADER *fh);
void SwapEndian(char *p, int size);
int SwapEndianFHeader(F_HEADER *h);
int InitDBase(void);
FILE *OpenFile(const cfac_t *cfac, const char *fn, F_HEADER *fhdr);
int CloseFile(FILE *f, F_HEADER *fhdr);
int InitFile(FILE *f, F_HEADER *fhdr, void *rhdr);
int DeinitFile(FILE *f, F_HEADER *fhdr);
int PrintTable(const cfac_t *cfac, char *ifn, char *ofn, int v);
int MemENTable(const cfac_t *cfac, char *fn);
int WriteENRecord(FILE *f, EN_RECORD *r);
int WriteENFRecord(FILE *f, ENF_RECORD *r);
int PrintENTable(FILE *f1, FILE *f2, int v, int swp);
int PrintENFTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianENHeader(EN_HEADER *h);
int SwapEndianENRecord(EN_RECORD *r);
int SwapEndianENFHeader(ENF_HEADER *h);
int SwapEndianENFRecord(ENF_RECORD *r);
int WriteTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx);
int PrintTRTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianTRHeader(TR_HEADER *h);
int SwapEndianTRRecord(TR_RECORD *r);
int WriteTRFRecord(FILE *f, TRF_RECORD *r);
int PrintTRFTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianTRFHeader(TRF_HEADER *h);
int SwapEndianTRFRecord(TRF_RECORD *r);
int WriteCERecord(FILE *f, CE_RECORD *r);
int PrintCETable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianCEHeader(CE_HEADER *h);
int SwapEndianCERecord(CE_RECORD *r);
int WriteCEFRecord(FILE *f, CEF_RECORD *r);
int PrintCEFTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianCEFHeader(CEF_HEADER *h);
int SwapEndianCEFRecord(CEF_RECORD *r);
int WriteCEMFRecord(FILE *f, CEMF_RECORD *r);
int PrintCEMFTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianCEMFHeader(CEMF_HEADER *h);
int SwapEndianCEMFRecord(CEMF_RECORD *r);
int WriteRRRecord(FILE *f, RR_RECORD *r);
int PrintRRTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianRRHeader(RR_HEADER *h);
int SwapEndianRRRecord(RR_RECORD *r);
int WriteAIRecord(FILE *f, AI_RECORD *r);
int PrintAITable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianAIHeader(AI_HEADER *h);
int SwapEndianAIRecord(AI_RECORD *r);
int WriteAIMRecord(FILE *f, AIM_RECORD *r);
int PrintAIMTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianAIMHeader(AIM_HEADER *h);
int SwapEndianAIMRecord(AIM_RECORD *r);
int WriteCIRecord(FILE *f, CI_RECORD *r);
int PrintCITable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianCIHeader(CI_HEADER *h);
int SwapEndianCIRecord(CI_RECORD *r);
int WriteCIMRecord(FILE *f, CIM_RECORD *r);
int PrintCIMTable(FILE *f1, FILE *f2, int v, int swp);
int SwapEndianCIMHeader(CIM_HEADER *h);
int SwapEndianCIMRecord(CIM_RECORD *r);
int AppendTable(const cfac_t *cfac, char *fn);
int JoinTable(const cfac_t *cfac, char *fn1, char *fn2, char *fn);
int FindLevelByName(char *fn, int nele, char *nc, char *cnr, char *cr);
int ISearch(int i, int n, int *ia);

int StoreInit(const cfac_t *cfac,
    const char *fn, int reset, sqlite3 **db, unsigned long *sid);
int StoreTable(const cfac_t *cfac,
    sqlite3 *db, unsigned long int sid, const char *ifn);
int StoreENTable(const cfac_t *cfac, sqlite3 *db,
    unsigned long int sid, FILE *fp, int swp);
int StoreTRTable(const cfac_t *cfac, sqlite3 *db,
    unsigned long int sid, FILE *fp, int swp);
int StoreCETable(const cfac_t *cfac, sqlite3 *db,
    unsigned long int sid, FILE *fp, int swp);
int StoreRRTable(const cfac_t *cfac,
    sqlite3 *db, unsigned long int sid, FILE *fp, int swp);
int StoreAITable(const cfac_t *cfac, sqlite3 *db,
    unsigned long int sid, FILE *fp, int swp);
int StoreCITable(const cfac_t *cfac, sqlite3 *db,
    unsigned long int sid, FILE *fp, int swp);
int StoreClose(const cfac_t *cfac, sqlite3 *db,
    unsigned long int sid, const char *cmdline);

#endif
