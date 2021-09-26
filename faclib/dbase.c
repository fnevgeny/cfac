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
#include <string.h>
#include <time.h>
#include <math.h>

#include "cfacP.h"
#include "consts.h"
#include "dbase.h"
#include "structure.h"
#include "transition.h"

static int version_read[NDB];
static F_HEADER fheader[NDB];
static EN_HEADER en_header;
static ENF_HEADER enf_header;
static TR_HEADER tr_header;
static TRF_HEADER trf_header;
static CE_HEADER ce_header;
static CEF_HEADER cef_header;
static CEMF_HEADER cemf_header;
static RR_HEADER rr_header;
static AI_HEADER ai_header;
static AIM_HEADER aim_header;
static CI_HEADER ci_header;
static CIM_HEADER cim_header;

static EN_SRECORD *mem_en_table = NULL;
static int mem_en_table_size = 0;
static EN_SRECORD *mem_enf_table = NULL;
static int mem_enf_table_size = 0;
static int iground;
static int iuta = 0;

#define _WSF0(sv, f) {                                  \
    n = fwrite(&(sv), sizeof(sv), 1, f);                \
    if (n != 1) return 0;                               \
    m += sizeof(sv);                                    \
  }while(0)
#define _RSF0(sv, f) {                                  \
    n = fread(&(sv), sizeof(sv), 1, f);                 \
    if (n != 1) return 0;                               \
    m += sizeof(sv);                                    \
  }while(0)
#define _WSF1(sv, s, k, f) {                            \
    n = fwrite(sv, s, k, f);                            \
    if ((n) != (k)) return 0;                           \
    m += (s)*(k);                                       \
  }while(0)
#define _RSF1(sv, s, k, f) {                            \
    n = fread(sv, s, k, f);                             \
    if ((n) != (k)) return 0;                           \
    m += (s)*(k);                                       \
  }while(0)
#define WSF0(sv) _WSF0(sv, f)
#define WSF1(sv, s, k) _WSF1(sv, s, k, f)
#define RSF0(sv) _RSF0(sv, f)
#define RSF1(sv, s, k) _RSF1(sv, s, k, f)

void SetVersionRead(int t, int v) {
  version_read[t-1] = v;
}

int CheckEndian(F_HEADER *fh) {
  unsigned short t = 0x01;
  char *p;

  if (fh) {
    if (fh->version > 0 || fh->sversion >= 7) {
      return (int) (fh->symbol[3]);
    }
  }

  p = (char *) &t;
  p += sizeof(unsigned short)-1;
  if ((unsigned short) (*p) == 1) return 1;
  else return 0;
}

void SwapEndian(char *p, int size) {
  int t1, t2;
  char tmp;

  t1 = 0;
  t2 = size-1;
  while (t2 > t1) {
    tmp = p[t1];
    p[t1] = p[t2];
    p[t2] = tmp;
    t1++;
    t2--;
  }
}

int SwapEndianFHeader(F_HEADER *h) {
  SwapEndian((char *) &(h->tsession), sizeof(long int));
  SwapEndian((char *) &(h->version), sizeof(int));
  SwapEndian((char *) &(h->sversion), sizeof(int));
  SwapEndian((char *) &(h->ssversion), sizeof(int));
  SwapEndian((char *) &(h->type), sizeof(int));
  SwapEndian((char *) &(h->atom), sizeof(float));
  SwapEndian((char *) &(h->nblocks), sizeof(int));
  return 0;
}

int SwapEndianENHeader(EN_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->nlevels), sizeof(int));
  return 0;
}

int SwapEndianENFHeader(ENF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->nlevels), sizeof(int));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
  return 0;
}

int SwapEndianENRecord(EN_RECORD *r) {
  SwapEndian((char *) &(r->p), sizeof(short));
  SwapEndian((char *) &(r->j), sizeof(short));
  SwapEndian((char *) &(r->ilev), sizeof(int));
  SwapEndian((char *) &(r->ibase), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(double));
  return 0;
}

int SwapEndianENFRecord(ENF_RECORD *r) {
  SwapEndian((char *) &(r->ilev), sizeof(int));
  SwapEndian((char *) &(r->energy), sizeof(double));
  SwapEndian((char *) &(r->pbasis), sizeof(int));
  return 0;
}

int SwapEndianTRHeader(TR_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->gauge), sizeof(int));
  SwapEndian((char *) &(h->mode), sizeof(int));
  SwapEndian((char *) &(h->multipole), sizeof(int));
  return 0;
}

int SwapEndianTRRecord(TR_RECORD *r) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->rme), sizeof(float));
  return 0;
}

int SwapEndianTRFHeader(TRF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->gauge), sizeof(int));
  SwapEndian((char *) &(h->mode), sizeof(int));
  SwapEndian((char *) &(h->multipole), sizeof(int));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
  return 0;
}

int SwapEndianTRFRecord(TRF_RECORD *r) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  return 0;
}

int SwapEndianCEHeader(CE_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  SwapEndian((char *) &(h->nparams), sizeof(int));
  SwapEndian((char *) &(h->pw_type), sizeof(int));
  SwapEndian((char *) &(h->msub), sizeof(int));
  SwapEndian((char *) &(h->te0), sizeof(float));
  return 0;
}

int SwapEndianCERecord(CE_RECORD *r) {
  int m;

  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->nsub), sizeof(int));
  SwapEndian((char *) &(r->bethe), sizeof(float));
  for (m = 0; m < 2; m++) {
    SwapEndian((char *) &(r->born[m]), sizeof(float));
  }
  return 0;
}

int SwapEndianCEFHeader(CEF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->te0), sizeof(float));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
  return 0;
}

int SwapEndianCEFRecord(CEF_RECORD *r) {
  int m;

  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  SwapEndian((char *) &(r->bethe), sizeof(float));
  for (m = 0; m < 2; m++) {
    SwapEndian((char *) &(r->born[m]), sizeof(float));
  }
  return 0;
}

int SwapEndianCEMFHeader(CEMF_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->n_thetagrid), sizeof(int));
  SwapEndian((char *) &(h->n_phigrid), sizeof(int));
  SwapEndian((char *) &(h->te0), sizeof(float));
  SwapEndian((char *) &(h->efield), sizeof(double));
  SwapEndian((char *) &(h->bfield), sizeof(double));
  SwapEndian((char *) &(h->fangle), sizeof(double));
  return 0;
}

int SwapEndianCEMFRecord(CEMF_RECORD *r) {
  SwapEndian((char *) &(r->lower), sizeof(int));
  SwapEndian((char *) &(r->upper), sizeof(int));
  return 0;
}

int SwapEndianRRHeader(RR_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->qk_mode), sizeof(int));
  SwapEndian((char *) &(h->multipole), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  SwapEndian((char *) &(h->nparams), sizeof(int));
  return 0;
}

int SwapEndianRRRecord(RR_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->kl), sizeof(int));
  return 0;
}

int SwapEndianAIHeader(AI_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->emin), sizeof(float));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  return 0;
}

int SwapEndianAIMHeader(AIM_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->emin), sizeof(float));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  return 0;
}

int SwapEndianAIRecord(AI_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->rate), sizeof(float));
  return 0;
}

int SwapEndianAIMRecord(AIM_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->nsub), sizeof(int));
  return 0;
}

int SwapEndianCIHeader(CI_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->qk_mode), sizeof(int));
  SwapEndian((char *) &(h->n_tegrid), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  SwapEndian((char *) &(h->nparams), sizeof(int));
  SwapEndian((char *) &(h->pw_type), sizeof(int));
  return 0;
}

int SwapEndianCIRecord(CI_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->kl), sizeof(int));
  return 0;
}

int SwapEndianCIMHeader(CIM_HEADER *h) {
  SwapEndian((char *) &(h->position), sizeof(long int));
  SwapEndian((char *) &(h->length), sizeof(long int));
  SwapEndian((char *) &(h->nele), sizeof(int));
  SwapEndian((char *) &(h->ntransitions), sizeof(int));
  SwapEndian((char *) &(h->n_egrid), sizeof(int));
  SwapEndian((char *) &(h->egrid_type), sizeof(int));
  SwapEndian((char *) &(h->n_usr), sizeof(int));
  SwapEndian((char *) &(h->usr_egrid_type), sizeof(int));
  return 0;
}

int SwapEndianCIMRecord(CIM_RECORD *r) {
  SwapEndian((char *) &(r->b), sizeof(int));
  SwapEndian((char *) &(r->f), sizeof(int));
  SwapEndian((char *) &(r->nsub), sizeof(int));
  return 0;
}

void CEMF2CEFHeader(CEMF_HEADER *mh, CEF_HEADER *h) {
  h->position = mh->position;
  h->length = mh->length;
  h->nele = mh->nele;
  h->ntransitions = mh->ntransitions;
  h->n_tegrid = mh->n_tegrid;
  h->n_egrid = mh->n_egrid;
  h->te0 = mh->te0;
  h->efield = mh->efield;
  h->bfield = mh->bfield;
  h->fangle = mh->fangle;
  h->tegrid = mh->tegrid;
  h->egrid = mh->egrid;
}

void CEMF2CEFRecord(CEMF_RECORD *mr, CEF_RECORD *r, CEMF_HEADER *mh,
                    int ith, int iph) {
  int k;

  r->lower = mr->lower;
  r->upper = mr->upper;
  k = ith*mh->n_phigrid + iph;
  r->bethe = mr->bethe[k];
  r->born[0] = mr->born[k];
  r->born[1] = mr->born[mh->n_phigrid*mh->n_thetagrid];
  r->strength = &(mr->strength[k*mh->n_egrid]);
}

static void FreeMemENTable(void) {
  if (mem_en_table) {
    free(mem_en_table);
    mem_en_table = NULL;
    mem_en_table_size = 0;
  }
}

int InitDBase(void) {
  int i;
  for (i = 0; i < NDB; i++) {
    fheader[i].tsession = (long int) time(0);
    fheader[i].version = CFAC_VERSION;
    fheader[i].sversion = CFAC_SUBVERSION;
    fheader[i].ssversion = CFAC_SUBSUBVERSION;
    fheader[i].symbol[2] = '\0';
    fheader[i].symbol[3] = (char) (CheckEndian(NULL));
    fheader[i].type = 0;
    fheader[i].atom = 0;
    fheader[i].nblocks = 0;
  }

  FreeMemENTable();

  mem_enf_table = NULL;
  mem_enf_table_size = 0;

  iground = 0;

  return 0;
}

int WriteFHeader(FILE *f, F_HEADER *fh) {
  int n, m = 0;

  WSF0(fh->tsession);
  WSF0(fh->version);
  WSF0(fh->sversion);
  WSF0(fh->ssversion);
  WSF0(fh->type);
  WSF0(fh->atom);
  WSF1(fh->symbol, sizeof(char), 4);
  WSF0(fh->nblocks);

  return m;
}

int ReadFHeader(FILE *f, F_HEADER *fh, int *swp) {
  int n, m = 0;

  RSF0(fh->tsession);
  RSF0(fh->version);
  RSF0(fh->sversion);
  RSF0(fh->ssversion);
  RSF0(fh->type);
  RSF0(fh->atom);
  RSF1(fh->symbol, sizeof(char), 4);
  RSF0(fh->nblocks);

  *swp = 0;

  if (CheckEndian(fh) != (int) (fheader[0].symbol[3])) {
    *swp = 1;
    SwapEndianFHeader(fh);
  }

  SetVersionRead(fh->type, fh->version*100+fh->sversion*10+fh->ssversion);

  return m;
}

int WriteENHeader(FILE *f, EN_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->nlevels);

  return m;
}

int WriteENFHeader(FILE *f, ENF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->nlevels);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);

  return m;
}

int WriteTRHeader(FILE *f, TR_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->gauge);
  WSF0(h->mode);
  WSF0(h->multipole);

  return m;
}

int WriteTRFHeader(FILE *f, TRF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->gauge);
  WSF0(h->mode);
  WSF0(h->multipole);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);

  return m;
}

int WriteCEHeader(FILE *f, CE_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF0(h->nparams);
  WSF0(h->pw_type);
  WSF0(h->msub);
  WSF0(h->te0);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);

  return m;
}

int WriteCEFHeader(FILE *f, CEF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->te0);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);

  return m;
}

int WriteCEMFHeader(FILE *f, CEMF_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->n_thetagrid);
  WSF0(h->n_phigrid);
  WSF0(h->te0);
  WSF0(h->efield);
  WSF0(h->bfield);
  WSF0(h->fangle);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->thetagrid, sizeof(double), h->n_thetagrid);
  WSF1(h->phigrid, sizeof(double), h->n_phigrid);

  return m;
}

int WriteRRHeader(FILE *f, RR_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->qk_mode);
  WSF0(h->multipole);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF0(h->nparams);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);

  return m;
}

int WriteAIHeader(FILE *f, AI_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->emin);
  WSF0(h->n_egrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);

  return m;
}

int WriteAIMHeader(FILE *f, AIM_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->emin);
  WSF0(h->n_egrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);

  return m;
}

int WriteCIHeader(FILE *f, CI_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->qk_mode);
  WSF0(h->n_tegrid);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF0(h->nparams);
  WSF0(h->pw_type);
  WSF1(h->tegrid, sizeof(double), h->n_tegrid);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);

  return m;
}

int WriteCIMHeader(FILE *f, CIM_HEADER *h) {
  int n, m = 0;

  WSF0(h->position);
  WSF0(h->length);
  WSF0(h->nele);
  WSF0(h->ntransitions);
  WSF0(h->n_egrid);
  WSF0(h->egrid_type);
  WSF0(h->n_usr);
  WSF0(h->usr_egrid_type);
  WSF1(h->egrid, sizeof(double), h->n_egrid);
  WSF1(h->usr_egrid, sizeof(double), h->n_usr);

  return m;
}

int WriteENRecord(FILE *f, EN_RECORD *r) {
  int n, m = 0, len;

  if (en_header.length == 0) {
    fheader[DB_EN-1].nblocks++;
    n = WriteENHeader(f, &en_header);
  }

  WSF0(r->p);
  WSF0(r->j);
  WSF0(r->ilev);
  WSF0(r->ibase);
  WSF0(r->energy);

  /* make sure the strings are zero-padded */
  len = strlen(r->ncomplex);
  memset(r->ncomplex + len, 0, LNCOMPLEX - len);
  len = strlen(r->sname);
  memset(r->sname + len, 0, LSNAME - len);
  len = strlen(r->name);
  memset(r->name + len, 0, LNAME - len);

  WSF1(r->ncomplex, sizeof(char), LNCOMPLEX);
  WSF1(r->sname, sizeof(char), LSNAME);
  WSF1(r->name, sizeof(char), LNAME);

  en_header.nlevels += 1;
  en_header.length += m;

  return m;
}

int WriteENFRecord(FILE *f, ENF_RECORD *r) {
  int n, m = 0;

  if (enf_header.length == 0) {
    fheader[DB_ENF-1].nblocks++;
    n = WriteENFHeader(f, &enf_header);
  }

  WSF0(r->ilev);
  WSF0(r->energy);
  WSF0(r->pbasis);

  enf_header.nlevels += 1;
  enf_header.length += m;

  return m;
}

int WriteTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx) {
  int n, m = 0;

  if (tr_header.length == 0) {
    fheader[DB_TR-1].nblocks++;
    n = WriteTRHeader(f, &tr_header);
  }

  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->rme);

  WSF0(rx->de);
  WSF0(rx->sdev);

  tr_header.ntransitions += 1;
  tr_header.length += m;

  return m;
}

int WriteTRFRecord(FILE *f, TRF_RECORD *r) {
  int n, m = 0;

  if (trf_header.length == 0) {
    fheader[DB_TRF-1].nblocks++;
    n = WriteTRFHeader(f, &trf_header);
  }

  WSF0(r->lower);
  WSF0(r->upper);
  WSF1(r->strength, sizeof(float), 2*abs(trf_header.multipole)+1);

  trf_header.ntransitions += 1;
  trf_header.length += m;

  return m;
}

int WriteCERecord(FILE *f, CE_RECORD *r) {
  int n;
  int m0, m = 0;

  if (ce_header.length == 0) {
    fheader[DB_CE-1].nblocks++;
    n = WriteCEHeader(f, &ce_header);
  }

  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->nsub);
  WSF0(r->bethe);
  WSF1(r->born, sizeof(float), 2);

  if (ce_header.msub) {
    m0 = r->nsub;
  } else m0 = 0;
  if (m0) {
    WSF1(r->params, sizeof(float), m0);
  }
  m0 = ce_header.n_usr * r->nsub;
  WSF1(r->strength, sizeof(float), m0);

  ce_header.ntransitions += 1;
  ce_header.length += m;

  return m;
}

int WriteCEFRecord(FILE *f, CEF_RECORD *r) {
  int n;
  int m0, m = 0;

  if (cef_header.length == 0) {
    fheader[DB_CEF-1].nblocks++;
    n = WriteCEFHeader(f, &cef_header);
  }

  WSF0(r->lower);
  WSF0(r->upper);
  WSF0(r->bethe);
  WSF1(r->born, sizeof(float), 2);

  m0 = cef_header.n_egrid;
  WSF1(r->strength, sizeof(float), m0);

  cef_header.ntransitions += 1;
  cef_header.length += m;

  return m;
}

int WriteCEMFRecord(FILE *f, CEMF_RECORD *r) {
  int n;
  int m0, m = 0;

  if (cemf_header.length == 0) {
    fheader[DB_CEMF-1].nblocks++;
    n = WriteCEMFHeader(f, &cemf_header);
  }

  WSF0(r->lower);
  WSF0(r->upper);

  m0 = cemf_header.n_thetagrid * cemf_header.n_phigrid;
  WSF1(r->bethe, sizeof(float), m0);
  WSF1(r->born, sizeof(float), m0+1);

  m0 = cemf_header.n_egrid * m0;
  WSF1(r->strength, sizeof(float), m0);

  cemf_header.ntransitions += 1;
  cemf_header.length += m;

  return m;
}

int WriteRRRecord(FILE *f, RR_RECORD *r) {
  int n;
  int m = 0, m0;

  if (rr_header.length == 0) {
    fheader[DB_RR-1].nblocks++;
    n = WriteRRHeader(f, &rr_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->kl);

  if (rr_header.qk_mode == QK_FIT) {
    m0 = rr_header.nparams;
    WSF1(r->params, sizeof(float), m0);
  }
  m0 = rr_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);

  rr_header.ntransitions += 1;
  rr_header.length += m;

  return m;
}

int WriteAIRecord(FILE *f, AI_RECORD *r) {
  int n, m = 0;

  if (ai_header.length == 0) {
    fheader[DB_AI-1].nblocks++;
    WriteAIHeader(f, &ai_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->rate);

  ai_header.ntransitions += 1;
  ai_header.length += m;

  return m;
}

int WriteAIMRecord(FILE *f, AIM_RECORD *r) {
  int n, m = 0;

  if (aim_header.length == 0) {
    fheader[DB_AIM-1].nblocks++;
    WriteAIMHeader(f, &aim_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->nsub);
  WSF1(r->rate, sizeof(float), r->nsub);

  aim_header.ntransitions += 1;
  aim_header.length += m;

  return m;
}

int WriteCIRecord(FILE *f, CI_RECORD *r) {
  int n;
  int m = 0, m0;

  if (ci_header.length == 0) {
    fheader[DB_CI-1].nblocks++;
    WriteCIHeader(f, &ci_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->kl);
  m0 = ci_header.nparams;
  WSF1(r->params, sizeof(float), m0);
  m0 = ci_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);

  ci_header.ntransitions += 1;
  ci_header.length += m;

  return m;
}

int WriteCIMRecord(FILE *f, CIM_RECORD *r) {
  int n;
  int m = 0, m0;

  if (cim_header.length == 0) {
    fheader[DB_CIM-1].nblocks++;
    WriteCIMHeader(f, &cim_header);
  }

  WSF0(r->b);
  WSF0(r->f);
  WSF0(r->nsub);
  m0 = r->nsub*cim_header.n_usr;
  WSF1(r->strength, sizeof(float), m0);

  cim_header.ntransitions += 1;
  cim_header.length += m;

  return m;
}

int ReadENHeader(FILE *f, EN_HEADER *h, int swp) {
  int n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->nlevels);

  if (swp) SwapEndianENHeader(h);

  return m;
}

int ReadENFHeader(FILE *f, ENF_HEADER *h, int swp) {
  int n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->nlevels);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);

  if (swp) SwapEndianENFHeader(h);

  return m;
}

int ReadENRecord(FILE *f, EN_RECORD *r, int swp) {
  int n, m = 0;

  RSF0(r->p);
  RSF0(r->j);
  RSF0(r->ilev);
  RSF0(r->ibase);
  RSF0(r->energy);
  RSF1(r->ncomplex, sizeof(char), LNCOMPLEX);
  RSF1(r->sname, sizeof(char), LSNAME);
  RSF1(r->name, sizeof(char), LNAME);

  if (swp) SwapEndianENRecord(r);

  return m;
}

int ReadENFRecord(FILE *f, ENF_RECORD *r, int swp) {
  int n, m = 0;

  RSF0(r->ilev);
  RSF0(r->energy);
  RSF0(r->pbasis);

  if (swp) SwapEndianENFRecord(r);

  return m;
}

int ReadTRHeader(FILE *f, TR_HEADER *h, int swp) {
  int n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->gauge);
  RSF0(h->mode);
  RSF0(h->multipole);

  if (swp) SwapEndianTRHeader(h);

  if (h->length/h->ntransitions > SIZE_TR_RECORD) {
    iuta = 1;
  } else {
    iuta = 0;
  }

  return m;
}

int ReadTRFHeader(FILE *f, TRF_HEADER *h, int swp) {
  int n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->gauge);
  RSF0(h->mode);
  RSF0(h->multipole);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);
  if (swp) SwapEndianTRFHeader(h);

  return m;
}

int ReadTRRecord(FILE *f, TR_RECORD *r, TR_EXTRA *rx, int swp) {
  int n, m = 0;

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->rme);

  if (iuta) {
    RSF0(rx->de);
    RSF0(rx->sdev);
  }

  if (swp) SwapEndianTRRecord(r);

  return m;
}

int ReadTRFRecord(FILE *f, TRF_RECORD *r, int swp, TRF_HEADER *h) {
  int n, m = 0, i, nq;

  RSF0(r->lower);
  RSF0(r->upper);
  nq = 2*abs(h->multipole) + 1;
  r->strength = (float *) malloc(sizeof(float)*nq);
  RSF1(r->strength, sizeof(float), nq);

  if (swp) {
    SwapEndianTRFRecord(r);
    for (i = 0; i < nq; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadCEHeader(FILE *f, CE_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);
  RSF0(h->nparams);
  RSF0(h->pw_type);
  RSF0(h->msub);
  RSF0(h->te0);

  if (swp) SwapEndianCEHeader(h);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCEFHeader(FILE *f, CEF_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->te0);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);

  if (swp) SwapEndianCEFHeader(h);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCEMFHeader(FILE *f, CEMF_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->n_thetagrid);
  RSF0(h->n_phigrid);
  RSF0(h->te0);
  RSF0(h->efield);
  RSF0(h->bfield);
  RSF0(h->fangle);

  if (swp) SwapEndianCEMFHeader(h);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->thetagrid = (double *) malloc(sizeof(double)*h->n_thetagrid);
  RSF1(h->thetagrid, sizeof(double), h->n_thetagrid);

  h->phigrid = (double *) malloc(sizeof(double)*h->n_phigrid);
  RSF1(h->phigrid, sizeof(double), h->n_phigrid);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_thetagrid; i++) {
      SwapEndian((char *) &(h->thetagrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_phigrid; i++) {
      SwapEndian((char *) &(h->phigrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCERecord(FILE *f, CE_RECORD *r, int swp, CE_HEADER *h) {
  int i, n, m = 0, m0;

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->nsub);
  RSF0(r->bethe);
  RSF1(r->born, sizeof(float), 2);

  if (swp) SwapEndianCERecord(r);

  if (h->msub) {
    m0 = r->nsub;
  } else m0 = 0;
  r->params = NULL;
  if (m0) {
    r->params = (float *) malloc(sizeof(float)*m0);
    RSF1(r->params, sizeof(float), m0);
    if (swp) {
      for (i = 0; i < m0; i++) {
        SwapEndian((char *) &(r->params[i]), sizeof(float));
      }
    }
  }

  m0 = h->n_usr * r->nsub;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadCEFRecord(FILE *f, CEF_RECORD *r, int swp, CEF_HEADER *h) {
  int i, n, m = 0, m0;

  RSF0(r->lower);
  RSF0(r->upper);
  RSF0(r->bethe);
  RSF1(r->born, sizeof(float), 2);

  if (swp) SwapEndianCEFRecord(r);

  m0 = h->n_egrid;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadCEMFRecord(FILE *f, CEMF_RECORD *r, int swp, CEMF_HEADER *h) {
  int i, n, m = 0, m0;

  RSF0(r->lower);
  RSF0(r->upper);
  if (swp) SwapEndianCEMFRecord(r);

  m0 = h->n_thetagrid * h->n_phigrid;
  r->bethe = (float *) malloc(sizeof(float)*m0);
  RSF1(r->bethe, sizeof(float), m0);
  r->born = (float *) malloc(sizeof(float)*(m0+1));
  RSF1(r->born, sizeof(float), m0+1);

  m0 = h->n_egrid*m0;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
    m0 = h->n_thetagrid * h->n_phigrid;
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->bethe[i]), sizeof(float));
    }
    for (i = 0; i <= m0; i++) {
      SwapEndian((char *) &(r->born[i]), sizeof(float));
    }
  }

  return m;
}

int ReadRRHeader(FILE *f, RR_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->qk_mode);
  RSF0(h->multipole);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);
  RSF0(h->nparams);

  if (swp) SwapEndianRRHeader(h);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadRRRecord(FILE *f, RR_RECORD *r, int swp, RR_HEADER *h) {
  int i, n, m = 0, m0;

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->kl);

  if (swp) SwapEndianRRRecord(r);

  if (h->qk_mode == QK_FIT) {
    m0 = h->nparams;
    r->params = (float *) malloc(sizeof(float)*m0);
    RSF1(r->params, sizeof(float), m0);
    if (swp) {
      for (i = 0; i < m0; i++) {
        SwapEndian((char *) &(r->params[i]), sizeof(float));
      }
    }
  }
  m0 = h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadAIHeader(FILE *f, AI_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->emin);
  RSF0(h->n_egrid);

  if (swp) SwapEndianAIHeader(h);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadAIMHeader(FILE *f, AIM_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->emin);
  RSF0(h->n_egrid);

  if (swp) SwapEndianAIMHeader(h);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadAIRecord(FILE *f, AI_RECORD *r, int swp) {
  int n, m = 0;

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->rate);

  if (swp) SwapEndianAIRecord(r);

  return m;
}

int ReadAIMRecord(FILE *f, AIM_RECORD *r, int swp) {
  int n, i, m = 0;

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->nsub);
  if (swp) {
    SwapEndianAIMRecord(r);
  }

  r->rate = (float *) malloc(sizeof(float)*r->nsub);
  RSF1(r->rate, sizeof(float), r->nsub);
  if (swp) {
    for (i = 0; i < r->nsub; i++) {
      SwapEndian((char *) &(r->rate[i]), sizeof(float));
    }
  }
  return m;
}

int ReadCIHeader(FILE *f, CI_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->qk_mode);
  RSF0(h->n_tegrid);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);
  RSF0(h->nparams);
  RSF0(h->pw_type);

  if (swp) SwapEndianCIHeader(h);

  h->tegrid = (double *) malloc(sizeof(double)*h->n_tegrid);
  RSF1(h->tegrid, sizeof(double), h->n_tegrid);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);

  if (swp) {
    for (i = 0; i < h->n_tegrid; i++) {
      SwapEndian((char *) &(h->tegrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCIRecord(FILE *f, CI_RECORD *r, int swp, CI_HEADER *h) {
  int i, n, m = 0, m0;

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->kl);

  if (swp) SwapEndianCIRecord(r);

  m0 = h->nparams;
  r->params = (float *) malloc(sizeof(float)*m0);
  RSF1(r->params, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->params[i]), sizeof(float));
    }
  }

  m0 = h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

int ReadCIMHeader(FILE *f, CIM_HEADER *h, int swp) {
  int i, n, m = 0;

  RSF0(h->position);
  RSF0(h->length);
  RSF0(h->nele);
  RSF0(h->ntransitions);
  RSF0(h->n_egrid);
  RSF0(h->egrid_type);
  RSF0(h->n_usr);
  RSF0(h->usr_egrid_type);

  if (swp) SwapEndianCIMHeader(h);

  h->egrid = (double *) malloc(sizeof(double)*h->n_egrid);
  RSF1(h->egrid, sizeof(double), h->n_egrid);

  h->usr_egrid = (double *) malloc(sizeof(double)*h->n_usr);
  RSF1(h->usr_egrid, sizeof(double), h->n_usr);

  if (swp) {
    for (i = 0; i < h->n_egrid; i++) {
      SwapEndian((char *) &(h->egrid[i]), sizeof(double));
    }
    for (i = 0; i < h->n_usr; i++) {
      SwapEndian((char *) &(h->usr_egrid[i]), sizeof(double));
    }
  }

  return m;
}

int ReadCIMRecord(FILE *f, CIM_RECORD *r, int swp, CIM_HEADER *h) {
  int i, n, m = 0, m0;

  RSF0(r->b);
  RSF0(r->f);
  RSF0(r->nsub);

  if (swp) SwapEndianCIMRecord(r);

  m0 = r->nsub*h->n_usr;
  r->strength = (float *) malloc(sizeof(float)*m0);
  RSF1(r->strength, sizeof(float), m0);
  if (swp) {
    for (i = 0; i < m0; i++) {
      SwapEndian((char *) &(r->strength[i]), sizeof(float));
    }
  }

  return m;
}

FILE *OpenFile(const cfac_t *cfac, const char *fn, F_HEADER *fhdr) {
  int ihdr;
  FILE *f;

  ihdr = fhdr->type - 1;

  f = fopen(fn, "r+b");
  if (f == NULL) {
    if (fheader[ihdr].nblocks > 0) {
      cfac_errmsg(cfac, "A single file for one DB type must be used in one session.\n");
      return NULL;
    }
    f = fopen(fn, "wb");
    if (f == NULL) {
      cfac_errmsg(cfac, "cannot open file %s\n", fn);
      return NULL;
    }
  } else {
    if (fheader[ihdr].nblocks == 0) {
      fclose(f);
      f = fopen(fn, "wb");
    }
  }

  if (f == NULL) {
    cfac_errmsg(cfac, "cannot open file %s\n", fn);
    return NULL;
  }

  fheader[ihdr].type = fhdr->type;
  memcpy(fheader[ihdr].symbol, fhdr->symbol, 2);
  fheader[ihdr].symbol[2] = '\0';
  fheader[ihdr].atom = fhdr->atom;
  WriteFHeader(f, &(fheader[ihdr]));

  return f;
}

int CloseFile(FILE *f, F_HEADER *fhdr) {
  int ihdr;

  ihdr = fhdr->type-1;
  fseek(f, 0, SEEK_SET);
  fheader[ihdr].type = fhdr->type;
  WriteFHeader(f, &(fheader[ihdr]));

  fclose(f);
  return 0;
}

int InitFile(FILE *f, F_HEADER *fhdr, void *rhdr) {
  EN_HEADER *en_hdr;
  ENF_HEADER *enf_hdr;
  TR_HEADER *tr_hdr;
  TRF_HEADER *trf_hdr;
  CE_HEADER *ce_hdr;
  CEF_HEADER *cef_hdr;
  CEMF_HEADER *cemf_hdr;
  RR_HEADER *rr_hdr;
  AI_HEADER *ai_hdr;
  AIM_HEADER *aim_hdr;
  CI_HEADER *ci_hdr;
  CIM_HEADER *cim_hdr;
  long int p;

  if (f == NULL) return 0;

  fseek(f, 0, SEEK_END);
  p = ftell(f);

  switch (fhdr->type) {
  case DB_EN:
    en_hdr = (EN_HEADER *) rhdr;
    en_header.position = p;
    en_header.length = 0;
    en_header.nele = en_hdr->nele;
    en_header.nlevels = 0;
    break;
  case DB_TR:
    tr_hdr = (TR_HEADER *) rhdr;
    memcpy(&tr_header, tr_hdr, sizeof(TR_HEADER));
    tr_header.position = p;
    tr_header.length = 0;
    tr_header.ntransitions = 0;
    break;
  case DB_CE:
    ce_hdr = (CE_HEADER *) rhdr;
    memcpy(&ce_header, ce_hdr, sizeof(CE_HEADER));
    ce_header.position = p;
    ce_header.length = 0;
    ce_header.ntransitions = 0;
    break;
  case DB_RR:
    rr_hdr = (RR_HEADER *) rhdr;
    memcpy(&rr_header, rr_hdr, sizeof(RR_HEADER));
    rr_header.position = p;
    rr_header.length = 0;
    rr_header.ntransitions = 0;
    break;
  case DB_AI:
    ai_hdr = (AI_HEADER *) rhdr;
    memcpy(&ai_header, ai_hdr, sizeof(AI_HEADER));
    ai_header.position = p;
    ai_header.length = 0;
    ai_header.ntransitions = 0;
    break;
  case DB_CI:
    ci_hdr = (CI_HEADER *) rhdr;
    memcpy(&ci_header, ci_hdr, sizeof(CI_HEADER));
    ci_header.position = p;
    ci_header.length = 0;
    ci_header.ntransitions = 0;
    break;
  case DB_AIM:
    aim_hdr = (AIM_HEADER *) rhdr;
    memcpy(&aim_header, aim_hdr, sizeof(AIM_HEADER));
    aim_header.position = p;
    aim_header.length = 0;
    aim_header.ntransitions = 0;
    break;
  case DB_CIM:
    cim_hdr = (CIM_HEADER *) rhdr;
    memcpy(&cim_header, cim_hdr, sizeof(CIM_HEADER));
    cim_header.position = p;
    cim_header.length = 0;
    cim_header.ntransitions = 0;
    break;
  case DB_ENF:
    enf_hdr = (ENF_HEADER *) rhdr;
    memcpy(&enf_header, enf_hdr, sizeof(ENF_HEADER));
    enf_header.position = p;
    enf_header.length = 0;
    enf_header.nele = enf_hdr->nele;
    enf_header.nlevels = 0;
    break;
  case DB_TRF:
    trf_hdr = (TRF_HEADER *) rhdr;
    memcpy(&trf_header, trf_hdr, sizeof(TRF_HEADER));
    trf_header.position = p;
    trf_header.length = 0;
    trf_header.ntransitions = 0;
    break;
  case DB_CEF:
    cef_hdr = (CEF_HEADER *) rhdr;
    memcpy(&cef_header, cef_hdr, sizeof(CEF_HEADER));
    cef_header.position = p;
    cef_header.length = 0;
    cef_header.ntransitions = 0;
    break;
  case DB_CEMF:
    cemf_hdr = (CEMF_HEADER *) rhdr;
    memcpy(&cemf_header, cemf_hdr, sizeof(CEMF_HEADER));
    cemf_header.position = p;
    cemf_header.length = 0;
    cemf_header.ntransitions = 0;
    break;
  default:
    break;
  }

  return 0;
}

int DeinitFile(FILE *f, F_HEADER *fhdr) {
  int n;

  if (f == NULL || fhdr->type <= 0) return 0;

  switch (fhdr->type) {
  case DB_EN:
    fseek(f, en_header.position, SEEK_SET);
    if (en_header.length > 0) {
      n = WriteENHeader(f, &en_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_TR:
    fseek(f, tr_header.position, SEEK_SET);
    if (tr_header.length > 0) {
      n = WriteTRHeader(f, &tr_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_CE:
    fseek(f, ce_header.position, SEEK_SET);
    if (ce_header.length > 0) {
      n = WriteCEHeader(f, &ce_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_RR:
    fseek(f, rr_header.position, SEEK_SET);
    if (rr_header.length > 0) {
      n = WriteRRHeader(f, &rr_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_AI:
    fseek(f, ai_header.position, SEEK_SET);
    if (ai_header.length > 0) {
      n = WriteAIHeader(f, &ai_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_CI:
    fseek(f, ci_header.position, SEEK_SET);
    if (ci_header.length > 0) {
      n = WriteCIHeader(f, &ci_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_AIM:
    fseek(f, aim_header.position, SEEK_SET);
    if (aim_header.length > 0) {
      n = WriteAIMHeader(f, &aim_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_CIM:
    fseek(f, cim_header.position, SEEK_SET);
    if (cim_header.length > 0) {
      n = WriteCIMHeader(f, &cim_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_ENF:
    fseek(f, enf_header.position, SEEK_SET);
    if (enf_header.length > 0) {
      n = WriteENFHeader(f, &enf_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_TRF:
    fseek(f, trf_header.position, SEEK_SET);
    if (trf_header.length > 0) {
      n = WriteTRFHeader(f, &trf_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_CEF:
    fseek(f, cef_header.position, SEEK_SET);
    if (cef_header.length > 0) {
      n = WriteCEFHeader(f, &cef_header);
      if (!n) {
        return 1;
      }
    }
    break;
  case DB_CEMF:
    fseek(f, cemf_header.position, SEEK_SET);
    if (cemf_header.length > 0) {
      n = WriteCEMFHeader(f, &cemf_header);
      if (!n) {
        return 1;
      }
    }
    break;
  default:
    break;
  }
  return 0;
}

static int MemENFTable(const cfac_t *cfac, char *fn) {
  F_HEADER fh;
  ENF_HEADER h;
  ENF_RECORD r;
  FILE *f;
  int n, i, nlevels;
  int swp, sr;

  f = fopen(fn, "rb");
  if (f == NULL) return -1;

  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) return 0;

  sr = sizeof(r.ilev) + sizeof(r.energy) + sizeof(r.pbasis);

  if (mem_enf_table) free(mem_enf_table);

  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadENFHeader(f, &h, swp);
    if (n == 0) break;
    nlevels = h.nlevels;
    if (h.length > sr) {
      if (fseek(f, h.length-sr, SEEK_CUR) != 0) {
        cfac_errmsg(cfac, "Error parsing file %s!\n", fn);
        fclose(f);
        return -1;
      }
    }
    n = ReadENFRecord(f, &r, swp);
    if (r.ilev >= nlevels) nlevels = r.ilev+1;
  }

  if (nlevels <= 0) {
    cfac_errmsg(cfac, "No levels found in the DB file %s!\n", fn);
    fclose(f);
    return -1;
  }

  mem_enf_table = (EN_SRECORD *) malloc(sizeof(EN_SRECORD)*nlevels);
  if (!mem_enf_table) {
    return -1;
  }
  mem_enf_table_size = nlevels;

  fseek(f, SIZE_F_HEADER, SEEK_SET);
  while (1) {
    n = ReadENFHeader(f, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENFRecord(f, &r, swp);
      if (n == 0) break;
      mem_enf_table[r.ilev].energy = r.energy;
      DecodeBasisEB(r.pbasis, &mem_enf_table[r.ilev].p, &mem_enf_table[r.ilev].j);
    }
  }

  fclose(f);

  return 0;
}

int PrintTable(const cfac_t *cfac, char *ifn, char *ofn, int v) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;

  f1 = fopen(ifn, "rb");
  if (f1 == NULL) return -1;

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) return -1;

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    goto DONE;
  }

  if (v) {
    if (fh.type == DB_EN && mem_en_table == NULL) {
      MemENTable(cfac, ifn);
    } else
    if (fh.type == DB_ENF && mem_enf_table == NULL) {
      MemENFTable(cfac, ifn);
    }
  }

  if (v && !mem_en_table) {
    cfac_errmsg(cfac, "Energy level table has not been built in memory.\n");
    return -1;
  }
  if (v && fh.type >= DB_ENF && mem_enf_table == NULL) {
    cfac_errmsg(cfac, "Field dependent energy table has not been built in memory.\n");
    goto DONE;
  }

  fprintf(f2, "cFAC %d.%d.%d\n", fh.version, fh.sversion, fh.ssversion);
  fprintf(f2, "Endian\t= %d\n", (int) CheckEndian(&fh));
  fprintf(f2, "TSess\t= %lu\n", fh.tsession);
  fprintf(f2, "Type\t= %d\n", fh.type);
  fprintf(f2, "Verbose\t= %d\n", v);
  fprintf(f2, "%s Z\t= %5.1f\n", fh.symbol, fh.atom);
  fprintf(f2, "NBlocks\t= %d\n", fh.nblocks);
  fflush(f2);
  switch (fh.type) {
  case DB_EN:
    if (v) {
      fprintf(f2, "E0\t= %-d, %15.8E\n",
              iground, (mem_en_table[iground].energy * HARTREE_EV));
    }
    n = PrintENTable(f1, f2, v, swp);
    break;
  case DB_TR:
    n = PrintTRTable(f1, f2, v, swp);
    break;
  case DB_CE:
    n = PrintCETable(f1, f2, v, swp);
    break;
  case DB_RR:
    n = PrintRRTable(f1, f2, v, swp);
    break;
  case DB_AI:
    n = PrintAITable(f1, f2, v, swp);
    break;
  case DB_CI:
    n = PrintCITable(f1, f2, v, swp);
    break;
  case DB_AIM:
    n = PrintAIMTable(f1, f2, v, swp);
    break;
  case DB_CIM:
    n = PrintCIMTable(f1, f2, v, swp);
    break;
  case DB_ENF:
    n = PrintENFTable(f1, f2, v, swp);
    break;
  case DB_TRF:
    n = PrintTRFTable(f1, f2, v, swp);
    break;
  case DB_CEF:
    n = PrintCEFTable(f1, f2, v, swp);
    break;
  case DB_CEMF:
    n = PrintCEMFTable(f1, f2, v, swp);
    break;
  default:
    break;
  }

 DONE:
  fclose(f1);
  if (f2 != stdout) fclose(f2);
  else fflush(f2);

  return n;
}

int UTAFromENRecord(const EN_RECORD *r) {
  if (r->j < 0) return 1;
  else return 0;
}

int JFromENRecord(const EN_RECORD *r) {
  if (r->j < 0) return r->ibase;
  else return r->j;
}

int IBaseFromENRecord(const EN_RECORD *r) {
  if (r->j < 0) return -1;
  else return r->ibase;
}

int MemENTable(const cfac_t *cfac, char *fn) {
  F_HEADER fh;
  EN_HEADER h;
  EN_RECORD r;
  FILE *f;
  int n, i, nlevels;
  float e0;
  int swp, sr;

  f = fopen(fn, "rb");
  if (f == NULL) return -1;

  n = ReadFHeader(f, &fh, &swp);
  if (n == 0) return 0;

  if (fh.type != DB_EN) {
    fclose(f);
    if (fh.type == DB_ENF) {
      return MemENFTable(cfac, fn);
    } else {
      cfac_errmsg(cfac, "File type is neither DB_EN nor DB_ENF\n");
      return -1;
    }
  }

  if (version_read[DB_EN-1] < 109) sr = sizeof(EN_RECORD);
  else sr = SIZE_EN_RECORD;

  nlevels = 0;
  for (i = 0; i < fh.nblocks; i++) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    nlevels = h.nlevels;
    if (h.length > sr) {
      if (fseek(f, h.length-sr, SEEK_CUR) != 0) {
        cfac_errmsg(cfac, "Error parsing file %s!\n", fn);
        fclose(f);
        return -1;
      }
    }
    n = ReadENRecord(f, &r, swp);
    if (r.ilev >= nlevels) nlevels = r.ilev+1;
  }

  if (nlevels <= 0) {
    cfac_errmsg(cfac, "No levels found in the DB file %s!\n", fn);
    fclose(f);
    return -1;
  }

  mem_en_table = (EN_SRECORD *) malloc(sizeof(EN_SRECORD)*nlevels);
  if (!mem_en_table) {
    return -1;
  }
  mem_en_table_size = nlevels;

  e0 = 0.0;
  if (version_read[DB_EN-1] < 109) {
    fseek(f, sizeof(F_HEADER), SEEK_SET);
  } else {
    fseek(f, SIZE_F_HEADER, SEEK_SET);
  }
  while (1) {
    n = ReadENHeader(f, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENRecord(f, &r, swp);
      if (n == 0) break;
      if (r.energy < e0) {
        e0 = r.energy;
        iground = r.ilev;
      }
      mem_en_table[r.ilev].energy = r.energy;
      mem_en_table[r.ilev].p = r.p;
      mem_en_table[r.ilev].j = JFromENRecord(&r);
      mem_en_table[r.ilev].ibase = r.ibase;
    }
  }

  fclose(f);
  return 0;
}

int PrintENTable(FILE *f1, FILE *f2, int v, int swp) {
  EN_HEADER h;
  EN_RECORD r;
  int n, i;
  int nb;
  double e;
  int p, vnl, j, ibase;

  nb = 0;
  while (1) {
    n = ReadENHeader(f1, &h, swp);
    if (n == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NLEV\t= %d\n", h.nlevels);
    fprintf(f2, "  ILEV  IBASE    ENERGY       P   VNL   2J\n");
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENRecord(f1, &r, swp);
      if (n == 0) break;
      e = r.energy;
      if (v) {
        e -= mem_en_table[iground].energy;
        e *= HARTREE_EV;
      }
      if (r.p < 0) {
        p = 1;
        vnl = -r.p;
      } else {
        p = 0;
        vnl = r.p;
      }
      j = JFromENRecord(&r);
      ibase = IBaseFromENRecord(&r);
      fprintf(f2, "%6d %6d %15.8E %1d %5d %4d %-20s %-20s %-s\n",
              r.ilev, ibase, e, p, vnl, j, r.ncomplex, r.sname, r.name);
    }
    nb += 1;
  }

  return nb;
}

int PrintENFTable(FILE *f1, FILE *f2, int v, int swp) {
  ENF_HEADER h;
  ENF_RECORD r;
  int n, i;
  int nb, ilev, mlev;
  double e;

  nb = 0;
  while (1) {
    n = ReadENFHeader(f1, &h, swp);
    if (n == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NLEV\t= %d\n", h.nlevels);
    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);
    fprintf(f2, "%6s %22s %6s %4s\n", "ILEV", "ENERGY", "PBASIS", "M");
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENFRecord(f1, &r, swp);
      if (n == 0) break;
      e = r.energy;
      if (v) {
        e -= mem_en_table[iground].energy;
        e *= HARTREE_EV;
      }
      DecodeBasisEB(r.pbasis, &ilev, &mlev);
      fprintf(f2, "%6d %22.15E %6d %4d\n", r.ilev, e, ilev, mlev);
    }
    nb++;
  }

  return nb;
}

int PrintTRTable(FILE *f1, FILE *f2, int v, int swp) {
  TR_HEADER h;
  TR_RECORD r;
  TR_EXTRA rx;
  int n, i;
  int nb;
  double e, a, gf;

  nb = 0;

  while (1) {
    n = ReadTRHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "MULTIP\t= %d\n", (int)h.multipole);
    fprintf(f2, "GAUGE\t= %d\n", (int)h.gauge);
    fprintf(f2, "MODE\t= %d\n", (int)h.mode);

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadTRRecord(f1, &r, &rx, swp);
      if (n == 0) break;
      if (v) {
        e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
        e += rx.de;

        gf = OscillatorStrength(h.multipole, e, r.rme, &a);
        a /= (mem_en_table[r.upper].j + 1.0);
          fprintf(f2, "%5d %4d %5d %4d %13.6E %11.4E %13.6E %13.6E %13.6E\n",
                  r.upper, mem_en_table[r.upper].j,
                  r.lower, mem_en_table[r.lower].j,
                  (e*HARTREE_EV),
                  (rx.sdev*HARTREE_EV), gf, a*RATE_AU, r.rme);
      } else {
        e = rx.de;
        fprintf(f2, "%5d %5d %13.6E %11.4E %13.6E\n",
              r.upper, r.lower, e, rx.sdev, r.rme);
      }
    }
    nb += 1;
  }

  return nb;
}

int PrintTRFTable(FILE *f1, FILE *f2, int v, int swp) {
  TRF_HEADER h;
  TRF_RECORD r;
  int n, i, j;
  int nb, nq;
  double e, a, gf, ta;

  nb = 0;

  while (1) {
    n = ReadTRFHeader(f1, &h, swp);
    if (n == 0) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "MULTIP\t= %d\n", (int)h.multipole);
    fprintf(f2, "GAUGE\t= %d\n", (int)h.gauge);
    fprintf(f2, "MODE\t= %d\n", (int)h.mode);
    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);
    nq = 2*abs(h.multipole)+1;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadTRFRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (v) {
        e = mem_enf_table[r.upper].energy - mem_enf_table[r.lower].energy;
        ta = 0.0;
        for (j = 0; j < nq; j++) {
          gf = OscillatorStrength(h.multipole, e, (double)(r.strength[j]), &a);
          a *= RATE_AU;
          ta += a;
          fprintf(f2, "%6d %6d %3d %6d %6d %3d %2d %13.6E %13.6E %13.6E %13.6E %13.6E\n",
                  r.upper, mem_enf_table[r.upper].p, mem_enf_table[r.upper].j,
                  r.lower, mem_enf_table[r.lower].p, mem_enf_table[r.lower].j,
                  j-abs(h.multipole),(e*HARTREE_EV), gf, a, r.strength[j], ta);
        }
      } else {
        for (j = 0; j < nq; j++) {
          fprintf(f2, "%6d %6d %13.6E\n",
                  r.upper, r.lower, r.strength[j]);
        }
      }
      free(r.strength);
    }
    nb += 1;
  }

  return nb;
}

int PrintCETable(FILE *f1, FILE *f2, int v, int swp) {
  CE_HEADER h;
  CE_RECORD r;
  int n, i, t;
  int nb;
  int k, p2;
  float a, e = 0.0;

  nb = 0;
  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "MSUB\t= %d\n", h.msub);
    fprintf(f2, "PWTYPE\t= %d\n", h.pw_type);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);

    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (v) {
        e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
        fprintf(f2, "%6d %2d %6d %2d %11.4E %d\n",
                r.lower, mem_en_table[r.lower].j,
                r.upper, mem_en_table[r.upper].j,
                e*HARTREE_EV, r.nsub);
        fprintf(f2, "%11.4E %11.4E %11.4E\n",
                r.bethe, r.born[0], r.born[1]*HARTREE_EV);
      } else {
        fprintf(f2, "%6d %6d %d\n",
                r.lower, r.upper, r.nsub);
        fprintf(f2, "%11.4E %11.4E %11.4E\n",
                r.bethe, r.born[0], r.born[1]);
      }

      p2 = 0;
      for (k = 0; k < r.nsub; k++) {
        if (h.msub) {
          fprintf(f2, "%11.4E\n", r.params[k]);
        }
        for (t = 0; t < h.n_usr; t++) {
          if (v) {
            a = h.usr_egrid[t];
            if (h.usr_egrid_type == 1) a += e;
            a *= 2.0*(1.0 + 0.5*FINE_STRUCTURE_CONST2 * a);
            a = M_PI * AREA_AU20/a;
            if (!h.msub) a /= (mem_en_table[r.lower].j+1.0);
            a *= r.strength[p2];
            fprintf(f2, "%11.4E %11.4E %11.4E\n",
                    h.usr_egrid[t]*HARTREE_EV,
                    r.strength[p2], a);
          } else {
            fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[p2]);
          }
          p2++;
        }
        if (k < r.nsub-1) {
          fprintf(f2, "--------------------------------------------\n");
        }
      }
      fflush(f2);
      if (h.msub) free(r.params);
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    nb += 1;
  }

  return nb;
}

int PrintCEFTable(FILE *f1, FILE *f2, int v, int swp) {
  CEF_HEADER h;
  CEF_RECORD r;
  int n, i, t;
  int nb;
  float a, e = 0.0;

  nb = 0;

  while (1) {
    n = ReadCEFHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);

    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }

    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCEFRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (v) {
        e = mem_enf_table[r.upper].energy - mem_enf_table[r.lower].energy;
        fprintf(f2, "%6d %6d %3d %6d %2d %3d %11.4E\n",
                r.lower, mem_enf_table[r.lower].p, mem_enf_table[r.lower].j,
                r.upper, mem_enf_table[r.upper].p, mem_enf_table[r.upper].j,
                e*HARTREE_EV);
        fprintf(f2, "%11.4E %11.4E %11.4E\n",
                r.bethe, r.born[0], r.born[1]*HARTREE_EV);
      } else {
        fprintf(f2, "%6d %6d\n",
                r.lower, r.upper);
        fprintf(f2, "%11.4E %11.4E %11.4E\n",
                r.bethe, r.born[0], r.born[1]);
      }

      for (t = 0; t < h.n_egrid; t++) {
        if (v) {
          a = h.egrid[t];
          a += e;
          a *= 2.0*(1.0 + 0.5*FINE_STRUCTURE_CONST2 * a);
          a = M_PI * AREA_AU20/a;
          a *= r.strength[t];
            fprintf(f2, "%11.4E %11.4E %11.4E\n",
                    h.egrid[t]*HARTREE_EV,
                    r.strength[t], a);
        } else {
          fprintf(f2, "%11.4E %11.4E\n", h.egrid[t], r.strength[t]);
        }
      }
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    nb += 1;
  }

  return nb;
}

int PrintCEMFTable(FILE *f1, FILE *f2, int v, int swp) {
  CEMF_HEADER h;
  CEMF_RECORD r;
  int n, i, t;
  int nb;
  int k, na, ith, iph;
  float a, e = 0.0;

  nb = 0;

  while (1) {
    n = ReadCEMFHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);

    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "TE0\t= %15.8E\n", h.te0 * HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }

    fprintf(f2, "NTHETA\t= %d\n", h.n_thetagrid);
    for (i = 0; i < h.n_thetagrid; i++) {
      fprintf(f2, "\t %15.8E\n", h.thetagrid[i]);
    }
    fprintf(f2, "NPHI\t= %d\n", h.n_phigrid);
    for (i = 0; i < h.n_phigrid; i++) {
      fprintf(f2, "\t %15.8E\n", h.phigrid[i]);
    }

    fprintf(f2, "EFIELD\t= %15.8E\n", h.efield);
    fprintf(f2, "BFIELD\t= %15.8E\n", h.bfield);
    fprintf(f2, "FANGLE\t= %15.8E\n", h.fangle);

    na = h.n_thetagrid * h.n_phigrid;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCEMFRecord(f1, &r, swp, &h);
      if (n == 0) break;
      k = 0;
      for (ith = 0; ith < h.n_thetagrid; ith++) {
        for (iph = 0; iph < h.n_phigrid; iph++) {
          if (v) {
            e = mem_enf_table[r.upper].energy - mem_enf_table[r.lower].energy;
            fprintf(f2, "%6d %6d %3d %6d %2d %3d %11.4E\n",
                    r.lower, mem_enf_table[r.lower].p, mem_enf_table[r.lower].j,
                    r.upper, mem_enf_table[r.upper].p, mem_enf_table[r.upper].j,
                    e*HARTREE_EV);
            fprintf(f2, "%11.4E %11.4E %11.4E %11.4E %11.4E\n",
                    h.thetagrid[ith]*180.0/M_PI, h.phigrid[iph]*180.0/M_PI,
                    r.bethe[k], r.born[k], r.born[na]*HARTREE_EV);
          } else {
            fprintf(f2, "%6d %6d\n",
                    r.lower, r.upper);
            fprintf(f2, "%11.4E %11.4E %11.4E\n",
                    r.bethe[k], r.born[k], r.born[na]);
          }
          for (t = 0; t < h.n_egrid; t++) {
            if (v) {
              a = h.egrid[t];
              a += e;
              a *= 2.0*(1.0 + 0.5*FINE_STRUCTURE_CONST2 * a);
              a = M_PI * AREA_AU20/a;
              a *= r.strength[t+h.n_egrid*k];
              fprintf(f2, "%11.4E %11.4E %11.4E\n",
                      h.egrid[t]*HARTREE_EV,
                      r.strength[t+h.n_egrid*k], a);
            } else {
              fprintf(f2, "%11.4E %11.4E\n", h.egrid[t], r.strength[t+h.n_egrid*k]);
            }
          }
          k++;
        }
      }
      free(r.strength);
      free(r.bethe);
      free(r.born);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.thetagrid);
    free(h.phigrid);
    nb += 1;
  }

  return nb;
}

int PrintRRTable(FILE *f1, FILE *f2, int v, int swp) {
  RR_HEADER h;
  RR_RECORD r;
  int n, i, t;
  int nb;
  float e = 0.0, eph, ee, phi, rr;

  nb = 0;
  while (1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "MULTIP\t= %d\n", h.multipole);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);
    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;

      if (v) {
        e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
        fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
                r.b, mem_en_table[r.b].j,
                r.f, mem_en_table[r.f].j,
                (e*HARTREE_EV), r.kl);
      } else {
        fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.kl);
      }

      if (h.qk_mode == QK_FIT) {
        for (t = 0; t < h.nparams; t++) {
          if (v && t == h.nparams-1) {
            fprintf(f2, "%11.4E ", r.params[t]*HARTREE_EV);
          } else {
            fprintf(f2, "%11.4E ", r.params[t]);
          }
        }
        fprintf(f2, "\n");
      }

      for (t = 0; t < h.n_usr; t++) {
        if (v) {
          if (h.usr_egrid_type == 0) {
            eph = h.usr_egrid[t];
            ee = eph - e;
          } else {
            ee = h.usr_egrid[t];
            eph = ee + e;
          }
          phi = FINE_STRUCTURE_CONST2*ee;
          phi = 2.0*M_PI*FINE_STRUCTURE_CONST*r.strength[t]*AREA_AU20;
          rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
          rr /= 1.0+0.5*FINE_STRUCTURE_CONST2*ee;
          phi /= (mem_en_table[r.b].j + 1.0);
          rr /= (mem_en_table[r.f].j + 1.0);
          fprintf(f2, "%11.4E %11.4E %11.4E %11.4E\n",
                  h.usr_egrid[t]*HARTREE_EV, rr, phi, r.strength[t]);
        } else {
          fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[t]);
        }
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }

    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);

    nb++;
  }

  return nb;
}

int PrintAITable(FILE *f1, FILE *f2, int v, int swp) {
  AI_HEADER h;
  AI_RECORD r;
  int n, i;
  int nb;
  float e, sdr, er;

  nb = 0;

  while (1) {
    n = ReadAIHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "EMIN\t= %15.8E\n", h.emin*HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadAIRecord(f1, &r, swp);
      if (n == 0) break;
      if (v) {
        e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
        if (e < 0) er = e - h.emin;
        else er = e;
        sdr = 0.5*(mem_en_table[r.b].j + 1.0);
        sdr *= M_PI*M_PI*r.rate/(er*(mem_en_table[r.f].j + 1.0));
        sdr *= AREA_AU20*HARTREE_EV;
        fprintf(f2, "%6d %2d %6d %2d %11.4E %11.4E %11.4E\n",
                r.b, mem_en_table[r.b].j,
                r.f, mem_en_table[r.f].j,
                e*HARTREE_EV, (RATE_AU*r.rate), sdr);
      } else {
        fprintf(f2, "%6d %6d %15.8E\n", r.b, r.f, r.rate);
      }
    }

    free(h.egrid);
    nb++;
  }

  return nb;
}

int PrintAIMTable(FILE *f1, FILE *f2, int v, int swp) {
  AIM_HEADER h;
  AIM_RECORD r;
  int n, i, m;
  int nb;
  float e;
  double u = AREA_AU20*HARTREE_EV;

  nb = 0;
  while (1) {
    n = ReadAIMHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "EMIN\t= %15.8E\n", h.emin*HARTREE_EV);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);

    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadAIMRecord(f1, &r, swp);
      if (n == 0) break;
      if (v) {
        e = mem_en_table[r.b].energy - mem_en_table[r.f].energy;
        fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
                r.b, mem_en_table[r.b].j,
                r.f, mem_en_table[r.f].j,
                e*HARTREE_EV, r.nsub);
        for (m = 0; m < r.nsub; m += 2) {
          fprintf(f2, "%11.4E %11.4E\n",
                  r.rate[m]*RATE_AU, r.rate[m+1]*u);
        }
      } else {
        fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.nsub);
        for (m = 0; m < r.nsub; m += 2) {
          fprintf(f2, "%11.4E %11.4E\n", r.rate[m], r.rate[m+1]);
        }
      }
      free(r.rate);
    }

    free(h.egrid);
    nb++;
  }

  return nb;
}

int PrintCITable(FILE *f1, FILE *f2, int v, int swp) {
  CI_HEADER h;
  CI_RECORD r;
  int n, i, t;
  int nb;
  float e = 0.0, a;

  nb = 0;

  while (1) {
    n = ReadCIHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "QKMODE\t= %d\n", h.qk_mode);
    fprintf(f2, "NPARAMS\t= %d\n", h.nparams);
    fprintf(f2, "PWTYPE\t= %d\n", h.pw_type);
    fprintf(f2, "NTEGRID\t= %d\n", h.n_tegrid);
    for (i = 0; i < h.n_tegrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.tegrid[i]);
      }
    }
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIRecord(f1, &r, swp, &h);
      if (n == 0) break;

      if (v) {
        e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
        fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
                r.b, mem_en_table[r.b].j,
                r.f, mem_en_table[r.f].j,
                e*HARTREE_EV, r.kl);
      } else {
        fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.kl);
      }

      for (t = 0; t < h.nparams; t++) {
        fprintf(f2, "%11.4E ", r.params[t]);
      }
      fprintf(f2, "\n");
      for (t = 0; t < h.n_usr; t++) {
        if (v) {
          a = h.usr_egrid[t];
          if (h.usr_egrid_type == 1) a += e;
          a *= 1.0 + 0.5*FINE_STRUCTURE_CONST2*a;
          a = AREA_AU20/(2.0*a*(mem_en_table[r.b].j + 1.0));
          a *= r.strength[t];
          fprintf(f2, "%11.4E %11.4E %11.4E\n",
                  h.usr_egrid[t]*HARTREE_EV, r.strength[t], a);
        } else {
          fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[t]);
        }
      }
      free(r.params);
      free(r.strength);
    }

    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);

    nb++;
  }

  return nb;
}

int PrintCIMTable(FILE *f1, FILE *f2, int v, int swp) {
  CIM_HEADER h;
  CIM_RECORD r;
  int n, i, t, q;
  int nb, k;
  float e, a;

  nb = 0;

  while (1) {
    n = ReadCIMHeader(f1, &h, swp);
    if (n == 0) break;

    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRANS\t= %d\n", h.ntransitions);
    fprintf(f2, "ETYPE\t= %d\n", h.egrid_type);
    fprintf(f2, "NEGRID\t= %d\n", h.n_egrid);
    for (i = 0; i < h.n_egrid; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.egrid[i]);
      }
    }
    fprintf(f2, "UTYPE\t= %d\n", h.usr_egrid_type);
    fprintf(f2, "NUSR\t= %d\n", h.n_usr);
    for (i = 0; i < h.n_usr; i++) {
      if (v) {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]*HARTREE_EV);
      } else {
        fprintf(f2, "\t %15.8E\n", h.usr_egrid[i]);
      }
    }

    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIMRecord(f1, &r, swp, &h);
      if (n == 0) break;

      if (v) {
        e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
        fprintf(f2, "%6d %2d %6d %2d %11.4E %2d\n",
                r.b, mem_en_table[r.b].j,
                r.f, mem_en_table[r.f].j,
                e*HARTREE_EV, r.nsub);
      } else {
        fprintf(f2, "%6d %6d %2d\n", r.b, r.f, r.nsub);
      }

      if (v) {
        q = 0;
        for (k = 0; k < r.nsub; k ++) {
          for (t = 0; t < h.n_usr; t++) {
            a = h.usr_egrid[t];
            if (h.usr_egrid_type == 1) a += e;
            a *= 1.0 + 0.5*FINE_STRUCTURE_CONST2*a;
            a = AREA_AU20/(2.0*a);
            a *= r.strength[q];
            fprintf(f2, "%11.4E %11.4E %11.4E\n",
                    h.usr_egrid[t]*HARTREE_EV, r.strength[q], a);
            q++;
          }
          if (k < r.nsub-1) {
            fprintf(f2, "--------------------------------------------\n");
          }
        }
      } else {
        q = 0;
        for (k = 0; k < r.nsub; k++) {
          for (t = 0; t < h.n_usr; t++) {
            fprintf(f2, "%11.4E %11.4E\n", h.usr_egrid[t], r.strength[q]);
            q++;
          }
        }
      }
      free(r.strength);
    }

    free(h.egrid);
    free(h.usr_egrid);

    nb++;
  }

  return nb;
}

int AppendTable(const cfac_t *cfac, char *fn) {
  F_HEADER fh;
  FILE *f;
  int n, swp;

  f = fopen(fn, "rb");
  if (f == NULL) return -1;
  n = ReadFHeader(f, &fh, &swp);
  if (!n) {
    return 1;
  }
  if (swp) {
    cfac_errmsg(cfac, "File %s is in different byte-order\n", fn);
    fclose(f);
    return -1;
  }
  memcpy(&(fheader[fh.type-1]), &fh, sizeof(F_HEADER));
  fclose(f);

  return 0;
}

int JoinTable(const cfac_t *cfac, char *fn1, char *fn2, char *fn) {
  F_HEADER fh1, fh2;
  FILE *f1, *f2, *f;
  int n, swp1, swp2;
#define NBUF 8192
  char buf[NBUF];

  f1 = fopen(fn1, "rb");
  if (f1 == NULL) return -1;
  f2 = fopen(fn2, "rb");
  if (f2 == NULL) return -1;

  n = ReadFHeader(f1, &fh1, &swp1);
  if (n == 0) {
    fclose(f1);
    fclose(f2);
    return 0;
  }
  n = ReadFHeader(f2, &fh2, &swp2);
  if (n == 0) {
    fclose(f1);
    fclose(f2);
    return 0;
  }
  if (swp1 != swp2) {
    cfac_errmsg(cfac, "Files %s and %s have different byte-order\n", fn1, fn2);
    return -1;
  }
  if (fh1.type != fh2.type) {
    cfac_errmsg(cfac, "Files %s and %s are of different type\n", fn1, fn2);
    return -1;
  }
  if (fh1.atom != fh2.atom) {
    cfac_errmsg(cfac, "Files %s and %s are for different element\n", fn1, fn2);
    return -1;
  }

  f = fopen(fn, "w");
  if (f == NULL) return -1;
  fh1.nblocks += fh2.nblocks;

  WriteFHeader(f, &fh1);
  while (1) {
    n = fread(buf, 1, NBUF, f1);
    if (n > 0) {
      if (n > fwrite(buf, 1, n, f)) {
        cfac_errmsg(cfac, "write error\n");
        return -1;
      }
    }
    if (n < NBUF) break;
  }
  while (1) {
    n = fread(buf, 1, NBUF, f2);
    if (n > 0) {
      if (n > fwrite(buf, 1, n, f)) {
        cfac_errmsg(cfac, "write error\n");
        return -1;
      }
    }
    if (n < NBUF) break;
  }

  fclose(f1);
  fclose(f2);
  fclose(f);

  return 0;
#undef NBUF
}

int SaveLevels(const cfac_t *cfac, const char *fn, int start, int n) {
  EN_HEADER en_hdr;
  F_HEADER fhdr;
  FILE *f = NULL;
  int k;
  int nele0;

  nele0 = -1;

  if (n < 0) {
    n = cfac->n_levels - start;
  }

  fhdr.type = DB_EN;
  strcpy(fhdr.symbol, cfac_get_atomic_symbol(cfac));
  fhdr.atom = cfac_get_atomic_number(cfac);

  f = OpenFile(cfac, fn, &fhdr);
  if (!f) {
    return -1;
  }

  for (k = 0; k < n; k++) {
    STATE sp, *s;
    LEVEL *lev;
    EN_RECORD r;
    int p, j, nele, vnl, ibase;
    char sname[LEVEL_NAME_LEN];
    char nc[LEVEL_NAME_LEN];
    int size;

    int i = start + k;
    lev = GetLevel(cfac, i);

    if (lev->uta) {
      sp.kgroup = lev->uta_cfg_g;
      sp.kcfg = lev->uta_g_cfg;
      sp.kstate = 0;
      s = &sp;

      p = lev->uta_p;
      j = -1;

      ibase = lev->uta_g - 1;
    } else {
      SYMMETRY *sym;

      sym = GetSymmetry(cfac, lev->pj);
      s = GetSymmetryState(sym, lev->pb);

      DecodePJ(lev->pj, &p, &j);

      ibase = lev->ibase;
    }

    r.j = j;
    r.ilev = i;
    r.ibase = ibase;
    r.energy = lev->energy;

    nele = ConstructLevelName(cfac, s, NULL, sname, nc, &vnl);
    if (p == 0) {
      r.p = vnl;
    } else {
      r.p = -vnl;
    }
    if (LNAME > strlen(lev->name)) {
        size = strlen(lev->name) + 1;
    } else {
        size = LNAME;
    }
    memcpy(r.name, lev->name, size);
    memcpy(r.sname, sname, LSNAME);
    memcpy(r.ncomplex, nc, LNCOMPLEX);
    r.name[LNAME-1] = '\0';
    r.sname[LSNAME-1] = '\0';
    r.ncomplex[LNCOMPLEX-1] = '\0';

    if (nele != nele0) {
      if (nele0 >= 0) {
        DeinitFile(f, &fhdr);
      }
      nele0 = nele;
      en_hdr.nele = nele;
      InitFile(f, &fhdr, &en_hdr);
    }
    WriteENRecord(f, &r);
  }

  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

  return 0;
}

int ISearch(int i, int n, int *ia) {
  int k;

  for (k = 0; k < n; k++) {
    if (ia[k] == i) {
      return k;
    }
  }

  return -1;
}


static int tr_sink(const cfac_t *cfac,
    const cfac_rtrans_data_t *rtdata, void *udata)
{
    TR_RECORD r;
    TR_EXTRA rx;
    FILE *f = udata;

    r.upper = rtdata->ii;
    r.lower = rtdata->fi;
    r.rme   = rtdata->rme;

    rx.de   = rtdata->uta_de;
    rx.sdev = rtdata->uta_sd;

    WriteTRRecord(f, &r, &rx);

    return 0;
}

int SaveTransition(cfac_t *cfac,
    unsigned nlow, unsigned *low, unsigned nup, unsigned *up,
    const char *fn, int mpole) {
  int mode;
  TR_HEADER tr_hdr;
  F_HEADER fhdr;
  FILE *f;

  if (nlow <= 0 || nup <= 0) {
    return -1;
  }

  if (mpole == 1) { /* always FR for M1 transitions */
    mode = M_FR;
  } else {
    mode = GetTransitionMode(cfac);
  }

  fhdr.type = DB_TR;
  strcpy(fhdr.symbol, cfac_get_atomic_symbol(cfac));
  fhdr.atom = cfac_get_atomic_number(cfac);

  tr_hdr.nele = GetNumElectrons(cfac, low[0]);
  tr_hdr.multipole = mpole;
  tr_hdr.gauge = GetTransitionGauge(cfac);
  tr_hdr.mode = mode;

  f = OpenFile(cfac, fn, &fhdr);
  if (!f) {
    return -1;
  }
  InitFile(f, &fhdr, &tr_hdr);

  crac_calculate_rtrans(cfac, nlow, low, nup, up, mpole, mode, tr_sink, f);

  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

  ReinitRadial(cfac, 1);

  return 0;
}

int StoreTable(const cfac_t *cfac,
    sqlite3 *db, unsigned long int sid, const char *ifn)
{
    F_HEADER fh;
    FILE *fp;
    int n, swp;
    int retval = 0;

    fp = fopen(ifn, "rb");
    if (fp == NULL) {
        return -1;
    }

    n = ReadFHeader(fp, &fh, &swp);
    if (n == 0) {
        fclose(fp);
        return -1;
    }

    switch (fh.type) {
    case DB_EN:
        retval = StoreENTable(cfac, db, sid, fp, swp);
        break;
    case DB_TR:
        retval = StoreTRTable(cfac, db, sid, fp, swp);
        break;
    case DB_CE:
        retval = StoreCETable(cfac, db, sid, fp, swp);
        break;
    case DB_RR:
        retval = StoreRRTable(cfac, db, sid, fp, swp);
        break;
    case DB_AI:
        retval = StoreAITable(cfac, db, sid, fp, swp);
        break;
    case DB_CI:
        retval = StoreCITable(cfac, db, sid, fp, swp);
        break;
    default:
        break;
    }

    fclose(fp);

    return retval;
}
