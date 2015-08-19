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

/*************************************************************
  Implementation of the module "recouple".
  This module calculates the recoupling coefficients. 

  The main task is to determine which electrons are the 
  interacting ones, and calculate the reduced matrix elements
  of the operator Z and ZxZ0, 

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cfacP.h"
#include "angular.h"
#include "rcfp.h"
#include "recouple.h"

static void InitInteractDatum(void *p, int n) {
  INTERACT_DATUM *d;
  int i;

  d = (INTERACT_DATUM *) p;
  for (i = 0; i < n; i++) {
    d[i].n_shells = 0;
    d[i].bra = NULL;
  }
}

/* 
** FUNCTION:    FreeInteractDatum
** PURPOSE:     free memory of an INTERACT_DATUM struct.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
static void FreeInteractDatum(void *p) {
  INTERACT_DATUM *d;
  
  if (!p) return;
  d = (INTERACT_DATUM *) p;
  if (d->n_shells > 0) {
    free(d->bra);
    d->bra = NULL;
    d->n_shells = -1;
  }
}
  
/* 
** FUNCTION:    SetMaxRank
** PURPOSE:     set the maximum rank of the operators.
** INPUT:       {int k},
**              the maximum rank. it should be twice the 
**              actual value.
** RETURN:      {int}
**              always 0.
** SIDE EFFECT: the static cfac->recouple.max_rank is set to k.
** NOTE:        
*/
int SetMaxRank(cfac_t *cfac, int k) {
  cfac->recouple.max_rank = k;
  return 0;
}

/* 
** FUNCTION:    GetMaxRank
** PURPOSE:     retrieve the maximum rank.
** INPUT:       
** RETURN:      {int},
**              the maximum rank.
** SIDE EFFECT: 
** NOTE:        
*/
int GetMaxRank(const cfac_t *cfac) {
  return cfac->recouple.max_rank;
}

static double DecoupleShellRecursive(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
			      int n_interact, int *interact, int *rank) {
  double coeff, a;
  int Jbra, Jket, j1bra, j1ket, j2bra, j2ket;
  int k1, k2, k;

  coeff = 0.0;
  if (n_interact > 0) {
    if (n_shells >= n_interact && interact[0] < n_shells) {
      if (n_shells == 1) return 1.0; /* coeff. is 1, if one shell left */
      Jbra = (bra[0]).totalJ;
      Jket = (ket[0]).totalJ;
      j1bra = (bra[1]).totalJ;
      j1ket = (ket[1]).totalJ;
      j2bra = (bra[0]).shellJ;
      j2ket = (ket[0]).shellJ;

      k = rank[0];
      if (interact[0] < n_shells-1) {
	k2 = 0;
	k1 = rank[0];
      } else {
	k2 = rank[1];
	k1 = (n_interact == 1) ? 0: rank[2];
      }
      
      /* it is a 9j symbol, although in most cases it reduces to a
	 6j symbol or even a number. shall distinguish them if maximum
	 efficiency is needed here */
      coeff = (sqrt((Jbra+1.0)*(Jket+1.0)*(rank[0]+1.0)) *
	       W9j(j1bra, j2bra, Jbra,
		   j1ket, j2ket, Jket,
		   k1, k2, k));
      if (fabs(coeff) < EPS30) return 0.0;
      if (interact[0] < n_shells - 1) {
	/* if the current shell is not an interacting one, the two shells
	   in bra and ket state should be identical. and the reduced matrix
	   element of the identity operator should be included */
	coeff *= sqrt(j2bra + 1);
      } else {
	/* otherwise proceed to the next interacting one */
	n_interact--;
	interact++;
	rank += 2;
      }
      /* strip the outmost shell and call DecoupleShell recursively */
      n_shells--;
      bra++;
      ket++;
      a = DecoupleShellRecursive(n_shells, bra, ket,
				 n_interact, interact, rank);
      coeff *= a;
    }
  } else {
    coeff = sqrt(bra[0].totalJ + 1.0);
  }
  return coeff;
}

/* 
** MACRO:       IsOrder
** PURPOSE:     check the ordering of the shells.
** INPUT:       {int order[4]},
**              the shell indexes to be checked.
**              {int a, b, c, d},
**              the ordering required.
** RETURN:      {int},
**              0: if the order is not of the required type.
**              1: otherwise.
** SIDE EFFECT: 
** NOTE:        
*/
#define IsOrder(order, a, b, c, d) (((order)[0] == (a)) && \
				    ((order)[1] == (b)) && \
				    ((order)[2] == (c)) && \
				    ((order)[3] == (d)))

/* 
** FUNCTION:    DecoupleShell
** PURPOSE:     decouple the operators, so that their reduced 
**              matrix elements are expressed in terms of those 
**              involve individual shells.
** INPUT:       {int n_shells},
**              number of shells in the bra and ket states.
**              {SHELL_STATE *bra, *ket},
**              the shell states. if the bra and ket have different 
**              shell structures, they must be padded with empty shells
**              to make them identical. this is done in GetInteract.
**              {int n_interact},
**              number of interacting shells.
**              {int *interact},
**              the indexes of the interacting shells.
**              {int *rank},
**              the ranks of the operators.
** RETURN:      {double},
**              the decoupling coefficient.
** SIDE EFFECT: 
** NOTE:        
*/
static double DecoupleShell(int n_shells, SHELL_STATE *bra, SHELL_STATE *ket,
		     int n_interact, int *interact, int *rank) {
  int i, j, k;
  double coeff;

  /* check the delta function for non-interacting shells first */
  i = 0;
  for (j = 0; j < n_interact; j++) {
    k = n_shells - interact[j] - 1;
    for (; i < k; i++) {
      if (bra[i].shellJ != ket[i].shellJ ||
	  bra[i].nu != ket[i].nu ||
	  bra[i].Nr != ket[i].Nr) return 0.0;
    }
    i++;
  }
  
  for (; i < n_shells; i++) {
    if (bra[i].shellJ != ket[i].shellJ ||
	bra[i].totalJ != ket[i].totalJ ||
	bra[i].nu != ket[i].nu ||
	bra[i].Nr != ket[i].Nr) return 0.0;
  }

  coeff = DecoupleShellRecursive(n_shells, bra, ket, n_interact, interact, rank);

  return coeff;
}
  
/* 
** FUNCTION:    AngularZ
** PURPOSE:     calculate the reduced matrix element of Z operator.
** INPUT:       {double **coeff},
**              a pointer to a double array, which holds the 
**              results on exit.
**              {int **k},
**              a pointer to an int array, which holds all possible 
**              ranks of the operator.
**              {int n_shells},
**              number of the shells in the bra and ket states.
**              {SHELL_STATE *bra, *ket},
**              the bra and ket states.
**              {INTERACT_SHELL *s1, *s2},
**              two interacting shells involved in the operator.
** RETURN:      {int},
**              0: error occured.
**              1: succeeded.
** SIDE EFFECT: 
** NOTE:        For the definition of the operator see 
**              Bar-Shalom et al. Phys. Rev A. 38, 1773. 

**              if the rank of the operator is passed in, the total 
**              number of ranks and ranks should be set in nk, and **kk, 
**              and the storage for the coeff be provided. 
**              otherwise all possible ranks are determined and 
**              storage for the coeff. and kk allocated in this routine. 

**              Note that the reduced matrix element calculated here 
**              does not include the overall phase factor that arises 
**              from interchanging of electron shells, which depends 
**              on the occupation numbers of the config.
**              However, the phase factor arise from the interchanging 
**              of operators are included. This is also the same for 
**              the next routine, AngularZxZ0
*/   
int AngularZ(double **coeff, int **kk, int nk,
	     int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
	     INTERACT_SHELL *s1, INTERACT_SHELL *s2, int max_rank){
  SHELL_STATE st1, st2;
  int rank[4];
  int interact[2];
  int n_interact;
  int kmin, kmax, k, m;
  double coeff1, coeff2;
  RCFP_STATE rcfp_bra, rcfp_ket;

  if (nk == 0) {
    k = (bra[0]).totalJ;
    m = (ket[0]).totalJ;
    kmin = Max(abs(s1->j - s2->j), abs(k - m));
    kmax = Min(s1->j + s2->j, k + m);
    kmax = Min(kmax, max_rank);
    if (IsOdd(kmin)) kmin++; 
    if (kmax < kmin) return -1;
    nk = (kmax - kmin)/2 + 1;

    (*kk) = malloc(sizeof(int) * nk);

    if (!(*kk)) {
      return -1;
    }
    
    m = 0;
    for (k = kmin; k <= kmax; k += 2) {
      (*kk)[m++] = k;
    }

    (*coeff) = malloc(sizeof(double) * nk);
    if (!(*coeff)) return -1;  
  }

  if (s1->index == s2->index) { /* interacting shells are the same */
    n_interact = 1;
    interact[0] = s1->index;
    for (m = 0; m < nk; m++) {
      rank[0] = (*kk)[m];
      rank[1] = (*kk)[m];
      (*coeff)[m] = DecoupleShell(n_shells, bra, ket, 
				  n_interact, interact, rank);
    }
    st1 = bra[n_shells - s1->index -1];
    st2 = ket[n_shells - s1->index -1];
    rcfp_bra.state = RCFPTermIndex(s1->j, (st1).nu, 
				   (st1).Nr, (st1).shellJ);
    rcfp_bra.nq = s1->nq_bra;
    rcfp_bra.subshellMQ = rcfp_bra.nq - (s1->j + 1)/2;

    rcfp_ket.state = RCFPTermIndex(s1->j, (st2).nu, 
				   (st2).Nr, (st2).shellJ);
    rcfp_ket.nq = s1->nq_ket;
    rcfp_ket.subshellMQ = rcfp_ket.nq - (s1->j + 1)/2;
    for (k = 0; k < nk; k++) {
      (*coeff)[k] *= ReducedW(&rcfp_bra, &rcfp_ket, (*kk)[k]/2, 1, -1);
    }

  } else { /* two different interacting shells */
    if (s1->index < s2->index) {
      n_interact = 2;
      interact[0] = s2->index;
      interact[1] = s1->index;
      rank[1] = s2->j;
      rank[2] = s1->j;
      rank[3] = s1->j;
    } else {
      n_interact = 2;
      interact[0] = s1->index;
      interact[1] = s2->index;
      rank[1] = s1->j;
      rank[2] = s2->j;
      rank[3] = s2->j;
    }

    for (m = 0; m < nk; m++) {
      rank[0] = (*kk)[m];
      (*coeff)[m] = DecoupleShell(n_shells, bra, ket, 
				  n_interact, interact, rank);  
    }
 
    st1 = bra[n_shells - s1->index -1];
    st2 = ket[n_shells - s1->index -1];
    rcfp_bra.state = RCFPTermIndex(s1->j, (st1).nu, 
				   (st1).Nr, (st1).shellJ);
    rcfp_bra.nq = s1->nq_bra;
    rcfp_bra.subshellMQ = rcfp_bra.nq - (s1->j + 1)/2;

    rcfp_ket.state = RCFPTermIndex(s1->j, (st2).nu, 
				   (st2).Nr, (st2).shellJ);
    rcfp_ket.nq = s1->nq_ket;
    rcfp_ket.subshellMQ = rcfp_ket.nq - (s1->j + 1)/2;  
    coeff1 = ReducedA(&rcfp_bra, &rcfp_ket, 1);

    st1 = bra[n_shells - s2->index -1];
    st2 = ket[n_shells - s2->index -1];
    rcfp_bra.state = RCFPTermIndex(s2->j, (st1).nu, 
				   (st1).Nr, (st1).shellJ);
    rcfp_bra.nq = s2->nq_bra;
    rcfp_bra.subshellMQ = rcfp_bra.nq - (s2->j + 1)/2;
    rcfp_ket.state = RCFPTermIndex(s2->j, (st2).nu, 
				   (st2).Nr, (st2).shellJ);
    rcfp_ket.nq = s2->nq_ket;
    rcfp_ket.subshellMQ = rcfp_ket.nq - (s2->j + 1)/2;  
    coeff2 = ReducedA(&rcfp_bra, &rcfp_ket, -1);
    for (m = 0; m < nk; m++) {
      if (s1->index > s2->index && IsEven((s1->j + s2->j - (*kk)[m])/2)) {
	(*coeff)[m] = -(*coeff)[m];
      }
      (*coeff)[m] *= coeff1*coeff2;
    }
  }
 
  for (k = 0; k < nk; k++) {
    (*coeff)[k] /= -sqrt((*kk)[k] + 1);
  }

  return nk;
}

/* 
** FUNCTION:    SumCoeff
** PURPOSE:     perform the summation due to the exchange of
**              two operators.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
static void SumCoeff(double *coeff,  int *kk,  int nk,  int p,
	      double *coeff1, int *kk1, int nk1, int p1,
	      int phase, int j1, int j2, int j3, int j4) {
  int i, j;
  double x;
  
  for (i = 0; i < nk; i++) {
    coeff[i] = 0.0;
    for (j = 0; j < nk1; j++) {
      if (fabs(coeff1[j]) > 0.0) {
	x = W6j(j1, j2, kk[i], j3, j4, kk1[j]);
	if (fabs(x) > 0.0) {
	  x *= sqrt(kk1[j] + 1.0);
	  if (p1 && IsOdd(kk1[j]/2)) x = -x;
	  coeff[i] += x * coeff1[j];
	}
      }
    }
    if (fabs(coeff[i]) > 0.0) {
      coeff[i] *= sqrt(kk[i]+1.0);
      if (p) {
	if (IsOdd(phase + kk[i]/2)) coeff[i] = -coeff[i];
      } else {
	  if (IsOdd(phase)) coeff[i] = -coeff[i];
      }
    }
  }
}

/* 
** FUNCTION:    SortShell
** PURPOSE:     sort the interacting shells.
**              calculate the phase of the interchange.
** INPUT:       {int ns},
**              number of shells.
**              {INTERACT_SHELL *s},
**              shell indexes to be sorted.
**              {int *order},
**              the order returned.
** RETURN:      {int},
**              the phase.
** SIDE EFFECT: 
** NOTE:        
*/
static int SortShell(int ns, INTERACT_SHELL *s, int *order) {
  int i, j, k;
  int phase;

  phase = 0;
  for (j = 0; j < ns-1; j++){
    for (i = ns-1; i > j; i--) {
      if (s[order[i]].index < s[order[i-1]].index) {
	k = order[i];
	order[i] = order[i-1];
	order[i-1] = k;
	phase++;
      }
    }
  }
  return phase;
}

/*
** FUNCTION:    AngularZxZ0
** PURPOSE:     calculate the reduced matrix element of 
**              (Z \dot Z) operator.
** INPUT:       {double **coeff},
**              a pointer to a double array, which holds the 
**              results on exit.
**              {int **k},
**              a pointer to an int array, which holds all possible 
**              ranks of the operator.
**              {int n_shells},
**              number of the shells in the bra and ket states.
**              {SHELL_STATE *bra, *ket},
**              the bra and ket states.
**              {INTERACT_SHELL s[4]},
**              4 interacting shells involved in the operator.
** RETURN:      {int},
**              0: error occured.
**              1: succeeded.
** SIDE EFFECT: 
** NOTE:        see the notes of the routine AngularZ. 
**              Note that the order of the interacting shells is 
**              different from that in the radial part. 
**              e.g. if the order in the slater integral is R(ab, cd), 
**              then the order passing into this routine is a, c, b, d 
*/
int AngularZxZ0(double **coeff, int **kk, int nk,
		int n_shells, SHELL_STATE *bra, SHELL_STATE *ket, 
		INTERACT_SHELL *s, int max_rank) {
  
  SHELL_STATE st1, st2;
  int rank[8];
  int interact[4];
  int nops[4]={0,0,0,0};
  int n_interact;
  int order[4] = {0, 1, 2, 3};
  int kmin, kmax, k, m;
  int qm[4] = {1, -1, 1, -1};
  int recouple_operators = 0;
  int phase;
  double *coeff1;
  int *kk1, nk1;
  int i, j;
  double r, a, b;  
  RCFP_STATE rcfp_bra[4], rcfp_ket[4];

  /* if nk is non positive, allocate the memory for the coeff. and kk */
  if (nk <= 0) {
    kmin = Max(abs(s[0].j-s[1].j), abs(s[2].j-s[3].j));
    kmax = Min(s[0].j+s[1].j, s[2].j+s[3].j);
    kmax = Min(kmax, max_rank);
    if (IsOdd(kmin)) kmin++; 
    if (kmax < kmin) return -1;
    nk = (kmax - kmin)/2 + 1;
    (*kk) = malloc(sizeof(int)*nk);
    if (!(*kk)) {
      return -1;
    }
    m = 0;
    for (k = kmin; k <= kmax; k += 2) {
      (*kk)[m++] = k;
    }

    (*coeff) = malloc(sizeof(double)*nk);
    if (!(*coeff)) return -1;
  }

  coeff1 = *coeff;
  kk1 = *kk;
  nk1 = nk;

  /* sort the order of interacting shells */
  phase = SortShell(4, s, order);
  n_interact = 0;
  j = -1;
  for (i = 3; i >= 0; i--) {
    if (i == 3 || s[order[i]].index != s[order[i+1]].index) {
      n_interact++;
      j++;
      interact[j] = s[order[i]].index;
      nops[j] = 1;
      st1 = bra[n_shells - interact[j] - 1];
      st2 = ket[n_shells - interact[j] - 1];
      rcfp_bra[j].state = RCFPTermIndex(s[order[i]].j, (st1).nu,
					(st1).Nr, (st1).shellJ);
      rcfp_bra[j].nq = s[order[i]].nq_bra;
      rcfp_bra[j].subshellMQ = rcfp_bra[j].nq - (s[order[i]].j + 1)/2;
      rcfp_ket[j].state = RCFPTermIndex(s[order[i]].j, (st2).nu,
					(st2).Nr, (st2).shellJ);
      rcfp_ket[j].nq = s[order[i]].nq_ket;
      rcfp_ket[j].subshellMQ = rcfp_ket[j].nq - (s[order[i]].j + 1)/2;
    } else {
      nops[j]++;
    }
  }
  

  rank[0] = 0;
  
  if (nops[0] == 4) { /* 4 identical interacting shells */
    for (m = 0; m < nk1; m++){
      rank[1] = 0;      
      a = DecoupleShell(n_shells, bra, ket, 
			n_interact, interact, rank); 
      b = ReducedWxW0(rcfp_bra, rcfp_ket, kk1[m]/2, 1, -1, 1, -1);
      coeff1[m] = a*b;
    }
    recouple_operators = 0;
  } else if (nops[0] == 3 && nops[1] == 1) { /* 1 x 3 */
    for (m = 0; m < nk1; m++){
      rank[1] = s[order[0]].j;
      rank[2] = rank[1];
      rank[3] = rank[1];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
      if (order[0] == 0 || order[0] == 1) { 
	coeff1[m] *= ReducedAxW(rcfp_bra, rcfp_ket, kk1[m]/2, rank[1],
				qm[order[1]], qm[order[2]], qm[order[3]]);
	coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[0]]);
      } else {     
	coeff1[m] *= ReducedWxA(rcfp_bra, rcfp_ket, kk1[m]/2, rank[1],
				qm[order[1]], qm[order[2]], qm[order[3]]);
	coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[0]]);
      }
    }
    if (order[0] == 0 || order[0] == 1) recouple_operators = 1;
    else recouple_operators = 2;
  } else if (nops[0] == 1 && nops[1] == 3) { /* 3 x 1 */
    for (m = 0; m < nk1; m++){
      rank[1] = s[order[3]].j;
      rank[2] = rank[1];
      rank[3] = rank[1];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
      if (order[3] == 0 || order[3] == 1) {
	coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
	coeff1[m] *= ReducedAxW(rcfp_bra+1, rcfp_ket+1, kk1[m]/2, rank[1],
				qm[order[0]], qm[order[1]], qm[order[2]]);
      } else {
	coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
	coeff1[m] *= ReducedWxA(rcfp_bra+1, rcfp_ket+1, kk1[m]/2, rank[1],
				qm[order[0]], qm[order[1]], qm[order[2]]);
      }
    }
    if (order[3] == 0 || order[3] == 1) recouple_operators = 3;
    else recouple_operators = 4;  
  } else if (nops[0] == 1 && nops[1] == 2) { /* 1 x 2 x 1 */
    if (!(order[1] == 0 && order[2] == 1) &&
	!(order[1] == 2 && order[2] == 3)) {
      m = 0;
      kmin = abs(s[order[0]].j - s[order[3]].j);
      kmax = Min(2*s[order[1]].j, s[order[0]].j+s[order[3]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }      
    for (m = 0; m < nk1; m++){
      rank[1] = s[order[3]].j;
      rank[2] = s[order[3]].j;
      rank[3] = kk1[m];
      rank[4] = s[order[0]].j;
      rank[5] = rank[4];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
      coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
      coeff1[m] *= ReducedW(rcfp_bra+1, rcfp_ket+1, kk1[m]/2,
			    qm[order[1]], qm[order[2]]);
      coeff1[m] *= ReducedA(rcfp_bra+2, rcfp_ket+2, qm[order[0]]);
    }
    recouple_operators = 5;
  } else if (nops[0] == 2 && nops[1] == 2) { /* 2 x 2 */
    if (!(order[0] == 0 && order[1] == 1) &&
	!(order[0] == 2 && order[1] == 3)) {
      m = 0;
      kmin = 0;
      kmax = 2 * Min(s[order[0]].j, s[order[2]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }

    for (m = 0; m < nk1; m++){
      rank[1] = kk1[m];
      rank[2] = rank[1];
      rank[3] = rank[1];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 

      r = ReducedW(rcfp_bra, rcfp_ket, kk1[m]/2, 
		   qm[order[2]], qm[order[3]]);
      coeff1[m] *= r;

      r = ReducedW(rcfp_bra+1, rcfp_ket+1, kk1[m]/2,
		   qm[order[0]], qm[order[1]]);
      coeff1[m] *= r;

    }
    recouple_operators = 6;
  } else if (nops[0] ==  2 && nops[1] == 1) { /* 1 x 1 x 2 */
    if (!(order[2] == 0 && order[3] == 1) &&
	!(order[2] == 2 && order[3] == 3)) {
      m = 0;
      kmin = abs(s[order[0]].j - s[order[1]].j);
      kmax = Min(2*s[order[3]].j, s[order[0]].j+s[order[1]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }
    for (m = 0; m < nk1; m++) {
      rank[1] = kk1[m];
      rank[2] = kk1[m];
      rank[3] = s[order[1]].j;
      rank[4] = s[order[0]].j;
      rank[5] = rank[4];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 
      coeff1[m] *= ReducedW(rcfp_bra, rcfp_ket, kk1[m]/2,
			    qm[order[2]], qm[order[3]]);
      coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[1]]);
      coeff1[m] *= ReducedA(rcfp_bra+2, rcfp_ket+2, qm[order[0]]);
    }
    recouple_operators = 7;
  } else if (nops[0] == 1 && nops[1] == 1 && nops[2] == 2) { /* 2 x 1 x 1 */
    if (!(order[0] == 0 && order[1] == 1) &&
	!(order[0] == 2 && order[1] == 3)) {
      m = 0;
      kmin = abs(s[order[2]].j - s[order[3]].j);
      kmax = Min(2*s[order[0]].j, s[order[2]].j+s[order[3]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }
    for (m = 0; m < nk1; m++) {
      rank[1] = s[order[3]].j;
      rank[2] = s[order[3]].j;
      rank[3] = s[order[2]].j;
      rank[4] = kk1[m];
      rank[5] = rank[4];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank); 

      r = ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
      coeff1[m] *= r;

      r = ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[2]]);
      coeff1[m] *= r;

      r = ReducedW(rcfp_bra+2, rcfp_ket+2, kk1[m]/2,
		   qm[order[0]], qm[order[1]]);
      coeff1[m] *= r;

    }
    recouple_operators = 8;
  } else { /* 1 x 1 x 1 x 1 */
    i = order[0] + order[1];
    if (i != 1 && i != 5) {
      m = 0;
      kmin = Max(abs(s[order[0]].j - s[order[1]].j), 
		 abs(s[order[2]].j - s[order[3]].j));
      kmax = Min(s[order[0]].j + s[order[1]].j,
		 s[order[2]].j + s[order[3]].j);
      if (IsOdd(kmin)) kmin++; 
      if (kmax < kmin) return -1;
      nk1 = (kmax - kmin)/2 + 1;
      coeff1 = malloc(sizeof(double)*nk1);
      kk1 = malloc(sizeof(int)*nk1);
      if (!kk1 || !coeff1) return -1;
      for (k = kmin; k <= kmax; k += 2) {
	kk1[m++] = k;
      }
    }
    for (m = 0; m < nk1; m++) {
      rank[1] = s[order[3]].j;
      rank[2] = rank[1];
      rank[3] = s[order[2]].j;
      rank[4] = kk1[m];
      rank[5] = s[order[1]].j;
      rank[6] = s[order[0]].j;
      rank[7] = rank[6];
      coeff1[m] = DecoupleShell(n_shells, bra, ket, 
				n_interact, interact, rank);
      coeff1[m] *= ReducedA(rcfp_bra, rcfp_ket, qm[order[3]]);
      coeff1[m] *= ReducedA(rcfp_bra+1, rcfp_ket+1, qm[order[2]]);
      coeff1[m] *= ReducedA(rcfp_bra+2, rcfp_ket+2, qm[order[1]]);
      coeff1[m] *= ReducedA(rcfp_bra+3, rcfp_ket+3, qm[order[0]]);
    }
    recouple_operators = 9;
  }

  /* now adjust the phase factor, and do the recoupling if necessary */
  switch (recouple_operators) {
  case 0:
    break;
  case 1:
    if (IsOrder(order, 0, 1, 2, 3)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 1, 0, 2, 3)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[0].j+s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    }
    break;
  case 2:
    if (IsOrder(order, 2, 0, 1, 3)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[3].j-s[2].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 3, 0, 1, 2)) {
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    }
    break;
  case 3:
    if (IsOrder(order, 1, 2, 3, 0)) {
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 2, 3, 1)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[0].j-s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    }
    break;
  case 4:
    if (IsOrder(order, 0, 1, 2, 3)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 1, 3, 2)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[2].j+s[3].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    }
    break;
  case 5:
    if (IsOrder(order, 2, 0, 1, 3)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[2].j-s[3].j+kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 3, 0, 1, 2)) {     
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 2, 3, 1)) {
      for (m = 0; m < nk1; m++)
	if (IsOdd(phase + (s[0].j-s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 1, 2, 3, 0)) { 
      if (IsEven(phase)) {
	for (m = 0; m < nk1; m++)
	  coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 1, 0, 2, 3)) {
	phase += (s[2].j + s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 3, 0, 2, 1)) {
	phase += (s[1].j-s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 1, 0, 3, 2)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 2, 0, 3, 1)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 1, 2, 3)) {
	phase += (s[0].j+s[1].j+s[2].j+s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 1, 2, 0)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 1, 3, 2)) {
	phase += (s[0].j+s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 1, 3, 0)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      }
      free(coeff1);
      free(kk1);
    }
    break;
  case 6:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++) coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2) ||
		 IsOrder(order, 1, 2, 0, 3)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      }
      free(coeff1);
      free(kk1);
    }
    break;      
  case 7:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      /* do nothing in this case */
    } else if (IsOrder(order, 1, 0, 2, 3)) {
      for (m = 0; m < nk1; m++) {
	if (IsOdd(phase + (s[0].j+s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 3, 2, 0, 1)) {
      for (m = 0; m < nk1; m++) {
	if (IsOdd(phase + (s[2].j+s[3].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 0, 1, 3)) {
	phase += (s[0].j-s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 3, 1, 0, 2)) {
	phase += (s[2].j - s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 0, 1, 2)) {
	phase += (s[0].j+s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase, 
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 1, 2, 0, 3)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 2, 1, 0, 3)) {
	phase += 1;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      }
      free(coeff1);
      free(kk1);
    }
    break; 
  case 8:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      if (IsOdd(phase)) {
	for (m = 0; m < nk1; m++) coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 0, 1, 3, 2)) {

      for (m = 0; m < nk1; m++) {
	if (IsOdd(phase + (s[2].j+s[3].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else if (IsOrder(order, 2, 3, 1, 0)) {

      for (m = 0; m < nk1; m++) {
	if (IsOdd(phase + (s[0].j+s[1].j-kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
      }
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j+s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 2, 3, 1)) {
	phase += (s[2].j-s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 1, 3, 2, 0)) {
	phase += (s[0].j - s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 3, 2, 1)) {
	phase += 1;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase, 
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 1, 2, 0, 3)) {
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 1, 2, 3, 0)) {
	phase += (s[0].j+s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      }
      free(coeff1);
      free(kk1);
    }
    break; 
  case 9:
    if (IsOrder(order, 0, 1, 2, 3) ||
	IsOrder(order, 2, 3, 0, 1)) {
      if (IsOdd(phase)) 
	for (m = 0; m < nk1; m++) coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 0, 1, 3, 2) ||
	       IsOrder(order, 3, 2, 0, 1)) {
      for (m = 0; m < nk1; m++) 
	if (IsOdd(phase + (s[2].j + s[3].j - kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 1, 0, 2, 3) ||
	       IsOrder(order, 2, 3, 1, 0)) {
      for (m = 0; m < nk1; m++) 
	if (IsOdd(phase + (s[0].j + s[1].j - kk1[m])/2)) 
	  coeff1[m] = -coeff1[m];
    } else if (IsOrder(order, 1, 0, 3, 2) ||
	       IsOrder(order, 3, 2, 1, 0)) {
      for (m = 0; m < nk1; m++) 
	if (IsOdd(phase+(s[0].j+s[1].j+s[2].j+s[3].j)/2)) 
	  coeff1[m] = -coeff1[m];
    } else {
      if (IsOrder(order, 0, 2, 1, 3) ||
	  IsOrder(order, 1, 3, 0, 2)) {
	phase += (s[1].j + s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 2, 3, 1) ||
		 IsOrder(order, 3, 1, 0, 2)) {
	phase += (s[2].j - s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 0, 1, 3) ||
		 IsOrder(order, 1, 3, 2, 0)) {
	phase += (s[0].j - s[1].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 2, 0, 3, 1) ||
		 IsOrder(order, 3, 1, 2, 0)) {
	phase += (s[0].j + s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 1, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[3].j, s[2].j);
      } else if (IsOrder(order, 0, 3, 1, 2) ||
		 IsOrder(order, 1, 2, 0, 3)) {
	phase += (s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 0, 3, 2, 1) ||
		 IsOrder(order, 2, 1, 0, 3)) {
	phase += 1;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 0, 1, 2) ||
		 IsOrder(order, 1, 2, 3, 0)) {
	phase += (s[0].j + s[3].j + s[1].j - s[2].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 0, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      } else if (IsOrder(order, 3, 0, 2, 1) ||
		 IsOrder(order, 2, 1, 3, 0)) {
	phase += (s[0].j - s[3].j)/2;
	SumCoeff((*coeff), (*kk), nk, 0, coeff1, kk1, nk1, 1, phase,
		 s[0].j, s[1].j, s[2].j, s[3].j);
      }      
      free(coeff1);
      free(kk1);
    }
    break;
  }

  /* adjust the prefactor */
  for (m = 0; m < nk; m++) {
    /* this factor arise from the definition of Z^k and the 
       scalar product. */
    (*coeff)[m] /= sqrt((*kk)[m] + 1);
    if (IsOdd((*kk)[m]/2)) (*coeff)[m] = -(*coeff)[m];
  }

  return nk;
}

/* Analyze the structure of the configuration of bra and ket to
   determine if they can interact, get the interacting shells that
   must interact, and determine the phase factor of the recoupling
   coeff. The latter does not include the phase resulting from the
   reordering of operators; it is calculated in AngularZxZ0 and
   AngularZ0 */
static int InteractingShells(const CONFIG *cbra, const CONFIG *cket,
                      INTERACT_DATUM **idatum,
		      const SHELL_STATE *csf_i, const SHELL_STATE *csf_j,
		      SHELL_STATE **sbra, SHELL_STATE **sket) {
  int i, j, k, m, pb, pk;
  int n, kl, jj, nq, nq_plus, nq_minus, qd;
  SHELL *bra;
  INTERACT_SHELL *s;
  int interaction[4] = {-1, -1, -1, -1};
  int n_shells;
  int jmin, jmax;

  /* allocate memory. use the sum of two configs to avoid counting the
     exact size in advance */
  (*idatum)->bra = malloc(sizeof(SHELL)*(cbra->n_shells + cket->n_shells));
  if ((*idatum)->bra == NULL) {
    printf("error allocating idatam->bra.\n");
    exit(1);
  }
  
  bra = (*idatum)->bra;
  s = (*idatum)->s;
  for (i = 0; i < 4; i++) {
    s[i].index = -1;
  }
  
  i = 0; 
  j = 0;
  m = 0;
  pb = 0;
  pk = 1;
  nq_plus = 0;
  nq_minus = 0;
  while (1) {
    if (i >= cbra->n_shells) {
      if (j < cket->n_shells) {
	k = -1;
      } else {
	break;
      }
    } else {
      if (j >= cket->n_shells) {
	k = 1;
      } else {
	k = CompareShell(&cbra->shells[i], &cket->shells[j]);
      }
    }
    
    if (k > 0) { /* bra has a shell (i-th) that does not exist in ket */
      if (nq_plus >= 2) break;
      
      if (cbra->shells[i].nq > 0) {
	bra[m] = cbra->shells[i];
	UnpackShell(&bra[m], &n, &kl, &jj, &nq);
	
        nq_plus += nq;
	if (nq_plus > 2) break;
	
        s[pb].index = m;
	s[pb].n = n;
	s[pb].kappa = bra[m].kappa;
	s[pb].j = jj;
	s[pb].kl = kl;
	s[pb].nq_bra = nq;
	s[pb].nq_ket = 0;
	
        pb += 2;
	
        if (nq == 2) {
	  s[pb] = s[0];
	} else {
	  interaction[pb-2] = m;
	}
	
        m++;
      }
      
      i++;
    } else
    if (k < 0) { /* ket has a shell (j-th) that does not exist in bra */
      if (nq_minus >= 2) break;
      
      if (cket->shells[j].nq > 0) {
	bra[m] = cket->shells[j];
	UnpackShell(&bra[m], &n, &kl, &jj, &nq);
	
        nq_minus += nq;
	if (nq_minus > 2) break;
	
        s[pk].index = m;
	s[pk].n = n;
	s[pk].kappa = bra[m].kappa;
	s[pk].j = jj;
	s[pk].kl = kl;
	s[pk].nq_bra = 0;
	bra[m].nq = 0;
	s[pk].nq_ket = nq;
	
        pk += 2;
	
        if (nq == 2) {
	  s[pk] = s[1];
	} else {
	  interaction[pk-2] = m;
	}
        
	m++;
      }
      
      j++;
    } else { /* both bra and ket has the shell */
      bra[m] = cbra->shells[i];
      qd = cbra->shells[i].nq - cket->shells[j].nq;
      if (qd > 0) { /* bra has more electrons in the shell */
	if (nq_plus >= 2) break;
	
        UnpackShell(bra+m, &n, &kl, &jj, &nq);
	
        nq_plus += qd;
	if (nq_plus > 2) break;
	
        s[pb].index = m;
	s[pb].n = n;
	s[pb].kappa = bra[m].kappa;
	s[pb].j = jj;
	s[pb].kl = kl;
	s[pb].nq_bra = nq;
	s[pb].nq_ket = cket->shells[j].nq;
	
        pb += 2;
	
        if (qd == 2) {
	  s[pb] = s[0];
	} else {
	  interaction[pb-2] = m;
	}
      } else if (qd < 0) { /* ket has more electrons in the shell */
	if (nq_minus >= 2) break;
	
        UnpackShell(bra+m, &n, &kl, &jj, &nq);
	
        nq_minus -= qd;
	if (nq_minus > 2) break;
	
        s[pk].index = m;
	s[pk].n = n;
	s[pk].kappa = bra[m].kappa;
	s[pk].j = jj;
	s[pk].kl = kl;
	s[pk].nq_bra = nq;
	s[pk].nq_ket = cket->shells[j].nq;
	
        pk += 2;
	
        if (qd == -2) {
	  s[pk] = s[1];
	} else {
	  interaction[pk-2] = m;
	}
      }
      
      i++;
      j++;
      m++;
    }
  }
  
  n_shells = m;
  if (nq_plus != nq_minus ||
      nq_plus > 2 || nq_minus > 2 ||
      i < cbra->n_shells || j < cket->n_shells) {
    free(bra);
    bra = NULL;
    n_shells = -1;
    goto END;
  }

  if (!csf_i || !csf_j || !sbra || !sket) goto END;

  /* determine the phase factor */
  (*idatum)->phase = 0;
  if (interaction[0] >= 0) {
    if (interaction[0] < interaction[1]) {
      jmin = interaction[0] + 1;
      jmax = interaction[1];
    } else {
      jmin = interaction[1] + 1;
      jmax = interaction[0];
    }
    for (j = jmin; j <= jmax; j++) {
      (*idatum)->phase += bra[j].nq;
    }
  } else {
    (*idatum)->phase += 1;
  }
  
  if (interaction[2] >= 0) {
    if (interaction[2] < interaction[3]) {
      jmin = interaction[2] + 1;
      jmax = interaction[3];
    } else {
      jmin = interaction[3] + 1;
      jmax = interaction[2];
    }
    for (j = jmin; j <= jmax; j++) {
      (*idatum)->phase += bra[j].nq;
    }
  } else {
    (*idatum)->phase += 1;
  }

  if (n_shells > 0) {
    i = 0;
    j = 0;
      
    (*sbra) = calloc(n_shells, sizeof(SHELL_STATE));
    (*sket) = calloc(n_shells, sizeof(SHELL_STATE));
    for (m = 0; m < n_shells; m++) {
      if (i < cbra->n_shells) {
	if (bra[m].n == cbra->shells[i].n &&
	    bra[m].kappa == cbra->shells[i].kappa) {
	  (*sbra)[m] = csf_i[i];
	  i++;
	} else {
	  (*sbra)[m].totalJ = csf_i[i].totalJ;
	}
      }
      if (j < cket->n_shells) {
	if (bra[m].n == cket->shells[j].n &&
	    bra[m].kappa == cket->shells[j].kappa) {
	  (*sket)[m] = csf_j[j];
	  j++;
	} else {
	  (*sket)[m].totalJ = csf_j[j].totalJ;
	}	
      }
    }
  }

 END:
  if (n_shells == 0) n_shells = -1;
  (*idatum)->n_shells = n_shells;
  if (n_shells > 0) {
    (*idatum)->bra = realloc(bra, sizeof(SHELL)*n_shells);
    /* adjust the index so that it counts from inner shells */
    for (i = 0; i < 4; i++) {
      if (s[i].index >= 0)
	s[i].index = n_shells - 1 - s[i].index;
    }
  }

  return n_shells;
}

/* 
** FUNCTION:    GetInteract
** PURPOSE:     determing which shells can be interacting.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int GetInteract(cfac_t *cfac, INTERACT_DATUM **idatum,
		SHELL_STATE **sbra, 
		SHELL_STATE **sket, 
		int kgi, int kgj,
		int kci, int kcj, 
		int ki, int kj, int ifb) {
  int i, j, m;
  CONFIG *ci, *cj, cip;
  SHELL_STATE *csf_i, *csf_j, *csf_ip;
  SHELL *bra;
  int n_shells;
  int index[4];

  ci = GetConfigFromGroup(cfac, kgi, kci);
  cj = GetConfigFromGroup(cfac, kgj, kcj);
  if (ci->n_csfs > 0) {
    csf_i = ci->csfs + ki;
    csf_j = cj->csfs + kj;
  } else {
    csf_i = NULL;
    csf_j = NULL;
  }
  if (ci->n_shells <= 0 || cj->n_shells <= 0) return -1;
  if (abs(ci->n_shells+ifb - cj->n_shells) > 2) return -1;

  if (csf_i != NULL) {
    n_shells = -1;
    /* check if this is a repeated call,
     * if not, search in the array.
     */
    if (*idatum == NULL) {
      index[0] = kgi;
      index[1] = kgj;
      index[2] = kci;
      index[3] = kcj;
      (*idatum) = (INTERACT_DATUM *) MultiSet(cfac->recouple.int_shells, index, 
					      NULL);
    }
    if ((*idatum)->n_shells < 0) return -1;
  } else {
    (*idatum) = malloc(sizeof(INTERACT_DATUM));
    (*idatum)->n_shells = 0;
  }
  if ((*idatum)->n_shells > 0) {
    n_shells = (*idatum)->n_shells;
    bra = (*idatum)->bra;
    i = 0;
    j = 0;
    (*sbra) = calloc(n_shells, sizeof(SHELL_STATE));
    (*sket) = calloc(n_shells, sizeof(SHELL_STATE));
    for (m = 0; m < n_shells; m++) {
      if (i < ci->n_shells) {
	if (bra[m].n == ci->shells[i].n &&
	    bra[m].kappa == ci->shells[i].kappa) {
	  (*sbra)[m] = csf_i[i];
	  i++;
	} else {
	  (*sbra)[m].totalJ = csf_i[i].totalJ;
	}
      }
      if (j < cj->n_shells) {
	if (bra[m].n == cj->shells[j].n &&
	    bra[m].kappa == cj->shells[j].kappa) {
	  (*sket)[m] = csf_j[j];
	  j++;
	} else {
	  (*sket)[m].totalJ = csf_j[j].totalJ;
	}	
      }
    }
    /* the quantum number for the free electron must be reset */
    if (ifb) {
      (*sbra)[0].shellJ = 1;
      (*sbra)[0].totalJ = 1;
      (*sbra)[0].nu = 1;
      (*sbra)[0].Nr = 0;
    }
  } else {
    if (ifb) {
      cip.n_shells = ci->n_shells + 1;
      cip.shells = malloc(sizeof(SHELL)*cip.n_shells);
      cip.shells[0].n = 9999;
      cip.shells[0].nq = 1;
      cip.shells[0].kappa = -1;
      memcpy(cip.shells+1, ci->shells, sizeof(SHELL)*ci->n_shells);
      if (csf_i) {
	cip.csfs = malloc(sizeof(SHELL_STATE)*cip.n_shells);
	cip.n_csfs = 1;
	csf_ip = cip.csfs;
	csf_ip[0].shellJ = 1;
	csf_ip[0].totalJ = 1;
	csf_ip[0].nu = 1;
	csf_ip[0].Nr = 0;
	memcpy(csf_ip+1, csf_i, sizeof(SHELL_STATE)*ci->n_shells);
      } else {
	cip.n_csfs = 0;
	csf_ip = NULL;
      }
      n_shells = InteractingShells(&cip, cj, idatum, csf_ip, csf_j, sbra, sket);
      free(cip.shells);
      if (csf_i) {
	free(cip.csfs);
      }
    } else {
      n_shells = InteractingShells(ci, cj, idatum, csf_i, csf_j, sbra, sket);
    }
  }

  if (n_shells < 0 && csf_i == NULL) {
    free((*idatum));
    idatum = NULL;
  }

  return n_shells;
}

/* 
** FUNCTION:    cfac_init_recouple
** PURPOSE:     Initialize the module "recouple"
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int cfac_init_recouple(cfac_t *cfac) {
    int blocks[4] = {10, 10, 64, 64};
    int ndim = 4;
    cfac->recouple.max_rank = MAXRANK;
    cfac->recouple.int_shells = malloc(sizeof(MULTI));
    if (!cfac->recouple.int_shells) {
        return -1;
    }

    return MultiInit(cfac->recouple.int_shells, sizeof(INTERACT_DATUM),
        ndim, blocks, FreeInteractDatum, InitInteractDatum);
}

void cfac_free_recouple(cfac_t *cfac)
{
    MultiFree(cfac->recouple.int_shells);
    free(cfac->recouple.int_shells);
}

/* 
** FUNCTION:    ReinitRecouple
** PURPOSE:     Reinitialize the module "recouple"
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
int ReinitRecouple(cfac_t *cfac) {
  MultiFreeData(cfac->recouple.int_shells);  
  return 0;
}
