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

#ifndef _RECOUPLE_H_
#define _RECOUPLE_H_ 1

/*************************************************************
  Header for module "recouple".
  This module calculates the recoupling coefficients.

  The main task is to determine which electrons are the
  interacting ones, and calculate the reduced matrix elements
  of the operator Z and ZxZ0,

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

#include "config.h"

/*
** STRUCT:      INTERACT_SHELL
** PURPOSE:     stores information of an itneracting shell.
** FIELDS:      {int index},
**              the index of the shell within the SHELL_STATE.
**              {int n, j, kl, kappa},
**              the principle quantum number, the angular momentum,
**              the orbital angular momentum, and the relativistic
**              angular quantum number of the shell.
**              {int nq_bra, nq_ket},
**              the occupation numbers int the bra and the ket states.
** NOTE:
*/
typedef struct _INTERACT_SHELL_ {
  int index;
  int n;
  int j;
  int kl;
  int kappa;
  int nq_bra;
  int nq_ket;
} INTERACT_SHELL;

/*
** STRUCT:      INTERACT_DATUM
** PURPOSE:     the information about interacting shells to be
**              saved for later use.
** FIELDS:      {SHELL *bra},
**              the shell structure of the bra state.
**              {INTERACT_SHELL s[4]},
**              the indexes of all interacting shells,
**              s[0] and s[2] are shells from the bra state,
**              s[1] and s[3] are shells from the ket state.
**              {short n_shells},
**              number of shells in the bra state.
**              {short phase},
**              the phase resulting from the decoupling that depends
**              on the shell structure of the states.
** NOTE:
*/
typedef struct _INTERACT_DATUM_ {
  SHELL *bra;
  INTERACT_SHELL s[4];
  short n_shells;
  short phase;
} INTERACT_DATUM;

/* the coeff of type (Z^k dot Z^k) */
int AngularZxZ0(double **coeff, int **kk, int nk,
                int n_shells, SHELL_STATE *bra, SHELL_STATE *ket,
                INTERACT_SHELL *s, int max_rank);

/* the coeff of type Z^k */
int AngularZ(double **coeff, int **kk, int nk,
             int n_shells, SHELL_STATE *bra, SHELL_STATE *ket,
             INTERACT_SHELL *s1, INTERACT_SHELL *s2, int max_rank);

int GetInteract(cfac_t *cfac, INTERACT_DATUM **idatum,
                SHELL_STATE **sbra,
                SHELL_STATE **sket,
                int kgi, int kgj,
                int kci, int kcj,
                int ki, int kj, int bf);

int SetMaxRank(cfac_t *cfac, int k);
int GetMaxRank(const cfac_t *cfac);

int ReinitRecouple(cfac_t *cfac);

#endif
