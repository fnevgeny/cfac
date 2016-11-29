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

#include <stdio.h>

#include <gsl/gsl_ieee_utils.h>

#include "cfacP.h"
#include "global.h"
#include "excitation.h"
#include "ionization.h"
#include "recombination.h"
#include "dbase.h"
#include "init.h"

cfac_t *cfac = NULL;

int Info(void) {
  printf("==========================================\n");
  printf("cFAC-%d.%d.%d http://github.com/fnevgeny/cfac\n\n",
    CFAC_VERSION, CFAC_SUBVERSION, CFAC_SUBSUBVERSION);
  printf("Based on the Flexible Atomic Code (FAC)\n");
  printf("by Ming Feng Gu\n\n");
  printf("Maintained by Evgeny Stambulchik\n");
  printf("==========================================\n");
  return 0;
}

int InitFac() {
  gsl_ieee_env_setup();

  cfac = cfac_new();
  if (!cfac) {
    printf("Initialization failed\n");
    return -1;
  }

  InitDBase();
  InitExcitation();
  InitRecombination();
  InitIonization(cfac);
  
  return 0;
}
