#include <stdio.h>

#include <gsl/gsl_ieee_utils.h>

#include "cfacP.h"
#include "global.h"
#include "consts.h"
#include "coulomb.h"
#include "recouple.h"
#include "angular.h"
#include "radial.h"
#include "excitation.h"
#include "ionization.h"
#include "recombination.h"
#include "dbase.h"
#include "init.h"

cfac_t *cfac = NULL;

#if FAC_DEBUG
  FILE *debug_log = NULL;
#endif

#ifdef PERFORM_STATISTICS
  FILE *perform_log = NULL;
#endif

int Info(void) {
  printf("========================================\n");
  printf("The Flexible Atomic Code (FAC)\n");
  printf("Version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
  printf("Bugs and suggestions, please contact:\n");
  printf("Ming Feng Gu, mfgu@ssl.berkeley.edu\n");
  printf("========================================\n");
  return 0;
}

int InitFac() {
  int ierr;

#if FAC_DEBUG
  debug_log = fopen("debug.log", "w");
#endif

#ifdef PERFORM_STATISTICS
  perform_log = fopen("perform.log", "w");
#endif

  cfac = cfac_new();

  InitCoulomb();
  InitAngular();
  InitRecouple();

  ierr = InitRadial();
  if (ierr < 0) {
    printf("initialize failed in InitRadial\n");
    return ierr;
  }
  
  InitDBase();
  InitExcitation();
  InitRecombination();
  InitIonization();
  
  gsl_ieee_env_setup();

  return 0;
}
