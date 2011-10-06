#include <stdio.h>

#include <gsl/gsl_ieee_utils.h>

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

int InitFac(void) {
  int ierr;

#if FAC_DEBUG
  debug_log = fopen("debug.log", "w");
#endif

#ifdef PERFORM_STATISTICS
  perform_log = fopen("perform.log", "w");
#endif

  ierr = InitConfig();
  if ( ierr < 0) {
    printf("initialize failed in InitConfig\n");
    return ierr;
  }

  InitCoulomb();
  InitAngular();
  InitRecouple();

  ierr = InitRadial();
  if (ierr < 0) {
    printf("initialize failed in InitRadial\n");
    return ierr;
  }
  
  InitDBase();
  InitStructure();
  InitExcitation();
  InitRecombination();
  InitIonization();
  
  gsl_ieee_env_setup();

  return 0;
}

int ReinitFac(int m_config, int m_recouple, int m_radial,
	      int m_dbase, int m_structure, int m_excitation,
	      int m_recombination, int m_ionization) {
  ReinitExcitation(m_excitation);
  ReinitRecombination(m_recombination);
  ReinitIonization(m_ionization);
  ReinitRecouple(m_recouple);
  ReinitRadial(m_radial);
  ReinitDBase(m_dbase);
  ReinitStructure(m_structure);
  ReinitConfig(m_config);

  return 0;
}
