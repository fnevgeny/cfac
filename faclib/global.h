#ifndef _GLOBAL_H_
#define _GLOBAL_H_ 1

#include "cfac.h"

/*************************************************************
  Header file defining some global constants, macros.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

/*
** VARIABLE:    DEBUG_RECOUPLE, DEBUG_STRUCTURE, FAC_DEBUG
** TYPE:        macro constants
** PURPOSE:     debugging flags.
** NOTE:        
*/
#define DEBUG_RECOUPLE  10
#define DEBUG_STRUCTURE 20
#define FAC_DEBUG 0
#if FAC_DEBUG
/*
** VARIABLE:    debug_log
** TYPE:        global.
** PURPOSE:     file handler for the output of debug information.
** NOTE:        it is only defined when FAC_DEBUG is true.
*/
extern FILE *debug_log;
#endif

/*
** VARIABLE:    PERFORM_STATISTICS
** TYPE:        macro constant.
** PURPOSE:     indicate whether the profiling 
**              information should be compiled in.
** NOTE:        normally, should be commented out.
*/

/* #define PERFORM_STATISTICS 1 */

#ifdef PERFORM_STATISTICS
# include <stdio.h>
# include <time.h>
extern FILE *perform_log;
#endif

extern cfac_t *cfac;

#endif

