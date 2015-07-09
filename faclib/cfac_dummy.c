#include <stdlib.h>
#include <stdio.h>

void cfac_dummy_dcoul(double a1, double a2, int a3, double a4, double *a5,
		 double *a6, double *a7, double *a8, int *a9)
{
    fputs("This version of cFAC was compiled without COULCC module\n", stderr);
    abort();
}

void cfac_dummy_cmultip(double a1, double a2, double a3, double a4, double a5,
		  int a6, int a7, int a8, double *a9, int *a10)
{
    fputs("This version of cFAC was compiled without CMULTIP module\n", stderr);
    abort();
}
