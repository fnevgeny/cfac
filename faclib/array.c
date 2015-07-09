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
  Implementation of module "array"
  
  This module implements a variable length one- and 
  multi-dimensional array.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

#include <stdlib.h>
#include <string.h>

#include "array.h"

void InitIntData(void *p, int n) {
  memset(p, 0, sizeof(int)*n);
}

void InitDoubleData(void *p, int n) {
  memset(p, 0, sizeof(double)*n);
}

void InitPointerData(void *p, int n) {
  memset(p, 0, sizeof(void *)*n);
}

void InitArrayData(void *p, int n) {
  memset(p, 0, sizeof(ARRAY)*n);
}

/* 
** FUNCTION:    ArrayInit
** PURPOSE:     initialize the one-dimensional array.
** INPUT:       {ARRAY *a},
**              pointer to the array to be initialized.
**              {int esize},
**              size of the elements in bytes.
**              {int block},
**              number of elements in one block.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/
int ArrayInit(ARRAY *a, int esize, int block,
    ARRAY_ELEM_FREE FreeElem, ARRAY_DATA_INIT InitData) {
  a->esize = esize;
  a->block = block;
  a->bsize = ((int)esize)*((int)block);
  a->dim = 0;
  a->data = NULL;
  a->FreeElem = FreeElem;
  a->InitData = InitData;
  return 0;
}

/* 
** FUNCTION:    ArrayGet
** PURPOSE:     retrieve the i-th element of the array.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {int i},
**              index of the element.
** RETURN:      {void *},
**              pointer to the element. 
**              NULL, if does not exist.
** SIDE EFFECT: 
** NOTE:        
*/
void *ArrayGet(ARRAY *a, int i) {
  DATA *p;
  if (i < 0 || i >= a->dim) return NULL;
  p = a->data;
  while (i >= a->block) {
    p = p->next;
    i -= a->block;
  }
  if (p->dptr) {
    return ((char *) p->dptr) + i*(a->esize);
  } else {
    return NULL;
  }
}

/* 
** FUNCTION:    ArraySet
** PURPOSE:     set the i-th element.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {int i},
**              index of the element.
**              {void *d},
**              pointer to the data to be copied.
**              {void (*InitData)(void *, int)},
**              a function to be called to initialize the data
**              when first created.
** RETURN:      {void *},
**              pointer to the element.
** SIDE EFFECT: 
** NOTE:        if d == NULL, this function simply retrieve the
**              i-th element. if the element does not exist,
**              an empty one is created.
*/
void *ArraySet(ARRAY *a, int i, const void *d) {
  void *pt;
  char *ct;
  DATA *p;
 
  if (a->dim == 0) {
    a->data = malloc(sizeof(DATA));
    a->data->dptr = malloc(a->bsize);
    if (a->InitData) a->InitData(a->data->dptr, a->block);
    a->data->next = NULL;
  }
  p = a->data;
  if (a->dim <= i) a->dim = i+1;
  while (i >= a->block) {
    if (!(p->next)) {
      p->next = malloc(sizeof(DATA));
      p->next->dptr = NULL;
      p->next->next = NULL;
    }
    p = p->next;
    i -= a->block;
  }

  if (!(p->dptr)) {
    p->dptr = malloc(a->bsize);
    if (a->InitData) a->InitData(p->dptr, a->block);
  }
  
  ct = p->dptr;
  for (; i > 0; i--) {
    ct += a->esize;
  }
  pt = ct;

  if (d) memcpy(pt, d, a->esize);
  return pt;
}

/* 
** FUNCTION:    ArrayAppend
** PURPOSE:     append an element to the array
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {void *d},
**              data to be appened.
** RETURN:      {void *},
**              pointer to the appended element.
** SIDE EFFECT: 
** NOTE:        
*/
void *ArrayAppend(ARRAY *a, const void *d) {
  int i;
  if (!a) {
    return NULL;
  }
  i = a->dim;
  return ArraySet(a, i, d);
}

/* 
** FUNCTION:    ArrayFreeData
** PURPOSE:     free the data stored in the array.
** INPUT:       {DATA *p},
**              pointer to the data to be freed
**              {int esize},
**              size of the element in bytes.
**              {int block},
**              number of elements in one block.
**              {void (*FreeElem)(void *)},
**              a function called before freeing the data.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        this function calls itself recursively.
*/
int ArrayFreeData(ARRAY *a, DATA *p) {
  char *pt;
  int i;

  if (p == NULL) return 0;

  if (p->next) {
    ArrayFreeData(a, p->next);
  }
    
  if (a->FreeElem && p->dptr) {
    pt = p->dptr;
    for (i = 0; i < a->block; i++) {
      a->FreeElem(pt);
      pt += a->esize;
    }
  }
  if (p->dptr) {
    free(p->dptr);
  }
  free(p);

  return 0;
}

/* 
** FUNCTION:    ArrayFree
** PURPOSE:     deinitialize the array.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {void (*FreeElem)(void *)},
**              a function called before freeing each element.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/    
int ArrayFree(ARRAY *a) {
  if (!a || !a->dim) return 0;
  ArrayFreeData(a, a->data);
  a->dim = 0;
  a->data = NULL;
  return 0;
}

/* 
** FUNCTION:    ArrayTrim
** PURPOSE:     Trim the tail of an array to a given length.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {int n},
**              length of the final array.
**              {void (*FreeElem)(void *)},
**              a function called before freeing each element.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        if the length of array is <= n, nothing happens.
*/    
int ArrayTrim(ARRAY *a, int n) {
  DATA *p;
  void *pt;
  int i;

  if (!a) return 0;
  if (a->dim <= n) return 0;
  
  if (n == 0) {
    ArrayFree(a);
    return 0;
  }

  i = n;
  p = a->data;
  while (i >= a->block) {
    p = p->next;
    i -= a->block;
  }

  if (i == 0) {
    ArrayFreeData(a, p);
    p = NULL;
  } else {
    if (p->next) {
      ArrayFreeData(a, p->next);
      p->next = NULL;
    }
    if (p->dptr && a->FreeElem) {
      pt = ((char *) p->dptr) + i*(a->esize);
      for (; i < a->block; i++) {
	a->FreeElem(pt);
	pt = ((char *) pt) + a->esize;
      }
    }
  }

  a->dim = n;

  return 0;
}         

/* hash table size 2^NHASH */
typedef unsigned long int ub4;
typedef struct _MDATA_ {
  int *index;
  void *data;
} MDATA;

void InitMDataData(void *p, int n) {
  MDATA *d;
  int i;
  
  d = (MDATA *) p;
  for (i = 0; i < n; i++) {
    d[i].index = NULL;
    d[i].data = NULL;
  }
}

#define HashSize(n) ((ub4)1<<(((n)/2)+16))
#define HashMask(n) (HashSize(n)-1)
#define Mix(a, b, c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

static int Hash2(int *id, ub4 length, ub4 initval, int n) {
  register ub4 a, b, c, len, *k;
  ub4 kd[32], i;

  k = kd;
  for (i = 0; i < length; i++) k[i] = id[i];

  len = length;  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;           /* the previous hash value */
  
  /*---------------------------------------- handle most of the key */
  while (len >= 3) {
    a += k[0];
    b += k[1];
    c += k[2];
    Mix(a,b,c);
    k += 3; len -= 3;
  }

  /*-------------------------------------- handle the last 2 ub4's */
  c += length;
  switch(len) {
    /* c is reserved for the length */
  case 2 : b+=k[1];
  case 1 : a+=k[0];
    /* case 0: nothing left to add */
  }
  Mix(a,b,c);
  /*-------------------------------------------- report the result */
  return (int) (c & HashMask(n));
}

int NMultiInit(MULTI *ma, int esize, int ndim, int *block,
    ARRAY_ELEM_FREE FreeElem, ARRAY_DATA_INIT InitData) {
  int i, n;

  ma->maxelem = -1;
  ma->numelem = 0;
  ma->ndim = ndim;
  ma->isize = sizeof(int)*ndim;
  ma->esize = esize;
  ma->block = (unsigned short *) malloc(sizeof(unsigned short)*ndim);
  n = HashSize(ma->ndim);

  ma->array = (ARRAY *) malloc(sizeof(ARRAY)*n);
  for (i = 0; i < n; i++) {
    ArrayInit(&(ma->array[i]), sizeof(MDATA), 10, FreeElem, InitData);
  }

  return 0;
}

void *NMultiGet(MULTI *ma, int *k) {
  ARRAY *a;
  MDATA *pt;
  DATA *p;
  int i, j, m, h;

  h = Hash2(k, ma->ndim, 0, ma->ndim);
  a = &(ma->array[h]);
  p = a->data;
  i = a->dim;
  j = 0;
  while (p) {
    pt = (MDATA *) p->dptr;
    for (m = 0; m < a->block && j < i; j++, m++) {
      if (memcmp(pt->index, k, ma->isize) == 0) {
	return pt->data;
      }
      pt++;
    }
    p = p->next;
  }

  return NULL;
}

void *NMultiSet(MULTI *ma, int *k, void *d) {
  int i, j, m = 0, h;
  MDATA *pt = NULL;
  ARRAY *a;
  DATA *p, *p0 = NULL;

  if (ma->maxelem > 0 && ma->numelem >= ma->maxelem) {
    NMultiFreeData(ma);
    ma->numelem = 0;
  }
  h = Hash2(k, ma->ndim, 0, ma->ndim);
  a = &(ma->array[h]);
  if (a->dim == 0) {
    a->data = (DATA *) malloc(sizeof(DATA));
    a->data->dptr = malloc(a->bsize);
    InitMDataData(a->data->dptr, a->block);
    a->data->next = NULL;
    pt = (MDATA *) a->data->dptr;
  } else {
    p = a->data;
    i = a->dim;
    j = 0;
    while (p) {
      pt = (MDATA *) p->dptr;
      for (m = 0; m < a->block && j < i; j++, m++) {
	if (memcmp(pt->index, k, ma->isize) == 0) {
	  if (d) {
	    memcpy(pt->data, d, ma->esize);
	  }
	  return pt->data;
	}
	pt++;
      }
      p0 = p;
      p = p->next;
    }  
    if (m == a->block) {
      p0->next = (DATA *) malloc(sizeof(DATA));
      p = p0->next;
      p->dptr = malloc(a->bsize);
      InitMDataData(p->dptr, a->block);
      p->next = NULL;
      pt = (MDATA *) p->dptr;
    }
  }

  ma->numelem++;
  pt->index = malloc(ma->isize);
  memcpy(pt->index, k, ma->isize);
  pt->data = malloc(ma->esize);
  if (a && a->InitData) a->InitData(pt->data, 1);
  if (d) memcpy(pt->data, d, ma->esize);
  (a->dim)++;

  return pt->data;
}

static int NMultiArrayFreeData(DATA *p, int esize, int block, 
			       ARRAY_ELEM_FREE FreeElem) { 
  MDATA *pt;
  int i;
  
  if (p->next) {
    NMultiArrayFreeData(p->next, esize, block, FreeElem);
  }

  if (p->dptr) {
    pt = p->dptr;
    for (i = 0; i < block; i++) {
      free(pt->index);
      if (FreeElem && pt->data) FreeElem(pt->data);
      free(pt->data);
      pt++;
    }
    free(p->dptr);
  }
  if (p) {
    free(p);
  }
  p = NULL;
  return 0;
}
    
int NMultiFreeDataOnly(ARRAY *a) {
  if (!a) return 0;
  if (a->dim == 0) return 0;
  NMultiArrayFreeData(a->data, a->esize, a->block, a->FreeElem);
  a->dim = 0;
  a->data = NULL;
  return 0;
}

int NMultiFreeData(MULTI *ma) {
  ARRAY *a;
  int i, n;

  if (!ma) return 0;
  
  n = HashSize(ma->ndim);
  for (i = 0; i < n; i++) {
    a = &(ma->array[i]);
    NMultiFreeDataOnly(a);
  }
  return 0;
}

int NMultiFree(MULTI *ma) {
  if (!ma) return 0;
  if (ma->ndim <= 0) return 0;
  NMultiFreeData(ma);
  free(ma->array);
  ma->array = NULL;
  free(ma->block);
  ma->block = NULL;
  ma->ndim = 0;
  return 0;
}
