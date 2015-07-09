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

#ifndef _ARRAY_H_
#define _ARRAY_H_ 1

/*************************************************************
  Header of module "array"
  
  This module implements a variable length one- and 
  multi-dimensional array.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

/* 
<** The following format is used for documenting the source **>
*/

/* documenting a struct */
/*
** STRUCT:      
** PURPOSE:     
** FIELDS:      
** NOTE:        
*/

/* documenting a function */
/* 
** FUNCTION:    
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a macro function */
/* 
** MACRO:       
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a global, static varialbe or a macro constant */
/*
** VARIABLE:    
** TYPE:        
** PURPOSE:     
** NOTE:        
*/

#define MultiInit NMultiInit
#define MultiGet NMultiGet
#define MultiSet NMultiSet
#define MultiFreeData NMultiFreeData
#define MultiFree NMultiFree


/*
** STRUCT:      DATA
** PURPOSE:     data structure of the array elements
** FIELDS:      {void *dptr},
**              pointer to a block of data.
**              {DATA *next},
**              pointer to the next block.
** NOTE:        
*/
typedef struct _DATA_ {
  void *dptr;
  struct _DATA_ *next;
} DATA;

/*
** STRUCT:      ARRAY
** PURPOSE:     a one-dimensional array.
** FIELDS:      {short esize},
**              the size of each element in bytes.
**              {short block},
**              number of elements in each block.
**              {int dim},
**              the size of the array.
** NOTE:        
*/

typedef struct _ARRAY_ ARRAY;

typedef void (*ARRAY_ELEM_FREE)(void *elem);
typedef void (*ARRAY_DATA_INIT)(void *elem, int);

struct _ARRAY_ {
  unsigned short esize;
  unsigned short block;
  int bsize;
  int dim;
  DATA  *data;
  
  ARRAY_ELEM_FREE FreeElem;
  ARRAY_DATA_INIT InitData;
};

/*
** STRUCT:      MULTI
** PURPOSE:     a multi-dimensional array.
** FIELDS:      {short ndim},
**              the dimension of the array.
**              {short esize},
**              size of each array element in bytes.
**              {short *block},
**              number of elements in each block for each dimension.
**              {ARRAY *array},
**              the multi-dimensional array is implemented as array 
**              of arrays. 
** NOTE:        
*/
typedef struct _MULTI_ {
  int numelem, maxelem;
  unsigned short ndim;
  unsigned short isize;
  unsigned short esize;
  unsigned short *block;
  ARRAY *array;
} MULTI;

int   ArrayInit(ARRAY *a, int esize, int block,
    ARRAY_ELEM_FREE FreeElem, ARRAY_DATA_INIT InitData);
void *ArrayGet(ARRAY *a, int i);
void *ArraySet(ARRAY *a, int i, const void *d);
void *ArrayAppend(ARRAY *a, const void *d);
int   ArrayTrim(ARRAY *a, int n);
int   ArrayFree(ARRAY *a);
int   ArrayFreeData(ARRAY *a, DATA *p);

/*
** the following set of funcitons are a different implementation
** for the MULTI array,
*/
int   NMultiInit(MULTI *ma, int esize, int ndim, int *block,
    ARRAY_ELEM_FREE FreeElem, ARRAY_DATA_INIT InitData);
void *NMultiGet(MULTI *ma, int *k);
void *NMultiSet(MULTI *ma, int *k, void *d);
int   NMultiFree(MULTI *ma);
int   NMultiFreeDataOnly(ARRAY *a);
int   NMultiFreeData(MULTI *ma);

void  InitIntData(void *p, int n);
void  InitDoubleData(void *p, int n);
void  InitPointerData(void *p, int n);
void  InitArrayData(void *p, int n);
void  InitMDataData(void *p, int n);

#endif
