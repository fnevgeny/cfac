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

#ifndef _STOKEN_H_
#define _STOKEN_H_ 1

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "array.h"
#include "parser.h"

#define MAXNARGS 128
#define MAXLINELENGTH 4000
#define MAXMETHODNAME 64

#define ERR_LINEUNTERMINATED (-1)
#define ERR_LINETOOLONG      (-2)
#define ERR_NOVARIABLE       (-3)
#define ERR_ARGSTOOMANY      (-4)
#define ERR_SYNTAX           (-5)
#define ERR_EVAL             (-6)

#define NUMBER  (0)
#define STRING  (1)
#define TUPLE   (2)
#define LIST    (3)
#define KEYWORD (4)

#define COMMENT       ('#')
#define CONTINUE      ('\\')

typedef struct _METHOD_ {
  char name[MAXMETHODNAME];
  int (*func)(int argc, char *argv[], int argt[], ARRAY *variables);
} METHOD;

typedef struct _VARIABLE_ {
  char *name;
  char *value;
  int type;
} VARIABLE;

typedef struct _STATEMENT_ {
  int nline;
  int imethod;
  int argc;
  char *argv[MAXNARGS];
  int argt[MAXNARGS];
} STATEMENT;

VARIABLE *VariableExists(char *name, ARRAY *variables);
int DecodeArgs(char *s, char *argv[], int argt[], ARRAY *variables);
int GetLine(FILE *f, char *line, int *nlines);
int GetValidLine(FILE *f, char *line, int *nlines);
int MethodIndex(char *name, METHOD *methods);
int TokenizeLine(int nline, char *line, METHOD *methods,
                 ARRAY *statements, ARRAY *variables);
int EvalFile(FILE *f, int exebyline, METHOD *methods, FILE *ferr);
int EvalStatement(STATEMENT *st, METHOD *methods, ARRAY *variables);
void ErrorOcurred(int ierr, int loc, FILE *ferr);

#endif
