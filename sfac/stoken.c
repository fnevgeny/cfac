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

#include "stoken.h"

VARIABLE *VariableExists(char *name, ARRAY *variables) {
  int i;

  for (i = 0; i < variables->dim; i++) {
    VARIABLE *v = ArrayGet(variables, i);
    if (strcmp(v->name, name) == 0) return v;
  }

  return NULL;
}

int DecodeArgs(char *s, char *argv[], int argt[], ARRAY *variables) {
  VARIABLE *v;
  int r, i;
  char token[MAXLINELENGTH];
  int brkpos;
  int next;
  int quotepos;
  int n;

  SetParserWhite(" \t");
  SetParserBreak(",=");
  SetParserQuote("\"'([", "\"')]");
  SetParserEscape('\0');

  next = 0;
  r = Parse(token, MAXLINELENGTH, s, &next, &brkpos, &quotepos);

  i = 0;
  while (1) {
    if (r > 0) break;
    if (r < 0) return ERR_SYNTAX;
    if (quotepos >= 0) {
      if (quotepos == 0 || quotepos == 1) {
        argt[i] = STRING;
      } else if (quotepos == 2) {
        argt[i] = TUPLE;
      } else if (quotepos == 3) {
        argt[i] = LIST;
      }
      n = strlen(token);
      argv[i] = (char *) malloc(n+1);
      strcpy(argv[i], token);
    } else {
      if (brkpos == 1) {
        argt[i] = KEYWORD;
        n = strlen(token);
        argv[i] = (char *) malloc(n+1);
        strcpy(argv[i], token);
      } else {
        if (token[0] == '$') {
          v = VariableExists(&(token[1]), variables);
          if (v == NULL) {
            return ERR_NOVARIABLE;
          }
          n = strlen(v->value);
          argv[i] = (char *) malloc(n+1);
          strcpy(argv[i], v->value);
          argt[i] = v->type;
        } else {
          n = strlen(token);
          if (n == 0) break;
          argv[i] = malloc(n+1);
          strcpy(argv[i], token);
          argt[i] = NUMBER;
        }
      }
    }
    i++;
    if (i == MAXNARGS) {
      return ERR_ARGSTOOMANY;
    }

    r = Parse(token, MAXLINELENGTH, s, &next, &brkpos, &quotepos);
  }

  return i;
}

int GetLine(FILE *f, char *line, int *nlines) {
  if (f == stdin) {
    fprintf(stdout, ">>> ");
  }
  if (fgets(line, MAXLINELENGTH, f) == NULL) return -1;
  (*nlines)++;
  return 0;
}

int GetValidLine(FILE *f, char *line, int *nlines) {
  int n, m;
  char buf[MAXLINELENGTH];
  char r;

  while (1) {
    if (GetLine(f, line, nlines) < 0) return 0;
    r = StrTrim(line, COMMENT);
    if (r == EOF || (line[0] && line[0] != COMMENT)) break;
  }
  n = strlen(line);
  if (n == 0) return n;
  if (line[n-1] == CONTINUE) {
    if (r == EOF) return ERR_LINEUNTERMINATED;
    line[n-1] = '\0';
    r = StrTrim(line, COMMENT);
    n = strlen(line);
    while (1) {
      while (1) {
        if (GetLine(f, buf, nlines) < 0) return 0;
        r = StrTrim(buf, COMMENT);
        if (r == EOF || (buf[0] && buf[0] != COMMENT)) break;
      }
      m = strlen(buf);
      if (m == 0) return ERR_LINEUNTERMINATED;
      n = n+m;
      if (n >= MAXLINELENGTH) return ERR_LINETOOLONG;
      strcat(line, buf);
      if (buf[m-1] != CONTINUE) break;
      line[n-1] = '\0';
      r = StrTrim(line, COMMENT);
      n = strlen(line);
      if (r == EOF) return ERR_LINEUNTERMINATED;
    }
  }
  return n;
}

int MethodIndex(char *name, METHOD *methods) {
  int i;

  i = 0;
  while (methods[i].func) {
    if (strncmp(methods[i].name, name, MAXMETHODNAME) == 0) {
      return i;
    }
    i++;
  }
  return -1;
}

int TokenizeLine(int nline, char *line, METHOD *methods,
                 ARRAY *statements, ARRAY *variables) {
  STATEMENT *fs;
  VARIABLE *v, *w;
  int i, r, n;
  char token[MAXLINELENGTH];
  char *t;
  int brkpos;
  int quotepos;
  int next;

  SetParserWhite(" \t");
  SetParserBreak(",=");
  SetParserQuote("\"'([", "\"')]");
  SetParserEscape('\0');

  next = 0;
  r = Parse(token, MAXLINELENGTH, line, &next, &brkpos, &quotepos);
  if (r) return ERR_SYNTAX;

  t = token;
  i = MethodIndex(t, methods);
  if (i >= 0) {
    fs = (STATEMENT *) ArraySet(statements, statements->dim, NULL);
    fs->nline = nline;
    fs->imethod = i;
    r = Parse(token, MAXLINELENGTH, line, &next, &brkpos, &quotepos);
    if (r) {
      fs->argc = 0;
    } else {
      if (quotepos != 2) return ERR_SYNTAX;
      fs->argc = DecodeArgs(token, fs->argv, fs->argt, variables);
      if (fs->argc < 0) return fs->argc;
    }
    return 1;
  } else {
    if (!isalpha(token[0])) return ERR_SYNTAX;
    if (quotepos >= 0) return ERR_SYNTAX;
    if (brkpos != 1) return ERR_SYNTAX;
    v = VariableExists(token, variables);
    if (v) {
      free(v->value);
    } else {
      v = (VARIABLE *) ArraySet(variables, variables->dim, NULL);
      n = strlen(token);
      v->name = (char *) malloc(n+1);
      strcpy(v->name, token);
    }

    SetParserBreak(",;");
    r = Parse(token, MAXLINELENGTH, line, &next, &brkpos, &quotepos);

    if (token[0] == '$') {
      w = VariableExists(&(token[1]), variables);
      if (w == NULL) {
        return ERR_NOVARIABLE;
      }
      n = strlen(w->value);
      v->value = (char *) malloc(n+1);
      strcpy(v->value, w->value);
      v->type = w->type;
    } else {
      n = strlen(token);
      v->value = (char *) malloc(n+1);
      strcpy(v->value, token);
      if (quotepos == 0 || quotepos == 1) {
        v->type = STRING;
      } else if (quotepos == 2) {
        v->type = TUPLE;
      } else if (quotepos == 3) {
        v->type = LIST;
      } else {
        v->type = NUMBER;
      }
    }
    return 0;
  }
}

static void InitStatementData(void *p, int n) {
  STATEMENT *st;
  int i;

  st = (STATEMENT *) p;
  for (i = 0; i < n; i++) {
    st[i].argc = 0;
  }
}

static void FreeStatementData(void *p) {
  STATEMENT *st;
  int i;

  st = (STATEMENT *) p;
  for (i = 0; i < st->argc; i++) {
    free(st->argv[i]);
  }
}

static void InitVariableData(void *p, int n) {
  VARIABLE *v;
  int i;

  v = (VARIABLE *) p;
  for (i = 0; i < n; i++) {
    v[i].name = NULL;
    v[i].value = NULL;
  }
}

static void FreeVariableData(void *p) {
  VARIABLE *v;

  v = (VARIABLE *) p;
  free(v->name);
  free(v->value);
}

int EvalFile(FILE *f, int exebyline, METHOD *methods, FILE *ferr) {
  ARRAY statements;
  ARRAY variables;
  STATEMENT *st;
  char buf[MAXLINELENGTH];
  int i, nlines;
  int ierr;

  ArrayInit(&statements,
    sizeof(STATEMENT), 1024, FreeStatementData, InitStatementData);
  ArrayInit(&variables,
    sizeof(VARIABLE), 1024, FreeVariableData, InitVariableData);

  nlines = 0;
  while (1) {
    i = GetValidLine(f, buf, &nlines);
    if (i == 0) break;
    if (i < 0) ErrorOcurred(i, nlines, ferr);
    i = TokenizeLine(nlines, buf, methods, &statements, &variables);
    if (i < 0) {
      ErrorOcurred(i, nlines, ferr);
      if (!exebyline) return -1;
    }
    if (exebyline && i > 0) {
      st = (STATEMENT *) ArrayGet(&statements, statements.dim-1);
      ierr = EvalStatement(st, methods, &variables);
      if (ierr < 0) {
        ErrorOcurred(ERR_EVAL, st->nline, ferr);
      }
    }
  }

  if (!exebyline) {
    for (i = 0; i < statements.dim; i++) {
      st = (STATEMENT *) ArrayGet(&statements, i);
      ierr = EvalStatement(st, methods, &variables);
      if (ierr < 0) {
        ErrorOcurred(ERR_EVAL, st->nline, ferr);
        return -1;
      }
    }
  }

  ArrayFree(&statements);
  ArrayFree(&variables);

  return 0;
}

int EvalStatement(STATEMENT *st, METHOD *methods, ARRAY *variables) {
  return methods[st->imethod].func(st->argc, st->argv, st->argt, variables);
}

void ErrorOcurred(int ierr, int loc, FILE *ferr) {
  fprintf(ferr, "Error at Line %d: ", loc);
  switch (ierr) {
  case ERR_LINEUNTERMINATED:
    fprintf(ferr, "Line Unterminated\n");
    break;
  case ERR_LINETOOLONG:
    fprintf(ferr, "Line Too Long, MaxLength: %d\n", MAXLINELENGTH);
    break;
  case ERR_NOVARIABLE:
    fprintf(ferr, "Variable does not exist\n");
    break;
  case ERR_ARGSTOOMANY:
    fprintf(ferr, "Arguments Too Many, Max: %d\n", MAXNARGS);
    break;
  case ERR_SYNTAX:
    fprintf(ferr, "Syntax Error\n");
    break;
  case ERR_EVAL:
    fprintf(ferr, "Evaluation Error\n");
    break;
  default:
    break;
  }
}

