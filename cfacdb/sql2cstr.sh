#!/bin/sh

sqlfile=$1

strname=`basename $1 .sql`

echo "static const char *${strname}[] = {"
cat $sqlfile | tr -d \\n | tr \; \\n |
sed '
s/\ \ */\ /g
s/^/"/
s/$/",/
'

echo "NULL};"
