#!/bin/sh

echo "static const char *schema_str[] = {"
tr -d \\n\\r | tr \; \\n |
sed '
s/\ \ */\ /g
s/^/"/
s/$/",/
' "$@"

echo "NULL};"
