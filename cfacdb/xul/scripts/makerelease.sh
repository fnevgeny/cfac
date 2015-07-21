#!/bin/sh

XPI=cfacdb.xpi
OUTDIR=releases

BASE_URL=http://plasma-gate.weizmann.ac.il/cfacdb/xul

KEYFILE=keys/keyfile.pem
UHURA=/usr/local/mxtools/uhura

VERSION=`grep em:version install.rdf | tr '<>' '\t' | cut -f3`

VERXPI=$OUTDIR/cfacdb-${VERSION}.xpi

cp $XPI $VERXPI
ln -sf cfacdb-${VERSION}.xpi $OUTDIR/$XPI

$UHURA -o $OUTDIR/update.rdf -k $KEYFILE $VERXPI ${BASE_URL}/$VERXPI
