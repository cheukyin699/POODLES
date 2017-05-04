#!/bin/bash

#usage execPoodleS.sh LIST FASTA_DIR

INSTALLDIR=/data/share/GenSET_root/programs/attribute_programs/POODLE-S

cat $1 | while read list

do
mkdir $INSTALLDIR/usrData/$list
./safeFasta.perl $INSTALLDIR/$2/$list.fst > $INSTALLDIR/usrData/$list/query.fst
cp $INSTALLDIR/execPoodleS_SRC $INSTALLDIR/usrData/$list/exec.sh
chmod 755 $INSTALLDIR/usrData/$list/exec.sh

cd $INSTALLDIR/usrData/$list
./exec.sh query.fst 0
cat ./caspFmt/psitmp.txt > $INSTALLDIR/results/$list
cd ../../
rm -rf $INSTALLDIR/usrData/$list
done
