#!/bin/bash
#
# Rerun a phantom simulation
#
read -p 'Filename: ' jobname
length_job=${#jobname}
if [ $length_job == 0 ]; then
   exit 2
fi

read -p 'Dumpfile: ' dumpfile
length_dump=${#dumpfile}
if [ $length_dump == 0 ]; then
   starting_dump=$jobname'_00000.tmp'
else
   starting_dump=$dumpfile
fi

rm phantom*
if [ $length_dump == 0 ]; then
   rm $jobname*
else
   rm $jobname.*
fi
make cleanall

make; make setup

if [[ -f 'phantom' ]] || [[ -f 'phantomsetup' ]]; then
	./phantomsetup $jobname
	./phantomsetup $jobname
	if [[ -f $jobname'.in' ]]; then
      starting_dump_default=$jobname'_00000.tmp'
      sed -i "s/$starting_dump_default/$starting_dump/g" $jobname'.in'
	   ./phantom $jobname.in
	fi
fi
