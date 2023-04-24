#!/bin/bash
#
# Rerun a phantom simulation
#
read -p 'Filename: ' jobname
length=${#jobname}

if [ $length == 0 ]; then
   exit 2
fi

rm p*
rm $jobname*
make cleanall

make; make setup

if [[ -f 'phantom' ]] || [[ -f 'phantomsetup' ]]; then
	./phantomsetup $jobname
	./phantomsetup $jobname
	if [[ -f $jobname'.in' ]]; then
	./phantom $jobname.in
	fi
fi
