#!/bin/bash
#
# Rerun a phantom simulation
#
read -p 'Filename: ' jobname

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
