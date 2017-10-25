
#!/bin/bash

syst=279
echo $syst
export TTBAR=/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/
cd $TTBAR

alias src="source rcSetup.sh"

$src
cd FitPackage_allinOne_Apr26_allCF

python FitSystScript_lephad1excl2incl3W_single.py 279 279
#python FitSystScript_lephad2incl3W_single.py 279 279
