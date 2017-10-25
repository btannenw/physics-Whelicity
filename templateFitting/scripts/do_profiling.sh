# user specific aliases and functions                                                                                                                               

#. /cvmfs/atlas.cern.ch/repo/sw/atlas-gcc/432/x86_64/setup.sh
#. /cvmfs/atlas.cern.ch/repo/sw/software/17.4.0/sw/lcg/app/releases/ROOT/5.30.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh

#export CURRY=`pwd`
#cd $CURRY

#./runGoeProfiling mu testinput CaliCurvesFraction BkgFit 3 1 0
#./runGoeProfiling mu Validation CaliCurvesNuisanceParam BkgFit 1 1 15 0
#./runGoeProfiling mu testinput CaliCurvesNuisanceParam BkgFit 3 1 0
#./runGoeProfiling mu testinput FluctuateNuisanceParam BkgFit 3 1 0
#./runGoeProfiling el testinput Datafit BkgFit 1 0 0
#./runGoeProfiling mu testinput Datafit BkgNui 1 1 0
#./runGoeProfiling mu testinput Datafit_Single BkgNui 1 1 0
./runGoeProfiling mu testinput Datafit_Single BkgFit 1 1 0
