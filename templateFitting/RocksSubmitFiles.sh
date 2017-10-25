#!/bin/bash

for ((i = 0; i < 100; i++)) 
do 
 echo "#!/bin/bash" > submit_profiling_${i}.sh
 #echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_Philipp/' >> submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/pstolte/project_promo1/trunk/' >> submit_profiling_${i}.sh
 echo 'mkdir -p $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
 echo 'cd $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp -r $PROFILECODE/* .' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo '. /cvmfs/atlas.cern.ch/repo/sw/atlas-gcc/432/x86_64/setup.sh' >> submit_profiling_${i}.sh
 echo '. /cvmfs/atlas.cern.ch/repo/sw/software/x86_64-slc5-gcc43-opt/17.4.0/sw/lcg/app/releases/ROOT/5.30.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh' >> submit_profiling_${i}.sh
 echo 'export BATINSTALLDIR=/home/aknue/GoeSelection_11_branch_0208/BAT-0.9.1/' >> submit_profiling_${i}.sh
 echo 'export KLFITTERINSTALLDIR=/home/aknue/GoeSelection_11_branch_0208/KLFitter-00-05-12/' >> submit_profiling_${i}.sh
 echo 'export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH:$BATINSTALLDIR/lib:$KLFITTERINSTALLDIR' >> submit_profiling_${i}.sh
 echo 'export CPATH=$CPATH:$BATINSTALLDIR/include' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'env' >> submit_profiling_${i}.sh
 echo 'ls -l' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 #echo './runGoeProfiling mu testinput CaliCurvesFraction BkgFit 10 0 ' $i >> submit_profiling_${i}.sh
 #echo './runGoeProfiling mu testinput CaliCurvesNuisanceParam BkgFit 10 1 ' $i >> submit_profiling_${i}.sh
 echo './runGoeProfiling mu testinput FluctuateNuisanceParam BkgFit 10 1 ' $i >> submit_profiling_${i}.sh
 #echo './runGoeProfiling mu testinput Datafit BkgFit 10 0 ' $i >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp *.root $PROFILECODE' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
#  cmd="qsub -q atlasS submit_profiling_${i}.sh"
  #echo $cmd 
  `echo $cmd` 

  sleep 1
done
 
