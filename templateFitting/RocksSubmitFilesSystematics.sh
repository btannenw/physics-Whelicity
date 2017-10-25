#!/bin/bash

for ((i = 0; i < 35; i++)) 
do 
 syst=Nominal
 if [ $i = 0 ]; then syst=jes
 fi
 if [ $i = 1 ]; then syst=jer
 fi
 if [ $i = 2 ]; then syst=jeff
 fi
 if [ $i = 3 ]; then syst=cellout
 fi
 if [ $i = 4 ]; then syst=softjet
 fi
 if [ $i = 5 ]; then syst=pileup
 fi
 if [ $i = 6 ]; then syst=Shower
 fi
 if [ $i = 7 ]; then syst=MCGen1
 fi
 if [ $i = 8 ]; then syst=MCGen2
 fi
 if [ $i = 9 ]; then syst=ISR_FSR
 fi
 if [ $i = 10 ]; then syst=CR_Perugia
 fi
 if [ $i = 11 ]; then syst=muid
 fi 
 if [ $i = 12 ]; then syst=mums
 fi 
 if [ $i = 13 ]; then syst=btag
 fi 
 if [ $i = 14 ]; then syst=ctag
 fi
 if [ $i = 15 ]; then syst=mistag
 fi
 if [ $i = 16 ]; then syst=musc
 fi
 if [ $i = 17 ]; then syst=muonid
 fi
 if [ $i = 18 ]; then syst=muonreco
 fi
 if [ $i = 19 ]; then syst=muontrig
 fi
 if [ $i = 20 ]; then syst=jvfsf
 fi
 if [ $i = 21 ]; then syst=UE
 fi
 if [ $i = 22 ]; then syst=TopMass
 fi
 if [ $i = 23 ]; then syst=Calibration
 fi
 if [ $i = 24 ]; then syst=Datafit
 fi
 if [ $i = 25 ]; then syst=mu_QCDshape
 fi
 if [ $i = 26 ]; then syst=wjetsptjmin10
 fi
 if [ $i = 27 ]; then syst=wjets_iqop3
 fi
 if [ $i = 28 ]; then syst=wjets_bb4
 fi
 if [ $i = 29 ]; then syst=wjets_bb5
 fi
 if [ $i = 30 ]; then syst=wjets_c4
 fi
 if [ $i = 31 ]; then syst=wjets_c5
 fi
 if [ $i = 32 ]; then syst=wjets_bbcc
 fi
 if [ $i = 33 ]; then syst=wjets_bbccc
 fi
 if [ $i = 34 ]; then syst=TemplateStat
 fi
 echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate mu 5000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 2
done
 

for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 35; i++)) 
do 
 syst=Nominal
 if [ $i = 0 ]; then syst=jes
 fi
 if [ $i = 1 ]; then syst=jer
 fi
 if [ $i = 2 ]; then syst=jeff
 fi
 if [ $i = 3 ]; then syst=cellout
 fi
 if [ $i = 4 ]; then syst=softjet
 fi
 if [ $i = 5 ]; then syst=pileup
 fi
 if [ $i = 6 ]; then syst=Shower
 fi
 if [ $i = 7 ]; then syst=MCGen1
 fi
 if [ $i = 8 ]; then syst=MCGen2
 fi
 if [ $i = 9 ]; then syst=ISR_FSR
 fi
 if [ $i = 10 ]; then syst=CR_Perugia
 fi
 if [ $i = 11 ]; then syst=eer
 fi 
 if [ $i = 12 ]; then syst=ees
 fi 
 if [ $i = 13 ]; then syst=btag
 fi 
 if [ $i = 14 ]; then syst=ctag
 fi
 if [ $i = 15 ]; then syst=mistag
 fi
 if [ $i = 16 ]; then syst=electronid
 fi
 if [ $i = 17 ]; then syst=electronreco
 fi
 if [ $i = 18 ]; then syst=electrontrig
 fi
 if [ $i = 19 ]; then syst=jvfsf
 fi
 if [ $i = 20 ]; then syst=UE
 fi
 if [ $i = 21 ]; then syst=TopMass
 fi
 if [ $i = 22 ]; then syst=Calibration
 fi
 if [ $i = 23 ]; then syst=Datafit
 fi
 if [ $i = 24 ]; then syst=el_QCD_FAKE
 fi
 if [ $i = 25 ]; then syst=el_QCD_REAL
 fi
 if [ $i = 26 ]; then syst=wjetsptjmin10
 fi
 if [ $i = 27 ]; then syst=wjets_iqop3
 fi
 if [ $i = 28 ]; then syst=wjets_bb4
 fi
 if [ $i = 29 ]; then syst=wjets_bb5
 fi
 if [ $i = 30 ]; then syst=wjets_c4
 fi
 if [ $i = 31 ]; then syst=wjets_c5
 fi
 if [ $i = 32 ]; then syst=wjets_bbcc
 fi
 if [ $i = 33 ]; then syst=wjets_bbccc
 fi
 if [ $i = 34 ]; then syst=TemplateStat
 fi

  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el 5000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 2
done
 
for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 41; i++)) 
do 
 syst=Nominal
 if [ $i = 0 ]; then syst=jes
 fi
 if [ $i = 1 ]; then syst=jer
 fi
 if [ $i = 2 ]; then syst=jeff
 fi
 if [ $i = 3 ]; then syst=cellout
 fi
 if [ $i = 4 ]; then syst=softjet
 fi
 if [ $i = 5 ]; then syst=pileup
 fi
 if [ $i = 6 ]; then syst=Shower
 fi
 if [ $i = 7 ]; then syst=MCGen1
 fi
 if [ $i = 8 ]; then syst=MCGen2
 fi
 if [ $i = 9 ]; then syst=ISR_FSR
 fi
 if [ $i = 10 ]; then syst=CR_Perugia
 fi
 if [ $i = 11 ]; then syst=eer
 fi 
 if [ $i = 12 ]; then syst=ees
 fi 
 if [ $i = 13 ]; then syst=btag
 fi 
 if [ $i = 14 ]; then syst=ctag
 fi
 if [ $i = 15 ]; then syst=mistag
 fi
 if [ $i = 16 ]; then syst=electronid
 fi
 if [ $i = 17 ]; then syst=electronreco
 fi
 if [ $i = 18 ]; then syst=electrontrig
 fi
 if [ $i = 19 ]; then syst=jvfsf
 fi
 if [ $i = 20 ]; then syst=UE 
 fi    
 if [ $i = 21 ]; then syst=muid
 fi 
 if [ $i = 22 ]; then syst=mums
 fi 
 if [ $i = 23 ]; then syst=musc
 fi
 if [ $i = 24 ]; then syst=muonid
 fi
 if [ $i = 25 ]; then syst=muonreco
 fi
 if [ $i = 26 ]; then syst=muontrig
 fi
 if [ $i = 27 ]; then syst=TopMass
 fi
 if [ $i = 28 ]; then syst=Calibration
 fi
 if [ $i = 29 ]; then syst=Datafit
 fi
 if [ $i = 30 ]; then syst=mu_QCDshape
 fi
 if [ $i = 31 ]; then syst=el_QCD_FAKE
 fi
 if [ $i = 32 ]; then syst=el_QCD_REAL
 fi
 if [ $i = 33 ]; then syst=wjetsptjmin10
 fi
 if [ $i = 34 ]; then syst=wjets_iqop3
 fi
 if [ $i = 35 ]; then syst=wjets_bb4
 fi
 if [ $i = 36 ]; then syst=wjets_bb5
 fi
 if [ $i = 37 ]; then syst=wjets_c4
 fi
 if [ $i = 38 ]; then syst=wjets_c5
 fi
 if [ $i = 39 ]; then syst=wjets_bbcc
 fi
 if [ $i = 40 ]; then syst=wjets_bbccc
 fi

  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el_mu 5000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 2
done
 

for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 53; i++))
do
 syst=PDF_CT10_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 2
done

for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 41; i++))
do
 syst=PDF_MSTW_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd`

  sleep 2
done

for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done

for ((i = 0; i < 100; i++))
do
 syst=PDF_NNPDF_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 2
done

for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 53; i++))
do
 syst=PDF_CT10_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate mu 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 2
done


for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 41; i++))
do
 syst=PDF_MSTW_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate mu 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd`

  sleep 2
done


for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 100; i++))
do
 syst=PDF_NNPDF_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate mu 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 1
done


for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done


for ((i = 0; i < 53; i++))
do
 syst=PDF_CT10_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el_mu 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 1
done

for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done

for ((i = 0; i < 41; i++))
do
 syst=PDF_MSTW_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el_mu 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd`

  sleep 1
done

for ((i = 0; i < 35; i++))
do
  echo $i
  sleep 5
done

for ((i = 0; i < 100; i++))
do
 syst=PDF_NNPDF_$i
 echo $syst
  echo "#!/bin/bash" > submit_profiling_${i}.sh
 echo 'export PROFILECODE=/home/aknue/GoeProfilingPackage_26012013/' >> submit_profiling_${i}.sh
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
 echo './runEvaluateExternalUpdate el_mu 2000 ' ${syst} ' 3D' >> submit_profiling_${i}.sh
 echo '#--------------------------' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/.' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.txt $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalSystematicsOutput_el_mu_3D/*.html $PROFILECODE/Results/ExternalSystematicsOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.eps $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cp ExternalCalibrationOutput_el_mu_3D/*.root $PROFILECODE/Results/ExternalCalibrationOutput_el_mu_3D/' >> submit_profiling_${i}.sh
 echo 'cd ../' >> submit_profiling_${i}.sh
 echo 'rm -rf $TMPDIR/$PBS_JOBID' >> submit_profiling_${i}.sh
  cmd="qsub -q shorttime submit_profiling_${i}.sh"
  `echo $cmd` 

  sleep 1
done

