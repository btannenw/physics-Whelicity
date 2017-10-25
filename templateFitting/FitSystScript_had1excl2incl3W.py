#!/usr/bin/python

import os, subprocess, sys

systMin = int(sys.argv[1])
systMax = int(sys.argv[2])

# Set general tags
elMu='el_mu_bTag'
nPseudo='5000'
nDim='3D'
newSample='true'
calibMode='single'
btag='1excl2incl'
lephad='had'
fixfloat='float'
wBkg='3W'


#Loop over sets in systList.sh and submit job for each
systList = open ('systConfig_Sig.txt','r')
readVars=[]

i=0
for systVar in systList:
    
    if '#' not in systVar:
        if i >= systMin and i <= systMax:
            systVar = systVar.split(',')[0]
                            
            #print systVar.split('\n')[0], i
            if systVar not in readVars:
                readVars.append(systVar)
                #print "./runEvaluateExternalUpdate {0} {1} {2} {3} {4} {5} {6}".format(elMu, nPseudo, systVar.split('\n')[0], nDim, newSample, calibMode, btag)
                os.system("./runEvaluateExternalUpdate {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(elMu, nPseudo, systVar.split('\n')[0], nDim, newSample, calibMode, btag, lephad, fixfloat, wBkg))

        #iterate counter
        i=i+1

#./runEvaluateExternalUpdate el 5000 jes_EffectiveNP_Detector1 3D false single 2incl
os.system('rm temp.root')
