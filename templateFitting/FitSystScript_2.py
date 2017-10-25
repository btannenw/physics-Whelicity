#!/usr/bin/python

import os, subprocess, sys

systMin = int(sys.argv[1])
systMax = int(sys.argv[2])

# Set general tags
btag='2incl'
nDim='3D'
nPseudo='100'
elMu='el_mu_lephad'
calibMode='single'
newSample='true'

#Loop over sets in systList.sh and submit job for each
systList = open ('systConfigFull.txt','r')
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
                os.system("./runEvaluateExternalUpdate {0} {1} {2} {3} {4} {5} {6}".format(elMu, nPseudo, systVar.split('\n')[0], nDim, newSample, calibMode, btag))

        #iterate counter
        i=i+1

#./runEvaluateExternalUpdate el 5000 jes_EffectiveNP_Detector1 3D false single 2incl
os.system('rm temp.root')
