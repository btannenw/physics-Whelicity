#!/usr/bin/python

import os, subprocess, sys

systMin = int(sys.argv[1])
systMax = int(sys.argv[2])

def submitLine ( region, syst ):
    """something about something"""
    os.system("cp batchsubmit_{0}.sh batchsubmit_temp.sh".format(region))
    os.system("sed -i 's/$1/{0}/g' batchsubmit_temp.sh".format(syst))
    os.system("bsub -q 8nh -J {0}_{1} < ./batchsubmit_temp.sh".format(region, syst))

i=systMin
while i < systMax:
    submitLine("lep",i)
    submitLine("had",i)
    submitLine("lephad",i)

    i=i+1
