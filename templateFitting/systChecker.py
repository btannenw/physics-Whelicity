#!/usr/bin/python

import os, sys

os.system("cat ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_lep/SystematicOutput_el_mu_bTag.txt | sort > sort_lep.txt")
os.system("cat ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had/SystematicOutput_el_mu_bTag.txt | sort > sort_had.txt")
os.system("cat ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/SystematicOutput_el_mu_lephad_bTag.txt | sort > sort_lephad.txt")

f_lep=open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr13/sort_lep.txt','r')
f_had=open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr13/sort_had.txt','r')
f_lephad=open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr13/sort_lephad.txt','r')

lep=[]
had=[]
lephad=[]

for line in f_lep:
    #print line.split('\n')[0].split('\t')[0]
    variation = line.split('\n')[0].split('\t')[0]
    if variation not in lep:
        lep.append(variation)

for line in f_had:
    #print line.split('\n')[0].split('\t')[0]
    variation = line.split('\n')[0].split('\t')[0]
    if variation not in had:
        had.append(variation)
    else:
        print variation,"extra, had"

for line in f_lephad:
    #print line.split('\n')[0].split('\t')[0]
    variation = line.split('\n')[0].split('\t')[0]
    if variation not in lephad:
        lephad.append(variation)
    else:
        print variation,"extra, lephad"


print len(lep),"in lep file"
print len(had),"in had file"
print len(lephad),"in lephad file"

for item in lep:
    if item not in had:
        print item,"not in had!"
    if item not in lephad:
        print item,"not in lephad!"


for item in lephad:
    if item not in lep:
        print item,"not in lep!"
