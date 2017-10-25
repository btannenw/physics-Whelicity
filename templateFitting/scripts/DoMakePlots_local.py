#! /usr/bin/python

import subprocess,sys
import os,string
import os.path

fPath  = os.getcwd()+"/"+sys.argv[1]

#fFolderList = []
#fFolderList.append("Output_Sample0")

#Determine number of nuisance parameters
nui_num = 1
nui_num_str= str(nui_num)
print 'Number of nuisance parameters : '+nui_num_str

#Run DoMakePlots

for n in range(0,1):
    fFolderList = []
    n_str=str(n)
    fFolder = "Output_Sample_mu_"+n_str
    cmd_run = os.getcwd()+"/runProducePlots;"+fFolder+";mu;"+n_str+";"+nui_num_str+";mode;"+fFolder+"_Plots"
    print cmd_run.split(";")
    subprocess.call(cmd_run.split(";"), cwd=fPath)
  
