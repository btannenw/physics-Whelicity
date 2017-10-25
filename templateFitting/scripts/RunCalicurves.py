#! /usr/bin/python

import subprocess,sys
import os,string
import os.path

import DoPullDistributions

#variables to be defined by the user
fPath      = (os.getcwd()).replace("scripts", "")
#fRocksPath = "/rocks/aknue/GoeProfilingPackage_Philipp/"
fRocksPath = "/rocks/pstolte/project_promo1/trunk/"

fNuiVar    = 1         # Nuisance parameter that is varied
#fCaliType  = "CaliFraction" # CaliNui, FracNui, CaliFraction or CaliDatafit
#fCaliType  = "CaliDatafit" # CaliNui, FracNui, CaliFraction or CaliDatafit
#fCaliType  = "CaliNui" # CaliNui, FracNui, CaliFraction or CaliDatafit
fCaliType  = "FracNui" # CaliNui, FracNui, CaliFraction or CaliDatafit
fCaliBool  = True      # if False, make only pull distributions
fCaliPull  = True     # if False: already have pull distributions and just want to work on CaliCurves

#variables to be determined
fNuiNum = 0
is_mu   = False
is_el   = False

#copy all root files from rocks to pcatlas
if os.path.exists(fRocksPath):
    print fRocksPath," does exist!!!"
else:
    print "ERROR: ",fRocksPath," does NOT exist!!!"

print "Get all *.root files from ",fRocksPath
RootFileListMu = []
RootFileListEl = []
RootFileListMu = DoPullDistributions.GetRootFileList(fRocksPath, "mu")
RootFileListEl = DoPullDistributions.GetRootFileList(fRocksPath, "el")
print "Copy files from ",fRocksPath," to pcatlas"
for MuFile in RootFileListMu:
    cmd_cp         = "cp;"+MuFile+";"+fPath+"/"
    subprocess.call(cmd_cp.split(";"), cwd=fPath)
for ElFile in RootFileListEl:
    cmd_cp         = "cp;"+ElFile+";"+fPath+"/"
    subprocess.call(cmd_cp.split(";"), cwd=fPath)
            

#Determine number of nuisance parameters
print "Determine number of nuisance parameters"
for num in range(0,99):
    num_str=str(num)	
    if os.path.exists(fPath+"/Output_mu_F0_0.7_NrNui"+num_str+"_NuiVar"+str(fNuiVar)+"_p0_0.root"):
        fNuiNum = num
        is_mu   = True
    elif os.path.exists(fPath+"/Output_el_F0_0.7_NrNui"+num_str+"_NuiVar"+str(fNuiVar)+"_p0_0.root"):
        fNuiNum = num
        is_el   = True
    elif os.path.exists(fPath+"/DatafitPseudo_Output_mu_NrNui"+num_str+"_NuiVar"+str(fNuiVar)+"_p0_0.root"):
        fNuiNum = num
        is_mu   = True
    elif os.path.exists(fPath+"/DatafitPseudo_Output_el_NrNui"+num_str+"_NuiVar"+str(fNuiVar)+"_p0_0.root"):
        fNuiNum = num
        is_el   = True
    elif os.path.exists(fPath+"/OutputFluctuate_mu_F0_0.7_NrNui"+num_str+"_NuiVar"+str(fNuiVar)+"_0.root"):
        fNuiNum = num
        is_mu   = True
    elif os.path.exists(fPath+"/OutputFluctuate_el_F0_0.7_NrNui"+num_str+"_NuiVar"+str(fNuiVar)+"_0.root"):
        fNuiNum = num
        is_el   = True
nui_num_str= str(fNuiNum)




print 'Number of nuisance parameters : '+nui_num_str

if is_mu:
    fChannel = "mu"
    if fCaliType == "CaliNui":
        OutputFolder = "Output_Calicurve_F0_0.7_NrNui"+str(fNuiNum)+"_NuiVar"+str(fNuiVar)
        DoPullDistributions.MakeFolder(OutputFolder, fPath)
        OutputFolder = fPath+"/"+OutputFolder
        print "Make pull distributions for channel mu, fNui is varied"
        if fCaliPull:
            DoPullDistributions.MakePullNui(fNuiNum, fNuiVar, "mu", fPath)
        if fCaliBool:
            cmd_cali = fPath+"/runProduceCalibration;"+fPath+";"+OutputFolder+";"+fChannel+";CaliNui;"+str(fNuiNum)+";"+str(fNuiVar)
            subprocess.call(cmd_cali.split(";"), cwd=fPath)        
    elif fCaliType == "CaliFraction":
        OutputFolder = "Output_Calicurve_Fraction_NrNui"+str(fNuiNum)+"_NuiVar"+str(fNuiVar)
        if fCaliPull:
            DoPullDistributions.MakeFolder(OutputFolder, fPath)
        OutputFolder = fPath+"/"+OutputFolder                
        print "Make pull distributions for channel mu, F0 is varied"
        if fCaliPull:
            DoPullDistributions.MakePullFraction(fNuiNum, fNuiVar, "mu", fPath)
        if fCaliBool:
            cmd_cali = fPath+"/runProduceCalibration;"+fPath+";"+OutputFolder+";"+fChannel+";CaliFraction;"+str(fNuiNum)+";"+str(fNuiVar)
            print cmd_cali.split(";")
            subprocess.call(cmd_cali.split(";"), cwd=fPath)
    elif fCaliType == "CaliDatafit":
        OutputFolder = "DatafitPseudo_Output_Calicurve_Fraction_NrNui"+str(fNuiNum)+"_NuiVar"+str(fNuiVar)
        #if fCaliPull:
        #    DoPullDistributions.MakeFolder(OutputFolder, fPath)
        OutputFolder = fPath+"/"+OutputFolder                
        print "Make pull datafit distributions for channel mu"
        if fCaliPull:
            DoPullDistributions.MakePullDatafit(fNuiNum, fNuiVar, "mu", fPath)
        #if fCaliBool:
        #    cmd_cali = fPath+"/runProduceCalibration;"+fPath+";"+OutputFolder+";"+fChannel+";CaliDatafit;"+str(fNuiNum)+";"+str(fNuiVar)
        #    print cmd_cali.split(";")
        #    subprocess.call(cmd_cali.split(";"), cwd=fPath)
    elif fCaliType == "FracNui":
        OutputFolder = "OutputFluctuate_Fraction_NrNui"+str(fNuiNum)+"_NuiVar"+str(fNuiVar)
        #if fCaliPull:
        #    DoPullDistributions.MakeFolder(OutputFolder, fPath)
        OutputFolder = fPath+"/"+OutputFolder                
        print "Make fluctuate pull distributions for channel mu"
        if fCaliPull:
            DoPullDistributions.MakePullFluctuate(fNuiNum, fNuiVar, "mu", fPath)
    else:
        print "This is not a valid fCaliType !!!"







