#! /usr/bin/python

import subprocess,sys
import os,string
import os.path

def MakeFolder(Folder, Path):
    print Folder
    print Path
    if not os.path.exists(Path+"/"+Folder):
        cmd_mkdir = "mkdir;"+Path+"/"+Folder
        subprocess.call(cmd_mkdir.split(";"), cwd=Path)
        print cmd_mkdir.split(";")


def GetRootFileList(Path, Channel):
    fRootFileList = []
    fHelpFileList = os.listdir(Path)
    for HelpFile in fHelpFileList:
        if ".root" in HelpFile:
            if Channel in HelpFile:
                fRootFileList.append(Path+"/"+HelpFile)
    return fRootFileList


def SortFilesIntoFolders(RootFileList, Label, Folder, Path):
    for RootFile in RootFileList:
        if Label+"_" in RootFile:
            cmd_mv = "mv;"+RootFile+";"+Path+"/"+Folder+"/"
            #print cmd_mv.split(";")
            subprocess.call(cmd_mv.split(";"), cwd=Path)
    

def MakePullNui(fNrNui, fNuiVar, fChannel, fPath):
    fHelpString = []
    fHelpString.append("NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_m1.5")
    fHelpString.append("NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_m1")
    fHelpString.append("NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_m0.5")
    fHelpString.append("NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fHelpString.append("NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0.5")
    fHelpString.append("NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p1")
    fHelpString.append("NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p1.5")
    fRootFileList = GetRootFileList(fPath, fChannel)
    for fString in fHelpString:
        SampleFolder = "Output_"+fChannel+"_F0_0.7_"+fString
        OutputFolder = "Output_"+fChannel+"_F0_0.7_"+fString+"_Plots"
        MakeFolder(SampleFolder, fPath)
        MakeFolder(OutputFolder, fPath)
        SortFilesIntoFolders(fRootFileList, fString, SampleFolder, fPath)
        cmd_run = fPath+"/runProducePlots;"+SampleFolder+";"+fChannel+";0.7;"+str(fNrNui)+";mode;"+OutputFolder
        print cmd_run.split(";")
        subprocess.call(cmd_run.split(";"), cwd=fPath)
        

def MakePullFraction(fNrNui, fNuiVar, fChannel, fPath):
    fHelpString = []
    fHelpString.append("_F0_0.4_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fHelpString.append("_F0_0.5_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fHelpString.append("_F0_0.6_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fHelpString.append("_F0_0.7_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fHelpString.append("_F0_0.8_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fHelpString.append("_F0_0.9_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fHelpString.append("_F0_1_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fF0 = []
    fF0.append(str(0.4))
    fF0.append(str(0.5))
    fF0.append(str(0.6))
    fF0.append(str(0.7))
    fF0.append(str(0.8))
    fF0.append(str(0.9))
    fF0.append(str(1))
    fRootFileList = GetRootFileList(fPath, fChannel)
    counter = 0
    for fString in fHelpString:
        SampleFolder = "Output_"+fChannel+fString
        OutputFolder = "Output_"+fChannel+fString+"_Plots"
        MakeFolder(SampleFolder, fPath)
        MakeFolder(OutputFolder, fPath)
        SortFilesIntoFolders(fRootFileList, fString, SampleFolder, fPath)
        cmd_run = fPath+"/runProducePlots;"+SampleFolder+";"+fChannel+";"+fF0[counter]+";"+str(fNrNui)+";mode;"+OutputFolder
        print cmd_run.split(";")
        subprocess.call(cmd_run.split(";"), cwd=fPath)
        counter = counter+1
        
def MakePullDatafit(fNrNui, fNuiVar, fChannel, fPath):	
    fHelpString = []
    fHelpString.append("_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar)+"_p0")
    fF0 = []
    fF0.append(str(0)) #dummy value...
    fRootFileList = GetRootFileList(fPath, fChannel)
    counter = 0
    #for fString in fHelpString:
    SampleFolder = "DatafitPseudo_Output_"+fChannel+fHelpString[counter]
    OutputFolder = "DatafitPseudo_Output_"+fChannel+fHelpString[counter]+"_Plots"
    DatafitFolder = "Datafit_Output"
    MakeFolder(SampleFolder, fPath)
    MakeFolder(OutputFolder, fPath)
    MakeFolder(DatafitFolder, fPath)
    SortFilesIntoFolders(fRootFileList, fHelpString[counter], SampleFolder, fPath)
    SortFilesIntoFolders(fRootFileList, "Datafit_Output", DatafitFolder, fPath)
    cmd_run = fPath+"/runProducePlots;"+SampleFolder+";"+fChannel+";"+fF0[counter]+";"+str(fNrNui)+";Datafit;"+OutputFolder
    print cmd_run.split(";")
    subprocess.call(cmd_run.split(";"), cwd=fPath)
    
def MakePullFluctuate(fNrNui, fNuiVar, fChannel, fPath):	
    fHelpString = []
    fHelpString.append("_NrNui"+str(fNrNui)+"_NuiVar"+str(fNuiVar) )
    fF0 = []
    fF0.append(str(0)) #dummy value...
    fRootFileList = GetRootFileList(fPath, fChannel)
    counter = 0
    #for fString in fHelpString:
    SampleFolder = "OutputFluctuate_"+fChannel+fHelpString[counter]
    OutputFolder = "OutputFluctuate_"+fChannel+fHelpString[counter]+"_Plots"
    FluctuateFolder = "OutputFluctuate"
    MakeFolder(SampleFolder, fPath)
    MakeFolder(OutputFolder, fPath)
    MakeFolder(FluctuateFolder, fPath)
    SortFilesIntoFolders(fRootFileList, fHelpString[counter], SampleFolder, fPath)
    SortFilesIntoFolders(fRootFileList, "OutputFluctuate", FluctuateFolder, fPath)
    cmd_run = fPath+"/runProducePlots;"+SampleFolder+";"+fChannel+";"+fF0[counter]+";"+str(fNrNui)+";Datafit;"+OutputFolder
    print cmd_run.split(";")
    subprocess.call(cmd_run.split(";"), cwd=fPath)
