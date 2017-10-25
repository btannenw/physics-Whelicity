#! /usr/bin/python

import subprocess,sys
import os,string
import os.path

fPath  = os.getcwd()

fRocksPath = fPath  #fPath.replace("/net/ph2/auto/home", "/rocks")

#Determine number of nuisance parameters
nui_num = 0
for num in range(0,99):
    num_str=str(num)	
    if os.path.exists(fRocksPath+"/Output_Sample_mu_0_nui"+num_str+"_0.root"):
      nui_num = num
    elif os.path.exists(fRocksPath+"/Output_Sample_el_0_nui"+num_str+"_0.root"):
      nui_num = num
nui_num_str= str(nui_num)
print 'Number of nuisance parameters : '+nui_num_str

#Run DoMakePlots

for n in range(0,8):
    fFolderList = []
    n_str=str(n)
    if os.path.exists(fRocksPath+"/Output_Sample_mu_"+n_str+"_nui"+nui_num_str +"_1.root") ==True:
      fFolderList.append("Output_Sample_mu_"+n_str)
      for fFolder in fFolderList:
        #Make Folder for OutputFiles
        cmd_folder     = "mkdir;"+fPath+"/"+fFolder
        subprocess.call(cmd_folder.split(";"),     cwd=fPath)
        fSubFolderList = os.listdir(fRocksPath)
        for fSubFolder in fSubFolderList:
          if fFolder in fSubFolder:
            cmd_cp         = "cp;"+fRocksPath+"/"+fSubFolder+";"+fPath+"/"+fFolder+"/"
            subprocess.call(cmd_cp.split(";"), cwd=fPath)
        cmd_run = "./runProducePlots;"+fFolder+";mu;"+n_str+";"+nui_num_str+";mode;"+fFolder+"_Plots"
        subprocess.call(cmd_run.split(";"), cwd=fPath)
    elif os.path.exists(fRocksPath+"/Output_Sample_el_"+n_str+"_nui"+nui_num_str +"_1.root") ==True:	
      fFolderList.append("Output_Sample_el_"+n_str)
      for fFolder in fFolderList:
        #Make Folder for OutputFiles
        cmd_folder     = "mkdir;"+fPath+"/"+fFolder
        subprocess.call(cmd_folder.split(";"),     cwd=fPath)
        fSubFolderList = os.listdir(fRocksPath)
        for fSubFolder in fSubFolderList:
          if fFolder in fSubFolder:
            cmd_cp         = "cp;"+fRocksPath+"/"+fSubFolder+";"+fPath+"/"+fFolder+"/"
            subprocess.call(cmd_cp.split(";"), cwd=fPath)
        cmd_run = "./runProducePlots;"+fFolder+";el;"+n_str+";"+nui_num_str+";mode;"+fFolder+"_Plots"
        subprocess.call(cmd_run.split(";"), cwd=fPath)
    else:
      print 'File with number '+n_str+' does not exist'



for n in range(0,8):
    fFolderList2 = []
    #channel_str
    n_str=str(n)
    if os.path.exists(fRocksPath+"/Output_Sample_mu_"+n_str+"_nui"+nui_num_str +"_1.root") ==True:
      fFolderList2.append("Output_Sample_mu_"+n_str)
      for fFolder in fFolderList2:
        #Make Folder for OutputFiles
        cmd_folder     = "mkdir;"+fPath+"/"+fFolder
        subprocess.call(cmd_folder.split(";"),     cwd=fPath)
        fSubFolderList = os.listdir(fRocksPath)
        for fSubFolder in fSubFolderList:
          if fFolder in fSubFolder:
            cmd_cp         = "cp;"+fRocksPath+"/"+fSubFolder+";"+fPath+"/"+fFolder+"/"
            subprocess.call(cmd_cp.split(";"), cwd=fPath)
      channel_str="mu"
      ncounter = n
    elif os.path.exists(fRocksPath+"/Output_Sample_el_"+n_str+"_nui"+nui_num_str +"_1.root") ==True:	
      fFolderList2.append("Output_Sample_el_"+n_str)
      for fFolder in fFolderList2:
        #Make Folder for OutputFiles
        cmd_folder     = "mkdir;"+fPath+"/"+fFolder
        subprocess.call(cmd_folder.split(";"),     cwd=fPath)
        fSubFolderList = os.listdir(fRocksPath)
        for fSubFolder in fSubFolderList:
          if fFolder in fSubFolder:
            cmd_cp         = "cp;"+fRocksPath+"/"+fSubFolder+";"+fPath+"/"+fFolder+"/"
            subprocess.call(cmd_cp.split(";"), cwd=fPath)
      channel_str="el"
      ncounter = n
ncounter_str=str(ncounter)      
print 'File with highest number : '+ncounter_str
cmd_run = "./runProducePlotsCal;"+fFolder+";"+channel_str+";"+ncounter_str+";"+nui_num_str+";mode;Output_"+channel_str+"_Cal_Plots"
subprocess.call(cmd_run.split(";"), cwd=fPath)


#try:
#    some_object
#except NameError:
#    do_something()
#else:
#    do_something_else()
