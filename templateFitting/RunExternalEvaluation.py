#! /usr/bin/python

import subprocess,sys
import os,string
import os.path

fPath = os.getcwd()

fSyst = []
#fSyst.append("Shower")
#fSyst.append("MCGen1")
#fSyst.append("MCGen2")
#fSyst.append("ISR_FSR")
#fSyst.append("mu_UE")
#fSyst.append("CR_Perugia")
#fSyst.append("mu_TopMass")

#fSyst.append("PDF_CT10_0")
#fSyst.append("PDF_MSTW_0")
#fSyst.append("PDF_NNPDF_0")
#fSyst.append("el_QCD_FAKE")
#fSyst.append("el_QCD_REAL")

#fSyst.append("wjets_bb4")
#fSyst.append("el_QCDshape")
fSyst.append("PrimaryVertex")
#fSyst.append("Datafit")
#fSyst.append("jes")
#fSyst.append("jer")
#fSyst.append("Calibration")
#fSyst.append("TemplateStat")
#fSyst.append("jeff")
#fSyst.append("muid")
#fSyst.append("mums")
#fSyst.append("btag")
#fSyst.append("ctag")
#fSyst.append("mistag")
#fSyst.append("musc")
#fSyst.append("muonid")
#fSyst.append("muonreco")
#fSyst.append("muontrig")
#fSyst.append("jvfsf")
#fSyst.append("cellout")
#fSyst.append("softjet")
#fSyst.append("pileup")

#fSyst.append("mu_QCDShape")

for Syst in fSyst:
    cmd = fPath+"/runEvaluateExternalPileup;el_mu;0;"+Syst+";3D;2"
    print cmd.split(";")
    subprocess.call(cmd.split(";"), cwd=fPath)
    cmd = fPath+"/runEvaluateExternalPileup;el_mu;0;"+Syst+";3D;3"
    print cmd.split(";")
    subprocess.call(cmd.split(";"), cwd=fPath)
    cmd = fPath+"/runEvaluateExternalPileup;el_mu;0;"+Syst+";3D;4"
    print cmd.split(";")
    subprocess.call(cmd.split(";"), cwd=fPath)
    cmd = fPath+"/runEvaluateExternalPileup;el_mu;0;"+Syst+";3D;5"
    print cmd.split(";")
    subprocess.call(cmd.split(";"), cwd=fPath)
    cmd = fPath+"/runEvaluateExternalPileup;el_mu;0;"+Syst+";3D;6"
    print cmd.split(";")
    subprocess.call(cmd.split(";"), cwd=fPath)
    cmd = fPath+"/runEvaluateExternalPileup;el_mu;0;"+Syst+";3D;7"
    print cmd.split(";")
    subprocess.call(cmd.split(";"), cwd=fPath)
    cmd = fPath+"/runEvaluateExternalPileup;el_mu;0;PrimaryVertexEval;3D;6"
    print cmd.split(";")
    subprocess.call(cmd.split(";"), cwd=fPath)
            
    
