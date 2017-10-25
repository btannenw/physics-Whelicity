#!/usr/bin/python

import os, sys, operator

#sys.argv.append( '-b-' )

import ROOT

if len(sys.argv) == 1 :
    print "NO OPTION!!! ABORT"
    exit(0)

region = sys.argv[1]
#fraction = sys.argv[2]

filedir = ''
outname = ''
if region == "lep":
    #filedir = 'ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep'
    #outname = filedir+"/"+"SystematicOutput_el_mu.script.txt"
    #filedir = 'ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_lep'
    filedir = '../syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_lep' # aug 6, 2016
    outname = filedir+"/"+"SystematicOutput_el_mu_bTag.script.txt"
elif region == "had":
    #filedir = 'ExternalSystematicsOutput_el_mu_2incl_3D_3W_had'
    #outname = filedir+"/"+"SystematicOutput_el_mu.script.txt"
    #filedir = 'ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had' 
    filedir = '../syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had' # aug 6, 2016
    outname = filedir+"/"+"SystematicOutput_el_mu_bTag.script.txt"
elif region == "lephad":
    #filedir = 'ExternalSystematicsOutput_el_mu_lephad_2incl_3D_3W_lephad'
    #outname = filedir+"/"+"SystematicOutput_el_mu_lephad.script.txt"
    filedir = 'ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad'
    filedir = '../syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_lephad_2incl_3D_3W_lephad' # aug 6, 2016
    outname = filedir+"/"+"SystematicOutput_el_mu_lephad_bTag.script.txt"
elif region == "lep_2incl":
    #filedir = 'ExternalSystematicsOutput_el_mu_2incl_3D_3W_had'
    #outname = filedir+"/"+"SystematicOutput_el_mu.script.txt
    #filedir = '../syst_outputs_05Jul/ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep'
    filedir = '../syst_outputs_lep_2incl_aug02/ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep/'
    outname = filedir+"/"+"SystematicOutput_el_mu_lep.script.txt"
elif region == "had_2incl":
    #filedir = 'ExternalSystematicsOutput_el_mu_2incl_3D_3W_had'
    #outname = filedir+"/"+"SystematicOutput_el_mu.script.txt
    #filedir = '../syst_outputs_05Jul/ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep'
    filedir = '../syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_2incl_3D_3W_had/'
    outname = filedir+"/"+"SystematicOutput_el_mu_had.script.txt"
elif region == "lephad_2incl":
    #filedir = 'ExternalSystematicsOutput_el_mu_2incl_3D_3W_had'
    #outname = filedir+"/"+"SystematicOutput_el_mu.script.txt"
    filedir = '../syst_outputs_05Jul/ExternalSystematicsOutput_el_mu_lephad_2incl_3D_3W_lephad' 
    outname = filedir+"/"+"SystematicOutput_el_mu_lephad.script.txt"
elif region == "lephad_bTag":
    #filedir = 'ExternalSystematicsOutput_el_mu_2incl_3D_3W_had'
    #outname = filedir+"/"+"SystematicOutput_el_mu.script.txt"
    #filedir = '../syst_outputs_05Jul/JERnp/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/' 
    
    # Aug 03, 2016- inclusion of 17% single top systematic
    filedir = '/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_outputs_8ch_03Aug/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad'
    outname = filedir+"/"+"SystematicOutput_el_mu_lephad_bTag.script.txt"
else:
    print "nope. Problem with input name: "+region
    print "ABORT"
    exit(0)


dumptxt = open(outname,'w')
print outname
nom=''
nomF0 = 0
nomFL = 0
nomFR = 0

nomF0_afii = 0
nomFL_afii = 0
nomFR_afii = 0

psF0_herwig = 0
psFL_herwig = 0
psFR_herwig = 0
psF0_pythia = 0
psFL_pythia = 0
psFR_pythia = 0


meF0_nlo = 0
meFL_nlo = 0
meFR_nlo = 0
meF0_powheg = 0
meFL_powheg = 0
meFR_powheg = 0


for subdir, dirs, files in os.walk(filedir):
    flist = []
    fdict = {}
    for file in files:

        if ".root" in file and "Nominal" in file and 'AFII' not in file:
            nom = str(filedir+"/"+file)
            
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            nomF0 = fgaus0.GetParameter(1)
            nomFL = fgausL.GetParameter(1)
            nomFR = fgausR.GetParameter(1)

        if ".root" in file and "Nominal" in file and 'AFII' in file:
            nom = str(filedir+"/"+file)
            
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            nomF0_afii = fgaus0.GetParameter(1)
            nomFL_afii = fgausL.GetParameter(1)
            nomFR_afii = fgausR.GetParameter(1)

        if ".root" in file and 'Powheg+fHerwig(PS)' in file:
            nom = str(filedir+"/"+file)
            print file
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            psF0_herwig = fgaus0.GetParameter(1)
            psFL_herwig = fgausL.GetParameter(1)
            psFR_herwig = fgausR.GetParameter(1)

        if ".root" in file and 'Powheg+Pythia6(PS)' in file:
            nom = str(filedir+"/"+file)
            print file
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            psF0_pythia = fgaus0.GetParameter(1)
            psFL_pythia = fgausL.GetParameter(1)
            psFR_pythia = fgausR.GetParameter(1)

        if ".root" in file and 'NLO' in file:
            nom = str(filedir+"/"+file)
            print file
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            meF0_nlo = fgaus0.GetParameter(1)
            meFL_nlo = fgausL.GetParameter(1)
            meFR_nlo = fgausR.GetParameter(1)

        if ".root" in file and 'Powheg+fHerwig' in file and '(PS)' not in file:
            nom = str(filedir+"/"+file)
            print file
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            meF0_powheg = fgaus0.GetParameter(1)
            meFL_powheg = fgausL.GetParameter(1)
            meFR_powheg = fgausR.GetParameter(1)


for subdir, dirs, files in os.walk(filedir):

    for file in files:
        filepath = subdir + os.sep + file

        if filepath.endswith(".root") and ("CT10" not in file and "MSTW" not in file and "NNPDF" not in file and "Mass" not in file and "Nominal" not in file):
            #print file
            name = file.split('Syst_')[1].split("_el_mu")[0]

            f = ROOT.TFile.Open(filepath, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")

            if 'JER10' not in name:
                if 'hdamp' not in name and 'PS' not in name and 'NLO' not in name and 'Powheg+fHerwig' not in file:
                    dumptxt.write(name+"\t& "+str(nomF0-fgaus0.GetParameter(1))+"\t& "+str(nomFL-fgausL.GetParameter(1))+"\t& "+str(nomFR-fgausR.GetParameter(1))+" \ q\n")
                elif 'hdamp' in name:
                    dumptxt.write(name+"\t& "+str(nomF0_afii-fgaus0.GetParameter(1))+"\t& "+str(nomFL_afii-fgausL.GetParameter(1))+"\t& "+str(nomFR_afii-fgausR.GetParameter(1))+" \ q\n")
                elif 'Herwig(PS)' in name:
                    dumptxt.write('PartonShower'+"\t& "+str(psF0_herwig-psF0_pythia)+"\t& "+str(psFL_herwig-psFL_pythia)+"\t& "+str(psFR_herwig-psFR_pythia)+" \ q\n")
                elif 'NLO' in name:                    
                    dumptxt.write('MEgenerator'+"\t& "+str(meF0_nlo-meF0_powheg)+"\t& "+str(meFL_nlo-meFL_powheg)+"\t& "+str(meFR_nlo-meFR_powheg)+" \ q\n")
