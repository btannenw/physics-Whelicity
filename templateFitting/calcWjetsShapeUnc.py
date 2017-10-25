#!/usr/bin/python

import os, sys, operator

#sys.argv.append( '-b-' )

import ROOT


def getMean(variedDir, nominalDir, region):

    means = []
    filename = ''
    if region == "lephad":
        filename = '/Syst_Nominal_el_mu_lephad_bTag.root'
    else:
        filename = '/Syst_Nominal_el_mu_bTag.root'

    nom = nominalDir + filename 
    var = variedDir + filename
    fgaus0 = ROOT.TF1("fgaus0","gaus")
    fgausL = ROOT.TF1("fgausL","gaus")
    fgausR = ROOT.TF1("fgausR","gaus")
   
    fnom = ROOT.TFile.Open(nom, "read")
    tnom = fnom.Get("EnsembleTree")
    tnom.Fit("fgaus0","F0","","Q")
    tnom.Fit("fgausL","FL","","Q")
    tnom.Fit("fgausR","FR","","Q")
    nomF0 = fgaus0.GetParameter(1)
    nomFL = fgausL.GetParameter(1)
    nomFR = fgausR.GetParameter(1)

    fvar = ROOT.TFile.Open(var, "read")
    tvar = fvar.Get("EnsembleTree")
    tvar.Fit("fgaus0","F0","","Q")
    tvar.Fit("fgausL","FL","","Q")
    tvar.Fit("fgausR","FR","","Q")
    varF0 = fgaus0.GetParameter(1)
    varFL = fgausL.GetParameter(1)
    varFR = fgausR.GetParameter(1)

    #print "nomF0:",nomF0,"nomFL:",nomFL,"nomFR:",nomFR,"tnom.GetEntries():",tnom.GetEntries()
    #print "varF0:",varF0,"varFL:",varFL,"varFR:",varFR,"tvar.GetEntries():",tvar.GetEntries()
    
    means.append(nomF0-varF0)
    means.append(nomFL-varFL)
    means.append(nomFR-varFR)

    return means

if len(sys.argv) != 2 :
    print "NO OPTION!!! ABORT"
    exit(0)

region = sys.argv[1]
#fixed = sys.argv[2]

if region != "lep" and region != "had" and region != "lephad" :
    print "region",region, "not recognized. ABORT"
    exit(0)

#if fixed != "light" and fixed != "c" and fixed != "bbcc" :
#    print "fixed",fixed, "not recognized. ABORT"
#    exit(0)


fileDir_light = ''
fileDir_c = ''
fileDir_bbcc = ''
#outname = ''

topDir_light  = "../FitPackage_allinOne_May4_CFlight_fixOtherW/"
topDir_c      = "../FitPackage_allinOne_May4_CFc_fixOtherW/"
topDir_bbcc   = "../FitPackage_allinOne_May4_CFbbcc_fixOtherW/"
regionDir     = ''

if region == "lep":
    regionDir = 'ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_lep'
elif region == "had":
    regionDir = 'ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had'
elif region == "lephad":
    regionDir = 'ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad'
else:
    print "nope. Problem with input name: "+region
    print "ABORT"
    exit(0)

fileDir_light = topDir_light + regionDir
fileDir_c     = topDir_c + regionDir
fileDir_bbcc  = topDir_bbcc + regionDir

lightMeans = getMean(fileDir_light, regionDir, region)
cMeans = getMean(fileDir_c, regionDir, region)
bbccMeans = getMean(fileDir_bbcc, regionDir, region)

envelope = lightMeans + cMeans + bbccMeans

print "LightMeans", lightMeans
print "cMeans", cMeans
print "bbccMeans", bbccMeans

print "Envelope = ",max(envelope),"-",min(envelope),"=",max(envelope)-min(envelope)
