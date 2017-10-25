#!/usr/bin/python

import os, sys, operator

#sys.argv.append( '-b-' )

import ROOT

if len(sys.argv) ==1 :
    print "NO OPTION!!! ABORT"
    exit(0)

region = sys.argv[1]

filedir = ''
if region == "lep":
    filedir = '../ExternalSystematicsOutput_el_mu_bTag_2incl_3D_3W_lep'
elif region == "had":
    #filedir = '../ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had'
    filedir = '../../FitPackage_allinOne_Apr26_allCF/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had'
elif region == "lephad":
    filedir = '../ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad'
else:
    print "nope. Problem with input name: "+region
    print "ABORT"
    exit(0)

ct10_F0 = open ('CT10_F0_'+region+'.txt','w')
ct10_FL = open ('CT10_FL_'+region+'.txt','w')
ct10_FR = open ('CT10_FR_'+region+'.txt','w')
mstw_F0 = open ('MSTW_F0_'+region+'.txt','w')
mstw_FL = open ('MSTW_FL_'+region+'.txt','w')
mstw_FR = open ('MSTW_FR_'+region+'.txt','w')
nnpdf_F0 = open ('NNPDF_F0_'+region+'.txt','w')
nnpdf_FL = open ('NNPDF_FL_'+region+'.txt','w')
nnpdf_FR = open ('NNPDF_FR_'+region+'.txt','w')


#for subdir, dirs, files in os.walk(filedir+"_pdfSysts_Mar11Templates"):
for subdir, dirs, files in os.walk(filedir):
    flist = []
    fdict = {}
    for file in files:
        if ".root" in file and ("CT10" in file or "MSTW" in file or "NNPDF" in file):
            flist.append(str(file))
            #print file, str(file).split('/')[0].split('_')[2]
            fdict[str(file)]=int(str(file).split('/')[0].split('_')[2])
            
    flist.sort()
    sorted_fdict = sorted(fdict.items(), key=operator.itemgetter(1))
    print sorted_fdict
    print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$"

    for entry in sorted_fdict:
        #print os.path.join(subdir, file)
        file = entry[0]
        filepath = subdir + os.sep + file

        if filepath.endswith(".root"):
            #print filepath
            f = ROOT.TFile.Open(filepath, "read")
            t = f.Get("EnsembleTree")
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            #print "F0:",fgaus0.GetParameter(1)
            
            if "CT10" in filepath:
                #print filepath
                ct10_F0.write("{0}\n".format(fgaus0.GetParameter(1)))
                ct10_FL.write("{0}\n".format(fgausL.GetParameter(1)))
                ct10_FR.write("{0}\n".format(fgausR.GetParameter(1)))

            if "MSTW" in filepath:
                print filepath
                mstw_F0.write("{0}\n".format(fgaus0.GetParameter(1)))
                mstw_FL.write("{0}\n".format(fgausL.GetParameter(1)))
                mstw_FR.write("{0}\n".format(fgausR.GetParameter(1)))

            if "NNPDF" in filepath:
                print filepath
                nnpdf_F0.write("{0}\n".format(fgaus0.GetParameter(1)))
                nnpdf_FL.write("{0}\n".format(fgausL.GetParameter(1)))
                nnpdf_FR.write("{0}\n".format(fgausR.GetParameter(1)))



            

            
