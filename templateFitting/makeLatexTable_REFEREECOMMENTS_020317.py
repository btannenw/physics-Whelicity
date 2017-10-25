#!/usr/bin/python

import os, sys, math, operator
import ROOT

#f_lephadBTag = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_outputs_8ch_03Aug/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/SystematicOutput_el_mu_lephad_bTag.script.txt') # Aug 3, 2016- inclusion of 17% single top unc

# Mohammad Test (2incl only) Jul 5, 2016
#f_hadBTag = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_outputs_05Jul/ExternalSystematicsOutput_el_mu_lephad_2incl_3D_3W_lephad/SystematicOutput_el_mu_lephad.script.txt')
#f_lepBTag = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_outputs_05Jul/ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep/SystematicOutput_el_mu_lep.script.txt')

# NORMAL FILES TO USE (JUL 5, 2016)
#f_lephadBTag=open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/FitPackage_allinOne_Apr26_allCF/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/SystematicOutput_el_mu_lephad_bTag.script.txt')

# Aug 08 Files
f_lep2incl    = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_outputs_lep_2incl_aug02/ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep/SystematicOutput_el_mu_lep.script.txt', 'r')
f_had2incl    = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_2incl_3D_3W_had/SystematicOutput_el_mu_had.script.txt', 'r')
f_lephad2incl = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_lephad_2incl_3D_3W_lephad/SystematicOutput_el_mu_lephad_bTag.script.txt', 'r')
f_lepBTag     = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_lep/SystematicOutput_el_mu_bTag.script.txt', 'r')
f_hadBTag     = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_output_06Aug_2016/ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_had/SystematicOutput_el_mu_bTag.script.txt', 'r')
f_lephadBTag  = open('/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_outputs_8ch_03Aug/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/SystematicOutput_el_mu_lephad_bTag.script.txt', 'r')


d_lep2incl={}
d_had2incl={}
d_lephad2incl={}
d_lepBTag={}
d_hadBTag={}
d_lephadBTag={}

statBkgUnc = [0,0,0,0]
PDFunc = [0,0,0,0]
TopMassUnc = [0,0,0,0]
TemplateStatUnc = [0,0,0,0]

tablelist={}

def makeFullSensitivityTable ( whichFraction, lepDict, hadDict, lephadDict, method=''):
   """This function takes one dictionary as input and outputs the correctly formatted final table with sources broken down by type and categorized byby 0/L/R fraction passed from whichFraction. Additionally provides % effect of giving systematic w.r.t. expected value of helicity fraction"""

   #print "Q#$%#%@$: ",sumMethod
   fraction = -1
   nomFrac = 0
   #statBkgUnc = 0;
   #PDFunc = 0;
   #TopMassUnc = 0;
   if whichFraction == 'F0':
      fraction = 0
      nomFrac = 0.697
      #statBkgUnc = [0.010, 0.012, 0.007]
      #PDFunc = [0.003, 0.001, 0.002]
      #TopMassUnc = [0.002, 0.003, 0.001]
      #TemplateStatUnc = [0.007, 0.007, 0.005]

   if whichFraction == 'FL':
      fraction = 1
      nomFrac = 0.301
      #statBkgUnc = [0.007, 0.022, 0.005]
      #PDFunc = [0.004, 0.002, 0.003]
      #TopMassUnc = [0.005, 0.009, 0.003]
      #TemplateStatUnc = [0.005, 0.016, 0.003]

   if whichFraction == 'FR':
      fraction = 2
      nomFrac = 0.002
      #statBkgUnc = [0.005, 0.026, 0.004]
      #PDFunc = [0.001, 0.002, 0.002]
      #TopMassUnc = [0.003, 0.006, 0.004]
      #TemplateStatUnc = [0.003, 0.016, 0.003]

   print "$$$$$$$$$$$$$", statBkgUnc[0], statBkgUnc[1], statBkgUnc[2]
   lepTotals = calcTotals (lepDict, whichFraction, statBkgUnc[0], PDFunc[0], TopMassUnc[0], TemplateStatUnc[0])
   hadTotals = calcTotals (hadDict, whichFraction, statBkgUnc[1], PDFunc[1], TopMassUnc[1], TemplateStatUnc[1])
   lephadTotals = calcTotals (lephadDict, whichFraction, statBkgUnc[2], PDFunc[2], TopMassUnc[2], TemplateStatUnc[2])

   print "\n### el+mu 2incl {0} ###".format(whichFraction)
   print "\\begin{table}[h!]"
   print "\centering"
   print "\\begin{tabular}{lcccc}"
   print "\hline\hline"
   if whichFraction == 'F0':
      print "\multicolumn{5}{c}{\\fo}\\\\"
   elif whichFraction == 'FL':
      print "\multicolumn{5}{c}{\\fl}\\\\"
   elif whichFraction == 'FR':
      print "\multicolumn{5}{c}{\\fr}\\\\"
   print "\hline"
   print "Systematic uncertainty & N$_{syst}$ & Lep 2incl & Lep+Had 2incl & Lep+Had 1excl+2incl \\\\\\hline"
   print "\multicolumn{5}{c}{Reconstructed Objects} \\\\\\hline"
   if method == 'alt':
      
      printThreePlusMinusLine('Muon', 6, lephadTotals[0], lepTotals[1], lepTotals[2], hadTotals[1], hadTotals[2], lephadTotals[1], lephadTotals[2])
      printThreePlusMinusLine('Electron', 5, lephadTotals[3], lepTotals[4], lepTotals[5], hadTotals[4], hadTotals[5], lephadTotals[4], lephadTotals[5])
      printThreePlusMinusLine('JES', 26, lephadTotals[6], lepTotals[7], lepTotals[8], hadTotals[7], hadTotals[8], lephadTotals[7], lephadTotals[8])
      printThreePlusMinusLine('JER', 11, lephadTotals[9], lepTotals[10], lepTotals[11], hadTotals[10], hadTotals[11], lephadTotals[10], lephadTotals[11])
      printThreePlusMinusLine('JVF', 1, lephadTotals[12], lepTotals[13], lepTotals[14], hadTotals[13], hadTotals[14], lephadTotals[13], lephadTotals[14])
      #printThreePlusMinusLine('JEff', 1, lephadTotals[15], lepTotals[16], lepTotals[17], hadTotals[16], hadTotals[17], lephadTotals[16], lephadTotals[17])
      printThreePlusMinusLine('\\bt tagging', 3, 3, lepTotals[18], lepTotals[19], hadTotals[18], hadTotals[19], lephadTotals[18], lephadTotals[19])
      #printThreePlusMinusLine('MET', 2, lephadTotals[20], lepTotals[21], lepTotals[22], hadTotals[21], hadTotals[22], lephadTotals[21], lephadTotals[22])
      print "\n\\hline\\hline"
      printThreePlusMinusLine('Sum of Reco Objects', ' ',' ', lepTotals[23], lepTotals[24], hadTotals[23], hadTotals[24], lephadTotals[23], lephadTotals[24])
      print "\n\\hline"
      print "\multicolumn{5}{c}{Modeling} \\\\\\hline"
      #printThreePlusMinusLine('Radiation', 1, n_rad, sum_up_rad, sum_down_rad, nomFrac, lephadTotals[43])
      printThreePlusMinusLine('Radiation, radLo', 1, lephadTotals[25], lepTotals[26], lepTotals[28], hadTotals[26], hadTotals[28], lephadTotals[26], lephadTotals[28])
      printThreePlusMinusLine('Radiation, radHi', 1, lephadTotals[25], lepTotals[27], lepTotals[28], hadTotals[27], hadTotals[28], lephadTotals[27], lephadTotals[28])
      printThreePlusMinusLine('Parton Shower', 1, lephadTotals[29], lepTotals[30], lepTotals[31], hadTotals[30], hadTotals[31], lephadTotals[30], lephadTotals[31])
      #printThreePlusMinusLine('Color Reconnection', 1, lephadTotals[32], lepTotals[33], lepTotals[34], hadTotals[33], hadTotals[34], lephadTotals[33], lephadTotals[34])
      #printThreePlusMinusLine('Underlying Event', 1, lephadTotals[35], lepTotals[36], lepTotals[37], hadTotals[36], hadTotals[37], lephadTotals[36], lephadTotals[37])
      printThreePlusMinusLine('ME Generator', 1, lephadTotals[38], lepTotals[39], lepTotals[40], hadTotals[39], hadTotals[40], lephadTotals[39], lephadTotals[40])
      printThreePlusMinusLine('PDF', 3, 3, lepTotals[45], lepTotals[45], hadTotals[45], hadTotals[45], lephadTotals[45], lephadTotals[45])
      printThreePlusMinusLine('Top Mass', 3, 3, lepTotals[46], lepTotals[46], hadTotals[46], hadTotals[46], lephadTotals[46], lephadTotals[46])
      print "\n\\hline\\hline"  
      printThreePlusMinusLine('Sum of Modeling', ' ', ' ', lepTotals[41], lepTotals[42], hadTotals[41], hadTotals[42], lephadTotals[41], lephadTotals[42])
      print "\n\\hline"
      print "\multicolumn{5}{c}{Method Uncertainty} \\\\\\hline"
      printThreePlusMinusLine('Template Statistics', 3, 3, lepTotals[47], lepTotals[47], hadTotals[47], hadTotals[47], lephadTotals[47], lephadTotals[47])
      print "\n\\hline\\hline"
      #printThreeTotalLine('Total Syst.', lepTotals[23] + lepTotals[41], lepTotals[24] + lepTotals[42], hadTotals[23] + hadTotals[41], hadTotals[24] + hadTotals[42], lephadTotals[23] + lephadTotals[41], lephadTotals[24] + lephadTotals[42])
      printThreePlusMinusLine('Total Syst.', ' ', ' ', lepTotals[23] + lepTotals[41] + lepTotals[47], lepTotals[24] + lepTotals[42] + lepTotals[47], hadTotals[23] + hadTotals[41] + hadTotals[47], hadTotals[24] + hadTotals[42] + hadTotals[47], lephadTotals[23] + lephadTotals[41] + lephadTotals[47], lephadTotals[24] + lephadTotals[42] + lephadTotals[47])
      print "Stat. + Bkg. & - & {0} & {1} & {2} \\\\\\hline".format(round(lepTotals[44],4), round(hadTotals[44],4), round(lephadTotals[44],4))
      print "\n\\hline\\hline"
      #print "Total &  & {0} & {1} & {2} \\\\\\hline\\hline".format(round(math.sqrt(lepTotals[43]),4), round(math.sqrt(hadTotals[43]),4), round(math.sqrt(lephadTotals[43]),4))

#   else:
#      print "Muon & 6({0}) & {1} & {2} \\\\\\hline".format(n_mu, round(math.sqrt(sum_mu),4), round(100*math.sqrt(sum_mu)/nomFrac,1))
#      print "Electron & 5({0}) & {1} & {2} \\\\\\hline".format(n_el, round(math.sqrt(sum_el),4), round(100*math.sqrt(sum_el)/nomFrac,1))
#      print "JES & 26({0}) & {1} & {2} \\\\\\hline".format(n_JES, round(math.sqrt(sum_JES),4), round(100*math.sqrt(sum_JES)/nomFrac,1))
#      print "JER & 1({0}) & {1} & {2} \\\\\\hline".format(n_JER, round(math.sqrt(sum_JER),4), round(100*math.sqrt(sum_JER)/nomFrac,1))
#      print "Jet Reco & 1({0}) & {1} & {2} \\\\\\hline".format(n_JEff, round(math.sqrt(sum_JEff),4), round(100*math.sqrt(sum_JEff)/nomFrac,1))
#      print "JVF & 1({0}) & {1} & {2} \\\\\\hline".format(n_JVF, round(math.sqrt(sum_JVF),4), round(100*math.sqrt(sum_JVF)/nomFrac,1))
#      print "\\bt tagging & 3({0}) & {1} & {2} \\\\\\hline".format(3, round(math.sqrt(sum_btag),4), round(100*math.sqrt(sum_btag)/nomFrac,1))
#      print "MET & 2({0}) & {1} & {2} \\\\\\hline".format(n_met, round(math.sqrt(sum_met),4), round(100*math.sqrt(sum_met)/nomFrac,1))
#      print "\n\\hline\\hline"
#      print "Sum & & {0} & {1} \\\\\\hline".format(round(math.sqrt(sum_el+sum_mu+sum_JES+sum_JER+sum_btag+sum_met),4), round(100*math.sqrt(sum_el+sum_mu+sum_JES+sum_JER+sum_btag+sum_met)/nomFrac,1))
#      print "\n\\hline"
#      print "\multicolumn{4}{c}{Modeling} \\\\\\hline"
#      print "Radiation & 1({0}) & {1} & {2} \\\\\\hline".format(n_rad, round(math.sqrt(sum_rad),4), round(100*math.sqrt(sum_rad)/nomFrac,1))
#      print "Parton Shower & 1({0}) & {1} & {2} \\\\\\hline".format(n_ps, round(math.sqrt(sum_ps),4), round(100*math.sqrt(sum_ps)/nomFrac,1))
#      print "Color Reconnection & 1({0}) & {1} & {2} \\\\\\hline".format(n_cr, round(math.sqrt(sum_cr),4), round(100*math.sqrt(sum_cr)/nomFrac,1))
#      print "Underlying Event & 1({0}) & {1} & {2} \\\\\\hline".format(n_ue, round(math.sqrt(sum_ue),4), round(100*math.sqrt(sum_ue)/nomFrac,1))
#      print "MC Generator & 1({0}) & {1} & {2} \\\\\\hline".format(n_mcgen, round(math.sqrt(sum_mcgen),4), round(100*math.sqrt(sum_mcgen)/nomFrac,1))
#      print "\n\\hline\\hline"  
#      print "Sum & & {0} & {1} \\\\\\hline".format(round(math.sqrt(sum_rad+sum_ps+sum_cr+sum_ue+sum_mcgen),4), round(100*math.sqrt(sum_rad+sum_ps+sum_cr+sum_ue+sum_mcgen)/nomFrac,1))
#      print "\n\\hline\\hline"
#      print "Total Syst. & - & {0} & {1} \\\\\\hline".format(round(math.sqrt(sum_el+sum_mu+sum_JES+sum_JER+sum_btag+sum_met+sum_rad+sum_ps+sum_cr+sum_ue+sum_mcgen),4), round(100*math.sqrt(sum_el+sum_mu+sum_JES+sum_JER+sum_btag+sum_met+sum_rad+sum_ps+sum_cr+sum_ue+sum_mcgen)/nomFrac,1))
#      print "Stat. + Bkg. & - & {0} & {1} \\\\\\hline".format(round(statBkgUnc,4), round(100*statBkgUnc/nomFrac,1))
#      print "\n\\hline\\hline"
#      print "Total &  & {0} & {1} \\\\\\hline\\hline".format(round(math.sqrt(sum_el+sum_mu+sum_JES+sum_JER+sum_btag+sum_met+sum_rad+sum_ps+sum_cr+sum_ue+sum_mcgen+statBkgUnc*statBkgUnc),4), round(100*round(math.sqrt(sum_el+sum_mu+sum_JES+sum_JER+sum_btag+sum_met+sum_rad+sum_ps+sum_cr+sum_ue+sum_mcgen+statBkgUnc*statBkgUnc),4)/nomFrac,1))
   print "\end{tabular}"
   print "\end{table}"
   return






def makeFullSensitivityTable4 ( whichFraction, lepDict, lep2Dict, hadDict, had2Dict, method=''):
   """This function takes one dictionary as input and outputs the correctly formatted final table with sources broken down by type and categorized byby 0/L/R fraction passed from whichFraction. Additionally provides % effect of giving systematic w.r.t. expected value of helicity fraction"""

   fraction = -1
   nomFrac = 0

   if whichFraction == 'F0':
      fraction = 0
      nomFrac = 0.697

   if whichFraction == 'FL':
      fraction = 1
      nomFrac = 0.301

   if whichFraction == 'FR':
      fraction = 2
      nomFrac = 0.002

   print "$$$$$$$$$$$$$", statBkgUnc[0], statBkgUnc[1], statBkgUnc[2], statBkgUnc[3]
   print "$$$$$$$$$$$$$", TemplateStatUnc[0], TemplateStatUnc[1], TemplateStatUnc[2], TemplateStatUnc[3]
   lepTotals = calcTotals  (lepDict,  whichFraction, statBkgUnc[0], PDFunc[0], TopMassUnc[0], TemplateStatUnc[0])
   lep2Totals = calcTotals (lep2Dict, whichFraction, statBkgUnc[1], PDFunc[1], TopMassUnc[1], TemplateStatUnc[1])
   hadTotals = calcTotals  (hadDict,  whichFraction, statBkgUnc[2], PDFunc[2], TopMassUnc[2], TemplateStatUnc[2])
   had2Totals = calcTotals (had2Dict, whichFraction, statBkgUnc[3], PDFunc[3], TopMassUnc[3], TemplateStatUnc[3])

   print "\n### el+mu 2incl {0} ###".format(whichFraction)
   print "\\begin{table}[h!]"
   print "\centering"
   print "\\begin{tabular}{lccccc}"
   print "\hline\hline"
   if whichFraction == 'F0':
      print "\multicolumn{6}{c}{\\fo}\\\\"
   elif whichFraction == 'FL':
      print "\multicolumn{6}{c}{\\fl}\\\\"
   elif whichFraction == 'FR':
      print "\multicolumn{6}{c}{\\fr}\\\\"
   print "\hline"
   print "Systematic uncertainty & N$_{syst}$ & Lep 2incl & Lep 1excl+2incl & Had 2incl & Had 1excl+2incl \\\\\\hline"
   print "\multicolumn{5}{c}{Reconstructed Objects} \\\\\\hline"
   if method == 'alt':
      
      printFourPlusMinusLine('Muon', 6, hadTotals[0], lepTotals[1], lepTotals[2], lep2Totals[1], lep2Totals[2], hadTotals[1], hadTotals[2], had2Totals[1], had2Totals[2])
      printFourPlusMinusLine('Electron', 5, hadTotals[3], lepTotals[4], lepTotals[5], lep2Totals[4], lep2Totals[5], hadTotals[4], hadTotals[5], had2Totals[4], had2Totals[5])
      printFourPlusMinusLine('JES', 26, hadTotals[6], lepTotals[7], lepTotals[8], lep2Totals[7], lep2Totals[8], hadTotals[7], hadTotals[8], had2Totals[7], had2Totals[8])
      printFourPlusMinusLine('JER', 11, hadTotals[9], lepTotals[10], lepTotals[11], lep2Totals[10], lep2Totals[11], hadTotals[10], hadTotals[11], had2Totals[10], had2Totals[11])
      printFourPlusMinusLine('JVF', 1, hadTotals[12], lepTotals[13], lepTotals[14], lep2Totals[13], lep2Totals[14], hadTotals[13], hadTotals[14], had2Totals[13], had2Totals[14])
      printFourPlusMinusLine('\\bt tagging', 3, 3, lepTotals[18], lepTotals[19], lep2Totals[18], lep2Totals[19], hadTotals[18], hadTotals[19], had2Totals[18], had2Totals[19])
      print "\n\\hline\\hline"
      printFourPlusMinusLine('Sum of Reco Objects', ' ',' ', lepTotals[23], lepTotals[24], lep2Totals[23], lep2Totals[24], hadTotals[23], hadTotals[24], had2Totals[23], had2Totals[24])
      print "\n\\hline"
      print "\multicolumn{5}{c}{Modeling} \\\\\\hline"
      printFourPlusMinusLine('Radiation, radLo', 1, hadTotals[25], lepTotals[26], lepTotals[28], lep2Totals[26], lep2Totals[28], hadTotals[26], hadTotals[28], had2Totals[26], had2Totals[28])
      printFourPlusMinusLine('Radiation, radHi', 1, hadTotals[25], lepTotals[27], lepTotals[28], lep2Totals[27], lep2Totals[28], hadTotals[27], hadTotals[28], had2Totals[27], had2Totals[28])
      printFourPlusMinusLine('Parton Shower', 1, hadTotals[29], lepTotals[30], lepTotals[31], lep2Totals[30], lep2Totals[31], hadTotals[30], hadTotals[31], had2Totals[30], had2Totals[31])
      printFourPlusMinusLine('ME Generator', 1, hadTotals[38], lepTotals[39], lepTotals[40], lep2Totals[39], lep2Totals[40], hadTotals[39], hadTotals[40], had2Totals[39], had2Totals[40])
      printFourPlusMinusLine('PDF', ' ', ' ', lepTotals[45], lepTotals[45], lep2Totals[45], lep2Totals[45], hadTotals[45], hadTotals[45], had2Totals[45], had2Totals[45])
      printFourPlusMinusLine('Top Mass', ' ', ' ', lepTotals[46], lepTotals[46], lep2Totals[46], lep2Totals[46], hadTotals[46], hadTotals[46], had2Totals[46], had2Totals[46])
      print "\n\\hline\\hline"  
      printFourPlusMinusLine('Sum of Modeling', ' ', ' ', lepTotals[41], lepTotals[42], lep2Totals[41], lep2Totals[42], hadTotals[41], hadTotals[42], had2Totals[41], had2Totals[42])
      print "\n\\hline"
      print "\multicolumn{5}{c}{Method Uncertainty} \\\\\\hline"
      printFourPlusMinusLine('Template Statistics', ' ', ' ', lepTotals[47], lepTotals[47], lep2Totals[47], lep2Totals[47], hadTotals[47], hadTotals[47], had2Totals[47], had2Totals[47])
      print "\n\\hline\\hline"
      printFourPlusMinusLine('Total Syst.', ' ', ' ', lepTotals[23] + lepTotals[41] + lepTotals[47], lepTotals[24] + lepTotals[42] + lepTotals[47], lep2Totals[23] + lep2Totals[41] + lep2Totals[47], lep2Totals[24] + lep2Totals[42] + lep2Totals[47], hadTotals[23] + hadTotals[41] + hadTotals[47], hadTotals[24] + hadTotals[42] + hadTotals[47], had2Totals[23] + had2Totals[41] + had2Totals[47], had2Totals[24] + had2Totals[42] + had2Totals[47])
      print "Stat. + Bkg. & - & {0} & {1} & {2} & {3} \\\\\\hline".format(round(lepTotals[44],4), round(hadTotals[44],4), round(hadTotals[44],4), round(had2Totals[44],4))
      print "\n\\hline\\hline"

   print "\end{tabular}"
   print "\end{table}"
   return





def calcTotals(inDict, whichFraction, statUnc_t, PDFUnc_t, TopMassUnc_t, TemplateStatUnc_t):
   totals=[]

   #print "Q#$%#%@$: ",sumMethod
   fraction = -1
   nomFrac = 0
   #statBkgUnc = 0;
   if whichFraction == 'F0':
      fraction = 0
      nomFrac = 0.697
      #statBkgUnc = 0.0078
   if whichFraction == 'FL':
      fraction = 1
      nomFrac = 0.301
      #statBkgUnc = 0.0050
   if whichFraction == 'FR':
      fraction = 2
      nomFrac = 0.002
      #statBkgUnc = 0.0041

   sum_el=0
   sum_mu=0
   sum_JES=0
   sum_JER=0
   sum_JEff=0
   sum_JVF=0
   sum_btag=0
   sum_met=0
   sum_ps=0
   sum_ue=0
   sum_mcgen=0
   sum_cr=0
   sum_rad=0

   sum_up_el=0
   sum_up_mu=0
   sum_up_JES=0
   sum_up_JER=0
   sum_up_JEff=0
   sum_up_JVF=0
   sum_up_btag=0
   sum_up_met=0
   sum_up_ps=0
   sum_up_ue=0
   sum_up_mcgen=0
   sum_up_cr=0
   sum_up_rad=0

   sum_down_el=0
   sum_down_mu=0
   sum_down_JES=0
   sum_down_JER=0
   sum_down_JEff=0
   sum_down_JVF=0
   sum_down_btag=0
   sum_down_met=0
   sum_down_ps=0
   sum_down_ue=0
   sum_down_mcgen=0
   sum_down_cr=0
   sum_down_rad=0

   sum_hi_rad=0
   sum_lo_rad=0

   n_el=0
   n_mu=0
   n_JES=0
   n_JER=0
   n_JEff=0
   n_JVF=0
   n_btag=0
   n_met=0
   n_ps=0
   n_ue=0
   n_mcgen=0
   n_cr=0
   n_rad=0

   for key, value in sorted(inDict.items()):
      stepUnc = [0, 0, 0, 0]
      radLoHi = [0, 0]
      plusMinusUnc = [0, 0]
      #print key
      for updown, value2 in sorted(inDict[key].items(), reverse=True):
         #print key, updown
         
         # for the plus/minus groupings
         #if float(inDict[key][updown][fraction]) > 0:
         if '_up' in updown or 'Up' in updown or 'UP' in updown:
            #print inDict[key]
            plusMinusUnc[0] = plusMinusUnc[0] + float(inDict[key][updown][fraction])**2
         else:
            plusMinusUnc[1] = plusMinusUnc[1] + float(inDict[key][updown][fraction])**2

         # for the 'use-the-largest' algorithm
         if 'jes' in updown and '_up' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'BTAG' in updown and 'Up' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'ELE' in updown and 'UP' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'MUON' in updown and 'UP' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
            
         if 'jes' in updown and '_down' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'BTAG' in updown and 'Down' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'ELE' in updown and 'DOWN' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'MUON' in updown and 'DOWN' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         #if 'JER' in updown:
         #   stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'jer_NP0' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'jer_Noise' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'jer_DataMC' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         #if 'JEff' in updown:
         #   stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'jvf' in updown:
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         #if 'NLO' in updown:
         if key == 'ME Generator':
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         #if '(PS)' in updown:
         if key == 'Parton Shower':
            stepUnc[0] = stepUnc[0] + float(inDict[key][updown][fraction])**2
         if 'mu=' in updown:
            print 'stepUnc[0]: ',stepUnc[0],'current val: ',inDict[key][updown][fraction]
            if stepUnc[0] == 0:
               stepUnc[0] = float(inDict[key][updown][fraction])**2
            elif math.sqrt(stepUnc[0]) < abs(float(inDict[key][updown][fraction])):
               stepUnc[0] = float(inDict[key][updown][fraction])**2

         if key == 'Radiation':
            #print 'key: ',key,'\tupdown: ',updown, '\tvalues: ',inDict[key][updown][fraction]
            #print 'stepUnc[0]: ',stepUnc[0],'current val: ',inDict[key][updown][fraction]
            if 'mu=0.' in updown:
               radLoHi[0] = float(inDict[key][updown][fraction])
            else:
               radLoHi[1] = float(inDict[key][updown][fraction])

      #print key, round(math.sqrt(stepUnc[0]),4)
      #print 'stepUnc[0]: ',stepUnc[0],'\t plusMinusUnc[0]: ',plusMinusUnc[0],'\t plusMinusUnc[1]: ',plusMinusUnc[1]

      if "jes_" in key:
         sum_JES = sum_JES + stepUnc[0]
         sum_up_JES = sum_up_JES + plusMinusUnc[0]
         sum_down_JES = sum_down_JES + plusMinusUnc[1]
         n_JES = n_JES+1
      if "ELE" in key:
         sum_el = sum_el + stepUnc[0]
         sum_up_el = sum_up_el + plusMinusUnc[0]
         sum_down_el = sum_down_el + plusMinusUnc[1]
         n_el = n_el + 1
      if "MUON" in key:
         sum_mu = sum_mu + stepUnc[0]
         sum_up_mu = sum_up_mu + plusMinusUnc[0]
         sum_down_mu = sum_down_mu + plusMinusUnc[1]
         n_mu = n_mu +1
      if "BTAG" in key:
         sum_btag = sum_btag + stepUnc[0]
         sum_up_btag = sum_up_btag + plusMinusUnc[0]
         sum_down_btag = sum_down_btag + plusMinusUnc[1]
         n_btag = n_btag + 1
      if "JER" in key or 'jer_Noise' in key or 'jer_DataMC' in key or 'jer' in key:
         sum_JER = sum_JER + stepUnc[0]
         sum_up_JER = sum_up_JER + plusMinusUnc[0]
         sum_down_JER = sum_down_JER + plusMinusUnc[1]
         n_JER = n_JER +1
      #if "JEff" in key:
      #   #print "HHHEEAAAAAY"
      #   sum_JEff = sum_JEff + stepUnc[0]
      #   sum_up_JEff = sum_up_JEff + plusMinusUnc[0]
      #   sum_down_JEff = sum_down_JEff + plusMinusUnc[1]
      #   n_JEff = n_JEff +1
      if "jvf" in key:
         sum_JVF = sum_JVF + stepUnc[0]
         sum_up_JVF = sum_up_JVF + plusMinusUnc[0]
         sum_down_JVF = sum_down_JVF + plusMinusUnc[1]
         n_JVF = n_JVF +1
      if "Radiation" in key:
         sum_rad = sum_rad + stepUnc[0]
         sum_lo_rad = sum_lo_rad + radLoHi[0]
         sum_hi_rad = sum_hi_rad + radLoHi[1]      
         sum_up_rad = sum_up_rad + plusMinusUnc[0]
         sum_down_rad = sum_down_rad + plusMinusUnc[1]
         n_rad = n_rad +1
      if "ME Generator" in key:
         sum_mcgen = sum_mcgen + stepUnc[0]
         sum_up_mcgen = sum_up_mcgen + plusMinusUnc[0]
         sum_down_mcgen = sum_down_mcgen + plusMinusUnc[1]
         n_mcgen = n_mcgen +1
      if "Parton Shower" in key:
         sum_ps = sum_ps + stepUnc[0]
         sum_up_ps = sum_up_ps + plusMinusUnc[0]
         sum_down_ps = sum_down_ps + plusMinusUnc[1]
         n_ps = n_ps +1

   sum_keep_rad=0
   # A. Radiation Prescription 1: Take largest as full uncertainty
   if abs(sum_lo_rad) > abs(sum_hi_rad):
      sum_keep_rad = float(sum_lo_rad)**2
   else:
      sum_keep_rad = float(sum_hi_rad)**2
   # B. Radiation Prescription 2: Take uncertainty as 0.5*|hi - lo|
   #sum_keep_rad = float( 0.5 * abs(sum_hi_rad-sum_lo_rad) )**2   

   # C. Fully symmetrize mcgen and parton shower
   if abs(sum_up_mcgen) > abs(sum_down_mcgen):
      sum_down_mcgen = sum_up_mcgen
   else:
      sum_up_mcgen = sum_down_mcgen

   if abs(sum_up_ps) > abs(sum_down_ps):
      sum_down_ps = sum_up_ps
   else:
      sum_up_ps = sum_down_ps
   # D. Fully symmetrize one-sided uncertainties (JEff, JER, XX)
   #if abs(sum_up_JER) > abs(sum_down_JER):
   #   sum_down_JER = sum_up_JER
   #else:
   #   sum_up_JER = sum_down_JER

   if abs(sum_up_JEff) > abs(sum_down_JEff):
      sum_down_JEff = sum_up_JEff
   else:
      sum_up_JEff = sum_down_JEff


   sum_mod_up = sum_keep_rad + sum_up_ps + sum_up_cr + sum_up_ue + sum_up_mcgen + PDFUnc_t**2 + TopMassUnc_t**2
   sum_mod_down = sum_keep_rad + sum_down_ps + sum_down_cr + sum_down_ue + sum_down_mcgen + PDFUnc_t**2 + TopMassUnc_t**2
   sum_det_up = sum_up_mu + sum_up_el + sum_up_JES + sum_up_JER + sum_up_JVF + sum_up_btag + sum_up_met
   sum_det_down = sum_down_mu + sum_down_el + sum_down_JES + sum_down_JER + sum_down_JVF + sum_down_btag + sum_down_met

   totErr=sum_det_up + sum_det_down + sum_mod_up + sum_mod_down + statUnc_t**2 + TemplateStatUnc_t**2

   totals.append(n_mu)
   totals.append(sum_up_mu)
   totals.append(sum_down_mu)
   totals.append(n_el)
   totals.append(sum_up_el)
   totals.append(sum_down_el)
   totals.append(n_JES)
   totals.append(sum_up_JES)
   totals.append(sum_down_JES)
   totals.append(n_JER)
   totals.append(sum_up_JER)
   totals.append(sum_down_JER)
   totals.append(n_JVF)
   totals.append(sum_up_JVF)
   totals.append(sum_down_JVF)
   totals.append(n_JEff)
   totals.append(sum_up_JEff)
   totals.append(sum_down_JEff)
   totals.append(sum_up_btag)
   totals.append(sum_down_btag)
   totals.append(n_met)
   totals.append(sum_up_met)
   totals.append(sum_down_met)
   totals.append(sum_det_up) #23
   totals.append(sum_det_down) #24
   totals.append(n_rad)
   totals.append(sum_lo_rad)
   totals.append(sum_hi_rad)
   totals.append(sum_keep_rad)
   totals.append(n_ps)
   totals.append(sum_up_ps)
   totals.append(sum_down_ps)
   totals.append(n_cr)
   totals.append(sum_up_cr)
   totals.append(sum_down_cr)
   totals.append(n_ue)
   totals.append(sum_up_ue)
   totals.append(sum_down_ue)
   totals.append(n_mcgen)
   totals.append(sum_up_mcgen)
   totals.append(sum_down_mcgen)
   totals.append(sum_mod_up) #41
   totals.append(sum_mod_down) #42
   totals.append(totErr)
   totals.append(statUnc_t) #44
   totals.append(PDFUnc_t**2)
   totals.append(TopMassUnc_t**2)
   totals.append(TemplateStatUnc_t**2) # 47
  

   return totals

def printPlusMinusLine(systCat, nTot, nKept, varUp, varDown, nomFrac, totErr):
   
   #if (nTot == 0):
   #   print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& \multirow{2}{*}{-} & \multirow{2}{*}{-}\\\\'
   #   return

   sVars=''
   if nTot==' ':
      sVars='-'
   else:
      sVars='{0}({1})'.format(nTot,nKept)

   if 'Radiation' not in systCat:
      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& +{0} & {1}\\\\'.format(round(math.sqrt(varUp),4), round(100*varUp/totErr,1))
      print '                      &                       & -{0} & {1}\\\\\\hline'.format(round(math.sqrt(varDown),4), round(100*varDown/totErr,1))
   else:
      relUnc=''
      print 'VarDown: ',math.sqrt(varDown),'\tVarUp: ',varUp
      if abs(float(math.sqrt(varDown))) == abs(float(varUp)):
         relUnc = round(100*abs(varUp**2)/totErr,1)
      else:
         relUnc = '*'
      if 'radLo' in systCat:
         print '\multirow{2}{*}{Radiation} & radLo & '+'{0} & {1}\\\\'.format(round(varUp,4), relUnc)
      if 'radHi' in systCat:
         print '                           & radHi & '+'{0} & {1}\\\\'.format(round(varUp,4), relUnc)

def printThreeTotalLine(systCat, lepUp, lepDown, hadUp, hadDown, lephadUp, lephadDown):
   
   threshold = 0.015
   
   if abs(lepUp - lepDown) < threshold and abs(hadUp- hadDown) < threshold and abs(lephadUp- lephadDown) < threshold:
      if lepUp - lepDown > 0:
         lepBig = round(math.sqrt(lepUp),4)
      else:
         lepBig = round(math.sqrt(lepDown),4)
      if hadUp - hadDown > 0:
         hadBig = round(math.sqrt(hadUp),4)
      else:
         hadBig = round(math.sqrt(hadDown),4)
      if lephadUp - lephadDown > 0:
         lephadBig = round(math.sqrt(lephadUp),4)
      else:
         lephadBig = round(math.sqrt(lephadDown),4)

      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'}'+'  &  & \multirow{2}{*}{'+'{0}'.format(lepBig)+'} & \multirow{2}{*}{'+'{0}'.format(hadBig)+'}  & \multirow{2}{*}{'+'{0}'.format(lephadBig)+'}  \\\\  \\\\ \\hline'
   else:
      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& {0} & {1} & {2}\\\\'.format( lepUp, hadUp, lephadUp)
      print '                      &                       & {0} & {1} & {2}\\\\\\hline'.format(lepDown, hadDown, lephadDown)

  
def printThreePlusMinusLine(systCat, nTot, nKept, lepUp, lepDown, hadUp, hadDown, lephadUp, lephadDown):
   
   #if (nTot == 0):
   #   print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& \multirow{2}{*}{-} & \multirow{2}{*}{-}\\\\'
   #   return

   sVars=''
   if nTot==' ':
      sVars='-'
   else:
      sVars='{0}({1})'.format(nTot,nKept)

   # empty uncertainty
   if nKept == 0:
      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& {0} & {1} & {2}\\\\'.format('-', '-', '-')
      print '                      &                       & {0} & {1} & {2}\\\\\\hline'.format('-', '-', '-')
   #elif systCat == 'Total Syst.':
   #   if abs(lepUp - lepDown) < 0.015 && abs(hadpUp - hadDown) <0.015 && abs(lephadUp - lephadDown) <0.015:
   #      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& {0} & {1} & {2}\\\\'.format( lepUp, hadUp, lephadUp)
   #   print '                      &                       & {0} & {1} & {2}\\\\\\hline'.format(lepDown, hadDown, lephadDown)
   elif 'Radiation' not in systCat:
      lepUp = '+'+str(round(math.sqrt(lepUp),4)) if abs(round(math.sqrt(lepUp),4)) != 0.0 else '<0.0001'
      hadUp = '+'+str(round(math.sqrt(hadUp),4)) if abs(round(math.sqrt(hadUp),4)) != 0.0 else '<0.0001'
      lephadUp = '+'+str(round(math.sqrt(lephadUp),4)) if abs(round(math.sqrt(lephadUp),4)) != 0.0 else '<0.0001'
      lepDown = '-'+str(round(math.sqrt(lepDown),4)) if abs(round(math.sqrt(lepDown),4)) != 0.0 else '<0.0001'
      hadDown = '-'+str(round(math.sqrt(hadDown),4)) if abs(round(math.sqrt(hadDown),4)) != 0.0 else '<0.0001'
      lephadDown = '-'+str(round(math.sqrt(lephadDown),4)) if abs(round(math.sqrt(lephadDown),4)) != 0.0 else '<0.0001'
      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& {0} & {1} & {2}\\\\'.format( lepUp, hadUp, lephadUp)
      print '                      &                       & {0} & {1} & {2}\\\\\\hline'.format(lepDown, hadDown, lephadDown)
   else:
      relUnc=''
      #print 'VarDown: ',math.sqrt(varDown),'\tVarUp: ',varUp
      #if abs(float(math.sqrt(varDown))) == abs(float(varUp)):
      #   relUnc = round(100*abs(varUp**2)/totErr,1)
      #else:
      #   relUnc = '*'
      if 'radLo' in systCat:
         print '\multirow{2}{*}{Radiation} & radLo & '+'{0} & {1} & {2}\\\\'.format(round(lepUp,4), round(hadUp,4), round(lephadUp,4))
      if 'radHi' in systCat:
         print '                           & radHi & '+'{0} & {1} & {2}\\\\ \\hline'.format(round(lepUp,4), round(hadUp,4), round(lephadUp,4))





def printFourPlusMinusLine(systCat, nTot, nKept, lepUp, lepDown, lep2Up, lep2Down, hadUp, hadDown, had2Up, had2Down):
   
   sVars=''
   if nTot==' ':
      sVars='-'
   else:
      sVars='{0}({1})'.format(nTot,nKept)

   # empty uncertainty
   if nKept == 0:
      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& {0} & {1} & {2} & {3}\\\\'.format('-', '-', '-', '-')
      print '                      &                       & {0} & {1} & {2} & {3}\\\\\\hline'.format('-', '-', '-', '-')

   elif 'Radiation' not in systCat:
      lepUp = '+'+str(round(math.sqrt(lepUp),4)) if abs(round(math.sqrt(lepUp),4)) != 0.0 else '<0.0001'
      lep2Up = '+'+str(round(math.sqrt(lep2Up),4)) if abs(round(math.sqrt(lep2Up),4)) != 0.0 else '<0.0001'
      hadUp = '+'+str(round(math.sqrt(hadUp),4)) if abs(round(math.sqrt(hadUp),4)) != 0.0 else '<0.0001'
      had2Up = '+'+str(round(math.sqrt(had2Up),4)) if abs(round(math.sqrt(had2Up),4)) != 0.0 else '<0.0001'
      lepDown = '-'+str(round(math.sqrt(lepDown),4)) if abs(round(math.sqrt(lepDown),4)) != 0.0 else '<0.0001'
      lep2Down = '-'+str(round(math.sqrt(lep2Down),4)) if abs(round(math.sqrt(lep2Down),4)) != 0.0 else '<0.0001'
      hadDown = '-'+str(round(math.sqrt(hadDown),4)) if abs(round(math.sqrt(hadDown),4)) != 0.0 else '<0.0001'
      had2Down = '-'+str(round(math.sqrt(had2Down),4)) if abs(round(math.sqrt(had2Down),4)) != 0.0 else '<0.0001'
      print '\multirow{2}{*}{'+'{0}'.format(systCat)+'} & \multirow{2}{*}{'+'{0}'.format(sVars)+'} '+'& {0} & {1} & {2} & {3} \\\\'.format( lepUp, lep2Up, hadUp, had2Up)
      print '                      &                       & {0} & {1} & {2} & {3}\\\\\\hline'.format(lepDown, lep2Down, hadDown, had2Down)
   else:
      relUnc=''
      if 'radLo' in systCat:
         print '\multirow{2}{*}{Radiation} & radLo & '+'{0} & {1} & {2} & {3}\\\\'.format(round(lepUp,4), round(lep2Up,4), round(hadUp,4), round(had2Up,4))
      if 'radHi' in systCat:
         print '                           & radHi & '+'{0} & {1} & {2} & {3}\\\\ \\hline'.format(round(lepUp,4), round(lep2Up,4), round(hadUp,4), round(had2Up,4))




def makeFullTable( whichFraction, lepBTagDict, hadBTagDict, lephadBTagDict):
   """This function takes three dictionaries as input and prints out a table corresponding to the 0/L/R fraction passed from whichFraction"""

   fraction = -1
   nomFrac = 0
   #statBkgUnc = 0;
   #PDFunc = 0;
   #TopMassUnc = 0;
   if whichFraction == 'F0':
      fraction = 0
      nomFrac = 0.697
      #statBkgUnc = [0.010, 0.012, 0.007]
      #PDFunc = [0.003, 0.001, 0.002]
      #TopMassUnc = [0.002, 0.003, 0.001]
      #TemplateStatUnc = [0.007, 0.007, 0.005]

   if whichFraction == 'FL':
      fraction = 1
      nomFrac = 0.301
      #statBkgUnc = [0.007, 0.022, 0.005]
      #PDFunc = [0.004, 0.002, 0.003]
      #TopMassUnc = [0.005, 0.009, 0.003]
      #TemplateStatUnc = [0.005, 0.016, 0.003]

   if whichFraction == 'FR':
      fraction = 2
      nomFrac = 0.002
      #statBkgUnc = [0.005, 0.026, 0.004]
      #PDFunc = [0.001, 0.002, 0.002]
      #TopMassUnc = [0.003, 0.006, 0.004]
      #TemplateStatUnc = [0.003, 0.016, 0.003]

   print '$$$$$', len(PDFunc)
   print len(PDFunc), len(TopMassUnc), len(TemplateStatUnc)
   lepTotals = calcTotals (lepBTagDict, whichFraction, statBkgUnc[0], PDFunc[0], TopMassUnc[0], TemplateStatUnc[0])
   hadTotals = calcTotals (hadBTagDict, whichFraction, statBkgUnc[1], PDFunc[1], TopMassUnc[1], TemplateStatUnc[1])
   lephadTotals = calcTotals (lephadBTagDict, whichFraction, statBkgUnc[2], PDFunc[2], TopMassUnc[2], TemplateStatUnc[2])


   #fraction = -1
   #if whichFraction == 'F0':
   #   fraction = 0
   #if whichFraction == 'FL':
   #   fraction = 1
   #if whichFraction == 'FR':
   #   fraction = 2
      
   print "\n### el+mu 2incl {0} ###".format(whichFraction)
   print "\\begin{table}[h!]"
   print "\centering"
   print "\\begin{tabular}{lcccc}"
   print "\hline\hline"
   if whichFraction == 'F0':
      print "\multicolumn{5}{c}{\\fo}\\\\\hline"
   elif whichFraction == 'FL':
      print "\multicolumn{5}{c}{\\fl}\\\\\hline" 
   elif whichFraction == 'FR':
      print "\multicolumn{5}{c}{\\fr}\\\\\hline"
   
   
   print "Systematic uncertainty & Up/Down & Leptonic 2incl & Had 1excl+2incl & Lep+Had 1excl+2incl \\\\\\hline"
   # ~~~~~~~ Modelling systematics first ~~~~~~
   print '\multicolumn{5}{c}{Modeling} \\\\ \\hline'
   #print len(lephadBTagDict.items())
   for key, value in sorted(lephadBTagDict.items()):
      #print len(lephadBTagDict[key]), key, lephadBTagDict[key]
      #print key, mu2incl[key]
      for updown, value2 in sorted(lephadBTagDict[key].items(), reverse=True):
         # single-sided uncertainties
         #if 'Parton Shower' in updown or 'ME Generator' in updown or 'Herwig' in updown:
         if 'PartonShower' in updown or 'MEgenerator' in updown or 'Herwig' in updown:
            #print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'  &  & \multirow{2}{*}{'+'{0}'.format(round(float(lephadBTagDict[key][updown][fraction]),4))+'}      \\\\ \\hline'
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'  &  & \multirow{2}{*}{'+'{0}'.format(round(float(lepBTagDict[key][updown][fraction]),4))+'} & \multirow{2}{*}{'+'{0}'.format(round(float(hadBTagDict[key][updown][fraction]),4))+'} & \multirow{2}{*}{'+'{0}'.format(round(float(lephadBTagDict[key][updown][fraction]),4))+'}   \\\\ \\\\ \\hline'
         #radiation
         if 'mu=2' in updown:
            #print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'      & radHi   &     {0} & {1} & {2}     \\\\'.format(round(float(lepBTagDict[key][updown][fraction]),4), round(float(hadBTagDict[key][updown][fraction]),4), round(float(lephadBTagDict[key][updown][fraction]),4))
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'      & radHi   &     {0}     &     {1}     &     {2}     \\\\'.format(round(float(lepBTagDict[key][updown][fraction]),4), round(float(hadBTagDict[key][updown][fraction]),4), round(float(lephadBTagDict[key][updown][fraction]),4))
         if 'mu=0' in updown:
            #print '                          & radLo &     {0}         \\\\ \\hline'.format(round(float(lephadBTagDict[key][updown][fraction]),4))
            print '                          & radLo &     {0}     &     {1}     &     {2}         \\\\ \\hline'.format(round(float(lepBTagDict[key][updown][fraction]),4), round(float(hadBTagDict[key][updown][fraction]),4), round(float(lephadBTagDict[key][updown][fraction]),4))
   

   # ~~~~~~~ Object systematics second ~~~~~~
   print '\multicolumn{5}{c}{Reconstructed Objects} \\\\ \\hline'
   for key, value in sorted(lephadBTagDict.items()):
      #print len(lephadBTagDict[key]), key, lephadBTagDict[key]
      #print key, mu2incl[key]
      for updown, value2 in sorted(lephadBTagDict[key].items(), reverse=True):
         #if 'jer_NP8' not in key:
         if '_up' in updown or 'Up' in updown or 'UP' in updown or 'scaleup' in updown:
            #print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'      & up   &     {0}      \\\\'.format(round(float(lephadBTagDict[key][updown][fraction]),4))
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'      & up   &     {0}     &     {1}     &     {2}      \\\\'.format(round(float(lepBTagDict[key][updown][fraction]),4), round(float(hadBTagDict[key][updown][fraction]),4), round(float(lephadBTagDict[key][updown][fraction]),4))
         if '_down' in updown or 'Down' in updown or 'DOWN' in updown or 'scaledown' in updown:
            #print '                                       & down &     {0}           \\\\ \\hline'.format(round(float(lephadBTagDict[key][updown][fraction]),4))
            print '                                       & down &     {0}     &     {1}     &     {2}       \\\\ \\hline'.format(round(float(lepBTagDict[key][updown][fraction]),4), round(float(hadBTagDict[key][updown][fraction]),4), round(float(lephadBTagDict[key][updown][fraction]),4))
         # single-sided uncertainties
         #if 'JER' in updown or 'JEff' in updown or 'jer_DataMC' in updown or 'jer_Noise' in updown:
         #if 'JER' in updown or 'jer_DataMC' in updown or 'jer_Noise' in updown:
         #if 'JEff' in updown or 'jer_DataMC' in updown or 'jer_Noise' in updown:
         if 'jer_DataMC' in updown or 'jer_Noise' in updown:
            #print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'  &  & \multirow{2}{*}{'+'{0}'.format(round(float(lephadBTagDict[key][updown][fraction]),4))+'}      \\\\ \\hline'
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'  &  & \multirow{2}{*}{'+'{0}'.format(round(float(lepBTagDict[key][updown][fraction]),4))+'} & \multirow{2}{*}{'+'{0}'.format(round(float(hadBTagDict[key][updown][fraction]),4))+'}  & \multirow{2}{*}{'+'{0}'.format(round(float(lephadBTagDict[key][updown][fraction]),4))+'}  \\\\  \\\\ \\hline'
   
   #print len(mu2dict)
   print "\n\\hline\\hline"
   printThreePlusMinusLine('JER', 11, lephadTotals[9], lepTotals[10], lepTotals[11], hadTotals[10], hadTotals[11], lephadTotals[10], lephadTotals[11])
   print "\n\\hline\\hline"
   printThreePlusMinusLine('Total Syst.', ' ', ' ', lepTotals[23] + lepTotals[41] + lepTotals[47], lepTotals[24] + lepTotals[42] + lepTotals[47], hadTotals[23] + hadTotals[41] + hadTotals[47], hadTotals[24] + hadTotals[42] + hadTotals[47], lephadTotals[23] + lephadTotals[41] + lephadTotals[47], lephadTotals[24] + lephadTotals[42] + lephadTotals[47])
   
   print "\end{tabular}"
   print "\end{table}"
                                                                                                       
   return


def makeSumTable( lephadBTagDict):
   """This function takes one dictionary as input and prints out a table showing the aggregate sum of the uncertainty between 0/L/R as well as each individual component"""

   print "\n### el+mu 1excl+2incl ###"
   print "\\begin{table}[h!]"
   print "\centering"
   print "\\begin{tabular}{lccccc}"
   print "\hline\hline"
   print "\multicolumn{6}{c}{Leptonic+Hadronic}\\\\\hline"
     
   
   print "Systematic uncertainty & Up/Down & \\fo & \\fl & \\fr & Sum \\\\\\hline"
   # ~~~~~~~ Modelling systematics first ~~~~~~
   print '\multicolumn{6}{c}{Modeling} \\\\ \\hline'
   #print len(lephadBTagDict.items())
   for key, value in sorted(lephadBTagDict.items()):
      #print len(lephadBTagDict[key]), key, lephadBTagDict[key]
      #print key, mu2incl[key]
      for updown, value2 in sorted(lephadBTagDict[key].items(), reverse=True):
         f0 = round(float(lephadBTagDict[key][updown][0]),4) if abs(round(float(lephadBTagDict[key][updown][0]),4)) != 0.0 else '<0.0001'
         fl = round(float(lephadBTagDict[key][updown][1]),4) if abs(round(float(lephadBTagDict[key][updown][1]),4)) != 0.0 else '<0.0001'
         fr = round(float(lephadBTagDict[key][updown][2]),4) if abs(round(float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'
         sumf = round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4) if  abs(round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'     
         # single-sided uncertainties
         if 'Parton Shower' in updown or 'ME Generator' in updown or 'Herwig' in updown:
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'  &  & \multirow{2}{*}{'+'{0}'.format(f0)+'} & \multirow{2}{*}{'+'{0}'.format(fl)+'} & \multirow{2}{*}{'+'{0}'.format(fr)+'} & \multirow{2}{*}{'+'{0}'.format(sumf)+'}    \\\\ \\hline'
         #radiation
         if 'mu=2' in updown:
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'      & radHi   &     {0}     &     {1}     &     {2}   &   {3}  \\\\'.format(f0, fl, fr, sumf)
         if 'mu=0' in updown:
            print '                          & radLo &     {0}     &     {1}     &     {2}    &   {3}  \\\\'.format(f0, fl, fr, sumf)
   

   # ~~~~~~~ Object systematics second ~~~~~~
   print '\multicolumn{6}{c}{Reconstructed Objects} \\\\ \\hline'
   for key, value in sorted(lephadBTagDict.items()):
      #print len(lephadBTagDict[key]), key, lephadBTagDict[key]
      #print key, mu2incl[key]
      for updown, value2 in sorted(lephadBTagDict[key].items(), reverse=True):
         if '_up' in updown or 'Up' in updown or 'UP' in updown or 'scaleup' in updown:
            f0 = round(float(lephadBTagDict[key][updown][0]),4) if abs(round(float(lephadBTagDict[key][updown][0]),4)) != 0.0 else '<0.0001'
            fl = round(float(lephadBTagDict[key][updown][1]),4) if abs(round(float(lephadBTagDict[key][updown][1]),4)) != 0.0 else '<0.0001'
            fr = round(float(lephadBTagDict[key][updown][2]),4) if abs(round(float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'
            sumf = round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4) if  abs(round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'     
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'      & up   &     {0}     &     {1}     &     {2}  & {3}      \\\\'.format(f0, fl, fr, sumf)
         if '_down' in updown or 'Down' in updown or 'DOWN' in updown or 'scaledown' in updown:
            f0 = round(float(lephadBTagDict[key][updown][0]),4) if abs(round(float(lephadBTagDict[key][updown][0]),4)) != 0.0 else '<0.0001'
            fl = round(float(lephadBTagDict[key][updown][1]),4) if abs(round(float(lephadBTagDict[key][updown][1]),4)) != 0.0 else '<0.0001'
            fr = round(float(lephadBTagDict[key][updown][2]),4) if abs(round(float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'
            sumf = round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4) if  abs(round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'     
            print '                                       & down &     {0}     &     {1}     &     {2}   &   {3}    \\\\ \\hline'.format(f0, fl, fr, sumf)
         # single-sided uncertainties
         #if 'JER' in updown or 'JEff' in updown:
         if 'JER' in updown :
            f0 = round(float(lephadBTagDict[key][updown][0]),4) if abs(round(float(lephadBTagDict[key][updown][0]),4)) != 0.0 else '<0.0001'
            fl = round(float(lephadBTagDict[key][updown][1]),4) if abs(round(float(lephadBTagDict[key][updown][1]),4)) != 0.0 else '<0.0001'
            fr = round(float(lephadBTagDict[key][updown][2]),4) if abs(round(float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'
            sumf = round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4) if abs(round(float(lephadBTagDict[key][updown][0]) + float(lephadBTagDict[key][updown][1]) + float(lephadBTagDict[key][updown][2]),4)) != 0.0 else '<0.0001'
            print '\multirow{2}{*}{'+'{0}'.format(key.replace('_','\\_'))+'}'+'  &  & \multirow{2}{*}{'+'{0}'.format(f0)+'} & \multirow{2}{*}{'+'{0}'.format(fl)+'}  & \multirow{2}{*}{'+'{0}'.format(fr)+'}  & \multirow{2}{*}{'+'{0}'.format(sumf) +'}  \\\\ \\hline'
   
   #print len(mu2dict)
   print "\end{tabular}"
   print "\end{table}"
                                                                                                       
   return



def loopFile( outDict, inFile):
   """This function reads an input file and fills a nested dictionary { 'sysVar: {upVar: [f0, fL, fR] }, {downVar: [f0, fL, fR] } }"""

   i = 0
   for variation in inFile:
       histo    = variation.split('&')
       upDict   = {}
       downDict = {}

       
       if '_up' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           if len(upName.split('_up')) > 1:
               #downName = upName.split('_up')[0]+'_down'+upName.split('_up')[1]
               fullName = upName.split('_up')[0]+upName.split('_up')[1]
           else:
               #downName = upName.split('_up')[0]+'_down'
               fullName = upName.split('_up')[0]
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

#print 'i: ', i, fullName, upName, up_f0, up_fL, up_fR  

       if 'Up' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           if len(upName.split('Up')) > 1:
               #downName = upName.split('_up')[0]+'_down'+upName.split('_up')[1]
               fullName = upName.split('Up')[0]+upName.split('Up')[1]
           else:
               #downName = upName.split('_up')[0]+'_down'
               fullName = upName.split('Up')[0]
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       if '_UP' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           if len(upName.split('_UP')) > 1:
               #downName = upName.split('_up')[0]+'_down'+upName.split('_up')[1]
               fullName = upName.split('_UP')[0]+upName.split('_UP')[1]
           else:
               #downName = upName.split('_up')[0]+'_down'
               fullName = upName.split('_UP')[0]
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       if 'scaleup' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           if len(upName.split('up')) > 1:
               #downName = upName.split('_up')[0]+'_down'+upName.split('_up')[1]
               fullName = upName.split('up')[0]+upName.split('up')[1]
           else:
               #downName = upName.split('_up')[0]+'_down'
               fullName = upName.split('up')[0]
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       if 'mu=2' in variation:
           upName   = variation.split('&')[0].split('\t')[0]
           fullName = 'Radiation'
           
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)


       if '_down' in variation:
           downName   = variation.split('&')[0].split('\t')[0]
           if len(downName.split('_down')) > 1:
               fullName = downName.split('_down')[0]+downName.split('_down')[1]
           else:
               fullName = downName.split('_down')[0]

           down_f0 = variation.split('&')[1].split('\t')[0]
           down_fL = variation.split('&')[2].split('\t')[0]
           down_fR = variation.split('&')[3].split(' \ q')[0]

           downList = [down_f0, down_fL, down_fR]
           downDict = { downName : downList }
           if fullName not in outDict.keys():
               outDict[fullName]=downDict
           else:
               outDict[fullName].update(downDict)

       if 'Down' in variation:
           downName   = variation.split('&')[0].split('\t')[0]
           if len(downName.split('Down')) > 1:
               fullName = downName.split('Down')[0]+downName.split('Down')[1]
           else:
               fullName = downName.split('Down')[0]

           down_f0 = variation.split('&')[1].split('\t')[0]
           down_fL = variation.split('&')[2].split('\t')[0]
           down_fR = variation.split('&')[3].split(' \ q')[0]

           downList = [down_f0, down_fL, down_fR]
           downDict = { downName : downList }
           if fullName not in outDict.keys():
               outDict[fullName]=downDict
           else:
               outDict[fullName].update(downDict)

       if '_DOWN' in variation:
           downName   = variation.split('&')[0].split('\t')[0]
           if len(downName.split('_DOWN')) > 1:
               fullName = downName.split('_DOWN')[0]+downName.split('_DOWN')[1]
           else:
               fullName = downName.split('_DOWN')[0]

           down_f0 = variation.split('&')[1].split('\t')[0]
           down_fL = variation.split('&')[2].split('\t')[0]
           down_fR = variation.split('&')[3].split(' \ q')[0]

           downList = [down_f0, down_fL, down_fR]
           downDict = { downName : downList }
           if fullName not in outDict.keys():
               outDict[fullName]=downDict
           else:
               outDict[fullName].update(downDict)

       if 'scaledown' in variation:
           downName   = variation.split('&')[0].split('\t')[0]
           if len(downName.split('down')) > 1:
               fullName = downName.split('down')[0]+downName.split('down')[1]
           else:
               fullName = downName.split('down')[0]

           down_f0 = variation.split('&')[1].split('\t')[0]
           down_fL = variation.split('&')[2].split('\t')[0]
           down_fR = variation.split('&')[3].split(' \ q')[0]

           downList = [down_f0, down_fL, down_fR]
           downDict = { downName : downList }
           if fullName not in outDict.keys():
               outDict[fullName]=downDict
           else:
               outDict[fullName].update(downDict)

       if 'mu=0.5' in variation:
           downName = variation.split('&')[0].split('\t')[0]
           fullName = 'Radiation'

           down_f0 = variation.split('&')[1].split('\t')[0]
           down_fL = variation.split('&')[2].split('\t')[0]
           down_fR = variation.split('&')[3].split(' \ q')[0]

           downList = [down_f0, down_fL, down_fR]
           downDict = { downName : downList }
           if fullName not in outDict.keys():
               outDict[fullName]=downDict
           else:
               outDict[fullName].update(downDict)
 
       # single sided uncertainties....
       if 'JER' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           fullName = upName
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       if 'JEff' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           fullName = upName
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       if 'jer_DataMC' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           fullName = upName
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       if 'jer_Noise' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           fullName = upName
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       if 'PartonShower' in variation:
           upName = variation.split('&')[0].split('\t')[0]
           fullName = 'Parton Shower'
               
           i= i + 1
           up_f0  = variation.split('&')[1].split('\t')[0]
           up_fL  = variation.split('&')[2].split('\t')[0]
           up_fR  = variation.split('&')[3].split(' \ q')[0]
           upList = [up_f0, up_fL, up_fR]
           upDict = { upName : upList }
           if fullName not in outDict.keys():
               outDict[fullName]=upDict
           else:
               outDict[fullName].update(upDict)

       #if 'MC@NLO+fHerwig' in variation:
       if 'MEgenerator' in variation:
          upName = variation.split('&')[0].split('\t')[0]
          fullName = 'ME Generator'
          
          i= i + 1
          up_f0  = variation.split('&')[1].split('\t')[0]
          up_fL  = variation.split('&')[2].split('\t')[0]
          up_fR  = variation.split('&')[3].split(' \ q')[0]
          upList = [up_f0, up_fL, up_fR]
          upDict = { upName : upList }
          if fullName not in outDict.keys():
             outDict[fullName]=upDict
          else:
             outDict[fullName].update(upDict)
             
   return 
          

def makeCovarianceMatrix( lephadBTagDict):
   """This function takes one dictionary as input and creates the covariance matrix for all systematic effects as devised by Mark Own"""

   m = 3
   covSystMatrix = [[0 for x in range(m)] for y in range(m)]
   covStatMatrix = [[0 for x in range(m)] for y in range(m)]
   covMatrix = [[0 for x in range(m)] for y in range(m)]

   for key, value in sorted(lephadBTagDict.items()):
      #print len(lephadBTagDict[key]), key, lephadBTagDict[key]
      #print key, mu2incl[key]
      for updown, value2 in sorted(lephadBTagDict[key].items(), reverse=True):
         #f0 = round(float(lephadBTagDict[key][updown][0]),4) 
         #fl = round(float(lephadBTagDict[key][updown][1]),4) 
         #fr = round(float(lephadBTagDict[key][updown][2]),4) 
         f0 = float(lephadBTagDict[key][updown][0]) 
         fl = float(lephadBTagDict[key][updown][1]) 
         fr = float(lephadBTagDict[key][updown][2])

         #print updown
         #print f0*f0,'\t',f0*fl,'\t',f0*fr
         #print fl*f0,'\t',fl*fl,'\t',fl*fr
         #print fr*f0,'\t',fr*fl,'\t',fr*fr
         
         covSystMatrix[0][0] = covSystMatrix[0][0] + f0*f0
         covSystMatrix[0][1] = covSystMatrix[0][1] + f0*fl
         covSystMatrix[0][2] = covSystMatrix[0][2] + f0*fr

         covSystMatrix[1][0] = covSystMatrix[1][0] + fl*f0
         covSystMatrix[1][1] = covSystMatrix[1][1] + fl*fl
         covSystMatrix[1][2] = covSystMatrix[1][2] + fl*fr

         covSystMatrix[2][0] = covSystMatrix[2][0] + fr*f0
         covSystMatrix[2][1] = covSystMatrix[2][1] + fr*fl
         covSystMatrix[2][2] = covSystMatrix[2][2] + fr*fr

   # *** 2. Hardcode stat matrix
   s_f0 = 6.95208e-03 # 8-channel
   s_fl = 4.67178e-03 # 8-channel
   s_fr = 3.68269e-03 # 8-channel

   covStatMatrix[0][0] = s_f0*s_f0 
   covStatMatrix[0][1] = s_f0*s_fl
   covStatMatrix[0][2] = s_f0*s_fr
   
   covStatMatrix[1][0] = s_fl*s_f0
   covStatMatrix[1][1] = s_fl*s_fl
   covStatMatrix[1][2] = s_fl*s_fr
   
   covStatMatrix[2][0] = s_fr*s_f0
   covStatMatrix[2][1] = s_fr*s_fl
   covStatMatrix[2][2] = s_fr*s_fr
   
   # *** 3. Print Results
   print "Total Systematic COV"
   printMatrix(covSystMatrix)

   print "Total Statistical COV"
   printMatrix(covStatMatrix)

   # *** 4. Calculate and print total COV matrix
   covMatrix[0][0] = covStatMatrix[0][0] + covSystMatrix[0][0]
   covMatrix[0][1] = covStatMatrix[0][1] + covSystMatrix[0][1]
   covMatrix[0][2] = covStatMatrix[0][2] + covSystMatrix[0][2]
   covMatrix[1][0] = covStatMatrix[1][0] + covSystMatrix[1][0]
   covMatrix[1][1] = covStatMatrix[1][1] + covSystMatrix[1][1]
   covMatrix[1][2] = covStatMatrix[1][2] + covSystMatrix[1][2]
   covMatrix[2][0] = covStatMatrix[2][0] + covSystMatrix[2][0]
   covMatrix[2][1] = covStatMatrix[2][1] + covSystMatrix[2][1]
   covMatrix[2][2] = covStatMatrix[2][2] + covSystMatrix[2][2]

   print "Total COV"
   printMatrix(covMatrix)

   
   return

def makeCorrelationMatrixTXT( systInfile ):
   """This function takes txt file as input and creates correlation matrix for all systematic effects as devised by Mark Owen"""
   lep2incl = True
   m = 3
   corrSystMatrix = [[0 for x in range(m)] for y in range(m)]
   corrStatMatrix = [[0 for x in range(m)] for y in range(m)]
   corrMatrix = [[0 for x in range(m)] for y in range(m)]

   txtfile = open(systInfile,'r')
   
   rad0 = 0
   radL = 0
   radR = 0
   radTot = 0

   for line in txtfile:
      #if 'loCR' not in line and 'mpiHi' not in line and 'Parton' not in line and 'ME' not in line and 'hdamp' not in line:
      if 'loCR' not in line and 'mpiHi' not in line and 'hdamp' not in line:
         #print line.split('&')[0]
      #print line.split('&')[1].split('\t')[0], line.split('&')[2].split('\t')[0], line.split('&')[3].split('\\')[0]
         sig0 = float(line.split('&')[1].split('\t')[0])
         sigL = float(line.split('&')[2].split('\t')[0])
         sigR = float(line.split('&')[3].split('\\')[0])
         
      #print nom.split('//')[1], '\t sig0:', round(sig0,4), '\t sigL:', round(sigL,4), '\t sigR:', round(sigR,4)
         f_0L = sig0*sigL
         f_0R = sig0*sigR
         f_LR = sigL*sigR
         
         corrSystMatrix[0][0] = corrSystMatrix[0][0] + sig0*sig0
         corrSystMatrix[0][1] = corrSystMatrix[0][1] + f_0L
         corrSystMatrix[0][2] = corrSystMatrix[0][2] + f_0R
         
         corrSystMatrix[1][0] = corrSystMatrix[1][0] + f_0L
         corrSystMatrix[1][1] = corrSystMatrix[1][1] + sigL*sigL
         corrSystMatrix[1][2] = corrSystMatrix[1][2] + f_LR
         
         corrSystMatrix[2][0] = corrSystMatrix[2][0] + f_0R
         corrSystMatrix[2][1] = corrSystMatrix[2][1] + f_LR
         corrSystMatrix[2][2] = corrSystMatrix[2][2] + sigR*sigR
      
      elif 'hdamp' in line:
         print line.split('&')[0]
         print "1: ",line
         sig0 = float(line.split('&')[1].split('\t')[0])
         sigL = float(line.split('&')[2].split('\t')[0])
         sigR = float(line.split('&')[3].split('\\')[0])
         if rad0 == 0 and radL == 0 and radR == 0:
            rad0 = sig0
            radL = sigL
            radR = sigR
            radTot = math.sqrt(rad0*rad0 + radL*radL +radR*radR)
         else:
            print "2: ",line
            tempTot = math.sqrt(rad0*rad0 + radL*radL +radR*radR)
            if tempTot > radTot:
               rad0 = sig0
               radL = sigL
               radR = sigR

            f_0L = rad0*radL
            f_0R = rad0*radR
            f_LR = radL*radR
            
            corrSystMatrix[0][0] = corrSystMatrix[0][0] + rad0*rad0
            corrSystMatrix[0][1] = corrSystMatrix[0][1] + f_0L
            corrSystMatrix[0][2] = corrSystMatrix[0][2] + f_0R
            
            corrSystMatrix[1][0] = corrSystMatrix[1][0] + f_0L
            corrSystMatrix[1][1] = corrSystMatrix[1][1] + radL*radL
            corrSystMatrix[1][2] = corrSystMatrix[1][2] + f_LR
            
            corrSystMatrix[2][0] = corrSystMatrix[2][0] + f_0R
            corrSystMatrix[2][1] = corrSystMatrix[2][1] + f_LR
            corrSystMatrix[2][2] = corrSystMatrix[2][2] + radR*radR
         
   # *** 2. Hardcode stat matrix
   #s_f0 = 6.95208e-03 # 8-channel
   #s_fl = 4.67178e-03 # 8-channel
   #s_fr = 3.68269e-03 # 8-channel
         
   s_f0fl = -2.94883618515706658e-05 # 8-channel
   s_f0fr = -2.10904470438445183e-05 # 8-channel
   s_flfr = 7.10647872315009303e-06 # 8-channel
   
   if lep2incl:
      s_f0fl = -9.16664135974798455e-05 # lep2incl
      s_f0fr = -6.05296831697142972e-05 # lep2incl
      s_flfr = 2.83165308800971708e-05 # lep2incl
      
   corrStatMatrix[0][0] = s_f0fl*s_f0fr/s_flfr
   corrStatMatrix[0][1] = s_f0fl
   corrStatMatrix[0][2] = s_f0fr
   
   corrStatMatrix[1][0] = s_f0fl
   corrStatMatrix[1][1] = s_f0fl*s_flfr/s_f0fr
   corrStatMatrix[1][2] = s_flfr
   
   corrStatMatrix[2][0] = s_f0fr
   corrStatMatrix[2][1] = s_flfr
   corrStatMatrix[2][2] = s_f0fr*s_flfr/s_f0fl
   
   # *** 3. Hardcode addition of top mass and template stat unc
   #  ** a. Top Mass Unc
   mt_f0 = 0.001 # 8-channel
   mt_fl = 0.003 # 8-channel
   mt_fr = 0.004 # 8-channel
   
   if lep2incl:
      mt_f0 = 0.002 # lep2incl
      mt_fl = 0.005 # lep2incl
      mt_fr = 0.003 # lep2incl
   
   corrSystMatrix[0][0] = corrSystMatrix[0][0] + mt_f0*mt_f0
   corrSystMatrix[0][1] = corrSystMatrix[0][1] + mt_f0*mt_fl
   corrSystMatrix[0][2] = corrSystMatrix[0][2] + mt_f0*mt_fr
   
   corrSystMatrix[1][0] = corrSystMatrix[1][0] + mt_fl*mt_f0
   corrSystMatrix[1][1] = corrSystMatrix[1][1] + mt_fl*mt_fl
   corrSystMatrix[1][2] = corrSystMatrix[1][2] + mt_fl*mt_fr
   
   corrSystMatrix[2][0] = corrSystMatrix[2][0] + mt_fr*mt_f0
   corrSystMatrix[2][1] = corrSystMatrix[2][1] + mt_fr*mt_fl
   corrSystMatrix[2][2] = corrSystMatrix[2][2] + mt_fr*mt_fr
   
   #  ** b. Top PDF Unc
   pdf_f0 = 0.002 # 8-channel
   pdf_fl = 0.003 # 8-channel
   pdf_fr = 0.002 # 8-channel
   
   if lep2incl:
      pdf_f0 = 0.003 # lep2incl
      pdf_fl = 0.004 # lep2incl
      pdf_fr = 0.001 # lep2incl
      
   corrSystMatrix[0][0] = corrSystMatrix[0][0] + pdf_f0*pdf_f0
   corrSystMatrix[0][1] = corrSystMatrix[0][1] + pdf_f0*pdf_fl
   corrSystMatrix[0][2] = corrSystMatrix[0][2] + pdf_f0*pdf_fr
   
   corrSystMatrix[1][0] = corrSystMatrix[1][0] + pdf_fl*pdf_f0
   corrSystMatrix[1][1] = corrSystMatrix[1][1] + pdf_fl*pdf_fl
   corrSystMatrix[1][2] = corrSystMatrix[1][2] + pdf_fl*pdf_fr
   
   corrSystMatrix[2][0] = corrSystMatrix[2][0] + pdf_fr*pdf_f0
   corrSystMatrix[2][1] = corrSystMatrix[2][1] + pdf_fr*pdf_fl
   corrSystMatrix[2][2] = corrSystMatrix[2][2] + pdf_fr*pdf_fr
   
   
   #  ** c. Template Statistics Unc
   ts_f0 = 0.005 # 8-channel
   ts_fl = 0.003 # 8-channel
   ts_fr = 0.003 # 8-channel
   
   if lep2incl:
      ts_f0 = 0.009 # lep2incl
      ts_fl = 0.006 # lep2incl
      ts_fr = 0.004 # lep2incl
      
   corrSystMatrix[0][0] = corrSystMatrix[0][0] + ts_f0*ts_f0
   corrSystMatrix[0][1] = corrSystMatrix[0][1] + ts_f0*ts_fl
   corrSystMatrix[0][2] = corrSystMatrix[0][2] + ts_f0*ts_fr
   
   corrSystMatrix[1][0] = corrSystMatrix[1][0] + ts_fl*ts_f0
   corrSystMatrix[1][1] = corrSystMatrix[1][1] + ts_fl*ts_fl
   corrSystMatrix[1][2] = corrSystMatrix[1][2] + ts_fl*ts_fr
   
   corrSystMatrix[2][0] = corrSystMatrix[2][0] + ts_fr*ts_f0
   corrSystMatrix[2][1] = corrSystMatrix[2][1] + ts_fr*ts_fl
   corrSystMatrix[2][2] = corrSystMatrix[2][2] + ts_fr*ts_fr
   
   # *** 4. Print Results
   print "Total Systematic CORR"
   printMatrix(corrSystMatrix)
   
   print "Total Statistical CORR"
   printMatrix(corrStatMatrix)
   
   # *** 5. Calculate and print total CORR matrix
   corrMatrix[0][0] = corrStatMatrix[0][0] + corrSystMatrix[0][0]
   corrMatrix[0][1] = corrStatMatrix[0][1] + corrSystMatrix[0][1]
   corrMatrix[0][2] = corrStatMatrix[0][2] + corrSystMatrix[0][2]
   corrMatrix[1][0] = corrStatMatrix[1][0] + corrSystMatrix[1][0]
   corrMatrix[1][1] = corrStatMatrix[1][1] + corrSystMatrix[1][1]
   corrMatrix[1][2] = corrStatMatrix[1][2] + corrSystMatrix[1][2]
   corrMatrix[2][0] = corrStatMatrix[2][0] + corrSystMatrix[2][0]
   corrMatrix[2][1] = corrStatMatrix[2][1] + corrSystMatrix[2][1]
   corrMatrix[2][2] = corrStatMatrix[2][2] + corrSystMatrix[2][2]
   
   print "Total CORR"
   printMatrix(corrMatrix)

   print " #### Correlation Coeffs ####"
   rho_0L = corrMatrix[0][1] / ( math.sqrt(corrMatrix[0][0]) * math.sqrt(corrMatrix[1][1]) )
   rho_0R = corrMatrix[0][2] / ( math.sqrt(corrMatrix[0][0]) * math.sqrt(corrMatrix[2][2]) )
   rho_LR = corrMatrix[1][2] / ( math.sqrt(corrMatrix[1][1]) * math.sqrt(corrMatrix[2][2]) )
   print "rho(F0, FL) = ", rho_0L
   print "rho(F0, FR) = ", rho_0R
   print "rho(FL, FR) = ", rho_LR
   
def makeCorrelationMatrix( lephadBTagDict):
   """This function takes one dictionary as input and creates the correlation matrix for all systematic effects as devised by Mark Own"""

   lep2incl = True

   #filedir = 'ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/'
   #filedir = '../syst_outputs_05Jul/JERnp/ExternalSystematicsOutput_el_mu_lephad_bTag_1excl2incl_3D_3W_lephad/'
   #filedir = 'ExternalSystematicsOutput_el_mu_bTag_1excl2incl_3D_3W_lep/'
   filedir = 'ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep/'

   m = 3
   corrSystMatrix = [[0 for x in range(m)] for y in range(m)]
   corrStatMatrix = [[0 for x in range(m)] for y in range(m)]
   corrMatrix = [[0 for x in range(m)] for y in range(m)]

   h_c0R = ROOT.TH1D("h_c0R","h_c0R",20,-0.85,-0.5)
   h_c0L = ROOT.TH1D("h_c0L","h_c0L",20,-0.90,-0.85)
   h_cLR = ROOT.TH1D("h_cLR","h_cLR",20,.2,.5)
   h_c0R.SetXTitle("c0R")
   h_c0L.SetXTitle("c0L")
   h_cLR.SetXTitle("cLR")

   c0 = ROOT.TCanvas("c0","c0",800,800)
   cL = ROOT.TCanvas("cL","cL",800,800)
   cR = ROOT.TCanvas("cR","cR",800,800)
   h_d0R = ROOT.TH1D("h_rho_0R","h_rho_0R",20,-0.1,0.1)
   h_d0L = ROOT.TH1D("h_rho_0L","h_rho_0L",20,-0.1,0.1)
   h_dLR = ROOT.TH1D("h_rho_LR","h_rho_LR",20,-0.1,0.1)

   
   for subdir, dirs, files in os.walk(filedir):
      flist = []
      fdict = {}

      for file in files:
         if 'Nominal' in file and 'root' in file:

            nom = str(filedir+"/"+file)
            #print nom
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")

            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            
            if 'AFII' not in nom:
               nomF0 = fgaus0.GetParameter(1)
               nomFL = fgausL.GetParameter(1)
               nomFR = fgausR.GetParameter(1)
            else:
               afiiF0 = fgaus0.GetParameter(1)
               afiiFL = fgausL.GetParameter(1)
               afiiFR = fgausR.GetParameter(1)

      chk0 = 0
      chkL = 0
      chkR = 0
      for file in files:
         if 'root' in file and ("CT10" not in file and "MSTW" not in file and "NNPDF" not in file and "Mass" not in file and "Nominal" not in file and 'AFII' not in file and 'mpiHi' not in file and 'loCR' not in file and 'Powheg+Pythia_' not in file and '(PS)' not in file):

            nom = str(filedir+"/"+file)
            print nom
            f = ROOT.TFile.Open(nom, "read")
            t = f.Get("EnsembleTree")
            
            ## *** Method C. Follow Mark's methodology, rho(f0,fL) = sig0*sigL
            sig0 = 0
            sigL = 0
            sigR = 0
            fgaus0 = ROOT.TF1("fgaus0","gaus")
            fgausL = ROOT.TF1("fgausL","gaus")
            fgausR = ROOT.TF1("fgausR","gaus")
            t.Fit("fgaus0","F0","","Q")
            t.Fit("fgausL","FL","","Q")
            t.Fit("fgausR","FR","","Q")
            if 'mu=' in nom:
               #print nom, 'AFII sub'
               if '2.0' in nom:
                  sigL = afiiFL - fgausL.GetParameter(1)
               if '0.5' in nom:
                  sig0 = afiiF0 - fgaus0.GetParameter(1)
                  sigR = afiiFR - fgausR.GetParameter(1)
            #print 'sig0:', sig0, 'sigL:', sigL, 'sigR:', sigR

            else:
               
               if 'PS' in nom:
                  if 'Pythia' in nom:
                     sig0 = nomF0 - fgaus0.GetParameter(1)
                  else:
                     sigL = nomFL - fgausL.GetParameter(1)
                     sigR = nomFR - fgausR.GetParameter(1)
                  print 'sig0:', sig0, 'sigL:', sigL, 'sigR:', sigR
               else:
                  sig0 = nomF0 - fgaus0.GetParameter(1)
                  sigL = nomFL - fgausL.GetParameter(1)
                  sigR = nomFR - fgausR.GetParameter(1)
               

            #print nom.split('//')[1], '\t sig0:', round(sig0,4), '\t sigL:', round(sigL,4), '\t sigR:', round(sigR,4)
            f_0L = sig0*sigL
            f_0R = sig0*sigR
            f_LR = sigL*sigR

            #print file, "C_0L:", round(f_0L,6), "C_0R:", round(f_0R,6), "C_LR:", round(f_LR,6) 
            #print file, "sig0:", round(sig0,4), "sigL:", round(sigL,4), "sigR:", round(sigR,4), "sum: ", round(sig0+sigL+sigR,4) 
            
            if 'BTAG' in nom:
               chk0 = chk0 + sig0*sig0
               chkL = chkL + sigL*sigL
               chkR = chkR + sigR*sigR


            corrSystMatrix[0][0] = corrSystMatrix[0][0] + sig0*sig0
            corrSystMatrix[0][1] = corrSystMatrix[0][1] + f_0L
            corrSystMatrix[0][2] = corrSystMatrix[0][2] + f_0R
            
            corrSystMatrix[1][0] = corrSystMatrix[1][0] + f_0L
            corrSystMatrix[1][1] = corrSystMatrix[1][1] + sigL*sigL
            corrSystMatrix[1][2] = corrSystMatrix[1][2] + f_LR
            
            corrSystMatrix[2][0] = corrSystMatrix[2][0] + f_0R
            corrSystMatrix[2][1] = corrSystMatrix[2][1] + f_LR
            corrSystMatrix[2][2] = corrSystMatrix[2][2] + sigR*sigR
            
   # *** 2. Hardcode stat matrix
   #s_f0 = 6.95208e-03 # 8-channel
   #s_fl = 4.67178e-03 # 8-channel
   #s_fr = 3.68269e-03 # 8-channel

   s_f0fl = -2.94883618515706658e-05 # 8-channel
   s_f0fr = -2.10904470438445183e-05 # 8-channel
   s_flfr = 7.10647872315009303e-06 # 8-channel

   if lep2incl:
      s_f0fl = -9.16664135974798455e-05 # lep2incl
      s_f0fr = -6.05296831697142972e-05 # lep2incl
      s_flfr = 2.83165308800971708e-05 # lep2incl
   
   corrStatMatrix[0][0] = s_f0fl*s_f0fr/s_flfr
   corrStatMatrix[0][1] = s_f0fl
   corrStatMatrix[0][2] = s_f0fr
   
   corrStatMatrix[1][0] = s_f0fl
   corrStatMatrix[1][1] = s_f0fl*s_flfr/s_f0fr
   corrStatMatrix[1][2] = s_flfr
   
   corrStatMatrix[2][0] = s_f0fr
   corrStatMatrix[2][1] = s_flfr
   corrStatMatrix[2][2] = s_f0fr*s_flfr/s_f0fl
   
   # *** 3. Hardcode addition of top mass and template stat unc
   #  ** a. Top Mass Unc
   mt_f0 = 0.001 # 8-channel
   mt_fl = 0.003 # 8-channel
   mt_fr = 0.004 # 8-channel

   if lep2incl:
      mt_f0 = 0.002 # lep2incl
      mt_fl = 0.005 # lep2incl
      mt_fr = 0.003 # lep2incl
   
   corrSystMatrix[0][0] = corrSystMatrix[0][0] + mt_f0*mt_f0
   corrSystMatrix[0][1] = corrSystMatrix[0][1] + mt_f0*mt_fl
   corrSystMatrix[0][2] = corrSystMatrix[0][2] + mt_f0*mt_fr
   
   corrSystMatrix[1][0] = corrSystMatrix[1][0] + mt_fl*mt_f0
   corrSystMatrix[1][1] = corrSystMatrix[1][1] + mt_fl*mt_fl
   corrSystMatrix[1][2] = corrSystMatrix[1][2] + mt_fl*mt_fr
   
   corrSystMatrix[2][0] = corrSystMatrix[2][0] + mt_fr*mt_f0
   corrSystMatrix[2][1] = corrSystMatrix[2][1] + mt_fr*mt_fl
   corrSystMatrix[2][2] = corrSystMatrix[2][2] + mt_fr*mt_fr

   #  ** b. Top PDF Unc
   pdf_f0 = 0.002 # 8-channel
   pdf_fl = 0.003 # 8-channel
   pdf_fr = 0.002 # 8-channel

   if lep2incl:
      pdf_f0 = 0.003 # lep2incl
      pdf_fl = 0.004 # lep2incl
      pdf_fr = 0.001 # lep2incl
   
   corrSystMatrix[0][0] = corrSystMatrix[0][0] + pdf_f0*pdf_f0
   corrSystMatrix[0][1] = corrSystMatrix[0][1] + pdf_f0*pdf_fl
   corrSystMatrix[0][2] = corrSystMatrix[0][2] + pdf_f0*pdf_fr
   
   corrSystMatrix[1][0] = corrSystMatrix[1][0] + pdf_fl*pdf_f0
   corrSystMatrix[1][1] = corrSystMatrix[1][1] + pdf_fl*pdf_fl
   corrSystMatrix[1][2] = corrSystMatrix[1][2] + pdf_fl*pdf_fr
   
   corrSystMatrix[2][0] = corrSystMatrix[2][0] + pdf_fr*pdf_f0
   corrSystMatrix[2][1] = corrSystMatrix[2][1] + pdf_fr*pdf_fl
   corrSystMatrix[2][2] = corrSystMatrix[2][2] + pdf_fr*pdf_fr


   #  ** c. Template Statistics Unc
   ts_f0 = 0.005 # 8-channel
   ts_fl = 0.003 # 8-channel
   ts_fr = 0.003 # 8-channel

   if lep2incl:
      ts_f0 = 0.009 # lep2incl
      ts_fl = 0.006 # lep2incl
      ts_fr = 0.004 # lep2incl

   corrSystMatrix[0][0] = corrSystMatrix[0][0] + ts_f0*ts_f0
   corrSystMatrix[0][1] = corrSystMatrix[0][1] + ts_f0*ts_fl
   corrSystMatrix[0][2] = corrSystMatrix[0][2] + ts_f0*ts_fr
   
   corrSystMatrix[1][0] = corrSystMatrix[1][0] + ts_fl*ts_f0
   corrSystMatrix[1][1] = corrSystMatrix[1][1] + ts_fl*ts_fl
   corrSystMatrix[1][2] = corrSystMatrix[1][2] + ts_fl*ts_fr
   
   corrSystMatrix[2][0] = corrSystMatrix[2][0] + ts_fr*ts_f0
   corrSystMatrix[2][1] = corrSystMatrix[2][1] + ts_fr*ts_fl
   corrSystMatrix[2][2] = corrSystMatrix[2][2] + ts_fr*ts_fr

   # *** 4. Print Results
   print "Total Systematic CORR"
   printMatrix(corrSystMatrix)
   
   print "Total Statistical CORR"
   printMatrix(corrStatMatrix)
   
   # *** 5. Calculate and print total CORR matrix
   corrMatrix[0][0] = corrStatMatrix[0][0] + corrSystMatrix[0][0]
   corrMatrix[0][1] = corrStatMatrix[0][1] + corrSystMatrix[0][1]
   corrMatrix[0][2] = corrStatMatrix[0][2] + corrSystMatrix[0][2]
   corrMatrix[1][0] = corrStatMatrix[1][0] + corrSystMatrix[1][0]
   corrMatrix[1][1] = corrStatMatrix[1][1] + corrSystMatrix[1][1]
   corrMatrix[1][2] = corrStatMatrix[1][2] + corrSystMatrix[1][2]
   corrMatrix[2][0] = corrStatMatrix[2][0] + corrSystMatrix[2][0]
   corrMatrix[2][1] = corrStatMatrix[2][1] + corrSystMatrix[2][1]
   corrMatrix[2][2] = corrStatMatrix[2][2] + corrSystMatrix[2][2]
   
   print "Total CORR"
   printMatrix(corrMatrix)

   print " #### Correlation Coeffs ####"
   rho_0L = corrMatrix[0][1] / ( math.sqrt(corrMatrix[0][0]) * math.sqrt(corrMatrix[1][1]) )
   rho_0R = corrMatrix[0][2] / ( math.sqrt(corrMatrix[0][0]) * math.sqrt(corrMatrix[2][2]) )
   rho_LR = corrMatrix[1][2] / ( math.sqrt(corrMatrix[1][1]) * math.sqrt(corrMatrix[2][2]) )
   print "rho(F0, FL) = ", rho_0L
   print "rho(F0, FR) = ", rho_0R
   print "rho(FL, FR) = ", rho_LR

   #print 'chk0:',chk0, '\t chkL:', chkL,  '\t chkR:', chkR

   c1 = ROOT.TCanvas("c1","c1",800,800)
   c1.Divide(3,1)
   c1.cd(1)
   h_c0L.Draw()
   c1.cd(2)
   h_c0R.Draw()
   c1.cd(3)
   h_cLR.Draw()
   c1.Print("correlation1D_systs.png")
   
   c0.cd()
   h_d0L.SetXTitle("#Delta#rho(F_{0},F_{L})")
   h_d0L.SetYTitle("Entries / Bin")
   h_d0L.Draw()
   c0.Print("diffCorrelationFactor_cf0fL.png")

   cL.cd()
   h_d0R.SetXTitle("#Delta#rho(F_{0},F_{R})")
   h_d0R.SetYTitle("Entries / Bin")
   h_d0R.Draw()
   cL.Print("diffCorrelationFactor_cf0fR.png")

   cR.cd()
   h_dLR.SetXTitle("#Delta#rho(yF_{L},F_{R})")
   h_dLR.SetYTitle("Entries / Bin")
   h_dLR.Draw()
   cR.Print("diffCorrelationFactor_cfLfR.png")

   return
         
def printMatrix (matrix):
   """print matrix. self-explanatory"""
            
   tot = matrix[0][0]

   print round(matrix[0][0],5),'\t',round(matrix[0][1],5),'\t',round(matrix[0][2],5)
   print round(matrix[1][0],5),'\t',round(matrix[1][1],5),'\t',round(matrix[1][2],5)
   print round(matrix[2][0],5),'\t',round(matrix[2][1],5),'\t',round(matrix[2][2],5)
   
def callUncertainties(region, i, whichFraction):
   """set hard-coded uncertainties for top mass, top pdf, template stat, and stat+bkg unc"""

   fraction = -1
   if whichFraction == 'F0':
      fraction = 0
   if whichFraction == 'FL':
      fraction = 1
   if whichFraction == 'FR':
      fraction = 2

   #statBkgUnc = [0, 0, 0]
   #PDFunc = [0, 0, 0]
   #TopMassUnc = [0, 0, 0]
   #TemplateStatUnc = [0, 0, 0]
 
   #lep 2incl
   statBkg_lep2incl         = [0.012, 0.008, 0.006]
   PDF_lep2incl             = [0.003, 0.004, 0.001]
   TopMass_lep2incl         = [0.002, 0.005, 0.003]
   TemplateStat_lep2incl    = [0.009, 0.006, 0.004]
   #had 2incl
   statBkg_had2incl         = [0.012, 0.022, 0.023]
   PDF_had2incl             = [0.002, 0.003, 0.002]
   TopMass_had2incl         = [0.001, 0.008, 0.007]
   TemplateStat_had2incl    = [0.009, 0.017, 0.017]
   #lep+had 2incl
   statBkg_lephad2incl      = [0.008, 0.006, 0.004]
   PDF_lephad2incl          = [0.003, 0.004, 0.001]
   TopMass_lephad2incl      = [0.001, 0.003, 0.004]
   TemplateStat_lephad2incl = [0.006, 0.004, 0.004]
   #lep 1excl+2incl
   statBkg_lepBTag          = [0.010, 0.006, 0.005]
   PDF_lepBTag              = [0.003, 0.004, 0.001]
   TopMass_lepBTag          = [0.002, 0.005, 0.003]
   TemplateStat_lepBTag     = [0.008, 0.005, 0.004] 
   #had 1excl+2incl
   statBkg_hadBTag          = [0.010, 0.021, 0.022]
   PDF_hadBTag              = [0.001, 0.002, 0.002]
   TopMass_hadBTag          = [0.003, 0.0095, 0.0068]
   TemplateStat_hadBTag     = [0.008, 0.016, 0.016]

   #lep+had 1excl+2incl
   statBkg_lephadBTag       = [0.007, 0.005, 0.004]
   PDF_lephadBTag           = [0.002, 0.003, 0.002]
   TopMass_lephadBTag       = [0.001, 0.003, 0.004]
   TemplateStat_lephadBTag  = [0.005, 0.004, 0.003]

   if region=="lep2incl":
      statBkgUnc[i]      = statBkg_lep2incl[fraction]
      PDFunc[i]          = PDF_lep2incl[fraction]
      TopMassUnc[i]      = TopMass_lep2incl[fraction]
      TemplateStatUnc[i] = TemplateStat_lep2incl[fraction] 
   if region=="had2incl":
      statBkgUnc[i]      = statBkg_had2incl[fraction]
      PDFunc[i]          = PDF_had2incl[fraction]
      TopMassUnc[i]      = TopMass_had2incl[fraction]
      TemplateStatUnc[i] = TemplateStat_had2incl[fraction] 
   if region=="lephad2incl":
      statBkgUnc[i]      = statBkg_lephad2incl[fraction]
      PDFunc[i]          = PDF_lephad2incl[fraction]
      TopMassUnc[i]      = TopMass_lephad2incl[fraction]
      TemplateStatUnc[i] = TemplateStat_lephad2incl[fraction] 
   if region=="lepBTag":
      statBkgUnc[i]      = statBkg_lepBTag[fraction]
      PDFunc[i]          = PDF_lepBTag[fraction]
      TopMassUnc[i]      = TopMass_lepBTag[fraction]
      TemplateStatUnc[i] = TemplateStat_lepBTag[fraction] 
   if region=="hadBTag":
      statBkgUnc[i]      = statBkg_hadBTag[fraction]
      PDFunc[i]          = PDF_hadBTag[fraction]
      TopMassUnc[i]      = TopMass_hadBTag[fraction]
      TemplateStatUnc[i] = TemplateStat_hadBTag[fraction] 
   if region=="lephadBTag":
      statBkgUnc[i]      = statBkg_lephadBTag[fraction]
      PDFunc[i]          = PDF_lephadBTag[fraction]
      TopMassUnc[i]      = TopMass_lephadBTag[fraction]
      TemplateStatUnc[i] = TemplateStat_lephadBTag[fraction] 

   #print "$$$$$$$$$$$$$", statBkgUnc[0], statBkgUnc[1], statBkgUnc[2], whichFraction

def setUncertainties(reg0, reg1, reg2, fraction):
   """ simplify call to setUncertainties --> only need one line instead of three in main"""
   callUncertainties(reg0, 0, fraction) 
   callUncertainties(reg1, 1, fraction) 
   callUncertainties(reg2, 2, fraction) 

def setUncertainties4(reg0, reg1, reg2, reg3, fraction):
   """ simplify call to setUncertainties --> only need one line instead of three in main"""
   callUncertainties(reg0, 0, fraction) 
   callUncertainties(reg1, 1, fraction) 
   callUncertainties(reg2, 2, fraction) 
   callUncertainties(reg3, 3, fraction) 


def makePaperTable( lepDict, hadDict ):
   r0='lep2incl'
   r1='had2incl'
   r2='hadBTag'

   setUncertainties(r0, r1, r2, 'F0') 
   lepTotals_f0 = calcTotals (lepDict, 'F0', statBkgUnc[0], PDFunc[0], TopMassUnc[0], TemplateStatUnc[0])
   hadTotals_f0 = calcTotals (hadDict, 'F0', statBkgUnc[1], PDFunc[1], TopMassUnc[1], TemplateStatUnc[1])
   setUncertainties(r0, r1, r2, 'FL') 
   lepTotals_fl = calcTotals (lepDict, 'FL', statBkgUnc[0], PDFunc[0], TopMassUnc[0], TemplateStatUnc[0])
   hadTotals_fl = calcTotals (hadDict, 'FL', statBkgUnc[1], PDFunc[1], TopMassUnc[1], TemplateStatUnc[1])
   setUncertainties(r0, r1, r2, 'FR') 
   lepTotals_fr = calcTotals (lepDict, 'FR', statBkgUnc[0], PDFunc[0], TopMassUnc[0], TemplateStatUnc[0])
   hadTotals_fr = calcTotals (hadDict, 'FR', statBkgUnc[1], PDFunc[1], TopMassUnc[1], TemplateStatUnc[1])


   print "\n### el+mu 2incl {0} ###".format(whichFraction)
   print "\\begin{table}[h!]"
   print "\centering"
   print "\\begin{tabular}{lcccc}"
   print "\hline\hline"
   if whichFraction == 'F0':
      print "\multicolumn{5}{c}{\\fo}\\\\"
   elif whichFraction == 'FL':
      print "\multicolumn{5}{c}{\\fl}\\\\"
   elif whichFraction == 'FR':
      print "\multicolumn{5}{c}{\\fr}\\\\"
   print "\hline"
   print "Systematic uncertainty & N$_{syst}$ & Lep 2incl & Lep+Had 2incl & Lep+Had 1excl+2incl \\\\\\hline"
   print "\multicolumn{5}{c}{Reconstructed Objects} \\\\\\hline"
   if method == 'alt':
      
      printThreePlusMinusLine('Muon', 6, lephadTotals[0], lepTotals[1], lepTotals[2], hadTotals[1], hadTotals[2], lephadTotals[1], lephadTotals[2])
      printThreePlusMinusLine('Electron', 5, lephadTotals[3], lepTotals[4], lepTotals[5], hadTotals[4], hadTotals[5], lephadTotals[4], lephadTotals[5])

### MAIN CODE ###

# *** -1. Graveyard of un-used functions
#makeTable( 'F0', el2incl, mu2incl, elmu2incl, elmuBTag)
#calcTotalUnc( 'F0', el2incl, mu2incl, elmu2incl, elmuBTag)
#getTopTen ( 'F0', el2incl, mu2incl, elmu2incl, elmuBTag)
#makeFinalTable( 'F0', elmu2incl_lep, elmu2incl)
#makeFinalSensitivityTable('F0', elmuBTag)
#makeFinalSensitivityTable('F0', lephadBTag, 'alt')


# *** 0. Loop over files and build dictionaries
loopFile( d_lephadBTag, f_lephadBTag)
loopFile( d_lepBTag, f_lepBTag)
loopFile( d_hadBTag, f_hadBTag)
loopFile( d_lep2incl, f_lep2incl)
loopFile( d_had2incl, f_had2incl)
loopFile( d_lephad2incl, f_lephad2incl)

# *** 1. Set Template Stat, Top Mass, PDF, and stat+bkg Unc

#setUncertainties('adsf', 0)

# *** 2. Dump table with full up/down information for all systematic variations
r0='lep2incl'
r1='hadBTag'
#r2='had2incl'
#r3='hadBTag'
#r0='lepBTag'
#r2='hadBTag'
r2='lephadBTag'


setUncertainties(r0, r1, r2, 'F0') 
#makeFullTable( 'F0', d_lepBTag, d_hadBTag, d_lephadBTag)
makeFullTable( 'F0', d_lep2incl, d_hadBTag, d_lephadBTag)
#makeFullTable( 'F0', d_lep2incl, d_had2incl, d_hadBTag)
#makeFullTable( 'F0', d_lep2incl, d_had2incl, d_lepBTag)
setUncertainties(r0, r1, r2, 'FL') 
#makeFullTable( 'FL', d_lep2incl, d_had2incl, d_lepBTag)
makeFullTable( 'FL', d_lep2incl, d_hadBTag, d_lephadBTag)
#makeFullTable( 'FL', d_lep2incl, d_had2incl, d_hadBTag)
setUncertainties(r0, r1, r2, 'FR') 
#makeFullTable( 'FR', d_lep2incl, d_had2incl, d_hadBTag)
#makeFullTable( 'FR', d_lep2incl, d_had2incl, d_lepBTag)
makeFullTable( 'FR', d_lep2incl, d_hadBTag, d_lephadBTag)

# *** 3. Dump table all kept systematics summarized by category


setUncertainties(r0, r1, r2, 'F0') 
makeFullSensitivityTable('F0', d_lep2incl, d_hadBTag, d_lephadBTag, 'alt')
#setUncertainties(r0, r1, r2, 'F0') 
#makeFullSensitivityTable('F0', d_lep2incl, d_hadBTag, d_lephadBTag, 'alt')
#makeFullSensitivityTable('F0', d_lep2incl, d_had2incl, d_hadBTag, 'alt')
setUncertainties(r0, r1, r2, 'FL') 
makeFullSensitivityTable('FL', d_lep2incl, d_hadBTag, d_lephadBTag, 'alt')
#makeFullSensitivityTable('FL', d_lep2incl, d_had2incl, d_lephadBTag, 'alt')
setUncertainties(r0, r1, r2, 'FR') 
makeFullSensitivityTable('FR', d_lep2incl, d_hadBTag, d_hadBTag, 'alt')

#setUncertainties4(r0, r1, r2, r3, 'F0') 
#makeFullSensitivityTable4('F0', d_lep2incl, d_lepBTag, d_had2incl, d_hadBTag, 'alt')
#setUncertainties4(r0, r1, r2, r3, 'FL') 
#makeFullSensitivityTable4('FL', d_lep2incl, d_lepBTag, d_had2incl, d_hadBTag, 'alt')
#setUncertainties4(r0, r1, r2, r3, 'FR') 
#makeFullSensitivityTable4('FR', d_lep2incl, d_lepBTag, d_had2incl, d_hadBTag, 'alt')

# *** 4. dump table with full up/down information for lephad --> should show sum is zero (within pseudoExperiment errors)
#makeSumTable( d_lephadBTag)

# *** 5. Calculate full covariance matrix
#makeCorrelationMatrix( d_lephadBTag)
#makeCorrelationMatrix( d_lep2incl)
#makeCorrelationMatrix( d_lepBTag)
makeCorrelationMatrixTXT( '/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/syst_outputs_lep_2incl_aug02/ExternalSystematicsOutput_el_mu_2incl_3D_3W_lep/SystematicOutput_el_mu_lep.script.txt')

# *** 6. Make paper table
#makePaperTable(d_lep2incl, d_had2incl)

sys.exit()



#  LocalWords:  fL
