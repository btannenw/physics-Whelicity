

#       WriteParameterStatus("EvaluateExternal", "Input Channel: "               + string(argv[1]));
#       WriteParameterStatus("EvaluateExternal", "Number of PE: "                + string(argv[2]));
#       WriteParameterStatus("EvaluateExternal", "Systematic in eval: "          + string(argv[3]));
#       WriteParameterStatus("EvaluateExternal", "Fitting mode: "                + string(argv[4]));

# ------------- How to use ------------------
# InputChannel: either el or mu, or el_mu_BTag
# Number of PE: Pseudoexperiments: 0-5000
# Systematic: Datafit , Calibration , TemplateStat , jes , muid , Nominal,
#			  mums, ees, eer, cellout, softjet, pileup, btag, ctag, mistag, muonid, muonreco, muontrig, electronid, jvfsf
#			  electrontrig, electronreco, jer, jeff, musc, ISR_FSR, UE, TopMass, CR_Perugia, MCGen1, MCGen2, 
# Fitting mode: 2D or 3D: number of signal templates
# NewSample: true or false
# CalibrationMode: dual or single (single: use same signal for PseudoData and systematicPD)
# 1excl, 1incl, 2incl, 1excl2incl 
#./runEvaluateExternalUpdate mu 1000 Datafit 3D #BkgFit

#./runEvaluateExternalUpdate el_mu_BTag 1000 Calibration 3D true single 1excl2incl #BkgFit
./runEvaluateExternalUpdate el 1000 Calibration 3D true single 2incl #BkgFit

#./runEvaluateExternalUpdate mu 5000 Calibration 3D true single 1incl #BkgFit 
