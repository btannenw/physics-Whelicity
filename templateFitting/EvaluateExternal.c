#include "StatusLogbook.h"
#include "ProfilingClass.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
 
     if(argc != 4){

       WriteErrorStatus("EvaluateExternal", "Number of input variables is wrong!!!");
       WriteErrorStatus("EvaluateExternal", "EXIT   ");

       return 1;

     }
     else{
       
       WriteParameterStatus("EvaluateExternal", "Input Channel: "               + string(argv[1]));
       WriteParameterStatus("EvaluateExternal", "Number of PE: "                + string(argv[2]));
       WriteParameterStatus("EvaluateExternal", "Systematic in eval: "          + string(argv[3]));

     }
     
     std::string InputChannel   = argv[1];
     int NumberOfPE             = TString(argv[2]).Atoi();
     std::string Systematic     = argv[3];

     int NumberOfBins  =   15;
     double lower_edge = -1.0;
     double upper_edge =  1.0;

     //Globally defined histograms, needed for fit function
     hist      = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
     hist_sum  = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
     hist_help = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);

     ProfilingClass fProf = ProfilingClass();

     // add signal parameter with param name and hist name
     fProf.AddSignalParameter("N0",  "CosTheta_F0", "F0eff");
     fProf.AddSignalParameter("NL",  "CosTheta_FL", "FLeff");
     fProf.AddSignalParameter("NR",  "CosTheta_FR", "FReff");

     // add bkg parameter with param name and hist name
     fProf.AddBackgroundParameter("Wjets",  "Wjets",  "WjetsNorm",  "WjetsUnc");
     fProf.AddBackgroundParameter("QCD",    "QCD",    "QCDNorm",    "QCDUnc");
     fProf.AddBackgroundParameter("RemBkg", "RemBkg", "RemBkgNorm", "RemBkgUnc");

     // only need to specify data hist name
     //     fProf.SetDataHist("Data");
     
     std::string fInputFolder  = "root_files_05122012_4incl";
     std::string fOutputFolder = "ExternalSystematicsOutput";

     std::string OutputTxtFile = fInputFolder+"/SystematicOutput.txt";

     fProf.SetOutputFolder(fOutputFolder);
     fProf.SetOutputTxtFile(OutputTxtFile);

     fProf.SetNominalInputFile(fInputFolder+"/Templates_" + InputChannel + ".root");

     if(Systematic == "jes"){

       fProf.AddSystematicPseudoData("Nominal", fInputFolder+"/Templates_" + InputChannel + ".root",         "PseudoData");
       fProf.AddSystematicPseudoData("JESup",   fInputFolder+"/Templates_" +InputChannel+ "_jesu_1k.root", "PseudoData");
       fProf.AddSystematicPseudoData("JESdown", fInputFolder+"/Templates_" +InputChannel+ "_jesd_1k.root", "PseudoData");
      
       fProf.SetEvaluationMode("LargestDiff");
       
     }
     else if(Systematic == "muid"){

       fProf.AddSystematicPseudoData("Nominal",  fInputFolder+"/Templates_Combined_" + InputChannel + ".root",          "PseudoData");
       fProf.AddSystematicPseudoData("MUIDup",   fInputFolder+"/Templates_muidu_4inclJets_1tagInTags_MCatNLO_mu.root", "PseudoData");
       fProf.AddSystematicPseudoData("MUIDdown", fInputFolder+"/Templates_muidd_4inclJets_1tagInTags_MCatNLO_mu.root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "mums"){
       
       fProf.AddSystematicPseudoData("Nominal",  fInputFolder+"/Templates_Combined_" + InputChannel + ".root",          "PseudoData");
       fProf.AddSystematicPseudoData("MUMSup",   fInputFolder+"/Templates_"+Systematic+"u_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("MUMSdown", fInputFolder+"/Templates_"+Systematic+"d_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "ees"){

       fProf.AddSystematicPseudoData("Nominal",  fInputFolder+"/Templates_Combined_" + InputChannel + ".root",          "PseudoData");
       fProf.AddSystematicPseudoData("EESup",   fInputFolder+"/Templates_"+Systematic+"u_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("EESdown", fInputFolder+"/Templates_"+Systematic+"d_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "eer"){

       fProf.AddSystematicPseudoData("Nominal",  fInputFolder+"/Templates_Combined_" + InputChannel + ".root",          "PseudoData");
       fProf.AddSystematicPseudoData("EERup",   fInputFolder+"/Templates_"+Systematic+"u_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("EERdown", fInputFolder+"/Templates_"+Systematic+"d_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "cellout"){

       fProf.AddSystematicPseudoData("Nominal",     fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("CellOutup",   fInputFolder+"/Templates_"+Systematic+"u_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("CellOutdown", fInputFolder+"/Templates_"+Systematic+"d_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "softjet"){

       fProf.AddSystematicPseudoData("Nominal",     fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("Softjetup",   fInputFolder+"/Templates_"+Systematic+"u_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("Softjetdown", fInputFolder+"/Templates_"+Systematic+"d_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "pileup"){

       fProf.AddSystematicPseudoData("Nominal",     fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("Pileupup",    fInputFolder+"/Templates_"+Systematic+"u_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("Pileupdown",  fInputFolder+"/Templates_"+Systematic+"d_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "btag"){

       fProf.AddSystematicPseudoData("Nominal",  fInputFolder+"/Templates_Combined_" + InputChannel + ".root",          "PseudoData");
       fProf.AddSystematicPseudoData("BTAGup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("BTAGdown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "ctag"){

       fProf.AddSystematicPseudoData("Nominal",  fInputFolder+"/Templates_Combined_" + InputChannel + ".root",         "PseudoData");
       fProf.AddSystematicPseudoData("CTAGup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("CTAGdown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "mistag"){

       fProf.AddSystematicPseudoData("Nominal",    fInputFolder+"/Templates_Combined_" + InputChannel + ".root",            "PseudoData");
       fProf.AddSystematicPseudoData("MISTAGup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("MISTAGdown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");


       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "muonid"){

       fProf.AddSystematicPseudoData("Nominal",    fInputFolder+"/Templates_Combined_" + InputChannel + ".root",            "PseudoData");
       fProf.AddSystematicPseudoData("MuonIdup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("MuonIddown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "muonreco"){

       fProf.AddSystematicPseudoData("Nominal",      fInputFolder+"/Templates_Combined_" + InputChannel + ".root",              "PseudoData");
       fProf.AddSystematicPseudoData("MuonRecoup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("MuonRecodown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       
       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "muontrig"){

       fProf.AddSystematicPseudoData("Nominal",      fInputFolder+"/Templates_Combined_" + InputChannel + ".root",              "PseudoData");
       fProf.AddSystematicPseudoData("MuonTrigup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("MuonTrigdown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "electronid"){

       fProf.AddSystematicPseudoData("Nominal",        fInputFolder+"/Templates_Combined_" + InputChannel + ".root",            "PseudoData");
       fProf.AddSystematicPseudoData("ElectronIdup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("ElectronIddown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "jvfsf"){

       fProf.AddSystematicPseudoData("Nominal",     fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("JVFSFup",     fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("JVFSFdown",   fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "electrontrig"){

       fProf.AddSystematicPseudoData("Nominal",      fInputFolder+"/Templates_Combined_" + InputChannel + ".root",              "PseudoData");
       fProf.AddSystematicPseudoData("ElectronTrigup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("ElectronTrigdown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "electronreco"){

       fProf.AddSystematicPseudoData("Nominal",      fInputFolder+"/Templates_Combined_" + InputChannel + ".root",              "PseudoData");
       fProf.AddSystematicPseudoData("ElectronRecoup",   fInputFolder+"/Templates_"+Systematic+"up_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
       fProf.AddSystematicPseudoData("ElectronRecodown", fInputFolder+"/Templates_"+Systematic+"down_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");

       fProf.SetEvaluationMode("LargestDiff");

     }
     else if(Systematic == "jer"){

       fProf.AddSystematicPseudoData("Nominal",    fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("BTAGup",   fInputFolder+"/Templates_"+Systematic+"_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root", "PseudoData");
  
       fProf.SetEvaluationMode("FullDiff");

     }
     else if(Systematic == "jeff"){

       fProf.AddSystematicPseudoData("Nominal",    fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("JEFF",       fInputFolder+"/Templates_"+Systematic+"_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root",     "PseudoData");

       fProf.SetEvaluationMode("FullDiff");

     }
     else if(Systematic == "musc"){

       fProf.AddSystematicPseudoData("Nominal",    fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("MUSC",       fInputFolder+"/Templates_"+Systematic+"_4inclJets_1tagInTags_MCatNLO_"+InputChannel+".root",     "PseudoData");

       fProf.SetEvaluationMode("FullDiff");

     }
     else if(Systematic == "ISR_FSR"){

       fProf.AddSystematicPseudoData("MorePS",     fInputFolder+"/Templates_MorePS_4inclJets_1tagInTags_MorePS_" + InputChannel + ".root",    "PseudoData");
       fProf.AddSystematicPseudoData("LessPS",     fInputFolder+"/Templates_LessPS_4inclJets_1tagInTags_LessPS_" + InputChannel + ".root",    "PseudoData");

       fProf.SetEvaluationMode("HalfDiff");

     }
     else if(Systematic == "mu_UE"){

       fProf.AddSystematicPseudoData("MoreUE",     fInputFolder+"/Templates_TTbar_MoreUE_SystType_" + InputChannel + ".root",    "PseudoData");
       fProf.AddSystematicPseudoData("LessUE",     fInputFolder+"/Templates_TTbar_LessUE_SystType_" + InputChannel + ".root",    "PseudoData");

       fProf.SetEvaluationMode("FullDiff");

     }
     else if(Systematic == "mu_TopMass"){

       fProf.AddSystematicPseudoData("Mtop167", fInputFolder+"/Templates_TTbar_TopMass167_SystType_" + InputChannel + ".root", "PseudoData");
       fProf.AddSystematicPseudoData("Mtop170", fInputFolder+"/Templates_TTbar_TopMass170_SystType_" + InputChannel + ".root", "PseudoData");
       fProf.AddSystematicPseudoData("Mtop175", fInputFolder+"/Templates_TTbar_TopMass175_SystType_" + InputChannel + ".root", "PseudoData");
       fProf.AddSystematicPseudoData("Mtop177", fInputFolder+"/Templates_TTbar_TopMass177_SystType_" + InputChannel + ".root", "PseudoData");

       fProf.SetEvaluationMode("Dependence");

     }
     else if(Systematic == "CR_Perugia"){
       
       fProf.AddSystematicPseudoData("Perugia",      fInputFolder+"/Templates_CR_Perugia_4inclJets_1tagInTags_CR_Perugia_" + InputChannel + ".root",     "PseudoData");
       fProf.AddSystematicPseudoData("Perugia_NoCR", fInputFolder+"/Templates_CR_Perugia_NoCR_4inclJets_1tagInTags_CR_Perugia_NoCR_" + InputChannel + ".root", "PseudoData");

       fProf.SetEvaluationMode("FullDiff");

     }
     else if(Systematic == "mu_CR_TuneAPro"){

       fProf.AddSystematicPseudoData("TuneACRPro",  fInputFolder+"/Templates_TTbar_CR_TuneACRPro_SystType_" + InputChannel + ".root", "PseudoData");
       fProf.AddSystematicPseudoData("TuneAPro",    fInputFolder+"/Templates_TTbar_CR_TuneAPro_SystType_" + InputChannel + ".root",    "PseudoData");

       fProf.SetEvaluationMode("FullDiff");

     }
     else if(Systematic == "MCGen1"){

       fProf.AddSystematicPseudoData("MCatNLO_AF",  fInputFolder+"/Templates_MCGen_MCatNLO_AF_4inclJets_1tagInTags_MCGen_MCatNLO_AF_" + InputChannel + ".root",    "PseudoData");
       fProf.AddSystematicPseudoData("Powheg_AF",   fInputFolder+"/Templates_MCGen_PowhegJimmy_AF_4inclJets_1tagInTags_MCGen_PowhegJimmy_AF_" + InputChannel + ".root", "PseudoData");

       fProf.SetEvaluationMode("FullDiff");
       
     }
     else if(Systematic == "MCGen2"){

       fProf.AddSystematicPseudoData("MCatNLO_FS",  fInputFolder+"/Templates_MCGen_MCatNLO_FS_4inclJets_1tagInTags_MCGen_MCatNLO_FS_" + InputChannel + ".root",    "PseudoData");
       fProf.AddSystematicPseudoData("AlpgenFS",    fInputFolder+"/Templates_MCGen_Alpgen_FS_4inclJets_1tagInTags_MCGen_Alpgen_FS_" + InputChannel + ".root", "PseudoData");

       fProf.SetEvaluationMode("FullDiff");

     }

     else if(Systematic == "Shower"){


       fProf.AddSystematicPseudoData("Powheg_Py_P2011C", fInputFolder+"/Templates_PS_PowhegPythia_P2011C_AF_4inclJets_1tagInTags_PS_PowhegPythia_P2011C_AF_" + InputChannel + ".root",    "PseudoData");
       fProf.AddSystematicPseudoData("Powheg_Jimmy",     fInputFolder+"/Templates_MCGen_PowhegJimmy_AF_4inclJets_1tagInTags_MCGen_PowhegJimmy_AF_" + InputChannel + ".root",         "PseudoData");

       fProf.SetEvaluationMode("FullDiff");

     }
     else if(Systematic == "mu_QCDShape"){
       
       fProf.AddSystematicPseudoData("Nominal", fInputFolder+"/Templates_Combined_" + InputChannel + ".root",             "PseudoData");
       fProf.AddSystematicPseudoData("QCD_MMA", fInputFolder+"/Templates_TTbar_QCDMMA_SystType_" + InputChannel + ".root",   "PseudoData");
       fProf.AddSystematicPseudoData("QCD_MMB", fInputFolder+"/Templates_TTbar_QCDMMB_SystType_" + InputChannel + ".root",   "PseudoData");

     fProf.SetEvaluationMode("LargestDiff");
     
     }

     // lumi in pb
     fProf.SetLumi(4655.74);

     // xsec in pb
     fProf.SetXsec(89.2);

     fProf.SetInputChannel(InputChannel);
     fProf.SetNumberOfPE(NumberOfPE);
     fProf.SetNumberOfBins(NumberOfBins);
     fProf.SetLowerEdge(lower_edge);
     fProf.SetUpperEdge(upper_edge);
    
     fProf.ReadNominalInfo();

     WriteInfoStatus("EvaluateExternal", "Do Validation");

     fProf.DoExternalSystematicEvaluation(Systematic);

     WriteInfoStatus("EvaluateExternal", "finalize");

     WriteInfoStatus("EvaluateExternal", "finished");

    return 0;

}
