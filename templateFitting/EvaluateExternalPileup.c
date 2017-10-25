#include "StatusLogbook.h"
#include "ProfilingClass.h"
#include "plots.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
 
     if(argc != 6){
    
       WriteErrorStatus("EvaluateExternal", "Number of input variables is wrong!!!");
       WriteErrorStatus("EvaluateExternal", "EXIT   ");
       
       return 1;
       
     } 
     else{
    
       WriteParameterStatus("EvaluateExternal", "Input Channel: "               + string(argv[1]));
       WriteParameterStatus("EvaluateExternal", "Number of PE: "                + string(argv[2]));
       WriteParameterStatus("EvaluateExternal", "Systematic in eval: "          + string(argv[3]));
       WriteParameterStatus("EvaluateExternal", "Fitting mode: "                + string(argv[4]));

     }
     
     std::string InputChannel   = argv[1];
     int NumberOfPE             = TString(argv[2]).Atoi();
     std::string Systematic     = argv[3];
     std::string FittingMode    = argv[4];

     int i = TString(argv[5]).Atoi();
     
     if(Systematic == "PrimaryVertex" || Systematic == "Mu"){

	 int NumberOfBins  =   15;
	 double lower_edge = -1.0;
	 double upper_edge =  1.0;
	 
	 std::string fInputFolder  = "TemplateFiles_Pileup_new";
	 
	 std::string TemplateFile;
	 
	 std::stringstream oss;
	 oss << i;
	
	 if(Systematic == "PrimaryVertex")
	   TemplateFile = fInputFolder+"/TotalTemplatesFile_4inclJets_1inclTags_" + InputChannel + "_CosTheta_Primary_"+oss.str()+".root";
	 else
	   TemplateFile = fInputFolder+"/TotalTemplatesFile_4inclJets_1inclTags_" + InputChannel + "_CosTheta_Mu_"+oss.str()+".root";

	 std::cout << TemplateFile.c_str() << std::endl;
	 
	 if(InputChannel == "el_mu"){
	   
	   NumberOfBins =   30;
	   lower_edge   = -2.0;
	   upper_edge   =  2.0;
	   
	 }
	 
	 std::vector<std::string> VarNames;
	 std::vector<std::string> HistoNames;
	 
	 //Globally defined histograms, needed for fit function
	 hist      = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
	 hist_sum  = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
	 hist_help = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
	 
	 ProfilingClass fProf = ProfilingClass();
	 
	 if(InputChannel == "el" || InputChannel == "mu"){
	   
	   if(InputChannel == "el")
	     fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "#geq 1 b-tags");
	   else if(InputChannel == "mu")
	     fProf.SetSubplotLabels("#mu+jets",   "#geq 4 jets", "#geq 1 b-tags");
	   
	   // add signal parameter with param name and hist name
	   
	   fProf.AddSignalParameter("N0",  "CosTheta_F0", "F0eff");
	   fProf.AddSignalParameter("NL",  "CosTheta_FL", "FLeff");
	   if(FittingMode == "3D")
	     fProf.AddSignalParameter("NR",  "CosTheta_FR", "FReff");
	   
	   // add bkg parameter with param name and hist name
	   fProf.AddBackgroundParameter("Wjets",  "Wjets",  "WjetsNorm",  "WjetsUnc");
	   fProf.AddBackgroundParameter("QCD",    "QCD",    "QCDNorm",    "QCDUnc");
	   fProf.AddBackgroundParameter("RemBkg", "RemBkg", "RemBkgNorm", "RemBkgUnc");
	   
	   HistoNames.push_back("CosTheta_F0");
	   HistoNames.push_back("CosTheta_FL");
	   if(FittingMode == "3D")
	     HistoNames.push_back("CosTheta_FR");
	   
	   HistoNames.push_back("Wjets");
	   HistoNames.push_back("QCD");
	   HistoNames.push_back("RemBkg");
	   
	   VarNames.push_back("F0eff");
	   VarNames.push_back("FLeff");
	   if(FittingMode == "3D")
	     VarNames.push_back("FReff");
	   VarNames.push_back("WjetsNorm");
	   VarNames.push_back("QCDNorm");
	   VarNames.push_back("RemBkgNorm");
	   
	 }
	 else if (InputChannel == "el_mu"){
	   
	   NumberOfBins =   30;
	   lower_edge   = -2.0;
	   upper_edge   =  2.0;
	   
	   fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "#geq 1 b-tags");
	   fProf.SetSubplotLabels("#mu+jets", "#geq 4 jets", "#geq 1 b-tags");
	   
	   // add signal parameter with param name and hist name
	   fProf.AddSignalParameter("N0",  "CosTheta_F0_allCombined", "F0eff_allCombined");
	   fProf.AddSignalParameter("NL",  "CosTheta_FL_allCombined", "FLeff_allCombined");
	   if(FittingMode == "3D")
	     fProf.AddSignalParameter("NR",  "CosTheta_FR_allCombined", "FReff_allCombined");
	   
	   // add bkg parameter with param name and hist name
	   fProf.AddBackgroundParameter("Wjets_4incl_1incl_el",  "Wjets_4incl_1incl_el",  "WjetsNorm_4incl_1incl_el",  "WjetsUnc_4incl_1incl_el");
	   fProf.AddBackgroundParameter("Wjets_4incl_1incl_mu",  "Wjets_4incl_1incl_mu",  "WjetsNorm_4incl_1incl_mu",  "WjetsUnc_4incl_1incl_mu");
	   fProf.AddBackgroundParameter("QCD_4incl_1incl_el",    "QCD_4incl_1incl_el",    "QCDNorm_4incl_1incl_el",    "QCDUnc_4incl_1incl_el");
	   fProf.AddBackgroundParameter("QCD_4incl_1incl_mu",    "QCD_4incl_1incl_mu",    "QCDNorm_4incl_1incl_mu",    "QCDUnc_4incl_1incl_mu");
	   fProf.AddBackgroundParameter("RemBkg_4incl_1incl_el", "RemBkg_4incl_1incl_el", "RemBkgNorm_4incl_1incl_el", "RemBkgUnc_4incl_1incl_el");
	   fProf.AddBackgroundParameter("RemBkg_4incl_1incl_mu", "RemBkg_4incl_1incl_mu", "RemBkgNorm_4incl_1incl_mu", "RemBkgUnc_4incl_1incl_mu");
	   
	   HistoNames.push_back("CosTheta_F0_allCombined");
	   HistoNames.push_back("CosTheta_FL_allCombined");
	   if(FittingMode == "3D")
	     HistoNames.push_back("CosTheta_FR_allCombined");
	   HistoNames.push_back("Wjets_4incl_1incl_el");
	   HistoNames.push_back("Wjets_4incl_1incl_mu");
	   HistoNames.push_back("QCD_4incl_1incl_el");
	   HistoNames.push_back("QCD_4incl_1incl_mu");
	   HistoNames.push_back("RemBkg_4incl_1incl_el");
	   HistoNames.push_back("RemBkg_4incl_1incl_mu");
	   
	   VarNames.push_back("F0eff_allCombined");
	   VarNames.push_back("FLeff_allCombined");
	   if(FittingMode == "3D")
	     VarNames.push_back("FReff_allCombined");
	   VarNames.push_back("WjetsNorm_4incl_1incl_el");
	   VarNames.push_back("WjetsNorm_4incl_1incl_mu");
	   VarNames.push_back("QCDNorm_4incl_1incl_el");
	   VarNames.push_back("QCDNorm_4incl_1incl_mu");
	   VarNames.push_back("RemBkgNorm_4incl_1incl_el");
	   VarNames.push_back("RemBkgNorm_4incl_1incl_mu");
	   
	 }
	 
	 std::string fOutputFolder;
	 std::string fOutputFolderCali;

	 if(Systematic == "PrimaryVertex")
	   fOutputFolder     = "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PrimaryVertex_"+oss.str();
	 else
	   fOutputFolder     = "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/Mu_"+oss.str();

	 std::string OutputTxtFile = fOutputFolder+"/SystematicOutput_"+InputChannel+".txt";
	 
	 fProf.SetOutputFolder(fOutputFolder);
	 
	 fProf.SetOutputTxtFile(OutputTxtFile);
	 
	 fProf.SetNominalInputFile(TemplateFile);
	 
	 fProf.SetTemplateFittingMode(FittingMode);
	 
	 // lumi in pb
	 fProf.SetLumi(4655.74);
	 
	 // xsec in pb
	 fProf.SetXsec(89.2);

	 fProf.SetDataHist("Data");
	 
	 fProf.AddSystematicPseudoData("Data", TemplateFile, "Data");
	 
	 fProf.SetFractions(0.698, 0.301, 0.00041);
	 
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
	 
	 fProf.DoExternalSystematicEvaluation("Datafit");
	 
	 WriteInfoStatus("EvaluateExternal", "finalize");
	 
	 WriteInfoStatus("EvaluateExternal", "finished");
	 
	 //       }
       
     }
     else if(Systematic == "PrimaryVertexEval"){
       
       TGraphErrors *fGraphF0 = new TGraphErrors(6);
       TGraphErrors *fGraphFL = new TGraphErrors(6);
       TGraphErrors *fGraphFR = new TGraphErrors(6);

       for(int k = 2; k < 8; k++){
       
	 std::stringstream oss;
	 oss << k;

	 TFile *fNewFile = new TFile(("ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PrimaryVertex_"+oss.str()+"/DatafitCorrectedUnc_"+InputChannel+".root").c_str(), "READ");
     
	 TTree *tree  = (TTree*) fNewFile -> Get("EnsembleTree");

	 double F0,     FL,     FR;
	 double F0err1, FLerr1, FRerr1;
	 double  F0err,  FLerr,  FRerr;

	 tree -> SetBranchAddress("F0",     &F0);
	 tree -> SetBranchAddress("FL",     &FL);
	 tree -> SetBranchAddress("FR",     &FR);

	 tree -> SetBranchAddress("F0_err", &F0err);
	 tree -> SetBranchAddress("FL_err", &FLerr);
	 tree -> SetBranchAddress("FR_err", &FRerr);
	 
	 tree -> GetEntry(0);

	 std::cout << F0    << "\t" << FL    << "\t" << FR    << std::endl;
	 std::cout << F0err << "\t" << FLerr << "\t" << FRerr << std::endl;
	 std::cout << "\t"  << std::endl;

	 fGraphF0 -> SetPoint(k-2, k, F0);
	 fGraphF0 -> SetPointError(k-2, 0, F0err);

	 fGraphFL -> SetPoint(k-2, k, FL);
         fGraphFL -> SetPointError(k-2, 0, FLerr);

	 fGraphFR -> SetPoint(k-2, k, FR);
         fGraphFR -> SetPointError(k-2, 0, FRerr);

	 fNewFile -> Close();


       }

       PlotPileupDistribution (fGraphF0, "No. primary vertices", "F_{0}", "BLA", 0.5, 9.5, 0.00, 1.30, InputChannel, "", "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PileupDependence_F0_"+InputChannel+"_"+FittingMode+".eps");       
       PlotPileupDistribution (fGraphFL, "No. primary vertices", "F_{L}", "BLA", 0.5, 9.5, -0.20, 0.80, InputChannel, "", "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PileupDependence_FL_"+InputChannel+"_"+FittingMode+".eps");
   
       if(FittingMode == "3D")
	 PlotPileupDistribution (fGraphFR, "No. primary vertices", "F_{R}", "BLA", 0.5, 9.5, -0.30, 0.30, InputChannel, "", "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PileupDependence_FR_"+InputChannel+"_"+FittingMode+".eps");
       
       return 0;
      
     }
     else if(Systematic == "MuEval"){
       
       TGraphErrors *fGraphF0 = new TGraphErrors(6);
       TGraphErrors *fGraphFL = new TGraphErrors(6);
       TGraphErrors *fGraphFR = new TGraphErrors(6);

       for(int k = 4; k < 10; k++){
	 
	 std::stringstream oss;
	 oss << k;
	 
	 TFile *fNewFile = new TFile(("ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/Mu_"+oss.str()+"/DatafitCorrectedUnc_"+InputChannel+".root").c_str(), "READ");
	 
	 TTree *tree  = (TTree*) fNewFile -> Get("EnsembleTree");
	 
	 double F0,     FL,     FR;
	 double F0err1, FLerr1, FRerr1;
	 double  F0err,  FLerr,  FRerr;
	 
	 tree -> SetBranchAddress("F0",        &F0);
	 tree -> SetBranchAddress("FL",        &FL);
	 tree -> SetBranchAddress("FR",        &FR);
	 
	 tree -> SetBranchAddress("F0_err", &F0err);
	 tree -> SetBranchAddress("FL_err", &FLerr);
	 tree -> SetBranchAddress("FR_err", &FRerr);
	 
	 tree -> GetEntry(0);

	 std::cout << F0    << "\t" << FL    << "\t" << FR    << std::endl;
	 std::cout << F0err << "\t" << FLerr << "\t" << FRerr << std::endl;
	 std::cout << "\t"  << std::endl;

	 fGraphF0 -> SetPoint(k-4, 2*k-4, F0);
	 fGraphF0 -> SetPointError(k-4, 0, F0err);

	 fGraphFL -> SetPoint(k-4, 2*k-4, FL);
         fGraphFL -> SetPointError(k-4, 0, FLerr);

	 fGraphFR -> SetPoint(k-4, 2*k-4, FR);
         fGraphFR -> SetPointError(k-4, 0, FRerr);

	 fNewFile -> Close();


       }

       PlotPileupDistribution (fGraphF0, "<#mu>", "F_{0}", "BLA", 3.5, 15.5, 0.00, 1.30, InputChannel, "", "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PileupDependence_Mu_F0_"+InputChannel+"_"+FittingMode+".eps");       
       PlotPileupDistribution (fGraphFL, "<#mu>", "F_{L}", "BLA", 3.5, 15.5, -0.20, 0.80, InputChannel, "", "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PileupDependence_Mu_FL_"+InputChannel+"_"+FittingMode+".eps");
       
       if(FittingMode == "3D")
	 PlotPileupDistribution (fGraphFR, "<#mu>", "F_{R}", "BLA", 3.5, 15.5, -0.30, 0.30, InputChannel, "", "ExternalPileupOutput_"+InputChannel+"_"+FittingMode+"/PileupDependence_Mu_FR_"+InputChannel+"_"+FittingMode+".eps");
       
       return 0;
      
     }

 
}
