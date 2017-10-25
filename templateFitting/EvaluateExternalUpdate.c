#include "StatusLogbook.h"
#include "ProfilingClass.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <sys/stat.h>

using namespace std;

std::vector<std::string> ReadInputFiles(const char * filename);
void runUpDownSystematicVariation( ProfilingClass &fProf, std::string TemplateFile, std::string Systematic, double m_f0, double m_fL, double m_fR, std::string pseudoData, std::string wjetsMode);
void runSingleSystematicVariation( ProfilingClass &fProf, std::string TemplateFile, std::string Systematic, double m_f0, double m_fL, double m_fR, std::string wjetsMode);
void runSinglePseudoData( ProfilingClass &fProf, std::string TemplateFile, std::string Systematic, double m_f0, double m_fL, double m_fR, std::string wjetsMode);
std::string add3Wtag (std::string &s, const std::string &toReplace, const std::string &replaceWith);

int main(int argc, char *argv[])
{

  std::cout<<"n_args: "<<argc<<std::endl;
     if(argc != 11){

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
     std::string FittingMode    = argv[4]; // 2D , 3D
     std::string NewSample    = argv[5]; // true , false
     std::string CalibrationMode    = argv[6]; // single , dual
     std::string tagMode   = argv[7]; // 1excl ,1incl, 2incl
     std::string angleMode   = argv[8]; // lep, had, lephad
     std::string fixfloat_bkg   = argv[9]; // fix, float
     std::string wjetsMode   = argv[10]; // 1W, 3W

     std::string TemplatesMainDir   = "new_02Aug2016_LHcut";
     //std::string TemplatesMainDir   = "new_11Mar2016_LHcut_CFerr";
     //std::string TemplatesMainDir   = "new_11Mar2016_LHcut_2CFerr";
     
     if(!(string(tagMode) == "1excl" || string(tagMode) == "1incl" || string(tagMode) == "2incl" ||string(tagMode)=="1excl2incl" ))
     {
       WriteErrorStatus("EvaluateExternal", "Unknown b-tag. Use 1excl or 2incl !!!");
       WriteErrorStatus("EvaluateExternal", "EXIT   ");
       return 1;
     }

     
     if(!(string(angleMode) == "lep" || string(angleMode) == "had" || string(angleMode) == "lephad"))
     {
       WriteErrorStatus("EvaluateExternal", "Unknown angle mode. lep or had are lephad allowed !!!");
       WriteErrorStatus("EvaluateExternal", "EXIT   ");
       return 1;
     }

     if(!(string(fixfloat_bkg) == "fix" || string(fixfloat_bkg) == "float"))
     {
       WriteErrorStatus("EvaluateExternal", "Unknown bkg. normalisation treatment mode. fix or float are allowed !!!");
       WriteErrorStatus("EvaluateExternal", "EXIT   ");
       return 1;
     }

     if(!(string(wjetsMode) == "1W" || string(wjetsMode) == "3W"))
     {
       WriteErrorStatus("EvaluateExternal", "Unknown wjets mode. 1W or 3W 2incl !!!");
       WriteErrorStatus("EvaluateExternal", "EXIT   ");
       return 1;
     }

//     int NumberOfBins = 21;
     int NumberOfBins = 15;
     double lower_edge = -1.0;
     double upper_edge =  1.0;
     
     bool useNewSample = true; // for 110404
     std::string KLFoption="KLF5jOPT"; // KLF5jets , KLF5jOPT

     string TemplatesDir,TemplatesSubDir;
     string uncSuffix;
     bool fixedBkg = false;
     if(fixfloat_bkg == "fix") fixedBkg = true;
     

     double m_fL, m_f0, m_fR;

     double fL_el= 0.302763, f0_el= 0.695922, fR_el= 0.00131471, fL_mu= 0.302914, f0_mu= 0.695735, fR_mu= 0.0013513;
     
     if(InputChannel == "el" || InputChannel == "el_BTag"){
        m_fL=fL_el ; m_f0=f0_el; m_fR=fR_el;
     }
     else if (InputChannel == "mu" || InputChannel == "mu_BTag"){
        m_fL=fL_mu ; m_f0=f0_mu; m_fR=fR_mu;
     }
     else{ // using the avarage for the leptonic combined channels
        m_fL= 0.5*(fL_mu + fL_el) ; m_f0=0.5*(f0_mu + f0_el); m_fR=0.5*(fR_mu + fR_el);
     }

     if(angleMode=="lep")
      TemplatesSubDir = "/Leptonic";
    else if(angleMode=="had")
      TemplatesSubDir = "/Hadronic";
    else
      TemplatesSubDir = "/LepHad";

     if(Systematic =="Calibration")
       //TemplatesDir = "/afs/cern.ch/work/m/mkareem/Wpol_mc/templateMaker/TemplateFiles/systTemplates/systTemplates_110404/"+TemplatesMainDir+TemplatesSubDir;
       TemplatesDir = "/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/new_02Aug2016_LHcut"+TemplatesSubDir;
     else
       //TemplatesDir = "/afs/cern.ch/work/m/mkareem/Wpol_mc/templateMaker/TemplateFiles/systTemplates/systTemplates_110404/"+TemplatesMainDir+TemplatesSubDir; 
       TemplatesDir = "/afs/cern.ch/work/b/btannenw/ttbar_Whelicity5/AnalysisTop-1.12.0/new_02Aug2016_LHcut"+TemplatesSubDir;
  


     if(NewSample=="false") {useNewSample=false; std::cout<< "### Used sample: 117050" <<std::endl;}
     else std::cout<< "### Used sample: 110404" <<std::endl;

     std::cout<< "### CalibrationMode: " << CalibrationMode <<std::endl;
     std::cout<< "### tagMode: " << tagMode <<std::endl;
     std::cout<< "### angleMode: " << angleMode <<std::endl;
     std::cout<< "### fixfloat_bkg: " << fixfloat_bkg <<std::endl;
     std::cout<< "### wjetsMode: " << wjetsMode <<std::endl;

     std::string SignalSampleNumber;
     std::string OtherSignalSampleNumber;
     

     if(useNewSample)
      {
        SignalSampleNumber = "110404";
        OtherSignalSampleNumber = "117050";
      }
     else
     {
        SignalSampleNumber = "117050";
        OtherSignalSampleNumber = "110404";
     }
    
     if(CalibrationMode=="single")
       OtherSignalSampleNumber = SignalSampleNumber;

     // 25 Feb. 2016 bt: Commented out since pdf uncertainties handled differently now. left structure to minimize changes
     bool IsPDF = false;
     /*
     // check if Systematic is part of PDF
     size_t found;
     found = Systematic.find("PDF");
     if(found !=string::npos)
       IsPDF = true;
     */

     std::string FitType; // BkgSum, WSum, RemBkgSum,QCDSum,All
     FitType = "All";
     //FitType = "RemBkgSum";

     

     std::string TemplateFile;
     std::string oldSignalTemplateFile;

     //TemplateFile = "Templates_Combined_"+InputChannel+"_Powheg_2.root"; // 7 TeV
          
     if(useNewSample){
       if(Systematic =="Calibration")
	 TemplateFile = TemplatesDir+"/Templates_110404_syst_"+tagMode+"_"+KLFoption+"_"+InputChannel+".root";
       else
        TemplateFile = TemplatesDir+"/Templates_110404_syst_"+tagMode+"_"+KLFoption+"_"+InputChannel+".root";
     }
     else{
       if(Systematic =="Calibration")
	 TemplateFile = TemplatesDir+"/Templates_117050_syst_"+tagMode+"_"+KLFoption+"_"+InputChannel+".root";
       else
        TemplateFile = TemplatesDir+"/Templates_117050_syst_"+tagMode+"_"+KLFoption+"_"+InputChannel+".root";
     }
     if(InputChannel == "el_mu"){ // 2 channels
       
       NumberOfBins =   NumberOfBins*2;
       lower_edge   = -2.0;
       upper_edge   =  2.0;

     }
     else if(InputChannel == "el_BTag" || InputChannel == "mu_BTag"){ // 2 channels

       NumberOfBins =   NumberOfBins*2;
       lower_edge   = -2.0;
       upper_edge   =  2.0;

     }
     else if(InputChannel == "el_mu_bTag"){ // 4 channels

       NumberOfBins =   NumberOfBins*4;
       lower_edge   = -4.0;
       upper_edge   =  4.0;

     }
     else if(InputChannel == "el_mu_lephad"){ // 4 channels

       NumberOfBins =   NumberOfBins*4;
       lower_edge   = -4.0;
       upper_edge   =  4.0;

     }

     else if(InputChannel == "el_mu_lephad_bTag"){ // 8 channels

       NumberOfBins =   NumberOfBins*8;
       lower_edge   = -8.0;
       upper_edge   =  8.0;

     }

     std::vector<std::string> VarNames;
     std::vector<std::string> HistoNames;

     //Globally defined histograms, needed for fit function
     hist      = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
     hist_sum  = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);
     hist_help = TH1D("hist", "hist", NumberOfBins, lower_edge, upper_edge);

     ProfilingClass fProf = ProfilingClass();
//---------------------------------------------------------------------------------------------------------


// Option for running the bkg as fixed or floating in the fitter (17.01.2016 by MJK)
if(fixedBkg) uncSuffix = "Fixed";
else uncSuffix ="";
//------------------------------------------
     if(InputChannel == "el" || InputChannel == "mu"){
       if(InputChannel == "el"){
          if(tagMode=="2incl") fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "#geq 2>= b-tags");
          else if (tagMode == "1excl") fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "1 b-tags");
       }
        
       else if(InputChannel == "mu"){
          if(tagMode=="2incl") fProf.SetSubplotLabels("#mu+jets",   "#geq 4 jets", "#geq 2>= b-tags");
          else if (tagMode == "1excl") fProf.SetSubplotLabels("#mu+jets",   "#geq 4 jets", "1 b-tags");
       }
       // add signal parameter with param name and hist name

       
       fProf.AddSignalParameter("N0",  "CosTheta_F0", "F0eff");
       fProf.AddSignalParameter("NL",  "CosTheta_FL", "FLeff");
       if(FittingMode == "3D")
	       fProf.AddSignalParameter("NR",  "CosTheta_FR", "FReff");
       
       
       
       // add bkg parameter with param name and hist name
       if(wjetsMode=="1W")
	       fProf.AddBackgroundParameter("Wjets",  "Wjets",  "WjetsNorm",  "WjetsUnc"+uncSuffix);
       else if(wjetsMode=="3W"){
	       fProf.AddBackgroundParameter("WLight", "WLight", "WLightNorm",  "WLightUnc"+uncSuffix);
       	       fProf.AddBackgroundParameter("Wc",     "Wc",     "WcNorm",     "WcUnc"+uncSuffix);
	       fProf.AddBackgroundParameter("Wbbcc",  "Wbbcc",  "WbbccNorm",  "WbbccUnc"+uncSuffix);
       }
       fProf.AddBackgroundParameter("QCD",    "QCD",    "QCDNorm",    "QCDUnc"+uncSuffix);
       fProf.AddBackgroundParameter("RemBkg", "RemBkg", "RemBkgNorm", "RemBkgUnc"+uncSuffix);
       
      
       HistoNames.push_back("CosTheta_F0");
       HistoNames.push_back("CosTheta_FL");
       if(FittingMode == "3D")
	       HistoNames.push_back("CosTheta_FR");
      
      
       if(wjetsMode=="1W")
	         HistoNames.push_back("Wjets");
       else if (wjetsMode=="3W"){
	       HistoNames.push_back("WLight");
	       HistoNames.push_back("Wc");
	       HistoNames.push_back("Wbbcc"); 
       }

       HistoNames.push_back("QCD");
       HistoNames.push_back("RemBkg");

       
      
       VarNames.push_back("F0eff");
       VarNames.push_back("FLeff");
       if(FittingMode == "3D")
	       VarNames.push_back("FReff");
      

       if(wjetsMode=="1W")
	       VarNames.push_back("WjetsNorm");
       else if (wjetsMode=="3W"){
	       VarNames.push_back("WLightNorm");
	       VarNames.push_back("WcNorm");
	       VarNames.push_back("WbbccNorm"); 
       }
       VarNames.push_back("QCDNorm");
       VarNames.push_back("RemBkgNorm");
       
     }
     //-------------------------------------------------------------------------------------------------------
     else if (InputChannel == "el_mu"){

       fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "#geq 2>= b-tags");
       fProf.SetSubplotLabels("#mu+jets", "#geq 4 jets", "#geq 2>= b-tags");

       // add signal parameter with param name and hist name
       fProf.AddSignalParameter("N0",  "CosTheta_F0_allCombined", "F0eff_allCombined");
       fProf.AddSignalParameter("NL",  "CosTheta_FL_allCombined", "FLeff_allCombined");
       if(FittingMode == "3D")
         fProf.AddSignalParameter("NR",  "CosTheta_FR_allCombined", "FReff_allCombined");

      if(tagMode=="1excl"){

       fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "#geq 1 b-tags");
       fProf.SetSubplotLabels("#mu+jets", "#geq 4 jets", "#geq 1 b-tags");

       if(FitType == "All"){
       // add bkg parameter with param name and hist name
	 if(wjetsMode=="1W"){
	   fProf.AddBackgroundParameter("Wjets_4incl_1excl_el",  "Wjets_4incl_1excl_el",  "WjetsNorm_4incl_1excl_el",  "WjetsUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("Wjets_4incl_1excl_mu",  "Wjets_4incl_1excl_mu",  "WjetsNorm_4incl_1excl_mu",  "WjetsUnc"+uncSuffix+"_4incl_1excl_mu");
	 }
	 else if(wjetsMode=="3W"){
	   fProf.AddBackgroundParameter("WLight_4incl_1excl_el",  "WLight_4incl_1excl_el",  "WLightNorm_4incl_1excl_el",  "WLightUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("WLight_4incl_1excl_mu",  "WLight_4incl_1excl_mu",  "WLightNorm_4incl_1excl_mu",  "WLightUnc"+uncSuffix+"_4incl_1excl_mu");
	   fProf.AddBackgroundParameter("Wc_4incl_1excl_el",  "Wc_4incl_1excl_el",  "WcNorm_4incl_1excl_el",  "WcUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("Wc_4incl_1excl_mu",  "Wc_4incl_1excl_mu",  "WcNorm_4incl_1excl_mu",  "WcUnc"+uncSuffix+"_4incl_1excl_mu");
	   fProf.AddBackgroundParameter("Wbbcc_4incl_1excl_el",  "Wbbcc_4incl_1excl_el",  "WbbccNorm_4incl_1excl_el",  "WbbccUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("Wbbcc_4incl_1excl_mu",  "Wbbcc_4incl_1excl_mu",  "WbbccNorm_4incl_1excl_mu",  "WbbccUnc"+uncSuffix+"_4incl_1excl_mu");
	 } 
	 fProf.AddBackgroundParameter("QCD_4incl_1excl_el",  "QCD_4incl_1excl_el",  "QCDNorm_4incl_1excl_el",  "QCDUnc"+uncSuffix+"_4incl_1excl_el");
	 fProf.AddBackgroundParameter("QCD_4incl_1excl_mu",  "QCD_4incl_1excl_mu",  "QCDNorm_4incl_1excl_mu",  "QCDUnc"+uncSuffix+"_4incl_1excl_mu"); 
	 fProf.AddBackgroundParameter("RemBkg_4incl_1excl_el",  "RemBkg_4incl_1excl_el",  "RemBkgNorm_4incl_1excl_el",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_el");
	 fProf.AddBackgroundParameter("RemBkg_4incl_1excl_mu",  "RemBkg_4incl_1excl_mu",  "RemBkgNorm_4incl_1excl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_mu");
       
        }

       HistoNames.push_back("CosTheta_F0_allCombined");
       HistoNames.push_back("CosTheta_FL_allCombined");
       if(FittingMode == "3D")
	 HistoNames.push_back("CosTheta_FR_allCombined");
       if(wjetsMode=="1W"){
	 HistoNames.push_back("Wjets_4incl_1excl_el");
	 HistoNames.push_back("Wjets_4incl_1excl_mu");
       }
       else if(wjetsMode=="3W"){
	 HistoNames.push_back("WLight_4incl_1excl_el");
	 HistoNames.push_back("WLight_4incl_1excl_mu");
	 HistoNames.push_back("Wc_4incl_1excl_el");
	 HistoNames.push_back("Wc_4incl_1excl_mu");
	 HistoNames.push_back("Wbbcc_4incl_1excl_el");
	 HistoNames.push_back("Wbbcc_4incl_1excl_mu");
       }
       HistoNames.push_back("QCD_4incl_1excl_el");
       HistoNames.push_back("QCD_4incl_1excl_mu");
       HistoNames.push_back("RemBkg_4incl_1excl_el");
       HistoNames.push_back("RemBkg_4incl_1excl_mu");

       VarNames.push_back("F0eff_allCombined");
       VarNames.push_back("FLeff_allCombined");
       if(FittingMode == "3D")
	       VarNames.push_back("FReff_allCombined");
       
       if(wjetsMode=="1W"){
	 VarNames.push_back("WjetsNorm_4incl_1excl_el");
	 VarNames.push_back("WjetsNorm_4incl_1excl_mu");
       }
       else if (wjetsMode=="3W"){
	 VarNames.push_back("WLightNorm_4incl_1excl_el");
	 VarNames.push_back("WLightNorm_4incl_1excl_mu");
	 VarNames.push_back("WcNorm_4incl_1excl_el");
	 VarNames.push_back("WcNorm_4incl_1excl_mu");
	 VarNames.push_back("WbbccNorm_4incl_1excl_el");
	 VarNames.push_back("WbbccNorm_4incl_1excl_mu");
       }
       VarNames.push_back("QCDNorm_4incl_1excl_el");
       VarNames.push_back("QCDNorm_4incl_1excl_mu");
       VarNames.push_back("RemBkgNorm_4incl_1excl_el");
       VarNames.push_back("RemBkgNorm_4incl_1excl_mu");

     }
     else if (tagMode == "2incl"){
       if(FitType == "All"){
	 // add bkg parameter with param name and hist name
	 if(wjetsMode=="1W"){
	   // fProf.AddBackgroundParameter("Wjets_4incl_2incl_el",  "Wjets_4incl_2incl_el",  "WjetsNorm_4incl_2incl_el",  "WjetsUnc"+uncSuffix+"_4incl_2incl_el");
	   // fProf.AddBackgroundParameter("Wjets_4incl_2incl_mu",  "Wjets_4incl_2incl_mu",  "WjetsNorm_4incl_2incl_mu",  "WjetsUnc"+uncSuffix+"_4incl_2incl_mu");
    fProf.AddBackgroundParameter("Wjets_4incl_2incl",  "Wjets_4incl_2incl",  "WjetsNorm_4incl_2incl",  "WjetsUnc"+uncSuffix+"_4incl_2incl");
	 }
	 else if(wjetsMode=="3W"){
	   // fProf.AddBackgroundParameter("WLight_4incl_2incl_el",  "WLight_4incl_2incl_el",  "WLightNorm_4incl_2incl_el",  "WLightUnc"+uncSuffix+"_4incl_2incl_el");
	   // fProf.AddBackgroundParameter("WLight_4incl_2incl_mu",  "WLight_4incl_2incl_mu",  "WLightNorm_4incl_2incl_mu",  "WLightUnc"+uncSuffix+"_4incl_2incl_mu");
    fProf.AddBackgroundParameter("WLight_4incl_2incl",  "WLight_4incl_2incl",  "WLightNorm_4incl_2incl",  "WLightUnc"+uncSuffix+"_4incl_2incl");
	   // fProf.AddBackgroundParameter("Wc_4incl_2incl_el",  "Wc_4incl_2incl_el",  "WcNorm_4incl_2incl_el",  "WcUnc"+uncSuffix+"_4incl_2incl_el");
	   // fProf.AddBackgroundParameter("Wc_4incl_2incl_mu",  "Wc_4incl_2incl_mu",  "WcNorm_4incl_2incl_mu",  "WcUnc"+uncSuffix+"_4incl_2incl_mu");
    fProf.AddBackgroundParameter("Wc_4incl_2incl",  "Wc_4incl_2incl",  "WcNorm_4incl_2incl",  "WcUnc"+uncSuffix+"_4incl_2incl");
	   // fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_el",  "Wbbcc_4incl_2incl_el",  "WbbccNorm_4incl_2incl_el",  "WbbccUnc"+uncSuffix+"_4incl_2incl_el");
	   // fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_mu",  "Wbbcc_4incl_2incl_mu",  "WbbccNorm_4incl_2incl_mu",  "WbbccUnc"+uncSuffix+"_4incl_2incl_mu");
    fProf.AddBackgroundParameter("Wbbcc_4incl_2incl",  "Wbbcc_4incl_2incl",  "WbbccNorm_4incl_2incl",  "WbbccUnc"+uncSuffix+"_4incl_2incl");
	 }
	 fProf.AddBackgroundParameter("QCD_4incl_2incl_el",  "QCD_4incl_2incl_el",  "QCDNorm_4incl_2incl_el",  "QCDUnc"+uncSuffix+"_4incl_2incl_el");
	 fProf.AddBackgroundParameter("QCD_4incl_2incl_mu",  "QCD_4incl_2incl_mu",  "QCDNorm_4incl_2incl_mu",  "QCDUnc"+uncSuffix+"_4incl_2incl_mu"); 
	 // fProf.AddBackgroundParameter("RemBkg_4incl_2incl_el",  "RemBkg_4incl_2incl_el",  "RemBkgNorm_4incl_2incl_el",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_el");
	 // fProf.AddBackgroundParameter("RemBkg_4incl_2incl_mu",  "RemBkg_4incl_2incl_mu",  "RemBkgNorm_4incl_2incl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_mu");
   fProf.AddBackgroundParameter("RemBkg_4incl_2incl",  "RemBkg_4incl_2incl",  "RemBkgNorm_4incl_2incl",  "RemBkgUnc"+uncSuffix+"_4incl_2incl");
	 
       }
       
       HistoNames.push_back("CosTheta_F0_allCombined");
       HistoNames.push_back("CosTheta_FL_allCombined");
       if(FittingMode == "3D")
         HistoNames.push_back("CosTheta_FR_allCombined");
       
       if(wjetsMode=="1W"){
	       // HistoNames.push_back("Wjets_4incl_2incl_el");
	       // HistoNames.push_back("Wjets_4incl_2incl_mu");
        HistoNames.push_back("Wjets_4incl_2incl");
       }
       else if(wjetsMode=="3W"){
	       // HistoNames.push_back("WLight_4incl_2incl_el");
	       // HistoNames.push_back("WLight_4incl_2incl_mu");
        HistoNames.push_back("WLight_4incl_2incl");
	       // HistoNames.push_back("Wc_4incl_2incl_el");
	       // HistoNames.push_back("Wc_4incl_2incl_mu");
        HistoNames.push_back("Wc_4incl_2incl");
	       // HistoNames.push_back("Wbbcc_4incl_2incl_el");
	       // HistoNames.push_back("Wbbcc_4incl_2incl_mu");
        HistoNames.push_back("Wbbcc_4incl_2incl");
       }
       HistoNames.push_back("QCD_4incl_2incl_el");
       HistoNames.push_back("QCD_4incl_2incl_mu");
       // HistoNames.push_back("RemBkg_4incl_2incl_el");
       // HistoNames.push_back("RemBkg_4incl_2incl_mu");
       HistoNames.push_back("RemBkg_4incl_2incl");

       VarNames.push_back("F0eff_allCombined");
       VarNames.push_back("FLeff_allCombined");
       if(FittingMode == "3D")
         VarNames.push_back("FReff_allCombined");
       
       if(wjetsMode=="1W"){
	 // VarNames.push_back("WjetsNorm_4incl_2incl_el");
	 // VarNames.push_back("WjetsNorm_4incl_2incl_mu");
        VarNames.push_back("WjetsNorm_4incl_2incl");
       }
       if(wjetsMode=="3W"){
	       // VarNames.push_back("WLightNorm_4incl_2incl_el");
	       // VarNames.push_back("WLightNorm_4incl_2incl_mu");
        VarNames.push_back("WLightNorm_4incl_2incl");
	       // VarNames.push_back("WcNorm_4incl_2incl_el");
	       // VarNames.push_back("WcNorm_4incl_2incl_mu");
        VarNames.push_back("WcNorm_4incl_2incl");
	       // VarNames.push_back("WbbccNorm_4incl_2incl_el");
	       // VarNames.push_back("WbbccNorm_4incl_2incl_mu");
        VarNames.push_back("WbbccNorm_4incl_2incl");
       }

       VarNames.push_back("QCDNorm_4incl_2incl_el");
       VarNames.push_back("QCDNorm_4incl_2incl_mu");
       // VarNames.push_back("RemBkgNorm_4incl_2incl_el");
       // VarNames.push_back("RemBkgNorm_4incl_2incl_mu");
       VarNames.push_back("RemBkgNorm_4incl_2incl");

     }
   }
   //---------------------------------------------------------------------  

     else if (InputChannel == "el_BTag" || InputChannel == "mu_BTag"){

       if(InputChannel == "el_BTag"){
	       //fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "0 b-tags");
	       fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "1 b-tags");
	       fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "#geq 2 b-tags");
       }
       else{
	       //fProf.SetSubplotLabels("#mu+jets",   "#geq 4 jets", "0 b-tags");
         fProf.SetSubplotLabels("#mu+jets",   "#geq 4 jets", "1 b-tags");
         fProf.SetSubplotLabels("#mu+jets",   "#geq 4 jets", "#geq 2 b-tags");
       }

       // add signal parameter with param name and hist name
       fProf.AddSignalParameter("N0",  "CosTheta_F0_allCombined", "F0eff_allCombined");
       fProf.AddSignalParameter("NL",  "CosTheta_FL_allCombined", "FLeff_allCombined");
       if(FittingMode == "3D")
         fProf.AddSignalParameter("NR",  "CosTheta_FR_allCombined", "FReff_allCombined");

       // add bkg parameter with param name and hist name

       if(InputChannel == "el_BTag"){

	 if(FitType == "All"){

	   //fProf.AddBackgroundParameter("Wjets_4incl_0excl_el",   "Wjets_4incl_0excl_el",   "WjetsNorm_4incl_0excl_el",   "WjetsUnc"+uncSuffix+"_4incl_0excl_el");
	   if(wjetsMode=="1W"){
	     fProf.AddBackgroundParameter("Wjets_4incl_1excl_el",   "Wjets_4incl_1excl_el",   "WjetsNorm_4incl_1excl_el",   "WjetsUnc"+uncSuffix+"_4incl_1excl_el");
	     fProf.AddBackgroundParameter("Wjets_4incl_2incl_el",   "Wjets_4incl_2incl_el",   "WjetsNorm_4incl_2incl_el",   "WjetsUnc"+uncSuffix+"_4incl_2incl_el"); 
	   }
	   else if(wjetsMode=="3W"){
	     fProf.AddBackgroundParameter("WLight_4incl_1excl_el",   "WLight_4incl_1excl_el",   "WLightNorm_4incl_1excl_el",   "WLightUnc"+uncSuffix+"_4incl_1excl_el");
	     fProf.AddBackgroundParameter("WLight_4incl_2incl_el",   "WLight_4incl_2incl_el",   "WLightNorm_4incl_2incl_el",   "WLightUnc"+uncSuffix+"_4incl_2incl_el"); 
	     fProf.AddBackgroundParameter("Wc_4incl_1excl_el",   "Wc_4incl_1excl_el",   "WcNorm_4incl_1excl_el",   "WcUnc"+uncSuffix+"_4incl_1excl_el");
	     fProf.AddBackgroundParameter("Wc_4incl_2incl_el",   "Wc_4incl_2incl_el",   "WcNorm_4incl_2incl_el",   "WcUnc"+uncSuffix+"_4incl_2incl_el"); 
	     fProf.AddBackgroundParameter("Wbbcc_4incl_1excl_el",   "Wbbcc_4incl_1excl_el",   "WbbccNorm_4incl_1excl_el",   "WbbccUnc"+uncSuffix+"_4incl_1excl_el");
	     fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_el",   "Wbbcc_4incl_2incl_el",   "WbbccNorm_4incl_2incl_el",   "WbbccUnc"+uncSuffix+"_4incl_2incl_el"); 
	   }
	   //fProf.AddBackgroundParameter("QCD_4incl_0excl_el",     "QCD_4incl_0excl_el",     "QCDNorm_4incl_0excl_el",     "QCDUnc"+uncSuffix+"_4incl_0excl_el");
	   fProf.AddBackgroundParameter("QCD_4incl_1excl_el",     "QCD_4incl_1excl_el",     "QCDNorm_4incl_1excl_el",     "QCDUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("QCD_4incl_2incl_el",     "QCD_4incl_2incl_el",     "QCDNorm_4incl_2incl_el",     "QCDUnc"+uncSuffix+"_4incl_2incl_el");
	   //fProf.AddBackgroundParameter("RemBkg_4incl_0excl_el",  "RemBkg_4incl_0excl_el",  "RemBkgNorm_4incl_0excl_el",  "RemBkgUnc"+uncSuffix+"_4incl_0excl_el");
	   fProf.AddBackgroundParameter("RemBkg_4incl_1excl_el",  "RemBkg_4incl_1excl_el",  "RemBkgNorm_4incl_1excl_el",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("RemBkg_4incl_2incl_el",  "RemBkg_4incl_2incl_el",  "RemBkgNorm_4incl_2incl_el",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_el");
	 }
	 else if(FitType == "WSum"){
	   
	   fProf.AddBackgroundParameter("Wjets_allCombined",   "Wjets_allCombined",   "WjetsNorm_allCombined",   "WjetsUnc"+uncSuffix+"_allCombined");
	   
	   fProf.AddBackgroundParameter("QCD_4incl_0excl_el",     "QCD_4incl_0excl_el",     "QCDNorm_4incl_0excl_el",     "QCDUnc"+uncSuffix+"_4incl_0excl_el");
	   fProf.AddBackgroundParameter("QCD_4incl_1excl_el",     "QCD_4incl_1excl_el",     "QCDNorm_4incl_1excl_el",     "QCDUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("QCD_4incl_2incl_el",     "QCD_4incl_2incl_el",     "QCDNorm_4incl_2incl_el",     "QCDUnc"+uncSuffix+"_4incl_2incl_el");
	   fProf.AddBackgroundParameter("RemBkg_4incl_0excl_el",  "RemBkg_4incl_0excl_el",  "RemBkgNorm_4incl_0excl_el",  "RemBkgUnc"+uncSuffix+"_4incl_0excl_el");
	   fProf.AddBackgroundParameter("RemBkg_4incl_1excl_el",  "RemBkg_4incl_1excl_el",  "RemBkgNorm_4incl_1excl_el",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_el");
	   fProf.AddBackgroundParameter("RemBkg_4incl_2incl_el",  "RemBkg_4incl_2incl_el",  "RemBkgNorm_4incl_2incl_el",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_el");
	 
	 }
	 else if(FitType == "BkgSum"){

	   fProf.AddBackgroundParameter("Wjets_allCombined",  "Wjets_allCombined",  "WjetsNorm_allCombined",   "WjetsUnc"+uncSuffix+"_allCombined");
	   fProf.AddBackgroundParameter("QCD_allCombined",    "QCD_allCombined",    "QCDNorm_allCombined",     "QCDUnc"+uncSuffix+"_allCombined");
	   fProf.AddBackgroundParameter("RemBkg_allCombined", "RemBkg_allCombined", "RemBkgNorm_allCombined",  "RemBkgUnc"+uncSuffix+"_allCombined");

	 }
	 else if(FitType == "RemBkgSum"){

           fProf.AddBackgroundParameter("Wjets_4incl_0excl_el",   "Wjets_4incl_0excl_el",   "WjetsNorm_4incl_0excl_el",   "WjetsUnc"+uncSuffix+"_4incl_0excl_el");
           fProf.AddBackgroundParameter("Wjets_4incl_1excl_el",   "Wjets_4incl_1excl_el",   "WjetsNorm_4incl_1excl_el",   "WjetsUnc"+uncSuffix+"_4incl_1excl_el");
           fProf.AddBackgroundParameter("Wjets_4incl_2incl_el",   "Wjets_4incl_2incl_el",   "WjetsNorm_4incl_2incl_el",   "WjetsUnc"+uncSuffix+"_4incl_2incl_el");
           fProf.AddBackgroundParameter("QCD_4incl_0excl_el",     "QCD_4incl_0excl_el",     "QCDNorm_4incl_0excl_el",     "QCDUnc"+uncSuffix+"_4incl_0excl_el");
           fProf.AddBackgroundParameter("QCD_4incl_1excl_el",     "QCD_4incl_1excl_el",     "QCDNorm_4incl_1excl_el",     "QCDUnc"+uncSuffix+"_4incl_1excl_el");
           fProf.AddBackgroundParameter("QCD_4incl_2incl_el",     "QCD_4incl_2incl_el",     "QCDNorm_4incl_2incl_el",     "QCDUnc"+uncSuffix+"_4incl_2incl_el");

           fProf.AddBackgroundParameter("RemBkg_allCombined",     "RemBkg_allCombined",     "RemBkgNorm_allCombined",     "RemBkgUnc"+uncSuffix+"_allCombined");

         }

       }
       else{ // if mu_Btag
	 
	 if(FitType == "All"){

	   //fProf.AddBackgroundParameter("Wjets_4incl_0excl_mu",   "Wjets_4incl_0excl_mu",   "WjetsNorm_4incl_0excl_mu",   "WjetsUnc"+uncSuffix+"_4incl_0excl_mu");
	   if(wjetsMode=="1W"){
	     fProf.AddBackgroundParameter("Wjets_4incl_1excl_mu",   "Wjets_4incl_1excl_mu",   "WjetsNorm_4incl_1excl_mu",   "WjetsUnc"+uncSuffix+"_4incl_1excl_mu");
	     fProf.AddBackgroundParameter("Wjets_4incl_2incl_mu",   "Wjets_4incl_2incl_mu",   "WjetsNorm_4incl_2incl_mu",   "WjetsUnc"+uncSuffix+"_4incl_2incl_mu"); 
	   }
	   else if(wjetsMode=="3W"){
	     fProf.AddBackgroundParameter("WLight_4incl_1excl_mu",   "WLight_4incl_1excl_mu",   "WLightNorm_4incl_1excl_mu",   "WLightUnc"+uncSuffix+"_4incl_1excl_mu");
	     fProf.AddBackgroundParameter("WLight_4incl_2incl_mu",   "WLight_4incl_2incl_mu",   "WLightNorm_4incl_2incl_mu",   "WLightUnc"+uncSuffix+"_4incl_2incl_mu"); 
	     fProf.AddBackgroundParameter("Wc_4incl_1excl_mu",   "Wc_4incl_1excl_mu",   "WcNorm_4incl_1excl_mu",   "WcUnc"+uncSuffix+"_4incl_1excl_mu");
	     fProf.AddBackgroundParameter("Wc_4incl_2incl_mu",   "Wc_4incl_2incl_mu",   "WcNorm_4incl_2incl_mu",   "WcUnc"+uncSuffix+"_4incl_2incl_mu"); 
	     fProf.AddBackgroundParameter("Wbbcc_4incl_1excl_mu",   "Wbbcc_4incl_1excl_mu",   "WbbccNorm_4incl_1excl_mu",   "WbbccUnc"+uncSuffix+"_4incl_1excl_mu");
	     fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_mu",   "Wbbcc_4incl_2incl_mu",   "WbbccNorm_4incl_2incl_mu",   "WbbccUnc"+uncSuffix+"_4incl_2incl_mu"); 
	   }
	   //fProf.AddBackgroundParameter("QCD_4incl_0excl_mu",     "QCD_4incl_0excl_mu",     "QCDNorm_4incl_0excl_mu",     "QCDUnc"+uncSuffix+"_4incl_0excl_mu");
	   fProf.AddBackgroundParameter("QCD_4incl_1excl_mu",     "QCD_4incl_1excl_mu",     "QCDNorm_4incl_1excl_mu",     "QCDUnc"+uncSuffix+"_4incl_1excl_mu");
	   fProf.AddBackgroundParameter("QCD_4incl_2incl_mu",     "QCD_4incl_2incl_mu",     "QCDNorm_4incl_2incl_mu",     "QCDUnc"+uncSuffix+"_4incl_2incl_mu");
	   //fProf.AddBackgroundParameter("RemBkg_4incl_0excl_mu",  "RemBkg_4incl_0excl_mu",  "RemBkgNorm_4incl_0excl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_0excl_mu");
	   fProf.AddBackgroundParameter("RemBkg_4incl_1excl_mu",  "RemBkg_4incl_1excl_mu",  "RemBkgNorm_4incl_1excl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_mu");
	   fProf.AddBackgroundParameter("RemBkg_4incl_2incl_mu",  "RemBkg_4incl_2incl_mu",  "RemBkgNorm_4incl_2incl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_mu");
	   
	 }
	 else if(FitType == "WSum"){
	   
	   fProf.AddBackgroundParameter("Wjets_allCombined",   "Wjets_allCombined",   "WjetsNorm_allCombined",   "WjetsUnc"+uncSuffix+"_allCombined");
	   
	   fProf.AddBackgroundParameter("QCD_4incl_0excl_mu",     "QCD_4incl_0excl_mu",     "QCDNorm_4incl_0excl_mu",     "QCDUnc"+uncSuffix+"_4incl_0excl_mu");
	   fProf.AddBackgroundParameter("QCD_4incl_1excl_mu",     "QCD_4incl_1excl_mu",     "QCDNorm_4incl_1excl_mu",     "QCDUnc"+uncSuffix+"_4incl_1excl_mu");
	   fProf.AddBackgroundParameter("QCD_4incl_2incl_mu",     "QCD_4incl_2incl_mu",     "QCDNorm_4incl_2incl_mu",     "QCDUnc"+uncSuffix+"_4incl_2incl_mu");
	   fProf.AddBackgroundParameter("RemBkg_4incl_0excl_mu",  "RemBkg_4incl_0excl_mu",  "RemBkgNorm_4incl_0excl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_0excl_mu");
	   fProf.AddBackgroundParameter("RemBkg_4incl_1excl_mu",  "RemBkg_4incl_1excl_mu",  "RemBkgNorm_4incl_1excl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_mu");
	   fProf.AddBackgroundParameter("RemBkg_4incl_2incl_mu",  "RemBkg_4incl_2incl_mu",  "RemBkgNorm_4incl_2incl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_mu");
	   
	 }
	 else if(FitType == "BkgSum"){
	   
           fProf.AddBackgroundParameter("Wjets_allCombined",  "Wjets_allCombined",  "WjetsNorm_allCombined",   "WjetsUnc"+uncSuffix+"_allCombined");
           fProf.AddBackgroundParameter("QCD_allCombined",    "QCD_allCombined",    "QCDNorm_allCombined",     "QCDUnc"+uncSuffix+"_allCombined");
           fProf.AddBackgroundParameter("RemBkg_allCombined", "RemBkg_allCombined", "RemBkgNorm_allCombined",  "RemBkgUnc"+uncSuffix+"_allCombined");
	   
         }
	 else if(FitType == "RemBkgSum"){

           fProf.AddBackgroundParameter("Wjets_4incl_0excl_mu",   "Wjets_4incl_0excl_mu",   "WjetsNorm_4incl_0excl_mu",   "WjetsUnc"+uncSuffix+"_4incl_0excl_mu");
           fProf.AddBackgroundParameter("Wjets_4incl_1excl_mu",   "Wjets_4incl_1excl_mu",   "WjetsNorm_4incl_1excl_mu",   "WjetsUnc"+uncSuffix+"_4incl_1excl_mu");
           fProf.AddBackgroundParameter("Wjets_4incl_2incl_mu",   "Wjets_4incl_2incl_mu",   "WjetsNorm_4incl_2incl_mu",   "WjetsUnc"+uncSuffix+"_4incl_2incl_mu");
           fProf.AddBackgroundParameter("QCD_4incl_0excl_mu",     "QCD_4incl_0excl_mu",     "QCDNorm_4incl_0excl_mu",     "QCDUnc"+uncSuffix+"_4incl_0excl_mu");
           fProf.AddBackgroundParameter("QCD_4incl_1excl_mu",     "QCD_4incl_1excl_mu",     "QCDNorm_4incl_1excl_mu",     "QCDUnc"+uncSuffix+"_4incl_1excl_mu");
           fProf.AddBackgroundParameter("QCD_4incl_2incl_mu",     "QCD_4incl_2incl_mu",     "QCDNorm_4incl_2incl_mu",     "QCDUnc"+uncSuffix+"_4incl_2incl_mu");
         
	   fProf.AddBackgroundParameter("RemBkg_allCombined",     "RemBkg_allCombined",     "RemBkgNorm_allCombined",     "RemBkgUnc"+uncSuffix+"_allCombined");
	   
         }

	 

       } // if mu_Btag

       HistoNames.push_back("CosTheta_F0_allCombined");
       HistoNames.push_back("CosTheta_FL_allCombined");
       if(FittingMode == "3D")
         HistoNames.push_back("CosTheta_FR_allCombined");

       if(InputChannel == "el_BTag"){

	 if(FitType == "All"){
	   //HistoNames.push_back("Wjets_4incl_0excl_el");
	   if(wjetsMode=="1W"){
	     HistoNames.push_back("Wjets_4incl_1excl_el");
	     HistoNames.push_back("Wjets_4incl_2incl_el"); 
	   }
	   if(wjetsMode=="3W"){
	     HistoNames.push_back("WLight_4incl_1excl_el");
	     HistoNames.push_back("WLight_4incl_2incl_el"); 
	     HistoNames.push_back("Wc_4incl_1excl_el");
	     HistoNames.push_back("Wc_4incl_2incl_el"); 
	     HistoNames.push_back("Wbbcc_4incl_1excl_el");
	     HistoNames.push_back("Wbbcc_4incl_2incl_el"); 
	   }
	   //HistoNames.push_back("QCD_4incl_0excl_el");
	   HistoNames.push_back("QCD_4incl_1excl_el");
	   HistoNames.push_back("QCD_4incl_2incl_el");
	   //HistoNames.push_back("RemBkg_4incl_0excl_el");
	   HistoNames.push_back("RemBkg_4incl_1excl_el");
	   HistoNames.push_back("RemBkg_4incl_2incl_el");
	   
	 }
	 else if(FitType == "WSum"){
	   
	   HistoNames.push_back("Wjets_allCombined");
	   
	   HistoNames.push_back("QCD_4incl_0excl_el");
	   HistoNames.push_back("QCD_4incl_1excl_el");
	   HistoNames.push_back("QCD_4incl_2incl_el");
	   HistoNames.push_back("RemBkg_4incl_0excl_el");
	   HistoNames.push_back("RemBkg_4incl_1excl_el");
	   HistoNames.push_back("RemBkg_4incl_2incl_el");
	 
	 }
	 else if(FitType == "BkgSum"){
	 
	   HistoNames.push_back("Wjets_allCombined");
	   HistoNames.push_back("QCD_allCombined");
	   HistoNames.push_back("RemBkg_allCombined");

	 }
	 else if(FitType == "RemBkgSum"){

           HistoNames.push_back("Wjets_4incl_0excl_el");
           HistoNames.push_back("Wjets_4incl_1excl_el");
           HistoNames.push_back("Wjets_4incl_2incl_el");

           HistoNames.push_back("QCD_4incl_0excl_el");
           HistoNames.push_back("QCD_4incl_1excl_el");
           HistoNames.push_back("QCD_4incl_2incl_el");

           HistoNames.push_back("RemBkg_allCombined");

         }

       }
       else{ // if mu_BTag
	 
	 if(FitType == "All"){
	   
           //HistoNames.push_back("Wjets_4incl_0excl_mu");
           if(wjetsMode=="1W"){
	     HistoNames.push_back("Wjets_4incl_1excl_mu");
	     HistoNames.push_back("Wjets_4incl_2incl_mu");
	   }
           else if(wjetsMode=="3W"){
	     HistoNames.push_back("WLight_4incl_1excl_mu");
	     HistoNames.push_back("WLight_4incl_2incl_mu");
	     HistoNames.push_back("Wc_4incl_1excl_mu");
	     HistoNames.push_back("Wc_4incl_2incl_mu");
	     HistoNames.push_back("Wbbcc_4incl_1excl_mu");
	     HistoNames.push_back("Wbbcc_4incl_2incl_mu");
	   }           
	   //HistoNames.push_back("QCD_4incl_0excl_mu");
           HistoNames.push_back("QCD_4incl_1excl_mu");
           HistoNames.push_back("QCD_4incl_2incl_mu");
           //HistoNames.push_back("RemBkg_4incl_0excl_mu");
           HistoNames.push_back("RemBkg_4incl_1excl_mu");
           HistoNames.push_back("RemBkg_4incl_2incl_mu");

         }
         else if(FitType == "WSum"){

           HistoNames.push_back("Wjets_allCombined");

           HistoNames.push_back("QCD_4incl_0excl_mu");
           HistoNames.push_back("QCD_4incl_1excl_mu");
           HistoNames.push_back("QCD_4incl_2incl_mu");
           HistoNames.push_back("RemBkg_4incl_0excl_mu");
           HistoNames.push_back("RemBkg_4incl_1excl_mu");
           HistoNames.push_back("RemBkg_4incl_2incl_mu");

         }
         else if(FitType == "BkgSum"){

           HistoNames.push_back("Wjets_allCombined");
           HistoNames.push_back("QCD_allCombined");
           HistoNames.push_back("RemBkg_allCombined");

         }
	 else if(FitType == "RemBkgSum"){

	   HistoNames.push_back("Wjets_4incl_0excl_mu");
           HistoNames.push_back("Wjets_4incl_1excl_mu");
           HistoNames.push_back("Wjets_4incl_2incl_mu");

           HistoNames.push_back("QCD_4incl_0excl_mu");
           HistoNames.push_back("QCD_4incl_1excl_mu");
           HistoNames.push_back("QCD_4incl_2incl_mu");
         
	         HistoNames.push_back("RemBkg_allCombined");

         }
       }

       VarNames.push_back("F0eff_allCombined");
       VarNames.push_back("FLeff_allCombined");
       if(FittingMode == "3D")
         VarNames.push_back("FReff_allCombined");
       
       if(InputChannel == "el_BTag"){
	 
	 if(FitType == "All"){

	   //VarNames.push_back("WjetsNorm_4incl_0excl_el");
	   if(wjetsMode=="1W"){
	     VarNames.push_back("WjetsNorm_4incl_1excl_el");
	     VarNames.push_back("WjetsNorm_4incl_2incl_el"); 
	   }
	   else if(wjetsMode=="3W"){
	     VarNames.push_back("WLightNorm_4incl_1excl_el");
	     VarNames.push_back("WLightNorm_4incl_2incl_el"); 
	     VarNames.push_back("WcNorm_4incl_1excl_el");
	     VarNames.push_back("WcNorm_4incl_2incl_el"); 
	     VarNames.push_back("WbbccNorm_4incl_1excl_el");
	     VarNames.push_back("WbbccNorm_4incl_2incl_el"); 
	   }	
	   //VarNames.push_back("QCDNorm_4incl_0excl_el");
	   VarNames.push_back("QCDNorm_4incl_1excl_el");
	   VarNames.push_back("QCDNorm_4incl_2incl_el");
	   //VarNames.push_back("RemBkgNorm_4incl_0excl_el");
	   VarNames.push_back("RemBkgNorm_4incl_1excl_el");
	   VarNames.push_back("RemBkgNorm_4incl_2incl_el");

	 }
	 else if(FitType == "WSum"){

	   VarNames.push_back("WjetsNorm_allCombined");
	   
	   VarNames.push_back("QCDNorm_4incl_0excl_el");
	   VarNames.push_back("QCDNorm_4incl_1excl_el");
	   VarNames.push_back("QCDNorm_4incl_2incl_el");
	   VarNames.push_back("RemBkgNorm_4incl_0excl_el");
	   VarNames.push_back("RemBkgNorm_4incl_1excl_el");
	   VarNames.push_back("RemBkgNorm_4incl_2incl_el");
	 
	 }
	 else if(FitType == "BkgSum"){
	 
	   VarNames.push_back("WjetsNorm_allCombined");
	   VarNames.push_back("QCDNorm_allCombined");
	   VarNames.push_back("RemBkgNorm_allCombined");

	 }
	 else if(FitType == "RemBkgSum"){

           VarNames.push_back("WjetsNorm_4incl_0excl_el");
           VarNames.push_back("WjetsNorm_4incl_1excl_el");
           VarNames.push_back("WjetsNorm_4incl_2incl_el");
           VarNames.push_back("QCDNorm_4incl_0excl_el");
           VarNames.push_back("QCDNorm_4incl_1excl_el");
           VarNames.push_back("QCDNorm_4incl_2incl_el");
          
	   VarNames.push_back("RemBkgNorm_allCombined");

         }


       }
       else{ // if mu_BTag
	 
	 if(FitType == "All"){

           //VarNames.push_back("WjetsNorm_4incl_0excl_mu");
           if(wjetsMode=="1W"){
	     VarNames.push_back("WjetsNorm_4incl_1excl_mu");
	     VarNames.push_back("WjetsNorm_4incl_2incl_mu");
	   }
           if(wjetsMode=="3W"){
	     VarNames.push_back("WLightNorm_4incl_1excl_mu");
	     VarNames.push_back("WLightNorm_4incl_2incl_mu");
	     VarNames.push_back("WcNorm_4incl_1excl_mu");
	     VarNames.push_back("WcNorm_4incl_2incl_mu");
	     VarNames.push_back("WbbccNorm_4incl_1excl_mu");
	     VarNames.push_back("WbbccNorm_4incl_2incl_mu");
	   } 
	   //VarNames.push_back("QCDNorm_4incl_0excl_mu");
           VarNames.push_back("QCDNorm_4incl_1excl_mu");
           VarNames.push_back("QCDNorm_4incl_2incl_mu");
           //VarNames.push_back("RemBkgNorm_4incl_0excl_mu");
           VarNames.push_back("RemBkgNorm_4incl_1excl_mu");
           VarNames.push_back("RemBkgNorm_4incl_2incl_mu");

         }
         else if(FitType == "WSum"){

           VarNames.push_back("WjetsNorm_allCombined");

           VarNames.push_back("QCDNorm_4incl_0excl_mu");
           VarNames.push_back("QCDNorm_4incl_1excl_mu");
           VarNames.push_back("QCDNorm_4incl_2incl_mu");
           VarNames.push_back("RemBkgNorm_4incl_0excl_mu");
           VarNames.push_back("RemBkgNorm_4incl_1excl_mu");
           VarNames.push_back("RemBkgNorm_4incl_2incl_mu");

         }
         else if(FitType == "BkgSum"){

           VarNames.push_back("WjetsNorm_allCombined");
           VarNames.push_back("QCDNorm_allCombined");
           VarNames.push_back("RemBkgNorm_allCombined");

         }
	 else if(FitType == "RemBkgSum"){

           VarNames.push_back("WjetsNorm_4incl_0excl_mu");
           VarNames.push_back("WjetsNorm_4incl_1excl_mu");
           VarNames.push_back("WjetsNorm_4incl_2incl_mu");
           VarNames.push_back("QCDNorm_4incl_0excl_mu");
           VarNames.push_back("QCDNorm_4incl_1excl_mu");
           VarNames.push_back("QCDNorm_4incl_2incl_mu");

           VarNames.push_back("RemBkgNorm_allCombined");

         }
	 
	 
       }
       
     }
     //--------------------------------------------------------------------------------
     else if (InputChannel == "el_mu_bTag"){ // 4-channels

       fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "1 b-tags");
       fProf.SetSubplotLabels("e+jets",   "#geq 4 jets", "#geq 2 b-tags");
       fProf.SetSubplotLabels("#mu+jets", "#geq 4 jets", "1 b-tags");
       fProf.SetSubplotLabels("#mu+jets", "#geq 4 jets", "#geq 2 b-tags");

       // add signal parameter with param name and hist name
       fProf.AddSignalParameter("N0",  "CosTheta_F0_allCombined", "F0eff_allCombined");
       fProf.AddSignalParameter("NL",  "CosTheta_FL_allCombined", "FLeff_allCombined");
       if(FittingMode == "3D")
         fProf.AddSignalParameter("NR",  "CosTheta_FR_allCombined", "FReff_allCombined");
       
       // add bkg parameter with param name and hist name
       if(wjetsMode=="1W"){
	 
	       // fProf.AddBackgroundParameter("Wjets_4incl_1excl_el",   "Wjets_4incl_1excl_el",   "WjetsNorm_4incl_1excl_el",   "WjetsUnc"+uncSuffix+"_4incl_1excl_el");
	       // fProf.AddBackgroundParameter("Wjets_4incl_2incl_el",   "Wjets_4incl_2incl_el",   "WjetsNorm_4incl_2incl_el",   "WjetsUnc"+uncSuffix+"_4incl_2incl_el");
	       // fProf.AddBackgroundParameter("Wjets_4incl_1excl_mu",   "Wjets_4incl_1excl_mu",   "WjetsNorm_4incl_1excl_mu",   "WjetsUnc"+uncSuffix+"_4incl_1excl_mu");
	       // fProf.AddBackgroundParameter("Wjets_4incl_2incl_mu",   "Wjets_4incl_2incl_mu",   "WjetsNorm_4incl_2incl_mu",   "WjetsUnc"+uncSuffix+"_4incl_2incl_mu");
        fProf.AddBackgroundParameter("Wjets_4incl",   "Wjets_4incl",   "WjetsNorm_4incl",   "WjetsUnc"+uncSuffix+"_4incl");
       }
       else if(wjetsMode=="3W"){
	       // fProf.AddBackgroundParameter("WLight_4incl_1excl_el",   "WLight_4incl_1excl_el",   "WLightNorm_4incl_1excl_el",   "WLightUnc"+uncSuffix+"_4incl_1excl_el");
	       // fProf.AddBackgroundParameter("WLight_4incl_2incl_el",   "WLight_4incl_2incl_el",   "WLightNorm_4incl_2incl_el",   "WLightUnc"+uncSuffix+"_4incl_2incl_el");
	       // fProf.AddBackgroundParameter("WLight_4incl_1excl_mu",   "WLight_4incl_1excl_mu",   "WLightNorm_4incl_1excl_mu",   "WLightUnc"+uncSuffix+"_4incl_1excl_mu");
	       // fProf.AddBackgroundParameter("WLight_4incl_2incl_mu",   "WLight_4incl_2incl_mu",   "WLightNorm_4incl_2incl_mu",   "WLightUnc"+uncSuffix+"_4incl_2incl_mu");
	       // fProf.AddBackgroundParameter("Wc_4incl_1excl_el",   "Wc_4incl_1excl_el",   "WcNorm_4incl_1excl_el",   "WcUnc"+uncSuffix+"_4incl_1excl_el");
	       // fProf.AddBackgroundParameter("Wc_4incl_2incl_el",   "Wc_4incl_2incl_el",   "WcNorm_4incl_2incl_el",   "WcUnc"+uncSuffix+"_4incl_2incl_el");
	       // fProf.AddBackgroundParameter("Wc_4incl_1excl_mu",   "Wc_4incl_1excl_mu",   "WcNorm_4incl_1excl_mu",   "WcUnc"+uncSuffix+"_4incl_1excl_mu");
	       // fProf.AddBackgroundParameter("Wc_4incl_2incl_mu",   "Wc_4incl_2incl_mu",   "WcNorm_4incl_2incl_mu",   "WcUnc"+uncSuffix+"_4incl_2incl_mu");
	       // fProf.AddBackgroundParameter("Wbbcc_4incl_1excl_el",   "Wbbcc_4incl_1excl_el",   "WbbccNorm_4incl_1excl_el",   "WbbccUnc"+uncSuffix+"_4incl_1excl_el");
	       // fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_el",   "Wbbcc_4incl_2incl_el",   "WbbccNorm_4incl_2incl_el",   "WbbccUnc"+uncSuffix+"_4incl_2incl_el");
	       // fProf.AddBackgroundParameter("Wbbcc_4incl_1excl_mu",   "Wbbcc_4incl_1excl_mu",   "WbbccNorm_4incl_1excl_mu",   "WbbccUnc"+uncSuffix+"_4incl_1excl_mu");
	       // fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_mu",   "Wbbcc_4incl_2incl_mu",   "WbbccNorm_4incl_2incl_mu",   "WbbccUnc"+uncSuffix+"_4incl_2incl_mu");

         // fProf.AddBackgroundParameter("WLight_4incl_el",   "WLight_4incl_el",   "WLightNorm_4incl_el",   "WLightUnc"+uncSuffix+"_4incl_el");
         // fProf.AddBackgroundParameter("WLight_4incl_mu",   "WLight_4incl_mu",   "WLightNorm_4incl_mu",   "WLightUnc"+uncSuffix+"_4incl_mu");
        fProf.AddBackgroundParameter("WLight_4incl",   "WLight_4incl",   "WLightNorm_4incl",   "WLightUnc"+uncSuffix+"_4incl");
         // fProf.AddBackgroundParameter("Wc_4incl_el",   "Wc_4incl_el",   "WcNorm_4incl_el",   "WcUnc"+uncSuffix+"_4incl_el");
         // fProf.AddBackgroundParameter("Wc_4incl_mu",   "Wc_4incl_mu",   "WcNorm_4incl_mu",   "WcUnc"+uncSuffix+"_4incl_mu");
        fProf.AddBackgroundParameter("Wc_4incl",   "Wc_4incl",   "WcNorm_4incl",   "WcUnc"+uncSuffix+"_4incl");
         // fProf.AddBackgroundParameter("Wbbcc_4incl_el",   "Wbbcc_4incl_el",   "WbbccNorm_4incl_el",   "WbbccUnc"+uncSuffix+"_4incl_el");
         // fProf.AddBackgroundParameter("Wbbcc_4incl_mu",   "Wbbcc_4incl_mu",   "WbbccNorm_4incl_mu",   "WbbccUnc"+uncSuffix+"_4incl_mu");
        fProf.AddBackgroundParameter("Wbbcc_4incl",   "Wbbcc_4incl",   "WbbccNorm_4incl",   "WbbccUnc"+uncSuffix+"_4incl");
         
       }
       
       fProf.AddBackgroundParameter("QCD_4incl_1excl_el",     "QCD_4incl_1excl_el",     "QCDNorm_4incl_1excl_el",     "QCDUnc"+uncSuffix+"_4incl_1excl_el");
       fProf.AddBackgroundParameter("QCD_4incl_2incl_el",     "QCD_4incl_2incl_el",     "QCDNorm_4incl_2incl_el",     "QCDUnc"+uncSuffix+"_4incl_2incl_el");
       fProf.AddBackgroundParameter("QCD_4incl_1excl_mu",     "QCD_4incl_1excl_mu",     "QCDNorm_4incl_1excl_mu",     "QCDUnc"+uncSuffix+"_4incl_1excl_mu");
       fProf.AddBackgroundParameter("QCD_4incl_2incl_mu",     "QCD_4incl_2incl_mu",     "QCDNorm_4incl_2incl_mu",     "QCDUnc"+uncSuffix+"_4incl_2incl_mu");
       
       // fProf.AddBackgroundParameter("RemBkg_4incl_1excl_el",  "RemBkg_4incl_1excl_el",  "RemBkgNorm_4incl_1excl_el",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_el");
       // fProf.AddBackgroundParameter("RemBkg_4incl_1excl_mu",  "RemBkg_4incl_1excl_mu",  "RemBkgNorm_4incl_1excl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_1excl_mu");
       // fProf.AddBackgroundParameter("RemBkg_4incl_2incl_el",  "RemBkg_4incl_2incl_el",  "RemBkgNorm_4incl_2incl_el",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_el");
       // fProf.AddBackgroundParameter("RemBkg_4incl_2incl_mu",  "RemBkg_4incl_2incl_mu",  "RemBkgNorm_4incl_2incl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_2incl_mu");
       
       // fProf.AddBackgroundParameter("RemBkg_4incl_el",  "RemBkg_4incl_el",  "RemBkgNorm_4incl_el",  "RemBkgUnc"+uncSuffix+"_4incl_el");
       // fProf.AddBackgroundParameter("RemBkg_4incl_mu",  "RemBkg_4incl_mu",  "RemBkgNorm_4incl_mu",  "RemBkgUnc"+uncSuffix+"_4incl_mu");
       fProf.AddBackgroundParameter("RemBkg_4incl",  "RemBkg_4incl",  "RemBkgNorm_4incl",  "RemBkgUnc"+uncSuffix+"_4incl");


       HistoNames.push_back("CosTheta_F0_allCombined");
       HistoNames.push_back("CosTheta_FL_allCombined");
       if(FittingMode == "3D")
         HistoNames.push_back("CosTheta_FR_allCombined");

       
       if(wjetsMode=="1W"){
	       // HistoNames.push_back("Wjets_4incl_1excl_el");
	       // HistoNames.push_back("Wjets_4incl_2incl_el");
	       // HistoNames.push_back("Wjets_4incl_1excl_mu");
	       // HistoNames.push_back("Wjets_4incl_2incl_mu");
        HistoNames.push_back("Wjets_4incl");
       }
       else if(wjetsMode=="3W"){
	         //HistoNames.push_back("WLight_4incl_1excl_el");
	         //HistoNames.push_back("WLight_4incl_2incl_el");
	         //HistoNames.push_back("WLight_4incl_1excl_mu");
	         //HistoNames.push_back("WLight_4incl_2incl_mu");
	         //HistoNames.push_back("Wc_4incl_1excl_el");
	         //HistoNames.push_back("Wc_4incl_2incl_el");
	         //HistoNames.push_back("Wc_4incl_1excl_mu");
	         //HistoNames.push_back("Wc_4incl_2incl_mu");	 
	         //HistoNames.push_back("Wbbcc_4incl_1excl_el");
	         //HistoNames.push_back("Wbbcc_4incl_2incl_el");
	         //HistoNames.push_back("Wbbcc_4incl_1excl_mu");
	         //HistoNames.push_back("Wbbcc_4incl_2incl_mu");

        // HistoNames.push_back("WLight_4incl_el");
        // HistoNames.push_back("WLight_4incl_mu");
        HistoNames.push_back("WLight_4incl");
        // HistoNames.push_back("Wc_4incl_el");
        // HistoNames.push_back("Wc_4incl_mu");
        HistoNames.push_back("Wc_4incl");
        // HistoNames.push_back("Wbbcc_4incl_el");
        // HistoNames.push_back("Wbbcc_4incl_mu");
        HistoNames.push_back("Wbbcc_4incl");
       }    
       
       HistoNames.push_back("QCD_4incl_1excl_el");
       HistoNames.push_back("QCD_4incl_2incl_el");
       HistoNames.push_back("QCD_4incl_1excl_mu");
       HistoNames.push_back("QCD_4incl_2incl_mu");
       
       // HistoNames.push_back("RemBkg_4incl_1excl_el");
       // HistoNames.push_back("RemBkg_4incl_2incl_el");
       // HistoNames.push_back("RemBkg_4incl_1excl_mu");
       // HistoNames.push_back("RemBkg_4incl_2incl_mu");
       
       // HistoNames.push_back("RemBkg_4incl_el");
       // HistoNames.push_back("RemBkg_4incl_mu");
       HistoNames.push_back("RemBkg_4incl");
       
       
       VarNames.push_back("F0eff_allCombined");
       VarNames.push_back("FLeff_allCombined");
       if(FittingMode == "3D")
         VarNames.push_back("FReff_allCombined");
       
       
       if(wjetsMode=="1W"){
	     // VarNames.push_back("WjetsNorm_4incl_1excl_el");
	     // VarNames.push_back("WjetsNorm_4incl_2incl_el");
	     // VarNames.push_back("WjetsNorm_4incl_1excl_mu");
	     // VarNames.push_back("WjetsNorm_4incl_2incl_mu");
        VarNames.push_back("WjetsNorm_4incl");
       }
       else if(wjetsMode=="3W"){
	     // VarNames.push_back("WLightNorm_4incl_1excl_el");
	     // VarNames.push_back("WLightNorm_4incl_2incl_el");
	     // VarNames.push_back("WLightNorm_4incl_1excl_mu");
	     // VarNames.push_back("WLightNorm_4incl_2incl_mu");
	     // VarNames.push_back("WcNorm_4incl_1excl_el");
	     // VarNames.push_back("WcNorm_4incl_2incl_el");
	     // VarNames.push_back("WcNorm_4incl_1excl_mu");
	     // VarNames.push_back("WcNorm_4incl_2incl_mu");
	     // VarNames.push_back("WbbccNorm_4incl_1excl_el");
	     // VarNames.push_back("WbbccNorm_4incl_2incl_el");
	     // VarNames.push_back("WbbccNorm_4incl_1excl_mu");
	     // VarNames.push_back("WbbccNorm_4incl_2incl_mu");

         // VarNames.push_back("WLightNorm_4incl_el");
         // VarNames.push_back("WLightNorm_4incl_mu");
        VarNames.push_back("WLightNorm_4incl");
         // VarNames.push_back("WcNorm_4incl_el");
         // VarNames.push_back("WcNorm_4incl_mu");
        VarNames.push_back("WcNorm_4incl");
         // VarNames.push_back("WbbccNorm_4incl_el");
         // VarNames.push_back("WbbccNorm_4incl_mu");
        VarNames.push_back("WbbccNorm_4incl");
       }
       VarNames.push_back("QCDNorm_4incl_1excl_el");
       VarNames.push_back("QCDNorm_4incl_2incl_el");
       VarNames.push_back("QCDNorm_4incl_1excl_mu");
       VarNames.push_back("QCDNorm_4incl_2incl_mu");
       // VarNames.push_back("RemBkgNorm_4incl_1excl_el");
       // VarNames.push_back("RemBkgNorm_4incl_2incl_el");
       // VarNames.push_back("RemBkgNorm_4incl_1excl_mu");
       // VarNames.push_back("RemBkgNorm_4incl_2incl_mu");
       
       // VarNames.push_back("RemBkgNorm_4incl_el");
       // VarNames.push_back("RemBkgNorm_4incl_mu");
       VarNames.push_back("RemBkgNorm_4incl");
	 
     }
     //-------------------------------------------------------------------------------
      // Lep-had combination
     else if (InputChannel == "el_mu_lephad"){ // 4-channels

       if(tagMode=="2incl"){

       fProf.SetSubplotLabels("e+jets(lep)",   "#geq 4 jets", "#geq 2 b-tags");
       fProf.SetSubplotLabels("#mu+jets(lep)", "#geq 4 jets", "#geq 2 b-tags");
       fProf.SetSubplotLabels("e+jets(had)",   "#geq 4 jets", "#geq 2 b-tags");
       fProf.SetSubplotLabels("#mu+jets(had)", "#geq 4 jets", "#geq 2 b-tags");
       

       // add signal parameter with param name and hist name
       fProf.AddSignalParameter("N0",  "CosTheta_F0_allCombined", "F0eff_allCombined");
       fProf.AddSignalParameter("NL",  "CosTheta_FL_allCombined", "FLeff_allCombined");
       if(FittingMode == "3D")
         fProf.AddSignalParameter("NR",  "CosTheta_FR_allCombined", "FReff_allCombined");
       
    // add bkg parameter with param name and hist name   
	 if(wjetsMode=="1W"){
	   // fProf.AddBackgroundParameter("Wjets_4incl_2incl_el",   "Wjets_4incl_2incl_el",   "WjetsNorm_4incl_2incl_el",   "WjetsUnc"+uncSuffix+"_4incl_2incl_el" );
	   // fProf.AddBackgroundParameter("Wjets_4incl_2incl_mu",   "Wjets_4incl_2incl_mu",   "WjetsNorm_4incl_2incl_mu",   "WjetsUnc"+uncSuffix+"_4incl_2incl_mu" );
    fProf.AddBackgroundParameter("Wjets_4incl_2incl",   "Wjets_4incl_2incl",   "WjetsNorm_4incl_2incl",   "WjetsUnc"+uncSuffix+"_4incl_2incl" );
	 }
	 else if(wjetsMode=="3W"){
	   // fProf.AddBackgroundParameter("WLight_4incl_2incl_el",   "WLight_4incl_2incl_el",   "WLightNorm_4incl_2incl_el",   "WLightUnc"+uncSuffix+"_4incl_2incl_el" );
	   // fProf.AddBackgroundParameter("WLight_4incl_2incl_mu",   "WLight_4incl_2incl_mu",   "WLightNorm_4incl_2incl_mu",   "WLightUnc"+uncSuffix+"_4incl_2incl_mu" );
    fProf.AddBackgroundParameter("WLight_4incl_2incl",   "WLight_4incl_2incl",   "WLightNorm_4incl_2incl",   "WLightUnc"+uncSuffix+"_4incl_2incl" );
	   // fProf.AddBackgroundParameter("Wc_4incl_2incl_el",   "Wc_4incl_2incl_el",   "WcNorm_4incl_2incl_el",   "WcUnc"+uncSuffix+"_4incl_2incl_el" );
	   // fProf.AddBackgroundParameter("Wc_4incl_2incl_mu",   "Wc_4incl_2incl_mu",   "WcNorm_4incl_2incl_mu",   "WcUnc"+uncSuffix+"_4incl_2incl_mu" );
    fProf.AddBackgroundParameter("Wc_4incl_2incl",   "Wc_4incl_2incl",   "WcNorm_4incl_2incl",   "WcUnc"+uncSuffix+"_4incl_2incl" );
	   // fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_el",   "Wbbcc_4incl_2incl_el",   "WbbccNorm_4incl_2incl_el",   "WbbccUnc"+uncSuffix+"_4incl_2incl_el" );
	   // fProf.AddBackgroundParameter("Wbbcc_4incl_2incl_mu",   "Wbbcc_4incl_2incl_mu",   "WbbccNorm_4incl_2incl_mu",   "WbbccUnc"+uncSuffix+"_4incl_2incl_mu" );
    fProf.AddBackgroundParameter("Wbbcc_4incl_2incl",   "Wbbcc_4incl_2incl",   "WbbccNorm_4incl_2incl",   "WbbccUnc"+uncSuffix+"_4incl_2incl" );

	 }
	 fProf.AddBackgroundParameter("QCD_4incl_2incl_el",   "QCD_4incl_2incl_el",   "QCDNorm_4incl_2incl_el",   "QCDUnc"+uncSuffix+"_4incl_2incl_el");
	 fProf.AddBackgroundParameter("QCD_4incl_2incl_mu",   "QCD_4incl_2incl_mu",   "QCDNorm_4incl_2incl_mu",   "QCDUnc"+uncSuffix+"_4incl_2incl_mu");
	 
	 // fProf.AddBackgroundParameter("RemBkg_4incl_2incl_el",   "RemBkg_4incl_2incl_el",   "RemBkgNorm_4incl_2incl_el",   "RemBkgUnc"+uncSuffix+"_4incl_2incl_el");
	 // fProf.AddBackgroundParameter("RemBkg_4incl_2incl_mu",   "RemBkg_4incl_2incl_mu",   "RemBkgNorm_4incl_2incl_mu",   "RemBkgUnc"+uncSuffix+"_4incl_2incl_mu");
   fProf.AddBackgroundParameter("RemBkg_4incl_2incl",   "RemBkg_4incl_2incl",   "RemBkgNorm_4incl_2incl",   "RemBkgUnc"+uncSuffix+"_4incl_2incl");
	 
	 HistoNames.push_back("CosTheta_F0_allCombined");
	 HistoNames.push_back("CosTheta_FL_allCombined");
	 if(FittingMode == "3D")
	   HistoNames.push_back("CosTheta_FR_allCombined");
	 
	 if(wjetsMode=="1W"){
	   // HistoNames.push_back("Wjets_4incl_2incl_el");
	   // HistoNames.push_back("Wjets_4incl_2incl_mu");
    HistoNames.push_back("Wjets_4incl_2incl");
	 }
	 else if(wjetsMode=="3W"){
	   // HistoNames.push_back("WLight_4incl_2incl_el");
	   // HistoNames.push_back("WLight_4incl_2incl_mu");
    HistoNames.push_back("WLight_4incl_2incl");
	   // HistoNames.push_back("Wc_4incl_2incl_el");
	   // HistoNames.push_back("Wc_4incl_2incl_mu");
    HistoNames.push_back("Wc_4incl_2incl");
	   // HistoNames.push_back("Wbbcc_4incl_2incl_el");
	   // HistoNames.push_back("Wbbcc_4incl_2incl_mu");
    HistoNames.push_back("Wbbcc_4incl_2incl");
	 }
	 HistoNames.push_back("QCD_4incl_2incl_el");
	 HistoNames.push_back("QCD_4incl_2incl_mu");
	
	 // HistoNames.push_back("RemBkg_4incl_2incl_el");
	 // HistoNames.push_back("RemBkg_4incl_2incl_mu");
   HistoNames.push_back("RemBkg_4incl_2incl");
	  
	 VarNames.push_back("F0eff_allCombined");
	 VarNames.push_back("FLeff_allCombined");
	 if(FittingMode == "3D")
         VarNames.push_back("FReff_allCombined");
       
	 if(wjetsMode=="1W"){
	   // VarNames.push_back("WjetsNorm_4incl_2incl_el");
	   // VarNames.push_back("WjetsNorm_4incl_2incl_mu");
    VarNames.push_back("WjetsNorm_4incl_2incl");
	 
	 }
	 else if(wjetsMode=="3W"){
	   // VarNames.push_back("WLightNorm_4incl_2incl_el");
	   // VarNames.push_back("WLightNorm_4incl_2incl_mu");
    VarNames.push_back("WLightNorm_4incl_2incl");
	   // VarNames.push_back("WcNorm_4incl_2incl_el");
	   // VarNames.push_back("WcNorm_4incl_2incl_mu");
    VarNames.push_back("WcNorm_4incl_2incl");
	   // VarNames.push_back("WbbccNorm_4incl_2incl_el");
	   // VarNames.push_back("WbbccNorm_4incl_2incl_mu");
    VarNames.push_back("WbbccNorm_4incl_2incl");
	 }
	 VarNames.push_back("QCDNorm_4incl_2incl_el");
	 VarNames.push_back("QCDNorm_4incl_2incl_mu");
	 
	 
	 // VarNames.push_back("RemBkgNorm_4incl_2incl_el");
	 // VarNames.push_back("RemBkgNorm_4incl_2incl_mu");
   VarNames.push_back("RemBkgNorm_4incl_2incl");
	 
       
  } // 2incl only (no need to have 1 excl combination)
       
} // end el_mu_lephad channel

//--------------------------------------------------------------------------

else if (InputChannel == "el_mu_lephad_bTag"){ // 8-channels

  fProf.SetSubplotLabels("e+jets(lep)",   "#geq 4 jets", "#eq 1 b-tags");
  fProf.SetSubplotLabels("e+jets(had)",   "#geq 4 jets", "#eq 1 b-tags");
  fProf.SetSubplotLabels("#mu+jets(lep)", "#geq 4 jets", "#eq 1 b-tags");
  fProf.SetSubplotLabels("#mu+jets(had)", "#geq 4 jets", "#eq 1 b-tags");
  fProf.SetSubplotLabels("e+jets(lep)",   "#geq 4 jets", "#geq 2 b-tags");
  fProf.SetSubplotLabels("e+jets(had)",   "#geq 4 jets", "#geq 2 b-tags");
  fProf.SetSubplotLabels("#mu+jets(lep)", "#geq 4 jets", "#geq 2 b-tags");
  fProf.SetSubplotLabels("#mu+jets(had)", "#geq 4 jets", "#geq 2 b-tags");


  // add signal parameter with param name and hist name
  fProf.AddSignalParameter("N0",  "CosTheta_F0_allCombined", "F0eff_allCombined");
  fProf.AddSignalParameter("NL",  "CosTheta_FL_allCombined", "FLeff_allCombined");
  if(FittingMode == "3D")
    fProf.AddSignalParameter("NR",  "CosTheta_FR_allCombined", "FReff_allCombined");
       
  if(FitType == "All"){
    // add bkg parameter with param name and hist name
    if(wjetsMode=="1W"){
      // fProf.AddBackgroundParameter("Wjets_4incl_1excl_el",   "Wjets_4incl_1excl_el",   "WjetsNorm_4incl_1excl_el",   "WjetsUnc"+uncSuffix+"_4incl_1excl_el" );
      // fProf.AddBackgroundParameter("Wjets_4incl_2incl_el",   "Wjets_4incl_2incl_el",   "WjetsNorm_4incl_2incl_el",   "WjetsUnc"+uncSuffix+"_4incl_2incl_el" );
      // fProf.AddBackgroundParameter("Wjets_4incl_1excl_mu",   "Wjets_4incl_1excl_mu",   "WjetsNorm_4incl_1excl_mu",   "WjetsUnc"+uncSuffix+"_4incl_1excl_mu" );
      // fProf.AddBackgroundParameter("Wjets_4incl_2incl_mu",   "Wjets_4incl_2incl_mu",   "WjetsNorm_4incl_2incl_mu",   "WjetsUnc"+uncSuffix+"_4incl_2incl_mu" );
      fProf.AddBackgroundParameter("Wjets_4incl",   "Wjets_4incl",   "WjetsNorm_4incl",   "WjetsUnc"+uncSuffix+"_4incl" );
    }
    else if(wjetsMode=="3W"){
      // fProf.AddBackgroundParameter("WLight_4incl_el",   "WLight_4incl_el",   "WLightNorm_4incl_el",   "WLightUnc"+uncSuffix+"_4incl_el" );
      // fProf.AddBackgroundParameter("WLight_4incl_mu",   "WLight_4incl_mu",   "WLightNorm_4incl_mu",   "WLightUnc"+uncSuffix+"_4incl_mu" );
      fProf.AddBackgroundParameter("WLight_4incl",   "WLight_4incl",   "WLightNorm_4incl",   "WLightUnc"+uncSuffix+"_4incl" );
      // fProf.AddBackgroundParameter("Wc_4incl_el",   "Wc_4incl_el",   "WcNorm_4incl_el",   "WcUnc"+uncSuffix+"_4incl_el" );
      // fProf.AddBackgroundParameter("Wc_4incl_mu",   "Wc_4incl_mu",   "WcNorm_4incl_mu",   "WcUnc"+uncSuffix+"_4incl_mu" );
      fProf.AddBackgroundParameter("Wc_4incl",   "Wc_4incl",   "WcNorm_4incl",   "WcUnc"+uncSuffix+"_4incl" );
      // fProf.AddBackgroundParameter("Wbbcc_4incl_el",   "Wbbcc_4incl_el",   "WbbccNorm_4incl_el",   "WbbccUnc"+uncSuffix+"_4incl_el" );
      // fProf.AddBackgroundParameter("Wbbcc_4incl_mu",   "Wbbcc_4incl_mu",   "WbbccNorm_4incl_mu",   "WbbccUnc"+uncSuffix+"_4incl_mu" );
      fProf.AddBackgroundParameter("Wbbcc_4incl",   "Wbbcc_4incl",   "WbbccNorm_4incl",   "WbbccUnc"+uncSuffix+"_4incl" );
    }
    fProf.AddBackgroundParameter("QCD_4incl_1excl_el",   "QCD_4incl_1excl_el",   "QCDNorm_4incl_1excl_el",   "QCDUnc"+uncSuffix+"_4incl_1excl_el");
    fProf.AddBackgroundParameter("QCD_4incl_1excl_mu",   "QCD_4incl_1excl_mu",   "QCDNorm_4incl_1excl_mu",   "QCDUnc"+uncSuffix+"_4incl_1excl_mu");
    fProf.AddBackgroundParameter("QCD_4incl_2incl_el",   "QCD_4incl_2incl_el",   "QCDNorm_4incl_2incl_el",   "QCDUnc"+uncSuffix+"_4incl_2incl_el");
    fProf.AddBackgroundParameter("QCD_4incl_2incl_mu",   "QCD_4incl_2incl_mu",   "QCDNorm_4incl_2incl_mu",   "QCDUnc"+uncSuffix+"_4incl_2incl_mu");
    

    // fProf.AddBackgroundParameter("RemBkg_4incl_el",   "RemBkg_4incl_el",   "RemBkgNorm_4incl_el",   "RemBkgUnc"+uncSuffix+"_4incl_el");
    // fProf.AddBackgroundParameter("RemBkg_4incl_mu",   "RemBkg_4incl_mu",   "RemBkgNorm_4incl_mu",   "RemBkgUnc"+uncSuffix+"_4incl_mu");
    fProf.AddBackgroundParameter("RemBkg_4incl",   "RemBkg_4incl",   "RemBkgNorm_4incl",   "RemBkgUnc"+uncSuffix+"_4incl");
  }
  
  if(FitType == "RemBkgSum"){
    // fProf.AddBackgroundParameter("RemBkg_4incl_el",   "RemBkg_4incl_el",   "RemBkgNorm_4incl_el",   "RemBkgUnc"+uncSuffix+"_4incl_el");
    // fProf.AddBackgroundParameter("RemBkg_4incl_mu",   "RemBkg_4incl_mu",   "RemBkgNorm_4incl_mu",   "RemBkgUnc"+uncSuffix+"_4incl_mu");
    fProf.AddBackgroundParameter("RemBkg_4incl",   "RemBkg_4incl",   "RemBkgNorm_4incl",   "RemBkgUnc"+uncSuffix+"_4incl");
  }

  HistoNames.push_back("CosTheta_F0_allCombined");
  HistoNames.push_back("CosTheta_FL_allCombined");
  if(FittingMode == "3D")
    HistoNames.push_back("CosTheta_FR_allCombined");

  if(FitType == "All"){
    if(wjetsMode=="1W"){
      HistoNames.push_back("Wjets_4incl_1excl_el");
      HistoNames.push_back("Wjets_4incl_1excl_mu");
      HistoNames.push_back("Wjets_4incl_2incl_el");
      HistoNames.push_back("Wjets_4incl_2incl_mu");
    }
    else if(wjetsMode=="3W"){
      // HistoNames.push_back("WLight_4incl_1excl_el");
      // HistoNames.push_back("WLight_4incl_1excl_mu");
      // HistoNames.push_back("Wc_4incl_1excl_el");
      // HistoNames.push_back("Wc_4incl_1excl_mu");
      // HistoNames.push_back("Wbbcc_4incl_1excl_el");
      // HistoNames.push_back("Wbbcc_4incl_1excl_mu");
      // HistoNames.push_back("WLight_4incl_2incl_el");
      // HistoNames.push_back("WLight_4incl_2incl_mu");
      // HistoNames.push_back("Wc_4incl_2incl_el");
      // HistoNames.push_back("Wc_4incl_2incl_mu");
      // HistoNames.push_back("Wbbcc_4incl_2incl_el");
      // HistoNames.push_back("Wbbcc_4incl_2incl_mu");
      
      // HistoNames.push_back("WLight_4incl_el");
      // HistoNames.push_back("WLight_4incl_mu");
      HistoNames.push_back("WLight_4incl");
      // HistoNames.push_back("Wc_4incl_el");
      // HistoNames.push_back("Wc_4incl_mu");
      HistoNames.push_back("Wc_4incl");
      // HistoNames.push_back("Wbbcc_4incl_el");
      // HistoNames.push_back("Wbbcc_4incl_mu");
      HistoNames.push_back("Wbbcc_4incl");
    }
    
    HistoNames.push_back("QCD_4incl_1excl_el");
    HistoNames.push_back("QCD_4incl_1excl_mu");
    HistoNames.push_back("QCD_4incl_2incl_el");
    HistoNames.push_back("QCD_4incl_2incl_mu");
    
    // HistoNames.push_back("RemBkg_4incl_1excl_el");
    // HistoNames.push_back("RemBkg_4incl_1excl_mu");
    // HistoNames.push_back("RemBkg_4incl_2incl_el");
    // HistoNames.push_back("RemBkg_4incl_2incl_mu");
    
    // HistoNames.push_back("RemBkg_4incl_el");
    // HistoNames.push_back("RemBkg_4incl_mu");
    HistoNames.push_back("RemBkg_4incl");
  }

  if(FitType == "RemBkgSum"){
    // HistoNames.push_back("RemBkg_4incl_el");
    // HistoNames.push_back("RemBkg_4incl_mu");
    HistoNames.push_back("RemBkg_4incl");
  }

  VarNames.push_back("F0eff_allCombined");
  VarNames.push_back("FLeff_allCombined");
  if(FittingMode == "3D")
    VarNames.push_back("FReff_allCombined");

  if(FitType == "All"){
    if(wjetsMode=="1W"){
      // VarNames.push_back("WjetsNorm_4incl_1excl_el");
      // VarNames.push_back("WjetsNorm_4incl_1excl_mu");
      // VarNames.push_back("WjetsNorm_4incl_2incl_el");
      // VarNames.push_back("WjetsNorm_4incl_2incl_mu");
      VarNames.push_back("WjetsNorm_4incl");
    }
    else if(wjetsMode=="3W"){
      // VarNames.push_back("WLightNorm_4incl_1excl_el");
      // VarNames.push_back("WLightNorm_4incl_1excl_mu");
      // VarNames.push_back("WcNorm_4incl_1excl_el");
      // VarNames.push_back("WcNorm_4incl_1excl_mu");
      // VarNames.push_back("WbbccNorm_4incl_1excl_el");
      // VarNames.push_back("WbbccNorm_4incl_1excl_mu");
      
      // VarNames.push_back("WLightNorm_4incl_2incl_el");
      // VarNames.push_back("WLightNorm_4incl_2incl_mu");
      // VarNames.push_back("WcNorm_4incl_2incl_el");
      // VarNames.push_back("WcNorm_4incl_2incl_mu");
      // VarNames.push_back("WbbccNorm_4incl_2incl_el");
      // VarNames.push_back("WbbccNorm_4incl_2incl_mu");
      
      // VarNames.push_back("WLightNorm_4incl_el");
      // VarNames.push_back("WLightNorm_4incl_mu");
      VarNames.push_back("WLightNorm_4incl");
      // VarNames.push_back("WcNorm_4incl_el");
      // VarNames.push_back("WcNorm_4incl_mu");
      VarNames.push_back("WcNorm_4incl");
      // VarNames.push_back("WbbccNorm_4incl_el");
      // VarNames.push_back("WbbccNorm_4incl_mu");
      VarNames.push_back("WbbccNorm_4incl");
    }
    VarNames.push_back("QCDNorm_4incl_1excl_el");
    VarNames.push_back("QCDNorm_4incl_1excl_mu");
    VarNames.push_back("QCDNorm_4incl_2incl_el");
    VarNames.push_back("QCDNorm_4incl_2incl_mu");
    
    // VarNames.push_back("RemBkgNorm_4incl_1excl_el");
    // VarNames.push_back("RemBkgNorm_4incl_1excl_mu");
    // VarNames.push_back("RemBkgNorm_4incl_2incl_el");
    // VarNames.push_back("RemBkgNorm_4incl_2incl_mu");
    // VarNames.push_back("RemBkgNorm_4incl_el");
    // VarNames.push_back("RemBkgNorm_4incl_mu");
    VarNames.push_back("RemBkgNorm_4incl");
  }

  if(FitType == "RemBkgSum"){
    // VarNames.push_back("RemBkgNorm_4incl_el");
    // VarNames.push_back("RemBkgNorm_4incl_mu");
    VarNames.push_back("RemBkgNorm_4incl");
  }
} // end el_mu_lephad_bTag channel

     
//==========================================================
     
     // only need to specify data hist name
     fProf.SetDataHist("Data");
     
     std::string fOutputFolder;
     std::string fOutputFolderCali;
     std::string fOutputFolderData;
     struct stat sb;
     
     
     
     if(Systematic == "Calibration"){
       if(fixedBkg){
	       fOutputFolderCali = "ExternalCalibrationOutput_"+InputChannel+"_"+tagMode+"_"+FittingMode+"_SignalOnly_fixedBkg_"+wjetsMode+"_"+angleMode;
	 //---------------- Creating fOutputFolder for fixedBkg
	 if (stat(fOutputFolderCali.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	   std::cout<<"directory: "<< fOutputFolderCali << "/ exist ...!"<< std::endl; 
	 else{
	   const int dir_err = mkdir(fOutputFolderCali.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	   if(dir_err==-1)
	     std::cout<<"Error creating directory "<< fOutputFolderCali << "/"<< std::endl; 
	   else
	     std::cout<<"directory: "<< fOutputFolderCali << "/ created ...!"<< std::endl; 
	 }
       }
       
       else{
	 fOutputFolderCali = "ExternalCalibrationOutput_"+InputChannel+"_"+tagMode+"_"+FittingMode+"_SignalOnly_"+wjetsMode+"_"+angleMode;
	 //---------------- Creating fOutputFolder for floating Bkg
	 if (stat(fOutputFolderCali.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	   std::cout<<"directory: "<< fOutputFolderCali << "/ exist ...!"<< std::endl; 
	 else{
	   const int dir_err = mkdir(fOutputFolderCali.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	   if(dir_err==-1)
	     std::cout<<"Error creating directory "<< fOutputFolderCali << "/"<< std::endl; 
	   else
	     std::cout<<"directory: "<< fOutputFolderCali << "/ created ...!"<< std::endl; 
	 }
       } 
     }
     
     if(Systematic == "Datafit"){
       if(fixedBkg){
	 fOutputFolderData = "ExternalDatafitOutput_"+InputChannel+"_"+tagMode+"_"+FittingMode+"_fixedBkg_"+wjetsMode+"_"+angleMode;
	 if (stat(fOutputFolderData.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	   std::cout<<"directory: "<< fOutputFolderData << "/ exist ...!"<< std::endl; 
	 else{
	   const int dir_err = mkdir(fOutputFolderData.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	   if(dir_err==-1)
	     std::cout<<"Error creating directory "<< fOutputFolderData << "/"<< std::endl; 
	   else
	     std::cout<<"directory: "<< fOutputFolderData << "/ created ...!"<< std::endl; 
	 }
       }
       else
	 fOutputFolderData = "ExternalDatafitOutput_"+InputChannel+"_"+tagMode+"_"+FittingMode+"_"+wjetsMode+"_"+angleMode;
       if (stat(fOutputFolderData.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	 std::cout<<"directory: "<< fOutputFolderData << "/ exist ...!"<< std::endl; 
       else{
	 const int dir_err = mkdir(fOutputFolderData.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	 if(dir_err==-1)
	   std::cout<<"Error creating directory "<< fOutputFolderData << "/"<< std::endl; 
	 else
	   std::cout<<"directory: "<< fOutputFolderData << "/ created ...!"<< std::endl; 
       }
     }
     
     
     //---------------- Creating fOutputFolder for systematics
     if(Systematic != "Calibration" && Systematic != "Datafit")
       {
	 if(fixedBkg){
	   fOutputFolder = "ExternalSystematicsOutput_"+InputChannel+"_"+tagMode+"_"+FittingMode+"_fixedBkg_"+wjetsMode+"_"+angleMode;
	   
	   if (stat(fOutputFolder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	     std::cout<<"directory: "<< fOutputFolder << "/ exist ...!"<< std::endl; 
	   else{
	     const int dir_err = mkdir(fOutputFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	     if(dir_err==-1)
	       std::cout<<"Error creating directory "<< fOutputFolder << "/"<< std::endl; 
	     else
	       std::cout<<"directory: "<< fOutputFolder << "/ created ...!"<< std::endl; 
	   }
	 }
	 
	 else{
	   fOutputFolder = "ExternalSystematicsOutput_"+InputChannel+"_"+tagMode+"_"+FittingMode+"_"+wjetsMode+"_"+angleMode;
	   
	   if (stat(fOutputFolder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	     std::cout<<"directory: "<< fOutputFolder << "/ exist ...!"<< std::endl; 
	   else{
	     const int dir_err = mkdir(fOutputFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	     if(dir_err==-1)
	       std::cout<<"Error creating directory "<< fOutputFolder << "/"<< std::endl; 
	     else
	       std::cout<<"directory: "<< fOutputFolder << "/ created ...!"<< std::endl; 
	   }
	   
	 }
       }

     
     std::string OutputTxtFile = fOutputFolder+"/SystematicOutput_"+InputChannel+".txt";
     std::string pseudoData;
     if (wjetsMode=="1W")
       pseudoData="PseudoData";
     else if(wjetsMode=="3W")
       pseudoData="PseudoData3W";

     if(Systematic == "Calibration")
       fProf.SetOutputFolder(fOutputFolderCali);
     else if(Systematic == "Datafit")
       fProf.SetOutputFolder(fOutputFolderData);
     else
       fProf.SetOutputFolder(fOutputFolder); // for systematics

     fProf.SetOutputTxtFile(OutputTxtFile);
                                               
     fProf.SetNominalInputFile(TemplateFile);

     fProf.SetTemplateFittingMode(FittingMode);

     // lumi in pb
     //fProf.SetLumi(4655.74); // mjk
     fProf.SetLumi(20276.9); // mjk : //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TopData12

     // xsec in pb
     //fProf.SetXsec(96);
     
     if(useNewSample)
      fProf.SetXsec(114.48); // for 110404 sample
     else
      fProf.SetXsec(114.49); // for 117050 sample
     

     if(IsPDF){
       
       std::cout << "HIER!!! PDF TEST!!! " << Systematic.c_str() << std::endl;
       fProf.AddSystematicPseudoData(Systematic, TemplateFile, pseudoData+"_"+Systematic);
       fProf.SetFractions(0.698, 0.301, 0.00041);

     }
     else{

       if(Systematic == "Calibration"){
	 
	 std::cout << "Systematic == Calibration!" << std::endl;
	 
	 if(FittingMode == "3D"){
	   
	   //ProducePseudoData( F0,  FL,  FR, vector<string> HistoNames, vector<string> VarNames, FileName,  OutputFileName)
	   fProf.ProducePseudoData(0.4, 0.45,  0.15, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.4_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.5, 0.40,  0.10, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.5_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.6, 0.35,  0.05, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.6_"+InputChannel+".root"); 
	   fProf.ProducePseudoData(0.7, 0.30,  0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.7_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.8, 0.25, -0.05, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.8_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.9, 0.20, -0.10, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.9_"+InputChannel+".root");
	   fProf.ProducePseudoData(1.0, 0.15, -0.15, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_1.0_"+InputChannel+".root"); 
	   
	   //TEST
	   // fProf.ProducePseudoData( 0.55, 0.447 , 0.003,   HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.4_"+InputChannel+".root");
	   // fProf.ProducePseudoData( 0.60, 0.398 , 0.002,   HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.5_"+InputChannel+".root");
	   // fProf.ProducePseudoData( 0.65, 0.349 , 0.001,   HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.6_"+InputChannel+".root"); 
	   // fProf.ProducePseudoData( 0.70, 0.30 ,  0.0,   HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.7_"+InputChannel+".root");
	   // fProf.ProducePseudoData( 0.75, 0.251 , -0.001,    HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.8_"+InputChannel+".root");
	   // fProf.ProducePseudoData( 0.80, 0.202 , -0.002,    HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.9_"+InputChannel+".root");
	   // fProf.ProducePseudoData( 0.85, 0.153 , -0.003,    HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_1.0_"+InputChannel+".root"); 
	 }
	 else if(FittingMode == "2D"){
	   
	   fProf.ProducePseudoData(0.475, 0.525, 0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.475_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.550, 0.450, 0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.550_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.625, 0.375, 0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.625_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.700, 0.300, 0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.700_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.775, 0.225, 0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.775_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.850, 0.150, 0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.850_"+InputChannel+".root");
	   fProf.ProducePseudoData(0.925, 0.075, 0.00, HistoNames, VarNames, TemplateFile, fOutputFolderCali+"/CalibrationOutput_"+SignalSampleNumber+"_F0_0.925_"+InputChannel+".root");
	   
	 }
	 	 
	 if(FittingMode == "3D")
	   fProf.SetFractions(0.4, 0.45, 0.15);
      //fProf.SetFractions(0.55, 0.447 , 0.003);
   else
     fProf.SetFractions(0.475, 0.525, 0.00); // for 2D
	 
   if(FittingMode == "3D")
     fProf.AddSystematicPseudoData("F0=0.4",   fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.4_"+InputChannel+".root", "PseudoData");
   else
     fProf.AddSystematicPseudoData("F0=0.475", fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.475_"+InputChannel+".root", "PseudoData");

   if(FittingMode == "3D")
     fProf.SetFractions(0.5, 0.40, 0.10);
     //fProf.SetFractions(0.60, 0.398 , 0.002);
   else
     fProf.SetFractions(0.550, 0.450, 0.00);
  
   if(FittingMode == "3D")
     fProf.AddSystematicPseudoData("F0=0.5",   fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.5_"+InputChannel+".root", "PseudoData");
   else
     fProf.AddSystematicPseudoData("F0=0.550", fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.550_"+InputChannel+".root", "PseudoData");

   if(FittingMode == "3D")
     fProf.SetFractions(0.6, 0.35, 0.05);
    //fProf.SetFractions(0.65, 0.349 , 0.001);
   else
     fProf.SetFractions(0.625, 0.375, 0.00);

   if(FittingMode == "3D")
     fProf.AddSystematicPseudoData("F0=0.6",   fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.6_"+InputChannel+".root", "PseudoData"); 
   else
     fProf.AddSystematicPseudoData("F0=0.625", fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.625_"+InputChannel+".root", "PseudoData");
	 
   if(FittingMode == "3D")
     fProf.SetFractions(0.70,  0.30, 0.00);
    //fProf.SetFractions(0.70,  0.30, 0.0);
   else
     fProf.SetFractions(0.7, 0.30, 0.00);

   if(FittingMode == "3D")
     fProf.AddSystematicPseudoData("F0=0.7",   fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.7_"+InputChannel+".root", "PseudoData");
   else
     fProf.AddSystematicPseudoData("F0=0.700", fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.700_"+InputChannel+".root", "PseudoData");

   
   if(FittingMode == "3D")
    fProf.SetFractions(0.8, 0.25, -0.05);
    //fProf.SetFractions(0.75, 0.251 , -0.001);
   else
     fProf.SetFractions(0.775, 0.225,  0.00);

   if(FittingMode == "3D")
     fProf.AddSystematicPseudoData("F0=0.8",   fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.8_"+InputChannel+".root", "PseudoData");
   else
     fProf.AddSystematicPseudoData("F0=0.775", fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.775_"+InputChannel+".root", "PseudoData");

   if(FittingMode == "3D")
    fProf.SetFractions(0.9, 0.20, -0.10);
    //fProf.SetFractions(0.80, 0.202 , -0.002);
   else
     fProf.SetFractions(0.850, 0.150,  0.00);

   if(FittingMode == "3D")
     fProf.AddSystematicPseudoData("F0=0.9",   fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.9_"+InputChannel+".root", "PseudoData"); 
   else
     fProf.AddSystematicPseudoData("F0=0.850", fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.850_"+InputChannel+".root", "PseudoData");


   if(FittingMode == "3D")
     fProf.SetFractions(1.0, 0.15, -0.15);
   //fProf.SetFractions(0.85, 0.153 , -0.003);
   else
     fProf.SetFractions(0.925, 0.075,  0.00);

   if(FittingMode == "3D")
     fProf.AddSystematicPseudoData("F0=1.0",   fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_1.0_"+InputChannel+".root", "PseudoData"); 
   else
     fProf.AddSystematicPseudoData("F0=0.925", fOutputFolderCali+"/CalibrationOutput_"+OtherSignalSampleNumber+"_F0_0.925_"+InputChannel+".root", "PseudoData");
   

   fProf.SetEvaluationMode("LinearDependence");
	       
       }
       
       else if(Systematic == "TemplateStat"){
	 
	 fProf.AddSystematicPseudoData("TemplateStat",     TemplateFile, pseudoData);
	 
	 fProf.SetFractions(0.698, 0.301, 0.00041);
	 
	 
	 /*       fProf.AddSystematicPseudoData("TemplateStatFL",     TemplateFile, pseudoData);
		  if(FittingMode == "3D")
		  fProf.AddSystematicPseudoData("TemplateStatFR",     TemplateFile, pseudoData);
		  
		  if(InputChannel == "el" || InputChannel == "mu"){
		  
		  fProf.AddSystematicPseudoData("TemplateStatWjets",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatQCD",    TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatRemBkg", TemplateFile, pseudoData);
		  
		  }
		  else if(InputChannel == "el_mu"){
		  
		  fProf.AddSystematicPseudoData("TemplateStatWjets_el",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatWjets_mu",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatQCD_el",    TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatQCD_mu",    TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatRemBkg_el", TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatRemBkg_mu", TemplateFile, pseudoData);
		  
		  }
		  else if(InputChannel == "el_BTag" || InputChannel == "mu_BTag"){
		  
		  fProf.AddSystematicPseudoData("TemplateStatWjets_0excl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatWjets_1excl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatWjets_2incl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatQCD_0excl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatQCD_1excl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatQCD_2incl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatRemBkg_0excl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatRemBkg_1excl",  TemplateFile, pseudoData);
		  fProf.AddSystematicPseudoData("TemplateStatRemBkg_2incl",  TemplateFile, pseudoData);
		  
		  }
		  
		  
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  if(FittingMode == "3D")
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  
		  if(InputChannel == "el_mu"){
		  
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  
		  }
		  else if(InputChannel == "el_BTag" || InputChannel == "mu_BTag"){
		  
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  fProf.SetFractions(0.698, 0.301, 0.00041);
		  
		  } */
	 
	 fProf.SetEvaluationMode("FractionWidth");
	 
       }
       else if(Systematic == "Datafit"){
	 
	 fProf.AddSystematicPseudoData("Data", TemplateFile, "Data"); //(std::string Syst, std::string input_file, std::string HistoName)
	 //fProf.SetFractions(0.698, 0.301, 0.00041);
   fProf.SetFractions(m_f0, m_fL, m_fR); //--- changed by MJK (14 Mar. 2016)
	 
       }
       //-----------------------------------------------new syst (23 Nov. 2015 by MJK)       
       else if(Systematic == "ees") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "err") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "mu_idres") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "mu_msres") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "mu_scale") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       //else if(Systematic == "mu_scaleup") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       //else if(Systematic == "mu_scaledown") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "jer") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "TopMass") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "Wjets_iqopt3") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "Wjets_ptjmin10") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "jeff") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "jvf") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       /*else if(Systematic == "jer_NP0") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP4") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP5") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP6") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP7") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
	 else if(Systematic == "jer_NP8") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}*/
       else if(Systematic == "jes_BJESUncert") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}     
       else if(Systematic == "jes_EffectiveNP_Detector1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Detector2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Detector3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Mixed1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Mixed2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Mixed3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Mixed4") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Modelling1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Modelling2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Modelling3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Modelling4") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Statistical1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Statistical2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Statistical3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EffectiveNP_Statistical4") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EtaIntercalibration_Modelling") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_EtaIntercalibration_TotalStat") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_FlavourComp") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_FlavourResponse") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_MuOffsetTerm") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_NPVOffsetTerm") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_PileupPtTerm") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_PunchThrough") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_RhoTopology") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "jes_SingleParticle_HighPt") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "met_res_soft") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "met_sc_soft") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "ELE_RECO") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "ELE_ID") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "ELE_TRIGGER") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "MUON_RECO") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "MUON_ID") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "MUON_TRIGGER") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_bTag0") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_bTag1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_bTag2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_bTag3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_bTag4") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_bTag5") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_cTag0") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_cTag1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_cTag2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_cTag3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag0") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag1") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag2") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag3") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag4") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag5") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag6") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag7") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag8") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag9") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag10") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "BTAG_misTag11") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       
       //-------------------- Modeling syst. -------------------
       else if(Systematic == "Radiation") {runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);}
       else if(Systematic == "MCgenerator") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "PartonShower") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "ColorReconnection") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "UnderlyingEvent") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       //-------------------------------------------------------
       
       //---------------------- PDF syst. ----------------------
       else if(Systematic == "CT10_0") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_1") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_2") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_3") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_4") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_5") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_6") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_7") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_8") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_9") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_10") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_11") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_12") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_13") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_14") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_15") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_16") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_17") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_18") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_19") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_20") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_21") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_22") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_23") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_24") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_25") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "CT10_26") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_0") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_1") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_2") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_3") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_4") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_5") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_6") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_7") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_8") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_9") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_10") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_11") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_12") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_13") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_14") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_15") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_16") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_17") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_18") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_19") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "MSTW_20") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_0") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_1") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_2") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_3") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_4") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_5") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_6") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_7") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_8") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_9") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_10") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_11") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_12") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_13") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_14") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_15") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_16") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_17") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_18") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_19") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_20") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_21") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_22") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_23") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_24") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_25") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_26") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_27") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_28") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_29") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_30") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_31") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_32") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_33") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_34") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_35") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_36") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_37") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_38") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_39") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_40") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_41") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_42") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_43") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_44") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_45") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_46") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_47") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_48") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_49") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       else if(Systematic == "NNPDF_50") {runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}
       //---------------------- single syst. ----------------------
       else if(Systematic.find("single") !=std::string::npos) {runSinglePseudoData( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);}

       //------------ ptRew syst
       else if(Systematic == "hdamp"){
	 runSingleSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, wjetsMode);      
       }
       else if(Systematic == "ptRew-ISR-FSR"){
	 runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);      
       }
       else if(Systematic == "ptRew-PartonShower"){
	 runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);      
       }
       else if(Systematic == "ptRew-MCgen"){
	 runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);      
       }
       else if(Systematic == "ptRew-jer"){
	 runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);      
       }
       else if(Systematic == "ptRew-bjes"){
	 runUpDownSystematicVariation( fProf, TemplateFile, Systematic, m_f0, m_fL, m_fR, pseudoData, wjetsMode);      
       }
             
       //-------------------------------------------------------
       
     }
     
     // lumi in pb
     //fProf.SetLumi(4655.74);
     fProf.SetLumi(20276.9); // mjk: //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TopData12
     
     // xsec in pb
     //fProf.SetXsec(96);
     if(useNewSample)
      fProf.SetXsec(114.48); // for 110404 sample
     else
      fProf.SetXsec(114.49); // for 117050 sample

     
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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<std::string> ReadInputFiles(const char * filename)
{
  // define input file 
  std::ifstream inputfile; 

  // open file 
  inputfile.open(filename); 

  // check if file is open 
  if (!inputfile.is_open())
    {
      std::cout << "TemplateMaker::ReadInputFiles(). File \"" << filename << "\" not found." << std::endl;              
    }
  
  std::vector<std::string> filenameVec;

  std::string line;
  while ( std::getline(inputfile, line) ) {
    if ( !line.empty() ){
      std::size_t found =line.find_first_of("#");
      if(found==std::string::npos)
        filenameVec.push_back(line);
    }
  }
  std::cout << "... "<<filenameVec.size()<<" file(s) per lepton listed "<< "... "<<std::endl;

  return filenameVec;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void runUpDownSystematicVariation( ProfilingClass &fProf, std::string TemplateFile, std::string Systematic, double m_f0, double m_fL, double m_fR, std::string pseudoData, std::string wjetsMode)
{

  //-----------read information about syst histos from systConfig.txt
  std::string List_helper = "systConfig_Sig.txt";
  std::vector<std::string> syst_list = ReadInputFiles(List_helper.c_str());
  std::stringstream test("this_is_a_test_string");
  std::string segment;
  std::vector<std::string> seglist;
  
  std::string varName_up="";
  std::string varName_down="";
  std::string pseudoData_up="";
  std::string pseudoData_down="";
  
  //------------- read information and split strings according to need
  for(int i=0; i < syst_list.size() ; i++){
    
    std::stringstream test(syst_list.at(i));
    std::string segment;
    std::vector<std::string> seglist;
    
    while(std::getline(test, segment, ','))
      {
	seglist.push_back(segment);
      }
    //for(int j=0; j < seglist.size() ; j++) {cout<<"++ j: "<<j<<" | "<<seglist.at(j)<<endl;}
    
    // set names for histograms and pseudoData to read if requested systematic variation is found
    if(seglist.at(0)==Systematic){
      varName_up=seglist.at(1);
      pseudoData_up=seglist.at(2);
      varName_down=seglist.at(3);
      pseudoData_down=seglist.at(4);
    }
    
  }

  if(varName_up==""){
    cout<<"Systematic: "<<Systematic<<" NOT FOUND in systConfigFull.txt . Exiting...."<<endl;
    return;
  }
  
  //add tag for using 3W pseudoData
  if(wjetsMode == "3W"){
    pseudoData_up = add3Wtag(pseudoData_up, "PseudoData", "PseudoData3W");
    pseudoData_down = add3Wtag(pseudoData_down, "PseudoData", "PseudoData3W");
  }

  //----------- do other stuff
  //else if(Systematic == "jes_EffectiveNP_Modelling1"){
  fProf.AddSystematicPseudoData("Nominal", TemplateFile, pseudoData);
  fProf.AddSystematicPseudoData( varName_up, TemplateFile, pseudoData_up);
  fProf.AddSystematicPseudoData( varName_down, TemplateFile, pseudoData_down);
         
  fProf.SetFractions(m_f0, m_fL, m_fR);
  fProf.SetFractions(m_f0, m_fL, m_fR);
  fProf.SetFractions(m_f0, m_fL, m_fR);
  
  fProf.SetEvaluationMode("LargestDiff");
  
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void runSinglePseudoData( ProfilingClass &fProf, std::string TemplateFile, std::string Systematic, double m_f0, double m_fL, double m_fR, std::string wjetsMode)
{

  //-----------read information about syst histos from systConfig.txt
  std::string List_helper = "systConfigSingle.txt";
  std::vector<std::string> syst_list = ReadInputFiles(List_helper.c_str());
  std::stringstream test("this_is_a_test_string");
  std::string segment;
  std::vector<std::string> seglist;
  
  std::string varName="";
  std::string pseudoData_name="";
  
  //------------- read information and split strings according to need
  for(int i=0; i < syst_list.size() ; i++){
    
    std::stringstream test(syst_list.at(i));
    std::string segment;
    std::vector<std::string> seglist;
    
    while(std::getline(test, segment, ','))
      {
	seglist.push_back(segment);
      }
    //for(int j=0; j < seglist.size() ; j++) {cout<<"++ j: "<<j<<" | "<<seglist.at(j)<<endl;}
    
    // set names for histograms and pseudoData to read if requested systematic variation is found
    //std::cout<<"!!!! seglist.size() = "<<seglist.size()<<std::endl;
    if(seglist.at(0)==Systematic){
      varName=seglist.at(1);
      pseudoData_name=seglist.at(2);
    }
    
  }
  std::cout<<"%%%% seglist.size() = "<<seglist.size()<<std::endl;
  std::cout<<"%%%% varName = "<<varName<<std::endl;
  std::cout<<"%%%% pseudoData_name = "<<pseudoData_name<<std::endl;

  if(varName==""){
    cout<<"Systematic: "<<Systematic<<" NOT FOUND in systConfigSingle.txt . Exiting...."<<endl;
    return;
  }
  
  //add tag for using 3W pseudoData
  if(wjetsMode == "3W"){
    pseudoData_name = add3Wtag(pseudoData_name, "PseudoData", "PseudoData3W");
  }

  //----------- do other stuff
  //else if(Systematic == "jes_EffectiveNP_Modelling1"){
  std::cout<<"&&&& seglist.size() = "<<seglist.size()<<std::endl;
  fProf.AddSystematicPseudoData( varName, TemplateFile, pseudoData_name);
  std::cout<<"**** seglist.size() = "<<seglist.size()<<std::endl;
  fProf.SetFractions(m_f0, m_fL, m_fR);
  std::cout<<")))) seglist.size() = "<<seglist.size()<<std::endl;

  //fProf.SetEvaluationMode("LargestDiff");
  
}

void runSingleSystematicVariation( ProfilingClass &fProf, std::string TemplateFile, std::string Systematic, double m_f0, double m_fL, double m_fR, std::string wjetsMode)
{

  //-----------read information about syst histos from systConfig.txt
  std::string List_helper = "systConfig_Sig.txt";
  std::vector<std::string> syst_list = ReadInputFiles(List_helper.c_str());
  std::stringstream test("this_is_a_test_string");
  std::string segment;
  std::vector<std::string> seglist;
  
  std::string varName_up="";
  std::string varName_down="";
  std::string pseudoData_up="";
  std::string pseudoData_down="";
  
  //------------- read information and split strings according to need
  for(int i=0; i < syst_list.size() ; i++){
    
    std::stringstream test(syst_list.at(i));
    std::string segment;
    std::vector<std::string> seglist;
    
    while(std::getline(test, segment, ','))
      {
	seglist.push_back(segment);
      }
    //for(int j=0; j < seglist.size() ; j++) {cout<<"++ j: "<<j<<" | "<<seglist.at(j)<<endl;}
    
    // set names for histograms and pseudoData to read if requested systematic variation is found
    if(seglist.at(0)==Systematic){
      varName_up=seglist.at(1);
      pseudoData_up=seglist.at(2);
      varName_down=seglist.at(3);
      pseudoData_down=seglist.at(4);
    }
    
  }
  
  //add tag for using 3W pseudoData
  if(wjetsMode == "3W"){
    pseudoData_up = add3Wtag(pseudoData_up, "PseudoData", "PseudoData3W");
    pseudoData_down = add3Wtag(pseudoData_down, "PseudoData", "PseudoData3W");
  }
  
  if(varName_up==""){
    cout<<"Systematic: "<<Systematic<<" NOT FOUND in systConfigFull.txt . Exiting...."<<endl;
    return;
  }
  std::cout<<"varName_up: "<<varName_up<<std::endl;
  std::cout<<"varName_down: "<<varName_down<<std::endl;

  //----------- do other stuff
  //else if(Systematic == "jes_EffectiveNP_Modelling1"){
  fProf.AddSystematicPseudoData( varName_up, TemplateFile, pseudoData_up);
  fProf.AddSystematicPseudoData( varName_down, TemplateFile, pseudoData_down);
         
  fProf.SetFractions(m_f0, m_fL, m_fR);
  fProf.SetFractions(m_f0, m_fL, m_fR);
  
  fProf.SetEvaluationMode("LargestDiff");
  
}

std::string add3Wtag (std::string &s, const std::string &toReplace, const std::string &replaceWith)
{
  std::cout<<"before"<<std::endl;
  std::cout<<"s: "<<s<<"\ttoReplace: "<<toReplace<<"\treplaceWith: "<<replaceWith<<std::endl;
  return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
  std::cout<<"after"<<std::endl;
}
