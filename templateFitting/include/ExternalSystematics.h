#ifndef EXTERNALSYSTEMATICS_H
#define EXTERNALSYSTEMATICS_H

#include "TF1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TRandom3.h"

#include "TemplateInfo.h"
#include "SystematicInfo.h"
#include "ProfilingClass.h"

#include <map>

extern std::vector<double> fBkgNorm2;
extern std::vector<double> fBkgUnc2;
extern std::vector<TH1D> fParameterHist2;
extern std::vector<TH1D> fParameterHistTS;
extern std::vector<double> SelEff2;

class ExternalSystematics
{
 public:

  ExternalSystematics(std::string);      //constructor
  ~ExternalSystematics();     //destructor

  TGraphErrors MakeTGraphErrors(TH1D);

  void clean();
  void CallEvaluation(std::string, std::string);
  TH1D RandomPseudodata(TH1D);
  void ProduceRandomTemplates();

  std::pair<double,double> CalculateParameterRanges(int, bool, bool);

  void FillTemplateVector(std::vector<TemplateInfo::MySignalTemplate>, std::vector<TemplateInfo::MyBkgTemplate>);

  void SetBkgParameters(std::vector<TemplateInfo::MyBkgTemplate> Bkg){fBkgParameter = Bkg;};
  void SetEvaluationMode(std::string mode){fEvaluationMode = mode;};
  void SetSignalParameters(std::vector<TemplateInfo::MySignalTemplate> Signal){fSignalParameter = Signal;};
  void SetSystematicPseudoData(std::vector<SystematicInfo::MySystematicPD> syst){fSystematicPD = syst;};
  void SetSystematicLabel(std::string syst){fSystematicLabel = syst;};

  void SetFitParameters(std::string);
  void SetLumi(double lumi){fLumi = lumi;};
  void SetNumberOfPE(int pe){fNPseudoExp = pe;};
  void SetOutputFolder(std::string folder){fOutputFolder = folder;};
  void SetOutputSampleNumber(int nr){fNOutputSample = nr;};
  void SetValidationMode(std::string mode){fValidationMode = mode;};
  void SetXsec(double xsec){fXsec = xsec;};
  void SetKfactor(double kfactor){fKfactor = kfactor;}; // added by MJK
  //void SetEff(double F0eff, double FLeff, double FReff){fF0eff = F0eff; fFLeff = FLeff; fFReff = FReff; }; // added by MJK
  void SetChannel(std::string channel){fChannel = channel;};

  void SetSubPlotLabels(std::vector<std::string> Lepton, std::vector<std::string> Jet, std::vector<std::string> BTag){

    fLeptonLabel = Lepton;
    fJetBinLabel = Jet;
    fBTagLabel   = BTag;

  };
    
  void SetFractions(std::vector<std::vector<double> > fraction){
    fFraction = fraction;
  };

  void AnalyseDatafit();
  void AnalyseParameter(int, int, std::string);

  void Validation(TH1D, int, std::string);
  void Datafit(std::vector<double>, std::string, std::string);

  void PullDistributions(int, int, double, double, double, double, std::string);
  void CalDistribution(int, int, double, double, double, double, std::string);
  void CalPullDistribution(int, int, double, double, double, double, std::string);
  void ErrorDataDistributions(int, std::string);

  //void EvaluateSystematic(std::string, std::string);
  void EvaluateSystematic(std::string, std::string, std::vector<std::string> );

  void InitializeGraphs();
  void InitializeHistograms(int);
  void InitializeHistogramsAnalysis(int,  int);
  void InitializeHistogramsFractions();

  void AnalyseFractions(int, std::string);

  void ProduceAllPlots(int);

  void ProduceAllFractionPlots();

  TH1D PseudoData(std::vector<double>, std::vector<double>);

  TH1D PseudoDataBestFit(std::vector<double>, std::vector<double>);

  void PlotSystematicDistributions(std::vector<TH1D>, std::vector<std::string>, std::string);

  void createHTML();

  TTree * InitializeOutputTree();

  std::vector<std::string> fSystList;
  std::vector<std::string> fParameterNames;
  std::vector<std::string> fFractionNames;
  std::vector<std::vector<double> > fFraction;

 //protected:

  TH2D *fN0NL,*fN0NR,*fNRNL,*fN0N0,*fNLNL,*fNRNR;
  TH2D *fF0FL,*fF0FR,*fFRFL,*fF0F0,*fFLFL,*fFRFR;
  
  double fErrormatrix[3][3];
  double fCorrelationmatrix[10][10];

  std::vector<std::vector<TGraphErrors> > fGraphs;
  std::vector<std::vector<TGraphErrors> > fFractionGraphs;

  //std::vector<std::vector<double> > fFraction;

  std::vector<std::string> fLeptonLabel,fJetBinLabel,fBTagLabel;

  std::map<std::string,std::string> _htmls;

  std::string fEvaluationMode;

  //std::vector<std::string> fSystList;
  std::vector<double> fNominalScaleParameters;

  std::vector<TemplateInfo::MyBkgTemplate>    fBkgParameter;
  std::vector<TemplateInfo::MySignalTemplate> fSignalParameter;
  std::vector<SystematicInfo::MySystematicPD> fSystematicPD;

//  std::vector<std::string> fParameterNames;
  std::vector<double>      fParameterValues;
  std::vector<double>      fParameterValues_err;
  std::vector<double>      fParameterValues_nom;

//  std::vector<std::string> fFractionNames;
  
  std::vector<double>      fFractionValues;
  std::vector<double>      fFractionValues_err;
  std::vector<double>      fFractionValues_nom;
  std::vector<double>      fFractionValues_pull;
  std::vector<double>      fFractionValues_diff;

  double fVarCorrN0NL;
  double fVarCorrN0NR;
  double fVarCorrNRNL;

  int fNOutputSample;
  int fNPseudoExp;
  int fNNui;

  double fNGes;

  double fLumi, fXsec, fKfactor;
  double fF0eff, fFLeff, fFReff;

  std::vector<double> diffF0;
  std::vector<double> diffFL;
  std::vector<double> diffFR;
  std::vector<double> meanF0;
  std::vector<double> meanFL;
  std::vector<double> meanFR;

  std::string fOutputFolder;
  std::string fFitMethod;
  std::string fValidationMode;
  std::string fSystematicLabel;
  std::string fChannel;
  std::string fHTMLLabel;
  std::string fFittingMode;

  std::vector<TH1D>   fHistoF0;
  std::vector<TH1D>   fHistoFL;
  std::vector<TH1D>   fHistoFR;

  std::vector<std::vector<TH1D> > fHistoAll;
  std::vector<std::vector<TH1D> > fHistoAnalysis;
  std::vector<std::vector<std::vector<TH1D> > > fFractionAnalysis;

  std::vector<TH1D>   fHistNomVec;
  std::vector<TH1D>   fHistNormVec;
  std::vector<double> fEffVec;
  std::vector<double> fBkgExpVec;
  std::vector<double> fBkgUncVec;
  int fNScaleFactors;

  TRandom3* fRandom;
  TRandom3* gaussRandom;

  TH1D fDataHist;
  TH1D fFitErrorHist;

  TF1     *fit;
  TGraph  *graph;
  TMinuit *fitter_prof;

  // variables for output tree
  double out_N0,         out_NL,       out_NR;
  double out_Wjets,      out_QCD,      out_RemBkg;

};

#endif

