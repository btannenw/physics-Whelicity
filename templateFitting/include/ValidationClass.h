#ifndef VALIDATION_H
#define VALIDATION_H

#include "TF1.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TRandom3.h"

#include "TemplateInfo.h"
#include "SystematicInfo.h"
#include "interpolation.h"
#include "ProfilingClass.h"

//extern TH1D hist_ensemble;
extern TH1D hist_data;
extern std::vector<tinterpolation> inter_obj;
extern std::vector<double> fBkgNorm;
extern std::vector<double> fBkgUnc;
extern std::vector<TH1D> fParameterHist;
extern std::vector<double> SelEff;

extern std::string fBkgMode;

class tvalidation
{
 public:

  tvalidation();      //constructor
  ~tvalidation();       //destructor

  void clean();

  void CallValidation(std::string, int, std::string);

  TH1D GetPseudodataFit(double, double, double, std::vector<double>);
  TH1D GetPseudodata(double, double, double);
  TH1D RandomPseudodata(TH1D);

  void FillInterpolationVector();
  void FillTemplateVector(std::vector<TemplateInfo::MySignalTemplate>, std::vector<TemplateInfo::MyBkgTemplate>);

  void SetBkgParameters(std::vector<TemplateInfo::MyBkgTemplate> Bkg){fBkgParameter = Bkg;};
  void SetSignalParameters(std::vector<TemplateInfo::MySignalTemplate> Signal){fSignalParameter = Signal;};
  void SetNuisanceParameters(std::vector<std::string> Nui){fNuisanceParameter = Nui;};
  //  void SetInterpolationMethod(std::string method){fInterpolationMethod = method;};
  void SetFitParameters(std::string);
  void SetLumi(double lumi){fLumi = lumi;};
  void SetNumberOfPE(int pe){fNPseudoExp = pe;};
  void SetOutputSampleNumber(int nr){fNOutputSample = nr;};
  void SetValidationMode(std::string mode){fValidationMode = mode;};
  void SetXsec(double xsec){fXsec = xsec;};

  void Validation(int, int, std::vector<double>, double, double, double, double, std::string, std::string);
  void Datafit(std::vector<double>, std::string, std::string);

  std::vector<double> Validation_corr(int, int, std::vector<double>, double, double, double, double, std::string);

  void PullDistributions(int, int, double, double, double, double, std::string);
  void CalDistribution(int, int, double, double, double, double, std::string);
  void CalPullDistribution(int, int, double, double, double, double, std::string);
  void ErrorDataDistributions(int, std::string);

  double error_profile_PE (TH1D, double, int, std::vector<double> &, std::vector<double> &, int, int, double, double, double, int);

  TTree * InitializeOutputTree();

 protected:

  std::vector<std::string> fNuisanceParameter;
  std::vector<double> fNominalScaleParameters;

  std::vector<TemplateInfo::MyBkgTemplate>    fBkgParameter;
  std::vector<TemplateInfo::MySignalTemplate> fSignalParameter;

  int fNOutputSample;
  int fNPseudoExp;
  int fNNui;

  double fLumi, fXsec;

  std::string fFitMethod;
  //  std::string fInterpolationMethod;
  std::string fValidationMode;

  TTree *fOutputTree;

  double out_F0,             out_FL,           out_FR;
  double out_N0,             out_NL,           out_NR;
  double out_Nges,           out_Nges_nom;
  double out_Wjets,          out_QCD,          out_RemBkg;
  double out_F0_nom,         out_FL_nom,       out_FR_nom;
  double out_N0_nom,         out_NL_nom,       out_NR_nom;
  double out_Wjets_nom,      out_QCD_nom,      out_RemBkg_nom;
  double out_F0_err,         out_FL_err,       out_FR_err;
  double out_N0_err,         out_NL_err,       out_NR_err;
  double out_Wjets_err,      out_QCD_err,      out_RemBkg_err;
  double out_N0_prof_err,    out_NL_prof_err,  out_NR_prof_err;
  double out_Wjets_prof_err, out_QCD_prof_err, out_RemBkg_prof_err;

  std::vector<double> out_nui;
  std::vector<double> out_nui_nom;
  std::vector<double> out_nui_err;
  std::vector<double> out_nui_prof_err;

  int out_MinuitStatus;
  int out_NumNuisancePar;
  std::vector<TH1D>   fHistNomVec;
  std::vector<TH1D>   fHistNormVec;
  std::vector<double> fEffVec;
  std::vector<double> fBkgExpVec;
  std::vector<double> fBkgUncVec;
  int fNScaleFactors;
  int fNNuiParameters;
  TRandom3* fRandom;
  TRandom3* gaussRandom;

  TF1     *fit;
  TGraph  *graph;
  TMinuit *fitter_prof;

  //  int fNBins;

};

#endif

