#ifndef PlotInterpolationCurves_H_
#define PlotInterpolationCurves_H_

#include "TemplateInfo.h"
#include "SystematicInfo.h"
//#include "ProfilingClass.h"

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"

#include <map>
#include <vector>
#include <string>

class PlotInterpolationCurves{

 public:

  PlotInterpolationCurves();
  PlotInterpolationCurves(std::vector<TemplateInfo::MySignalTemplate>, std::vector<TemplateInfo::MyBkgTemplate>, std::vector<SystematicInfo::MySystematic>);
  PlotInterpolationCurves(std::vector<TemplateInfo::MySignalTemplate>, std::vector<TemplateInfo::MyBkgTemplate>, int, std::string);

  virtual ~PlotInterpolationCurves();

  void AddFileLabel(std::string label)   {fLabel.push_back(label);};

  void AddTemplateFile(std::string file) {fTemplateFile.push_back(file);};

  void AddKValue(double val){fKValues.push_back(val);};

  void FillGraphs(std::vector<std::vector<TH1D> >, int);
  void FillVectorsBkgNorm(int);
  void FillVectorsNominal();
  void FillVectorsSystematics();

  TH1D *MakePointer(TH1D);

  std::vector<std::vector<TF1> > GetInterpolationFunctions(){return fFitFunc;};
  std::vector<std::vector<TF1> > GetInterpolationFunctions_QuadFit(){return fFitFunc_QuadFit;};

  void PrintGraphs(std::string, int, std::string);
    
  void SetLeptonType(std::string LeptonType){fLeptonType = LeptonType;};
  void SetJetBinLabel(std::string);
  void SetBtagLabel(std::string);
  void SetLumiLabel(std::string);
  void SetOutputFile(std::string output){fOutputFile = output;};
  void SetOutputFolder(std::string output){fOutputFolder = output;};
  //void SetNumberOfBins (int bins){fNBins = bins;};
  //void SetLowerEdge(double low){fLowerEdge = low;};
  //void SetUpperEdge(double up){fUpperEdge = up;};

  void SetSubPlotLabels(std::vector<std::string> Lepton, std::vector<std::string> Jet, std::vector<std::string> BTag){

    fLeptonLabel = Lepton;
    fJetBinLabel = Jet;
    fBTagLabel   = BTag;

  };

  void MakePlots(std::string);
  void DrawFinalPlot(std::vector<TH1D>, std::vector<std::string>, std::string, bool, bool);
  void DrawFinalPlotSingle(std::vector<TH1D>, std::vector<std::string>, std::string, bool, bool);

  void createHTML(std::string);
  void WriteEventYieldComp();

  TGraphErrors* GetRatioDistribution(TH1D*, TH1D*, int);
  TGraphErrors* GetRatioDistribution(std::vector<double>, std::vector<double>, std::vector<double>, TF1);
  TGraphErrors* GetRatioDistribution(std::vector<double>, std::vector<double>, std::vector<double>, TF1, TF1);

  double GetChi2Values(std::vector<double>, std::vector<double>, std::vector<double>, TF1);
  double GetChi2Values(std::vector<double>, std::vector<double>, std::vector<double>, TF1, TF1);


  TH1D *MakeHistogram(std::vector<double, std::allocator<double> >, TF1);

  void MakeChi2Distributions(std::string);

  void SetHtmlLable(std::string htmlLable) {fHTMLLabel=htmlLable;};

 protected:

  std::vector<std::string> fLeptonLabel,fJetBinLabel,fBTagLabel;

  std::map<std::string,std::string> _htmls;

  std::vector<std::string> fTemplateFile;

  std::vector<std::vector<TH1D> >   fHist;
  std::vector<std::vector<TF1> >    fFitFunc;
  std::vector<std::vector<TF1> >    fFitFunc_QuadFit;

  std::vector<std::vector<std::vector<double> > > fHistEntries;
  std::vector<std::vector<std::vector<double> > > fHistEntriesErr;

  std::vector<TemplateInfo::MyBkgTemplate>    fBkgTempl;
  std::vector<TemplateInfo::MySignalTemplate> fSignalTempl;
  std::vector<SystematicInfo::MySystematic >  fSystematics;

  std::vector<double> F0_lin, F0_quad, F0_lin_pos, F0_quadfit;
  std::vector<double> FL_lin, FL_quad, FL_lin_pos, FL_quadfit;
  std::vector<double> FR_lin, FR_quad, FR_lin_pos, FR_quadfit;
  std::vector<double> Wjets_lin, Wjets_quad, Wjets_lin_pos, Wjets_quadfit;
  std::vector<double> QCD_lin, QCD_quad, QCD_lin_pos, QCD_quadfit;
  std::vector<double> RemBkg_lin, RemBkg_quad, RemBkg_lin_pos, RemBkg_quadfit;

  std::vector<double> fKValues;

  std::vector<std::string> fLabel;
  std::vector<std::string> fHistLabels;

  std::string fLeptonType;
  //  std::string fJetBinLabel;
  //  std::string fBtagLabel;
  std::string fLumiLabel;
  std::string fOutputFile;
  std::string fOutputFolder;
  std::string fHTMLLabel;
  std::string fSystLabel;

  int fNKValues;
  int fNHistos;
  //int fNBins;
  double fXLow;
  double fXHigh;

};

#endif
