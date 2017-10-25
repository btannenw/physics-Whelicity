#ifndef PLOTS_H
#define PLOTS_H

#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"

void PlotDataFit(TH1D, TH1D, TH1D, TGraphErrors, TH1D, std::string, std::string);
void PlotDataFitSimple(TH1D, TH1D, TH1D, TGraphErrors, std::string, std::string);
void PlotDataFit (TH1D, TH1D, TH1D, TH1D, TGraphAsymmErrors*, TH1D, std::string, std::string); 
void PlotDistribution (TH1D*, double, std::string, std::string, std::string, std::string, TF1*, int, int, std::string);
void PlotDistribution_N (TH1D*, double, double, double, std::string, std::string, std::string, std::string, int, std::string, std::string);
void PlotErrorDistribution (TH1D*, double, std::string, std::string, double, std::string, std::string, int);
void PlotErrorDataComp (TH1D*, double, double, std::string, std::string, std::string, std::string, int);
void Plot2Dcorrelation_F (TH2D*, double, double, std::string, std::string, std::string, std::string);
void Plot2Dcorrelation_N (TH2D*, double, std::string, std::string, std::string, std::string);
void PlotPullDistribution (TH1D*, double, std::string, std::string, double, std::string, std::string, TF1*);
void PlotCalDistribution (TGraphErrors*, TGraphErrors*, std::string, std::string, std::string, double, double, double, double, std::string, int, std::string);
void PlotCalDistribution_N (TGraphErrors*, TGraphErrors*, std::string, std::string, std::string, double, double, std::string, std::string, std::string);
void PlotCalPullDistribution (TGraphErrors*, std::string, std::string, std::string, double, double, double, double, std::string, std::string, std::string);
void PlotCalRMSDistribution (TGraphErrors*, std::string, std::string, std::string, double, double, double, double, std::string, std::string, std::string);
void PlotPileupDistribution (TGraphErrors* , std::string, std::string, std::string, double, double, double, double, std::string, std::string, std::string);

#endif
