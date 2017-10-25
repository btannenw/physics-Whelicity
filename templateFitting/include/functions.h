#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "TemplateInfo.h"
#include "SystematicInfo.h"

extern std::vector<std::vector<TemplateInfo::SigInterpolObject> > fSignalInterpolation;
extern std::vector<std::vector<TemplateInfo::BkgInterpolObject> > fBkgInterpolation;

extern std::string fInterpolMethod;

extern std::vector<TFile*> input_files_down;
extern std::vector<TFile*> input_files_up;

extern int fNNui;
extern int fNBkg;
extern int fNSignal;

//void SetNumberOfBkgHist(int bkg){fNBkg = bkg;};
//void SetNumberOfSignalHist(int sig){fNSignal = sig;};

void SetInterpolationMethod(std::string);
std::string GetInterpolationMethod();

void AddBkgSystematic(std::vector<TemplateInfo::MyBkgTemplate>, std::vector<SystematicInfo::MySystematic>);
void AddSignalSystematic(std::vector<TemplateInfo::MySignalTemplate>, std::vector<SystematicInfo::MySystematic>);

TH1D* normalise(TH1D*);
TH1D normalise(TH1D);
TH1D normalise_fast(TH1D);
void CreateDummyFiles(std::vector<TH1D*> &);

//void eval_chisqrt(int &, double *, double &, double *, int);

	
#endif
