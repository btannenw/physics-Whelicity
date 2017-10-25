#include "PlotInterpolationCurves.h"
#include "AtlasStyle.h"
#include "StatusLogbook.h"
#include "functions.h"
#include "ProfilingClass.h"


#include <sys/stat.h>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TLine.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

using namespace std;

PlotInterpolationCurves::PlotInterpolationCurves()
{


}

//PlotInterpolationCurves object, pass information about templates and used systematics
PlotInterpolationCurves::PlotInterpolationCurves(std::vector<TemplateInfo::MySignalTemplate> SigTempl, std::vector<TemplateInfo::MyBkgTemplate> BkgTempl, std::vector<SystematicInfo::MySystematic> Syst)
{
 
  fNKValues = Syst.size()*2+1;
  fNHistos  = SigTempl.size() + BkgTempl.size();

  fSignalTempl = SigTempl;
  fBkgTempl    = BkgTempl;
  fSystematics = Syst;

  fHist.clear();
  fLabel.clear();
  fHistLabels.clear();
  fKValues.clear();

  fOutputFolder = fSystematics[0].Filename_up;
  fOutputFolder.replace(fOutputFolder.find("u_1k.root"), strlen("u_1k.root"), "");

  fSystLabel    = fSystematics[0].SystematicType;

  fHTMLLabel    = "index4jincl_"+fSystLabel+".html";

  gErrorIgnoreLevel = kError;

  //  createHTML("index");
  //  createHTML("index5jincl");
  //  createHTML("index4jincl");

  TH1::AddDirectory(kFALSE);
  
}

//PlotInterpolationCurves object, pass information about templates in case no additional systematic effects are considered
PlotInterpolationCurves::PlotInterpolationCurves(std::vector<TemplateInfo::MySignalTemplate> SigTempl, std::vector<TemplateInfo::MyBkgTemplate> BkgTempl, int NrBkg, std::string OutputFolder)
{

  fNKValues = 3;
  fNHistos  = SigTempl.size() + BkgTempl.size();

  fSignalTempl = SigTempl;
  fBkgTempl    = BkgTempl;
 
  fHist.clear();
  fLabel.clear();
  fHistLabels.clear();
  fKValues.clear();

  fSystLabel    = BkgTempl[NrBkg].ParamName+"Norm";

  fOutputFolder = OutputFolder+"_"+fSystLabel+"/";

  fHTMLLabel    = "index4jincl_"+fSystLabel+".html";

  gErrorIgnoreLevel = kError;

  TH1::AddDirectory(kFALSE);

}


PlotInterpolationCurves::~PlotInterpolationCurves()
{

}

//Fill vector fHist containing all input histograms if bkg is fitted via nuisance parameters
void PlotInterpolationCurves::FillVectorsBkgNorm(int NrBkg)
{

  std::vector<TH1D> HelpVec;

  fKValues.push_back(0.0);
  fKValues.push_back(1.0);
  fKValues.push_back(-1.0);
  
  fLabel.push_back("nominal");
  fLabel.push_back("uncert_up");
  fLabel.push_back("uncert_down");

  for(int k = 0; k < fKValues.size(); ++k){
    
    for(int iSig = 0; iSig < fSignalTempl.size(); ++iSig){
      
      std::string label = fSignalTempl[iSig].ParamName;
      
      fHistLabels.push_back(label);
           
      HelpVec.push_back(fSignalTempl[iSig].hist);
      
    }
    
    for(int iBkg = 0; iBkg < fBkgTempl.size(); ++iBkg){
      
      std::string label = fBkgTempl[iBkg].ParamName;
      
      fHistLabels.push_back(label);
      
      double rel_uncertainty = fBkgTempl[iBkg].Unc/fBkgTempl[iBkg].Exp;
      
      TH1D up   = fBkgTempl[iBkg].hist;    
      TH1D down = fBkgTempl[iBkg].hist;
      up.Sumw2();
      down.Sumw2();
      
      up.Scale(1.0+rel_uncertainty);   
      down.Scale(1.0-rel_uncertainty);

      fNBins = up.GetNbinsX();
      
      if(NrBkg == iBkg){

	if(k == 0)
	  HelpVec.push_back(fBkgTempl[iBkg].hist);
	else if(k == 1)
	  HelpVec.push_back(up);
	else if(k == 2){
	 
	  if(down.Integral() < 10e-5){
	   
	    for(int i = 1; i <= fNBins; ++i)
	      down.SetBinError(i, up.GetBinError(i));

	  }

	  HelpVec.push_back(down);
	
	}

	//	std::cout << iBkg << "\t" << NrBkg << "\t" << down.Integral() << "\t" << fBkgTempl[iBkg].hist.Integral() << "\t" << up.Integral() << std::endl;
	
      }
      else{
	
	HelpVec.push_back(fBkgTempl[iBkg].hist);
	
	//	std::cout << iBkg << "\t" << NrBkg << "\t" << HelpVec[iBkg].Integral() << std::endl;
	
      }
      
    }

    fHist.push_back(HelpVec);
    
    HelpVec.clear();
  
  }

}

//Fill vector fHist with nominal histograms of one certain systematic effect (see ProfilinClass.cxx, approx. l. 293) 
void PlotInterpolationCurves::FillVectorsNominal()
{

  std::vector<TH1D> HelpVec;

  fKValues.push_back(0);
  fLabel.push_back("Nominal");

  for(int iSig = 0; iSig < fSignalTempl.size(); ++iSig){

    std::string label = fSignalTempl[iSig].ParamName;

    fHistLabels.push_back(label);

    HelpVec.push_back(fSignalTempl[iSig].hist);

  }
  for(int iBkg = 0; iBkg < fBkgTempl.size(); ++iBkg){

    std::string label = fBkgTempl[iBkg].ParamName;

    fHistLabels.push_back(label);

    HelpVec.push_back(fBkgTempl[iBkg].hist);
      
  }

  fHist.push_back(HelpVec);

  HelpVec.clear();

}

//Fill vector fHist with information concerning one certain systematic effect
void PlotInterpolationCurves::FillVectorsSystematics()
{

  for(int iSyst = 0; iSyst < fSystematics.size(); ++iSyst){

    std::vector<TH1D> HelpVec;
    std::vector<TH1D> HelpVec2;

    std::vector<TH1D> CurrentVecUp   = fSystematics[iSyst].HistUp;
    std::vector<TH1D> CurrentVecDown = fSystematics[iSyst].HistDown;

    std::vector<std::string> Param = fSystematics[iSyst].ParamName;

    bool all_good = true;

    for(int iHist = 0; iHist < CurrentVecUp.size(); ++iHist){

      // check if the labels in the structure match the ones in the vector of Syst labels
      if(fHistLabels[iHist] != Param[iHist])
	all_good = false;
      
    }
    
    // only go on if all uncertainties are correctly assigned
    if(all_good){
      
      for(int iHist = 0; iHist < CurrentVecUp.size(); ++iHist){
	
	HelpVec.push_back(CurrentVecUp[iHist]);
	
      }
      fHist.push_back(HelpVec);

      //      HelpVec.clear();

      for(int iHist = 0; iHist < CurrentVecDown.size(); ++iHist){

        HelpVec2.push_back(CurrentVecDown[iHist]);

      }

      fHist.push_back(HelpVec2);
      
      fKValues.push_back(fSystematics[iSyst].kup);
      fKValues.push_back(fSystematics[iSyst].kdown);

      std::stringstream nr_up, nr_down;
      nr_up   << fSystematics[iSyst].kup;
      nr_down << fSystematics[iSyst].kdown;
      
      fLabel.push_back(fSystLabel+" ("+nr_up.str()+"#sigma)");
      fLabel.push_back(fSystLabel+" ("+nr_down.str()+"#sigma)");

    }
    else{

      WriteErrorStatus("PlotInterpolationCurves", "Systematics are not properly assigned!!!");

    }

    HelpVec.clear();
    HelpVec2.clear();
    
  }
  
}

//Function needed to make the different desired plots
void PlotInterpolationCurves::MakePlots(std::string fInterpolationMethod)
{
     cout << "fNBins " << fNBins << endl;
  std::vector<double> entries;
  std::vector<std::vector<double> > bins;

  for(int iK = 0; iK < fKValues.size(); ++iK)
    entries.push_back(-100.0);

  for(int iBin = 0; iBin < fNBins; ++iBin)
    bins.push_back(entries);

  for(int iHist = 0; iHist < fNHistos; ++iHist){
    fHistEntries.push_back(bins);
    fHistEntriesErr.push_back(bins);
  }

  createHTML("index4jincl");

  for(int iHist = 0; iHist < fNHistos; ++iHist){

    //Fill information of histograms in corresponding 3D-double vector; use fHistEntries
    this -> FillGraphs(fHist, iHist);

  }

  for(int iHist = 0; iHist < fNHistos; ++iHist){

    //Plot different interpolation curves in one canvas, calculate and fill fFitFunc
    this -> PrintGraphs(fHistLabels[iHist].c_str(), iHist, fInterpolationMethod);

  }

  std::string output = fOutputFolder;

  //Need different order in 2D- and 3D-vector, respectively.
  std::vector<std::vector<TH1D> > fHist_reverse;
  std::vector<TH1D> hist_helper;

  //for(int iHist = 0; iHist < fNHistos; ++iHist){
  for (int iHist = 0; iHist < fNHistos; ++iHist){

	hist_helper.clear();

	for(int iK = 0; iK < fHist.size(); ++iK){
	
		hist_helper.push_back(fHist[iK][iHist]);	
	}
  	fHist_reverse.push_back(hist_helper);
  }

  //std::string outputFile     = fOutputFolder+"/Interpolation_"+name+"_bin"+nr_bin.str()+".eps";
  std::string output_histos = fOutputFolder+"/Templates";

  for(int iHist = 0; iHist < fNHistos; ++iHist){

	//Draw plot which contains all different templates (including variations)
    this -> DrawFinalPlot(fHist_reverse[iHist], fLabel, output_histos+"_"+fHistLabels[iHist]+".png", true, false);
    
    
  }

  //Plot the different chi^2 values of all bins and all interpolation methods
  MakeChi2Distributions("N0");
  MakeChi2Distributions("NL");
  MakeChi2Distributions("NR");
  MakeChi2Distributions("Wjets");
  MakeChi2Distributions("QCD");
  MakeChi2Distributions("RemBkg");
 
}


//Plot chi^2  values of all bins and interpolation methods
void PlotInterpolationCurves::MakeChi2Distributions(std::string name)
{

  AtlasStyle();

  double val_lin[fNBins], val_quad[fNBins], val_lin_pos[fNBins], val_quadfit[fNBins];
 //double val_lin[15], val_quad[15], val_lin_pos[15], val_quadfit[15];

  std::vector<double> lin, lin_pos, quad, quadfit;

  if(name == "N0"){
    lin = F0_lin; lin_pos = F0_lin_pos; quad = F0_quad; quadfit = F0_quadfit;
  }
  if(name == "NL"){
    lin= FL_lin; lin_pos = FL_lin_pos; quad = FL_quad; quadfit= FL_quadfit;
  }
  if(name == "NR"){
    lin= FR_lin; lin_pos = FR_lin_pos; quad = FR_quad; quadfit= FR_quadfit;
  }
  if(name == "Wjets"){
    lin= Wjets_lin; lin_pos = Wjets_lin_pos; quad = Wjets_quad; quadfit= Wjets_quadfit;
  }
  if(name == "QCD"){
    lin= QCD_lin; lin_pos = QCD_lin_pos; quad = QCD_quad; quadfit= QCD_quadfit;
  }
  if(name == "RemBkg"){
    lin= RemBkg_lin; lin_pos = RemBkg_lin_pos; quad = RemBkg_quad; quadfit= RemBkg_quadfit;
  }

  double x[fNBins];
  //double x[15];  

  //Fill different chi^2 values
  for(int i = 1; i <= fNBins; ++i){
    val_lin[i-1]     = lin[i-1];
    val_lin_pos[i-1] = lin_pos[i-1];
    val_quad[i-1]    = quad[i-1];
    val_quadfit[i-1] = quadfit[i-1];
    x[i-1] = -1.0+(i-1)*2.0/15;
  }
  
  TCanvas *c0        = new TCanvas("", "", 720, 1000);

  TPad *pad1 = new TPad("pad1","pad1", 0.0, 0.800, 1.0, 1.000);
  TPad *pad2 = new TPad("pad2","pad2", 0.0, 0.600, 1.0, 0.800);
  TPad *pad3 = new TPad("pad3","pad3", 0.0, 0.400, 1.0, 0.600);
  TPad *pad4 = new TPad("pad4","pad4", 0.0, 0.175, 1.0, 0.400);

  pad1->SetBottomMargin(0.010);
  pad1->SetBorderMode(0);
  pad2->SetBottomMargin(0.010);
  pad2->SetBorderMode(0);
  pad3->SetBottomMargin(0.010);
  pad3->SetBorderMode(0);
  pad4->SetBottomMargin(0.200);

  pad1->SetTicks(1,1);
  pad2->SetTicks(1,1);
  pad3->SetTicks(1,1);
  pad4->SetTicks(1,1);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();

  
  pad1 -> cd();

  TGraph *gr   = new TGraph(fNBins, x, val_lin);
  //  gr->SetTitle(name.c_str());
  gr->SetMarkerColor(kCyan);
  gr->SetLineColor(kCyan);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(20);
  gr -> GetYaxis() -> SetTitle("Chi2: linear");
  gr -> GetXaxis() -> SetLabelSize(0.095);
  gr -> GetYaxis() -> SetLabelSize(0.095);
  gr -> GetXaxis() -> SetTitleSize(0.08);
  gr -> GetYaxis() -> SetTitleSize(0.08);
  gr -> GetYaxis() -> SetTitleOffset(0.40);
  gr -> GetYaxis() -> CenterTitle();
  gr->Draw("AP");

  pad2 -> cd();

  TGraph *gr2   = new TGraph(fNBins, x, val_quad);
  //  gr2->SetTitle(name.c_str());
  gr2->SetMarkerColor(kViolet);
  gr2->SetLineColor(kViolet);
  gr2->SetLineWidth(2);
  gr2->SetMarkerStyle(20);
  gr2->GetYaxis()->SetTitle("Chi2: quad. int.");
  gr2 -> GetXaxis() -> SetLabelSize(0.095);
  gr2 -> GetYaxis() -> SetLabelSize(0.095);
  gr2 -> GetXaxis() -> SetTitleSize(0.08);
  gr2 -> GetYaxis() -> SetTitleSize(0.08);
  gr2 -> GetYaxis() -> SetTitleOffset(0.40);
  gr2 -> GetYaxis() -> CenterTitle();

  gr2->Draw("AP");

  pad3 -> cd();

  TGraph *gr3   = new TGraph(fNBins, x, val_lin_pos);
  // gr3->SetTitle(name.c_str());
  gr3->SetMarkerColor(kGreen);
  gr3->SetLineColor(kGreen);
  gr3->SetLineWidth(2);
  gr3->SetMarkerStyle(20);
  gr3->GetYaxis()->SetTitle("Chi2: PW linear");
  gr3 -> GetXaxis() -> SetLabelSize(0.095);
  gr3 -> GetYaxis() -> SetLabelSize(0.095);
  gr3 -> GetXaxis() -> SetTitleSize(0.08);
  gr3 -> GetYaxis() -> SetTitleSize(0.08);
  gr3 -> GetYaxis() -> SetTitleOffset(0.40);
  gr3 -> GetYaxis() -> CenterTitle();
  gr3->Draw("AP");

  pad4 -> cd();

  TGraph *gr4   = new TGraph(fNBins, x, val_quadfit);
  // gr4->SetTitle(name.c_str());
  gr4->SetMarkerColor(kRed);
  gr4->SetLineColor(kRed);
  gr4->SetLineWidth(2);
  gr4->SetMarkerStyle(20);
  gr4->GetXaxis()->SetTitle("cos #theta*");
  gr4->GetYaxis()->SetTitle("Chi2: quad. fit");
  gr4 -> GetXaxis() -> SetLabelSize(0.085);
  gr4 -> GetYaxis() -> SetLabelSize(0.085);
  gr4 -> GetXaxis() -> SetTitleSize(0.11);
  gr4 -> GetYaxis() -> SetTitleSize(0.08);
  gr4 -> GetXaxis() -> SetTitleOffset(1.10);
  gr4 -> GetYaxis() -> SetTitleOffset(0.40);
  gr4 -> GetYaxis() -> CenterTitle();

  gr4->Draw("AP");

  c0->Draw();
  c0->Print((fOutputFolder+"/Chi2comparison_"+name+".png").c_str());
    

}


//Fill information of histograms in corresponding 3D-double vector
void PlotInterpolationCurves::FillGraphs(std::vector<std::vector<TH1D> > hist, int iHist)
{
  
  for(int iK = 0; iK < hist.size(); ++iK){

    for(int iBin = 0; iBin < hist[iK][iHist].GetNbinsX(); ++iBin){
      
      //fHistEntries defined in MakePlots, protected 3D-double vectors
      fHistEntries[iHist][iBin][iK]    = hist[iK][iHist].GetBinContent(iBin+1);
      if(iK == 0)
	fHistEntriesErr[iHist][iBin][iK] = 0.0001;
      else
	fHistEntriesErr[iHist][iBin][iK] = hist[iK][iHist].GetBinError(iBin+1);

    }
    
  }
  
}



//Plot different interpolation curves in one canvas, calculate and fill fFitFunc
void PlotInterpolationCurves::PrintGraphs(std::string name, int iHist, std::string fInterpolationMethod)
{

  SetAtlasStyle();

  std::vector<TF1> HelpFit, HelpFit_QuadFit;
  HelpFit.clear();
  HelpFit_QuadFit.clear();
  
  for(int iBin = 0; iBin < fNBins; ++iBin){

    int   n = fHistEntries[iHist][iBin].size();
    double x[n], xerr[n], y[n], yerr[n];

    for(int k = 0; k < n; ++k){
      
      y[k]    = fHistEntries[iHist][iBin][k];
      yerr[k] = fHistEntriesErr[iHist][iBin][k];

      x[k]    = fKValues[k];
      xerr[k] = 0.0; 

    }


    double add_term     =  y[0]; // value for k = 0
    double lin_term     = (y[1] - y[2])/2.0; // 1 -> up, 2 -> down 
    double lin_term_pos =  y[1] - y[0];
    double lin_term_neg =  y[0] - y[2];
    double quad_term    = (y[1] + y[2])/2.0 - y[0];
    double add_term_pos = y[1]; //hist_vec_F0[2*k+2]->GetBinContent(i);
    double add_term_neg = y[2]; //hist_vec_F0[2*k+1]->GetBinContent(i);
    double quad_term_pos = 1.5*y[1] + 0.5*y[2] - 2*y[0]; 
    double quad_term_neg = -1.5*y[2] - 0.5*y[1] + 2*y[0];  

    TF1 *fquadfit  = new TF1 ("fquadfit",  "[0] + x*[1] + x*x*[2]", -4.0,  4.0);
    TF1 *flin      = new TF1 ("flin",      "[0] + x*[1]",           -4.0,  4.0);
    TF1 *fquad     = new TF1 ("fquad",     "[0] + x*[1] + x*x*[2]", -1.0,  1.0);
    TF1 *flin_pos  = new TF1 ("flin_pos",  "[0] + x*[1]",            0.0,  4.0);
    TF1 *flin_neg  = new TF1 ("flin_neg",  "[0] + x*[1]",           -4.0,  0.0);
    TF1 *fquad_pos = new TF1 ("fquad_pos", "[0] + (x-1)*[1]",        1.0,  4.0);
    TF1 *fquad_neg = new TF1 ("fquad_neg", "[0] + (x+1)*[1]",       -4.0, -1.0);

    TF1 *flin_comb = new TF1 ("flin_comb"," (x>0)*([0] + x*[1]) + (x<0)*([0] + x*[2]) ", -4.0, 4.0);
    TF1 *fquad_comb= new TF1 ("fquad_comb"," (x>1)*([0] + (x-1)*[1]) + (TMath::Abs(x)<1)*([2] + x*[3] + x*x*[4]) + (x<(-1))*([5] + (x+1)*[6]) ", -4.0, 4.0);


    flin      -> SetParameter(0,add_term);
    flin      -> SetParameter(1,lin_term);
    fquad     -> SetParameter(0,add_term);
    fquad     -> SetParameter(1,lin_term);
    fquad     -> SetParameter(2,quad_term);
    flin_pos  -> SetParameter(0,add_term);
    flin_pos  -> SetParameter(1,lin_term_pos);
    flin_neg  -> SetParameter(0,add_term);
    flin_neg  -> SetParameter(1,lin_term_neg);
    fquad_pos -> SetParameter(0,add_term_pos);
    fquad_pos -> SetParameter(1,quad_term_pos);
    fquad_neg -> SetParameter(0,add_term_neg);
    fquad_neg -> SetParameter(1,quad_term_neg);

    flin_comb -> SetParameter(0,add_term);
    flin_comb -> SetParameter(1,lin_term_pos);
    flin_comb -> SetParameter(2,lin_term_neg);

    fquad_comb-> SetParameter(0,add_term_pos);
    fquad_comb-> SetParameter(1,quad_term_pos);
    fquad_comb-> SetParameter(2,add_term);
    fquad_comb-> SetParameter(3,lin_term);
    fquad_comb-> SetParameter(4,quad_term);
    fquad_comb-> SetParameter(5,add_term_neg);
    fquad_comb-> SetParameter(6,quad_term_neg);

    std::stringstream nr_bin;
    nr_bin << iBin + 1;

    std::string outputFile     = fOutputFolder+"/Interpolation_"+name+"_bin"+nr_bin.str()+".eps";
    std::string outputFile_png = fOutputFolder+"/Interpolation_"+name+"_bin"+nr_bin.str()+".png";

    TCanvas *c0        = new TCanvas("", "", 720, 1000);

    TPad *pad1 = new TPad("pad1","pad1", 0.0, 0.500, 1.0, 1.000);
    TPad *pad2 = new TPad("pad2","pad2", 0.0, 0.400, 1.0, 0.500);
    TPad *pad3 = new TPad("pad3","pad3", 0.0, 0.300, 1.0, 0.400);
    TPad *pad4 = new TPad("pad4","pad4", 0.0, 0.200, 1.0, 0.300);
    TPad *pad5 = new TPad("pad5","pad5", 0.0, 0.000, 1.0, 0.200);
    
    pad1->SetBottomMargin(0.010);
    pad1->SetBorderMode(0);
    pad2->SetBottomMargin(0.040);
    pad2->SetBorderMode(0);
    pad3->SetBottomMargin(0.040);
    pad3->SetBorderMode(0);
    pad4->SetBottomMargin(0.040);
    pad4->SetBorderMode(0);
    pad5->SetBottomMargin(0.500);
    // pad5->SetBorderMode(0);

    pad1->SetTicks(1,1);
    pad2->SetTicks(1,1);
    pad3->SetTicks(1,1);       
    pad4->SetTicks(1,1);
    pad5->SetTicks(1,1);
    
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();
    pad5->Draw();

    pad1->cd();

    TGraphErrors *gr   = new TGraphErrors(n,x,y,xerr,yerr);
    gr->SetTitle(name.c_str());
    gr->SetMarkerColor(4);
    gr->SetLineColor(4);
    gr->SetLineWidth(2);
    gr->SetMarkerStyle(20);
    gr->Draw("AP");

    gr->Fit("fquadfit", "QN");


    //Fill HelpFit according to flag:
    //QuadraticFit,QuadraticInterp,LinearInterp,PiecewiseLinear
    if(fInterpolationMethod == "QuadraticFit")       HelpFit.push_back(*fquadfit);
    else if (fInterpolationMethod == "LinearInterp")     HelpFit.push_back(*flin);
    else if (fInterpolationMethod == "PiecewiseLinear") HelpFit.push_back(*flin_comb);
    else if (fInterpolationMethod == "QuadraticInterp") HelpFit.push_back(*fquad_comb);
    else {
      WriteErrorStatus("PlotInterpolation", "No such option for the interpolation mode!");
      WriteErrorStatus("PlotInterpolation", "EXIT");
      exit(-1);;
    }

    //Fill always QuadFit function in correspondong vector
    HelpFit_QuadFit.push_back(*fquadfit);

    gr -> GetXaxis() -> SetTitle(("k("+name+", Bin"+nr_bin.str()+")").c_str());
    gr -> GetYaxis() -> SetTitle(("#scale[0.9]{Entries "+name+"}").c_str() );
    gr -> GetYaxis() -> SetTitleOffset(1.42);

    gr -> GetYaxis() -> SetLabelSize(0.038);

    TLatex l2a;
    l2a.SetTextAlign(9);
    l2a.SetTextFont(72);
    l2a.SetTextSize(0.04);
    l2a.SetNDC();
    l2a.DrawLatex(0.19, 0.880, "ATLAS");
    TLatex l3a;
    l3a.SetTextAlign(9);
    l3a.SetTextSize(0.04);
    l3a.SetNDC();
    l3a.DrawLatex(0.29, 0.880, "Work in progress");
    TLatex l4a;
    l4a.SetTextAlign(9);
    l4a.SetTextSize(0.041);
    l4a.SetNDC();
    l4a.DrawLatex(0.19, 0.82, "Interpolation studies per bin");

    // auch geklaut von Philipp

    flin->SetLineColor(kCyan);
    flin->SetLineWidth(2);
    flin->SetLineStyle(2);
    flin->Draw("SAME");

    fquad->SetLineColor(kViolet);
    fquad->SetLineWidth(2);
    fquad->SetLineStyle(5);
    fquad->Draw("SAME");

    fquad_pos->SetLineColor(kViolet);
    fquad_pos->SetLineWidth(2);
    fquad_pos->SetLineStyle(5);
    fquad_pos->Draw("same");

    fquad_neg->SetLineColor(kViolet);
    fquad_neg->SetLineWidth(2);
    fquad_neg->SetLineStyle(5);
    fquad_neg->Draw("same");

    flin_pos->SetLineColor(kGreen);
    flin_pos->SetLineWidth(2);
    flin_pos->SetLineStyle(4);
    flin_pos->Draw("SAME");

    flin_neg->SetLineColor(kGreen);
    flin_neg->SetLineWidth(2);
    flin_neg->SetLineStyle(4);
    flin_neg->Draw("SAME");

    fquadfit->SetLineColor(kRed);
    fquadfit->SetLineWidth(2);
    fquadfit->SetLineStyle(1);
    fquadfit->Draw("SAME");

    flin_comb->SetLineColor(kBlue+3);
    flin_comb->SetLineWidth(2);
    //flin_comb->Draw("SAME");

    fquad_comb->SetLineColor(kOrange);
    fquad_comb->SetLineWidth(2);
    //fquad_comb->Draw("SAME");

    TLegend *legend = new TLegend(0.60,0.05,0.9,0.32);
    //legend->AddEntry(hsig,"Data, #sqrt{s} = 7 TeV","lep");
    legend->SetTextFont(42);
    legend->AddEntry(gr,       "#scale[1.85]{Template values}","lep");       
    legend->AddEntry(flin,     "#scale[1.85]{Linear interpolation}",    "l");
    legend->AddEntry(fquad,    "#scale[1.85]{Quadratic interpolation}", "l");
    legend->AddEntry(flin_pos, "#scale[1.85]{Piecewise lin. interp.}",  "l");
    legend->AddEntry(fquadfit, "#scale[1.85]{Quadratic fit}",           "l");
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->Draw("SAME");

    TF1 *norm1 = new TF1("fa1","1", -4.0, 4.0);
    norm1 -> SetLineColor(kBlack);
    norm1 -> SetLineStyle(1);
    norm1 -> SetLineWidth(2);

    pad2 -> cd();

    TGraphErrors *comp_lin     = GetRatioDistribution(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *flin);
    TGraphErrors *comp_quad    = GetRatioDistribution(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *fquad);
    TGraphErrors *comp_lin_pos = GetRatioDistribution(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *flin_pos, *flin_neg);
    TGraphErrors *comp_quadfit = GetRatioDistribution(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *fquadfit);

    double chi2_comp_lin     = GetChi2Values(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *flin);
    double chi2_comp_quad    = GetChi2Values(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *fquad);
    double chi2_comp_lin_pos = GetChi2Values(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *flin_pos, *flin_neg);
    double chi2_comp_quadfit = GetChi2Values(fKValues, fHistEntries[iHist][iBin], fHistEntriesErr[iHist][iBin], *fquadfit);

    if(chi2_comp_lin     < 10e-4) chi2_comp_lin     = 0;
    if(chi2_comp_quad    < 10e-4) chi2_comp_quad    = 0;
    if(chi2_comp_lin_pos < 10e-4) chi2_comp_lin_pos = 0;
    if(chi2_comp_quadfit < 10e-4) chi2_comp_quadfit = 0;

    std::stringstream oss1,oss2,oss3,oss4;
    oss1 << setprecision(3) << chi2_comp_lin;
    oss2 << setprecision(3) << chi2_comp_quad;
    oss3 << setprecision(3) << chi2_comp_lin_pos;
    oss4 << setprecision(3) << chi2_comp_quadfit;

    if(name == "N0"){
      F0_lin.push_back(chi2_comp_lin);
      F0_quad.push_back(chi2_comp_quad);
      F0_lin_pos.push_back(chi2_comp_lin_pos);
      F0_quadfit.push_back(chi2_comp_quadfit);
    }
    if(name == "NL"){
      FL_lin.push_back(chi2_comp_lin);
      FL_quad.push_back(chi2_comp_quad);
      FL_lin_pos.push_back(chi2_comp_lin_pos);
      FL_quadfit.push_back(chi2_comp_quadfit);
    }
    if(name == "NR"){
      FR_lin.push_back(chi2_comp_lin);
      FR_quad.push_back(chi2_comp_quad);
      FR_lin_pos.push_back(chi2_comp_lin_pos);
      FR_quadfit.push_back(chi2_comp_quadfit);
    }
    if(name == "Wjets"){
      Wjets_lin.push_back(chi2_comp_lin);
      Wjets_quad.push_back(chi2_comp_quad);
      Wjets_lin_pos.push_back(chi2_comp_lin_pos);
      Wjets_quadfit.push_back(chi2_comp_quadfit);
    }
    if(name == "QCD"){
      QCD_lin.push_back(chi2_comp_lin);
      QCD_quad.push_back(chi2_comp_quad);
      QCD_lin_pos.push_back(chi2_comp_lin_pos);
      QCD_quadfit.push_back(chi2_comp_quadfit);
    }
    if(name == "RemBkg"){
      RemBkg_lin.push_back(chi2_comp_lin);
      RemBkg_quad.push_back(chi2_comp_quad);
      RemBkg_lin_pos.push_back(chi2_comp_lin_pos);
      RemBkg_quadfit.push_back(chi2_comp_quadfit);
    }



    //    std::cout << chi2_comp_lin << "\t" << chi2_comp_quadfit << std::endl;

    double max_help[4], min_help[4];

    max_help[0] = TMath::MaxElement(fKValues.size(), comp_lin     -> GetY());  
    min_help[0] = TMath::MinElement(fKValues.size(), comp_lin     -> GetY());
    max_help[1] = TMath::MaxElement(fKValues.size(), comp_quad    -> GetY());
    min_help[1] = TMath::MinElement(fKValues.size(), comp_quad    -> GetY());
    max_help[2] = TMath::MaxElement(fKValues.size(), comp_lin_pos -> GetY());
    min_help[2] = TMath::MinElement(fKValues.size(), comp_lin_pos -> GetY());
    max_help[3] = TMath::MaxElement(fKValues.size(), comp_quadfit -> GetY());
    min_help[3] = TMath::MinElement(fKValues.size(), comp_quadfit -> GetY());

    double max = TMath::MaxElement(4, max_help) + 0.025;
    double min = TMath::MinElement(4, min_help) - 0.025;

    comp_lin -> SetMarkerColor(kCyan);
    comp_lin -> SetLineColor(kCyan);
    comp_lin -> SetLineWidth(2);
    comp_lin -> SetMarkerStyle(20);
    comp_lin -> SetMarkerSize(1.2);
    comp_lin -> GetXaxis() -> SetLabelFont(42);
    comp_lin -> GetXaxis() -> SetLabelSize(0.15);
    comp_lin -> GetYaxis() -> SetLabelSize(0.18);
    comp_lin -> GetYaxis() -> SetLabelFont(42);
    comp_lin -> GetXaxis() -> SetLabelOffset(0.03);
    comp_lin -> GetYaxis() -> SetTitleFont(42);
    comp_lin -> GetXaxis() -> SetTitleFont(42);
    comp_lin -> GetYaxis() -> SetTitleSize(0.20);
    comp_lin -> GetYaxis() -> SetTitleOffset(0.35);
    comp_lin -> GetYaxis() -> SetTitle("Tem./Lin.");
    comp_lin -> GetYaxis() -> CenterTitle();
    comp_lin -> GetYaxis() -> SetNdivisions(505);
    comp_lin -> SetMaximum(max);
    comp_lin -> SetMinimum(min);
    comp_lin -> Draw("AP");

    std::string label1 = "#chi^{2}/ndf = "+oss1.str();
    TLatex l1;
    l1.SetTextAlign(9);
    l1.SetTextSize(0.17);
    l1.SetNDC();
    l1.DrawLatex(0.77, 0.710, label1.c_str());

    norm1 -> Draw("SAME");
    
    pad3 -> cd();

    comp_quad -> SetMarkerColor(kViolet);
    comp_quad -> SetLineColor(kViolet);
    comp_quad -> SetLineWidth(2);
    comp_quad -> SetMarkerStyle(20);
    comp_quad -> SetMarkerSize(1.2);
    comp_quad -> GetXaxis() -> SetLabelFont(42);
    comp_quad -> GetXaxis() -> SetLabelSize(0.15);
    comp_quad -> GetYaxis() -> SetLabelSize(0.18);
    comp_quad -> GetYaxis() -> SetLabelFont(42);
    comp_quad -> GetXaxis() -> SetLabelOffset(0.03);
    comp_quad -> GetYaxis() -> SetTitleFont(42);
    comp_quad -> GetXaxis() -> SetTitleFont(42);
    comp_quad -> GetYaxis() -> SetTitleSize(0.20);
    comp_quad -> GetYaxis() -> SetTitleOffset(0.35);
    comp_quad -> GetYaxis() -> SetTitle("Tem./Quad.");
    comp_quad -> GetYaxis() -> CenterTitle();
    comp_quad -> GetYaxis() -> SetNdivisions(505);
    comp_quad -> SetMaximum(max);
    comp_quad -> SetMinimum(min);
    comp_quad -> Draw("AP");
  
    std::string label2 = "#chi^{2}/ndf = "+oss2.str();
    TLatex l2;
    l2.SetTextAlign(9);
    l2.SetTextSize(0.17);
    l2.SetNDC();
    l2.DrawLatex(0.77, 0.710, label2.c_str());

    norm1 -> Draw("SAME");

    pad4 -> cd();

    comp_lin_pos -> SetMarkerColor(kGreen);
    comp_lin_pos -> SetLineColor(kGreen);
    comp_lin_pos -> SetLineWidth(2);
    comp_lin_pos -> SetMarkerStyle(20);
    comp_lin_pos -> SetMarkerSize(1.2);
    comp_lin_pos -> GetXaxis() -> SetLabelFont(42);
    comp_lin_pos -> GetXaxis() -> SetLabelSize(0.15);
    comp_lin_pos -> GetYaxis() -> SetLabelSize(0.18);
    comp_lin_pos -> GetYaxis() -> SetLabelFont(42);
    comp_lin_pos -> GetXaxis() -> SetLabelOffset(0.03);
    comp_lin_pos -> GetYaxis() -> SetTitleFont(42);
    comp_lin_pos -> GetXaxis() -> SetTitleFont(42);
    comp_lin_pos -> GetYaxis() -> SetTitleSize(0.20);
    comp_lin_pos -> GetYaxis() -> SetTitleOffset(0.35);
    comp_lin_pos -> GetYaxis() -> SetTitle("Tem./PW");
    comp_lin_pos -> GetYaxis() -> CenterTitle();
    comp_lin_pos -> GetYaxis() -> SetNdivisions(505);
    comp_lin_pos -> SetMaximum(max);
    comp_lin_pos -> SetMinimum(min);
    comp_lin_pos -> Draw("AP");

    std::string label3 = "#chi^{2}/ndf = "+oss3.str();
    TLatex l3;
    l3.SetTextAlign(9);
    l3.SetTextSize(0.17);
    l3.SetNDC();
    l3.DrawLatex(0.77, 0.710, label3.c_str());

    norm1 -> Draw("SAME");

    pad5 -> cd();
  
    comp_quadfit -> SetMarkerColor(kRed);
    comp_quadfit -> SetLineColor(kRed);
    comp_quadfit -> SetLineWidth(2);
    comp_quadfit -> SetMarkerStyle(20);
    comp_quadfit -> SetMarkerSize(1.2);
    comp_quadfit -> GetXaxis() -> SetLabelFont(42);
    comp_quadfit -> GetXaxis() -> SetLabelSize(0.10);
    comp_quadfit -> GetYaxis() -> SetLabelSize(0.09);
    comp_quadfit -> GetYaxis() -> SetLabelFont(42);
    comp_quadfit -> GetXaxis() -> SetLabelOffset(0.03);
    comp_quadfit -> GetYaxis() -> SetTitleFont(42);
    comp_quadfit -> GetXaxis() -> SetTitleFont(42);
    comp_quadfit -> GetYaxis() -> SetTitleSize(0.10);
    comp_quadfit -> GetYaxis() -> SetTitleOffset(0.70);
    comp_quadfit -> GetYaxis() -> SetTitle("Tem./Fit");
    comp_quadfit -> GetYaxis() -> CenterTitle();
    comp_quadfit -> GetYaxis() -> SetNdivisions(505);
    comp_quadfit -> SetMaximum(max);
    comp_quadfit -> SetMinimum(min);
    comp_quadfit -> Draw("AP");

    std::string label4 = "#chi^{2}/ndf = "+oss4.str();
    TLatex l4;
    l4.SetTextAlign(9);
    l4.SetTextSize(0.0875);
    l4.SetNDC();
    l4.DrawLatex(0.77, 0.795, label4.c_str());

    norm1 -> Draw("SAME");

    comp_quadfit -> GetXaxis() -> SetTitle(("k("+name+", Bin"+nr_bin.str()+")").c_str());
    comp_quadfit -> GetXaxis() -> SetTitleSize(0.10);
    comp_quadfit -> GetXaxis() -> SetTitleOffset(1.35);

    c0 -> Draw();
    c0 -> Print((outputFile).c_str());
    c0 -> Print((outputFile_png).c_str());

    std::string new_outputFile_png = outputFile_png.replace(outputFile_png.find(fOutputFolder), strlen(fOutputFolder.c_str()), "");

    new_outputFile_png = new_outputFile_png.replace(new_outputFile_png.find("/"), strlen("/"), "");

    std::ofstream* htmlfile = new std::ofstream();

    htmlfile->open((_htmls[fHTMLLabel]).c_str(),fstream::app);
    *htmlfile << "<tr>"
		<< "<td>"<< name <<"</td>";

    *htmlfile << "<td>";
    *htmlfile << "<p><a href=\"" << new_outputFile_png.c_str() << "\">";
    *htmlfile << "<img src=\""   << new_outputFile_png.c_str() << "\">";
    *htmlfile << "</td>" << std::endl;
   
    htmlfile->close();
    
    delete flin;
    delete fquad;
    delete flin_pos;
    delete flin_neg;
    delete fquadfit;
    delete gr;
    delete c0; 
 
  }

  fFitFunc.push_back(HelpFit);
  fFitFunc_QuadFit.push_back(HelpFit_QuadFit);

}



//Draw plot which contains all different templates (including variations) for one certain systematic effect
void PlotInterpolationCurves::DrawFinalPlot(std::vector<TH1D > histo0, std::vector<std::string> Label, std::string outputFile, bool ResidualFlag, bool norm)
{
  SetAtlasStyle();

  int NTotalBins = histo0[0].GetNbinsX();

  std::cout << NTotalBins << std::endl;

  std::vector<TH1D*> histo;

  //Convert TH1D into TH1D*
  for (int i = 0; i != histo0.size(); i++) {
    
    histo0[i].Sumw2();

    if(norm)    
      histo0[i].Scale(1.0/histo0[i].Integral());
    
    TH1D * hist_dummy = &histo0[i];
    histo.push_back( hist_dummy);
    
  }

  TCanvas *c0;

  TPad *pad1;
  TPad *pad2;
  TPad *pad3;

  if(!ResidualFlag){
    c0   = new TCanvas("", "", 720, 720);
    pad1 = new TPad("pad1","pad1", 0.0, 0.0, 1.0, 1.0);
    pad2 = new TPad("pad2","pad2", 0.0, 0.0, 1.0, 0.0);
  }
  else{
    c0   = new TCanvas("", "", 720, 900);
    pad1 = new TPad("pad1","pad1", 0.0, 0.400, 1.0, 1.000);
    pad2 = new TPad("pad2","pad2", 0.0, 0.250, 1.0, 0.400);
    pad3 = new TPad("pad3","pad3", 0.0, 0.000, 1.0, 0.250);

    pad1->SetBottomMargin(0.010);
    pad1->SetBorderMode(0);
    pad2->SetBottomMargin(0.010);
    pad2->SetBorderMode(0);
    pad3->SetBottomMargin(0.400);
    pad3->SetBorderMode(0);

    pad1->SetTicks(1,1);
    pad2->SetTicks(1,1);
    pad3->SetTicks(1,1);

    pad1->Draw();
    pad2->Draw();
    pad3->Draw();

  }


  double nmax = histo[1] -> GetBinContent(histo[1] -> GetMaximumBin())*1.55;

  int nbins = histo[1] -> GetNbinsX();

  double lower_edge  = histo[0] -> GetBinLowEdge(1);
  double bin_width   = histo[0] -> GetBinWidth(1);
  double number_bins = histo[0] -> GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;


  //  fTotalSum1 = new TH1D("Sum1", "", nbins, lower_edge, upper_edge);

  /*  fTotalSum4 -> Add(fTTbarSum4);
  fTotalSum4 -> Add(fWjetsSum);
  fTotalSum4 -> Add(fZjetsSum);
  fTotalSum4 -> Add(fDibosonSum);
  fTotalSum4 -> Add(fSingleTopSum);
  fTotalSum4 -> Add(fQCDSum); */


  pad1 -> cd();

  for(int iHist = 0; iHist < histo.size(); ++iHist){

    //histo[iHist] = normalise(histo[iHist]);   //only for shape plots!!!
    histo[iHist] -> SetFillColor(kWhite);
    histo[iHist] -> SetLineWidth(2.25);
    //histo[iHist] -> SetLineStyle(iHist+1);
  
    if(iHist==2){
      histo[iHist] -> SetLineColor(kCyan+2);
      histo[iHist] -> SetLineStyle(2);
    }
    if(iHist==1){
      histo[iHist] -> SetLineColor(kOrange+8);
      histo[iHist] -> SetLineStyle(3);
    }
    //    if(iHist!=4) histo[iHist] -> SetLineColor(iHist+1);
    //    else histo[iHist] -> SetLineColor(kViolet+5);
    histo[iHist] -> SetMarkerSize(0);
    if(iHist == 0)
      histo[iHist] -> Draw();
    else
      histo[iHist] -> Draw("SAME");

  }

  TLine line1 = TLine(0.0, histo[0] -> GetMinimum(), 0.0, nmax/2.0);
  line1.SetLineStyle(2);
  
  if(NTotalBins > 15)
    line1.Draw();
    
  //nmax = histo[1] -> GetBinContent(histo[1] -> GetMaximumBin())*1.2; //only for shape plots!!!

  int counter_zero = 0;

  for(int ibin = 1; ibin <= nbins; ++ibin)
    if(histo[0]  -> GetBinContent(ibin) == 0)
      counter_zero++;

  int ngoodbins = nbins - counter_zero;

  TGraphErrors *comp1 = new TGraphErrors();
  TGraphErrors *comp2 = new TGraphErrors();

  comp1 = GetRatioDistribution(histo[0], histo[1], ngoodbins);
  if(histo.size() > 2)
    comp2 = GetRatioDistribution(histo[0], histo[2], ngoodbins);

  comp1->GetXaxis()->SetLimits(lower_edge, upper_edge); 
  comp2->GetXaxis()->SetLimits(lower_edge, upper_edge);

  double max = 1.20;
  double min = 0.80;

  if(norm){

    max = 1.05;
    min = 0.95;
    
  }

  std::stringstream width;
  if(bin_width > 9.5){
    width << std::setprecision(2) << bin_width;
  }
  else
    width << std::setprecision(1) << bin_width;

  //  histo[0] -> GetXaxis() -> SetTitle("#scale[3.5]{cos #theta*}");

  if(norm)
    histo[0] -> GetYaxis() -> SetTitle("#scale[1.15]{Normalised}");
  else
    histo[0] -> GetYaxis() -> SetTitle("#scale[1.15]{Events}");

  histo[0]->GetXaxis()->SetLabelSize(0.04);
  histo[0]->GetYaxis()->SetLabelSize(0.05);
  if(!norm)
    histo[0]->SetMinimum(0.001);
  histo[0] -> GetYaxis() -> SetTitleOffset(1.5);
  //  histo[0] -> GetXaxis() -> SetTitleOffset(1.2);

  //define the legend...
  TLegend *fLegend = new TLegend(0.60, 0.65, 0.80, 0.90);
  for(unsigned int k = 0; k < histo.size(); ++k)
  if(norm)  
    fLegend -> AddEntry(histo[k],  Label[k].c_str(), "l");
  else
    fLegend -> AddEntry(histo[k],  (Label[k]+" ("+ std::to_string((int)histo[k]->Integral())+")").c_str(), "l");
  fLegend  -> SetFillColor(0);
  fLegend  -> SetLineColor(0);
  fLegend  -> SetBorderSize(0);
  fLegend  -> SetTextFont(72);
  fLegend  -> SetTextSize(0.035);
  fLegend  -> Draw();

  histo[0]  -> SetMaximum(nmax);
  //histo[0]  -> SetMinimum(0.0);

  //set labels....
  TLatex l1;
  l1.SetTextAlign(9);
  l1.SetTextSize(0.04);
  l1.SetNDC();
  l1.DrawLatex(0.21, 0.750, fLumiLabel.c_str());
  TLatex l2;
  l2.SetTextAlign(9);
  l2.SetTextFont(72);
  l2.SetTextSize(0.04);
  l2.SetNDC();
  l2.DrawLatex(0.19, 0.86, "ATLAS");
  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextSize(0.04);
  l3.SetNDC();
  l3.DrawLatex(0.295, 0.86, "Work in progress");
  //l3.DrawLatex(0.19, 0.81, "Templates used for interpolation");
  /*  TLatex l4;
  l4.SetTextAlign(9);
  l4.SetTextSize(0.042);
  l4.SetTextFont(72);
  l4.SetNDC();
  l4.DrawLatex(0.66, 0.85, (fBtagLabel+" "+fJetBinLabel).c_str()); */

  int NrSubPlots = fLeptonLabel.size();

  float SubPlotWidth     = 0.85/float(NrSubPlots);
  float SubPlotHalfWidth = SubPlotWidth/2.0;

  for(int iPlot = 0; iPlot < NrSubPlots; ++iPlot){

    //    std::cout << iPlot << "\t" << fLeptonLabel[iPlot].c_str() << std::endl;
    
    float SubPlotPosition = (1+2*iPlot)*SubPlotHalfWidth;

    if(iPlot == 0) SubPlotPosition += 0.10;
    if(iPlot == 1) SubPlotPosition += 0.05;

    std::cout << iPlot << "\t" << SubPlotWidth << "\t" << SubPlotHalfWidth << "\t" << SubPlotPosition << std::endl;
    
    TLatex l3;
    l3.SetTextAlign(9);
    l3.SetTextFont(42);
    l3.SetTextSize(0.037);
    l3.SetNDC();
    l3.DrawLatex(SubPlotPosition, 0.120, (fLeptonLabel[iPlot]).c_str());
    l3.DrawLatex(SubPlotPosition, 0.080, (fJetBinLabel[iPlot]).c_str());
    l3.DrawLatex(SubPlotPosition, 0.040, (fBTagLabel[iPlot]).c_str());

  }


  pad2->cd();

  TF1 *norm1 = new TF1("fa1","1", lower_edge, upper_edge);
  norm1 -> SetLineColor(kBlack);
  norm1 -> SetLineStyle(1);
  norm1 -> SetLineWidth(2);
  
  comp1 -> SetMarkerColor(kOrange+8);
  comp1 -> SetMarkerStyle(20);
  comp1 -> SetMarkerSize(1.2);
  comp1 -> SetLineColor(kOrange+8);
  comp1 -> SetLineWidth(2);
  comp1 -> GetXaxis() -> SetLabelFont(42);
  comp1 -> GetXaxis() -> SetLabelSize(0.165);
  comp1 -> GetYaxis() -> SetLabelSize(0.165);
  comp1 -> GetYaxis() -> SetLabelFont(42);
  comp1 -> GetXaxis() -> SetLabelOffset(0.03);
  comp1 -> GetYaxis() -> SetTitleFont(42);
  comp1 -> GetXaxis() -> SetTitleFont(42);
  comp1 -> GetYaxis() -> SetTitleSize(0.17);
  comp1 -> GetYaxis() -> SetTitleOffset(0.37);
  //  comp1 -> GetYaxis() -> SetTitle("1#sigma-Up/Nom.");
  comp1 -> GetYaxis() -> SetTitle("Up/Nom.");
  comp1 -> GetXaxis() -> SetNdivisions(505);
  comp1 -> GetYaxis() -> SetNdivisions(505);
  comp1 -> SetMaximum(max);
  comp1 -> SetMinimum(min);
  //------Take me out after check
  comp1 -> SetMaximum(1.03);
  comp1 -> SetMinimum(0.97);

  comp1 -> Draw("AP");
  norm1 -> Draw("Same");

  //  if(NTotalBins > 15){

    //std::cout << "Hier!!!   " << NTotalBins << "\t" << min << "\t" << max << std::endl;

    TLine line = TLine(0.0, 0.8, 0.0, 1.2);
    line.SetLineStyle(2);
    // line.Set

    if(NTotalBins > 15)
      line.Draw();

    //  }

  pad3->cd();

  //  int NBins = comp2 -> GetNbinsX();

  int counter = 0;

  //  comp2 -> GetXaxis() -> SetTitle("#scale[4.7]{cos #theta*}");
  //  comp2 -> GetXaxis()->SetTitleOffset(4.7);
  comp2 -> SetMarkerColor(kCyan+2);
  comp2 -> SetMarkerStyle(20);
  comp2 -> SetMarkerSize(1.2);
  comp2 -> SetLineColor(kCyan+2);
  comp2 -> SetLineWidth(2);
  comp2 -> GetXaxis() -> SetLabelFont(42);
  comp2 -> GetXaxis() -> SetLabelSize(0.12);
  comp2 -> GetYaxis() -> SetLabelSize(0.105);
  comp2 -> GetYaxis() -> SetLabelFont(42);
  //  comp2 -> GetXaxis() -> SetLabelOffset(0.03);
  comp2 -> GetYaxis() -> SetTitleFont(42);
  comp2 -> GetXaxis() -> SetTitleFont(42);
  comp2 -> GetYaxis() -> SetTitleSize(0.10);
  comp2 -> GetXaxis() -> SetTitleSize(0.15);
  comp2 -> GetYaxis() -> SetTitleOffset(0.6);
  comp2 -> GetXaxis() -> SetTitleOffset(1.3);
  //  comp2 -> GetYaxis() -> SetTitle("1#sigma-Down/Nom.");
  comp2 -> GetYaxis() -> SetTitle("Down/Nom.");
  comp2 -> GetXaxis() -> SetTitle("cos #theta*");
  comp2 -> GetYaxis() -> SetNdivisions(505);
  comp2 -> SetMaximum(max);
  comp2 -> SetMinimum(min);
  //------Take me out after check
  comp2 -> SetMaximum(1.03);
  comp2 -> SetMinimum(0.97);
  comp2 -> Draw("AP");
  norm1 -> Draw("Same");

  TGaxis *axis2 = new TGaxis(-2.0, min, 0.0, min, -1.0, 1.0, 505, "");
  axis2->SetTextAlign(9);
  axis2->SetLabelSize(0.12);
  axis2->SetLabelFont(42);
  axis2->SetLabelOffset(0.025);
  
  TGaxis *axis3 = new TGaxis(0.1, min, 2.0, min, -0.9, 1.0, 505, "");
  axis3->SetTextAlign(9);
  axis3->SetLabelSize(0.12);
  axis3->SetLabelFont(42);
  axis3->SetLabelOffset(0.025);
  
  if(NrSubPlots != 1){
    axis2->Draw();
    axis3->Draw();
  }

  TLine line2 = TLine(0.0, min, 0.0, max);
  line2.SetLineStyle(2);

  if(NTotalBins > 15)
    line2.Draw();

  c0 -> Draw();


  c0 -> Print((outputFile+".eps").c_str());
  c0 -> Print((outputFile+".png").c_str());

  //std::cout << outputFile.c_str() << std::endl;

  
  std::string tpng      = (outputFile+".png").c_str();
  size_t pos = tpng.find("/");
  tpng.erase(tpng.begin(), tpng.begin()+pos+1);

  std::string title = "Cos #theta*";

  std::ofstream* htmlfile = new std::ofstream();


  htmlfile->open(fHTMLLabel ,fstream::app);
  // if normal, print tab begin
  //if(strstr(title.c_str(), "norm") == NULL)
    *htmlfile << "<tr>"
              << "<td>"<< title <<"</td>";

  *htmlfile << "<td>";
  *htmlfile << "<p><a href=\"" << tpng.c_str() << "\">";
  *htmlfile << "<img src=\""   << tpng.c_str() << "\">";
  *htmlfile << "</td>" << std::endl;
  //if(strstr(title.c_str(), "norm") != NULL)
  //  *htmlfile << "</tr>"<< std::endl;
  htmlfile->close();


  delete fLegend;
  delete c0;

  fcloseall();

}

//Draw plot which contains all different templates (including variations) for one certain systematic effect: if only ONE systematic variation is considered
void PlotInterpolationCurves::DrawFinalPlotSingle(std::vector<TH1D > histo0, std::vector<std::string> Label, std::string outputFile, bool ResidualFlag, bool norm)
{
  SetAtlasStyle();

  int NTotalBins = histo0[0].GetNbinsX();

  std::cout << NTotalBins << std::endl;

  std::vector<TH1D*> histo;

  //Convert TH1D into TH1D*
  for (int i = 0; i != histo0.size(); i++) {
    
    histo0[i].Sumw2();

    if(norm)    
      histo0[i].Scale(1.0/histo0[i].Integral());
    
    TH1D * hist_dummy = &histo0[i];
    histo.push_back( hist_dummy);
    
  }

  TCanvas *c0;

  TPad *pad1;
  TPad *pad2;
  TPad *pad3;

  if(!ResidualFlag){
    c0   = new TCanvas("", "", 720, 720);
    pad1 = new TPad("pad1","pad1", 0.0, 0.0, 1.0, 1.0);
    pad2 = new TPad("pad2","pad2", 0.0, 0.0, 1.0, 0.0);
  }
  else{
    c0   = new TCanvas("", "", 720, 765);
    pad1 = new TPad("pad1","pad1", 0.0, 0.300, 1.0, 1.000);
    pad2 = new TPad("pad2","pad2", 0.0, 0.000, 1.0, 0.300);
    //    pad3 = new TPad("pad3","pad3", 0.0, 0.000, 1.0, 0.250);

    pad1->SetBottomMargin(0.010);
    pad1->SetBorderMode(0);
    pad2->SetBottomMargin(0.300);
    pad2->SetBorderMode(0);
    //    pad3->SetBottomMargin(0.300);
    //  pad3->SetBorderMode(0);

    pad1->SetTicks(1,1);
    pad2->SetTicks(1,1);
    // pad3->SetTicks(1,1);

    pad1->Draw();
    pad2->Draw();
    // pad3->Draw();

  }


  double nmax = histo[1] -> GetBinContent(histo[1] -> GetMaximumBin())*1.45;

  int nbins = histo[1] -> GetNbinsX();

  double lower_edge  = histo[0] -> GetBinLowEdge(1);
  double bin_width   = histo[0] -> GetBinWidth(1);
  double number_bins = histo[0] -> GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;


  //  fTotalSum1 = new TH1D("Sum1", "", nbins, lower_edge, upper_edge);

  /*  fTotalSum4 -> Add(fTTbarSum4);
  fTotalSum4 -> Add(fWjetsSum);
  fTotalSum4 -> Add(fZjetsSum);
  fTotalSum4 -> Add(fDibosonSum);
  fTotalSum4 -> Add(fSingleTopSum);
  fTotalSum4 -> Add(fQCDSum); */


  pad1 -> cd();

  for(int iHist = 0; iHist < histo.size(); ++iHist){

    //histo[iHist] = normalise(histo[iHist]);   //only for shape plots!!!
    histo[iHist] -> SetFillColor(kWhite);
    histo[iHist] -> SetLineWidth(2.25);
    //histo[iHist] -> SetLineStyle(iHist+1);
  
    if(iHist==2){
      histo[iHist] -> SetLineColor(kCyan+2);
      histo[iHist] -> SetLineStyle(2);
    }
    if(iHist==1){
      histo[iHist] -> SetLineColor(kOrange+8);
      histo[iHist] -> SetLineStyle(3);
    }
    //    if(iHist!=4) histo[iHist] -> SetLineColor(iHist+1);
    //    else histo[iHist] -> SetLineColor(kViolet+5);
    histo[iHist] -> SetMarkerSize(0);
    if(iHist == 0)
      histo[iHist] -> Draw();
    else
      histo[iHist] -> Draw("SAME");

  }

  TLine line1 = TLine(0.0, histo[0] -> GetMinimum(), 0.0, nmax/2.0);
  line1.SetLineStyle(2);
  
  if(NTotalBins > 15)
    line1.Draw();
    
  //nmax = histo[1] -> GetBinContent(histo[1] -> GetMaximumBin())*1.2; //only for shape plots!!!

  int counter_zero = 0;

  for(int ibin = 1; ibin <= nbins; ++ibin)
    if(histo[0]  -> GetBinContent(ibin) == 0)
      counter_zero++;

  int ngoodbins = nbins - counter_zero;

  TGraphErrors *comp1 = new TGraphErrors();
  TGraphErrors *comp2 = new TGraphErrors();

  comp1 = GetRatioDistribution(histo[0], histo[1], ngoodbins);
  if(histo.size() > 2)
    comp2 = GetRatioDistribution(histo[0], histo[2], ngoodbins);

  comp1->GetXaxis()->SetLimits(lower_edge, upper_edge); 
  comp2->GetXaxis()->SetLimits(lower_edge, upper_edge);

  double max = 1.15;
  double min = 0.85;

  if(norm){

    //    max = 1.1;
    //  min = 0.9;
    
  }

  std::stringstream width;
  if(bin_width > 9.5){
    width << std::setprecision(2) << bin_width;
  }
  else
    width << std::setprecision(1) << bin_width;

  //  histo[0] -> GetXaxis() -> SetTitle("#scale[3.5]{cos #theta*}");

  if(norm)
    histo[0] -> GetYaxis() -> SetTitle("#scale[1.45]{Normalised}");
  else
    histo[0] -> GetYaxis() -> SetTitle("#scale[1.45]{Events}");

  histo[0]->GetXaxis()->SetLabelSize(0.04);
  histo[0]->GetYaxis()->SetLabelSize(0.05);
  if(!norm){
    histo[0]->SetMinimum(0.001);
    histo[0] -> GetYaxis() -> SetTitleOffset(1.85);
  }
  else
    histo[0] -> GetYaxis() -> SetTitleOffset(1.85);

  //  histo[0] -> GetXaxis() -> SetTitleOffset(1.2);

  //define the legend...
  TLegend *fLegend = new TLegend(0.725, 0.65, 0.92, 0.92);
  for(unsigned int k = 0; k < histo.size(); ++k)
   fLegend -> AddEntry(histo[k],  Label[k].c_str(), "l");
  fLegend  -> SetFillColor(0);
  fLegend  -> SetLineColor(0);
  fLegend  -> SetBorderSize(0);
  fLegend  -> SetTextFont(72);
  fLegend  -> SetTextSize(0.0375);
  fLegend  -> Draw();

  histo[0]  -> SetMaximum(nmax);
  //histo[0]  -> SetMinimum(0.0);

  //set labels....
  TLatex l1;
  l1.SetTextAlign(9);
  l1.SetTextSize(0.04);
  l1.SetNDC();
  l1.DrawLatex(0.21, 0.750, fLumiLabel.c_str());
  TLatex l2;
  l2.SetTextAlign(9);
  l2.SetTextFont(72);
  l2.SetTextSize(0.04);
  l2.SetNDC();
  l2.DrawLatex(0.195, 0.865, "ATLAS");
  TLatex l3;
  l3.SetTextAlign(9);
  l3.SetTextSize(0.04);
  l3.SetNDC();
  l3.DrawLatex(0.300, 0.865, "Work in progress");
  //l3.DrawLatex(0.19, 0.81, "Templates used for interpolation");
  /*  TLatex l4;
  l4.SetTextAlign(9);
  l4.SetTextSize(0.042);
  l4.SetTextFont(72);
  l4.SetNDC();
  l4.DrawLatex(0.66, 0.85, (fBtagLabel+" "+fJetBinLabel).c_str()); */

  int NrSubPlots = fLeptonLabel.size();

  float SubPlotWidth     = 0.85/float(NrSubPlots);
  float SubPlotHalfWidth = SubPlotWidth/2.0;

  for(int iPlot = 0; iPlot < NrSubPlots; ++iPlot){

    //    std::cout << iPlot << "\t" << fLeptonLabel[iPlot].c_str() << std::endl;
    
    float SubPlotPosition = (1+2*iPlot)*SubPlotHalfWidth;

    if(iPlot == 0) SubPlotPosition += 0.10;
    if(iPlot == 1) SubPlotPosition += 0.05;

    if(NrSubPlots == 1) SubPlotPosition -= 0.02;

    std::cout << iPlot << "\t" << SubPlotWidth << "\t" << SubPlotHalfWidth << "\t" << SubPlotPosition << std::endl;
    
    TLatex l3;
    l3.SetTextAlign(9);
    l3.SetTextFont(42);
    l3.SetTextSize(0.037);
    l3.SetNDC();
    l3.DrawLatex(SubPlotPosition, 0.120, (fLeptonLabel[iPlot]).c_str());
    l3.DrawLatex(SubPlotPosition, 0.080, (fJetBinLabel[iPlot]).c_str());
    l3.DrawLatex(SubPlotPosition, 0.040, (fBTagLabel[iPlot]).c_str());

  }


  pad2->cd();

  TF1 *norm1 = new TF1("fa1","1", lower_edge, upper_edge);
  norm1 -> SetLineColor(kBlack);
  norm1 -> SetLineStyle(1);
  norm1 -> SetLineWidth(2);
  
  comp1 -> SetMarkerColor(kOrange+8);
  comp1 -> SetMarkerStyle(20);
  comp1 -> SetMarkerSize(1.2);
  comp1 -> SetLineColor(kOrange+8);
  comp1 -> SetLineWidth(2);
  comp1 -> GetXaxis() -> SetLabelFont(42);
  comp1 -> GetXaxis() -> SetLabelSize(0.12);
  comp1 -> GetYaxis() -> SetLabelSize(0.10);
  comp1 -> GetYaxis() -> SetLabelFont(42);
  comp1 -> GetXaxis() -> SetLabelOffset(0.03);
  comp1 -> GetYaxis() -> SetTitleFont(42);
  comp1 -> GetXaxis() -> SetTitleFont(42);
  comp1 -> GetYaxis() -> SetTitleSize(0.11); 
  comp1 -> GetXaxis() -> SetTitleSize(0.15);
  comp1 -> GetYaxis() -> SetTitleOffset(0.55);
  comp1 -> GetXaxis() -> SetTitleOffset(1.00);

  //  comp1 -> GetYaxis() -> SetTitle("1#sigma-Up/Nom.");
  comp1 -> GetYaxis() -> SetTitle("Var./Nom.");
  comp1 -> GetYaxis() -> CenterTitle();
  comp1 -> GetXaxis() -> SetNdivisions(505);
  comp1 -> GetYaxis() -> SetNdivisions(505);
  comp1 -> SetMaximum(max);
  comp1 -> SetMinimum(min);

  comp1 -> GetXaxis() -> SetTitle("cos #theta*");
  //comp1 -> SetMaximum(1.19);
  //comp1 -> SetMinimum(0.81);

  comp1 -> Draw("AP");
  norm1 -> Draw("Same");

  TLine line = TLine(0.0, 0.8, 0.0, 1.2);
  line.SetLineStyle(2);
   
  if(NTotalBins > 15)
    line.Draw();
 
  TGaxis *axis2 = new TGaxis(-2.0, min, 0.0, min, -1.0, 1.0, 505, "");
  axis2->SetTextAlign(9);
  axis2->SetLabelSize(0.12);
  axis2->SetLabelFont(42);
  axis2->SetLabelOffset(0.025);
  
  TGaxis *axis3 = new TGaxis(0.1, min, 2.0, min, -0.9, 1.0, 505, "");
  axis3->SetTextAlign(9);
  axis3->SetLabelSize(0.12);
  axis3->SetLabelFont(42);
  axis3->SetLabelOffset(0.025);

  if(NrSubPlots != 1){
    axis2->Draw();
    axis3->Draw();
  }
  
  TLine line2 = TLine(0.0, min, 0.0, max);
  line2.SetLineStyle(2);

  if(NTotalBins > 15)
    line2.Draw();

  c0 -> Draw();



  //  c0 -> Print((outputFile+".eps").c_str());
  //  c0 -> Print((outputFile+".gif").c_str());
  c0 -> Print((outputFile).c_str());

  //std::cout << outputFile.c_str() << std::endl;

  // geklaut von Anna
  //  std::string teps = (outputFile+".eps").c_str();
  std::string tpng      = (outputFile).c_str();
  //  std::string tgif = (outputFile+".gif").c_str();

  //  std::string teps = (outputFile+"_normalized.eps").c_str();
  //  std::string tpng_norm = (outputFile+"_normalized.png").c_str();
  //  std::string tgif = (outputFile+"_normalized.gif").c_str();

  //  r.replace( r.find( "Fuchs"), strlen("Fuchs"), "Hund");

  tpng.replace(tpng.find(fOutputFolder.c_str()), strlen(fOutputFolder.c_str()), "");
  // tpng.replace(tpng.find("/"),                   strlen("/"),                   "");

  std::string title = "Cos #theta*";

  std::ofstream* htmlfile = new std::ofstream();

  //  std::cout << "Open again: " << (_htmls[fHTMLLabel]).c_str() << std::endl;

  htmlfile->open((_htmls[fHTMLLabel]).c_str(),fstream::app);
  // if normal, print tab begin
  if(strstr(title.c_str(), "norm") == NULL)
    *htmlfile << "<tr>"
              << "<td>"<< title <<"</td>";

  *htmlfile << "<td>";
  *htmlfile << "<p><a href=\"" << tpng.c_str() << "\">";
  *htmlfile << "<img src=\""   << tpng.c_str() << "\">";
  *htmlfile << "</td>" << std::endl;
  if(strstr(title.c_str(), "norm") != NULL)
    *htmlfile << "</tr>"<< std::endl;
  htmlfile->close();


  delete fLegend;
  delete c0;

  fcloseall();

}


//Function usd to create a page including all important plots; called in MakePlots
void PlotInterpolationCurves::createHTML(std::string name){

  // create the page
  std::ofstream page;
  std::string pname = fOutputFolder+"/"+fHTMLLabel;

  struct stat buf;
  if(stat(pname.c_str(), &buf) == -1){

    page.open(pname.c_str());
    
    // main
    if(name == "index"){
      
      page << "<html><head><title> Comparison of different interpolations </title></head>" << std::endl;
      page << "<body>" << std::endl;
      page << "<h1> Comparison of different interpolations </h1>" << std::endl;
      page << "<a href=" << fHTMLLabel.c_str() << "> >=4 jets </a>"<< std::endl;
      
      page <<"<br> <br>" << std::endl;
    }
    page << "<table border = 1> <tr>"
	 << "<th> name </th>"
	 << "<th> variable </th>"
	 << "</tr>" << std::endl;
    _htmls[fHTMLLabel] = pname;
    page.close();

  }
  else{

    //    std::cout << pname.c_str() << "   already exists!" << std::endl;

  }

    return;

}

//Used in DrawFinalPlot to get histogram ratio
TGraphErrors *PlotInterpolationCurves::GetRatioDistribution(TH1D *DataHist, TH1D *SumHist, int ngoodbins)
{
  double x[ngoodbins],          xerr[ngoodbins];
  double y_new[ngoodbins], y_new_err[ngoodbins];
  double y_nom[ngoodbins], y_nom_err[ngoodbins];
  double quotient[ngoodbins], quotient_err[ngoodbins];

  int count_bins = 0;

  int nbins = SumHist -> GetNbinsX();

  for(int ibin = 1; ibin <= nbins; ++ibin){

    if(DataHist  -> GetBinContent(ibin) > 0){

      x[count_bins]             = SumHist  -> GetBinCenter(ibin);
      xerr[count_bins]          = 0.0;
      y_nom[count_bins]         = SumHist  -> GetBinContent(ibin);
      y_nom_err[count_bins]     = SumHist  -> GetBinError(ibin);
      y_new[count_bins]         = DataHist -> GetBinContent(ibin);
      y_new_err[count_bins]     = DataHist -> GetBinError(ibin);

      if(y_nom[count_bins] != 0 && y_new[count_bins] != 0){
	quotient[count_bins]      = y_nom[count_bins]/y_new[count_bins];
	quotient_err[count_bins]  = quotient[count_bins]*sqrt(pow(y_nom_err[count_bins]/y_nom[count_bins],2) + pow(y_new_err[count_bins]/y_new[count_bins],2));

	// if we assume 50% correlation...
	// quotient_err[count_bins]  = quotient[count_bins]*sqrt(pow(y_nom_err[count_bins]/y_nom[count_bins],2) + pow(y_new_err[count_bins]/y_new[count_bins],2) - 2.0*y_nom_err[count_bins]/y_nom[count_bins]*y_new_err[count_bins]/y_new[count_bins]*0.5);


      }

      count_bins++;

    }

  }

  TGraphErrors *gr_comp = new TGraphErrors(ngoodbins, x, quotient, xerr, quotient_err);
    // gr_comp->GetXaxis()->SetLimits(lower_edge, upper_edge);

  return gr_comp;
}


//Used in PrintGraphs to get chi2 values
double PlotInterpolationCurves::GetChi2Values(std::vector<double> k, std::vector<double> y, std::vector<double> yerr, TF1 fcn)
{
  
  TH1D *fTemplate;
  TH1D *fInterpol;

  if(k.size() == 7){
    fTemplate = new TH1D("", "", k.size(), -3.5, 3.5);
    fInterpol = new TH1D("", "", k.size(), -3.5, 3.5);
  }
  if(k.size() == 5){
    fTemplate = new TH1D("", "", k.size(), -2.5, 2.5);
    fInterpol = new TH1D("", "", k.size(), -2.5, 2.5);
  }
  if(k.size() == 3){
    fTemplate = new TH1D("", "", k.size(), -1.5, 1.5);
    fInterpol = new TH1D("", "", k.size(), -1.5, 1.5);
  }

  double chi2 = 0;

  for(int ibin = 1; ibin <= k.size(); ++ibin){

    fTemplate -> SetBinContent(ibin,  y[ibin-1]);
    // fTemplate -> SetBinError(ibin, yerr[ibin-1]);
    
    fInterpol -> SetBinContent(ibin, fcn.Eval(k[ibin-1]));

    chi2 += (y[ibin-1] - fcn.Eval(k[ibin-1]))*(y[ibin-1] - fcn.Eval(k[ibin-1]))/y[ibin-1];

  }
  //  double chi2 = fTemplate -> Chi2Test(fInterpol, "CHI2");

  return chi2/6;
}


//Used in PrintGraphs to get chi2 values
double PlotInterpolationCurves::GetChi2Values(std::vector<double> k, std::vector<double> y, std::vector<double> yerr, TF1 fcn_pos, TF1 fcn_neg)
{

  TH1D *fTemplate;
  TH1D *fInterpol;

  if(k.size() == 7){
    fTemplate = new TH1D("", "", k.size(), -3.5, 3.5);
    fInterpol = new TH1D("", "", k.size(), -3.5, 3.5);
  }
  if(k.size() == 5){
    fTemplate = new TH1D("", "", k.size(), -2.5, 2.5);
    fInterpol = new TH1D("", "", k.size(), -2.5, 2.5);
  }
  if(k.size() == 3){
    fTemplate = new TH1D("", "", k.size(), -1.5, 1.5);
    fInterpol = new TH1D("", "", k.size(), -1.5, 1.5);
  }

  double chi2 = 0;

  for(int ibin = 1; ibin <= k.size(); ++ibin){

    fTemplate -> SetBinContent(ibin,  y[ibin-1]);
    //    fTemplate -> SetBinError(ibin, yerr[ibin-1]);

    if(k[ibin-1] < 0){
      fInterpol -> SetBinContent(ibin, fcn_neg.Eval(k[ibin-1]));
      chi2 += (y[ibin-1] - fcn_neg.Eval(k[ibin-1]))*(y[ibin-1] - fcn_neg.Eval(k[ibin-1]))/y[ibin-1];

    }
    else{
      
      fInterpol -> SetBinContent(ibin, fcn_pos.Eval(k[ibin-1]));
      chi2 += (y[ibin-1] - fcn_pos.Eval(k[ibin-1]))*(y[ibin-1] - fcn_pos.Eval(k[ibin-1]))/y[ibin-1];


    }
  }

  //  double chi2 = fTemplate -> Chi2Test(fInterpol, "CHI2/NDF");

  return chi2/6;

}

//Used in PrintGraphs to get histogram ratio
TGraphErrors *PlotInterpolationCurves::GetRatioDistribution(std::vector<double> k, std::vector<double> y, std::vector<double> yerr, TF1 fcn)
{
  
  int ngoodbins = k.size();

  double x[ngoodbins],          xerr[ngoodbins];
  double y_new[ngoodbins], y_new_err[ngoodbins];
  double y_nom[ngoodbins], y_nom_err[ngoodbins];
  double quotient[ngoodbins], quotient_err[ngoodbins];

  int count_bins = 0;

  for(int i = 0; i < ngoodbins; ++i){

    x[count_bins]             = k[i];
    xerr[count_bins]          = 0.0;
    y_nom[count_bins]         = y[i];
    y_nom_err[count_bins]     = yerr[i];
    y_new[count_bins]         = fcn.Eval(k[i]);
    y_new_err[count_bins]     = 0.0;  // has to be calculated!!!
    
    
    if(y_nom[count_bins] != 0 && y_new[count_bins] != 0){
      
      quotient[count_bins]      = y_nom[count_bins]/y_new[count_bins]; 
      quotient_err[count_bins]  = quotient[count_bins]*sqrt(pow(y_nom_err[count_bins]/y_nom[count_bins],2) + pow(y_new_err[count_bins]/y_new[count_bins],2));
    
    }

    count_bins++;

  }

  TGraphErrors *gr_comp = new TGraphErrors(ngoodbins, x, quotient, xerr, quotient_err);

  return gr_comp;

}


//Used in PrintGraphs to get histogram ratio in case piecewise linear interpolation function considered
TGraphErrors *PlotInterpolationCurves::GetRatioDistribution(std::vector<double> k, std::vector<double> y, std::vector<double> yerr,  TF1 fcn_pos, TF1 fcn_neg)
{

  int ngoodbins = k.size();

  double x[ngoodbins],          xerr[ngoodbins];
  double y_new[ngoodbins], y_new_err[ngoodbins];
  double y_nom[ngoodbins], y_nom_err[ngoodbins];
  double quotient[ngoodbins], quotient_err[ngoodbins];

  int count_bins = 0;

  for(int i = 0; i < ngoodbins; ++i){

    if(k[i] < 0){

      x[count_bins]             = k[i];
      xerr[count_bins]          = 0.0;
      y_nom[count_bins]         = y[i];
      y_nom_err[count_bins]     = yerr[i];
      y_new[count_bins]         = fcn_neg.Eval(k[i]);
      y_new_err[count_bins]     = 0.0;  // has to be calculated!!!
      
      
      if(y_nom[count_bins] != 0 && y_new[count_bins] != 0){

	quotient[count_bins]      = y_nom[count_bins]/y_new[count_bins];
	quotient_err[count_bins]  = quotient[count_bins]*sqrt(pow(y_nom_err[count_bins]/y_nom[count_bins],2) + pow(y_new_err[count_bins]/y_new[count_bins],2));
	
      }

    }
    else{
      
      x[count_bins]             = k[i];
      xerr[count_bins]          = 0.0;
      y_nom[count_bins]         = y[i];
      y_nom_err[count_bins]     = yerr[i];
      y_new[count_bins]         = fcn_pos.Eval(k[i]);
      y_new_err[count_bins]     = 0.0;  // has to be calculated!!!

      
      if(y_nom[count_bins] != 0 && y_new[count_bins] != 0){

	quotient[count_bins]      = y_nom[count_bins]/y_new[count_bins];
	quotient_err[count_bins]  = quotient[count_bins]*sqrt(pow(y_nom_err[count_bins]/y_nom[count_bins],2) + pow(y_new_err[count_bins]/y_new[count_bins],2));

      }

    }

    count_bins++;

  }

  TGraphErrors *gr_comp = new TGraphErrors(ngoodbins, x, quotient, xerr, quotient_err);

  // gr_comp->GetXaxis()->SetLimits(lower_edge, upper_edge);

  return gr_comp;

}

TH1D *PlotInterpolationCurves::MakeHistogram(std::vector<double> k, TF1 fcn)
{
  //  TH1D fHelpHist = ;
}



TH1D *PlotInterpolationCurves::MakePointer(TH1D Histo)
{
  double lower_edge  = Histo.GetBinLowEdge(1);
  double bin_width   = Histo.GetBinWidth(1);
  double number_bins = Histo.GetNbinsX();
  double upper_edge  = lower_edge + number_bins*bin_width;

  TH1D *new_histo = new TH1D("", "", number_bins, lower_edge, upper_edge);

  for(int iBin = 1; iBin <= number_bins; ++iBin){

    double bin_content = Histo.GetBinContent(iBin);
    double bin_error   = Histo.GetBinError(iBin);
    
    new_histo -> SetBinContent(iBin, bin_content);
    new_histo -> SetBinError(iBin, bin_error);

  }

  return new_histo;

}
