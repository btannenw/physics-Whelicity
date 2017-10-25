#include <TFile.h>
#include <TH1.h>
#include <TObjArray.h>
#include <math.h>

#include "AtlasStyle.C"

void readHists(TFile* file, TObjArray* array, bool isUpVariation);
void set1DAddress(TObjArray* array, TFile* file, std::string variation, bool isUpVariation);
void set1DAddress_oneSided(TObjArray* array, TFile* file, std::string variation);
void set1DAddress_compareTwo(TObjArray* array, TFile* file, std::string variation1, std::string variation2);
//void set1DAddress_twoSample(TObjArray* array, TFile* file, std::string variation, bool isUpVariation);
void makeSum(TObjArray* array, TH1D*& hist);

void makeSystErrorBand()
{
  SetAtlasStyle();

  // read file
  TFile* fTemplates = new TFile("Templates_110404_syst_2incl_KLF5jOPT_el_mu.root","READ");

  // declare canvasses
  TCanvas* c1 = new TCanvas("c1","test", 800, 800);

  // declare arrays and sum histograms
  TObjArray* upVariations = new TObjArray();
  TObjArray* downVariations = new TObjArray();
  TH1D* h_sumUp(0);
  TH1D* h_sumDown(0);
  
  // fill arrays
  readHists(fTemplates, upVariations, true);
  readHists(fTemplates, downVariations, false);

  // create sum histograms
  makeSum(upVariations, h_sumUp);
  makeSum(downVariations, h_sumDown);
  
  // make negative histogram... negative, and divide by nominal
  h_sumDown->Scale(-1);
  TH1D* h_nom = (TH1D*)fTemplates->Get("PseudoData3W");
  h_sumDown->Divide(h_nom);
  h_sumUp->Divide(h_nom);
  
  // draw for human eye comparison
  c1->cd();
  
  h_sumUp->SetLineColor(kBlue);
  h_sumDown->SetLineColor(kRed);
  h_sumUp->SetMinimum( h_sumDown->GetMinimum()*1.1 );
  
  h_sumUp->Draw();
  h_sumDown->Draw("same");

}

void readHists(TFile* file, TObjArray* array, bool isUpVariation)
{
  // missing: Parton Shower, ME Gen
  //set1DAddress(array, file, "", isUpVariation);

  // ME generator variation
  set1DAddress_compareTwo(array, file, "PseudoData3W_105860", "PseudoData3W_105200");
  
  // Parton shower variation
  set1DAddress_compareTwo(array, file, "PseudoData3W_117050", "PseudoData3W_105860");

  // Radiation variation
  set1DAddress_oneSided(array, file, "PseudoData3W_110408"); // mu = 0.5
  /*
  // JEff variation
  set1DAddress_oneSided(array, file, "PseudoData3W_jeff");
  
  // one-sided JER variation
  set1DAddress_oneSided(array, file, "PseudoData3W_jer_DataMC_Difference");
  set1DAddress_oneSided(array, file, "PseudoData3W_jer_Noise_ForwardRegion");

  // electron variations
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_ELE_RECO_UP", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_ELE_ID_UP", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_ELE_TRIGGER_UP", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_ELE_RECO_DOWN", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_ELE_ID_DOWN", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_ELE_TRIGGER_DOWN", isUpVariation);
  
  // muon variations
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_MUON_RECO_UP", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_MUON_ID_UP", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_MUON_TRIGGER_UP", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_MUON_RECO_DOWN", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_MUON_ID_DOWN", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_MUON_TRIGGER_DOWN", isUpVariation);
  
  // JER variations
  set1DAddress(array, file, "PseudoData3W_jer_down_NP0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP4", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP5", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP6", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP7", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_down_NP8", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP4", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP5", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP6", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP7", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jer_up_NP8", isUpVariation);
  
  // JES variations
  set1DAddress(array, file, "PseudoData3W_jes_up_EffectiveNP_Modelling1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_up_EffectiveNP_Statistical1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_up_EtaIntercalibration_TotalStat", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_up_FlavourComp", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_up_FlavourResponse", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_up_RhoTopology", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_down_EffectiveNP_Modelling1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_down_EffectiveNP_Statistical1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_down_EtaIntercalibration_TotalStat", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_down_FlavourComp", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_down_FlavourResponse", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jes_down_RhoTopology", isUpVariation);
  
  // JVF variation
  set1DAddress(array, file, "PseudoData3W_jvf_up", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_jvf_down", isUpVariation);
  
  // b-tagging up variations
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarUp_0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarUp_1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarUp_2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarUp_3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarUp_4", isUpVariation);*/
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarUp_5", isUpVariation);
  /*set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarUp_0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarUp_1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarUp_2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarUp_3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_4", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_5", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_6", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_7", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_8", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_9", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_10", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarUp_11", isUpVariation);
  
  // b-tagging down variations
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarDown_0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarDown_1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarDown_2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarDown_3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarDown_4", isUpVariation);*/
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_bTagVarDown_5", isUpVariation);
  /*set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarDown_0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarDown_1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarDown_2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_cTagVarDown_3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_0", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_1", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_2", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_3", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_4", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_5", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_6", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_7", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_8", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_9", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_10", isUpVariation);
  set1DAddress(array, file, "PseudoData3W_nominal_scaleFactor_BTAG_misTagVarDown_11", isUpVariation);
  */
}

void set1DAddress(TObjArray* array, TFile* file, std::string variation, bool isUpVariation)
{
  TH1D* nom = (TH1D*)file->Get("PseudoData3W");
  TH1D* diff = (TH1D*)nom->Clone();
  TH1D* histo(0);
  //cout<<"----------------------------------"<<endl;
  if(file->GetListOfKeys()->Contains(variation.c_str()) ) //proceed if histogram in file
    {
      histo=(TH1D*)file->Get(variation.c_str()); //get histogram from input file
      diff->Add(histo, -1);
      diff->SetName(variation.c_str());
      
      if (isUpVariation && (variation.find("UP")!=string::npos || variation.find("Up")!=string::npos || variation.find("up")!=string::npos) ){ 
	diff->Draw();
	array->AddLast(diff); // add histo to TObjArray
	//cout<<"added "<<variation<<endl;
      }
      else if (!isUpVariation && (variation.find("DOWN")!=string::npos || variation.find("Down")!=string::npos || variation.find("down")!=string::npos) ){ 
	diff->Draw();
	array->AddLast(diff); // add histo to TObjArray
	//cout<<"added "<<variation<<endl;
      }
      else{
	if( !( (isUpVariation && (variation.find("DOWN")!=string::npos || variation.find("Down")!=string::npos || variation.find("down")!=string::npos)) ||
	       (!isUpVariation && (variation.find("UP")!=string::npos || variation.find("Up")!=string::npos || variation.find("up")!=string::npos)) ) )
	  cout<<"Histogram "<<variation<<" FOUND but not formatted in Up/Down correctly. isUpVariation: "<<isUpVariation<<endl;
      }
      
    }
  else{
    cout<<"Histogram "<<variation<<" NOT FOUND."<<endl;
  }

  diff->Print("all");
}


void set1DAddress_oneSided(TObjArray* array, TFile* file, std::string variation)
{
  TH1D* nom = (TH1D*)file->Get("PseudoData3W");
  if(variation.find("110408")!=string::npos) // if radiation variation, need to use afii nominal
    nom = (TH1D*)file->Get("PseudoData3W_110404_AFII");
  
  TH1D* diff = (TH1D*)nom->Clone();
  TH1D* histo(0);
  
  if(file->GetListOfKeys()->Contains(variation.c_str()) ) //proceed if histogram in file
    {
      histo=(TH1D*)file->Get(variation.c_str()); //get histogram from input file
      diff->Add(histo, -1);
      diff->SetName(variation.c_str());
      
      diff->Draw();
      array->AddLast(diff); // add histo to TObjArray
    }
  else{
    cout<<"One-sided histogram "<<variation<<" NOT FOUND."<<endl;
  }

  diff->Print("all");
  
}

void set1DAddress_compareTwo(TObjArray* array, TFile* file, std::string variation1, std::string variation2)
{
  TH1D* hist1(0);
  TH1D* hist2(0);
  
  if(file->GetListOfKeys()->Contains(variation1.c_str()) ) //proceed if histogram in file
    {
      hist1=(TH1D*)file->Get(variation1.c_str()); //get histogram from input file
      TH1D* diff = (TH1D*)hist1->Clone();
      if(file->GetListOfKeys()->Contains(variation2.c_str()) ) //proceed if histogram in file
	{
	  hist2=(TH1D*)file->Get(variation2.c_str()); //get histogram from input file
	  diff->Add(hist2, -1);
	  diff->SetName(variation1.c_str());
	  
	  diff->Draw();
	  array->AddLast(diff); // add histo to TObjArray
	  diff->Print("all");
	}
      else{
	cout<<"compareTwo histogram "<<variation1<<" NOT FOUND."<<endl;
      }
    }
  else{
    cout<<"compareTwo histogram "<<variation2<<" NOT FOUND."<<endl;
  }

}

void makeSum(TObjArray* array, TH1D*& hist)
{
  TObjArrayIter next(array);
  TObject* object;
  bool first = true;
  double binContent = 0;


  while( (object=next()) ){ //iterate over objects in new array, add to a_data
    TString type=object->ClassName();
    TString name=object->GetName();

    if(type=="TH1D"){
      TH1D* h1=(TH1D*)object;
      if (first){
	hist = h1;
	//for(int i = 0; i < hist->GetNbinsX() + 1; i++)
	//  hist->SetBinContent(i, h1->GetBinContent(i)*h1->GetBinContent(i));
	hist->Multiply(h1);
	first = false;
      }
      else{
	for(int i = 0; i < hist->GetNbinsX() + 1; i++){
	  //binContent = hist->GetBinContent(i)*hist->GetBinContent(i) + h1->GetBinContent(i)*h1->GetBinContent(i);
	  //hist -> SetBinContent(i, sqrt(binContent));
	  binContent = hist->GetBinContent(i) + h1->GetBinContent(i)*h1->GetBinContent(i);
       	  hist -> SetBinContent(i, binContent);
	 //hist -> SetBinContent(i, sqrt(binContent));
	}
      }
      
    }// end type == TH1D safety condition
  }// end while loop

  for(int i = 0; i < hist->GetNbinsX() + 1; i++)
    hist->SetBinContent(i, sqrt(hist->GetBinContent(i)));

  hist->Print("all");
}
