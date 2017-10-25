//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 23 01:05:30 2012 by ROOT version 5.30/02
// from TTree EnsembleTree/EnsembleTree
// found on file: Output_Sample0_0.root
//////////////////////////////////////////////////////////

#ifndef OutputTreeReader_h
#define OutputTreeReader_h

#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class OutputTreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        F0;
   Double_t        FL;
   Double_t        FR;
   Double_t        N0;
   Double_t        NL;
   Double_t        NR;
   Double_t        Nges;
   Double_t        Wjets;
   Double_t        QCD;
   Double_t        RemBkg;
   Double_t        F0_nom;
   Double_t        FL_nom;
   Double_t        FR_nom;
   Double_t        N0_nom;
   Double_t        NL_nom;
   Double_t        NR_nom;
   Double_t        Nges_nom;
   Double_t        Wjets_nom;
   Double_t        QCD_nom;
   Double_t        RemBkg_nom;
   Double_t        F0_err;
   Double_t        FL_err;
   Double_t        FR_err;
   Double_t        N0_err;
   Double_t        NL_err;
   Double_t        NR_err;
   Double_t        Wjets_err;
   Double_t        QCD_err;
   Double_t        RemBkg_err;
   Double_t        N0_prof_err;
   Double_t        NL_prof_err;
   Double_t        NR_prof_err;
   Double_t        Wjets_prof_err;
   Double_t        QCD_prof_err;
   Double_t        RemBkg_prof_err;
   Double_t        nui_0;
   Double_t        nui_nom_0;
   Double_t        nui_err_0;
   Double_t        nui_prof_err_0;
   Double_t        nui_1;
   Double_t        nui_nom_1;
   Double_t        nui_err_1;
   Double_t        nui_prof_err_1;
   Double_t        nui_2;
   Double_t        nui_nom_2;
   Double_t        nui_err_2;
   Double_t        nui_prof_err_2;
   Double_t        nui_3;
   Double_t        nui_nom_3;
   Double_t        nui_err_3;
   Double_t        nui_prof_err_3;
   Double_t        nui_4;
   Double_t        nui_nom_4;
   Double_t        nui_err_4;
   Double_t        nui_prof_err_4;
   Double_t        nui_5;
   Double_t        nui_nom_5;
   Double_t        nui_err_5;
   Double_t        nui_prof_err_5;
   Double_t        nui_6;
   Double_t        nui_nom_6;
   Double_t        nui_err_6;
   Double_t        nui_prof_err_6;
   Double_t        nui_7;
   Double_t        nui_nom_7;
   Double_t        nui_err_7;
   Double_t        nui_prof_err_7;
   Double_t        nui_8;
   Double_t        nui_nom_8;
   Double_t        nui_err_8;
   Double_t        nui_prof_err_8;
   Double_t        nui_9;
   Double_t        nui_nom_9;
   Double_t        nui_err_9;
   Double_t        nui_prof_err_9;
   Double_t        nui_10;
   Double_t        nui_nom_10;
   Double_t        nui_err_10;
   Double_t        nui_prof_err_10;
   Double_t        nui_11;
   Double_t        nui_nom_11;
   Double_t        nui_err_11;
   Double_t        nui_prof_err_11;
   Double_t        nui_12;
   Double_t        nui_nom_12;
   Double_t        nui_err_12;
   Double_t        nui_prof_err_12;
   Double_t        nui_13;
   Double_t        nui_nom_13;
   Double_t        nui_err_13;
   Double_t        nui_prof_err_13;
   Double_t        nui_14;
   Double_t        nui_nom_14;
   Double_t        nui_err_14;
   Double_t        nui_prof_err_14;
   Double_t        MinuitStatus;
   //Double_t        NumNuisancePar;
   //std::vector<Double_t> nui_value;//(num_nuisance);
   //std::vector<Double_t> nui_nom;//(num_nuisance);
   //std::vector<Double_t> nui_err;//(num_nuisance);
   //std::vector<Double_t> nui_prof_err;//(num_nuisance);
   Double_t test;

   // List of branches
   TBranch        *b_F0;   //!
   TBranch        *b_FL;   //!
   TBranch        *b_FR;   //!
   TBranch        *b_N0;   //!
   TBranch        *b_NL;   //!
   TBranch        *b_NR;   //!
   TBranch        *b_Nges;   //!
   TBranch        *b_Wjets;   //!
   TBranch        *b_QCD;   //!
   TBranch        *b_RemBkg;   //!
   TBranch        *b_F0_nom;   //!
   TBranch        *b_FL_nom;   //!
   TBranch        *b_FR_nom;   //!
   TBranch        *b_N0_nom;   //!
   TBranch        *b_NL_nom;   //!
   TBranch        *b_NR_nom;   //!
   TBranch        *b_Wjets_nom;   //!
   TBranch        *b_QCD_nom;   //!
   TBranch        *b_RemBkg_nom;   //!
   TBranch        *b_F0_err;   //!
   TBranch        *b_FL_err;   //!
   TBranch        *b_FR_err;   //!
   TBranch        *b_N0_err;   //!
   TBranch        *b_NL_err;   //!
   TBranch        *b_NR_err;   //!
   TBranch        *b_Wjets_err;   //!
   TBranch        *b_QCD_err;   //!
   TBranch        *b_RemBkg_err;   //!
   TBranch        *b_N0_prof_err;   //!
   TBranch        *b_NL_prof_err;   //!
   TBranch        *b_NR_prof_err;   //!
   TBranch        *b_Wjets_prof_err;   //!
   TBranch        *b_QCD_prof_err;   //!
   TBranch        *b_RemBkg_prof_err;   //!
   TBranch        *b_nui_0;   //!
   TBranch        *b_nui_nom_0;   //!
   TBranch        *b_nui_err_0;   //!
   TBranch        *b_nui_prof_err_0;   //!
   TBranch        *b_nui_1;   //!
   TBranch        *b_nui_nom_1;   //!
   TBranch        *b_nui_err_1;   //!
   TBranch        *b_nui_prof_err_1;   //!
   TBranch        *b_nui_2;   //!
   TBranch        *b_nui_nom_2;   //!
   TBranch        *b_nui_err_2;   //!
   TBranch        *b_nui_prof_err_2;   //!
   TBranch        *b_nui_3;   //!
   TBranch        *b_nui_nom_3;   //!
   TBranch        *b_nui_err_3;   //!
   TBranch        *b_nui_prof_err_3;   //!
   TBranch        *b_nui_4;   //!
   TBranch        *b_nui_nom_4;   //!
   TBranch        *b_nui_err_4;   //!
   TBranch        *b_nui_prof_err_4;   //!
   TBranch        *b_nui_5;   //!
   TBranch        *b_nui_nom_5;   //!
   TBranch        *b_nui_err_5;   //!
   TBranch        *b_nui_prof_err_5;   //!
   TBranch        *b_nui_6;   //!
   TBranch        *b_nui_nom_6;   //!
   TBranch        *b_nui_err_6;   //!
   TBranch        *b_nui_prof_err_6;   //!
   TBranch        *b_nui_7;   //!
   TBranch        *b_nui_nom_7;   //!
   TBranch        *b_nui_err_7;   //!
   TBranch        *b_nui_prof_err_7;   //!
   TBranch        *b_nui_8;   //!
   TBranch        *b_nui_nom_8;   //!
   TBranch        *b_nui_err_8;   //!
   TBranch        *b_nui_prof_err_8;   //!
   TBranch        *b_nui_9;   //!
   TBranch        *b_nui_nom_9;   //!
   TBranch        *b_nui_err_9;   //!
   TBranch        *b_nui_prof_err_9;   //!
   TBranch        *b_nui_10;   //!
   TBranch        *b_nui_nom_10;   //!
   TBranch        *b_nui_err_10;   //!
   TBranch        *b_nui_prof_err_10;   //!
   TBranch        *b_nui_11;   //!
   TBranch        *b_nui_nom_11;   //!
   TBranch        *b_nui_err_11;   //!
   TBranch        *b_nui_prof_err_11;   //!
   TBranch        *b_nui_12;   //!
   TBranch        *b_nui_nom_12;   //!
   TBranch        *b_nui_err_12;   //!
   TBranch        *b_nui_prof_err_12;   //!
   TBranch        *b_nui_13;   //!
   TBranch        *b_nui_nom_13;   //!
   TBranch        *b_nui_err_13;   //!
   TBranch        *b_nui_prof_err_13;   //!
   TBranch        *b_nui_14;   //!
   TBranch        *b_nui_nom_14;   //!
   TBranch        *b_nui_err_14;   //!
   TBranch        *b_nui_prof_err_14;   //!
   TBranch        *b_MinuitStatus;   //!
   //TBranch        *b_NumNuisancePar;   //!
   //std::vector<TBranch*> b_nui_value;//(num_nuisance);
   //std::vector<TBranch*> b_nui_nom;//(num_nuisance);
   //std::vector<TBranch*> b_nui_err;//(num_nuisance);
   //std::vector<TBranch*> b_nui_prof_err;//(num_nuisance);

   OutputTreeReader(TTree *tree=0, int num_of_nui_parameters=0 );
   virtual ~OutputTreeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, int num_of_nui_parameters);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef OutputTreeReader_cxx
OutputTreeReader::OutputTreeReader(TTree *tree, int num_of_nui_parameters)
{
/*    std::vector<Double_t> nui_value(num_of_nui_parameters);
    std::vector<Double_t> nui_nom(num_of_nui_parameters);
    std::vector<Double_t> nui_err(num_of_nui_parameters);
    std::vector<Double_t> nui_prof_err(num_of_nui_parameters);
    std::vector<TBranch*> b_nui_value(num_of_nui_parameters);
    std::vector<TBranch*> b_nui_nom(num_of_nui_parameters);
    std::vector<TBranch*> b_nui_err(num_of_nui_parameters);
    std::vector<TBranch*> b_nui_prof_err(num_of_nui_parameters);
*/


// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Output_Sample0_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Output_Sample0_0.root");
      }
      f->GetObject("EnsembleTree",tree);

   }
   Init(tree, num_of_nui_parameters);
}

OutputTreeReader::~OutputTreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t OutputTreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t OutputTreeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void OutputTreeReader::Init(TTree *tree, int num_of_nui_parameters )
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("F0", &F0, &b_F0);
   fChain->SetBranchAddress("FL", &FL, &b_FL);
   fChain->SetBranchAddress("FR", &FR, &b_FR);
   fChain->SetBranchAddress("N0", &N0, &b_N0);
   fChain->SetBranchAddress("NL", &NL, &b_NL);
   fChain->SetBranchAddress("NR", &NR, &b_NR);
   fChain->SetBranchAddress("Nges", &Nges, &b_Nges);
   fChain->SetBranchAddress("Wjets", &Wjets, &b_Wjets);
   fChain->SetBranchAddress("QCD", &QCD, &b_QCD);
   fChain->SetBranchAddress("RemBkg", &RemBkg, &b_RemBkg);
   fChain->SetBranchAddress("F0_nom", &F0_nom, &b_F0_nom);
   fChain->SetBranchAddress("FL_nom", &FL_nom, &b_FL_nom);
   fChain->SetBranchAddress("FR_nom", &FR_nom, &b_FR_nom);
   fChain->SetBranchAddress("N0_nom", &N0_nom, &b_N0_nom);
   fChain->SetBranchAddress("NL_nom", &NL_nom, &b_NL_nom);
   fChain->SetBranchAddress("NR_nom", &NR_nom, &b_NR_nom);
   fChain->SetBranchAddress("Wjets_nom", &Wjets_nom, &b_Wjets_nom);
   fChain->SetBranchAddress("QCD_nom", &QCD_nom, &b_QCD_nom);
   fChain->SetBranchAddress("RemBkg_nom", &RemBkg_nom, &b_RemBkg_nom);
   fChain->SetBranchAddress("F0_err", &F0_err, &b_F0_err);
   fChain->SetBranchAddress("FL_err", &FL_err, &b_FL_err);
   fChain->SetBranchAddress("FR_err", &FR_err, &b_FR_err);
   fChain->SetBranchAddress("N0_err", &N0_err, &b_N0_err);
   fChain->SetBranchAddress("NL_err", &NL_err, &b_NL_err);
   fChain->SetBranchAddress("NR_err", &NR_err, &b_NR_err);
   fChain->SetBranchAddress("Wjets_err", &Wjets_err, &b_Wjets_err);
   fChain->SetBranchAddress("QCD_err", &QCD_err, &b_QCD_err);
   fChain->SetBranchAddress("RemBkg_err", &RemBkg_err, &b_RemBkg_err);
   fChain->SetBranchAddress("N0_prof_err", &N0_prof_err, &b_N0_prof_err);
   fChain->SetBranchAddress("NL_prof_err", &NL_prof_err, &b_NL_prof_err);
   fChain->SetBranchAddress("NR_prof_err", &NR_prof_err, &b_NR_prof_err);
   fChain->SetBranchAddress("Wjets_prof_err", &Wjets_prof_err, &b_Wjets_prof_err);
   fChain->SetBranchAddress("QCD_prof_err", &QCD_prof_err, &b_QCD_prof_err);
   fChain->SetBranchAddress("RemBkg_prof_err", &RemBkg_prof_err, &b_RemBkg_prof_err);
   /*fChain->SetBranchAddress("nui_0", &nui_0, &b_nui_0);
   fChain->SetBranchAddress("nui_nom_0", &nui_nom_0, &b_nui_nom_0);
   fChain->SetBranchAddress("nui_err_0", &nui_err_0, &b_nui_err_0);
   fChain->SetBranchAddress("nui_prof_err_0", &nui_prof_err_0, &b_nui_prof_err_0);
   fChain->SetBranchAddress("nui_1", &nui_1, &b_nui_1);
   fChain->SetBranchAddress("nui_nom_1", &nui_nom_1, &b_nui_nom_1);
   fChain->SetBranchAddress("nui_err_1", &nui_err_1, &b_nui_err_1);
   fChain->SetBranchAddress("nui_prof_err_1", &nui_prof_err_1, &b_nui_prof_err_1);*/
   fChain->SetBranchAddress("MinuitStatus", &MinuitStatus, &b_MinuitStatus);
   //fChain->SetBranchAddress("NumNuisancePar", &NumNuisancePar, &b_NumNuisancePar);

   //Branch address for nuisance parameters:
  for (int i = 0; i != num_of_nui_parameters; ++i) 
  {
	//std::cout << "fChain for loop" << std::endl;

	if (i == 0) {
	fChain->SetBranchAddress("nui_0", &nui_0, &b_nui_0);
	fChain->SetBranchAddress("nui_nom_0", &nui_nom_0, &b_nui_nom_0);
   	fChain->SetBranchAddress("nui_err_0", &nui_err_0, &b_nui_err_0);
   	fChain->SetBranchAddress("nui_prof_err_0", &nui_prof_err_0, &b_nui_prof_err_0); }
	if (i == 1) {
   	fChain->SetBranchAddress("nui_1", &nui_1, &b_nui_1);
   	fChain->SetBranchAddress("nui_nom_1", &nui_nom_1, &b_nui_nom_1);
   	fChain->SetBranchAddress("nui_err_1", &nui_err_1, &b_nui_err_1);
   	fChain->SetBranchAddress("nui_prof_err_1", &nui_prof_err_1, &b_nui_prof_err_1); }
	if (i == 2) {
   	fChain->SetBranchAddress("nui_2", &nui_2, &b_nui_2);
   	fChain->SetBranchAddress("nui_nom_2", &nui_nom_2, &b_nui_nom_2);
   	fChain->SetBranchAddress("nui_err_2", &nui_err_2, &b_nui_err_2);
   	fChain->SetBranchAddress("nui_prof_err_2", &nui_prof_err_2, &b_nui_prof_err_2); }
	if (i == 3) {
   	fChain->SetBranchAddress("nui_3", &nui_3, &b_nui_3);
   	fChain->SetBranchAddress("nui_nom_3", &nui_nom_3, &b_nui_nom_3);
   	fChain->SetBranchAddress("nui_err_3", &nui_err_3, &b_nui_err_3);
   	fChain->SetBranchAddress("nui_prof_err_3", &nui_prof_err_3, &b_nui_prof_err_3); }
	if (i == 4) {
   	fChain->SetBranchAddress("nui_4", &nui_4, &b_nui_4);
   	fChain->SetBranchAddress("nui_nom_4", &nui_nom_4, &b_nui_nom_4);
   	fChain->SetBranchAddress("nui_err_4", &nui_err_4, &b_nui_err_4);
   	fChain->SetBranchAddress("nui_prof_err_4", &nui_prof_err_4, &b_nui_prof_err_4); }
	if (i == 5) {
   	fChain->SetBranchAddress("nui_5", &nui_5, &b_nui_5);
   	fChain->SetBranchAddress("nui_nom_5", &nui_nom_5, &b_nui_nom_5);
   	fChain->SetBranchAddress("nui_err_5", &nui_err_5, &b_nui_err_5);
   	fChain->SetBranchAddress("nui_prof_err_5", &nui_prof_err_5, &b_nui_prof_err_5); }
	if (i == 6) {
   	fChain->SetBranchAddress("nui_6", &nui_6, &b_nui_6);
   	fChain->SetBranchAddress("nui_nom_6", &nui_nom_6, &b_nui_nom_6);
   	fChain->SetBranchAddress("nui_err_6", &nui_err_6, &b_nui_err_6);
   	fChain->SetBranchAddress("nui_prof_err_6", &nui_prof_err_6, &b_nui_prof_err_6); }
	if (i == 7) {
   	fChain->SetBranchAddress("nui_7", &nui_7, &b_nui_7);
   	fChain->SetBranchAddress("nui_nom_7", &nui_nom_7, &b_nui_nom_7);
   	fChain->SetBranchAddress("nui_err_7", &nui_err_7, &b_nui_err_7);
   	fChain->SetBranchAddress("nui_prof_err_7", &nui_prof_err_7, &b_nui_prof_err_7); }
	if (i == 8) {
   	fChain->SetBranchAddress("nui_8", &nui_8, &b_nui_8);
   	fChain->SetBranchAddress("nui_nom_8", &nui_nom_8, &b_nui_nom_8);
   	fChain->SetBranchAddress("nui_err_8", &nui_err_8, &b_nui_err_8);
   	fChain->SetBranchAddress("nui_prof_err_8", &nui_prof_err_8, &b_nui_prof_err_8); }
	if (i == 9) {
   	fChain->SetBranchAddress("nui_9", &nui_9, &b_nui_9);
   	fChain->SetBranchAddress("nui_nom_9", &nui_nom_9, &b_nui_nom_9);
   	fChain->SetBranchAddress("nui_err_9", &nui_err_9, &b_nui_err_9);
   	fChain->SetBranchAddress("nui_prof_err_9", &nui_prof_err_9, &b_nui_prof_err_9); }
	if (i == 10) {
   	fChain->SetBranchAddress("nui_10", &nui_10, &b_nui_10);
   	fChain->SetBranchAddress("nui_nom_10", &nui_nom_10, &b_nui_nom_10);
   	fChain->SetBranchAddress("nui_err_10", &nui_err_10, &b_nui_err_10);
   	fChain->SetBranchAddress("nui_prof_err_10", &nui_prof_err_10, &b_nui_prof_err_10); }
	if (i == 11) {
   	fChain->SetBranchAddress("nui_11", &nui_11, &b_nui_11);
   	fChain->SetBranchAddress("nui_nom_11", &nui_nom_11, &b_nui_nom_11);
   	fChain->SetBranchAddress("nui_err_11", &nui_err_11, &b_nui_err_11);
   	fChain->SetBranchAddress("nui_prof_err_11", &nui_prof_err_11, &b_nui_prof_err_11); }
	if (i == 12) {
   	fChain->SetBranchAddress("nui_12", &nui_12, &b_nui_12);
   	fChain->SetBranchAddress("nui_nom_12", &nui_nom_12, &b_nui_nom_12);
   	fChain->SetBranchAddress("nui_err_12", &nui_err_12, &b_nui_err_12);
   	fChain->SetBranchAddress("nui_prof_err_12", &nui_prof_err_12, &b_nui_prof_err_12); }
	if (i == 13) {
   	fChain->SetBranchAddress("nui_13", &nui_13, &b_nui_13);
   	fChain->SetBranchAddress("nui_nom_13", &nui_nom_13, &b_nui_nom_13);
   	fChain->SetBranchAddress("nui_err_13", &nui_err_13, &b_nui_err_13);
   	fChain->SetBranchAddress("nui_prof_err_13", &nui_prof_err_13, &b_nui_prof_err_13); }
	if (i == 14) {
   	fChain->SetBranchAddress("nui_14", &nui_14, &b_nui_14);
   	fChain->SetBranchAddress("nui_nom_14", &nui_nom_14, &b_nui_nom_14);
   	fChain->SetBranchAddress("nui_err_14", &nui_err_14, &b_nui_err_14);
   	fChain->SetBranchAddress("nui_prof_err_14", &nui_prof_err_14, &b_nui_prof_err_14); }



	//fChain->SetBranchAddress(Form("nui_%d", i), &(nui_value[i]), &(b_nui_value[i]) );
	//fChain->SetBranchAddress(Form("nui_nom_%d", i), &(nui_nom[i]), &(b_nui_nom[i]) );
	//fChain->SetBranchAddress(Form("nui_err_%d", i), &(nui_err[i]), &(b_nui_err[i]) );
	//fChain->SetBranchAddress(Form("nui_prof_err_%d", i), &(nui_prof_err[i]), &(b_nui_prof_err[i] )  );
  }


   Notify();
}

Bool_t OutputTreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void OutputTreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t OutputTreeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef OutputTreeReader_cxx
