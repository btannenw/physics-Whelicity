
// compile with:  g++ -o PDFPlotter PDFPlotter.c `root-config --cflags --glibs`

//REF: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TopPdfUncertainty

#include <TFile.h>
#include <THStack.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TGraph.h>
#include <sys/stat.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
//#include <readparameters.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TVectorD.h>
#include <TGraph.h>
#include <TLine.h>
#include <TBox.h>

#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <cmath>


#include <math.h>
#include <vector>
#include <utility>
#include <string>
#include <AtlasStyle.h>

using namespace std;

std::string channel;
std::string frac;
template<typename T> T* vector2Array(std::vector<T> vec);
std::vector<double> ReadValues(const char * filename);
pair <double, pair <double,double > > hessian_uncertainty (vector <double> uncerts);
pair <double, pair <double, double > > hessian_asymmetric_uncertainty (vector <double> uncerts);
pair < double , pair <double, double > > RMS_NNPDF_uncertainty(vector<double> uncerts);
pair <double, double > envelope (double central, vector < pair < double , pair <double, double> > > pair_uncerts);
pair <double, double> calculate_envelope(double central, vector <double> uncert_hess, vector < double > uncert_asymm, vector < double > uncert_NNPDF );
bool doestringcontain(std::string str1, std::string str2); 

//void PDFPlotter(int argc, char **argv)
void PDFPlotter(const char *ch_CT10, const char *ch_MSTW, const char *ch_NNPDF)
{
  //SetAtlasStyle();
  // if (!(argc >= 4)) { // i.e. 2 parameters ...
  //   //Par1: data histo name and name for subdirectory
  //   std::cout << "at least 3 parameter needed: ./PDFPlotter [CT10 file] [MSTW file] [NNPDF file] [opt: HERAPDF file]" << std::endl;
  //   return 1;
  // } 

  std::vector <double> CT10 = ReadValues(ch_CT10);
  std::vector <double> MSTW = ReadValues(ch_MSTW);
  std::vector <double> NNPDF = ReadValues(ch_NNPDF);
  //std::vector <double> HERAPDF;
  //if (argc == 5) HERAPDF= ReadValues(argv[4]);

  std::string output_basebath(ch_CT10);
  output_basebath.replace(output_basebath.find("CT10"),4,"");
  output_basebath.replace(output_basebath.find(".txt"),4,"");
  output_basebath="./pdf_results_Whel_oct31/"+output_basebath;
 
  if (doestringcontain(output_basebath, "lephad")) channel="lephad";
  else if (doestringcontain(output_basebath, "had")) channel="had";
	else channel="lep";

 
  if (doestringcontain(output_basebath, "F0")) {
	output_basebath="./pdf_results_Whel_oct31/F0_"+channel;
	frac="F0";
  }
  else if (doestringcontain(output_basebath, "FL")) {
	output_basebath="./pdf_results_Whel_oct31/FL_"+channel;
	frac="FL";
  }
  else if (doestringcontain(output_basebath, "FR")) {
	output_basebath="./pdf_results_Whel_oct31/FR_"+channel;
	frac="FR";
  }
 
  int bins = CT10.size();
  if (bins < MSTW.size()) bins = MSTW.size();
  if (bins < NNPDF.size()) bins = NNPDF.size();

  pair < double , pair <double, double > > CT10_uncertainty 	= hessian_uncertainty (CT10);
  pair <double, pair <double, double > > MSTW_uncertainty 	= hessian_asymmetric_uncertainty (MSTW);
  pair < double , pair <double, double > > NNPDF_uncertainty	= RMS_NNPDF_uncertainty (NNPDF);

//pair < double , pair <double, double > > HERAPDF_uncertainty;
//if (argc == 5)	HERAPDF_uncertainty = RMS_NNPDF_uncertainty (HERAPDF);
 
std::cout<<"CT10: "<<CT10_uncertainty.first<<" + "<<CT10_uncertainty.second.first<<" - "<<CT10_uncertainty.second.second<<std::endl;
std::cout<<"MSTW: "<<MSTW_uncertainty.first<<" + "<<MSTW_uncertainty.second.first<<" - "<<MSTW_uncertainty.second.second<<std::endl;
std::cout<<"NNPDF: "<<NNPDF_uncertainty.first<<" + "<<NNPDF_uncertainty.second.first<<" - "<<NNPDF_uncertainty.second.second<<std::endl;
//if (argc == 5) std::cout<<"HERAPDF: "<<HERAPDF_uncertainty.first<<" + "<<HERAPDF_uncertainty.first<<" - "<<HERAPDF_uncertainty.second.second<<std::endl;

  TGraph * g_CT10 = new TGraph(bins);
  TGraph * g_MSTW = new TGraph(bins);
  TGraph * g_NNPDF = new TGraph(bins);
  //TGraph * g_HERAPDF = new TGraph(bins);
  g_CT10->SetTitle("");
  g_MSTW->SetTitle("");
  g_NNPDF->SetTitle("");


  for (int i = 0; i< CT10.size(); i++){
    if (i==0) 	g_CT10->SetPoint(i, i, CT10.at(i));
    else 	g_CT10->SetPoint(i, (i+1)/2, CT10.at(i));
  }

  for (int i = 0; i< MSTW.size(); i++){
    if (i==0) 	g_MSTW->SetPoint(i, i, MSTW.at(i));
    else 	g_MSTW->SetPoint(i, (i+1)/2, MSTW.at(i));
  }

  for (int i = 0; i< NNPDF.size(); i++){
    if (i==0) 	g_NNPDF->SetPoint(i, i, NNPDF.at(i));
    else 	g_NNPDF->SetPoint(i, (i+1)/2, NNPDF.at(i));
  }

  // for (int i = 0; i< HERAPDF.size(); i++){
  //   if (i==0) 	g_HERAPDF->SetPoint(i, i, HERAPDF.at(i));
  //   else 	g_HERAPDF->SetPoint(i, (i+1)/2, HERAPDF.at(i));
  // }

  TCanvas canni;
  g_CT10->SetMarkerStyle(20);
  g_MSTW->SetMarkerStyle(21);
  g_NNPDF->SetMarkerStyle(22);
  //g_HERAPDF->SetMarkerStyle(23);

  g_CT10->SetMarkerColor(kBlue);
  g_MSTW->SetMarkerColor(kRed);
  g_NNPDF->SetMarkerColor(kGreen);
  //g_HERAPDF->SetMarkerColor(kCyan);

  g_NNPDF->Draw("AP");
	TLegend * leg1 = new TLegend(0.65, 0.75, 0.85, 0.89);
  leg1->SetBorderSize(0);
	leg1->AddEntry(g_CT10, "CT10", "p");
	leg1->AddEntry(g_MSTW, "MSTW", "p");
	leg1->AddEntry(g_NNPDF, "NNPDF", "p");
	//if (argc == 5) leg1->AddEntry(g_HERAPDF, "HERAPDF", "p");
	leg1->Draw();
  
 if(frac=="F0"){
  	g_NNPDF->GetYaxis()->SetTitle("F_{0}");
	if(channel=="lephad")
		g_NNPDF->GetYaxis()->SetRangeUser(0.70, 0.72);
	else if(channel=="lep")
		g_NNPDF->GetYaxis()->SetRangeUser(0.69, 0.715);
	else if(channel=="had")
                g_NNPDF->GetYaxis()->SetRangeUser(0.725, 0.735);
 }
else if(frac=="FL"){
  	g_NNPDF->GetYaxis()->SetTitle("F_{L}");
	if(channel=="lephad")
                g_NNPDF->GetYaxis()->SetRangeUser(0.28, 0.31);
        else if(channel=="lep")
                g_NNPDF->GetYaxis()->SetRangeUser(0.285, 0.315);
        else if(channel=="had")
                g_NNPDF->GetYaxis()->SetRangeUser(0.32, 0.34);
 }
else if(frac=="FR"){
  	g_NNPDF->GetYaxis()->SetTitle("F_{R}");
	if(channel=="lephad")
                g_NNPDF->GetYaxis()->SetRangeUser(-0.0035, 0.004);
        else if(channel=="lep")
                g_NNPDF->GetYaxis()->SetRangeUser(0.001, 0.008);
        else if(channel=="had")
                 g_NNPDF->GetYaxis()->SetRangeUser(-0.085, -0.025);
 }
  
  g_NNPDF->GetYaxis()->SetTitleOffset(1.3);
  g_NNPDF->GetXaxis()->SetTitle("PDF set");
  
  TLine line_CT10(0, CT10_uncertainty.first, CT10.size()/2, CT10_uncertainty.first); 
  line_CT10.SetLineColor(kBlue);
  line_CT10.Draw();
  TBox box_CT10(0, CT10_uncertainty.first+CT10_uncertainty.second.first, CT10.size()/2, CT10_uncertainty.first-CT10_uncertainty.second.second);
  box_CT10.SetLineColor(kBlue);
  box_CT10.SetFillColor(kBlue);
  box_CT10.SetFillStyle(3004);
  box_CT10.Draw();

std::cout<<"DEBUG ! MSTW cen sec1 sec: "<<MSTW_uncertainty.first<<" / "<<MSTW_uncertainty.second.first<<" / "<<MSTW_uncertainty.second.second<<std::endl;
  TLine line_MSTW(0, MSTW_uncertainty.first, MSTW.size()/2, MSTW_uncertainty.first); 
  line_MSTW.SetLineColor(kRed);
  line_MSTW.Draw();
  TBox box_MSTW(0, MSTW_uncertainty.first+MSTW_uncertainty.second.first, MSTW.size()/2, MSTW_uncertainty.first-MSTW_uncertainty.second.second);
  box_MSTW.SetLineColor(kRed);
  box_MSTW.SetFillColor(kRed);
  box_MSTW.SetFillStyle(3004);
  box_MSTW.Draw();


  TLine line_NNPDF(0, NNPDF_uncertainty.first, NNPDF.size()/2, NNPDF_uncertainty.first); 
  line_NNPDF.SetLineColor(kGreen);
  line_NNPDF.Draw();
  TBox box_NNPDF(0, NNPDF_uncertainty.first+NNPDF_uncertainty.second.first, NNPDF.size()/2, NNPDF_uncertainty.first-NNPDF_uncertainty.second.second);
  box_NNPDF.SetLineColor(kGreen);
  box_NNPDF.SetFillColor(kGreen);
  box_NNPDF.SetFillStyle(3004);
  box_NNPDF.Draw();

  pair <double, double> overall_unc = calculate_envelope(0, CT10, MSTW, NNPDF);

  TLine line_total(0, overall_unc.first, bins/2., overall_unc.first); 
  line_total.SetLineColor(kMagenta);
  line_total.Draw();

  TLine line_total2(bins/2., overall_unc.first+overall_unc.second, bins/2., overall_unc.first-overall_unc.second); 
  line_total2.SetLineColor(kMagenta);
  line_total2.SetLineStyle(1);
  line_total2.Draw();
/*
  TBox box_total(0, total_uncertainty.first-total_uncertainty.second.first, total.size()/2, total_uncertainty.first+total_uncertainty.second.second);
  box_total.SetLineColor(kGreen);
  box_total.SetFillColor(kGreen);
  box_total.SetFillStyle(3004);
  box_total.Draw();
*/
/*
  double max_man, max_aut, min_man, min_aut;
  min_aut=overall_unc.first-overall_unc.second;
  max_aut=overall_unc.first+overall_unc.second;

  min_aut=overall_unc.first-overall_unc.second;
  max_aut=overall_unc.first+overall_unc.second;
*/
	TLatex texti, textAtlas, textCh;
	texti.SetNDC();
	texti.Draw();
 	texti.DrawLatex(0.65, 0.65, ((std::string)Form("#Delta PDF = %.4f", overall_unc.second)).c_str());

  textAtlas.SetNDC();
  textAtlas.Draw();
  textAtlas.SetTextFont(72);
  textAtlas.SetTextSize(0.035);
  textAtlas.DrawLatex(0.12, 0.85, "ATLAS");
  //textAtlas.DrawLatex(0.12, 0.81, "Work in progress");
  textAtlas.DrawLatex(0.12, 0.81, "Internal");
  
if(channel=="lep")
  textAtlas.DrawLatex(0.12, 0.73, "Leptonic (e+#mu) channel");
else if(channel=="had")
 textAtlas.DrawLatex(0.12, 0.73, "Hadronic (e+#mu) channel");
else if(channel == "lephad")
 textAtlas.DrawLatex(0.12, 0.73, "Combined channel");


  g_MSTW->Draw("P same");
  g_CT10->Draw("P same");

  canni.Print((output_basebath+".png").c_str());
  canni.Print((output_basebath+".pdf").c_str());
  canni.Print((output_basebath+".eps").c_str());


}

std::vector<double> ReadValues(const char * filename)
{
  // define input file 
  std::ifstream inputfile; 

  // open file 
  inputfile.open(filename); 

  // check if file is open 
  if (!inputfile.is_open())
    {
      std::cout << "ReadValues(). File \"" << filename << "\" not found." << std::endl;              
    }

  // reset parameters 
  std::vector<double> fracInVec; 
  std::string name = "empty";
	
	// Read in, in case input files are separeted with "," 
  // read a string via file since long string causes memory error in CINT when it is read via stdin
  std::string argStr;
  std::ifstream ifs(filename);
  std::getline(ifs,argStr);

	// split by ','
  for (size_t i=0,n; i <= argStr.length(); i=n+1)
    {
      n = argStr.find_first_of(',',i);
      if (n == std::string::npos && fracInVec.size()!=0) //end of file
				n = argStr.length();
			else if(n == std::string::npos) //either only one entry or separated by line
        break;
      std::string tmp = argStr.substr(i,n-i);
      fracInVec.push_back(TString(tmp).Atof());
    }
	
	//If there are no "," read in line by line
	if(fracInVec.size() ==0){
  	while(!inputfile.eof())
    	{
      	name = ""; 
      	inputfile >> name;
				if(name.size()==0)
					break;
      	fracInVec.push_back(TString(name).Atof()); 
    	}
	}	
	// close file 
  inputfile.close();    

  // no error 
  return fracInVec; 
}

template<typename T>
T* vector2Array(std::vector<T> vec)
{
        T* array = new T[vec.size()];
        for (int i = 0; i < vec.size(); ++i)
        {
                array[i] = vec[i];
        }
        return array;
}


pair <double, pair <double,double > > hessian_uncertainty (std::vector <double> uncerts) {
    if ( uncerts.size() % 2 != 1){
      cout << "uncerts not uneven" << uncerts.size() << endl;   
    }

    double sum_uncert=0.0;
    for (size_t n=1; n< (uncerts.size()-1) /2 ; n++) {
      //sum_uncert+=pow(fabs(uncerts[n*2]-uncerts[n*2+1]),2);
      sum_uncert+=pow(fabs(uncerts[n*2]-uncerts[n*2-1]),2);
    }

    double uncerty = 0.5 * sqrt(sum_uncert)/1.645;
    pair <double, pair <double,double> > result;
    pair <double, double > uncertainty;
    result.first=uncerts[0];
    uncertainty.first=uncerty;
    uncertainty.second=uncerty;
    result.second=uncertainty;
    return result;

  }

//def hessian_uncertainty(uncerts):
//  if len(uncerts) // 2 != 1:
//    print "uncerts not uneven",uncerts
//  sum_uncert=0.0
//  for n in [ idx for idx in range(len(uncerts)) if idx // 2 == 1  ]:
//    sum_uncert+=(abs(uncerts[n]-uncerts[n+1]))**2
//  uncertainty=0.5*sqrt(sum_uncert)
//  return uncerts[0],uncertainty

pair <double, pair <double, double > > hessian_asymmetric_uncertainty (std::vector <double> uncerts) {
   if ( uncerts.size() % 2 != 1){
     cout << "uncerts not uneven" << uncerts.size() << endl;   
   }

   double sum_uncert_up=0.0;
   double sum_uncert_down=0.0;

   std::cout<<"uncerts.size(); "<<uncerts.size()<<std::endl;
   double x0=uncerts[0];

   for (size_t n = 1; n<= uncerts.size()/2; n++) {
     double error1= uncerts[2*n]-x0;
     double error2= uncerts[2*n-1]-x0;
     sum_uncert_up+=pow( std::max(static_cast<double>(0), std::max(error1,error2)),2);
     sum_uncert_down+=pow(std::max(static_cast<double>(0), std::max(-error1,-error2)),2);
   }   

   pair <double, pair <double,double> > result;
   pair <double, double > uncertainty;
   uncertainty.first=sqrt(sum_uncert_up);
   uncertainty.second=sqrt(sum_uncert_down);

   result.first=uncerts[0];
   result.second=uncertainty;

   return result;
 }

//def hessian_asymmetric_uncertainty(uncerts):
//  if len(uncerts) // 2 != 1:
//    print "uncerts not uneven",uncerts
//  sum_uncert_up=0.0
//  sum_uncert_down=0.0
//  x0=uncerts[0]
//  for n in xrange(1,len(uncerts)):
//    error=uncerts[n]-x0
//    if error > 0:
//      sum_uncert_up+=error**2
//    else:
//      sum_uncert_down+=error**2
//  return uncerts[0],sqrt(sum_uncert_up),sqrt(sum_uncert_down)

pair < double , pair <double, double > > RMS_NNPDF_uncertainty(vector<double> uncerts) {
  if ( uncerts.size() % 2 != 1){
    cout << "uncerts not uneven" << uncerts.size() << endl;    
  }

  double average=0.0;
  for (size_t n=1; n<uncerts.size(); n++){
    average+=uncerts[n];
  }

  average=average / (uncerts.size()-1);
  double std=0.0;
  for (size_t n=1; n<uncerts.size(); n++){
    std+=pow(uncerts[n]-average,2);
  }

  std=sqrt(1.0/(uncerts.size()-1) * std);
  pair <double, pair <double,double> > result;
  pair <double, double > uncertainty;
  result.first=average;

  uncertainty.first=std;
  uncertainty.second=std;
  result.second=uncertainty;

  return result;
}

//def NNPDF_uncertainty(uncerts):
//  if len(uncerts) % 2 != 1:
//    print "uncerts not uneven",uncerts
//  average=sum( uncerts[1:]) / (len(uncerts)-1)
//  std=sqrt(1.0/(len(uncerts)-1.0) * sum ( [ (x-average)**2 for x in uncerts[1:] ] ) )
//  return average, std

pair <double, double > envelope (double central, vector < pair < double , pair <double, double> > > pair_uncerts){
  double ups=-1e10;
  double downs=1e10;
  for (size_t n=0; n< pair_uncerts.size() ; n++){
    double up=pair_uncerts[n].first+pair_uncerts[n].second.first;
    double down=pair_uncerts[n].first-pair_uncerts[n].second.second;
    if (up>ups) ups=up;
    if (down<downs) downs=down;
  }
  pair <double, double> result;
  result.first=(ups+downs)/2.0;
  result.second=(ups-downs)/2.0;
  return result;
}


//def envelope(central,pair_uncerts):
//  ups=-1e10
//  downs=1e10
//  for pair in pair_uncerts:
//    if len(pair )==2:
//      up=pair [0]+pair [1]
//      down=pair [0]-pair [1]
//    elif len(pair )==3:
//      up=pair [0]+pair [1]
//      down=pair [0]-pair [2]
//    else:
//      print "pair ",pair ,"not recognized"
//    if up > ups: ups=up
//    if down < downs : downs=down
//    print up,down, ups,downs
//  return central,(ups-downs)/2.0

pair <double, double> calculate_envelope(double central, vector <double> uncert_hess, vector < double > uncert_asymm, vector < double > uncert_NNPDF ) {
  vector < pair < double, pair < double, double > > > pairs;
  pairs.push_back(hessian_uncertainty(uncert_hess));
  pairs.push_back(hessian_asymmetric_uncertainty(uncert_asymm));
  pairs.push_back(RMS_NNPDF_uncertainty(uncert_NNPDF));
  return envelope(central, pairs);
}

//def calculate_envelope(central,uncert_hess,uncert_asymm,uncert_NNPDF):
//  all_uncerts=[hessian_uncertainty(uncert_hess),hessian_asymmetric_uncertainty(uncert_asymm),NNPDF_uncertainty(uncert_NNPDF)]
//  return envelope(central,all_uncerts)
  
bool doestringcontain(std::string str1, std::string str2){
 bool does = false;
 std::size_t found = str1.find(str2);
 if (found!=std::string::npos) does = true;
  return does;
}

//-----------------------------------

# ifndef __CINT__  // the following code will be invisible for the interpreter
int main(int argc, char **argv)
{
  
  if (!(argc >= 3)) {
      std::cout << "3 parameter needed: ./PDFPlotter [CT10 file] [MSTW file] [NNPDF file]" << std::endl;
      return 1;
  } 

  const char *ch_CT10 = argv[1];
  const char *ch_MSTW = argv[2];
  const char *ch_NNPDF = argv[3];
  //std::vector <double> HERAPDF;
  //if (argc == 5) HERAPDF= ReadValues(argv[4]);
  
  PDFPlotter(ch_CT10, ch_MSTW, ch_NNPDF);
}
# endif

//--------------------------------------
