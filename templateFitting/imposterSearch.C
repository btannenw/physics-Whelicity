#include <TFile.h>
#include <TH1.h>


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

void imposterSearch()
{
  std::vector<std::string> list_of_files = ReadInputFiles("ls1_sf.txt");


  for(int i=0;i<list_of_files.size();i++)
    //for(int i=0;i<10;i++)
    {
      if(list_of_files.at(i).find("nominal") != std::string::npos)
      {
        TFile* fF=new TFile(("/afs/cern.ch/work/m/mkareem/Wpol_mc/RunPlotFactory/OutputHistos_SFsyst_lhCut/NormalisedHistos_4incl_Nominal_1exclTags/"+list_of_files.at(i)).c_str(),"READ");
        //std::cout<<"file: "<<list_of_files.at(i) << std::endl;
      TH1D* h1=(TH1D*)fF->Get("NBtag");
      
      if(h1->GetBinContent(3)!=0)
      std::cout<<list_of_files.at(i).c_str()<<" IS WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      
      delete h1;
      fF->Close();
      delete fF;
      }
        
    }
}

