#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"

//#include "TProfile.h"
//#include "TAxis.h"

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "TGProgressBar.h"


#if defined (MAKECINT)
#pragma link C++ class vector<Long_t>+;   // in order to identify vector<Long_t> 
#endif


using namespace std;
// default read and write file from and to "In/OutPATH"/rootfiles/ "filename".   no "rootfiles" need in path!!!!  READ "_dcorrect" file !!!!!
void CombineTree(string InPATH="",string Infilename="",string OutPATH="",string Outfilename="", double in_tof_ref=0, double tof_ref_std=0, double t0=0){
      string inputfile = InPATH + "rootfiles/" + Infilename + ".root";
      TFile *fin = new TFile(inputfile.c_str());
      if(fin->IsZombie()){
      	cout<<"\e[1;33m Can not open file:  "<<inputfile<<"\e[0m"<<endl;
      	return;
      }

      TTree *intree = (TTree *)fin->Get("tree0");

       cout<<"\e[1;33m reading data from "<<inputfile<<"\e[0m"<<endl;

          Long64_t nevt = 0;
 //     Long64_t num = 0;
	      int nhits[10];
	//      int glo = 0;
	      vector <Int_t> *channel = new vector <Int_t>();
	      vector <Int_t> *edge = new vector <Int_t>();
	      vector <Int_t> *value = new vector <Int_t>();
	      vector <Double_t> *time = new vector <Double_t>();
	      vector <Double_t> *timec = new vector <Double_t>();
	      vector <Int_t> *sweeps = new vector <Int_t>();
	      vector <Long_t> * sweeps_global = new vector <Long_t>;
	      vector <Int_t> *glo = new vector <Int_t>();
	      vector <Int_t> *tag = new vector <Int_t>();
	      vector <Int_t> *lost = new vector <Int_t>();
	      double offset_newtree = 0;
	      
	      intree->SetBranchAddress("nevt", &nevt); 
	//      intree->SetBranchAddress("nevt", &num); 
	      intree->SetBranchAddress("nhits", nhits);
	      intree->SetBranchAddress("channel", &channel);
	      intree->SetBranchAddress("edge", &edge);
	      intree->SetBranchAddress("value", &value);
	      intree->SetBranchAddress("time", &time);
	      intree->SetBranchAddress("timec", &timec);
	      intree->SetBranchAddress("sweeps", &sweeps);
	      intree->SetBranchAddress("sweeps_global", &sweeps_global);
	      intree->SetBranchAddress("glo", &glo);
	      intree->SetBranchAddress("tag", &tag);
	      intree->SetBranchAddress("lost", &lost);
	      intree->SetBranchAddress("offset", &offset_newtree);

	  Long_t glo_sweep_old =0;  // get the global sweep histry
	  Long64_t nevt_old = 0;
	  TTree * outtree = NULL;
      string outputfile = OutPATH + "rootfiles/" + Outfilename + ".root";
      TFile *fout = new TFile(outputfile.c_str());
      if(fout->IsZombie()){
      	  cout<<"\e[1;33m Can not open file:  "<<outputfile<<" ; Create a new one!!"<<"\e[0m"<<endl;
      	  //delete fout;
      	  fout = new TFile(outputfile.c_str(),"UPDATE");
      	  outtree = new TTree("tree0","tree");
      	  outtree->Branch("nevt",&nevt);
		//  tree->Branch("nevt1",&num);
		  outtree->Branch("nhits",nhits,"nhits[10]/I");
		  outtree->Branch("channel",&channel);
		  outtree->Branch("edge",&edge);
		  outtree->Branch("value",&value); // in unit of 100ps
		  outtree->Branch("time",&time); // in unit of ns
		  outtree->Branch("timec",&timec); // in unit of ns
		  outtree->Branch("sweeps",&sweeps);
		  outtree->Branch("sweeps_global",&sweeps_global);
		  outtree->Branch("glo",&glo);
		  outtree->Branch("tag",&tag);
		  outtree->Branch("lost",&lost);
		  outtree->Branch("offset",&offset_newtree);
      }
      else{
      	    cout<<"\e[1;33m  open combined tree file:  "<<outputfile<<"\e[0m"<<endl;
      	    delete fout;
      	    fout = new TFile(outputfile.c_str(),"UPDATE");
      	    outtree = (TTree*) fout->Get("tree0");
	        outtree->SetBranchAddress("nevt", &nevt); 
	//      intree->SetBranchAddress("nevt", &num); 
	        outtree->SetBranchAddress("nhits", nhits);
	        outtree->SetBranchAddress("channel", &channel);
	        outtree->SetBranchAddress("edge", &edge);
	        outtree->SetBranchAddress("value", &value);
 	        outtree->SetBranchAddress("time", &time);
	        outtree->SetBranchAddress("timec", &timec);
	        outtree->SetBranchAddress("sweeps", &sweeps);
	        outtree->SetBranchAddress("sweeps_global", &sweeps_global);
	        outtree->SetBranchAddress("glo", &glo);
	        outtree->SetBranchAddress("tag", &tag);
	        outtree->SetBranchAddress("lost", &lost);
	        outtree->SetBranchAddress("offset", &offset_newtree);

	        Long64_t Nentry_oldtree = outtree->GetEntriesFast();
	        outtree->GetEntry(Nentry_oldtree-1);
	        glo_sweep_old = sweeps_global->at(0);
	        nevt_old = nevt;

	        channel->clear();
		    edge->clear();

		    value->clear();
		    time->clear();
		    timec->clear();
		    sweeps->clear();
		    sweeps_global->clear();
		    glo->clear();
		    tag ->clear();
		    lost->clear(); 
		    for(int i=0; i<10; i++) nhits[i] = 0;

      }

      double offset = (in_tof_ref -t0) / (tof_ref_std-t0);

      Long64_t nentries = intree->GetEntriesFast();
  	  cout<<"infile nentries = "<<nentries<<endl;

  	  int counter =0;

  	  for(Long64_t jentry=0; jentry<nentries; jentry++){
    		if(jentry == (int)counter*nentries/20){  // divide into 20 piese , each piese is equivalent to 5%
      		cout<< '\r' << "finish running : "<< counter*5 << "%" <<flush;
      		counter++;
     		}

     		intree->GetEntry(jentry);
     		nevt = nevt + nevt_old;
     		offset_newtree = offset_newtree * offset;

		    for(int j=0; j<((int)sweeps_global->size()); j++){
		      double current_sweep = sweeps_global->at(j);
		      sweeps_global->at(j) = current_sweep + glo_sweep_old;    // to make glo_sweep keep starting from the glo_sweep of last file
		      
		    }
     		     // calculate new timec
		    for(int j=0; j<((int)timec->size()); j++){
		      double itime = timec->at(j);
		      timec->at(j)=((itime-t0) / offset + t0);    // ratio method correct
		      // cout<<"offset = "<<offset<<", timec = "<<itime-offset<<endl;
		    }

		    outtree->Fill();
    
		    channel->clear();
		    edge->clear();

		    value->clear();
		    time->clear();
		    timec->clear();
		    sweeps->clear();
		    sweeps_global->clear();
		    glo->clear();
		    tag ->clear();
		    lost->clear(); 
		    for(int i=0; i<10; i++) nhits[i] = 0;
		}

		  cout<<endl;

	    cout<<"\e[1;33m writing tree to "<<outputfile<<"\e[0m"<<endl;
	    fout->cd();
	    outtree->Write();
	    fout->Close();
	    fin->Close();


}


//double tof_ref[]={6328290.1,6328283.5,6328283.9,6328281.4,6328286.602,6328282.5,6328283.9,6328282.7,6328287.6,6328283.4,6328289.1,
//	6328287.3,6328278.148,6328285,6328287.5,6328288.369,6328285.1,6328286.5,6328284.941,6328292.1,6328285.2,6328291.548,6328290.5,6328293.3};

double tof_ref[]={12442252.4431,12442244.4261,12442251.3291,12442251.7318,12442265.4788};

void autorun(){

		//int filenum[]={52,53,56,57,58,59,60,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79};
		int filenum[]={0,1,2,3,4,5};
		double t_0 = 130;
		for(int fileindex=0;fileindex<5;fileindex++){
			string infilename=Form("Ni_078Zn@690_107-S%02d_dcorrect",filenum[fileindex]);
			//string infilename=Form("Hi_088As@600_%03d_dcorrect",filenum[fileindex]);
			string outfilename="Ni_078Zn@690_107-Scan_dcorrect";
			string inpath= "/home/xian/winfoder/MT2021/211201MT_ana/Run107-Scan/";
			string outpath="/home/xian/winfoder/MT2021/211201MT_ana/Run107-Scan/";

			CombineTree(inpath,infilename,outpath,outfilename,tof_ref[fileindex],tof_ref[0],t_0);

		}



}

















