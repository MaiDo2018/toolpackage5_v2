#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#endif

void help(int choice=-1);
void welcome();
void rootlogon() {
/*
  //Using the RootTools
  //gROOT->ProcessLine(".L ~/ROOT/Macros/Utilities/RootTools.C+");

  //Base Style
  gROOT->SetStyle("Plain");
  // gROOT->SetStyle("Modern");
    // //gROOT->SetStyle("Classic");

  //Force Style
  gStyle->SetHistFillColor(7);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistLineColor(kBlue);
  gStyle->SetFuncColor(kRed);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetGridWidth(2);
  //gStyle->SetPadGridX(1);
  //gStyle->SetPadGridY(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(0);
  gStyle->SetStatX(0.9);  
  gStyle->SetStatY(0.9);  
  gStyle->SetPalette(1);
  gStyle->SetOptLogz(1);
  //  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1111111);
  gStyle->SetPadBorderMode(1);
  //gStyle->SetOptDate(1); */

  gStyle->SetOptStat(1001111);
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTitleFont(132,"XYZ");
 // gStyle->SetTitleFont(132,"");
  //gStyle->SetTextFont(132);
  gStyle->SetStatFont(132);
  /*
  TColor *dummyColor = new TColor();
  UInt_t Number = 4;
  Double_t Red[4]   = { 0.50, 0.00, 0.99, 0.99};
  Double_t Green[4] = { 0.99, 0.00, 0.00, 0.99};
  Double_t Blue[4]  = { 0.99, 0.99, 0.00, 0.00};
  Double_t Stops[4] = { 0.00, 0.30, 0.70, 0.99};
  dummyColor->CreateGradientColorTable(Number,Stops,Red,Green,Blue,100);

  gROOT->ProcessLine(".L Analysis.C+");
*/


   welcome();

}

void welcome(){
	std::cout<<"**************************************"<<std::endl;
	std::cout<<"* 		To analyse mass      *"<<std::endl;
	std::cout<<"*		   Run	      	     *"<<std::endl;
	std::cout<<"*		.L preview3.C+	     *"<<std::endl;
	std::cout<<"*   for beta-tof .L preview4.C+	     *"<<std::endl;
	std::cout<<"*   for more information => help()   *"<<std::endl;
	std::cout<<"**************************************"<<std::endl; 
}

void help(int choice){
     if(choice==-1){
	cout<<" there is some section:"<<endl;
	cout<<"0:header"<<endl;
	cout<<"1:define ion"<<endl;
	cout<<"2:display"<<endl;
	cout<<"3:fitting"<<endl;
	cout<<"4:mass searcher & mass calculate"<<endl;
	cout<<"5:setting of ref or X"<<endl;
	cout<<"6:Mark TOF and ROI"<<endl;
	cout<<"7:multicurves class"<<endl;
	cout<<"8:beta-tof"<<endl;
	cout<<"9:other"<<endl;
	cout<<"10:interactoin mode"<<endl;
	cout<<"11:general procedures"<<endl;
	return;
      }

	// reset,bgRed,bWhite,bBlue,bGreen,bYellow
	string color[]={"\e[0m","\e[1;47m","\e[1;37m","\e[1;36m","\e[1;32m","\e[1;32m"};

	const int Nkey1=12;
 	string keyword1[Nkey1]={"header:","define ion:","display:","fitting:","mass searcher & mass calculate:","setting of ref or X:","Mark TOF and ROI:","multicurves class:","beta-TOF:","others:","interaction mode:","general procedures:"}; // color[5]

	if(choice>=12){Error("choice number","should be less than 11");return;}

	const int Nkey2=4;
	string keyword2[Nkey2]={"Basic","function","Key","procedure"}; //color[3]

	string keyword3="warning"; //color[1]

	string keyword4="%"; //color[2]

	string keyword5="trick"; //color[4]

	bool GETKEY1=false;


	char buffer[1000];
	ifstream fin;
	fin.open("manual.txt",ios::in);
	if(!fin.is_open()){Error("open manual.txt",".txt file unavailable"); return;}


	while(fin.peek()!=EOF){
		fin.getline(buffer,1000);
		string information = buffer;
		bool IsPrint=false;

		if(!GETKEY1){
			  if(information.find(keyword1[choice].c_str())!=string::npos){
				GETKEY1=true;
				cout<<color[5].c_str()<<keyword1[choice].c_str()<<color[0].c_str()<<endl;
				continue;
			   }
			
		}
		else{
			for(int ikey2=0;ikey2<Nkey2;ikey2++){
				if(information.find(keyword2[ikey2].c_str())!=string::npos){
				   cout<<color[3].c_str()<<information.c_str()<<color[0].c_str()<<endl;
					IsPrint=true;
					break;
				}
			}

			if(information.find(keyword3.c_str())!=string::npos){
				if(!IsPrint){
					cout<<color[1].c_str()<<information.c_str()<<color[0].c_str()<<endl;
					IsPrint=true;
					continue;
				}
			}
			else if(information.find(keyword4.c_str())!=string::npos){
				if(!IsPrint){
					cout<<color[2].c_str()<<information.c_str()<<color[0].c_str()<<endl;
					IsPrint=true;
					continue;
				}
			}
			else if(information.find(keyword5.c_str())!=string::npos){
				if(!IsPrint){
					cout<<color[4].c_str()<<information.c_str()<<color[0].c_str()<<endl;
					IsPrint=true;
					continue;
				}
			}
			else if(information.find("<end>")!=string::npos){break;}
			else{
				if(!IsPrint){cout<<information.c_str()<<endl;}
			}
		}
	}//end of while

	fin.close();

}
