#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include "TMath.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"

using namespace std;

struct Iondva{
	char name[1000];
	char name_latex[1000];
	double mass;
	double mass_err;
	Iondva(){
		for(int i=0;i<1000;i++){
			name[i]='\0';
			name_latex[i]='\0';
		}
	}

};

TH1D* h[2];
THStack *hs=NULL;
Iondva **Ioncontainer=NULL;
TGraphErrors ** ge=NULL;
TMultiGraph *mg =NULL; 
TLegend *lg =NULL;
TCanvas * c1=NULL; 
//int numofmeasure_history=0;  // just for memory clean

double histo_max=0;
double histo_min=0;

void unzoom(){
		histo_max = h[0]->GetMaximum()*1.2;
		histo_min = h[1]->GetMinimum()*1.2;
		h[0]->SetMaximum(histo_max); h[1]->SetMaximum(histo_max);
		h[0]->SetMinimum(histo_min); h[1]->SetMinimum(histo_min);

}

void plotdva(string path = "../",string filename="without .txt"){
		if(strncmp(filename.c_str(),"without .txt",12) == 0){cout<<"\e[1;33m"<<"Please input path and filename!!!"<<"\e[0m"<<endl; return;}
		string infilename = path + filename +".txt";
		ifstream fin;
		fin.open(infilename.c_str(),ios::in);

		if(!fin){
			cerr<<"Can not open file: \""<<filename.c_str()<<"\" !!! abort!"<<endl;
			return;
		}

		char temc[1000];
		bool firstloop_file_done=false;  // no record in temc yet;
		int linenum=0;

		while(!fin.eof()){
			if(firstloop_file_done){linenum++;}
			fin.getline(temc,1000);
			firstloop_file_done=true;
		}

		int NumIons = 0;
		if(linenum>3){	NumIons = linenum-3;}
		else{cout<<"\e[1;33m"<<"find no information of ion, please check!!"<<"\e[0m"<<endl; return;}
		cout<<"num of Ions "<<NumIons<<endl;

		fin.clear();
		fin.seekg(0,ios::beg);
		fin.getline(temc,1000);
cout<<temc<<endl;
		int numofmeasure = -1;   // how many measurement results, not number of ions; compare result by different people; should be >=0 (0 means have ref only)
		fin>>numofmeasure;
		if(numofmeasure<0){ cout<<"no mass info of any ion!! break!"<<endl; return;} // ==0 => only have ref information; 1 means have results from 1 person
cout<<numofmeasure<<endl;
		Ioncontainer = new Iondva*[NumIons];
		for(int i=0;i<NumIons;i++){Ioncontainer[i] = new Iondva[numofmeasure+1];}
		//Iondva (*Ioncontainer)[Ncolumn] = new Iondva[NumIons][Ncolumn];  // int (*p)[column] = new int[row][column] !!!!

		fin.getline(temc,1000);  // >> ahead only read number ; there is a '\n' remain => thus have to read two time!!!!!!!
		fin.getline(temc,1000);
		int rowindex=0,colindex=0;
		cout<<temc<<endl;

		while(!fin.eof()){
			if(rowindex<NumIons){
				char name_ion[1000];
				char name_ion_latex[1000];
				fin>>name_ion>>name_ion_latex;
				for(colindex=0;colindex<numofmeasure+1;colindex++){
					strcpy(Ioncontainer[rowindex][colindex].name,name_ion);
					strcat(Ioncontainer[rowindex][colindex].name_latex,"  ");
					strcat(Ioncontainer[rowindex][colindex].name_latex,name_ion_latex);
					fin>>Ioncontainer[rowindex][colindex].mass>>Ioncontainer[rowindex][colindex].mass_err;
				}
			}
			else{break;}
			rowindex++;
		}// end of while loop to read data

		if(ge!=NULL){delete ge;delete mg;delete lg;}
		ge = new TGraphErrors*[numofmeasure];
		mg = new TMultiGraph("mg","mg");
		lg = new TLegend(0.75,0.7,0.9,0.9);

		for(int i=0;i<numofmeasure;i++){  // load data to TGraph  // column
			ge[i] = new TGraphErrors();
			for(int ionIndex=0;ionIndex<NumIons;ionIndex++){ // row
				ge[i]->SetPoint(ionIndex,ionIndex+1,Ioncontainer[ionIndex][i+1].mass-Ioncontainer[ionIndex][0].mass);
				ge[i]->SetPointError(ionIndex,0,Ioncontainer[ionIndex][i+1].mass_err);
			}
		}

		ge[0]->SetMarkerStyle(20); 
		ge[0]->GetXaxis()->SetNdivisions(0);
		ge[0]->GetXaxis()->SetRangeUser(0,NumIons);

		if(numofmeasure==2){
		ge[1]->SetMarkerStyle(24);
		ge[1]->SetMarkerColor(kRed);
		ge[1]->SetLineColor(kRed); 
		ge[1]->GetXaxis()->SetNdivisions(0);
		ge[1]->GetXaxis()->SetRangeUser(0,NumIons);}

		ge[0]->SetDrawOption("AP");
		if(numofmeasure==2)ge[1]->SetDrawOption("AP");
		mg->Add(ge[0],"pe");
		if(numofmeasure==2)mg->Add(ge[1],"pe");
		//mg->GetXaxis()->SetNdivisions(0);
		lg->AddEntry(ge[0],"MRTOF","lep");
		if(numofmeasure==2)lg->AddEntry(ge[1],"S.Giraud","lep");

		if(h[0]!=NULL){ delete h[0]; delete h[1];}
		h[0] = new TH1D("hu","",NumIons,0.5,NumIons+0.5);  // 0.5 for alignment
		h[1] = new TH1D("hd","",NumIons,0.5,NumIons+0.5);

		for(int i=0;i<NumIons;i++){
			h[0]->SetBinContent(i+1,Ioncontainer[i][0].mass_err);
			h[1]->SetBinContent(i+1,-Ioncontainer[i][0].mass_err);
			h[0]->GetXaxis()->SetBinLabel(i+1,Ioncontainer[i][0].name_latex);
			h[1]->GetXaxis()->SetBinLabel(i+1,Ioncontainer[i][0].name_latex);
		}

		histo_max = h[0]->GetMaximum()*1.2*2.1*20;
		histo_min = h[1]->GetMinimum()*1.2*2.1*20;
		h[0]->SetMaximum(histo_max); h[1]->SetMaximum(histo_max);
		h[0]->SetMinimum(histo_min); h[1]->SetMinimum(histo_min);

		//h[0]->GetYaxis()->SetTitle("ME_{MRTOF} - ME_{AME 2020} [keV]");
		//h[0]->GetYaxis()->CenterTitle();

		h[0]->SetFillColor(kGreen-10);h[0]->SetLineColor(kBlue-7);
		h[1]->SetFillColor(kGreen-10);h[1]->SetLineColor(kBlue-7);

		gStyle->SetOptStat(0);

		if(c1!=NULL)delete c1;
		c1 = new TCanvas("c1","c1",1200,800);
		c1->SetFillColor(0);
   		c1->SetBorderMode(0);
   		c1->SetBorderSize(2);
   		c1->SetFrameBorderMode(0);
   		c1->SetFrameBorderMode(0);

		hs = new THStack("hs","");
		//hs->GetYaxis()->SetTitle("ME_{MRTOF} - ME_{AME 2020} [keV]");
		hs->Add(h[0],"][");
		hs->Add(h[1],"][");
		hs->SetMaximum(histo_max);
		hs->SetMinimum(histo_min);
		//hs->GetYaxis()->SetTitle("ME_{MRTOF} - ME_{AME 2020} [keV]");// why Latex cause crash here??
		//hs->GetYaxis()->CenterTitle();
		hs->Draw("nostack");
		hs->GetYaxis()->SetTitle("ME_{MRTOF} - ME_{AME 2020} [keV]");// why Latex cause crash here??
		hs->GetYaxis()->CenterTitle();
		mg->Draw("P");
		lg->AddEntry(h[0],"AME 2020","lf");
		lg->Draw();

}
