/* ////////////////////////////////////////////////////////////////////////////////
*	For generate a summary file including all mass results
*	After using BatchMassCal.C
*	Folders like A2_83Ge  ; A2_74Ni; A2_38ArH  were created
*
*	make a txt file that contains the name of folders , ame_mass and ame_err in the 
*	same path as these folders
*
*	this code read the latest files in these folder and generate a summary text file
*//////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <time.h>
#include <algorithm>
#include "TROOT.h"
#include "TSystem.h"


using namespace std;

struct M_uamu{
	double mass;
	double mass_err;
};

struct M_keV{
	double mass;
	double mass_err;
};


string ExtractNewFileName(string PATH = "../A2_74Ni/", string ion_name = "A2_74Ni", string fin_keyword = "", string fout_keyword="");

void Autobatchplot(string PATH ="./", string listname = "including.txt"){

	char ionname_c[1000];
	string ionname;
	double mass_ame;
	double mass_err_ame;
	M_uamu mass_mean_by_m , mass_mean_by_rho;
	M_keV mass_mean_by_m_kev , mass_mean_by_rho_kev;

	FILE* fout =NULL;
	FILE* fin=NULL;

	listname = PATH + listname;

	fin = fopen(listname.c_str(),"r");

	if(fin==NULL){
			cout<<"\e[1;33m"<<"Error!! No list of ion is found!!!!"<<"\e[0m"<<endl;
			return;
	}

	struct tm *ptr_time;
	time_t time_now= time(0);
	ptr_time = localtime(&time_now);
	char time_str[20];
	strftime(time_str,20,"%Y%m%d%H%M%S",ptr_time);

	string foutname = Form("summary_auto_%s.txt",time_str);
	foutname = PATH + foutname;

	fout = fopen(foutname.c_str(),"w");

	fprintf(fout,"ionname\tmassAME\tmass_errAME\tmean_by_massUAMU\tmean_err_by_massUAMU\tmean_by_masskeV\tmean_err_by_masskeV\tmean_by_rhoUAMU\tmean_err_by_rhoUAMU\tmean_by_rhokeV\tmean_err_by_rhokeV\n");

	char file_1st_line[1000];

	fgets(file_1st_line,1000,fin);
	if(file_1st_line[0] == 'A') rewind(fin);
	

	while(!feof(fin)){
		fscanf(fin,"%s\t%lf\t%lf",ionname_c,&mass_ame,&mass_err_ame);

		if(!feof(fin)){
			string commandline = PATH + ionname_c + "/";
			if(gSystem->AccessPathName(commandline.c_str()) == true){
				cout<<"\e[1;31m"<<"Error!!!! folder "<<ionname_c<<" not exist!!!!"<<"\e[0m"<<endl;
				continue;
			}

			string getfilename = 	ExtractNewFileName(commandline,ionname_c,ionname_c,"read_new_");	
			if(strcmp(getfilename.c_str(),"fail") ==0) continue;

			batchplot(commandline,getfilename,mass_ame,mass_err_ame);


//read from plotnote (==> rho_i ==> mass_i ==> mass_mean)
			getfilename = 	ExtractNewFileName(commandline,ionname_c,"plotnote_","read_new_plotnote_");
			if(getfilename == "fail") continue;


			commandline = PATH + ionname_c + "/";
			commandline = commandline + getfilename + ".txt";
			if(gSystem->AccessPathName(commandline.c_str()) == true){
				cout<<"\e[1;31m"<<"Error!!!! plotenote "<<getfilename<<" not exist!!!!"<<"\e[0m"<<endl;
				continue;
			}



			FILE* fread =NULL;

			fread = fopen(commandline.c_str(),"r");

			char buffer[1000];
			
			while(!feof(fread)){
				fgets(buffer,1000,fread);
				if(!feof(fread)){
					string buffer_s = buffer;					
					if(buffer_s.find("Mass Mean[micro_amu]= ") != string::npos){
							string  remove_word = "Mass Mean[micro_amu]= ";
							char cc;
							buffer_s.erase(buffer_s.begin(),buffer_s.begin()+remove_word.size());
							istringstream is(buffer_s);
							is>>mass_mean_by_m.mass>>cc>>mass_mean_by_m.mass_err>>cc;
					}
					if(buffer_s.find("Mass Mean[ME in keV]= ") != string::npos){
							string  remove_word = "Mass Mean[ME in keV]= ";
							char cc;
							buffer_s.erase(buffer_s.begin(),buffer_s.begin()+remove_word.size());
							istringstream is(buffer_s);
							is>>mass_mean_by_m_kev.mass>>cc>>mass_mean_by_m_kev.mass_err>>cc;
					}
				}

			}
			
			fclose(fread);


//read from plotnote (==> rho_i ==> mass_i ==> mass_mean)

			commandline = PATH + ionname_c + "/";
			

			getfilename = 	ExtractNewFileName(commandline,ionname_c,"rho_","read_new_rho_");
			if(getfilename == "fail") continue;

			commandline = commandline +getfilename +".txt";
			if(gSystem->AccessPathName(commandline.c_str()) == true){
				cout<<"\e[1;31m"<<"Error!!!! file "<<getfilename<<" not exist!!!!"<<"\e[0m"<<endl;
				continue;
			}


			fread = fopen(commandline.c_str(),"r");
			while(!feof(fread)){
				fgets(buffer,1000,fread);
				if(!feof(fread)){
					if(strncmp(buffer,"combine",7) ==0){
						string buffer_s = buffer;
						string remove_word = "combine";
						buffer_s.erase(buffer_s.begin(),buffer_s.begin()+remove_word.size());
						istringstream is(buffer_s);
						double temm;
						is>>temm>>temm>>temm>>temm>>temm>>temm>>mass_mean_by_rho.mass;
						is>>mass_mean_by_rho.mass_err>>mass_mean_by_rho_kev.mass>>mass_mean_by_rho_kev.mass_err;
					}
				}

			}
			
			ionname = ionname_c;
			ionname.erase(ionname.begin(),ionname.begin()+ionname.find("_")+1);
			fprintf(fout,"%s\t%.4f\t%.4f\t%.4f\t%.4f",ionname.c_str(),mass_ame,mass_err_ame,mass_mean_by_m.mass,mass_mean_by_m.mass_err);
			fprintf(fout,"\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",mass_mean_by_m_kev.mass,mass_mean_by_m_kev.mass_err,mass_mean_by_rho.mass,mass_mean_by_rho.mass_err,
															mass_mean_by_rho_kev.mass,mass_mean_by_rho_kev.mass_err);


		}


	}//end of read ion list.txt

	fclose(fout);

}



string ExtractNewFileName(string PATH = "../A2_74Ni/", string ion_name = "A2_74Ni", string fin_keyword = "", string fout_keyword=""){

			string commandline = PATH;
			commandline = ".!ls -t "+ commandline + fin_keyword + "* > "  +  commandline + fout_keyword + ion_name + "_XX.txt";

			gROOT->ProcessLine(commandline.c_str());

			commandline = PATH;
			commandline = commandline + fout_keyword + ion_name + "_XX.txt";

			FILE* fmass_data =NULL;

			fmass_data = fopen(commandline.c_str(),"r");

			char buffer[1000];
			fgets(buffer,1000,fmass_data);

			fclose(fmass_data);

			string keyword = buffer;
			if(keyword.find(fin_keyword) == string::npos){
				cout<<"\e[1;31m"<<"Error!!!! file of "<<fin_keyword<<" not exist!!!!"<<"\e[0m"<<endl;
				return "fail";
			}

			
			while(keyword.find(fin_keyword,1) != string::npos){keyword.erase(keyword.begin(),keyword.begin()+keyword.find(fin_keyword,1)); }

			keyword.erase(keyword.begin()+keyword.find(".txt"),keyword.end());


		return keyword;


}