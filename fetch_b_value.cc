#include "TROOT.h"
#include "TRint.h"
#include"TSystem.h"
#include <stdlib.h>
#include <string>
#include <vector>


using namespace std;


void fetch_b_value(string PATH ="../upper_folder_of_Run/", int iRun_begin=0, int iRun_end=1){
			for(int irun=iRun_begin; (irun>=iRun_begin) && (irun<=iRun_end);irun++){
				string commandline = Form(".!echo Run%d >> fetch_b_value.txt",irun);
				gROOT->ProcessLine(commandline.c_str());
				commandline = Form("Run%d/LST/*@*.lst",irun);
				commandline = ".!grep 'arbraxbxsr:' "+  PATH + commandline +" >> fetch_b_value.txt";
				gROOT->ProcessLine(commandline.c_str());
			}

}