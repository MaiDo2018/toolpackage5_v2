#include"TROOT.h"
#include<iostream>
#include<stdlib.h>
#include"TSystem.h"
#include"TTree.h"
#include"TFile.h"
#include<fstream>
#include <vector>
#include <stdio.h>
#include"TMath.h"
#include"TTree.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TH1.h"
#include"TH2.h"
#include"TPaveStats.h"

#include"TStyle.h"
#include"TColor.h"
#include"TExec.h"
#include"TCutG.h"
#include"signal.h"

#include"RefineBeta2.h"

#if defined (MAKECINT)
#pragma link C++ class vector<Long_t>+;   // in order to identify vector<Long_t> 
#endif

using namespace std;

// Graph of beta-beta energy distribution and coincidence time difference distribution
TFile* f_beta_tree=NULL;
TTree* tbeta_beta=NULL;
TCanvas* c_beta_beta=NULL, *c_beta_coin_and_raw=NULL;
TH1D* histo_beta_deltaT=NULL; // time difference of coincidence
TH2D* histo_beta_EkeV=NULL;   // E distribution of beta coincidence
TH2D* histo_beta_EkeV_withcut=NULL; // E distribution with cut
TH2D* histo_beta_EkeV_CutAndRej=NULL; // E distribution with cut and reject
TH1D* histo_h1E=NULL;  // Silicon 1 ADC histo
TH1D* histo_h2E=NULL;  // Silicon2 ADC histo

TExec * exbcolor = new TExec("exbcolor","gStyle->SetPalette(kPastel);");
TExec * exdefaultcolor = new TExec("exdefaultcolor","gStyle->SetPalette(kBird);");

TH1D* h1_t_relate_raw=NULL,*h1_t_relate_coin=NULL,*h2_t_relate_raw=NULL,*h2_t_relate_coin=NULL;
TH2D* h1_t_relate_EkeV_raw=NULL, *h2_t_relate_EkeV_raw=NULL, *h1_t_relate_EkeV_coin=NULL, *h2_t_relate_EkeV_coin=NULL, *h_LogEratio_delta_t=NULL;

bool BetaFastMode=false; // load beta file basing on mrtof; align beta file to first trigger of mrtof
bool alphaMode=false; // Just take the first Si signal for beta-tof coincident measurement; Set B2 Info same as B1 in B1B2 coincident list 
bool ReadBetaOnline = true;
string OnlineBetaFilePath="/mnt/beta_driver/F11_MRTOF/240619MT/usb/";
bool BetaBetaFilter=false; // false-> keep all beta-beta coin. event in time gate; true-> keep only the beta-beta combination with time different approaching given time different value "b2tob1_standard" 

string MainFilePath; // the mother path of LST and rootfiles; for saving cut.root file

bool NeedToCompare=true;  // in BetaFastMode, if there is no new beta.lst file read, not need to compare again.


//%%%%%%%%%%%%%% save beta.lst as a tree %%%%%%%%%%%%%%%%%%
bool MakeBetaRawTree = false;
TTree* tbeta_Raw = NULL;
TFile* f_betaRaw_out = NULL;
//%%%%%%%%%%%%% end of saving beta.lst as a tree %%%%%%%%%%%%%


//%%%%%%%%%%%%% save global time of each sweeps in a tree %%%%%%%%%%%%%%
bool BetaSweepsTimeReady=false;
TTree* tbeta_sweeps_gtime = NULL;
TFile* f_sweeps_gtime_out = NULL;
//%%%%%%%%%%%%% end of saving global time of each sweep in a tree %%%%%%%%%



// beta energy calibration
//**********************************
double beta_slope[2]={0};
double beta_intercept[2]={0};

void GetBeta_E_Calibrate_para(){ // [0]=>Silicon1; [1]Silicon2
	printf("Silicon1: beta_slope[0]= %.2f ; beta_intercept[0]=%.2f\n",beta_slope[0],beta_intercept[0]);
	printf("Silicon2: beta_slope[1]= %.2f ; beta_intercept[1]=%.2f\n",beta_slope[1],beta_intercept[1]);
}

void SetBeta_E_Calibrate_para(double _slope1=1,double _intercept1=0,double _slope2=1,double _intercept2=0){ // [0]=>Silicon1; [1]Silicon2
	if(beta_slope[0]!=_slope1) beta_slope[0] =_slope1;
	if(beta_slope[1]!=_slope2) beta_slope[1] =_slope2;
	if(beta_intercept[0]!=_intercept1) beta_intercept[0] =_intercept1;
	if(beta_intercept[1]!=_intercept2) beta_intercept[1] =_intercept2;
	GetBeta_E_Calibrate_para();
}

double ADC2keV_1(int _adc){return beta_slope[0]*_adc + beta_intercept[0];}
double ADC2keV_2(int _adc){return beta_slope[1]*_adc + beta_intercept[1];}

double keV2ADC_1(double _keV){return (_keV-beta_intercept[0])/beta_slope[0];}
double keV2ADC_2(double _keV){return (_keV-beta_intercept[1])/beta_slope[1];}
//**************************************


//***************** Beta-beta and betaTof coincident condition *************
int ShapingTime[2]={500,500}; // unit in ns
int MCAThreshold[2]={500,110};

int ADC_max_ch=16384;
Long64_t gate_time=500;   // time window 500 ns
Long64_t Tgate_cent_offset=0;
int gate_adc_low_1=0,gate_adc_low_2=0;  // ADC1 and ADC2 low limit
int gate_adc_hi_1=ADC_max_ch,gate_adc_hi_2=ADC_max_ch;  // ADC1 and ADC2 High limit

void ShowCoinCondition(){
	printf("coincident condition:  gate time\t gate_adc_low_1 \t gate_adc_hi_1 \t gate_adc_low_2 \t gate_adc_hi_2 \t Tgate_cent_offset\n");
	cout<<"\e[1;37m"<<gate_time<<"\t"<<gate_adc_low_1<<"\t"<<gate_adc_hi_1<<"\t"<<gate_adc_low_2<<"\t"<<gate_adc_hi_2<<"\t"<<Tgate_cent_offset<<endl;
	printf("For changing condition:     ");
	cout<<"\e[1;37m"<<"SetCoinCondition("<<"\e[0m"<<endl;
	cout<<endl;
}
void SetCoinCondition(Long64_t _gate_time=500, int _gate_adc_low_1=0, int _gate_adc_hi_1=ADC_max_ch, int _gate_adc_low_2=0, int _gate_adc_hi_2=ADC_max_ch,Long64_t _Tgate_cent_offset=0){
	gate_time= _gate_time;   // time window 500 ns
 	gate_adc_low_1= _gate_adc_low_1;         // ADC1 low limit
 	gate_adc_hi_1= _gate_adc_hi_1; 
 	gate_adc_low_2= _gate_adc_low_2;         // ADC1 low limit
 	gate_adc_hi_2= _gate_adc_hi_2; 
 	Tgate_cent_offset = _Tgate_cent_offset;// only for better display
 	ShowCoinCondition();
}


Long64_t ModuleDelay = 50; // 50[ns] signl to ch4 has a 50 ns delay from synchronizer output !!!!!!!!!!!!!!!!!!!!! Becase of using module to convert NIM to TTL
Long64_t RIHalflive = 150 * 1000000;  // radioactive ion half live in ns
void ShowRIHalflife(){
	cout<<"\e[1;37m"<<"RI halflife:  "<<RIHalflive<<" [ns] !!!!"<<"\e[0m"<<endl;
	printf("Set halflife by:  SetRIHalflife()\n");
	cout<<endl;
}
void SetRIHalflife(double InMilliSecond = 150){
		RIHalflive = (Long64_t) InMilliSecond * 1000000;  // to [ns]
		ShowRIHalflife();
}
double Lifetime2Halflife(double Inlivetime_ms=150){
		double  output= Inlivetime_ms*TMath::Log(2);
		cout<<"Halflife T1/2 = "<< output<<" [ms] !!!"<<endl;
		return output;    // in [ms]
}
double Halflife2Lifetime(double Inhalflife_ms=150){
		double output= Inhalflife_ms/TMath::Log(2);
		cout<<"life time tao = "<<output<<" [ms] !!!"<<endl;
		return output;	// in [ms]
}
void SetRIHalflife_byLifeTime(double Inlivetime_ms=150){
	SetRIHalflife(Lifetime2Halflife(Inlivetime_ms));
}


Long64_t CorrectedBetaTime(Long64_t _BetaTime, int _ADC, int _channel){
	if(_channel!=1 && _channel!=2) return _BetaTime;

	return _BetaTime + TMath::Nint(TMath::Sqrt(TMath::Log(_ADC*1.0/MCAThreshold[_channel-1]) *2*ShapingTime[_channel-1]*ShapingTime[_channel-1]));
}


//**************************************************


//************ vasatile make cut for beta event ************************
// use it to make cut and to know the input time relative and Ekev is noise or not
// inside any cut means noise
// set GetCutStatus = true    to know whether there is any cut available or not
// by default, start a guide to make new cuts in c_beta_coin_and_raw Canvas
// input the time_relative and Ekev of an event from specific channel to know it is noise or not

int CutModifyCounter=0; // makecut and load cut will change this number; --> when number is change; allow to run "Compare_Beta_Beta()""
bool MakeBetaCut_Or_FindSignal(Long64_t _TimeRelative =-1000, double _EkeV =-1000, int WhichChannel=1, bool GetCutStatus = false, TFile* RootfilePtr=nullptr, int Clear0_Recreate1_Add2=1); 
// 0->just clear cuts, 1 ->clear and create new cut, 2->add cut


int syn_CH = 4; // beta-tof synchronizer channel
Long64_t SweepInterval = 50000040;

//vector <Long64_t> *nevt_v = new vector <Long64_t>();
vector<Long64_t> *sweeps_gclock_v = new vector<Long64_t>();
vector <Long64_t> *time_v = new vector <Long64_t>();
vector <Long64_t> *time_relative_v = new vector <Long64_t>();
vector<int> *adc_v = new vector<int>();
vector<int> *channel_v = new vector<int>();

Long64_t nevt_bta=0;  // how many lines or events in a file

string Path_beta_file="./";
string Current_beta_file="---"; // when Read_Beta_lst() read new bunch of beta files with "loadmore" = false; vector for storing raw data from file
										// is cleared; update the file name in this case.
										//
string Current_First_beta_file="!!!";  // only update by bath reading function in normal mode
string Current_First_beta_file_FastMode="!!!"; // only update by bath reading function in BetaFast mode
string ToLoad_beta_file="+++";

int PreSetADC1_low_Thres_ch=100;
int PreSetADC2_low_Thres_ch=100;

void LoadNewBetaFile(string filename){
	ToLoad_beta_file = filename;
	cout<<"\e[1;33m"<<"Beta file to be load: "<<filename.c_str()<<"\e[0m"<<endl;
}


double tof_b2l[]={16.5e6,17.e6,17.31e6,17.8e6};
double tof_b2h[]={10.151e6,10.653e6,12.433e6,12.935e6};

void ReadBinaryData(FILE* _fin, Long64_t & _time, int & _channel, int & _adc){
	unsigned char Onebyte[8];
	fread(Onebyte,sizeof(unsigned char),8,_fin);
	unsigned long long timedata=0;
	for(int i=0;i<6;i++){
		unsigned long long tembyte = Onebyte[i];
		tembyte = tembyte<<(8*(5-i));
		timedata = timedata | tembyte;
	}
	
	_time= (Long64_t) timedata*40;
	_channel = ((Onebyte[6] & 0xc0) >> 6) +1;
	unsigned int temadc = Onebyte[6];
	_adc = (int)((temadc & 0x3F) << 8) | Onebyte[7];

}


//%%%%%%%%%%%%%%%% identify Beta event from a new beta file that never read in RAM %%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%      to produce accumulated spectrum over different beta file in FastBetaMode           %%%%%%%%%%%%%%
//  avoid to load repetitive data to histogram
vector<Long64_t> FileStartTime;// store the start time of each file ; reference in Read_Beta_lst_batch();
bool IsBetaFromNewFile(Long64_t _B1time=0, Long64_t _B2time=0, bool UpdateRecord=false, bool ClearRecord = false){ // UpdateRecord -> update the index history of file that have been read, 
																																 // input _B1time and _B2time will be regarded as event already processed basing on the index history
																																// if _B1time and _B2time belong to the file read already in index history --> return false for reject
																																// otherwise --> true ,  allow to load in the histogram
		static bool * FileStatus = NULL;
		static vector<int> FileIndex_tmp;

		if(ClearRecord){
			if(FileStatus!=NULL){
				delete[] FileStatus;
				FileStatus = NULL;
			}

			FileIndex_tmp.clear();
			return true;
		}

		static int FileAmount=0;

		if(FileStatus == NULL){
			if(FileStartTime.size()>0){
				FileAmount = (int)FileStartTime.size();
				FileStatus = new bool[FileAmount];

				for(int indexFile=0; indexFile<FileAmount; indexFile++){
						FileStatus[indexFile]=true;  // true --> regarded as new file, data allowed to be fill in histogram
				}
			}
			else return true; // no record of FileStartTime;  
		}

		if(UpdateRecord && FileIndex_tmp.size()>0){// update file index record from FileIndex_tmp to FileStatus
			for(int i=0;i<(int)FileIndex_tmp.size();i++){ FileStatus[ FileIndex_tmp[i] ] = false;}
			FileIndex_tmp.clear();
			return true;
		}


		int CurrentIndexB1=-1, CurrentIndexB2=-1;

		for(int i=FileAmount-1;i>=0;i--){
			if(CurrentIndexB1 == -1){
				if(FileStartTime[i] <= _B1time){
					 CurrentIndexB1 = i;
					 bool NewB1Index=true;
					 for(int j=0;j<(int)FileIndex_tmp.size();j++){
					 	if(FileIndex_tmp[j] == CurrentIndexB1){ NewB1Index=false; break; }
					 }
					 if(NewB1Index) FileIndex_tmp.push_back(CurrentIndexB1);
				}
			}
				
			if(CurrentIndexB2 == -1){
				if(FileStartTime[i] <= _B2time){
					CurrentIndexB2 = i;
					bool NewB2Index=true;
					 for(int j=0;j<(int)FileIndex_tmp.size();j++){
					 	if(FileIndex_tmp[j] == CurrentIndexB2){ NewB2Index=false; break; }
					 }
					 if(NewB2Index) FileIndex_tmp.push_back(CurrentIndexB2);
				}
				
			}

			if(CurrentIndexB1 != -1 && CurrentIndexB2 != -1) break;

		}


		
		if(CurrentIndexB1 != -1 && CurrentIndexB2 != -1) return FileStatus[CurrentIndexB1] || FileStatus[CurrentIndexB2]; // either B1 or B2 in new file should be filled in histogram
		else if(CurrentIndexB1 != -1) return FileStatus[CurrentIndexB1];
		else if(CurrentIndexB2 != -1) return FileStatus[CurrentIndexB2];

		return true;

}



//********************** How to use Read_Beta_lst() **********************************************
// 1. Read a single file according to the PATH and filename: PATH + filename
// 2. To read a new data or a new batch of data, set loadmore = false.  Raw data already exist in vector will be cleared first.  To continue reading next file ==> set loadmore = true;
// 3. set ADC threshold for ch1 and ch2 in the unite of channel, default value is 100 channel while the full channel is 16383 ch.
// 4. In " BetaFastMode ", global time of all beta events will align to the first trigger (regard as the start time of beta acquisition --> t = 0) from MRTOF of that MRTOF run. 
//     so the TrigTime in [ns] is required !!!!!!!!!!!!! 
//		Read_Beta_lst_batch() determines the time of first trigger from the first beta.lst automatically. If this is no information of first trigger. 3.8 [s] is used by defalut
// 5. To read the last file for current batch. Set IsLastFile ==> true. Sweep and relative time are generated for data read already before return to the master process. If synchronization signal
//		is lost. The sweep and relative time of data beyond the last Ch4 event are corrected automatically
// 6. If in some test situation, No ch4 input. The sweeps and relative time information will be reconstructed according to " first trigger time = 3.8 s and each sweep = 25 ms"
//**********************************************************************************************

bool Read_Beta_lst(string PATH,string filename, bool loadmore=false,int ADC1_CH_low_thres=100, int ADC2_CH_low_thres=100, Long64_t TrigTime=0, bool IsLastFile=true){// full path to folder with .txt
	//filename should include ".txt"
	//set "loadmore" to false at the first beta list for loading new bunch of beta list (clear the memory of vector recoding all history events), the current filename is save to "Current_beta_file"
	//set "loadmore" to true for loading beta list at the same bunch (from the 2nd sub-list, keep the memory for the beta event have loaded)
	// ADC thres. in the unit of channel!!!!!
	// set first trigger time align beta file to the start (1st sweep)of MRTOF
	// default --> no need to set trigger time. Get global time of mrtof from beta file
	//IsLastFile when true --> generate sweep numbers and time_relative for beta event. For multi-file case. Set this to true when it is the last file to load.


	//*********************  using for 4ch reconstruct if this is no signal from 4ch *****************************
	static int Index_1st_file_current_batch = 1;  

	if(!loadmore){  // get the index of the first file in the current batch to read,  
			string filename_s = filename;
			string::size_type charpos=0;
			while(filename_s.find("_",charpos)!=string::npos){charpos = filename_s.find("_",charpos)+1;}
			filename_s.erase(filename_s.begin(),filename_s.begin()+charpos);
			charpos = filename_s.find(".lst");
			if(charpos != string::npos) filename_s.erase(filename_s.begin()+charpos,filename_s.end());
			Index_1st_file_current_batch = atoi(filename_s.c_str());
	}
	//********************* end of extracting the index of the first beta in current batch ***********************


	if(BetaFastMode && TrigTime==0) cout<<"\e[1;33m"<<"Warning!!! It is BetaFastMode.  Please set time of first trig!!!!"<<endl;

	FILE* fin =NULL;
	//string fpath_name = PATH + "LST/" + filename +".txt";
	string fpath_name = PATH + filename;
	cout<<fpath_name<<endl;
	fin = fopen(fpath_name.data(),"r");
	if(fin==NULL){
		cout<<"fail to open list file of beta!! break!!"<<endl;
		return false;
	}

	// check binary or txt file
char tryget[100] = {'\0'};
fgets(tryget,100,fin);
string trygetstr = tryget;
string::size_type first_comma = trygetstr.find(",");
string::size_type second_comma = trygetstr.find(",",first_comma+1);
bool IsTXTfile=true;

if(second_comma-first_comma == 2) {fclose(fin); fin = fopen(fpath_name.data(),"r"); IsTXTfile=true;} // txt data file
else {fclose(fin); fin = fopen(fpath_name.data(),"rb"); IsTXTfile=false;}

	static Long64_t nevt_bta_old=0;
	static Long64_t nevt_bta_old_old=0;  // for sort; to have a 10% overlap with the last file
	if(!loadmore){
		 nevt_bta=0;
		 nevt_bta_old=0;
		 nevt_bta_old_old=0;
	}

	static Long64_t sweeps_gclock=0;
	if(!loadmore){
		if(BetaFastMode){
			string filetem = filename;
			while(filetem.find("_") != string::npos)filetem.erase(filetem.begin(),filetem.begin()+filetem.find("_")+1);
			if(filetem.find(".lst") != string::npos) filetem.erase(filetem.begin()+filetem.find(".lst"),filetem.end());
			if(atoi(filetem.c_str())!=1){
				 sweeps_gclock=1; //no the first file ; set one to allow the Compare_Beta_Beta to handle very beginning data in the file
				 cout<<"set initial sweeps_gclock to 1"<<endl;
			}
			else sweeps_gclock=0; // for BetaFastMode but is the first beta file in the series.
		}
		else sweeps_gclock=0; // for normal mode
	} 

	Long64_t time =0;
	static Long64_t current_sweeps_time_gclock = 0;
	if(!loadmore) current_sweeps_time_gclock=0;
	static Long64_t sweeps_marker=0;
	if(!loadmore) sweeps_marker=0;
	static bool NosynSignal=false;
	if(!loadmore) NosynSignal=false;
	int adc=0;
	int channel=0;

	int num_nothing1=0, num_nothing2=0;

	//%%%%%%%%% for make sweeps global time tree %%%%%%%%%%%%
	static Long64_t Last_sweeps_gclock =0;
	if(!loadmore){
		Last_sweeps_gclock =0;
		if(!BetaSweepsTimeReady && tbeta_sweeps_gtime!=NULL && !BetaFastMode){
				tbeta_sweeps_gtime->Branch("sweeps_gtime",&current_sweeps_time_gclock,"sweeps_gtime/L");
				tbeta_sweeps_gtime->Fill();
		}
	}	
	//%%%%%%%% end of for making sweeps global time tree %%%%%%%%%%



/*	//vector <Long64_t> *nevt_v = new vector <Long64_t>();
	vector <Long64_t> *time_v = new vector <Long64_t>();
	vector<int> *adc_v = new vector<int>();
	vector<int> *channel_v = new vector<int>();
*/

	if(!loadmore){
		sweeps_gclock_v->clear(); sweeps_gclock_v->shrink_to_fit();
		time_v->clear(); time_v->shrink_to_fit();
		time_relative_v->clear(); time_relative_v->shrink_to_fit();
		adc_v->clear();  adc_v->shrink_to_fit();
		channel_v->clear();  channel_v->shrink_to_fit();
	}

	bool firstloop=false; // turn to true after finish first read
	Long64_t StartTimeCurrentFile=0;

	while(!feof(fin)){
		if(firstloop){
//if(channel==1 && adc<1000) goto skip;
//if(channel==2 && adc<500) goto skip;
/*
		if(!BetaFastMode){
					if(channel==syn_CH){
						 sweeps_gclock++;  // in mrtof.lst file: sweeps start from 1; NOT 0
						 current_sweeps_time_gclock = time;
					}


			if(time-StartTimeCurrentFile>5000000000 && current_sweeps_time_gclock == 0) {NosynSignal=true;}// no sweeps info in file (check the first 5 second); enable approximation

			if(NosynSignal){
				if(current_sweeps_time_gclock==0) cout<<"No syn signal from 4ch!!!!"<<endl;// just print once
				sweeps_gclock = (time / SweepInterval);
				current_sweeps_time_gclock = sweeps_gclock*SweepInterval; 
				sweeps_gclock= sweeps_gclock+1;// sweeps number starts from 1

			}
		}
*/

//double relative_time = time - current_sweeps_time_gclock;

//bool in_t_window1 = relative_time >=tof_b2l[0] && relative_time <=tof_b2l[0]+ 0.12e6;
//bool in_t_window2 = relative_time >=tof_b2l[1] && relative_time <=tof_b2l[1]+ 0.12e6;
//bool in_t_window3 = relative_time >=tof_b2l[2] && relative_time <=tof_b2l[2]+ 0.12e6;
//bool in_t_window4 = relative_time >=tof_b2l[3] && relative_time <=tof_b2l[3]+ 0.12e6;

if(channel==2 || channel==1){
	//if(in_t_window1 || in_t_window2 || in_t_window3 || in_t_window4) goto skip;
}

if(adc < ADC1_CH_low_thres && channel==1) goto skip;
if(adc < ADC2_CH_low_thres && channel==2) goto skip;

		//	sweeps_gclock_v->push_back(sweeps_gclock);
		//	time_relative_v->push_back(time - current_sweeps_time_gclock);
			nevt_bta++;
			time_v->push_back(CorrectedBetaTime(time,adc,channel)); //time 
			channel_v->push_back(channel);
			adc_v->push_back(adc);
skip:
			if(nevt_bta%1000==0)cout<<'\r'<<"nevt_bta= "<<nevt_bta<<flush;
			//if(nevt_bta==10)break;
//printf("%lld,%d,%d\n",time,channel,adc);
		}


		if(IsTXTfile)	fscanf(fin,"%lld,%d,%d",&time,&channel,&adc); // using this format in general
//fscanf(fin,"%d,%d,%d,%lld,%d",&num_nothing1,&channel,&adc,&time,&num_nothing2);// just for Niwase's data
//channel++;
//adc/=16.;
		else ReadBinaryData(fin,time,channel,adc);

		if(BetaFastMode){
				time=time-TrigTime;
			/*	if(channel==syn_CH){
						 sweeps_gclock++;  // in mrtof.lst file: sweeps start from 1; NOT 0
						 current_sweeps_time_gclock = time;
				}*/

				//sweeps_gclock = (time / SweepInterval);
				//if(sweeps_gclock-sweeps_marker==113)
				//current_sweeps_time_gclock = sweeps_gclock*SweepInterval+(160*sweeps_gclock/113)%SweepInterval; 
				//sweeps_gclock= sweeps_gclock+1-(160*sweeps_gclock/113)/SweepInterval;// sweeps number starts from 1
		}

	
		if(!firstloop)StartTimeCurrentFile = time;

		firstloop=true;
//if(nevt_bta==10) return 1;
	}
	cout<<'\r'<<"nevt_bta= "<<nevt_bta<<endl;
	fclose(fin);

/*
	for(Long64_t index=0;index<100;index++){
		sweeps_gclock = sweeps_gclock_v->at(index);
		time=time_v->at(index);
		adc=adc_v->at(index);
		channel = channel_v->at(index);
		printf("%lld,%lld,%d,%d\n",sweeps_gclock,time,channel,adc);
	}
*/

	/*if(histo_h1E!=NULL) delete histo_h1E;
		histo_h1E = new TH1D("histo_h1E"," Si1 and Si2 ADC ",512,0,1024); histo_h1E->SetLineColor(kBlack);
	if(histo_h2E!=NULL) delete histo_h2E;
		histo_h2E = new TH1D("histo_h2E"," Si1 and Si2 ADC ",512,0,1024); histo_h2E->SetLineColor(kRed);
	for(Long64_t index=0;index<nevt_bta;index++){
		if(channel_v->at(index)==1) histo_h1E->Fill(adc_v->at(index));
		if(channel_v->at(index)==2) histo_h2E->Fill(adc_v->at(index));
	}*/


	//%%%%%%%%%%%%%%%%%%%%% sort!!!! mess sequence found under high throughout rate
	Long64_t eventoffset=0;
	if(nevt_bta_old>0) eventoffset = nevt_bta_old_old + (Long64_t)((nevt_bta_old-nevt_bta_old_old)*0.9);

	Long64_t nevt_bta_backup = nevt_bta;

	if(MakeBetaRawTree) nevt_bta = (Long64_t)time_v->size();

	Long64_t* timecopy = new Long64_t[nevt_bta-eventoffset];
	int* adccopy = new int[nevt_bta-eventoffset];
	int* channelcopy = new int[nevt_bta-eventoffset];
	Long64_t* sequencelist = new Long64_t[nevt_bta-eventoffset];

	for(Long64_t index=eventoffset;index<nevt_bta;index++){
		timecopy[index-eventoffset] = time_v->at(index);
		adccopy[index-eventoffset] = adc_v->at(index);
		channelcopy[index-eventoffset] = channel_v->at(index);
	}

	TMath::Sort(nevt_bta-eventoffset,timecopy,sequencelist,kFALSE);

	for(Long64_t index=eventoffset;index<nevt_bta;index++){
		time_v->at(index) = timecopy[sequencelist[index-eventoffset]];
		adc_v->at(index) = adccopy[sequencelist[index-eventoffset]];
		channel_v->at(index) = channelcopy[sequencelist[index-eventoffset]];
/*
		if(channel_v->at(index)==syn_CH){
			sweeps_gclock++;
			current_sweeps_time_gclock = time_v->at(index);
		}

		sweeps_gclock_v->push_back(sweeps_gclock);
		time_relative_v->push_back(time_v->at(index)- current_sweeps_time_gclock);*/

	}

	nevt_bta = nevt_bta_backup;


	//**************** if make raw beta tree **************
	// load data -> sort -> calculate sweep -> fill to tree ->shrink vector
if(!MakeBetaRawTree){
	nevt_bta_old_old=nevt_bta_old;
	nevt_bta_old = nevt_bta; //update
}


	if(IsLastFile || MakeBetaRawTree){//%%%%%%%%%%%%%%%% for generating sweep number and time relative of beta event after geting a correct sequence
		//sweeps_gclock_v->clear();sweeps_gclock_v->shrink_to_fit();
		//time_relative_v->clear();time_relative_v->shrink_to_fit();
		//current_sweeps_time_gclock=0;

		bool Ch4_Reconstruct_required = true;

		for(unsigned long allindex=0;allindex<time_v->size();allindex++){
			if(channel_v->at(allindex)==syn_CH){
					 Ch4_Reconstruct_required = false;
					 break;
			}
			if(time_v->at(allindex) - time_v->at(0) > 4100000000) break; // just read the first 4.1 s maybe enough
		}


		if(Ch4_Reconstruct_required){ // Have to reconstruct ch4 signal
			if(BetaFastMode){
				for(unsigned long allindex=0;allindex<time_v->size();allindex++){
					if(time_v->at(allindex)<0){ // it is the first beta file in the series, using sweeps_glock =0, and current_sweeps_time_gclock =0
								sweeps_gclock_v->push_back(0);
								time_relative_v->push_back(time_v->at(allindex));
					}else{
								sweeps_gclock = time_v->at(allindex) / (SweepInterval/2);
								current_sweeps_time_gclock = sweeps_gclock * (SweepInterval/2);
								sweeps_gclock_v->push_back(sweeps_gclock+1);
								time_relative_v->push_back(time_v->at(allindex)- current_sweeps_time_gclock);
					}
				}

			}else{// normal mode

				unsigned long EndIndex=0;
				if(IsLastFile) EndIndex = time_v->size(); // except for "MakeBetaRawTree" mode, generate "sweeps_gclock" for all events in RAM; if this is the last file, generate "sweeps_gclock" for all events in RAM anyway.
				else EndIndex = nevt_bta_old;

				for(unsigned long allindex=0;allindex<EndIndex;allindex++){
					if(time_v->at(allindex)<3800000000){ // just assume the first trigger will come at 3800000000 ns like MCS at acquisition
							sweeps_gclock_v->push_back(0);
							time_relative_v->push_back(time_v->at(allindex));
					}else{
							sweeps_gclock = (time_v->at(allindex)-3800000000) / (SweepInterval/2);
							current_sweeps_time_gclock = sweeps_gclock * (SweepInterval/2) + 3800000000;
							sweeps_gclock_v->push_back(sweeps_gclock+1);
							time_relative_v->push_back(time_v->at(allindex)- current_sweeps_time_gclock);

							if(sweeps_gclock+1 > Last_sweeps_gclock && !BetaSweepsTimeReady && tbeta_sweeps_gtime!=NULL && !BetaFastMode){ tbeta_sweeps_gtime->Fill();  Last_sweeps_gclock = sweeps_gclock+1;} // fill sweeps global time tree

					}
				}
			}

		}else{ // No need to reconstruct ch4 

				if(!MakeBetaRawTree){
					if(time_v->at(0)-3800000000<0){
						sweeps_gclock = 0;
						current_sweeps_time_gclock=0;
					}else{
						sweeps_gclock = (time_v->at(0)-3800000000) / (SweepInterval/2);
						current_sweeps_time_gclock = sweeps_gclock * (SweepInterval/2) + 3800000000;
						sweeps_gclock++;  // sweeps number start from 1
					}
				}


				unsigned long EndIndex=0;

				if(IsLastFile) EndIndex = time_v->size(); // except for "MakeBetaRawTree" mode, generate "sweeps_gclock" for all events in RAM; if this is the last file, generate "sweeps_gclock" for all events in RAM anyway.
				else EndIndex = nevt_bta_old;


				for(unsigned long allindex=0;allindex<EndIndex;allindex++){
						if( (channel_v->at(allindex)==syn_CH && allindex>0) || (channel_v->at(allindex)==syn_CH  && (MakeBetaRawTree || !BetaFastMode) ) ){ 
																						// if not make beta raw tree,  index=0 case has been used for determining the global sweep number above. Add the condiction allindex > 0, avoid the 
																						// repetition of increase sweeps_gclock by one again when the index=0 is from channel 4
																						// if make beta raw tree, it must be under "NORMAL MODE", it will not use the above time_v->at(0) to determine "sweeps_glock". Do not skip the 
																						// index=0 judgement.
							sweeps_gclock++;
							current_sweeps_time_gclock = time_v->at(allindex);
						}

						sweeps_gclock_v->push_back(sweeps_gclock);
						time_relative_v->push_back(time_v->at(allindex)- current_sweeps_time_gclock);

						if(sweeps_gclock > Last_sweeps_gclock && !BetaSweepsTimeReady && tbeta_sweeps_gtime!=NULL && !BetaFastMode){ tbeta_sweeps_gtime->Fill();  Last_sweeps_gclock = sweeps_gclock;}


				}
		}



					//%%%%%%%%%%%%% save sweeps_gtime tree %%%%%%%%%%%%%%%%%%%5

					if(IsLastFile && !BetaSweepsTimeReady && tbeta_sweeps_gtime!=NULL && !BetaFastMode){
						f_sweeps_gtime_out->cd();
						tbeta_sweeps_gtime->Write();
						f_sweeps_gtime_out->Write();
						BetaSweepsTimeReady= true;
					}

					//%%%%% end of saving sweeps_gtime tree %%%%%%%%%%%%%%%%%%%



		//%%%%%%%%%%%%%% reconstruct the sweep trigger for  data coming after the stop of MCS  %%%%%%%%%%%%%%%%%%%%%%%%5

		if(IsLastFile){
				unsigned long LasttrigIndex=0;
				Long64_t LasttrigTime=0;
				Long64_t triInterval1=0;
				Long64_t triInterval2=0;
				Long64_t sweeps_lasttrig=0;

				for(unsigned long allindex=time_v->size()-1;;allindex--){
					if(channel_v->at(allindex) == syn_CH){
						if(LasttrigIndex==0){
							LasttrigIndex = allindex;
							LasttrigTime = time_v->at(allindex);
		sweeps_lasttrig = sweeps_gclock_v->at(allindex);
							continue;
						}
						else{
							Long64_t temInterval = LasttrigTime-time_v->at(allindex);
							if(temInterval%50000040>2000){
								triInterval1 =temInterval%50000040;  // actuall interval of every two sweeps. USB MCA resolution is 40 ns
								break;
							}
							
						}

					}

					if(allindex==0) break;

				}

				if(triInterval1>2000){
					triInterval2=50000040-triInterval1;
				}


				if(triInterval2>0){
		current_sweeps_time_gclock=LasttrigTime;
		sweeps_gclock=sweeps_lasttrig;
		int extra_sweep=0;
						for(unsigned long allindex=LasttrigIndex+1;allindex<time_v->size();allindex++){
								if(extra_sweep%2==0){
										if(time_v->at(allindex)-current_sweeps_time_gclock >triInterval2){
											current_sweeps_time_gclock+=triInterval2;
											sweeps_gclock++;
											extra_sweep++;
										}

								}
								else if(extra_sweep%2==1){
										if(time_v->at(allindex)-current_sweeps_time_gclock >triInterval1){
											current_sweeps_time_gclock+=triInterval1;
											sweeps_gclock++;
											extra_sweep++;
										}
								}

								sweeps_gclock_v->at(allindex)=sweeps_gclock;
								time_relative_v->at(allindex)=time_v->at(allindex)- current_sweeps_time_gclock;

						}


				} 
		}// this is the end for reconstruct the sweep for data coming after the stop of MCS

	}// the end of if(IsLastFile) to get or calculate the sweep informaiton for beta



	if(MakeBetaRawTree){ //%%%%%%%%%%%%%%%%%%%%%%%%%% fill beta raw data from vector to tree

		Long64_t time_tem =0;
		int channel_tem=0;
		int adc_tem=0;
		Long64_t sweeps_gclock_tem=0;
		Long64_t time_relative_tem=0;

		if(tbeta_Raw != NULL){
			tbeta_Raw->SetBranchAddress("time",&time_tem);
			tbeta_Raw->SetBranchAddress("channel",&channel_tem);
			tbeta_Raw->SetBranchAddress("adc",&adc_tem);
			tbeta_Raw->SetBranchAddress("sweeps_global",&sweeps_gclock_tem);
			tbeta_Raw->SetBranchAddress("time_relative",&time_relative_tem);
		}

		vector<Long64_t> time_tem_v;
		vector<int> adc_tem_v;
		vector<int> channel_tem_v;

		unsigned long EndIndex=0;
		if(IsLastFile) EndIndex = time_v->size(); // except for "MakeBetaRawTree" mode, generate "sweeps_gclock" for all events in RAM; if this is the last file, generate "sweeps_gclock" for all events in RAM anyway.
		else EndIndex = nevt_bta_old;


		if(EndIndex>0){
			for(unsigned long allindex=0;allindex<EndIndex;allindex++){ // fill upper section of vector
				time_tem = time_v->at(allindex);
				channel_tem = channel_v->at(allindex);
				adc_tem = adc_v->at(allindex);
				sweeps_gclock_tem = sweeps_gclock_v->at(allindex);
				time_relative_tem = time_relative_v->at(allindex);
				if(tbeta_Raw != NULL) tbeta_Raw->Fill();
			}

			for(unsigned long allindex=EndIndex;allindex<time_v->size();allindex++){ // backup lower section of vector
				time_tem_v.push_back(time_v->at(allindex));
				adc_tem_v.push_back(adc_v->at(allindex));
				channel_tem_v.push_back(channel_v->at(allindex));
			}

		}


		nevt_bta_old = time_v->size() - nevt_bta_old;

		if(time_tem_v.size()>0){
			*time_v = time_tem_v;
			*adc_v = adc_tem_v;
			*channel_v = channel_tem_v;
		}


		sweeps_gclock_v->clear();
		time_relative_v->clear();

		if(IsLastFile && tbeta_Raw != NULL){
			f_betaRaw_out->cd();
			tbeta_Raw->Write();
			f_betaRaw_out->Write();

		}


		if(tbeta_Raw != NULL) tbeta_Raw->ResetBranchAddresses();

	}//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of fill beta raw data from vector to tree



	if(!loadmore){
		if(filename.find("LST/")!=string::npos) filename.erase(filename.begin(),filename.begin()+filename.find("LST/")+4);
		Current_beta_file = filename;
	}

//temporary disconnect the auto update for a online path
//	Path_beta_file = PATH;

	if(!ReadBetaOnline) Path_beta_file = PATH;  // for BetaOnlineMode, Path_beta_file is set in 

	delete[] timecopy;
	delete []adccopy;
	delete[] channelcopy;
	delete[] sequencelist;


	return true;
}


bool Read_Beta_lst_batch(string PATH , string filename, int ADC1_CH_low_thres=100, int ADC2_CH_low_thres=100,Long64_t LoadTimeRangeL=-1, Long64_t LoadTimeRangeR=-1, bool ReloadBeta=true){
	// ADC thres. in the unit of CHANNEL; LoadTimeRange in the unit of ns
	// please set ReloadBeta= true; when changing the batch of file. For example, changing to another run;
	// ReloadBeta is for loading the list of beta file name and reloading info of STARTING TIME of each file. when false, not to extract the starting time of each file again

	if(MakeBetaRawTree){ // make tree is available only under NORMAL mode "!BetaFastMode"
		if(f_betaRaw_out != NULL){
			f_betaRaw_out->Close();
			tbeta_Raw=NULL;
			f_betaRaw_out=NULL;
		}

		if(BetaFastMode){
			cout<<endl;
			cout<<endl;
			cout<<"\e[1;33m MakeBetaRawTree option is invalid under \"BetaFastMode\". By default, BetaFastMode is set to false to continue.\e[0m"<<endl;
			BetaFastMode = false;
			sleep(2);
		}
	}

	BetaSweepsTimeReady=false;
	if(f_sweeps_gtime_out !=NULL){
		f_sweeps_gtime_out->Close();
		tbeta_sweeps_gtime = NULL;
		f_sweeps_gtime_out = NULL;
	}


	NeedToCompare=true;  // in BetaFastMode, if there is no new beta.lst file read, not need to compare again. default is "to compare beta-beta"

	if(!BetaFastMode) ReloadBeta=true; //

	if(BetaFastMode){
		if(LoadTimeRangeL>LoadTimeRangeR) swap(LoadTimeRangeL,LoadTimeRangeR);
	}

	if(BetaFastMode && (LoadTimeRangeL<0 || LoadTimeRangeR<0)){
		cout<<"\e[1;31m"<<"Warning!!! To use BetaFastMode, time range required!!! Abort!!!!"<<"\e[0m"<<endl;
		return false;
	}

	if(BetaFastMode && ReloadBeta == false){ // skip process of reading filename list
		if(Current_First_beta_file_FastMode != ToLoad_beta_file ){  // beta file is different from history --> should ReloadBeta
            if(ToLoad_beta_file=="+++"){
              cout<<"warning: no beta file to be load; please load file first"<<endl;
              cout<<"\e[1;33m"<<"beta filename with .lst"<<"\e[0m"<<endl;
              cout<<"using:  LoadNewBetaFile("<<endl;
              return false;
            }
            else{
            	cout<<"Maybe new beta file to load. ReloadBeta = true is necessary. try again"<<endl;
            	return false;
            }
        }else{ // same beta file --> Can skip get filename list process 
        		cout<<"Same batch of beta file, can skip get filename list process"<<endl;
        }

	}
 
 	static string PATH_tem=""; // if the ReloadBeta option is not chosen, will use the PATH_tem modified last time. Otherwise in BetaFastMode, correction to the PATH will be skipped expect for the first beta.lst encounter after press 'c'
 										// because the starting time of each beta.lst will not be read out again from the second TOF in the specified region by cursor.

	static int fstartindex_old=0;
	static int fendindex_old=0;
	static vector <string> filenamelist;
	static vector <int> fileIndexSequence;
	if(ReloadBeta){
		fstartindex_old =-1;
		fendindex_old =-1;
		filenamelist.clear();
		fileIndexSequence.clear();
	} 

	int * filenameIndex_in = NULL;
	int * filenameIndex_out = NULL;

	string beta_fine_name_head = filename;

	if(ReloadBeta){// generate a list of beta files and a list of correct index
				while(beta_fine_name_head.find("/")!= string::npos) beta_fine_name_head.erase(beta_fine_name_head.begin(),beta_fine_name_head.begin()+
																					beta_fine_name_head.find("/")+1);
				string::size_type charpos=0;
				while(beta_fine_name_head.find("_",charpos)!=string::npos){charpos = beta_fine_name_head.find("_",charpos)+1;} // from "L_1999_777_001.lst" to "L_1999_777_"
				beta_fine_name_head.erase(beta_fine_name_head.begin()+charpos,beta_fine_name_head.end());
			cout<<beta_fine_name_head.c_str()<<endl;

				if(PATH.find("/rootfiles/../")!=string::npos){
					MainFilePath = PATH; // BetaCut.root is save to mother path of "LST"  and "rootfiles"
					PATH_tem = PATH+"LST/";		// tbeta_beta tree is save to "rootfiles" by setting Path_beta_file = PATH in Read_Beta_lst(). Note: tbeta_beta tree is always save to Path_beta_file+"../rootfiles/"
				}else if(PATH.find("/rootfiles/")!=string::npos){
					MainFilePath = PATH + "../";  // BetaCut.root is save to mother path of "LST"  and "rootfiles"
					PATH_tem = PATH+"../LST/";  // tbeta_beta tree is save to "rootfiles"
				}else{// 
					MainFilePath = PATH; // scenario 1: using Beta_Beta_coin.cc independently and locally; Betacut.root is save to the same Path as that of beta.lst files; tbeta_beta tree is save to Path_beta_file (=PATH  in Read_Beta_lst())
					PATH_tem = PATH;
					if(ReadBetaOnline){ MainFilePath = Path_beta_file+"../";} // scenario 2: under BetaOnline mode, Path_beta_file is set to FilePath+"../LST/" (FilePath = "/xxx/xxx/rootfiles/") in Exec.h; Therefore 
																							//	BetaCut.root is save to mother path of "LST"; 	tbeta_beta tree is save to "rootfiles"		

																					//SCENARIO 3: if want to using independently but online; set Path_beta_file manually. tbeta_beta tree is save to Path_beta_file; BetaCut.root is save to Path_beta_file+"../" ;
				}

				string linux_command = "ls -tr " + PATH_tem + beta_fine_name_head + "*.lst";

				char getfilename[1000]={'\0'};

				FILE* pipe = popen(linux_command.c_str(), "r");
				if (!pipe) {
					std::cerr << "popen() failed!" << std::endl;
					return false;
				}

			try{
				while(!feof(pipe))
				{fgets(getfilename,1000,pipe);
					if(!feof(pipe))filenamelist.push_back(getfilename);
					else break;
				}
			}
			catch(...){
				pclose(pipe);
				throw;
			}


			pclose(pipe);

			filenameIndex_in= new int[filenamelist.size()];
			filenameIndex_out=new int[filenamelist.size()];

	
		for(int ifile=0;ifile<(int)filenamelist.size();ifile++){
			while(filenamelist[ifile].find("/")!= string::npos){
				filenamelist[ifile].erase(filenamelist[ifile].begin(),filenamelist[ifile].begin()+filenamelist[ifile].find("/")+1); // from "/home/xxx/xxx/L_123_001.lst" to "L_123_001.lst"
			}
			filenamelist[ifile].erase(filenamelist[ifile].begin()+filenamelist[ifile].find("lst")+3, filenamelist[ifile].end()); // clear any other characters after lst in case there are to have pure "L_123_001.lst"
			//cout<<filenamelist[ifile]<<endl;

			string filename_i_tem = filenamelist[ifile];
			filename_i_tem.erase(filename_i_tem.begin()+filename_i_tem.find(".lst"),filename_i_tem.end()); // from "L_123_001.lst" to "L_123_001"

			while(filename_i_tem.find("_")!=string::npos){
				filename_i_tem.erase(filename_i_tem.begin(),filename_i_tem.begin()+filename_i_tem.find("_")+1); // from "L_123_001" to "001"
			}

			filenameIndex_in[ifile]= atoi(filename_i_tem.c_str()); // from "001" to interger 1
			//PATH = PATH+"LST/";
		}

		TMath::Sort((int)filenamelist.size(),filenameIndex_in,filenameIndex_out,kFALSE);

	}// skip when no need to reload beta



	if(!BetaFastMode){// load all related beta file

		for(int ifile=0;ifile<(int)filenamelist.size();ifile++){
			cout<<filenamelist[filenameIndex_out[ifile]]<<endl;
			//PATH = PATH+"LST/"; // path has included LST/
//continue;
			if(ifile==0){
				//Current_First_beta_file = filenamelist[filenameIndex_out[ifile]]; // seems useless now
				bool IsLastFile=false;
				if(filenamelist.size()==1) IsLastFile=true;


				string ftem = filenamelist[filenameIndex_out[ifile]];
				if(ftem.find(".lst") != string::npos) ftem.erase(ftem.begin()+ftem.find(".lst"), ftem.end());
				string Path_BetaRaw = MainFilePath + "rootfiles/";

				//%%%%%%%%%%%%  check whether to make a tree for sweeps_gtime %%%%%%%%%%%%%%%%%%%%%%
				string f_sweeps_gtime_name = ftem + "_sweep_gtime.root";

				if(  gSystem->AccessPathName( (Path_BetaRaw+f_sweeps_gtime_name).c_str() )  ){// tree now exist
					cout<<"\e[1;33m Find no record of sweeps_gtime tree. Create a new one \e[0m"<<endl;
						BetaSweepsTimeReady=false;
						f_sweeps_gtime_out = new TFile( (Path_BetaRaw+f_sweeps_gtime_name).c_str(),"RECREATE");
						if(f_sweeps_gtime_out ==NULL){
							cout<<"\e[1;31m Error!! Fail to create root file to save sweeps_gtime tree!!!! Abort!!!\e[0m"<<endl;
							if(filenameIndex_in != NULL) delete[] filenameIndex_in;
							if(filenameIndex_out != NULL) delete[] filenameIndex_out;
							return false;

						}

						cout<<endl;
						cout<<"\e[1;33m Create "<<(Path_BetaRaw+f_sweeps_gtime_name).c_str()<<"   to save sweeps_gtime tree\e[0m"<<endl;
						f_sweeps_gtime_out->cd();
						tbeta_sweeps_gtime = new TTree("sweeps_gtime" ,"sweeps_gtime"); // set branch in Read_Beta_lst();

				}else{// sweeps_gtime tree exist already

						f_sweeps_gtime_out = new TFile( (Path_BetaRaw+f_sweeps_gtime_name).c_str(),"READ");
						if(f_sweeps_gtime_out ==NULL){
							cout<<"\e[1;31m Error!! Fail to read root file of sweeps_gtime tree!!!! Abort!!!\e[0m"<<endl;
							if(filenameIndex_in != NULL) delete[] filenameIndex_in;
							if(filenameIndex_out != NULL) delete[] filenameIndex_out;
							return false;

						}

						tbeta_sweeps_gtime = (TTree*)f_sweeps_gtime_out->Get("sweeps_gtime");
						BetaSweepsTimeReady=true;
						cout<<endl;
						cout<<"Get sweeps_gtime tree from history."<<endl<<endl;
				}
				//%%%%%%%%%%%%%%%%% finish checking whether to make a tree for sweeps_gtime %%%%%%%%%%%%%%%%



				//%%%%%%%%%%%%%%%% open or create beta raw tree %%%%%%%%%%%%%%%%
				if(MakeBetaRawTree){
					time_v->clear(); 				time_v->shrink_to_fit();
					channel_v->clear();			channel_v->shrink_to_fit();
					adc_v->clear();					adc_v->shrink_to_fit();
					sweeps_gclock_v->clear();	sweeps_gclock_v->shrink_to_fit();
					time_relative_v->clear();	time_relative_v->shrink_to_fit();


					string File_betaRaw_tree = Form("%s_Thres1_%d_Thres2_%d_RawBeta.root",ftem.c_str(),ADC1_CH_low_thres,ADC2_CH_low_thres);

					Path_BetaRaw = Path_BetaRaw + File_betaRaw_tree;

					f_betaRaw_out = new TFile(Path_BetaRaw.c_str(),"READ");
					if(f_betaRaw_out!=NULL) tbeta_Raw = (TTree*)f_betaRaw_out->Get("tbeta_Raw");

					if(f_betaRaw_out==NULL || tbeta_Raw==NULL){ // beta raw tree not exist
						f_betaRaw_out = new TFile(Path_BetaRaw.c_str(),"RECREATE");
						if(f_betaRaw_out ==NULL){
							cout<<"\e[1;31m Error!! Fail to create root file to save betaRaw tree!!!! Abort!!!\e[0m"<<endl;
							if(filenameIndex_in != NULL) delete[] filenameIndex_in;
							if(filenameIndex_out != NULL) delete[] filenameIndex_out;
							return false;

						}


						cout<<endl;
						cout<<"\e[1;33m Create "<<Path_BetaRaw<<"   to save betaRaw tree\e[0m"<<endl;
						f_betaRaw_out->cd();
						tbeta_Raw = new TTree("tbeta_Raw","tbeta_Raw");

						Long64_t time_tem =0;
						int channel_tem=0;
						int adc_tem=0;
						Long64_t sweeps_gclock_tem=0;
						Long64_t time_relative_tem=0;
		
						tbeta_Raw->Branch("time",&time_tem,"time/L");
						tbeta_Raw->Branch("channel",&channel_tem,"channel/I");
						tbeta_Raw->Branch("adc",&adc_tem,"adc/I");
						tbeta_Raw->Branch("sweeps_global",&sweeps_gclock_tem,"sweeps_global/L");
						tbeta_Raw->Branch("time_relative",&time_relative_tem,"time_relative/L");

					}else{// read history betaRaw tree successfully
						cout<<endl;
						cout<<"\e[1;33m betaRaw tree \""<<File_betaRaw_tree<<"\" exist. Load Beta Raw data from history\e[0m"<<endl;
						if(filenameIndex_in != NULL) delete[] filenameIndex_in;
						if(filenameIndex_out != NULL) delete[] filenameIndex_out;
						return true;
		
					}
				}
				//%%%%%%%%%%%%%%%% end of open or create beta raw tree %%%%%%%%%%%%%%%%%%%%

			 	if(!Read_Beta_lst(PATH_tem,filenamelist[filenameIndex_out[ifile]],false, ADC1_CH_low_thres, ADC2_CH_low_thres,0,IsLastFile)){// general mode does not need trigger time; align mrtof to bete file
			 		delete[] filenameIndex_in;
					delete[] filenameIndex_out;
			 		return false; 
			 	} 
			}
			else{
				bool IsLastFile=false;
				if(ifile==(int)filenamelist.size()-1) IsLastFile=true;
				 if(!Read_Beta_lst(PATH_tem,filenamelist[filenameIndex_out[ifile]],true, ADC1_CH_low_thres, ADC2_CH_low_thres,0,IsLastFile)){	
				 	delete[] filenameIndex_in;
					delete[] filenameIndex_out;
					return false;
				 } 
			}
		}
	}else{// load only specific beta file which cover the range of time window

		Long64_t gettime_tem;
		int getchannel_tem;
		int getadc_tem;
	//	static vector<Long64_t> FileStartTime;// store the start time of each file
		static Long64_t FirstTrigTime=3800000000;   //the first trigger come about 3.8 s after the start of usb MCA

		if(ReloadBeta){ //%%%%%%%%%%%%%%%%% ReadloadBeta in BetaFastMode %%%%%%%%%%%%%%%%%

			FileStartTime.clear(); FileStartTime.shrink_to_fit();
			FirstTrigTime=3800000000;

			for(int ifile=0;ifile<(int)filenamelist.size();ifile++){// extract the time of first event in each beta file
		if(ifile==0)Current_First_beta_file_FastMode = filenamelist[filenameIndex_out[0]];
	fileIndexSequence.push_back( filenameIndex_out[ifile] );
				string filePathandName = PATH_tem + filenamelist[filenameIndex_out[ifile]];	

//cout<<"Paht = "<<PATH<<endl;
//cout<<ifile<<"   "<<filenamelist[filenameIndex_out[ifile]]<<endl;

				FILE* FileHandle = fopen(filePathandName.c_str(),"r");
				if(FileHandle == NULL){
					cout<<"\e[1;31m"<<"Warning!!! Can not open file: "<<filePathandName<<"  Abort!!!!"<<"\e[0m"<<endl;
					if(filenameIndex_in != NULL) delete[] filenameIndex_in;
					if(filenameIndex_out != NULL) delete[] filenameIndex_out;
					return false;
				}

			// ////////////////Checking data format .txt or binary endian
				char tryget[100] = {'\0'};
				fgets(tryget,100,FileHandle);
				string trygetstr = tryget;
				string::size_type first_comma = trygetstr.find(",");
				string::size_type second_comma = trygetstr.find(",",first_comma+1);
				bool IsTXTfile=true;

				if(second_comma-first_comma == 2) {fclose(FileHandle); FileHandle = fopen(filePathandName.data(),"r"); IsTXTfile=true;} // txt data file
				else {fclose(FileHandle); FileHandle = fopen(filePathandName.data(),"rb"); IsTXTfile=false;}
				//////////////////////////////////////////////////////////////////////

				if(IsTXTfile)	fscanf(FileHandle,"%lld,%d,%d",&gettime_tem,&getchannel_tem,&getadc_tem);
				else	ReadBinaryData(FileHandle,gettime_tem,getchannel_tem,getadc_tem);
cout<<" file star time "<<ifile+1<<": "<<gettime_tem/1e9<<"  [s]"<<endl;
					FileStartTime.push_back(gettime_tem); //get the start time of current file
					if(ifile==0){
						for(int i=1;i<200000000 && !feof(FileHandle);i++){
							if(IsTXTfile)	fscanf(FileHandle,"%lld,%d,%d",&gettime_tem,&getchannel_tem,&getadc_tem);
							else ReadBinaryData(FileHandle,gettime_tem,getchannel_tem,getadc_tem);
							if(getchannel_tem==4){FirstTrigTime=gettime_tem;	break;	}
						}
						if(FirstTrigTime==3800000000){
							cout<<"\e[1;33m"<<"Warning!!! Find no signal from CH4!!! FirstTrigTime is set to 3.8 s as default!!!"<<"\e[0m"<<endl;
						}
cout<<"First trig in second = "<<FirstTrigTime/1e9<<"  , in ns = "<<FirstTrigTime<<endl;
					}// if find no signal from CH4, First trigger time keep as default value 3.8 s
					fclose(FileHandle);
			}

		}//%%%%%%%%%%% end of ReloadBeta process in BetaFastMode

		if(FileStartTime.size()==0){
			cout<<"\e[1;31m"<<"Warning!!! File starting time Info is not available. Please use option ReloadBeta=true for loading starting time Info of each file.  Abort!!!!"<<"\e[0m"<<endl;
			if(filenameIndex_in != NULL) delete[] filenameIndex_in;
			if(filenameIndex_out != NULL) delete[] filenameIndex_out;
			return false;
		}


		int fstartindex=0, fendindex=0;
		bool GetLowindex=false, GetHighindex=false;

		for(int index=(int)filenamelist.size()-1;index>=0;index--){ // from last file to the first file
			if( (FileStartTime[index]-FirstTrigTime+ModuleDelay) <LoadTimeRangeL && !GetLowindex) {fstartindex =index; GetLowindex = true;}
			if( (FileStartTime[index]-FirstTrigTime+ModuleDelay) <LoadTimeRangeR && !GetHighindex) {fendindex = index; GetHighindex = true;}
			if(GetLowindex && GetHighindex) break;
		}

//cout<<fstartindex<<"\t"<<fendindex<<endl;
		printf("\033[1;33m%03d.lst\033[0m, TwinL-filestart= %.3f, \033[1;33m%03d.lst\033[0m, TwinR-filestart= %.3f\n ",fstartindex+1,
			(LoadTimeRangeL-FileStartTime[fstartindex]+FirstTrigTime-ModuleDelay)/1e9,fendindex+1,(LoadTimeRangeR-FileStartTime[fendindex]+FirstTrigTime-ModuleDelay)/1e9);

		if(fstartindex < fstartindex_old || fendindex > fendindex_old || fstartindex_old==-1 || fendindex_old==-1){
			for(int index=fstartindex;index<=fendindex;index++){
				bool IsSuccess=true;
cout<<filenamelist[fileIndexSequence[index]]<<endl; //continue;
bool IsLastFile=false;
if(index==fendindex)IsLastFile=true;
				if(index==fstartindex) IsSuccess = Read_Beta_lst(PATH_tem,filenamelist[fileIndexSequence[index]],false, ADC1_CH_low_thres, ADC2_CH_low_thres,FirstTrigTime-ModuleDelay,IsLastFile);
				else IsSuccess = Read_Beta_lst(PATH_tem,filenamelist[fileIndexSequence[index]],true,ADC1_CH_low_thres, ADC2_CH_low_thres,FirstTrigTime-ModuleDelay,IsLastFile);
				// set trigger time(the last option) to align the start time of beta file to the start of mrtof
				if(!IsSuccess){
					cout<<"\e[1;33m"<<"error occur when loading beta file: "<< filenamelist[fileIndexSequence[index]] <<"basing on mrTOF!!"<<"\e[0m"<<endl;
					if(filenameIndex_in != NULL) delete[] filenameIndex_in;
					if(filenameIndex_out != NULL) delete[] filenameIndex_out;
					return false;
				}

			}

			fstartindex_old=fstartindex;
			fendindex_old=fendindex;

		}else{
			NeedToCompare=false;
			printf("Requried file is the same %03d to %03d,  skip loading.\n",fstartindex_old+1,fendindex_old+1);
		}



		if(filenameIndex_in != NULL) delete[] filenameIndex_in;
		if(filenameIndex_out != NULL) delete[] filenameIndex_out;
		return true;
		

	}// Yes or no BetaFastMode	




if(filenameIndex_in != NULL) delete[] filenameIndex_in;
if(filenameIndex_out != NULL) delete[] filenameIndex_out;


	return true;

}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  define in RefineBeta.h %%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct match_event{
	Long64_t B1_sweeps_gclock;
	Long64_t B1_time;
	Long64_t B1_time_relative;
	int B1_adc;
	Long64_t B2_sweeps_gclock;
	Long64_t B2_time;
	Long64_t B2_time_relative;
	int B2_adc;
	int delta_t;
	void clear(){
		B1_sweeps_gclock=-1;
		B1_time=-1;
		B1_time_relative=-1;
		B1_adc=-1;
		B2_sweeps_gclock=-1;
		B2_time=-1;
		B2_time_relative=-1;
		B2_adc=-1;
		delta_t=-1;
	}
};

struct hit_info{
	bool empty;   // container empty or not empty=> true; full=>false
	Long64_t sweeps_gclock;
	Long64_t time;
	Long64_t time_relative;
	int channel;
	int adc;
};
*/
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vector<match_event>Match_hit; // store matching pairs

void KeepNewMatching(match_event _NewMatch){ // keep only new matching pairs; skip same ones

	if(Match_hit.size()==0){ Match_hit.push_back(_NewMatch); return;}

	for(unsigned long index=Match_hit.size()-1;;index--){
		Long64_t time1 = Match_hit[index].B1_time;
		Long64_t time2 = Match_hit[index].B2_time;
		if(time1==_NewMatch.B1_time && time2==_NewMatch.B2_time) return; // skip; do not save repeated data
		if((_NewMatch.B1_time-time1) > 25e6 && (_NewMatch.B2_time-time2 >25e6)) break;
		if(index==0)break;
	}

	Match_hit.push_back(_NewMatch); // store new one;
	

}

//void Compare_Beta_Beta(Long64_t gate_time=500,int gate_adc_low=0, int gate_adc_hi=1000){
bool Compare_Beta_Beta(){

	static string file_old;
	static Long64_t gate_time_old=-1;   // for remember the condition of last run, skip "beta Compare" if setting no change
	static Long64_t Tgate_cent_offset_old=-1;
	static int gate_adc_low_old_1=-1, gate_adc_low_old_2=-1;         // for remember the condition of last run
	static int gate_adc_hi_old_1=-1,  gate_adc_hi_old_2=-1;  // adc gate are set in the unit of CHANNEL !! Can be set by SetCondiction(....)
	static int CutTag_old=0; // CutModifyCounter++; both Makecut and loadcut increase CutModifyCounter;
	static int PreSetADC1_low_Thres_ch_old=-1;
    static int PreSetADC2_low_Thres_ch_old=-1;
    static bool Last_BetaFastMode=false;
    static bool AlphaMode_old=false;
	bool WarningPrinted=false; // for printing warning of No beta noise cut for just once

	if(!BetaFastMode){// BetaFastMode will always update;
		if(Tgate_cent_offset_old == Tgate_cent_offset && gate_time_old==gate_time && gate_adc_low_old_1==gate_adc_low_1 && gate_adc_hi_old_1==gate_adc_hi_1 
			&& gate_adc_low_old_2==gate_adc_low_2 && gate_adc_hi_old_2==gate_adc_hi_2 && file_old==Current_beta_file && CutTag_old==CutModifyCounter
			&& PreSetADC1_low_Thres_ch_old == PreSetADC1_low_Thres_ch && PreSetADC2_low_Thres_ch_old == PreSetADC2_low_Thres_ch && AlphaMode_old == alphaMode && !Last_BetaFastMode){
			cout<<"\e[1;33m"<<"Beta-Beta coincident condition NO change! break"<<"\e[0m"<<endl;
			 return false;
		}
	}else{
		IsBetaFromNewFile(0,0,true); // update the list of files have been precessed preciously. Avoid repetitive fill the histogram in BetaFastMode
		if(!NeedToCompare){
			cout<<"\e[1;33m Same batch of beta.lst . Skip comparing process \e[0m"<<endl;
			return false;
		}
	}

	gate_time_old = gate_time ;
	Tgate_cent_offset_old = Tgate_cent_offset;
	gate_adc_low_old_1 = gate_adc_low_1;
	gate_adc_hi_old_1 = gate_adc_hi_1;
	gate_adc_low_old_2 = gate_adc_low_2;
	gate_adc_hi_old_2 = gate_adc_hi_2;
	file_old=Current_beta_file;
	CutTag_old=CutModifyCounter;
	PreSetADC1_low_Thres_ch_old = PreSetADC1_low_Thres_ch;
	PreSetADC2_low_Thres_ch_old = PreSetADC2_low_Thres_ch;
	Last_BetaFastMode = BetaFastMode;
	AlphaMode_old = alphaMode;


	Match_hit.clear();	Match_hit.shrink_to_fit();
	match_event OneMatch;  // for one pair matching hits
	OneMatch.clear();

	hit_info candidate_tem;
	candidate_tem.empty=true;
	Long64_t progress=0; // show progress

	bool InE1Gate=true;
	bool InE2Gate=true;


	//%%%%%%%%%%%%% for read beta raw tree %%%%%%%%%%%%%%%%%%
	Long64_t time_tem =0;
	int channel_tem=0;
	int adc_tem=0;
	Long64_t sweeps_gclock_tem=0;
	Long64_t time_relative_tem=0;

	if(MakeBetaRawTree && tbeta_Raw!=NULL){

		nevt_bta = tbeta_Raw->GetEntries();

		tbeta_Raw->SetBranchAddress("time",&time_tem);
		tbeta_Raw->SetBranchAddress("channel",&channel_tem);
		tbeta_Raw->SetBranchAddress("adc",&adc_tem);
		tbeta_Raw->SetBranchAddress("sweeps_global",&sweeps_gclock_tem);
		tbeta_Raw->SetBranchAddress("time_relative",&time_relative_tem);
	}
	//%%%%%%%%%%%%%%%% end of for reading beta raw tree %%%%%%%%%%%%%%%


	for(Long64_t line_index=0;line_index<nevt_bta;line_index++){
		progress++;
		if(progress%1000==0)cout<<"progress = "<<progress<<" / "<<nevt_bta<<"\r"<<flush;
		//int progress_ratio =TMath::Nint((double)(progress/nevt_bta)*100);
		//if(progress_ratio%5==0)cout<<"\r"<<"progress= "<<progress_ratio<<"%"<<flush;

			if(MakeBetaRawTree){
				tbeta_Raw->GetEntry(line_index);
					candidate_tem.time = time_tem;         
					candidate_tem.time_relative = time_relative_tem;               
					candidate_tem.adc = adc_tem;                         
					candidate_tem.channel = channel_tem;                      
					candidate_tem.sweeps_gclock = sweeps_gclock_tem;              

			}else{
					candidate_tem.time = time_v->at(line_index);
					candidate_tem.time_relative = time_relative_v->at(line_index);
					candidate_tem.adc = adc_v->at(line_index);
					candidate_tem.channel = channel_v->at(line_index);
					candidate_tem.sweeps_gclock = sweeps_gclock_v->at(line_index);
			}



		if(MakeBetaCut_Or_FindSignal(-1000,-1000,1,true)){// beta noise cuts are available
			bool IsSignal=false;
			if(candidate_tem.channel==1){
				 IsSignal=MakeBetaCut_Or_FindSignal(candidate_tem.time_relative,ADC2keV_1(candidate_tem.adc),candidate_tem.channel);
				 if(!IsSignal) continue;
			}
			else if(candidate_tem.channel==2){
				IsSignal=MakeBetaCut_Or_FindSignal(candidate_tem.time_relative,ADC2keV_2(candidate_tem.adc),candidate_tem.channel);
				if(!IsSignal) continue;
			} 

		}else{
			if(!WarningPrinted)cout<<"\e[1;31m"<<"No Beta noise cut available yet!!!!!"<<"\e[0m"<<endl;
			WarningPrinted = true;
		}

		if(alphaMode && candidate_tem.channel != 1) continue; // Select only the ch1 event

		InE1Gate=true;
		InE2Gate=true;

		if(candidate_tem.channel==1) InE1Gate = candidate_tem.adc>=gate_adc_low_1 && candidate_tem.adc<=gate_adc_hi_1 ;    // adc gate are set in the unit of CHANNEL !! Can be set by SetCondiction(....)
		if(candidate_tem.channel==2) InE2Gate = candidate_tem.adc>=gate_adc_low_2 && candidate_tem.adc<=gate_adc_hi_2 ;

		if( InE1Gate && InE2Gate && candidate_tem.channel !=syn_CH && candidate_tem.sweeps_gclock!=0){ 
		// energy gate event for first waiting candidate
		// eliminate the event of sweeps_gclock ==0 ==> MCS of TOF not start yet.

			if(alphaMode){
							OneMatch.B1_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B1_time = candidate_tem.time;
							OneMatch.B1_time_relative = candidate_tem.time_relative;
							OneMatch.B1_adc = candidate_tem.adc;
							OneMatch.B2_sweeps_gclock = OneMatch.B1_sweeps_gclock;
							OneMatch.B2_time = OneMatch.B1_time;
							OneMatch.B2_time_relative = OneMatch.B1_time_relative;
							OneMatch.B2_adc = OneMatch.B1_adc;
							OneMatch.delta_t=(int)(OneMatch.B2_time-OneMatch.B1_time);//delta_time;
							Match_hit.push_back(OneMatch);
							OneMatch.clear();
							continue;

			}


			for(Long64_t forward_i=line_index+1;forward_i<nevt_bta;forward_i++){ // seeking coincident forward(towards end of file)
					Long64_t tem_sweeps_gclock= 0;
					Long64_t tem_time = 0;
					Long64_t delta_time = 0;
					Long64_t tem_time_relative = 0;
					int tem_adc = 0;
					int tem_channel = 0;

					if(MakeBetaRawTree && tbeta_Raw!=NULL){

						tbeta_Raw->GetEntry(forward_i);

						 tem_sweeps_gclock= sweeps_gclock_tem;
						 tem_time = time_tem;
						 delta_time = TMath::Abs(tem_time-candidate_tem.time);
						 tem_time_relative = time_relative_tem;
						 tem_adc = adc_tem;
						 tem_channel = channel_tem;
					}else{
						 tem_sweeps_gclock= sweeps_gclock_v->at(forward_i);
						 tem_time = time_v->at(forward_i);
						 delta_time = TMath::Abs(tem_time-candidate_tem.time);
						 tem_time_relative = time_relative_v->at(forward_i);
						 tem_adc = adc_v->at(forward_i);
						 tem_channel = channel_v->at(forward_i);
					}


				if(tem_channel == candidate_tem.channel||tem_channel==syn_CH){// sweeps trig or same silicon
					if(delta_time>gate_time)break;
					else continue; 
				}
				else{
if(tem_channel==1) InE1Gate = tem_adc>=gate_adc_low_1 && tem_adc<=gate_adc_hi_1;
if(tem_channel==2) InE2Gate = tem_adc>=gate_adc_low_2 && tem_adc<=gate_adc_hi_2;


			//%%%%%%%%%%%%%%%5 reject noise by cut  %%%%%%%%%%%%%%%%%%%%%%%%%
			if(MakeBetaCut_Or_FindSignal(-1000,-1000,1,true)){// beta noise cuts are available
				bool IsSignal=false;
				if(tem_channel==1){
					 IsSignal=MakeBetaCut_Or_FindSignal(tem_time_relative,ADC2keV_1(tem_adc),tem_channel);
					 if(!IsSignal) continue;
				}
				else if(tem_channel==2){
					IsSignal=MakeBetaCut_Or_FindSignal(tem_time_relative,ADC2keV_2(tem_adc),tem_channel);
					if(!IsSignal) continue;
				} 

			}else{
				if(!WarningPrinted) cout<<"\e[1;31m"<<"No Beta noise cut available yet!!!!!"<<"\e[0m"<<endl;
				WarningPrinted = true;
			}
			//%%%%%%%%%%%%%%  end of rejecting noise by cut %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


					if(delta_time<=gate_time && InE1Gate && InE2Gate){// matching
						if(candidate_tem.channel==1){// load data to match vector
							OneMatch.B1_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B1_time = candidate_tem.time;
							OneMatch.B1_time_relative = candidate_tem.time_relative;
							OneMatch.B1_adc = candidate_tem.adc;
							OneMatch.B2_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B2_time = tem_time;
							OneMatch.B2_time_relative = tem_time_relative;
							OneMatch.B2_adc = tem_adc;
							OneMatch.delta_t=(int)(OneMatch.B2_time-OneMatch.B1_time);//delta_time;
						}
						else{
							OneMatch.B1_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B1_time = tem_time;
							OneMatch.B1_time_relative = tem_time_relative;
							OneMatch.B1_adc = tem_adc;
							OneMatch.B2_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B2_time = candidate_tem.time;
							OneMatch.B2_time_relative = candidate_tem.time_relative;
							OneMatch.B2_adc = candidate_tem.adc;
							OneMatch.delta_t=(int)(OneMatch.B2_time-OneMatch.B1_time);//delta_time;
						}
						KeepNewMatching(OneMatch);
						OneMatch.clear();
					}
					else{
						if(delta_time>gate_time)break;
						continue;
					}
				}

			}// loop forward


/*			for(int backward_i=line_index-1;backward_i>=0;backward_i--){ // seeking coincident backward(towards beginning of file)
				Long64_t tem_sweeps_gclock= sweeps_gclock_v->at(backward_i);
				Long64_t tem_time = time_v->at(backward_i);
				Long64_t delta_time = TMath::Abs(tem_time-candidate_tem.time);
				int tem_adc = adc_v->at(backward_i);
				int tem_channel = channel_v->at(backward_i);
				if(tem_channel == candidate_tem.channel||tem_channel==syn_CH){// sweeps trig or same silicon
					if(delta_time>gate_time)break;
					else continue; 
				}
				else{
if(tem_channel==1) InE1Gate = tem_adc>=gate_adc_low_1 && tem_adc<=gate_adc_hi_1;
if(tem_channel==2) InE2Gate = tem_adc>=gate_adc_low_2 && tem_adc<=gate_adc_hi_2;

					if(delta_time<=gate_time && InE1Gate && InE2Gate){// matching
						if(candidate_tem.channel==1){// load data to match vector
							OneMatch.B1_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B1_time = candidate_tem.time;
							OneMatch.B1_adc = candidate_tem.adc;
							OneMatch.B2_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B2_time = tem_time;
							OneMatch.B2_adc = tem_adc;
							OneMatch.delta_t=delta_time;
						}
						else{
							OneMatch.B1_sweeps_gclock = tem_sweeps_gclock;
							OneMatch.B1_time = tem_time;
							OneMatch.B1_adc = tem_adc;
							OneMatch.B2_sweeps_gclock = candidate_tem.sweeps_gclock;
							OneMatch.B2_time = candidate_tem.time;
							OneMatch.B2_adc = candidate_tem.adc;
							OneMatch.delta_t=delta_time;
						}
						KeepNewMatching(OneMatch);
						OneMatch.clear();
					}
					else{
						if(delta_time>gate_time)break;
						continue;
					}
				}

			}// loop backward   */


		}
		else{continue;}


	}// end of for loop line_index of file


	if(MakeBetaRawTree) tbeta_Raw->ResetBranchAddresses();


	/*for(unsigned long result_i=0;result_i<Match_hit.size();result_i++){
		Long64_t B1_t = Match_hit[result_i].B1_time;
		int B1_adc = Match_hit[result_i].B1_adc;
		Long64_t B2_t = Match_hit[result_i].B2_time;
		int B2_adc = Match_hit[result_i].B2_adc;
		printf("%lld,%d \t %lld,%d \t %d\n",B1_t,B1_adc,B2_t,B2_adc,Match_hit[result_i].delta_t);
	}*/
	cout<<endl;
	cout<<"total event before refine= "<<Match_hit.size()<<endl;
	cout<<endl;

	Match_hit.shrink_to_fit();

	if(!alphaMode && BetaBetaFilter){
		vector<match_event> Match_hit_refine;

		RefineMatchHit(Match_hit,Match_hit_refine);


		Match_hit.clear();Match_hit.shrink_to_fit();

		Match_hit = Match_hit_refine;

		Match_hit.shrink_to_fit();
	}

/*
cout<<endl;
	for(unsigned long result_i=0;result_i<Match_hit.size();result_i++){
		Long64_t B1_t = Match_hit[result_i].B1_time;
		int B1_adc = Match_hit[result_i].B1_adc;
		Long64_t B2_t = Match_hit[result_i].B2_time;
		int B2_adc = Match_hit[result_i].B2_adc;
		printf("%lld,%d \t %lld,%d \t %d\n",B1_t,B1_adc,B2_t,B2_adc,Match_hit[result_i].delta_t);
	}*/

	if(!alphaMode) cout<<"total event after refine= "<<Match_hit.size()<<endl;
	else cout<<"Alpha mode:  total event at CH1 = "<<Match_hit.size()<<endl;
	cout<<endl;

	ShowCoinCondition();


	if(BetaFastMode){
		for(Long64_t index=0;index<nevt_bta;index++){ // before coincident
			if(channel_v->at(index)==1) {
				if( IsBetaFromNewFile(time_v->at(index),-1) ){
					histo_h1E->Fill(adc_v->at(index));
					h1_t_relate_raw->Fill(time_relative_v->at(index));
					h1_t_relate_EkeV_raw->Fill( ADC2keV_1(adc_v->at(index)) , time_relative_v->at(index) );
				}
			}
			if(channel_v->at(index)==2){
				if( IsBetaFromNewFile(-1,time_v->at(index)) ){
					histo_h2E->Fill(adc_v->at(index));
					h2_t_relate_raw->Fill(time_relative_v->at(index));
					h2_t_relate_EkeV_raw->Fill( ADC2keV_2(adc_v->at(index)) , time_relative_v->at(index) );
				}
			} 
		}

		for(unsigned long index=0;index<Match_hit.size();index++){ // after coincident
			if( IsBetaFromNewFile(Match_hit[index].B1_time,Match_hit[index].B2_time) ){
				double E_B1 =  ADC2keV_1(Match_hit[index].B1_adc);
				double E_B2 =  ADC2keV_2(Match_hit[index].B2_adc);

				histo_beta_deltaT->Fill(Match_hit[index].delta_t);
				histo_beta_EkeV->Fill( E_B1,  E_B2);
				h1_t_relate_EkeV_coin->Fill( E_B1, Match_hit[index].B1_time_relative );
				h2_t_relate_EkeV_coin->Fill( E_B2, Match_hit[index].B2_time_relative );
				h1_t_relate_coin->Fill(Match_hit[index].B1_time_relative);
				h2_t_relate_coin->Fill(Match_hit[index].B2_time_relative);
				h_LogEratio_delta_t->Fill(Match_hit[index].delta_t,TMath::Log10(E_B2/E_B1));
			}

		}

	}




	return true;

}




Long64_t DownScale=100; // for beta raw histogram display; get data per DownScale events
bool StopFillBetaRawHisto = false;
void SetStopFillBetaRawHisto(int sig){
	StopFillBetaRawHisto =true;
}

void StopShowBeta_Beta(){
	StopFillBetaRawHisto =true;
}


bool Beta_Graph_Initialize(){//SetBeta_E_Calibrate_para(2,0,2.8570); //  remove when recalibrate	
	//if(c_beta_beta!=NULL) delete c_beta_beta;
	if(histo_beta_deltaT!=NULL) delete histo_beta_deltaT;
	if(histo_beta_EkeV!=NULL) delete histo_beta_EkeV;
	if(histo_beta_EkeV_withcut!=NULL) delete histo_beta_EkeV_withcut;
	if(histo_beta_EkeV_CutAndRej!=NULL) delete histo_beta_EkeV_CutAndRej;
	if(histo_h1E!=NULL) delete histo_h1E;
	if(histo_h2E!=NULL) delete histo_h2E;

if(h1_t_relate_raw!=NULL) delete h1_t_relate_raw;
if(h1_t_relate_coin!=NULL)delete h1_t_relate_coin;
if(h2_t_relate_raw!=NULL) delete h2_t_relate_raw;
if(h2_t_relate_coin!=NULL)delete h2_t_relate_coin;
if(h1_t_relate_EkeV_raw!=NULL)delete h1_t_relate_EkeV_raw;
if(h2_t_relate_EkeV_raw!=NULL)delete h2_t_relate_EkeV_raw;
if(h1_t_relate_EkeV_coin!=NULL)delete h1_t_relate_EkeV_coin;
if(h2_t_relate_EkeV_coin!=NULL)delete h2_t_relate_EkeV_coin;
if(h_LogEratio_delta_t != NULL) delete h_LogEratio_delta_t;



	if(tbeta_beta!=NULL) delete tbeta_beta;


if(c_beta_coin_and_raw==NULL || c_beta_coin_and_raw->GetCanvasImp()==nullptr){
	c_beta_coin_and_raw = new TCanvas("c_beta_coin_and_raw","c_beta_coin_and_raw",2000,2200);
	c_beta_coin_and_raw->SetCanvasSize(1800,2200);
	c_beta_coin_and_raw->Divide(2,5);

}

	if(c_beta_beta==NULL || c_beta_beta->GetCanvasImp()==nullptr){
		c_beta_beta = new TCanvas("c_beta_beta","c_beta_beta",2000,1200);
		c_beta_beta->Divide(2,2);
	}
	
	c_beta_beta->cd(1);

	if(f_beta_tree!=NULL){
		if(f_beta_tree->IsOpen()){f_beta_tree->Close(); delete f_beta_tree;}
		else{delete f_beta_tree;}
	}

	//************ create file to store beta-beta tree, Avoid crash!!!!  ****************
	string tem_path = Path_beta_file +"../rootfiles/";
	bool PathOK = gSystem->AccessPathName(tem_path.c_str());
	if(!PathOK){ // path exist
		tem_path += "tbeta_beta_tree.root";
		f_beta_tree = new TFile(tem_path.c_str(),"RECREATE");
	}
	else{
		tem_path = Path_beta_file + "tbeta_beta_tree.root";
		if(gSystem->AccessPathName(Path_beta_file.c_str())){
			cout<<" Path--> "<<Path_beta_file<<" is not exist. Not able to save tbeta_beta_tree.root!! Abort!!"<<endl;
			return false;
		}
		f_beta_tree = new TFile(tem_path.c_str(),"RECREATE");
	}

	if(f_beta_tree->IsOpen()){f_beta_tree->cd();}
	else{cout<<"\e[1;33m"<<"Fail to create file to store beta-beta tree!!!! Abort!!!"<<"\e[0m"<<endl; return false;}

	tbeta_beta = new TTree("tbeta_beta","tbeta_beta");

	histo_beta_deltaT = new TH1D("histo_beta_deltaT","Coincidence time difference(Si2-Si1)",gate_time/40,Tgate_cent_offset-gate_time,Tgate_cent_offset+gate_time);//gate_time/40,0,gate_time

	h_LogEratio_delta_t = new TH2D("h_LogEratio_delta_t","TMath::Log10(EkeV_2/EkeV_1):delta_t",gate_time/40,Tgate_cent_offset-gate_time,Tgate_cent_offset+gate_time,1000,-1.8,1.1);

	if(beta_slope[0]==0 && beta_slope[1]==0){
			SetBeta_E_Calibrate_para();
	}

	int ADC_dis_resolution = (int)((0.0028/10)*ADC_max_ch);//0.0098
	if(ADC_dis_resolution==0)ADC_dis_resolution=1;
	int Nbins_ADC=ADC_max_ch/ADC_dis_resolution;


	histo_beta_EkeV = new TH2D("histo_beta_EkeV","",Nbins_ADC,0.,ADC2keV_1(ADC_max_ch),Nbins_ADC,0.,ADC2keV_2(ADC_max_ch));
	//title "coincident Energy Si1 VS Si2" not set for easy to use cut
	histo_beta_EkeV->GetXaxis()->SetTitle("Si_1 [keV]"); histo_beta_EkeV->GetXaxis()->CenterTitle();
	histo_beta_EkeV->GetYaxis()->SetTitle("Si_2 [keV]"); histo_beta_EkeV->GetYaxis()->CenterTitle();

	histo_beta_EkeV_withcut = new TH2D("histo_beta_EkeV_withcut","EkeV_2:EkeV_1",Nbins_ADC,0.,ADC2keV_1(ADC_max_ch),Nbins_ADC,0.,ADC2keV_2(ADC_max_ch));  // show E-E 2D histo after cut
	histo_beta_EkeV_withcut->GetXaxis()->SetTitle("Si_1 [keV]"); histo_beta_EkeV_withcut->GetXaxis()->CenterTitle();
	histo_beta_EkeV_withcut->GetYaxis()->SetTitle("Si_2 [keV]"); histo_beta_EkeV_withcut->GetYaxis()->CenterTitle();

	histo_beta_EkeV_CutAndRej = new TH2D("histo_beta_EkeV_CutAndRej","EkeV_2:EkeV_1",Nbins_ADC,0.,ADC2keV_1(ADC_max_ch),Nbins_ADC,0.,ADC2keV_2(ADC_max_ch)); 
	histo_beta_EkeV_CutAndRej->GetXaxis()->SetTitle("Si_1 [keV]"); histo_beta_EkeV_CutAndRej->GetXaxis()->CenterTitle();
	histo_beta_EkeV_CutAndRej->GetYaxis()->SetTitle("Si_2 [keV]"); histo_beta_EkeV_CutAndRej->GetYaxis()->CenterTitle();


	h1_t_relate_EkeV_coin = new TH2D("h1_t_relate_EkeV_coin","time1_relative:EkeV_1",Nbins_ADC,0,ADC2keV_1(ADC_max_ch),25000,0,50e6);

	h2_t_relate_EkeV_coin = new TH2D("h2_t_relate_EkeV_coin","time2_relative:EkeV_2",Nbins_ADC,0,ADC2keV_2(ADC_max_ch),25000,0,50e6);

	//if(histo_h1E!=NULL) delete histo_h1E;
		histo_h1E = new TH1D("histo_h1E"," Si1 and Si2 ADC ",Nbins_ADC,0,ADC_max_ch); histo_h1E->SetLineColor(kBlack);
	//if(histo_h2E!=NULL) delete histo_h2E;
		histo_h2E = new TH1D("histo_h2E"," Si1 and Si2 ADC ",Nbins_ADC,0,ADC_max_ch); histo_h2E->SetLineColor(kRed);

h1_t_relate_raw = new TH1D("h1_t_relate_raw","time1_relative",25000,-1e6,50e6);
h1_t_relate_coin = new TH1D("h1_t_relate_coin","time1_relative",25000,-1e6,50e6);

h2_t_relate_raw = new TH1D("h2_t_relate_raw","time2_relative",25000,-1e6,50e6);
h2_t_relate_coin = new TH1D("h2_t_relate_coin","time2_relative",25000,-1e6,50e6);


h1_t_relate_EkeV_raw = new TH2D("h1_t_relate_EkeV_raw","time1_relative:EkeV_1",Nbins_ADC,0,ADC2keV_1(ADC_max_ch),2500,0,50e6);
h2_t_relate_EkeV_raw = new TH2D("h2_t_relate_EkeV_raw","time2_relative:EkeV_2",Nbins_ADC,0,ADC2keV_2(ADC_max_ch),2500,0,50e6);


if(!BetaFastMode){
	if(MakeBetaRawTree){
		if(tbeta_Raw!=NULL){

			Long64_t time_relative_tem=0;
			int channel_tem=0;
			int adc_tem=0;

			nevt_bta = tbeta_Raw->GetEntries();
			tbeta_Raw->SetBranchAddress("time_relative",&time_relative_tem);
			tbeta_Raw->SetBranchAddress("channel",&channel_tem);
			tbeta_Raw->SetBranchAddress("adc",&adc_tem);

			cout<<endl;
			cout<<"\e[1;33m Pending....... Maybe take a long time to plot raw beta histograms.  Using ctrl + c (otherwise use: StopShowBeta_Beta() )to stop and change DownScale factor to speed up! \e[0m"<<endl;

			StopFillBetaRawHisto = false; 

			if(DownScale>0){
					for(Long64_t index=0;index<nevt_bta;index++){
						signal(SIGINT,SetStopFillBetaRawHisto);
						gSystem->ProcessEvents();
						if(StopFillBetaRawHisto) break;

						if(index % 10000 ==0) cout<<"Process = "<<index<<" / "<<nevt_bta<<" ("<<index*1.0/nevt_bta*100<<"%) \r"<<flush;
						if(index % DownScale !=0) continue;
						tbeta_Raw->GetEntry(index);
						if(channel_tem==1){
							histo_h1E->Fill(adc_tem);
							h1_t_relate_raw->Fill(time_relative_tem);
							h1_t_relate_EkeV_raw->Fill( ADC2keV_1(adc_tem) , time_relative_tem );

						}

						if(channel_tem==2){
							histo_h2E->Fill(adc_tem);
							h2_t_relate_raw->Fill(time_relative_tem);
							h2_t_relate_EkeV_raw->Fill( ADC2keV_2(adc_tem) , time_relative_tem );
						}


					}
			}

			tbeta_Raw->ResetBranchAddresses();
			StopFillBetaRawHisto = false; // reset


			/*

			c_beta_beta->cd(1);
			tbeta_Raw->Draw("adc>>histo_h1E","channel==1");
			tbeta_Raw->Draw("adc>>histo_h2E","channel==2");

			c_beta_coin_and_raw->cd(6);
			tbeta_Raw->Draw("time_relative>>h1_t_relate_raw","channel==1");
			c_beta_coin_and_raw->cd(8);
			tbeta_Raw->Draw("time_relative>>h2_t_relate_raw","channel==2");
			c_beta_coin_and_raw->cd(2);
			tbeta_Raw->Draw("time_relative:ADC2keV_1(adc)>>h1_t_relate_EkeV_raw","channel==1","colz");
			c_beta_coin_and_raw->cd(4);
			tbeta_Raw->Draw("time_relative:ADC2keV_2(adc)>>h2_t_relate_EkeV_raw","channel==2","colz");
			*/

		}

	}else{
			cout<<endl;
			cout<<"\e[1;33m Pending....... Maybe take a long time to plot raw bata histograms.  Using ctrl + c (otherwise use: StopShowBeta_Beta() ) to stop and change DownScale factor to speed up! \e[0m"<<endl;

			StopFillBetaRawHisto = false; 

			if(DownScale>0){
				for(Long64_t index=0;index<nevt_bta;index++){
					gSystem->ProcessEvents();
					signal(SIGINT,SetStopFillBetaRawHisto);
					if(StopFillBetaRawHisto) break;

					if(index % 10000 ==0) cout<<"Process = "<<index<<" / "<<nevt_bta<<" ("<<index*1.0/nevt_bta*100<<"%) \r"<<flush;
					if(index % DownScale !=0) continue;

					if(channel_v->at(index)==1) {
						histo_h1E->Fill(adc_v->at(index));
						h1_t_relate_raw->Fill(time_relative_v->at(index));
						h1_t_relate_EkeV_raw->Fill( ADC2keV_1(adc_v->at(index)) , time_relative_v->at(index) );
					}
					if(channel_v->at(index)==2){
						histo_h2E->Fill(adc_v->at(index));
						h2_t_relate_raw->Fill(time_relative_v->at(index));
						h2_t_relate_EkeV_raw->Fill( ADC2keV_2(adc_v->at(index)) , time_relative_v->at(index) );
					} 
				}

			}

			StopFillBetaRawHisto = false; // reset
	}
}

	return true;

}


void PrintFrequentCommand(){
	printf("tbeta_beta->Draw(\"(TMath::Log(EkeV_2/EkeV_1)):time1_relative>>ht1(1000,0,30000000,100,-4,4)\",\"!cut1&&!cut2&&!cut3&&!cut4&&!cut5 && !cut6 && !cut7 && !cut8 && !cutpure && delta_t<600 && delta_t>-2000\",\"colz\") \n");
	printf("tbeta_beta->Draw(\"EkeV_2-EkeV_1:delta_t>>het(500,-30000,30000,200,-4000,4000)\",\"!cut1 && !cut2 && !cut3 && !cut4 && !cut5 && !cut6 && !cut7 && !cut8\",\"colz\") \n");
}

bool ShowBeta_Beta_Coin(){

	if(!BetaFastMode){
		if(!Beta_Graph_Initialize()){
			cout<<"\e[1;33m"<<"Fail to initialize the Beta_Beta_coin Graph beacaue of failure of creating file to save tbeta_beta tree!!"<<"\e[0m"<<endl;
			return false;
		}
	}

	c_beta_beta->cd(1);
	histo_h2E->Draw();
	c_beta_beta->cd(1)->Modified();
	c_beta_beta->cd(1)->Update();
	TPaveStats* st2 = (TPaveStats*) histo_h2E->FindObject("stats");
	//cout<<"hii ="<<st2<<endl;
	st2->SetY1NDC(0.73500001-0.21);
	st2->SetY2NDC(0.93500000-0.21);
	histo_h1E->Draw();
	histo_h2E->Draw("same");
	c_beta_beta->cd(1)->SetLogy();

	if(Match_hit.size()==0){ cout<<"\e[1;23m"<<"No coincidence in store"<<"\e[0m"<<endl; return false;}


	Long64_t time1=0; Long64_t time1_relative=0; Long64_t sweeps_gclock_1=0; int adc1=0; double EkeV_1=0;
	Long64_t time2=0; Long64_t time2_relative=0; Long64_t sweeps_gclock_2=0; int adc2=0; double EkeV_2=0;
	int delta_t=0;


	tbeta_beta->Branch("sweeps_gclock_1",&sweeps_gclock_1);
	tbeta_beta->Branch("time1",&time1);
	tbeta_beta->Branch("time1_relative",&time1_relative);
	tbeta_beta->Branch("adc1",&adc1);
	tbeta_beta->Branch("EkeV_1",&EkeV_1);
	tbeta_beta->Branch("sweeps_gclock_2",&sweeps_gclock_2);
	tbeta_beta->Branch("time2",&time2);
	tbeta_beta->Branch("time2_relative",&time2_relative);
	tbeta_beta->Branch("adc2",&adc2);
	tbeta_beta->Branch("EkeV_2",&EkeV_2);
	tbeta_beta->Branch("delta_t",&delta_t);

	for(unsigned long index=0;index< Match_hit.size();index++){
			// histo_beta_deltaT->Fill(Match_hit[index].delta_t);	
		sweeps_gclock_1 = Match_hit[index].B1_sweeps_gclock;
		time1 = Match_hit[index].B1_time;
		time1_relative = Match_hit[index].B1_time_relative;
		adc1 = Match_hit[index].B1_adc;
		EkeV_1 = ADC2keV_1(adc1);
		sweeps_gclock_2 = Match_hit[index].B2_sweeps_gclock;
		time2 = Match_hit[index].B2_time;
		time2_relative = Match_hit[index].B2_time_relative;
		adc2 = Match_hit[index].B2_adc;
		EkeV_2 = ADC2keV_2(adc2);
		delta_t = Match_hit[index].delta_t;
		tbeta_beta->Fill();
	}

	c_beta_beta->cd(3);
if(!BetaFastMode)	tbeta_beta->Draw("delta_t>>histo_beta_deltaT");
else histo_beta_deltaT->Draw();

	c_beta_beta->cd(2);
if(!BetaFastMode)	tbeta_beta->Draw("EkeV_2:EkeV_1>>histo_beta_EkeV","","colz");
else histo_beta_EkeV->Draw("colz");

	histo_beta_EkeV->SetTitle("EkeV_2:EkeV_1");
	exdefaultcolor->Draw();

	cout<<"beta-beta tree is created:  tbeta_beta"<<endl;
	cout<<endl;

  // Save the current global palette
    
const TArrayI &savedPalette = TColor::GetPalette();
Int_t nColors = savedPalette.GetSize();
Int_t* gColorArray_cp = new Int_t[nColors];
const Int_t* currentColorArray = savedPalette.GetArray();

for(int i=0;i<nColors;i++){
	gColorArray_cp[i] = currentColorArray[i];
}

//gStyle->SetPalette(kPastel);

//TExec * exbcolor = new TExec("exbcolor","gStyle->SetPalette(kPastel);");

	c_beta_coin_and_raw->cd(1);

if(!BetaFastMode)	tbeta_beta->Draw("time1_relative:EkeV_1>>h1_t_relate_EkeV_coin","","colz");
else h1_t_relate_EkeV_coin->Draw("colz");

	h1_t_relate_EkeV_coin->SetTitle("time1_relative:EkeV_1");
	h1_t_relate_EkeV_coin->GetXaxis()->SetTitle("Si1 Energy [keV]");
	h1_t_relate_EkeV_coin->GetXaxis()->CenterTitle();
	h1_t_relate_EkeV_coin->GetXaxis()->SetLabelSize(0.045);
	h1_t_relate_EkeV_coin->GetYaxis()->SetTitle("Si1 (w/ coin) relative time (time to sync. signal) [ns]");
	h1_t_relate_EkeV_coin->GetYaxis()->CenterTitle();
	h1_t_relate_EkeV_coin->GetYaxis()->SetLabelSize(0.045);
	h1_t_relate_EkeV_coin->GetXaxis()->SetTitleSize(0.045);
	h1_t_relate_EkeV_coin->GetYaxis()->SetTitleSize(0.045);
exbcolor->Draw();

	c_beta_coin_and_raw->cd(2);

	h1_t_relate_EkeV_raw->Draw("colz");
	h1_t_relate_EkeV_raw->SetTitle("time1_relative:EkeV_1");
	h1_t_relate_EkeV_raw->GetXaxis()->SetTitle("Si1 Energy [keV]");
	h1_t_relate_EkeV_raw->GetXaxis()->CenterTitle();
	h1_t_relate_EkeV_raw->GetXaxis()->SetLabelSize(0.045);
	h1_t_relate_EkeV_raw->GetYaxis()->SetTitle("Si1 (w/o coin) relative time (time to sync. signal) [ns]");
	h1_t_relate_EkeV_raw->GetYaxis()->CenterTitle();
	h1_t_relate_EkeV_raw->GetYaxis()->SetLabelSize(0.045);
	h1_t_relate_EkeV_raw->GetXaxis()->SetTitleSize(0.045);
	h1_t_relate_EkeV_raw->GetYaxis()->SetTitleSize(0.045);
exbcolor->Draw();

	c_beta_coin_and_raw->cd(3);

if(!BetaFastMode)	tbeta_beta->Draw("time2_relative:EkeV_2>>h2_t_relate_EkeV_coin","","colz");
else h2_t_relate_EkeV_coin->Draw("colz");
	h2_t_relate_EkeV_coin->SetTitle("time2_relative:EkeV_2");
	h2_t_relate_EkeV_coin->GetXaxis()->SetTitle("Si2 Energy [keV]");
	h2_t_relate_EkeV_coin->GetXaxis()->CenterTitle();
	h2_t_relate_EkeV_coin->GetXaxis()->SetLabelSize(0.045);
	h2_t_relate_EkeV_coin->GetYaxis()->SetTitle("Si2 (w/ coin) relative time (time to sync. signal) [ns]");
	h2_t_relate_EkeV_coin->GetYaxis()->CenterTitle();
	h2_t_relate_EkeV_coin->GetYaxis()->SetLabelSize(0.045);
	h2_t_relate_EkeV_coin->GetXaxis()->SetTitleSize(0.045);
	h2_t_relate_EkeV_coin->GetYaxis()->SetTitleSize(0.045);
exbcolor->Draw();

	c_beta_coin_and_raw->cd(4);

	h2_t_relate_EkeV_raw->Draw("colz");
	h2_t_relate_EkeV_raw->SetTitle("time2_relative:EkeV_2");
	h2_t_relate_EkeV_raw->GetXaxis()->SetTitle("Si2 Energy [keV]");
	h2_t_relate_EkeV_raw->GetXaxis()->CenterTitle();
	h2_t_relate_EkeV_raw->GetXaxis()->SetLabelSize(0.045);
	h2_t_relate_EkeV_raw->GetYaxis()->SetTitle("Si2 (w/o coin) relative time (time to sync. signal) [ns]");
	h2_t_relate_EkeV_raw->GetYaxis()->CenterTitle();
	h2_t_relate_EkeV_raw->GetYaxis()->SetLabelSize(0.045);
	h2_t_relate_EkeV_raw->GetXaxis()->SetTitleSize(0.045);
	h2_t_relate_EkeV_raw->GetYaxis()->SetTitleSize(0.045);
exbcolor->Draw();

	c_beta_coin_and_raw->cd(5);
if(!BetaFastMode)	tbeta_beta->Draw("time1_relative>>h1_t_relate_coin");
else h1_t_relate_coin->Draw();
	h1_t_relate_coin->GetXaxis()->SetTitle("Si1 w/ coin. relative time [ns]");
	h1_t_relate_coin->GetXaxis()->CenterTitle();
	h1_t_relate_coin->GetXaxis()->SetLabelSize(0.045);
	h1_t_relate_coin->GetXaxis()->SetTitleSize(0.045);



	c_beta_coin_and_raw->cd(6);
	h1_t_relate_raw->Draw();
	h1_t_relate_raw->GetXaxis()->SetTitle("Si1 w/o coin. relative time [ns]");
	h1_t_relate_raw->GetXaxis()->CenterTitle();
	h1_t_relate_raw->GetXaxis()->SetLabelSize(0.045);
	h1_t_relate_raw->GetXaxis()->SetTitleSize(0.045);
	

	c_beta_coin_and_raw->cd(7);
if(!BetaFastMode)	tbeta_beta->Draw("time2_relative>>h2_t_relate_coin");
else h2_t_relate_coin->Draw();
	h2_t_relate_coin->GetXaxis()->SetTitle("Si2 w/ coin. relative time [ns]");
	h2_t_relate_coin->GetXaxis()->CenterTitle();
	h2_t_relate_coin->GetXaxis()->SetLabelSize(0.045);
	h2_t_relate_coin->GetXaxis()->SetTitleSize(0.045);
	

	c_beta_coin_and_raw->cd(8);
	h2_t_relate_raw->Draw();
	h2_t_relate_raw->GetXaxis()->SetTitle("Si2 w/o coin. relative time [ns]");
	h2_t_relate_raw->GetXaxis()->CenterTitle();
	h2_t_relate_raw->GetXaxis()->SetLabelSize(0.045);
	h2_t_relate_raw->GetXaxis()->SetTitleSize(0.045);

	c_beta_coin_and_raw->cd(9);
if(!BetaFastMode)	tbeta_beta->Draw("TMath::Log10(EkeV_2/EkeV_1):delta_t>>h_LogEratio_delta_t","","colz");
else h_LogEratio_delta_t->Draw("colz");
	//h_LogEratio_delta_t->SetTitle("time2_relative:EkeV_2");
	h_LogEratio_delta_t->GetXaxis()->SetTitle("delta_t [ns]");
	h_LogEratio_delta_t->GetXaxis()->CenterTitle();
	h_LogEratio_delta_t->GetXaxis()->SetLabelSize(0.045);
	h_LogEratio_delta_t->GetYaxis()->SetTitle("Log10(E2/E1)");
	h_LogEratio_delta_t->GetYaxis()->CenterTitle();
	h_LogEratio_delta_t->GetYaxis()->SetLabelSize(0.045);
	h_LogEratio_delta_t->GetXaxis()->SetTitleSize(0.045);
	h_LogEratio_delta_t->GetYaxis()->SetTitleSize(0.045);
exbcolor->Draw();
	//gStyle->SetPalette(nColors, gColorArray_cp);

	c_beta_coin_and_raw->Modified();
	c_beta_coin_and_raw->Update();

	f_beta_tree->cd();
	histo_h1E->Write();
	histo_h2E->Write();
	tbeta_beta->Write();
	f_beta_tree->Write();
	tbeta_beta->ResetBranchAddresses();
PrintFrequentCommand();
	return true;

}


TCanvas* c_eject=NULL;
TGraph* Eject0[2] = {NULL,NULL};  // trigger leading edge ==>0; lagging edge ==>1
TGraph* Eject0_o[2] = {NULL,NULL};
TGraph* Eject1[2] = {NULL,NULL};
TGraph* Eject1_o[2] = {NULL,NULL};

void ShowEjectionBeta(double ejection0, double ejection1){  // ejection moment in ns

	ejection0 -= 2000; // mirror open 2us before daq delay
	ejection1 -=2000;

	if(c_eject==NULL){
		c_eject = new TCanvas("c_eject","Beta-Beta coincidence around ejection moment", 1200,800);
		c_eject->Divide(3,2);
	}

	if(Eject0[0]!=NULL){
	 	delete Eject0[0]; delete Eject0[1];
	 	delete Eject0_o[0]; delete Eject0_o[1];
	 	delete Eject1[0]; delete Eject1[1];
	 	delete Eject1_o[0]; delete Eject1_o[1];
	}

	double ejeW = 500 * 1000; // width of ejection trigger ==> 500 us to ns

	Eject0_o[0] = new TGraph(); Eject0_o[1] = new TGraph();  // point of eject
	Eject1_o[0] = new TGraph(); Eject1_o[1] = new TGraph();

	Eject0[0] = new TGraph(); Eject0[1] = new TGraph(); // points of beta-beta
	Eject1[0] = new TGraph(); Eject1[1] = new TGraph();

	for(int index=0;index<2;index++){
		Eject0_o[index]->SetPoint(0,ejection0+index*ejeW,ejection0+index*ejeW);
		if(index==0)Eject0_o[index]->SetPoint(1,ejection0+index*ejeW+3000,ejection0+index*ejeW+3000);
		else Eject0_o[index]->SetPoint(1,ejection0+index*ejeW+4000,ejection0+index*ejeW+4000);
		Eject0_o[index]->SetMarkerStyle(20);
		Eject0_o[index]->SetMarkerSize(1.1);

		if(index==0){
			Eject0_o[index]->SetTitle("time2_relative:time1_relative");
			Eject0_o[index]->GetXaxis()->SetTitle("Si_1 eje0 rising [ns]"); Eject0_o[index]->GetXaxis()->CenterTitle(); 
			Eject0_o[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-60000,ejection0+index*ejeW+60000); //-1000  +6000 original
			Eject0_o[index]->GetYaxis()->SetTitle("Si_2 eje0 rising [ns]"); Eject0_o[index]->GetYaxis()->CenterTitle();
			Eject0_o[index]->SetMaximum(ejection0+index*ejeW+60000);
			Eject0_o[index]->SetMinimum(ejection0+index*ejeW-60000);

			Eject0[index]->SetTitle("time2_relative:time1_relative");
			Eject0[index]->GetXaxis()->SetTitle("Si_1 eje0 rising [ns]"); Eject0[index]->GetXaxis()->CenterTitle();
			Eject0[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-60000,ejection0+index*ejeW+60000);
			Eject0[index]->GetYaxis()->SetTitle("Si_2 eje0 rising [ns]"); Eject0[index]->GetYaxis()->CenterTitle();
			Eject0[index]->SetMaximum(ejection0+index*ejeW+60000);
			Eject0[index]->SetMinimum(ejection0+index*ejeW-60000);
			Eject0[index]->SetMarkerStyle(31);
			Eject0[index]->SetMarkerSize(1.1);
			Eject0[index]->SetMarkerColor(kRed);
		}
		else{
			Eject0_o[index]->SetTitle("time2_relative:time1_relative");
			Eject0_o[index]->GetXaxis()->SetTitle("Si_1 eje0 dropping [ns]"); Eject0_o[index]->GetXaxis()->CenterTitle();
			Eject0_o[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-60000,ejection0+index*ejeW+60000);
			Eject0_o[index]->GetYaxis()->SetTitle("Si_2 eje0 dropping [ns]"); Eject0_o[index]->GetYaxis()->CenterTitle();
			Eject0_o[index]->SetMaximum(ejection0+index*ejeW+60000);
			Eject0_o[index]->SetMinimum(ejection0+index*ejeW-60000);

			Eject0[index]->SetTitle("time2_relative:time1_relative");
			Eject0[index]->GetXaxis()->SetTitle("Si_1 eje0 dropping [ns]"); Eject0[index]->GetXaxis()->CenterTitle();
			Eject0[index]->GetXaxis()->SetLimits(ejection0+index*ejeW-60000,ejection0+index*ejeW+60000);
			Eject0[index]->GetYaxis()->SetTitle("Si_2 eje0 dropping [ns]"); Eject0[index]->GetYaxis()->CenterTitle();
			Eject0[index]->SetMaximum(ejection0+index*ejeW+60000);
			Eject0[index]->SetMinimum(ejection0+index*ejeW-60000);
			Eject0[index]->SetMarkerStyle(31);
			Eject0[index]->SetMarkerSize(1.1);
			Eject0[index]->SetMarkerColor(kRed);
		}

		Eject1_o[index]->SetPoint(0,ejection1+index*ejeW,ejection1+index*ejeW);
		if(index==0)Eject1_o[index]->SetPoint(1,ejection1+index*ejeW+3000,ejection1+index*ejeW+3000);
		else Eject1_o[index]->SetPoint(1,ejection1+index*ejeW+4000,ejection1+index*ejeW+4000);
		Eject1_o[index]->SetMarkerStyle(20);
		Eject1_o[index]->SetMarkerSize(1.1);

		if(index==0){	
			Eject1_o[index]->SetTitle("time2_relative:time1_relative"); 
			Eject1_o[index]->GetXaxis()->SetTitle("Si_1 eje1 rising [ns]"); Eject1_o[index]->GetXaxis()->CenterTitle(); 
			Eject1_o[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-60000,ejection1+index*ejeW+60000);
			Eject1_o[index]->GetYaxis()->SetTitle("Si_2 eje1 rising [ns]"); Eject1_o[index]->GetYaxis()->CenterTitle();
			Eject1_o[index]->SetMaximum(ejection1+index*ejeW+60000);
			Eject1_o[index]->SetMinimum(ejection1+index*ejeW-60000);

			Eject1[index]->SetTitle("time2_relative:time1_relative");
			Eject1[index]->GetXaxis()->SetTitle("Si_1 eje1 rising [ns]"); Eject1[index]->GetXaxis()->CenterTitle();
			Eject1[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-60000,ejection1+index*ejeW+60000);
			Eject1[index]->GetYaxis()->SetTitle("Si_2 eje1 rising [ns]"); Eject1[index]->GetYaxis()->CenterTitle();
			Eject1[index]->SetMaximum(ejection1+index*ejeW+60000);
			Eject1[index]->SetMinimum(ejection1+index*ejeW-60000);
			Eject1[index]->SetMarkerStyle(31);
			Eject1[index]->SetMarkerSize(1.1);
			Eject1[index]->SetMarkerColor(kRed);
		}
		else{
			Eject1_o[index]->SetTitle("time2_relative:time1_relative");
			Eject1_o[index]->GetXaxis()->SetTitle("Si_1 eje1 dropping [ns]"); Eject1_o[index]->GetXaxis()->CenterTitle();
			Eject1_o[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-60000,ejection1+index*ejeW+60000);
			Eject1_o[index]->GetYaxis()->SetTitle("Si_2 eje1 dropping [ns]"); Eject1_o[index]->GetYaxis()->CenterTitle();
			Eject1_o[index]->SetMaximum(ejection1+index*ejeW+60000);
			Eject1_o[index]->SetMinimum(ejection1+index*ejeW-60000);

			Eject1[index]->SetTitle("time2_relative:time1_relative");
			Eject1[index]->GetXaxis()->SetTitle("Si_1 eje1 dropping [ns]"); Eject1[index]->GetXaxis()->CenterTitle();
			Eject1[index]->GetXaxis()->SetLimits(ejection1+index*ejeW-60000,ejection1+index*ejeW+60000);
			Eject1[index]->GetYaxis()->SetTitle("Si_2 eje1 dropping [ns]"); Eject1[index]->GetYaxis()->CenterTitle();
			Eject1[index]->SetMaximum(ejection1+index*ejeW+60000);
			Eject1[index]->SetMinimum(ejection1+index*ejeW-60000);
			Eject1[index]->SetMarkerStyle(31);
			Eject1[index]->SetMarkerSize(1.1);
			Eject1[index]->SetMarkerColor(kRed);
		}

	}// end of for format setting


	static vector <double> ejection0_si1_v[2];
	static vector <double> ejection0_si2_v[2];
	static vector <double> ejection1_si1_v[2];
	static vector <double> ejection1_si2_v[2];
	for(int iv=0;iv<2;iv++){
		ejection0_si1_v[iv].clear();
		ejection0_si2_v[iv].clear();
		ejection1_si1_v[iv].clear();
		ejection1_si2_v[iv].clear();
	}

	Long64_t si1_si2_avg=0;

	for(unsigned long index=0;index< Match_hit.size();index++){
			si1_si2_avg = (Match_hit[index].B1_time_relative + Match_hit[index].B2_time_relative)/2;
		for(int j=0;j<2;j++){
			if(si1_si2_avg > ejection0+j*ejeW-1000 && si1_si2_avg < ejection0+j*ejeW+6000){
				ejection0_si1_v[j].push_back( (double) Match_hit[index].B1_time_relative );
				ejection0_si2_v[j].push_back( (double) Match_hit[index].B2_time_relative );
				break;
			}
			if(si1_si2_avg > ejection1+j*ejeW-1000 && si1_si2_avg < ejection1+j*ejeW+6000){
				ejection1_si1_v[j].push_back( (double) Match_hit[index].B1_time_relative);
				ejection1_si2_v[j].push_back( (double) Match_hit[index].B2_time_relative);
				break;
			}
		}

	}

cout<<"eje0 up= "<<ejection0_si1_v[0].size()<<endl;
cout<<"eje0 down= "<<ejection0_si1_v[1].size()<<endl;
cout<<"eje1 up= "<<ejection1_si1_v[0].size()<<endl;
cout<<"eje1 down= "<<ejection1_si1_v[1].size()<<endl;

	for(int j=0;j<2;j++){
		for(unsigned long index=0;index< ejection0_si1_v[j].size();index++){
			Eject0[j]->SetPoint(index,ejection0_si1_v[j][index],ejection0_si2_v[j][index]);
		}
		for(unsigned long index=0;index< ejection1_si1_v[j].size();index++){
			Eject1[j]->SetPoint(index,ejection1_si1_v[j][index],ejection1_si2_v[j][index]);
		}
	}

	c_eject->cd(1)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(1)->Modified();
	Eject0_o[0]->Draw("AP");
	if(Eject0[0]->GetN())Eject0[0]->Draw("P");


	c_eject->cd(2)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(2)->Modified();
	Eject0_o[1]->Draw("AP");
	if(Eject0[1]->GetN())Eject0[1]->Draw("P");


	c_eject->cd(3)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(3)->Modified();
	Eject1_o[0]->Draw("AP");
	if(Eject1[0]->GetN())Eject1[0]->Draw("P");


	c_eject->cd(4)->SetMargin(0.120,0.052,0.1,0.1);
	c_eject->cd(4)->Modified();
	Eject1_o[1]->Draw("AP");
	if(Eject1[1]->GetN())Eject1[1]->Draw("P");


}




TCutG* BetaCut[15]={NULL};

//**************** make time relative vs. Ekev for beta in each channel **********************
// use it to make cut and to know the input time relative and Ekev is noise or not
// use MakeBetaCut_Or_FindSignal() bt default to start making new cuts
// inside any cut means noise
// set GetCutStatus = true    to know whether there is any cut available or not
// by default, start a guide to make new cuts in c_beta_coin_and_raw Canvas
// input the time_relative and Ekev of an event from specific channel to know it is noise or not ==> Is noise: return false; Signal: return true.

bool MakeBetaCut_Or_FindSignal(Long64_t _TimeRelative, double _EkeV, int WhichChannel, bool GetCutStatus, TFile* RootfilePtr, int Clear0_Recreate1_Add2){

	if(Clear0_Recreate1_Add2 >2 || Clear0_Recreate1_Add2<0){
		cout<<"Wrong action selection! Using 0 for clear cut; 1 for recreating cuts; 2 for adding cuts.  Abort!!"<<endl;
		return false;
	}

  	//*********** return the cut status ******************
	if(GetCutStatus){
		if(BetaCut[0]!=NULL) return true;
		else return false;
	}
	//**************************************************

	if(WhichChannel!=1 && WhichChannel!=2){
		cout<<"channel number limited to 1 or 2. Abort"<<endl;
		return false;
	}

	//********** ideentify an event is noise or not *************
	if(_TimeRelative != -1000 && _EkeV!=-1000){
		bool Isnoise=false;
		string energyName= Form("EkeV_%d",WhichChannel);
		for(int i=0;i<15;i++){
			if(BetaCut[i]!=NULL){
				string xname=BetaCut[i]->GetVarX();
				if(xname == energyName){
					Isnoise = BetaCut[i]->IsInside(_EkeV,_TimeRelative);
					if(Isnoise) break;
				}
				else continue;
			}
			else break;

		}

		return !Isnoise; // noise --> return false; signal --> return true;

	}

	
	if(RootfilePtr!= nullptr) RootfilePtr->cd();


//**************** clear cut exist and create new cut ********************
  static int cutindex=0;

    if(Clear0_Recreate1_Add2==2){// add cut
    	if(CutModifyCounter%2 ==1){// odd=> last active is made by loadcut; even -> last action is made by makecut
    		cutindex=0;
    		for(int i=0;i<15;i++){
    			if(BetaCut[i]!=NULL)cutindex++;
    		}
    	}
    } 
    else{
	    for(int i=0;i<15;i++){
			  	if(BetaCut[i]!=NULL){
			  	 	gROOT->GetListOfSpecials()->Remove((TObject*)BetaCut[i]);
			  	 	delete BetaCut[i];
			  	}
			  	BetaCut[i] = NULL;
		}

		TObject* obj;
		TIter iter(gROOT->GetListOfSpecials());
		while((obj=(TObject*)iter())){
		  		string cutname = obj->GetName();
		  		for(int i=0;i<15;i++){
		  			string cutname_candidate = Form("Beta_cut_%d",i+1);
		  			if(cutname == cutname_candidate){
		  				gROOT->GetListOfSpecials()->Remove(obj);
		  				break;
		  			}
		  		}
		}

		cutindex=0;
	}

  if(Clear0_Recreate1_Add2==0){
  		cout<<"\e[1;33m"<<"Beta cuts made by \"MakeBetaCut_Or_FindSignal()\" are all clear."<<endl;
  		if(CutModifyCounter%2==0)CutModifyCounter+=2;
  		else CutModifyCounter++;
  		return true;
  }

  if(c_beta_coin_and_raw == NULL){cout<<"canvas \"c_beta_coin_and_raw\" is not existed, ShowBeta_Beta_Coin() first!!"<<endl; return false;}
  c_beta_coin_and_raw->cd();

  cout << "You can draw a cut by clicking View->Toolbar->Graphical Cut" << endl;
  cout << "Double-click to close cut" << endl<<endl;


    while(1){
  		if(cutindex<15)cout<<"Add No. "<<cutindex+1<<" Beta noise cut? ('y' or 'n')"<<endl;
  		else cout<<"[1;33m"<<"Maxmimum number of cuts =15. Overwrite No. "<<cutindex%15+1<<" Beta noise cut? ('y' or 'n')"<<endl;

  		char yesorno='\0';
  		while(1){
  			cin>>yesorno;
  			if(yesorno=='y'|| yesorno=='n') break;
  			cout<<" 'y' or 'n' "<<endl;
  		}

  		if(yesorno=='n') break;
  		if(yesorno=='y'){
  			//remove the one will be overwrited first
  			if(cutindex>=15){
  			  	gROOT->GetListOfSpecials()->Remove((TObject*)BetaCut[cutindex%15]);
  			  	TObject* obj;
				TIter iter(gROOT->GetListOfSpecials());
				while((obj=(TObject*)iter())){
					string cutname = obj->GetName();
					string cutname_candidate = Form("Beta_cut_%d",cutindex%15+1);
					if(cutname == cutname_candidate){
						gROOT->GetListOfSpecials()->Remove(obj);
						break;
					}
				  		
				}
  	 			
  	 			delete BetaCut[cutindex%15];	
  	 			BetaCut[cutindex%15] = nullptr;
  			}

  			  cout<< "\e[1;33m"<<"Draw Time relative and Energy cut at c_beta_coin_and_raw (1 or 2 or 3 or 4)"<<"\e[0m"<<endl;
			  char cutname[15]={'\0'};
			  sprintf(cutname,"Beta_cut_%d",cutindex%15+1);
			  BetaCut[cutindex%15] = (TCutG*)c_beta_coin_and_raw->cd()->WaitPrimitive("CUTG");
			  BetaCut[cutindex%15]->SetName(cutname);
			  TObject* obj;
			  TIter iter(gROOT->GetListOfSpecials());
			  bool IsExistInList=false;
			  while((obj=(TObject*)iter())){
			  		string cutname = obj->GetName();
					string cutname_candidate = Form("Beta_cut_%d",cutindex%15+1);
					if(cutname == cutname_candidate){
						IsExistInList = true;
						break;
					}
				  		
			  }
			  if(!IsExistInList)gROOT->GetListOfSpecials()->Add((TObject *)BetaCut[cutindex%15]);			  
  		}

		cutindex++;

  	}
  	//%%% the end of loop for making new cuts 


  	//********* save cuts to file ******************************
  	TFile *cutf;
  
  	string cutf_name = MainFilePath + "cutfile_v2.root";
    cutf = new TFile(cutf_name.c_str(),"RECREATE");
    if(!cutf->IsOpen()){ cout<<"Can not recreate file.... Please check the path.... Abort!!!!"<<endl; return false;}
	  cutf->cd();
	  for(int i=0;i<15;i++){
	  		if(BetaCut[i]!=NULL)BetaCut[i]->Write(BetaCut[i]->GetName());
	  }
	  cutf->ls();
	  cutf->Close();
	  if(RootfilePtr != nullptr)RootfilePtr->cd();
	  //*********** cuts are save.

	if(CutModifyCounter%2==0)CutModifyCounter+=2; //CutModifyCounter even --> action by makecut; odd-> action by loadcut
  	else CutModifyCounter++;

	  return true;
}



bool LoadBetaCut(TFile* RootfilePtr= nullptr){

	string cutf_name = MainFilePath + "cutfile_v2.root";
  	TFile *fcut_in = new TFile(cutf_name.c_str(),"READ");
    if(!fcut_in->IsOpen()){ 
    	cout<<"Can not cut file is not existed!!!! Please check .... Abort!!!!"<<endl; 
    	if(RootfilePtr!=nullptr)RootfilePtr->cd();
    	return false;
    }

	if(RootfilePtr!=nullptr)RootfilePtr->cd();

	//%%%%%%%%%%%%%%%%% clear exist cut in memory
	for(int i=0;i<15;i++){
	  	if(BetaCut[i]!=NULL){
	  	 	gROOT->GetListOfSpecials()->Remove((TObject*)BetaCut[i]);
	  	 	delete BetaCut[i];
	  	}
	  	BetaCut[i] = NULL;
	}

  TObject* obj;
  TIter iter(gROOT->GetListOfSpecials());
  while((obj=(TObject*)iter())){
  		string cutname = obj->GetName();
  		for(int i=0;i<15;i++){
  			string cutname_candidate = Form("Beta_cut_%d",i+1);
  			if(cutname == cutname_candidate){
  				gROOT->GetListOfSpecials()->Remove(obj);
  				break;
  			}
  		}
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fcut_in->cd();
  int ncut=0;
  TIter nextobj(fcut_in->GetListOfKeys());
  while( (obj = (TObject *)nextobj()) ){
    if(fcut_in->Get(obj->GetName())->InheritsFrom("TCutG")){
    	if(ncut<15){
	      BetaCut[ncut] = (TCutG *)fcut_in->Get(obj->GetName());
	      ncut++;
  		}
    }
  }
/*
  for(int cutindex=0;cutindex<5;cutindex++){
  	if(BetaCut[cutindex]!=NULL){
  		gROOT->GetListOfSpecials()->Add( (TObject*)BetaCut[cutindex] );
  	}
  }*/

  fcut_in->ls();
  fcut_in->Close();
  if(RootfilePtr!=nullptr)RootfilePtr->cd();

  if(CutModifyCounter%2==1) CutModifyCounter+=2; //CutModifyCounter even --> action by makecut; odd-> action by loadcut
  else CutModifyCounter++;

  return true;
}


