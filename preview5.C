//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7
// main function input: string filename="",char origin_correct='o'/ 'c', double eje0=0, double eje1=0,double tof_ref=0
//canvas layout (1,2) 
//		     (3.4)
//get_para_draw(canvas_index) use as layout
//variables to change each time: m_ref ; bref; laps_ref;
//1.  histo_zoom_in_ref(int bins =0,double tof_ref=0,double halfwidth=0 double TAG=1); default; draw in cd(4)
//			defaule TAG 1 as reference; set TAG =0 can change to tag1 spectrum
//2.   fitquickly('r'); auto save peak centor to "tof_ref_cento"
//3.   histo_zoom_in_x(tag=0 ro 1); auto draw in cd(2)
//4.    fitquickly('x'); auto save peak centro to "tof_x_cento"
//5.    mass_calculator(Peak1orPeak2=1 or 2); 
//				auto read "tof_ref_cento" , "tof_x_cento" and laps_ref to calculate mass
//				1 for single peak curve case; 2 for double peak curve case
//6.   mass_calculator(double tof0,double tof_ref,int _laps_ref, double mass_ref, double b_ref); manually calculate any mass
//7. ErrCal_mass(double,double);  input: moq_x_err, _q_x 
//8. ErrCal_moq(double moq_x, double tofx, double tofx_err, double tofref, double tofref_err, double mass_ref, double mass_ref_err);
//9. PrintInfo()
//      Print information of reference: mass, q , laps, tof ; X Ion : tof
//10. ResetRef(double mass_ref=-1,double mass_ref_err=-1,int _qref=-1,int _laps_ref=-1,double tofref=-1,double tofref_err=-1);
//11. ResetRefMassTime(int Anum=0, const char* elename=" ", int _qref=-1,int _laps_ref=-1,double tofref=-1,double tofref_err=-1);
//               read in mass infor from AME2016.txt
//12. ResetX(int Peak1OrPeak2=-1 , double tofx=-1, double tofx_err=-1);   // -1 => don't set
//             Peak1 OR Peak2 control value set to tof_x_cento[0], tof_x_cento_err[0];
//13. MarkTof(double mass_xx , int _q_x, const char* IonName,int Tag=0,bool AddTrue=true);
//            tag==0 => show in tag0 spectrum; tag==1 => show in tag1 spectrum
//				addTrue ==> add marker or erase and refresh
//14. Rm_cal(TF1* inputline,int peaknum_index=1) input the line pointer and which peak (1 or 2) in case of double-peak curve to calculate Resoluting power
//15. SearchAME(int Anum,const char*element=" ",const char* filename="AME2016.txt")
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//10.	GenerateRejectList(tbeta_tof,mycut,&rejectlist);

#ifndef _PREVIEWER_
#define _PREVIEWER_


#include "TROOT.h"
#include "TRint.h"
#include"TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TAttAxis.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iomanip>
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TCutG.h"
#include "TObject.h"
#include "TCollection.h" // for TIter
#include "TFitResult.h"
#include<ctime>
#define _SETFUNC1_   // for gaus_exp fit function setting
#define _SETFUNC2_   // for exp_gaus_exp fit function setting

#ifdef _SETFUNC1_
#include "func1.h"   // define fit function
#endif

#ifdef _SETFUNC2_
#include "func2.h"   // define fit function
#endif

#include "Math/MinimizerOptions.h"// for ROOT::Math::MinimizerOptions...........
#include "Math/GenAlgoOptions.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/DataRange.h"
#include "Fit/UnBinData.h"
#include "Fit/Fitter.h"

#include "funcSample.h"
#include "funcN.h"

#include "chance.cc"

//#include "SearchAME.h"
//extern Double_t* SearchAME(const char* filename, int Anum,const char*element);

using namespace std;


#include "preview1.h"  // define SearchAME() and Resolving power calculator basing on different funcXX.h;  MUST be place after "using namespace std"

TCanvas *c1 = new TCanvas("c1","previewer",1200,800);

//&&&&&&&&&&&& for mass deviation plot drawing &&&&&&&&&&&&
TCanvas *c2 = NULL; // 
TGraphErrors *MassGraph=NULL;
TPad *padmass=NULL, *padbirge=NULL;
TLatex *text_delta_mass=NULL, *text_mass_mean=NULL, *text_birge=NULL;
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

bool canvas_initial = false;
bool initializer_stat = false;
TTree *intree;
TFile *fin = NULL;

string FilePath = "../";
string FileName = "file";


//&&&&&&&&&&&&& Reset beta tof detector efficiency measurement &&&&&&&&&&

bool Efficiency_histo_reset = false;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


//&&&&&&&&&&&&& default parameter for histo preview &&&&&&&&&&&&&
int sptrFW = 36000;  // full spectrum width (ns)
int binsfuspectrum = 10000;  // bins
int singleHalfwidth = 150;  // ns half spectrum width of ref histo
int bins_refS = 200;		// bins
int bins_nevt = 200;		// bins 2D histo X axis
string Usetag2D = "tag==1";   // 2D histogram selected by tag or not

int TagBit0 = 0; // using to reverse the eje time read from .lst to get proper preview
int TagBit1 = 1;
double EJE0=0;
double EJE1=0;
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



//&&&&&&&&&&& current active histo_X paras &&&&&&&&&&&&&&&&&&
// will be save in txt file
int active_tagX =0;
int active_histoX_Nbins = 0;
double active_histoX_RangeL = 0;
double active_histoX_RangeR =0;
string active_tree_name ="";

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



double numofevent=0;
long Maxsweeps=0;
double tof_ref_cento=0;        double  tof_ref_cento_err=0;       // contain tof input from labview calculator or fit tof result !!!!!!!!
double tof_x_cento[10]={0};     double  tof_x_cento_err[10]={0};   // contain tof fit result of x peak !!!!!!!!!!!!!
double tag1zoomL=0;	// zoom in spectrum left = tof_ref - 150
double tag1zoomR=0;   // zoom in spectrum rignt = tof_ref + 150


TH1D *h_xF = NULL;   // x ion full spectrum
TH1D *h_refF=NULL;   // ref ion full spectrum
TH1D *h_refS=NULL;   // ref ion single spectrum
TH2D *h_ref2D=NULL;  // ref ion 2D spectrum

/*TH1D *h_xF_tem = NULL;   // x ion full spectrum
TH1D *h_refF_tem=NULL;   // ref ion full spectrum
TH1D *h_refS_tem=NULL;   // ref ion single spectrum
TH2D *h_ref2D_tem=NULL;  // ref ion 2D spectrum
*/


TH1D *h_zoom_ref=NULL;
TH1D *h_zoom_x = NULL;
TH1D *h_zoom_x_shadow=NULL;

TH1D *histores = NULL; // residue histogram after original histo - value of fit curve

TH2D* h_beta_tof_2D=NULL; // 2D histogram of dk time Vs tof

//TLine *graphline;   // get manually draw line from canvas
TArrow *graphline;   // get manually draw line from canvas

//&&&&&&&&&&& Element TOF and ROI marker &&&&&&&&&&&&&&&
class TOFMarker{
	private:
			TLine *Pline;    // mark expected TOF of given mass
			TLatex *LabelXion;  // mark name of expected ion
	public:
			static int marker_amount;
			static bool marker_reset;
			double mass;
			int charge;
			int laps;
			double tof;		// the x position of the line
			double yHigh; // the high of the line
			int tag;
			char *name;

			void SetLine(double _x_tof,double _y_high){
				Pline->SetX1(_x_tof);
				Pline->SetX2(_x_tof);
				Pline->SetY1(0.);
				Pline->SetY2(_y_high);
				tof=_x_tof;
				yHigh = _y_high;
			}
			void SetName(string in_name){
				if(in_name.find("@")!=string::npos){in_name.erase(in_name.begin()+in_name.find("@"),in_name.end());}
				if(name !=NULL) delete[] name;
				name = new char[strlen(in_name.c_str())+1];
				strcpy(name,in_name.c_str());
			}
			void SetLabel(double _x,double _y, const char* _name){
				LabelXion->SetText(_x,_y,_name);
				LabelXion->SetTextAlign(21);
    			LabelXion->SetTextSize(0.07);
    			SetName(_name);
			}
			void SetMassAndCharge(double _mass,double _qx){
				mass = _mass;
				charge = _qx;
			}
			void SetColor(Int_t _kala){
				Pline->SetLineColor(_kala);
				LabelXion->SetTextColor(_kala);

			}
			void SetY(double _y){
				yHigh =_y;
				Pline->SetY2(yHigh);
				LabelXion->SetY(yHigh);
			}
			double GetX(){
				return Pline->GetX1();
			}
			void Copy(TOFMarker & inMarker){
				mass = inMarker.mass;
				charge = inMarker.charge;
				tof = inMarker.tof;
				yHigh = inMarker.yHigh;
				tag = inMarker.tag;
				laps = inMarker.laps;
				SetName(inMarker.name);
				SetLine(tof,yHigh);
				SetLabel(tof,yHigh,name);

			}
			void Clear(){Pline->Clear(); LabelXion->Clear();}
			void Draw(){Pline->Draw();LabelXion->Draw();}


			TOFMarker():name(NULL){
				Pline = new TLine;
				LabelXion = new TLatex;
				SetLine(0,0.1);
				tag = 0;
				laps = 0;
			}
			TOFMarker(TOFMarker & inMarker):name(NULL){
				TOFMarker();
				Copy(inMarker);
			}
			TOFMarker(double _x_tof,double _y_high, const char* _name="default"):name(NULL){
				Pline = new TLine;
				LabelXion = new TLatex;
				SetLine(_x_tof,_y_high);
				SetLabel(_x_tof,_y_high,_name);
				tag = 0;
				laps=0;
			}
			~TOFMarker(){
				delete Pline;
				delete LabelXion;
				if(name!= NULL)delete[] name;
			}

};

TOFMarker* marker_tof = new TOFMarker[40];
int TOFMarker::marker_amount = 0;
bool TOFMarker::marker_reset = true;


//*****************  ROI marker **************************
TLine * ROIL = new TLine[40];
TLine * ROIR = new TLine[40];
TLatex * ROI_label = new TLatex[40];
int ROI_tag[40];


double gXposition=0;  // to get tof of cursor by double click
double ROI_WIDTH =30;
int ROI_INDEX =1;    // begin from 1
bool ROI_initial=false;

void ROIClear(){
	for(int i=0;i<40;i++){
		ROIL[i].SetX1(1000);ROIL[i].SetX2(1000);
		ROIR[i].SetX1(1000);ROIR[i].SetX2(1000);
		ROI_label[i].SetX(1000);
		ROI_tag[i]=0;
	}

	ROI_INDEX =1;
	ROI_initial=false;
	c1->Modified();
	c1->Update();
}

void ROIadjY(TH1D* histo_in,double h_left,double h_right){
	string histo_name = histo_in->GetName();
	if(histo_name=="h_xF"){
		double LineHigh = histo_in->GetBinContent( histo_in->GetMaximumBin() ) * 1.01;
		if(LineHigh<1) LineHigh=1.01;  // 0 counts in histogram
		double TextHigh = LineHigh*0.5;
		if(ROI_INDEX<=40){
			for(int i=0;i<ROI_INDEX;i++){
				if(ROIL[i].GetX1()>10){ // not empty
					if(ROIL[i].GetX1()>=h_left && ROIR[i].GetX1()<=h_right){
						ROIL[i].SetY2(LineHigh);
						ROIR[i].SetY2(LineHigh);
						ROI_label[i].SetY(TextHigh);
					}
				}
			}
		}
		else{
			for(int i=0;i<40;i++){
				if(ROIL[i].GetX1()>10){ // not empty
					if(ROIL[i].GetX1()>=h_left && ROIR[i].GetX1()<=h_right){
						ROIL[i].SetY2(LineHigh);
						ROIR[i].SetY2(LineHigh);
						ROI_label[i].SetY(TextHigh);
					}
				}
			}
		}
	}
}



//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



//&&&&&&&&&& mass measurements Deviation plot &&&&&&&&&&
TLine massLowEdge;
TLine massHighEdge;
TLine ameLowEdge;
TLine ameHighEdge;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




//&&&&&&&&&&& predefined func for fitting &&&&&&&&&&&&&&
TF1 *gau;
TF1 *gaus_line;
TF1 *man1_func1; // 1 gaus_exp func
TF1 *man2_func1; // 1 gaus_exp func
TF1 *man3_func1; // 1 gaus_exp func
TF1 *ml_func1;   // 1 gaus_exp + linearbg func
TF1 *mexp_func1; // 1 gaus_exp + expbg func
TF1 *multi2_func1; // 2 gaus_exp func
TF1 *multi3_func1; // 3 gaus_exp func
TF1 *multi2_exp_func1; // 2 gaus_exp + 1 expbg func
TF1 *multi2_line_func1;// 2 gaus_exp + 1 linebg func

TF1 *man1_func2; // 1 gaus_exp func
TF1 *man2_func2; // 1 gaus_exp func
TF1 *man3_func2; // 1 gaus_exp func
TF1 *ml_func2;   // 1 gaus_exp + linearbg func
TF1 *mexp_func2; // 1 gaus_exp + expbg func
TF1 *multi2_func2; // 2 gaus_exp func
TF1 *multi3_func2; // 3 gaus_exp func
TF1 *multi2_exp_func2; // 2 gaus_exp + 1 expbg func
TF1 *multi2_line_func2;// 2 gaus_exp + 1 linebg func

TF1 *funcRb90=NULL;

TF1* tem_func=NULL; // get address of current fit function !!!!!!!!!!!!!

funcS* fs=NULL;
funcN* fext=NULL;

bool fixPara = false;  // set it to true if fix parameters to responding func manually

double tctail = 20; // default gaus_exp tail
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


//&&&&&&&&&&&& temper variable for fit func parameter &&&&&&&&&&&&
// get values from get_para_by_draw()
double tem_high=0;
double tem_cento=0;
double tem_sigma=0;
double tem_slop[2]={0};	// [0] for linear bg; [1] for exp bg
double tem_offset[2]={0};	// [0] for linear bg; [1] for exp bg
double fitrangeL=0;
double fitrangeR=0;
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&





//&&&&&&&&&&&&&&  variables for define and control general fit func &&&&&&&&&&&&&&&
int SetOneorTwo = 1;   // 0 = pure gaus;  1 = func1 (gauss + exp Right tail) ; 2 = func2 (gauss + exp Ritht and Left tail)
int NumOfPeaks = 1;    // at least 1 peak
int BackGroundCurve = 0;  // 0 = no background, 1 = linear; 2 = expotential
int MainPeakIndex =1;  // for multi-peak curve the most left one is 1; please use the strongest peak as main peak;
bool ParasLock= true; // for multi-peak curve case; all parameters of weaker peak are determined by corresponding parameters of main peak


int NumOfPars_TF_cal(int _SetOneorTwo , int _NumOfPeaks,int _whichbackground, bool _ParasLock){
	 int NumOfPars = 0;
	 int paras_offset = _SetOneorTwo+3;
	 if(_ParasLock){
	 	NumOfPars = 2 * (_NumOfPeaks-1)+ paras_offset;
	 }
	 else{
	 	NumOfPars = _NumOfPeaks * paras_offset;
	 }

	 if(_whichbackground !=0) NumOfPars+=2;

	 return NumOfPars;
}

// peakindex start from 0;   calculate the starting index of paras for the (NO. peakindex) peak in TF
int Paras_offset_cal(int peakindex,int _SetOneorTwo,int _MainPeakIndex,bool _ParasLock){
	if(_ParasLock){
		if(peakindex<=(_MainPeakIndex-1)) return peakindex*2; // before main peak
		else return (peakindex-1)*2 + (_SetOneorTwo+3);

	}
	else{// paras are independent for each peak
		return peakindex*(_SetOneorTwo+3);
	}

}

// form TF para to long setpar: parA_origin to parB_acceptor
void ParaDistributor(int _SetOneorTwo,int _NumOfPeaks,int _BackGroundCurve,int _MainPeakIndex, bool _ParasLock,double * parA_origin,double * parB_acceptor){
	// distribut paraA to paraB; number of paras in paraB >= paraA in multi-peak case
	
	int NumOfPars_B = (_SetOneorTwo+3)*_NumOfPeaks;
	if(_BackGroundCurve !=0)NumOfPars_B+=2;

	if(_ParasLock == false || _NumOfPeaks==1){// number of para_A =  number of para_B
		for(int i=0;i<NumOfPars_B;i++){parB_acceptor[i] = parA_origin[i];}
	}
	else{ // locked and multi-peak
		int para_offset = _SetOneorTwo+3;  // for para_B
		int NumOfPars_A = (_NumOfPeaks-1)*2 + (_SetOneorTwo+3);
		double main_para[4]={0}; // 0: sigma; 1: tailR; 2: tailL; 3: tailR_2;  for share between all peaks
		for(int i=0;i<(_SetOneorTwo+1);i++){
			main_para[i] = parA_origin[(_MainPeakIndex-1)*2+2+i];//[para_offset*(_MainPeakIndex-1) + 2+i];
		}

		int i=0;
		for( i=0;i<_NumOfPeaks;i++){
			
			if(i<=(_MainPeakIndex-1)){ // before main peak, only 2 paras for one peak in para_A
				parB_acceptor[para_offset*i] =	parA_origin[i*2];
				parB_acceptor[para_offset*i+1] = parA_origin[i*2+1];
				for(int j=0;j<(_SetOneorTwo+1);j++){
					parB_acceptor[para_offset*i+2+j] = main_para[j];
				}
			}
			else { // after the main peak
				parB_acceptor[para_offset*i] =	parA_origin[(i-1)*2+para_offset]; 
				parB_acceptor[para_offset*i+1] = parA_origin[(i-1)*2+para_offset+1];
				for(int j=0;j<(_SetOneorTwo+1);j++){
					parB_acceptor[para_offset*i+2+j] = main_para[j];
				}
			}
		}

		if(_BackGroundCurve !=0){
			parB_acceptor[para_offset*i] =	parA_origin[(i-1)*2+para_offset]; 
			parB_acceptor[para_offset*i+1] = parA_origin[(i-1)*2+para_offset+1];
		}
	}

}


void ParaContractor(int _SetOneorTwo,int _NumOfPeaks,int _BackGroundCurve,int _MainPeakIndex, bool _ParasLock,double * parA_origin,double * parB_acceptor){
		// number of paras in TF
	int NumOfPars_B = NumOfPars_TF_cal(_SetOneorTwo ,_NumOfPeaks,_BackGroundCurve,_ParasLock);
	if(_ParasLock == false || _NumOfPeaks==1){// number of para_A =  number of para_B
		for(int i=0;i<NumOfPars_B;i++){parB_acceptor[i] = parA_origin[i];}
	}
	else{//multi peak and paraslock
		int para_offset_A=0;
		int para_offset_B=0;
		int PeakIndex=0;
		for(PeakIndex=0;PeakIndex<_NumOfPeaks;PeakIndex++){
		 		para_offset_A = (_SetOneorTwo+3)*PeakIndex;
		 		para_offset_B = Paras_offset_cal(PeakIndex,_SetOneorTwo,_MainPeakIndex,_ParasLock);
		 	if(PeakIndex != (_MainPeakIndex-1)){ // not main peak case
		 		parB_acceptor[para_offset_B+0] = parA_origin[para_offset_A+0];
		 		parB_acceptor[para_offset_B+1] = parA_origin[para_offset_A+1];
		 	}
		 	else{// main peak case
		 		parB_acceptor[para_offset_B+0] = parA_origin[para_offset_A+0];
		 		parB_acceptor[para_offset_B+1] = parA_origin[para_offset_A+1];
		 		for(int i=0;i<(_SetOneorTwo+1);i++){
		 			parB_acceptor[para_offset_B+2+i] = parA_origin[para_offset_A+2+i]; // start from sigma => tctailR = tctailL
		 		}

		 	}
		}// after loop for peak number; PeakIndex point to Background

		if(_BackGroundCurve !=0){
				para_offset_A = (_SetOneorTwo+3)*PeakIndex;
		 		para_offset_B = Paras_offset_cal(PeakIndex,_SetOneorTwo,_MainPeakIndex,_ParasLock);
		 		parB_acceptor[para_offset_B+0] = parA_origin[para_offset_A+0];
		 		parB_acceptor[para_offset_B+1] = parA_origin[para_offset_A+1];
		}

	}

}


/*
double GeneralCurve(double *x, double *par){

		double YValue = 0;
		int PeakIndex=0;

		// calculate total number of parameter for all curve actually;
		int NumOfPars = (SetOneorTwo+3)*NumOfPeaks;
		if(BackGroundCurve !=0)NumOfPars+=2;

		//double * setpar = new double[NumOfPars];

		//distribute parameter input from outside; becasue sub peaks share same parameters with main peaks (sigma, taiL_1_Left, tail_1_R, tail_2_R)!!!


		if(SetOneorTwo ==0 || SetOneorTwo == 1){
			// multi peaks
			for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
				if(SetOneorTwo ==1) YValue = YValue + func1::fitfunc(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]);
				else YValue = YValue + func1::gas(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]);
			}
			// background
			if(BackGroundCurve == 0){ return YValue;}
			else if(BackGroundCurve == 1){ return YValue + func1::linearbg(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]); }
			else{return YValue + func1::expbg(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]);}

		}
		else if(SetOneorTwo ==2){

			for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
				YValue = YValue + func2::fitfunc(x,&par[5*PeakIndex]);
			}
			// background
			if(BackGroundCurve == 0){ return YValue;}
			else if(BackGroundCurve == 1){ return YValue + func2::linearbg(x,&par[5*PeakIndex]); }
			else{return YValue + func2::expbg(x,&par[5*PeakIndex]);}

		}
		else{
			cerr<<"fit function choice error: 0, 1 or 2"<<endl;
			return 0;
		}



}*/

double GeneralCurve(double *x, double *par){// new version

		double YValue = 0;
		int PeakIndex=0;

		// calculate total number of parameter for all curve actually;
		int NumOfPars = (SetOneorTwo+3)*NumOfPeaks;
		if(BackGroundCurve !=0)NumOfPars+=2;

		static double * setpar = NULL;
		if(setpar != NULL)delete[] setpar;
		setpar = new double[NumOfPars];

		//distribute parameter input from outside; becasue sub peaks share same parameters with main peaks (sigma, taiL_1_Left, tail_1_R, tail_2_R)!!!
		ParaDistributor(SetOneorTwo,NumOfPeaks,BackGroundCurve,MainPeakIndex,ParasLock,par,setpar);

		if(SetOneorTwo ==0 || SetOneorTwo == 1){
			// multi peaks
			for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
				if(SetOneorTwo ==1) YValue = YValue + func1::fitfunc(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]);
				else YValue = YValue + func1::gas(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]);
			}
			// background
			if(BackGroundCurve == 0){ return YValue;}
			else if(BackGroundCurve == 1){ return YValue + func1::linearbg(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]); }
			else{return YValue + func1::expbg(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]);}

		}
		else if(SetOneorTwo ==2){

			for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
				YValue = YValue + func2::fitfunc(x,&setpar[5*PeakIndex]);
			}
			// background
			if(BackGroundCurve == 0){ return YValue;}
			else if(BackGroundCurve == 1){ return YValue + func2::linearbg(x,&setpar[5*PeakIndex]); }
			else{return YValue + func2::expbg(x,&setpar[5*PeakIndex]);}

		}
		else{
			cerr<<"fit function choice error: 0, 1 or 2"<<endl;
			return 0;
		}


}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



//&&&&&&&&&&&&&& History container and Func spliter &&&&&&&&&&&&&&&7&&

//class DefGeneralFunc;
class SubTFPackage;

class DefGeneralFunc{
	friend class SubTFPackage;
	public:
		static int GloIndex;
		static double normalize_factor; // area of strongest peak!!!!!!!!!  for subline normalize
		char * name;
		char * Date;
		int SetOneorTwo;   // 0 = pure gaus;  1 = func1 (gauss + exp Right tail) ; 2 = func2 (gauss + exp Ritht and Left tail)
		int NumOfPeaks;    // at least 1 peak
		int BackGroundCurve;  // 0 = no background, 1 = linear; 2 = expotential
		int MainPeakIndex;
		bool ParasLock;
		double* TFParameters;
		TF1 * TFunc;
		
		int NumOfPars;  // how many paras for TF
		int TAG;
		double fitRangeL;
		double fitRangeR;
		int histo_nbins;
		double histo_L;
		double histo_R;
		double histo_bin_width;

		//&&&&&&&&&&&&&&   start of definition of class SubTGPackage  &&&&&&&&&&&&&&&&
		class SubTFPackage{
			friend class DefGeneralFunc;
			public:
				TF1 * subTF;
				TLatex *subTFLable[2]; // [0]: peak index; [1] intensity(peak area);
				TLine *subTFCenterLine;
				double subTFIntensity;
				int subTFIndex;
				Int_t kala[7];

				double Integral(double &LeftEdge,double &RightEdge,double _binwidth){
					
					double count_histo_sum= 0;
					for(double xposition=LeftEdge;xposition<RightEdge+_binwidth;xposition+=_binwidth){
						count_histo_sum += subTF->Eval(xposition);
					}
					return count_histo_sum;
				}

				SubTFPackage(const char* _name,int _SetOneorTwo, double *TFpars,double _binwidth){
						kala[0]=633;kala[1]=808;kala[2]=799;kala[3]=417;kala[4]=433;kala[5]=600;kala[6]=617;
						subTFLable[0]=NULL;subTFLable[1]=NULL;
						subTFCenterLine=NULL;
						subTFIndex = (DefGeneralFunc::GloIndex)++;

						double IntegralL=0;
						double IntegralR=0;
						switch (_SetOneorTwo){
							case 2:{
								subTF = new TF1(Form("%s_%02d",_name,subTFIndex),func2::fitfunc,0,25e6,5);
								subTF->SetParameters(TFpars);
								IntegralL=func2::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]-TFpars[2]*5,TFpars[1]-TFpars[2]); //5 times sigma to 1 sigma
								if(IntegralL== -1){ 
									IntegralL=func2::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]-1000,TFpars[1]-TFpars[2]); //-1000 ns to 1 sigma
								}
								IntegralR=func2::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]+TFpars[2],TFpars[1]+TFpars[2]*5); //5 times sigma to 1 sigma
								if(IntegralR== -1){ 
									IntegralR=func2::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]+TFpars[2],TFpars[1]+1000); //-1000 ns to 1 sigma
								}

								break;
							}
							default:{
								if(_SetOneorTwo==0){subTF = new TF1(Form("%s_%02d",_name,subTFIndex),"gaus",0,25e6);}
								else{subTF = new TF1(Form("%s_%02d",_name,subTFIndex),func1::fitfunc,0,25e6,4);}

								subTF->SetParameters(TFpars);
								IntegralL=func1::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]-TFpars[2]*5,TFpars[1]-TFpars[2]); //5 times sigma to 1 sigma
								if(IntegralL== -1){ 
									IntegralL=func1::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]-1000,TFpars[1]-TFpars[2]); //-1000 ns to 1 sigma
								}
								IntegralR=func1::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]+TFpars[2],TFpars[1]+TFpars[2]*5); //5 times sigma to 1 sigma
								if(IntegralR== -1){ 
									IntegralR=func1::find_x_to_y(subTF,TFpars[0]*0.01,TFpars[1]+TFpars[2],TFpars[1]+1000); //-1000 ns to 1 sigma
								}

								break;
							}
						}

						subTF->SetLineColor(kala[subTFIndex%7]);
						subTFIntensity = Integral(IntegralL,IntegralR,_binwidth);

						//&&&&&&&&&&&7 set lables &&&&&&&&&&&&
						for(int i=0;i<2;i++){  
							if(subTFLable[i]!=NULL){	subTFLable[i]->Clear(); delete subTFLable[i]; subTFLable[i]= new TLatex();}
							else{subTFLable[i]= new TLatex();}
							if(i==0){// index lable
								subTFLable[i]->SetText(TFpars[1],10,Form("%d",subTFIndex));
							}
							else{ // intensity lable
								subTFLable[i]->SetText(TFpars[1],TFpars[0]*1.1,Form("%.2f",subTFIntensity/DefGeneralFunc::normalize_factor*100));
							}
							subTFLable[i]->SetTextAlign(22);
							subTFLable[i]->SetTextSize(0.05);
						}

						if(subTFCenterLine !=NULL){subTFCenterLine->Clear();delete subTFCenterLine; subTFCenterLine=new TLine();}
						else{subTFCenterLine=new TLine();}
						subTFCenterLine->SetX1(TFpars[1]);subTFCenterLine->SetX2(TFpars[1]);
						subTFCenterLine->SetY1(0);subTFCenterLine->SetY2(TFpars[0]);

				}

				~SubTFPackage(){
					//cout<<"i am delete in SubTFPackage"<<endl;
					subTF->Clear();
					//subTF->Draw();
					subTF->Delete();
					delete subTFLable[0];
					delete subTFLable[1];
					delete subTFCenterLine;
					
				}

				void Draw(bool _showintens){
					subTF->Draw("same");
					subTFLable[0]->Draw();
					if(_showintens)subTFLable[1]->Draw();
					subTFCenterLine->Draw();
				}


		};
		//&&&&&&&&&&&&&&&&&&&&&& end of definition 	of class SubTFPackage &&&&&&&&&&&&&&



		SubTFPackage ** SubTF;
		//&&&&&&&&&&  control SubTFPackage from DefGeneralFunc &&&&&&&&&&&&&&&&&&&
		void ClearSubTF(){
			if(SubTF != NULL){
				for(int i=0;i<NumOfPeaks;i++){delete SubTF[i];} 
				delete[] SubTF; 
				SubTF= NULL;
			}
		}

		void BuildSubTF(){
			// clear SubTF memmmory first 
			ClearSubTF();
			SubTF = new SubTFPackage*[NumOfPeaks];

			int NumOfPars_all = (SetOneorTwo+3)*NumOfPeaks;
			if(BackGroundCurve !=0)NumOfPars_all+=2;

			double * setpar = new double[NumOfPars_all];
			::ParaDistributor(SetOneorTwo,NumOfPeaks,BackGroundCurve,MainPeakIndex,ParasLock,TFParameters,setpar);


			for(int i=0;i<NumOfPeaks;i++){
				int par_offset = (3+SetOneorTwo)*i;
				SubTF[i] = new SubTFPackage(name,SetOneorTwo,&setpar[par_offset],histo_bin_width);
			}

			delete[] setpar;

		}
		//!!!!!!!!!!  important: each time to draw; first clear SubTF,reset it to NULL; then generate object for SubTF again according to current 
		// curve structure;
		void DrawSubTF(bool showintens=true){
			BuildSubTF();
			if(SubTF !=NULL){// already built successfully; allow to draw
				for(int i=0;i<NumOfPeaks;i++){
					SubTF[i]->Draw(showintens);
				}
			}
			else{
				cout<<"\e[1;33m"<<"SubTF container is empty; No function can be draw"<<"\e[0m"<<endl;
			}
		}
		//&&&&&&&&&&&&&&&&&&& end of SubTFPackage controller &&&&&&&&&&&&&&&&&&&&&&&



		//&&&&&&&&&&&&&&&&&& member function of DefGeneralFunc &&&&&&&&&&&&&&&&&&&
		/*Double_t GeneralCurve(Double_t *x, Double_t *par){ // function to define TFunc;

				double YValue = 0;
				int PeakIndex=0;

				if(SetOneorTwo ==0 || SetOneorTwo == 1){
					// multi peaks
					for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
						if(SetOneorTwo ==1) YValue = YValue + func1::fitfunc(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]);
						else YValue = YValue + func1::gas(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]);
					}
					// background
					if(BackGroundCurve == 0){ return YValue;}
					else if(BackGroundCurve == 1){ return YValue + func1::linearbg(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]); }
					else{return YValue + func1::expbg(x,&par[ (4+(SetOneorTwo-1)) *PeakIndex]);}

				}
				else if(SetOneorTwo ==2){

					for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
						YValue = YValue + func2::fitfunc(x,&par[5*PeakIndex]);
					}
					// background
					if(BackGroundCurve == 0){ return YValue;}
					else if(BackGroundCurve == 1){ return YValue + func2::linearbg(x,&par[5*PeakIndex]); }
					else{return YValue + func2::expbg(x,&par[5*PeakIndex]);}

				}
				else{
					cerr<<"fit function choice error: 0, 1 or 2"<<endl;
					return 0;
				}
		}*/

		double GeneralCurve(double *x, double *par){// new version

				double YValue = 0;
				int PeakIndex=0;

				// calculate total number of parameter for all curve actually;
				int NumOfPars_all = (SetOneorTwo+3)*NumOfPeaks;
				if(BackGroundCurve !=0)NumOfPars_all+=2;

				static double * setpar = NULL;
				if(setpar!=NULL) delete[] setpar; 
				setpar = new double[NumOfPars_all];

				//distribute parameter input from outside; becasue sub peaks share same parameters with main peaks (sigma, taiL_1_Left, tail_1_R, tail_2_R)!!!
				::ParaDistributor(SetOneorTwo,NumOfPeaks,BackGroundCurve,MainPeakIndex,ParasLock,par,setpar);

				if(SetOneorTwo ==0 || SetOneorTwo == 1){
					// multi peaks
					for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
						if(SetOneorTwo ==1) YValue = YValue + func1::fitfunc(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]);
						else YValue = YValue + func1::gas(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]);
					}
					// background
					if(BackGroundCurve == 0){ return YValue;}
					else if(BackGroundCurve == 1){ return YValue + func1::linearbg(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]); }
					else{return YValue + func1::expbg(x,&setpar[ (4+(SetOneorTwo-1)) *PeakIndex]);}

				}
				else if(SetOneorTwo ==2){

					for(PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){ // each loop point to the start of next group of paras if available
						YValue = YValue + func2::fitfunc(x,&setpar[5*PeakIndex]);
					}
					// background
					if(BackGroundCurve == 0){ return YValue;}
					else if(BackGroundCurve == 1){ return YValue + func2::linearbg(x,&setpar[5*PeakIndex]); }
					else{return YValue + func2::expbg(x,&setpar[5*PeakIndex]);}

				}
				else{
					cerr<<"fit function choice error: 0, 1 or 2"<<endl;
					return 0;
				}


		}
	
		/*void cal_NumOfPars(){
				NumOfPars = (SetOneorTwo+3)*NumOfPeaks;
				if(BackGroundCurve !=0)NumOfPars+=2;
		}*/


		void SetTFInfo(int _SetOneorTwo=-1,int _NumOfPeaks=-1,int _BackGroundCurve=-1,bool _ParasLock=true,int _MainPeakIndex=-1,double* _parameters=NULL,double _fitRangeL=-1,double _fitRangeR=-1,int _histo_nbins=-1,double _histo_L=-1,double _histo_R=-1){
			bool SetStatus = false;  // to check whether the basic structure of the TF curve changed or not

			//!!!!!!!!!!!!!!!!!! 	IMPORTANT:   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//have to clear SubTF memory first, in case any change of curve structure !!!!!!!!!!!!!!!!!!!
			// because the clear precess rely on the current curve structure setting; if structure change happens before SubTF clear; the following 
			//clear porcess (every time to Draw subTF will invoke the build precess which will invoke clear process to clear old history first) will do 
			//based on the new curve structure, not the old one to clear old SubTF memory improperly; this will lead to protential
			// memory leak !!!!!!!!!!!!!!!!!!
			
			ClearSubTF();

			if(_SetOneorTwo !=-1){SetOneorTwo = _SetOneorTwo;SetStatus = true;}
			if(_NumOfPeaks!=-1){NumOfPeaks = _NumOfPeaks;SetStatus = true;}
			if(_BackGroundCurve!=-1){BackGroundCurve = _BackGroundCurve;SetStatus = true;}
			if(_ParasLock != ParasLock){ParasLock = _ParasLock; SetStatus = true;}
			if(_MainPeakIndex != -1){MainPeakIndex = _MainPeakIndex; SetStatus = true;}
			if(_parameters!=NULL){
				//cal_NumOfPars(); // renew number of parameters
				NumOfPars = ::NumOfPars_TF_cal(SetOneorTwo , NumOfPeaks,BackGroundCurve, ParasLock);
				if(TFParameters != NULL){delete[] TFParameters; TFParameters=NULL;}
				TFParameters = new double[NumOfPars];
				for(int i=0;i<NumOfPars;i++){TFParameters[i]=_parameters[i];}
			}
			if(_fitRangeL!=-1){fitRangeL = _fitRangeL;}
			if(_fitRangeR!=-1){fitRangeR = _fitRangeR;}
			if(_histo_nbins!=-1){histo_nbins = _histo_nbins;}
			if(_histo_L!=-1){histo_L = _histo_L;}
			if(_histo_R!=-1){histo_R = _histo_R;}

			if(SetStatus){ // structure of TF curve already changed
				if(TFunc != NULL)delete TFunc; 
				TFunc = new TF1(name,this,&DefGeneralFunc::GeneralCurve,0,25e6,NumOfPars,"1DefGeneralFunc","1GeneralCurve");
			}

			if(TFParameters != NULL && TFunc !=NULL){ // have availble parameter and TF curve exist; get paras to define TF curve;
				TFunc->SetParameters(TFParameters);
			}

			histo_bin_width = (histo_R - histo_L)/ histo_nbins;

			
		}


		void SetName(const char* _name){
			if(name != NULL)delete[]name;
			name = new char[strlen(_name)+1];
			strcpy(name,_name);
		}

		void SetDate(const char* _Date){
			if(Date != NULL)delete[]Date;
			Date = new char[strlen(_Date)+1];
			strcpy(Date,_Date);
		}

		DefGeneralFunc():SetOneorTwo(1),NumOfPeaks(1),BackGroundCurve(0),ParasLock(true),MainPeakIndex(1),name(NULL),TFParameters(NULL),SubTF(NULL),TFunc(NULL),Date(NULL){
		}
		DefGeneralFunc(const char* _name,int _tag=0,int _SetOneorTwo=1,int _NumOfPeaks=1,int _BackGroundCurve=0,int _MainPeakIndex=1, bool _ParasLock = true,double* _parameters=NULL,double _fitRangeL=0,double _fitRangeR=0,int _histo_nbins=0,double _histo_L=0,double _histo_R=0,const char* _Date=NULL){
			TFParameters = NULL;
			SubTF = NULL;
			TFunc = NULL;
			name = NULL;
			TAG = _tag;
			Date = NULL;
			SetName(_name);SetDate(_Date);
			SetTFInfo(_SetOneorTwo,_NumOfPeaks,_BackGroundCurve,_ParasLock,_MainPeakIndex,_parameters,_fitRangeL,_fitRangeR,_histo_nbins,_histo_L,_histo_R);
		}
		DefGeneralFunc(TF1* inputTF){
			if(inputTF != NULL){
				cout<<"Define GeneralFunc with current setting as tem_func!!!!"<<endl;
				TFParameters = NULL;
				SubTF = NULL;
				TFunc = NULL;
				TAG = active_tagX;
				name = NULL;Date = NULL;
				time_t now = time(0);
				SetDate(ctime(&now));
				SetName("current_fit_function");
				SetTFInfo(::SetOneorTwo,::NumOfPeaks,::BackGroundCurve,::ParasLock,::MainPeakIndex,inputTF->GetParameters(),::fitrangeL,::fitrangeR,::active_histoX_Nbins,::active_histoX_RangeL,::active_histoX_RangeR);
			}
			else{
				cout<<"Define empty GeneralFunc !!!!"<<endl;
				DefGeneralFunc("empty");
			}

		}

		~DefGeneralFunc(){
			if(TFParameters != NULL) delete TFParameters;
			if(TFunc != NULL)delete TFunc;
			if(name != NULL)delete[] name;
			if(Date != NULL)delete[] Date;
			ClearSubTF();
		}

		bool operator()(TF1* inputTF=NULL){
			if(inputTF != NULL){
				cout<<"Define GeneralFunc with current setting as tem_func!!!!"<<endl;
				//TFParameters = NULL;
				TAG = active_tagX;
				time_t now = time(0);
				SetDate(ctime(&now));
				SetName("current_fit_function");
				SetTFInfo(::SetOneorTwo,::NumOfPeaks,::BackGroundCurve,::ParasLock,::MainPeakIndex,inputTF->GetParameters(),::fitrangeL,::fitrangeR,::active_histoX_Nbins,::active_histoX_RangeL,::active_histoX_RangeR);
				return true;
			}
			else return false;
		}

		void SetTFPars2Container(){TFunc->GetParameters(TFParameters);} // copy paras in TFunc to TFParameters; while ->GetParameters() return address of container in TFunc
		void SetContainerPars2TF(){TFunc->SetParameters(TFParameters);}

		TF1* GetTFunc(){return TFunc;}

		void ScalePeaksHigh(double factor){ 
			// for easily match the height of peaks in current histo with different bins width
			// only change parameters in TFunc not those in container
			double* tem_paras = new double[NumOfPars];
			for(int i=0;i<NumOfPars;i++){tem_paras[i] = TFParameters[i];}

			int peak_para_offset=0;
			for(int PeakIndex=0;PeakIndex<NumOfPeaks;peak_para_offset = ::Paras_offset_cal(++PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock)){
				tem_paras[peak_para_offset]*=factor;
			} 

			TFunc->SetParameters(tem_paras);
			delete[] tem_paras;
		}

		void Draw(bool ShowSubTF=false){
			TFunc->Draw("same");
			if(ShowSubTF==true) DrawSubTF();

		}

};


int DefGeneralFunc::GloIndex=0;
double DefGeneralFunc::normalize_factor=100;
DefGeneralFunc *FHistory[30];
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&





//&&&&&&&&&&& variable for mass calculation &&&&&&&&&&&&
double t0 = 130;   double err_t0 = 10;// ns
const double m_ele = 548.58 ; const double err_me = 0.003;  // unit micro_u
double m_ref = 132905451.961;  double err_ref = 0.009;    // unit micro_u  133Cs
//const double m_ref = 84911789.738; const double err_ref = 0.006;    // unit micro_u  85Rb
//const double m_ref = 86909180.531; const double err_x =  0.006;	// unit micro_u 87Rb
//const double m_ref = 38963706.487; const double err_x =  0.005;  // unit micro_u 39K
int q_ref = 1;     // charge
double bref=33256.31;
double bref_err = 0; // negligible by default

int laps_ref = 600;
int laps_x = laps_ref;
int q_x = 1;
//const double KEV90 = 0.93149400380; const double err_KEV90 = 0.0004*1e-6; // 1 u = 931494.0038 KEV90, mass unit in micro u in program
const double KEV90 = 0.93149410242; const double err_KEV90 = 0.00028*1e-6;
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//&&&&&&&&&& variable for component mass &&&&&&&&&&&&&&&&&
double component_mass[20]={0};
double component_mass_err[20]={0};
char component_name[20][20];
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



//&&&&&&&&&&&   Mass recorder And Time recorder &&&&&&&&&&

double mass_recorded[20]={0};
double mass_recorded_err[20]={0};
const char * mass_recorded_message[20];

double time_recorded[20]={0};
double time_recorded_err[20]={0};
const char * time_recorded_message[20];

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


//&&&&&&&&&& ion information container &&&&&&&&&&&&&&&&&&
struct Ion{
	double mass;   // atomic mass, including all electron
	double mass_err;
	double tof;
	double tof_err;
	int charge;
	int counts;
	const char* name;
	string _name_tem;
	string _name;

	string MakeName_s(string formula){
		string iso[100];  // iostop name 12C
	    int weight[100];  // number of this iostop in a molecule
	    int compon = 0;

	    istringstream is(formula);
	    char buffer[100];
	    char *p;
	    while(is>>buffer){
		    p = strtok(buffer,";");
		    iso[compon] = p;
		    p = strtok(NULL,";");
		    if(p) weight[compon] = atoi(p);
		    else  weight[compon] = 1;
		    compon++;
  	    }

  	    int Anum[100]={0};
  	    for(int icompon=0; icompon<compon; icompon++){
    		string sA, El;
    		for(unsigned long int i=0; i<iso[icompon].size(); i++){
	     		 char c = iso[icompon][i];
	      		if(isdigit(c)) sA+=c;
	      		else El+=c;
	    	}

	    	Anum[icompon] = atoi(sA.c_str());
	    	iso[icompon] = El;
	    }

	    string str_return="";
	    for(int icompon=0; icompon<compon; icompon++){
	    	if(icompon<compon-1) str_return+= "^{" + to_string(Anum[icompon]) +"}"+ iso[icompon]+ to_string(weight[icompon])+":";
	    	else str_return+= "^{" + to_string(Anum[icompon]) +"}"+ iso[icompon]+ to_string(weight[icompon]);
	    }

	    return str_return;

	}

	Ion(double _mass=0, double _mass_err=0,int _charge=1,double _tof=0, double _tof_err=0, const char* _inname = NULL,int _counts=0){
		mass = _mass;
		mass_err = _mass_err;
		tof = _tof;
		tof_err= _tof_err;
		charge = _charge;
		name = _inname;
		counts = _counts;
	}

	Ion(string formula,int _charge=1,double _tof=0, double _tof_err=0,int _counts=0){
		double* getmass = SearchMolMass(formula);
		//formula.erase(remove(formula.begin(),formula.end(),';'), formula.end());
		_name_tem="";
		istringstream is(formula);
	  	char buffer[100];
	  	string sub_name;
	  	while(is>>buffer){
		    sub_name = buffer;
		    _name_tem = _name_tem + sub_name +":" ;
		}

		_name_tem[_name_tem.size()-1]='*';
		_name_tem.erase(remove(_name_tem.begin(),_name_tem.end(),'*'),_name_tem.end());

		_name = MakeName_s(formula);
		mass = getmass[0];
		mass_err = getmass[1];
		tof = _tof;
		tof_err = _tof_err;
		charge = _charge;
		name = _name_tem.c_str();
		counts = _counts;
	}

};
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



void LoadNewBetaFile(string); //for loading beta file
void autopreview(string PATH="../rootfiles/",char origin_correct='o',double _intime_ref=0);
void autoloadBetaFile();
void useTOFveto(double);
void clearance();
void initializer_func1();
void initializer_func2();
void initializer();
void get_para_by_draw(int canvas=1);

//*********************
void para_modifier(int FullWidth_ns = 36000,int binsFW = 10000,int refHalfWidth_ns = 150,int bins_refSingle=200,int bins_event=100,int whichtag2D=1);
//   width of full spectrum , bins of full spectrum , half width of ref spectrum from centro , bins of ref zoom in histo , bins of event in 2D
//*********************

void fitquickly_func1(char whichhisto='x',Option_t *Fitoption = "ME", bool FitAgain=false, double _RangeL=-1, double _RangeR=-1);
void fitquickly_func2(char whichhisto='x',Option_t *Fitoption = "ME", bool FitAgain=false, double _RangeL=-1, double _RangeR=-1);
void fitquickly(char whichhisto='x',int func1ORfunc2=1,Option_t *Fitoption = "ME", bool FitAgain=false, double _RangeL=-1, double _RangeR=-1,int N_loop=1);
//void fitgeneralcurve(char whichhisto='x',int func1ORfunc2=1, int NumPeaks2Fit = 1, int whichbackground =0,Option_t *Fitoption = "ME", bool FitAgain=false, double _RangeL=-1, double _RangeR=-1);
void fitgeneralcurve(char whichhisto='x',int func1ORfunc2=1, int NumPeaks2Fit=1, int whichbackground=0, int Index_MainPeak=1,bool _ParasLock=true,Option_t *Fitoption= "ME", bool FitAgain=false, double _RangeL=-1, double _RangeR=-1);
double ErrCal_mass(double,double);
double ErrCal_rho(double _rho, double tofx, double tofx_err, double tofref, double tofref_err);
double ErrCal_moq(double _rho, double _rho_err, double mass_ref, double mass_ref_err, int _q_ref);
//double ErrCal_moq(double moq_x, double tofx, double tofx_err, double tofref, double tofref_err, double mass_ref, double mass_ref_err);
double * mass_calculator(double tof0,double tof0_err,int _laps_x, int _q_x, double tof_ref,double tof_ref_err,int _laps_ref, double mass_ref, double mass_ref_err,
										double b_ref, double b_ref_err,int _q_ref, bool verbal=true);
void mass_calculator(int Peak1OrPeak2=1);
void mass_calculator(double tof0,double tof0_err,double tof_ref,double tof_ref_err,int _laps_ref, double mass_ref, double mass_ref_err,double b_ref,double b_ref_err,int _q_ref=1);
void ResetRef(double mass_ref=-1,double mass_ref_err=-1,int _qref=-1,int _laps_ref=-1,double b_ref=-1,double tofref=-1,double tofref_err=-1); // default -1 means don't set
void ResetRefMassTime(int Anum=0, const char* elename=" ", int _qref=-1,int _laps_ref=-1,double b_ref=-1,double tofref=-1,double tofref_err=-1); // default -1 means don't set
void ResetRefMassTime(string formula="12C;1 16O;2", int _qref=-1,int _laps_ref=-1,double b_ref=-1,double tofref=-1,double tofref_err=-1); // default -1 means don't set
void ResetX(int Peak1OrPeak2=-1 , double tofx=-1, double tofx_err=-1);  // default -1 means don't set

void MarkTof(double mass_xx , int _q_x, const char* IonName,int _laps_x=0,int Tag=0,bool renew_current=false,bool renew_all=false);
void MarkTof(int Anum,const char* EleName, int _q_x, int _laps_x=0, int Tag=0,bool renew_current=false,bool renew_all=false);
void MarkTof(string formula, int _q_x, int _laps_x=0,int Tag=0,bool renew_current=false,bool renew_all=false);

void SetROI(double TCento, double ROI_width=30, int ROI_index = 1, double TOffset=0,bool showcount=false,int Tag=0);
void SetMassROI(double _massxx,int _q_x=1, int _laps_x=-1, int ROI_index=1, double TOffset=0, double ROI_width=30, bool showcount=false,int Tag=0);
void SetEleROI(int Anum,const char *element, int _q_x =1, int _laps_x=-1, int ROI_index=1, double TOffset=0, double ROI_width=30,bool showcount=false,int Tag=0);
void SetMoleculeROI(string formula, int _q_x=1, int _laps_x=-1, int ROI_index=1, double TOffset=0, double ROI_width=30, bool showcount=false,int Tag=0);
void CreateComponent(int Compon_index,int Anum,const char* elename,int N_atom, bool addTrue=true);
void PrintComponentList();
void PrintMassTimeList();
void SaveMassResult(int result_index=1, const char * Note="", double mass_result=-1, double mass_result_err=-1);
void SaveTimeResult(int result_index=1, const char * Note="", double time_result=-1, double time_result_err=-1);
double* ExtractMassDeviate(double mass_A, double mass_err_A,double mass_B, double mass_err_B,bool saveMass=false,int result_index=1, const char * Note="");
void savefile_current_fitparas(const char* Note, bool NewIndex=false);
void savefile_history_fitparas(DefGeneralFunc * history, const char *Note=NULL,bool NewIndex=false);
void ReverseEjeTime();


bool dflag=true;
double veto_test=200;
double Tof_veto=0; //[ns] veto reject hits of the same sweeep in the window
int Tof_Channel2Read=1;
int TAG_hzoom=0;
double* sweeps_glo_array_hzoomx=nullptr;
double* Tof_array_hzoomx=nullptr;
long long Tof_array_hzoomx_entries=0;

// get TOF 
long long GetTof_wo_Twin(TTree* _intree = nullptr, string branches_name="time", string selection_str = "tag==1 && channel==1", 
						double** outV1=nullptr, double** outV2=nullptr){
		TTree* treecp =_intree;  //= _intree->CopyTree(selection_str.c_str());
		long long NEntriesGet =0;
		NEntriesGet = treecp->Draw(branches_name.c_str(),selection_str.c_str(),"goff");
		if(NEntriesGet == 0){
			cout<<"\e[1;31m"<<"Find no responding data, Please check!!"<<"\e[0m"<<endl;
			return NEntriesGet;
		} 

		treecp->SetEstimate(NEntriesGet);
		treecp->Draw(branches_name.c_str(),selection_str.c_str(),"goff");
		*outV1 = treecp->GetV1();
		*outV2 = treecp->GetV2();
		return NEntriesGet;


}


long long GetToF_w_Twin(TTree* _intree = nullptr,TObject* histo_obj = nullptr, string branches_name="time:sweeps_global", string GeneralSelection="" ,int _tag=1, int _channel=1,
						double veto_T_win=0,double** outV1=nullptr,double** outV2=nullptr, bool ShowCondiction=true){ 
	// define 1D or 2D histogdram histo_in ==> histo_obj = (TObject*) histo_in, hisot_in will be filled;
	//Automatically identify the input histogram object is 1D or 2D, if 1D, filling the histogram with tof; if 2D, fill the histogram with tof and global sweeps: tof => Y; global sweeps => X
	//define pointer double* gett1; double* gett2 ==> outV1=& gett1 ,  outV2=&gett2 to get output variables given in branches_name; access element by gett1[i] from outside;
	// V1 is tof ; V2 is sweeps_global as default
	//return number is the valid number of elements in array outV1 and outV2
	//when veto_T_win >0; branches_name must be "time:sweeps_global" or "timec:sweeps_global"
	// if tag!=1 and tag!=0, ==> do not use tag to select data

	if(veto_T_win<0){cout<<"veto_T_win must be >0; abort!!!"<<endl; return 0;}

	static double * outV1_inner=nullptr;
	static double * outV2_inner=nullptr;

	if(outV1_inner != nullptr){ delete[] outV1_inner; outV1_inner = nullptr; }
	if(outV2_inner != nullptr){ delete[] outV2_inner; outV2_inner = nullptr; }

	string sele_condiction="";
	if(_tag!=0 && _tag!=1) sele_condiction = Form("channel == %d",_channel);
	else sele_condiction = Form("tag == %d && channel == %d",_tag,_channel);

	if(GeneralSelection != "")sele_condiction = "(" + sele_condiction + ") && (" + GeneralSelection + ")";

	if(veto_T_win>0){
		if(branches_name.find("time") != 0 || branches_name.find("sweeps_global") == string::npos){// branch of time or timec is requried;
			cout<<"\e[1;31m"<<"When veto time window is provided; branches_name must be \"time:sweeps_global\" or \"timec:sweeps_global\""<<"\e[0m"<<endl;
			return 0; 
		}
	}

	if(ShowCondiction) cout<<"\e[1;33m"<<sele_condiction.c_str()<<"\e[0m"<<endl;

	double* getV1=nullptr, * getV2=nullptr;
	long long NEntries_raw = GetTof_wo_Twin(_intree,branches_name,sele_condiction,&getV1,&getV2);
//cout<<"raw entries = "<<NEntries_raw<<endl;

	if(NEntries_raw<1){
		cout<<"\e[1;31m"<<"find no tof data, Please check!!"<<"\e[0m"<<endl;
		return 0;
	}

	long long NEntries_filted=0;
	double Last_tof = 0;
	double Last_sweeps_glo=0;
	outV1_inner = new double[NEntries_raw];
	outV2_inner = new double[NEntries_raw];

	for(long long index=0;index<NEntries_raw;index++){
		if(veto_T_win>0){
			if(getV2[index] == Last_sweeps_glo){
				if(getV1[index]-Last_tof < veto_T_win) continue;
			}
		}
		
		outV1_inner[NEntries_filted] = getV1[index];
		outV2_inner[NEntries_filted++] = getV2[index];
		Last_sweeps_glo = getV2[index];
		Last_tof = getV1[index];
	}


	if(outV1 != nullptr) *outV1 = outV1_inner;
	if(outV2 != nullptr) *outV2 = outV2_inner;

	TH1D* hist_1D=nullptr; 
	TH2D* hist_2D=nullptr;
	if(histo_obj!=nullptr){
		if(histo_obj->InheritsFrom("TH2")) hist_2D = (TH2D*) histo_obj; // TH2 is also derived from TH1!!!!!
		else if(histo_obj->InheritsFrom("TH1")) hist_1D = (TH1D*) histo_obj;
		for(long long index=0;index<NEntries_filted;index++){
			double tof_tem = outV1_inner[index];
			double sweeps_glo_tem = outV2_inner[index];
			if(hist_1D != nullptr) hist_1D->Fill(tof_tem);
			if(hist_2D != nullptr) hist_2D->Fill(sweeps_glo_tem,tof_tem);  //(sweeps_global, time);
		}
	}

	return NEntries_filted;

}




void preview(string PATH="../rootfiles/",string filename="",char origin_correct='o', double eje0=0, double eje1=0,double tof_ref=0){

	FilePath  = PATH;
	FileName = filename;

	if((origin_correct!='o'&& origin_correct!='c') || eje0<0 || eje1<0 || tof_ref<0){
		cout<<"para error: origin_correct == 'o' or 'c'; eje0>0; eje1>0; tof_ref>0"<<endl;
		return;
	}

   if(eje0==0 && eje1==0 ){ // read para for spectrum from .lst file

   		double paracontainer[6]; // eje0; eje1; timer,massr, bref, qx

   		int goodornot = ReadParaLst(PATH,filename,paracontainer); // -1 for fail; 1 for good

   		if(goodornot==-1){cerr<<"\e[1;33m"<< "Parameters are incomplete in .LST . Parameters loading fail!!! Input paras manually to open file."<<"\e[0m"<<endl; return;}

   		eje0=paracontainer[TagBit0]*1e3; // time unit is micro s in list => change to ns in program
   		eje1=paracontainer[TagBit1]*1e3; // can be reversed by ReverseEjeTime() to get correct preview;
   		EJE0=eje0;   // actually this is just the DAQ delay; for real ejection moment = DAQ delay -2 us  ==> EJE0 - 2000
   		EJE1=eje1;

   		if(tof_ref == 0) tof_ref=paracontainer[2]*1e3;  // if not 0 , use inputed number;  read from .lst by default

   		double massfroml_func1st = paracontainer[3];

   		bref=paracontainer[4]*1000;   // br output in us =>change to [ns]

   		laps_ref = (int)(EJE1 / bref);
   		laps_x = laps_ref;

   		q_x = paracontainer[5];

   		double* massMatchAME = SearchEleByMass(massfroml_func1st*1e6);

   		if(massMatchAME[0] != -1){ // successfully open AME2016.txt
   				m_ref=massMatchAME[0];
   				err_ref =massMatchAME[1];
   		}
   		else cerr<<"\e[1;33m"<< "Fail to set mass of Ref automatically!!! Manually set please."<<"\e[0m"<<endl;

   		
   		DefGeneralFunc::GloIndex=0;
   }
   else{
   		double eje_tem[]={eje0,eje1};
	   	EJE0=eje_tem[TagBit0];
	   	EJE1=eje_tem[TagBit1];
	   	eje0 = EJE0;
	   	eje1 = EJE1;
   }


	clearance();

	if(initializer_stat==false){ initializer(); initializer_stat =true;}

	  tof_ref_cento = tof_ref;

      string inputfile;

      if(origin_correct == 'o') inputfile = PATH + filename + ".root";//"rootfiles/" + 
      else if(origin_correct == 'c') inputfile = PATH + filename + "_dcorrect.root";
      //else {cerr<<"\e[1;33m"<<"False: choice error 'o' for original ; 'c' for drift corrrection"<<"\e[0m"<<endl; return;}

      if(gSystem->AccessPathName(inputfile.c_str())) {cout<<"Path or file name not exist!!"<<endl;return;}

      fin = new TFile(inputfile.c_str(),"READ");
	if(!fin->IsOpen()){ cout<<"fail to open file: "<<filename<<".root"<<endl; return;}
	if(origin_correct=='o'){
      	intree = (TTree *)fin->Get("tree");
      	if(intree == NULL){
      		intree = (TTree *)fin->Get("tree0");
      		if(intree == NULL){
      			cout<<"\e[1;31m"<<"Find no tree!!!!!!! Error!!!!!!!"<<"\e[0m"<<endl;
      		}
      		else{
      			cout<<"\e[1;31m"<<"Find no original tree!!! Change to dcorrect tree!!!!"<<"\e[0m"<<endl;
      		}
      	}else{
      		cout<<"get oiginal tree"<<endl;
      	}
      	
	}
	else{
      	intree = (TTree *)fin->Get("tree0");
      	if(intree == NULL){
      			cout<<"\e[1;31m"<<"Find no tree!!!!!!! Error!!!!!!!"<<"\e[0m"<<endl;
      	}else{
      		cout<<"get dcorrect tree"<<endl;
      	}
      	
	}

	numofevent =(double) intree->GetEntriesFast();

	double tag0L = eje0 - 300;
	double tag0H = tag0L + sptrFW;
	double tag1L= eje1 - 300;
	double tag1H = tag1L + sptrFW;

	tag1zoomL = tof_ref - singleHalfwidth;
	tag1zoomR = tof_ref + singleHalfwidth;


	double Maxsweeps_tem=0;
	intree->Draw("sweeps_global","","goff");
	double* getsweeps = intree->GetV1();

	for(long isweeps=0;isweeps<intree->GetSelectedRows() % intree->GetEstimate();isweeps++){
		if(getsweeps[isweeps]>Maxsweeps_tem) Maxsweeps_tem =getsweeps[isweeps];
	}

	Maxsweeps=(long)Maxsweeps_tem;

	cout<<"File statistics: "<<endl;
	cout<<"numofevent = "<<numofevent<<endl;
	cout<<"Maxsweeps = "<<Maxsweeps<<"\t"<<"Time:  "<<Maxsweeps/40<<" [s]"<<endl;

	h_xF = new TH1D("h_xF","tag of ion X ",binsfuspectrum,tag0L,tag0H);
	h_refF = new TH1D("h_refF","tag of ion Reference ",binsfuspectrum,tag1L,tag1H);
	h_refS = new TH1D("h_refS","tag of ion Reference zoom in ",bins_refS,tag1zoomL,tag1zoomR);
	h_ref2D = new TH2D("h_ref2D","tag of ion Reference zoom in 2D",bins_nevt,0,Maxsweeps_tem,bins_refS,tag1zoomL,tag1zoomR);


if(!dflag)useTOFveto(veto_test);

	if(canvas_initial==false){ c1->Divide(2,2); canvas_initial = true;/*ctest->Divide(2,2)*/}

	if(origin_correct=='o'){
		double* getVar1=nullptr;
		double* getVar2=nullptr;
		long long getEntries=0;

		GetToF_w_Twin(intree,h_xF,"time:sweeps_global","",0,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);

		getEntries=GetToF_w_Twin(intree,h_refF,"time:sweeps_global","",1,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);

		for(long long index=0;index<getEntries;index++){
			h_refS->Fill(getVar1[index]);
			//h_ref2D->Fill(getVar2[index],getVar1[index]);
		}

		if(Usetag2D.find("0") != string::npos) GetToF_w_Twin(intree,h_ref2D,"time:sweeps_global","",0,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);
		else if(Usetag2D.find("1") != string::npos) GetToF_w_Twin(intree,h_ref2D,"time:sweeps_global","",1,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);
		else GetToF_w_Twin(intree,h_ref2D,"time:sweeps_global","",-1,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);


		c1->cd(1)->SetEditable(kTRUE);
			if(dflag){//intree->Draw("time>>h_xF","tag==0");
				h_xF->Draw();
			}
			else h_xF->Draw();
		c1->cd(3);
			if(dflag){//intree->Draw("time>>h_refF","tag==1");
				h_refF->Draw();
			}
			else h_refF->Draw();
		c1->cd(2)->SetEditable(kTRUE);
			if(dflag){//intree->Draw("time>>h_refS","tag==1");
				h_refS->Draw();
			}
			else h_refS->Draw();
		c1->cd(4);
		gStyle->SetPalette(kBird);
			if(dflag){//intree->Draw("time:sweeps_global>>h_ref2D",Usetag2D.c_str(),"colz"); //:nevt
				h_ref2D->Draw("colz");
			}
			else h_ref2D->Draw("colz");

/*	GetToF_w_Twin(intree,h_xF_tem,"time:sweeps_global","",0,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
	GetToF_w_Twin(intree,h_refF_tem,"time:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
	GetToF_w_Twin(intree,h_refS_tem,"time:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
	GetToF_w_Twin(intree,h_ref2D_tem,"time:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);

	
	ctest->cd(1);
	h_xF_tem->Draw();
	ctest->cd(3);
	h_refF_tem->Draw();
	ctest->cd(2);
	h_refS_tem->Draw();
	ctest->cd(4);
	h_ref2D_tem->Draw("colz");*/



	}else{
		double* getVar1=nullptr;
		double* getVar2=nullptr;
		long long getEntries=0;

		GetToF_w_Twin(intree,h_xF,"timec:sweeps_global","",0,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);

		getEntries=GetToF_w_Twin(intree,h_refF,"timec:sweeps_global","",1,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);

		for(long long index=0;index<getEntries;index++){
			h_refS->Fill(getVar1[index]);
			//h_ref2D->Fill(getVar2[index],getVar1[index]);
		}

		if(Usetag2D.find("0") != string::npos) GetToF_w_Twin(intree,h_ref2D,"timec:sweeps_global","",0,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);
		else if(Usetag2D.find("1") != string::npos) GetToF_w_Twin(intree,h_ref2D,"timec:sweeps_global","",1,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);
		else GetToF_w_Twin(intree,h_ref2D,"timec:sweeps_global","",-1,Tof_Channel2Read,Tof_veto,&getVar1,&getVar2,false);


		c1->cd(1)->SetEditable(kTRUE);
		if(dflag){//intree->Draw("timec>>h_xF","tag==0");
			h_xF->Draw();
		}
	else h_xF->Draw();
		c1->cd(3);
		if(dflag){//intree->Draw("timec>>h_refF","tag==1");
			h_refF->Draw();
		}
	else h_refF->Draw();
		c1->cd(2)->SetEditable(kTRUE);
		if(dflag){//intree->Draw("timec>>h_refS","tag==1");
			h_refS->Draw();
		}
	else h_refS->Draw();
		c1->cd(4);
		if(dflag){//intree->Draw("timec:sweeps_global>>h_ref2D",Usetag2D.c_str(),"colz"); //:nevt
			h_ref2D->Draw("colz");
		}
	else h_ref2D->Draw("colz");

/*	GetToF_w_Twin(intree,h_xF_tem,"timec:sweeps_global","",0,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
	GetToF_w_Twin(intree,h_refF_tem,"timec:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
	GetToF_w_Twin(intree,h_refS_tem,"timec:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
	GetToF_w_Twin(intree,h_ref2D_tem,"timec:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
	
	ctest->cd(1);
	h_xF_tem->Draw();
	ctest->cd(3);
	h_refF_tem->Draw();
	ctest->cd(2);
	h_refS_tem->Draw();
	ctest->cd(4);
	h_ref2D_tem->Draw("colz");*/


	}
	

	if(h_xF!=NULL){
		h_xF->GetXaxis()->SetTitle("tof[ns]");
		h_xF->GetXaxis()->SetTitleSize(0.05);
		h_xF->GetXaxis()->CenterTitle();
		h_xF->GetYaxis()->SetTitle(Form("counts/%.1f [ns]",(tag0H-tag0L)/binsfuspectrum));
		h_xF->GetYaxis()->SetTitleSize(0.05);
		h_xF->GetYaxis()->CenterTitle();
	}	
	if(h_refF!=NULL){
		h_refF->GetXaxis()->SetTitle("tof[ns]");
		h_refF->GetXaxis()->SetTitleSize(0.05);
		h_refF->GetXaxis()->CenterTitle();
		h_refF->GetYaxis()->SetTitle(Form("counts/%.1f [ns]",(tag1H-tag1L)/binsfuspectrum));
		h_refF->GetYaxis()->SetTitleSize(0.05);
		h_refF->GetYaxis()->CenterTitle();
	}	
	if(h_refS!=NULL){
		h_refS->GetXaxis()->SetTitle("tof[ns]");
		h_refS->GetXaxis()->SetTitleSize(0.05);
		h_refS->GetXaxis()->CenterTitle();
		h_refS->GetYaxis()->SetTitle(Form("counts/%.1f [ns]",(tag1zoomR-tag1zoomL)/bins_refS));
		h_refS->GetYaxis()->SetTitleSize(0.05);
		h_refS->GetYaxis()->CenterTitle();
	}
	if(h_ref2D!=NULL){
		h_ref2D->GetXaxis()->SetTitle("sweeps");
		h_ref2D->GetXaxis()->SetTitleSize(0.05);
		h_ref2D->GetXaxis()->CenterTitle();
		h_ref2D->GetYaxis()->SetTitle("tof[ns]");
		h_ref2D->GetYaxis()->SetTitleSize(0.05);
		h_ref2D->GetYaxis()->SetTitleOffset(1.1);
		h_ref2D->GetYaxis()->CenterTitle();
	}
	ROI_INDEX =1;
	ROI_initial = false; //reset ROI staus for interactive mode

	c1->Modified();
	c1->Update();
	c1->cd(1)->SetEditable(kFALSE);
	c1->cd(2)->SetEditable(kFALSE);

	autoloadBetaFile();

}


void clearance(){
	//c1->cd(1)->SetEditable(kTRUE);
	//c1->cd(2)->SetEditable(kTRUE);
	if(h_xF!=NULL){h_xF->Delete(); h_xF = NULL;}//h_xF->Clear();
	if(h_refF!=NULL){h_refF->Delete(); h_refF = NULL;}//h_refF->Clear();
	if(h_refS!=NULL){h_refS->Delete(); h_refS=NULL;}//h_refS->Clear();
	if(h_ref2D!=NULL){h_ref2D->Delete();h_ref2D = NULL;}//h_ref2D->Clear();

/*
	if(h_xF_tem!=NULL){h_xF_tem->Delete(); h_xF_tem = NULL;}//h_xF->Clear();
	if(h_refF_tem!=NULL){h_refF_tem->Delete(); h_refF_tem = NULL;}//h_refF->Clear();
	if(h_refS_tem!=NULL){h_refS_tem->Delete(); h_refS_tem=NULL;}//h_refS->Clear();
	if(h_ref2D_tem!=NULL){h_ref2D_tem->Delete();h_ref2D_tem = NULL;}//h_ref2D->Clear();  */


    // h_zoom_ref and h_zoom_x have links to tree file. when tree is clean, .root file is closed. These two histogram  are pointed to unknown space;
    // run "delete" in histo_zoom_in_ref/x() will cause crash.
    // Should be also "DELETED" when tree and file are "delete"!!!!!!!!!!!!!!
	if(h_zoom_ref!=NULL){ h_zoom_ref->Delete(); h_zoom_ref = NULL;} //h_zoom_ref->Clear();
	if(h_zoom_x!=NULL){h_zoom_x->Delete(); h_zoom_x =NULL;}//h_zoom_x->Clear();

	if(h_beta_tof_2D!=NULL){delete h_beta_tof_2D; h_beta_tof_2D=NULL;}

	if(marker_tof != NULL){
		for(int index=0;index<40;index++){
			marker_tof[index].Clear();
			
		}

		delete[] marker_tof;


		marker_tof = new TOFMarker[40];
		TOFMarker::marker_amount = 0;
		TOFMarker::marker_reset=true;

	}

	if(fs!=NULL){delete fs; fs=NULL;}
	fs = new funcS; // sample function
	
	if(fext!=NULL){delete fext; fext=NULL;}
	fext = new funcN;

	if(intree!=NULL){
		intree=NULL;
	}
	if(fin!=NULL){
		fin->Close();
		fin->Delete();
		fin = NULL;
	}

	Efficiency_histo_reset = false; // used in Exec.h


	//c1->cd(1)->Modified();
	//c2->cd(2)->Modified();
	//c1->Update();
}

void autopreview(string PATH,char origin_correct,double _intime_ref){//string PATH="../rootfiles/",char origin_correct='o'

	string commandline = PATH;
	//string::size_type pos=0;
	//pos=commandline.find("Run");
	//commandline.erase(commandline.begin(),commandline.begin()+pos+3);
	//pos = commandline.find("/");
	//commandline.erase(commandline.begin()+pos,commandline.end());
	
	commandline = ".! ls -t "+ PATH + "../LST/*@*.lst > " + PATH + "../LST/tof_filename.log";
	gROOT->ProcessLine(commandline.c_str());


	string tempathfile = PATH + "../LST/tof_filename.log";
	FILE* namelog =NULL;
	namelog = fopen(tempathfile.c_str(),"r");
	if(namelog==NULL){cout<<"\e[1;31m"<<" tof_filename.log not exist!!! "<<"\e[0m"<<endl; return;}

	char firstname[1000]{'\0'};
	string getname;
	bool fileexist = false;
	vector<string> getnamelist;

	while(!feof(namelog)){
		fscanf(namelog,"%s",firstname);
		getname = firstname;
		string::size_type pos=0;
		while((pos = getname.find("/")) != string::npos){ getname.erase(getname.begin(),getname.begin()+pos+1);}   // 1 is necessary!!!!

		if(getname.find("@") != string::npos && getname.find(".lst") != string::npos){
		 	fileexist = true; 
			getname.erase(getname.begin()+getname.find(".lst"),getname.end());
		 	getnamelist.push_back(getname);
		}

	}

	if(fileexist){
		if(getnamelist.size()>1){
			for(int i=0;i<(int)getnamelist.size();i++){
				for(int j=i+1;j<(int)getnamelist.size();j++){
					if(getnamelist[j]<getnamelist[i]){
						getname = getnamelist[i];
						getnamelist[i] = getnamelist[j];
						getnamelist[j] = getname;
					}
				}
			}
		}

		getname = getnamelist[0];

		if(_intime_ref>0) preview(PATH,getname,origin_correct,0,0,_intime_ref);
		else preview(PATH,getname,origin_correct);
	}
	else cout<<"\e[1;31m"<<"Found no TOF file. Manually run preview()"<<"\e[0m"<<endl;
	fclose(namelog);


}


void autoloadBetaFile(){
	//&&&&&&&&&&&&&&& load beta file automatically &&&&&&&&&&&&&&&
	
	string commandline = ".! ls -t "+ FilePath + "../LST/L*.lst > " + FilePath + "../LST/beta_filename.log";
	gROOT->ProcessLine(commandline.c_str());
	string tempathfile = FilePath + "../LST/beta_filename.log";
	FILE* namelog =NULL;
	namelog = fopen(tempathfile.c_str(),"r");
	if(namelog==NULL){cout<<"\e[1;31m"<<" beta_filename.log not exist!!! "<<"\e[0m"<<endl; return;}
	char firstname[1000]{'\0'};
	string firstBetaFile;
	vector<string> BfilenameList;
while(!feof(namelog)){
	fscanf(namelog,"%s",firstname);
	if(!feof(namelog)){
		string getname = firstname;
		string::size_type pos=0;
		while((pos = getname.find("/")) != string::npos){ getname.erase(getname.begin(),getname.begin()+pos+1);}   // 1 is necessary!!!!
		getname.erase( getname.begin()+getname.find(".lst")+4, getname.end());
	 	if(getname.find("L") != string::npos && getname.find(".lst") != string::npos) BfilenameList.push_back(getname);
	}

}
	//fscanf(namelog,"%s",firstname);
	//string getname = firstname;
	//string::size_type pos=0;
	//while((pos = getname.find("/")) != string::npos){ getname.erase(getname.begin(),getname.begin()+pos+1);}   // 1 is necessary!!!!
	
	if(BfilenameList.size() ==0) cout<<"\e[1;31m"<<"Found no beta file. Manually run LoadNewBetaFile(\"......lst\")"<<"\e[0m"<<endl;
	else{
		int * filenameIndex_in = new int[BfilenameList.size()];
		int * filenameIndex_out = new int[BfilenameList.size()];

		for(int ifile=0;ifile<(int)BfilenameList.size();ifile++){
			string filename_i_tem = BfilenameList[ifile];
			filename_i_tem.erase(filename_i_tem.begin()+filename_i_tem.find(".lst"),filename_i_tem.end()); // from "L_123_001.lst" to "L_123_001"

			while(filename_i_tem.find("_")!=string::npos){
				filename_i_tem.erase(filename_i_tem.begin(),filename_i_tem.begin()+filename_i_tem.find("_")+1); // from "L_123_001" to "001"
			}

			filenameIndex_in[ifile]= atoi(filename_i_tem.c_str()); // from "001" to interger 1
			//PATH = PATH+"LST/";
		}

		TMath::Sort((int)BfilenameList.size(),filenameIndex_in,filenameIndex_out,kFALSE);

		cout<<"Get a list of beta files:"<<endl;
		for(int ifile=0;ifile<(int)BfilenameList.size();ifile++){
			if(ifile%3==0) cout<<BfilenameList[filenameIndex_out[ifile]];
			else if(ifile%3==1) cout<<"\t"<<BfilenameList[filenameIndex_out[ifile]];
			else cout<<"\t"<<BfilenameList[filenameIndex_out[ifile]]<<endl;
		}

		cout<<endl;
		cout<<"Set the file: "<<BfilenameList[filenameIndex_out[0]]<<" to LoadNewBetaFile()"<<endl;
		LoadNewBetaFile(BfilenameList[filenameIndex_out[0]]);
	}

	fclose(namelog);

}



void useTOFveto(double Veto_window =200){
		vector< vector<Double_t> > eve_mr_time; // store tof  vector of each event in tree
		vector< vector <Long_t> > eve_mr_sweeps_global; // store sweeps_global vector of each event in tree
		vector< vector <Int_t> > eve_mr_tag;

		vector <Double_t> *mr_time = new vector <Double_t>();
		vector <Long_t> * mr_sweeps_global = new vector <Long_t>;
		vector <Int_t> *mr_tag = new vector <Int_t>();

		intree->ResetBranchAddresses();
		intree->SetBranchAddress("sweeps_global",&mr_sweeps_global);
		intree->SetBranchAddress("tag",&mr_tag);
		string treename = intree->GetName();
		if(treename == "tree") intree->SetBranchAddress("time",&mr_time);
		else{ intree->SetBranchAddress("timec",&mr_time);}

		Long64_t N_entries = intree->GetEntriesFast();

		for(Long64_t event_i=0;event_i<N_entries;event_i++){

			mr_time->clear();
			mr_sweeps_global->clear();
			mr_tag->clear();
			intree->GetEntry(event_i);

			eve_mr_time.push_back(*mr_time);
			eve_mr_sweeps_global.push_back(*mr_sweeps_global);
			eve_mr_tag.push_back(*mr_tag);

		}


		// combine data with same global sweeps
		for(Long64_t event_i=0;event_i<N_entries;event_i++){

			long glo_sweeps_i = eve_mr_sweeps_global[event_i][0];

			for(Long64_t event_i_before=0; event_i_before< event_i; event_i_before++){

				if(glo_sweeps_i == eve_mr_sweeps_global[event_i_before][0] ){ // find data with same global sweeps ==> combine

					for(unsigned long ihit=0;ihit<(unsigned long)eve_mr_sweeps_global[event_i].size();ihit++){
						eve_mr_sweeps_global[event_i_before].push_back(eve_mr_sweeps_global[event_i][ihit]);
						eve_mr_time[event_i_before].push_back(eve_mr_time[event_i][ihit]);
						eve_mr_tag[event_i_before].push_back(eve_mr_tag[event_i][ihit]);
					}
					eve_mr_sweeps_global[event_i][0] =-99; // marker for beening combined
					break; //==> next current event;
				} 

			}

		}

		// make sequence
		for(Long64_t event_i=0;event_i<N_entries;event_i++){
			if(eve_mr_sweeps_global[event_i][0]>0){// first global sweep ==1; reject glo sweeps ==-99
				vector<double> tof_to_fill;
				vector<int> tof_tag_to_fill;

				if(eve_mr_sweeps_global[event_i].size() >1 ){ // more than one hit 

						const unsigned long N_hits = eve_mr_time[event_i].size();
						double time_tem[N_hits];
						vector<int> tag_cp;

						for(unsigned long ihit=0;ihit<N_hits;ihit++){ // make an original copy of tof in one sweep
								time_tem[ihit]=eve_mr_time[event_i][ihit];
								tag_cp.push_back(eve_mr_tag[event_i][ihit]);
						}

						long long tof_sequence_list[N_hits];
						TMath::Sort((long long)N_hits,time_tem,tof_sequence_list,kFALSE);
						sort(eve_mr_time[event_i].begin(),eve_mr_time[event_i].end());

				
						for(unsigned long ihit=0;ihit<N_hits;ihit++){
							eve_mr_tag[event_i][ihit] = tag_cp[ tof_sequence_list[ihit] ];
						}

				}
				

				for(unsigned long ihit=0;ihit<(unsigned long)eve_mr_time[event_i].size();ihit++){
					if(ihit == 0){
						tof_to_fill.push_back(eve_mr_time[event_i][0]);
						tof_tag_to_fill.push_back(eve_mr_tag[event_i][0]);
					}
					else{
						if( (eve_mr_time[event_i][ihit]-tof_to_fill[tof_to_fill.size()-1]) > Veto_window){
							tof_to_fill.push_back(eve_mr_time[event_i][ihit]);
							tof_tag_to_fill.push_back(eve_mr_tag[event_i][ihit]);
						}
					}
				}
				
				for(unsigned long ihit=0;ihit<(unsigned long)tof_to_fill.size();ihit++){
					if(tof_tag_to_fill[ihit] ==0 ){
						h_xF->Fill(tof_to_fill[ihit]);
					}
					else{
						h_refF->Fill(tof_to_fill[ihit]);
						h_refS->Fill(tof_to_fill[ihit]);
						h_ref2D->Fill(eve_mr_sweeps_global[event_i][0],tof_to_fill[ihit]);
					}
				}
					
			}
		}


	intree->ResetBranchAddresses();

}



void para_modifier(int FullWidth_ns, int binsFW, int refHalfWidth_ns, int bins_refSingle, int bins_event, int whichtag2D){
	// width of full spectrum , bins of full spectrum , half width of ref spectrum from centro , bins of ref zoom in histo , bins of event in 2D
	sptrFW = FullWidth_ns;  // full spectrum width
	binsfuspectrum = binsFW;
	singleHalfwidth = refHalfWidth_ns;
	bins_refS = bins_refSingle;
	bins_nevt = bins_event;

	if(whichtag2D==0){ Usetag2D = "tag==0";}
	else if(whichtag2D==1){ Usetag2D = "tag==1";}
	else{ Usetag2D = "";}

}

void DisplaySetting(){
	printf("FullWidth_ns = %d\nbinsFW= %d\nrefHalfWidth_ns= %d\nbins_refSingle= %d\nbins_event= %d\nwhichtag2D= %s\n",sptrFW,binsfuspectrum,singleHalfwidth,bins_refS,bins_nevt,Usetag2D.c_str());
}




void histo_zoom_in_ref(int bins =0,double tof_ref=0,double halfwidth=0, double TAG=1){
		if(bins!=0){bins_refS = bins;}
		if(tof_ref!=0){tof_ref_cento = tof_ref;}
		if(halfwidth!=0){singleHalfwidth = halfwidth;}

		if(h_zoom_ref!=NULL){ h_zoom_ref->Delete(); h_zoom_ref = NULL;}//h_zoom_ref->Clear();


		h_zoom_ref = new TH1D("h_zoom_ref","Reference zoom in",bins_refS,tof_ref_cento-singleHalfwidth,tof_ref_cento+singleHalfwidth);
		c1->cd(4);

		cout<<"binwidth: "<<singleHalfwidth*2.0/bins_refS<<" [ns]"<<endl;

		h_zoom_ref->GetXaxis()->SetTitle("tof[ns]");
		h_zoom_ref->GetXaxis()->SetTitleSize(0.05);
		h_zoom_ref->GetXaxis()->CenterTitle();
		h_zoom_ref->GetYaxis()->SetTitle(Form("counts/%.1f [ns]",(double)singleHalfwidth*2/bins_refS));
		h_zoom_ref->GetYaxis()->SetTitleSize(0.05);
		h_zoom_ref->GetYaxis()->CenterTitle();

		string treename = intree->GetName();
		if(treename =="tree0"){
			if(TAG==1){// intree->Draw("timec>>h_zoom_ref","tag==1");
					GetToF_w_Twin(intree,h_zoom_ref,"timec:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
			}
			else if(TAG==0){ //intree->Draw("timec>>h_zoom_ref","tag==0");
					GetToF_w_Twin(intree,h_zoom_ref,"timec:sweeps_global","",0,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
			}
			else{
					GetToF_w_Twin(intree,h_zoom_ref,"timec:sweeps_global","",-1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
			}

		}
		else{
			if(TAG==1){ //intree->Draw("time>>h_zoom_ref","tag==1");
					GetToF_w_Twin(intree,h_zoom_ref,"time:sweeps_global","",1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
			}
			else if(TAG==0){ //intree->Draw("time>>h_zoom_ref","tag==0");
					GetToF_w_Twin(intree,h_zoom_ref,"time:sweeps_global","",0,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
			}
			else{
					GetToF_w_Twin(intree,h_zoom_ref,"time:sweeps_global","",-1,Tof_Channel2Read,Tof_veto,nullptr,nullptr,false);
			}
		}

		c1->cd(4);

		h_zoom_ref->Draw();

}



void histo_zoom_in_x(int tag=0,int bins=0,double histoL=0,double histoR=0){

	if(bins==0 || histoL==0 || histoR==0){
		cout<<"histo X parameter invalid"<<endl; 
		cout<<"get parameter from graph?? [y/n]"<<endl;
		char yesorno;
		do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
		if(yesorno=='y'){
			cout<<"pending..."<<endl;
			get_para_by_draw();
			bins = 100;
			histoL = fitrangeL;
			histoR = fitrangeR;
			printf("used: histo_zoom_in_x(%d,%d,%.1f,%.1f)\n",tag,bins,histoL,histoR);
		}
		else return;
	}

	if(h_zoom_x!=NULL){delete h_zoom_x; h_zoom_x=NULL;}//h_zoom_x->Clear();
	//if(h_zoom_x_shadow!=NULL){delete h_zoom_x_shadow;h_zoom_x_shadow=NULL;}

	h_zoom_x = new TH1D("h_zoom_x",Form("X ion tag%d zoom in",tag),bins,histoL,histoR);
	//h_zoom_x_shadow = new TH1D("h_zoom_x_shadow","X ion(tag0 or tag1) zoom in",bins,histoL,histoR);


	cout<<"binwidth: "<<(histoR-histoL)/bins<<" [ns]"<<endl;

		h_zoom_x->GetXaxis()->SetTitle("tof[ns]");
		h_zoom_x->GetXaxis()->SetTitleSize(0.05);
		h_zoom_x->GetXaxis()->CenterTitle();
		h_zoom_x->GetYaxis()->SetTitle(Form("counts/%.1f [ns]",(histoR-histoL)/bins));
		h_zoom_x->GetYaxis()->SetTitleSize(0.05);
		h_zoom_x->GetYaxis()->CenterTitle();

	// record current histoX setting, save to txt file when need;
	active_tagX = tag;
	active_histoX_Nbins = bins;
	active_histoX_RangeL = histoL;
	active_histoX_RangeR = histoR;
	active_tree_name = "intree";

	c1->cd(2)->SetEditable(kTRUE);
	string treename = intree->GetName();
	if(tag==1){
		if(treename=="tree0"){
			//intree->Draw("timec>>h_zoom_x","tag==1");
			Tof_array_hzoomx_entries = GetToF_w_Twin(intree,h_zoom_x,"timec:sweeps_global","",1,Tof_Channel2Read,Tof_veto,&Tof_array_hzoomx, &sweeps_glo_array_hzoomx,false);
		}
		else{// intree->Draw("time>>h_zoom_x","tag==1");
			Tof_array_hzoomx_entries = GetToF_w_Twin(intree,h_zoom_x,"time:sweeps_global","",1,Tof_Channel2Read,Tof_veto,&Tof_array_hzoomx, &sweeps_glo_array_hzoomx,false);
		}

	}
	else if(tag==0){
		if(treename=="tree0"){
			//intree->Draw("timec>>h_zoom_x","tag==0");
			Tof_array_hzoomx_entries = GetToF_w_Twin(intree,h_zoom_x,"timec:sweeps_global","",0,Tof_Channel2Read,Tof_veto,&Tof_array_hzoomx, &sweeps_glo_array_hzoomx,false);
		}	
		else{ //intree->Draw("time>>h_zoom_x","tag==0");
			Tof_array_hzoomx_entries = GetToF_w_Twin(intree,h_zoom_x,"time:sweeps_global","",0,Tof_Channel2Read,Tof_veto,&Tof_array_hzoomx, &sweeps_glo_array_hzoomx,false);
		}

	}
	else{
		if(treename=="tree0"){
			//intree->Draw("timec>>h_zoom_x","tag==0");
			Tof_array_hzoomx_entries = GetToF_w_Twin(intree,h_zoom_x,"timec:sweeps_global","",-1,Tof_Channel2Read,Tof_veto,&Tof_array_hzoomx, &sweeps_glo_array_hzoomx,false);
		}	
		else{ //intree->Draw("time>>h_zoom_x","tag==0");
			Tof_array_hzoomx_entries = GetToF_w_Twin(intree,h_zoom_x,"time:sweeps_global","",-1,Tof_Channel2Read,Tof_veto,&Tof_array_hzoomx, &sweeps_glo_array_hzoomx,false);
		}
	}

		h_zoom_x->Draw();

		c1->Update();
		c1->cd(2)->SetEditable(kFALSE);

}


void get_para_by_draw(int canvas){

    if(canvas<1||canvas>4) cout<<"canvas index should be [1,4]"<<endl;
    c1->cd(canvas)->SetEditable(kTRUE);
    graphline = (TArrow*)c1->WaitPrimitive("TArrow");

    double X1 =  graphline->GetX1(); double Y1 = graphline->GetY1();
    double X2 =  graphline->GetX2(); double Y2 = graphline->GetY2();

     printf("X1 = %.2f \t Y1 = %.2f \n",X1,Y1);
     printf("X2 = %.2f \t Y2 = %.2f \n",X2,Y2);
	cout<<"peak high = "<<TMath::Max(Y1,Y2)<<endl;
	cout<<"sigma = "<<TMath::Abs(X1-X2)<<endl;
	
	double slop = (Y1-Y2)/(X1-X2);
	double offset = -slop*X2+Y2;

	double exp_p0 = ( X2*TMath::Log(Y1) - X1*TMath::Log(Y2) ) / (X2-X1);
	double exp_p1 = ( TMath::Log(Y1) - TMath::Log(Y2) ) / (X1-X2);
	printf("linearbg: p0 + p1*x : p0= %f, p1 = %f\n",offset,slop);
	printf("expbg: exp(p0 + p1*x) : p0 = %f , p1 = %f\n",exp_p0,exp_p1);
	//graphline->Clear();
	graphline->Delete();

	cout<<"distance: X1-X2: "<<TMath::Abs(X1-X2)<<endl;

	tem_high=TMath::Max(Y1,Y2);
	tem_cento=X1;
	tem_sigma=TMath::Abs(X1-X2);
	tem_slop[0]=slop;		      // [0] for linear bg; [1] for exp bg
	tem_offset[0]=offset;
	tem_slop[1]=exp_p1;
	tem_offset[1]=exp_p0;

	fitrangeL=X1;
	fitrangeR=X2;

	if(canvas==1 || canvas==2)c1->cd(canvas)->SetEditable(kFALSE);
}

void printChi(TH1D* hin, char whichhist, double _fitRangeL, double _fitRangeR){
	int binL = hin->FindBin(_fitRangeL); int binR = hin->FindBin(_fitRangeR);
	int nfreepars = tem_func->GetNumberFreeParameters();
	int ninputpoints = 0;  // number of vaild input data ==> non-zero points
	for(int i=binL;i<=binR;i++){	if(hin->GetBinContent(i)>0) ninputpoints++;	} 
	int NDF = ninputpoints - nfreepars;

	if(NDF<0){ cout<<"\e[1;33m"<<"Too little input: output chi2 only"<<tem_func->GetChisquare()<<"\e[0m"<<endl; }
	else{
			double reduce_chisquare = tem_func->GetChisquare() / NDF;
			if(whichhist == 'r'){
				cout<<"reference: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< NDF <<" = "<< reduce_chisquare <<endl;
				cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,NDF) <<","
					 << TMath::ChisquareQuantile(0.95,NDF)<<"]"<<endl;					 
			}
			else{
				cout<<"X ion: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< NDF <<" = "<< reduce_chisquare <<endl;
				cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,NDF) <<","
					 << TMath::ChisquareQuantile(0.95,NDF)<<"]"<<endl;	
			}
	}

}


#include "FitHelper.h"


#ifdef _SETFUNC1_
void initializer_func1(){

	gau = new TF1("gau","gaus",0,20.e6);  gau->SetParameter(2,10);
	gaus_line = new TF1("gaus_line",func1::gaus_line,0,20.e6,5); gaus_line->SetParameter(2,10);  gaus_line->SetParameter(3,0); gaus_line->SetParameter(4,0);
	man1_func1 = new TF1("man1_func1",func1::fitfunc,0,20.e6,4);  man1_func1->SetParameter(2,10);  man1_func1->SetParameter(3,tctail);
					man1_func1->SetParLimits(3,0,1000);
	man2_func1 = new TF1("man2_func1",func1::fitfunc,0,20.e6,4);  man2_func1->SetParameter(2,10);  man2_func1->SetParameter(3,tctail);
					man2_func1->SetParLimits(3,0,1000);
	man3_func1 = new TF1("man3_func1",func1::fitfunc,0,20.e6,4);  man3_func1->SetParameter(2,10);  man3_func1->SetParameter(3,tctail);
					man3_func1->SetParLimits(3,0,1000);
	ml_func1 = new TF1("ml_func1",func1::fitfunc_line,0,20.e6,6); ml_func1->SetParameter(2,10);  ml_func1->SetParameter(3,tctail);
					ml_func1->SetParLimits(3,0,1000);
	mexp_func1 = new TF1("mexp_func1",func1::fitfunc_exp,0,20.e6,6); mexp_func1->SetParameter(2,10);  mexp_func1->SetParameter(3,tctail);
					mexp_func1->SetParLimits(3,0,1000);

	multi2_func1 = new TF1("multi2_func1",func1::multifunc2,0,20.e6,8);
		multi2_func1->SetParameter(2,10);  multi2_func1->SetParameter(3,tctail);
		multi2_func1->SetParameter(6,10);  multi2_func1->SetParameter(7,tctail);
		multi2_func1->SetParLimits(3,0,1000);	multi2_func1->SetParLimits(7,0,1000);

	multi3_func1 = new TF1("multi3_func1",func1::multifunc3,0,20.e6,12);
		multi3_func1->SetParameter(2,10);  multi3_func1->SetParameter(3,tctail);
		multi3_func1->SetParameter(6,10);  multi3_func1->SetParameter(7,tctail);
		multi3_func1->SetParameter(10,10);  multi3_func1->SetParameter(11,tctail);
		multi3_func1->SetParLimits(3,0,1000);	multi3_func1->SetParLimits(7,0,1000);	multi3_func1->SetParLimits(11,0,1000);

	multi2_exp_func1 = new TF1("multi2_exp_func1",func1::multifunc2_exp,0,20.e6,10);
		multi2_exp_func1->SetParameter(2,10);  multi2_exp_func1->SetParameter(3,tctail);
		multi2_exp_func1->SetParameter(6,10);  multi2_exp_func1->SetParameter(7,tctail);
		multi2_exp_func1->SetParLimits(3,0,1000);	multi2_exp_func1->SetParLimits(7,0,1000);

	multi2_line_func1 = new TF1("multi2_line_func1",func1::multifunc2_line,0,20.e6,10);
		multi2_line_func1->SetParameter(2,10);  multi2_line_func1->SetParameter(3,tctail);
		multi2_line_func1->SetParameter(6,10);  multi2_line_func1->SetParameter(7,tctail);
		multi2_line_func1->SetParLimits(3,0,1000);	multi2_line_func1->SetParLimits(7,0,1000);

}


void fitquickly_func1(char whichhisto, Option_t *Fitoption, bool FitAgain, double _RangeL, double _RangeR){

	int histochoice=2;
	if(whichhisto!='x')histochoice=4;
	static int choice = 0;

	if(FitAgain != true){
		//TF1* tem_func=NULL;
		cout<<"fit func choice: -1: gaus_linearbg; 0: pure gauss; 1. gaus_exp; 2. gaus_exp+linearbg; 3. gaus_exp+expbg; 4. double gaus_exp; 5.double gaus_exp+expbg; 6.double gaus_exp+linebg"<<endl;

		cin>>choice;
		if(choice<-1||choice>6){cout<<"no fit func match!! "<<endl;return;}

		switch(choice){
			case -1:{ NumOfPeaks=1;  SetOneorTwo=0; BackGroundCurve=1;
				char yesorno;
					cout<<"gaus_line is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = gaus_line;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false)	tem_func->SetParameter(2,tem_sigma);
					//tem_func->SetParameter(3,tctail);

					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(3,tem_offset[0]);
					tem_func->SetParameter(4,tem_slop[0]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(3,tem_p0);
						tem_func->SetParameter(4,tem_p1);
					}


					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;

			}
			case 0:{NumOfPeaks=1; SetOneorTwo=0; BackGroundCurve=0;
					char yesorno;
					cout<<"gauss is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = gau;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false)	tem_func->SetParameter(2,tem_sigma);
					//tem_func->SetParameter(3,tctail);

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;

			}
			case 1:{NumOfPeaks=1; BackGroundCurve=0;
					char yesorno;
					cout<<"gasu_exp is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = man1_func1;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
						}
					}
					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;

			}
			case 2:{NumOfPeaks=1; BackGroundCurve=1;
					char yesorno;
					cout<<"gasu_exp+linearbg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = ml_func1;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
						}
					}
					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(4,tem_offset[0]);
					tem_func->SetParameter(5,tem_slop[0]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(4,tem_p0);
						tem_func->SetParameter(5,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 3:{NumOfPeaks=1; BackGroundCurve=2;
					char yesorno;
					cout<<"gasu_exp+expbg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = mexp_func1;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
						}
					}
					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(4,tem_offset[1]);
					tem_func->SetParameter(5,tem_slop[1]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(4,tem_p0);
						tem_func->SetParameter(5,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 4:{NumOfPeaks=2; BackGroundCurve=0;
					char yesorno;
					cout<<"double gasu_exp is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = multi2_func1;  // set func use to fit

					cout<<"draw a line from top of peak1 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
						}
					}
					cout<<"draw a second line from top of peak2 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(4,tem_high);
					tem_func->SetParameter(5,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(6,tem_sigma);
						tem_func->SetParameter(7,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(7,tem_t0);
						}
					}
					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 5:{NumOfPeaks=2; BackGroundCurve=2;
					char yesorno;
					cout<<"double gasu_exp + expbg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = multi2_exp_func1;  // set func use to fit

					cout<<"draw a line from top of peak1 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
						}
					}
					cout<<"draw a second line from top of peak2 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(4,tem_high);
					tem_func->SetParameter(5,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(6,tem_sigma);
						tem_func->SetParameter(7,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(7,tem_t0);
						}
					}

					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(8,tem_offset[1]);
					tem_func->SetParameter(9,tem_slop[1]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(8,tem_p0);
						tem_func->SetParameter(9,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 6:{NumOfPeaks=2; BackGroundCurve=1;
					char yesorno;
					cout<<"double gasu_exp + linebg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = multi2_line_func1;  // set func use to fit

					cout<<"draw a line from top of peak1 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
						}
					}
					cout<<"draw a second line from top of peak2 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(4,tem_high);
					tem_func->SetParameter(5,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(6,tem_sigma);
						tem_func->SetParameter(7,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(7,tem_t0);
						}
					}

					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(8,tem_offset[0]);
					tem_func->SetParameter(9,tem_slop[0]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(8,tem_p0);
						tem_func->SetParameter(9,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}

		}
	}// end of Fitagain != true
	else{

		if(_RangeL != -1 && _RangeR != -1){fitrangeL = _RangeL ; fitrangeR = _RangeR;}
	}


	if(whichhisto=='x'){
		//if(!FitAgain)FitHelper(tem_func,h_zoom_x,fitrangeL,fitrangeR);
		h_zoom_x->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
	}
	else{
		//if(!FitAgain)FitHelper(tem_func,h_zoom_ref,fitrangeL,fitrangeR);
		h_zoom_ref->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
	}

	printf("width for fit: %.1f [ns]\n",fitrangeR-fitrangeL);
	printf("centro_1: %.4f(%.4f) \n",tem_func->GetParameter(1),tem_func->GetParError(1));
	if(choice==4 || choice==5 || choice==6){printf("centro_2: %.4f(%.4f) \n",tem_func->GetParameter(5),tem_func->GetParError(5));}

	// output fitting result to global variables
	if(whichhisto!='x'){// fit result for reference peak
			/*double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
			cout<<"reference: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<<reduce_chisquare<<endl;
			cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
				 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;	*/		

			printChi(h_zoom_ref,whichhisto,fitrangeL,fitrangeR);
			tof_ref_cento = tem_func->GetParameter(1);  tof_ref_cento_err = tem_func->GetParError(1);

	}
	else{// fit result for peak of interest
			/*double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
			cout<<"X ion: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
			cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
				 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;*/

			printChi(h_zoom_x,whichhisto,fitrangeL,fitrangeR); 
			tof_x_cento[0] = tem_func->GetParameter(1); 
			tof_x_cento_err[0] = tem_func->GetParError(1);
			if(choice==4 || choice==5 || choice==6){
				tof_x_cento[1] = tem_func->GetParameter(5); 
				tof_x_cento_err[1] = tem_func->GetParError(5);
			}

	}  


}

#endif



#ifdef _SETFUNC2_
void initializer_func2(){

	man1_func2 = new TF1("man1_func2",func2::fitfunc,0,20.e6,5);  man1_func2->SetParameter(2,10);  man1_func2->SetParameter(3,tctail);   man1_func2->SetParameter(4,tctail);
					man1_func2->SetParLimits(3,0,1000);man1_func2->SetParLimits(4,0,1000);
	man2_func2 = new TF1("man2_func2",func2::fitfunc,0,20.e6,5);  man2_func2->SetParameter(2,10);  man2_func2->SetParameter(3,tctail);   man2_func2->SetParameter(4,tctail);
					man2_func2->SetParLimits(3,0,1000);man2_func2->SetParLimits(4,0,1000);
	man3_func2 = new TF1("man3_func2",func2::fitfunc,0,20.e6,5);  man3_func2->SetParameter(2,10);  man3_func2->SetParameter(3,tctail);   man3_func2->SetParameter(4,tctail);
					man3_func2->SetParLimits(3,0,1000);man3_func2->SetParLimits(4,0,1000);
	ml_func2 = new TF1("ml_func2",func2::fitfunc_line,0,20.e6,7); ml_func2->SetParameter(2,10);  ml_func2->SetParameter(3,tctail);    ml_func2->SetParameter(4,tctail);
					ml_func2->SetParLimits(3,0,1000);ml_func2->SetParLimits(4,0,1000);
	mexp_func2 = new TF1("mexp_func2",func2::fitfunc_exp,0,20.e6,7); mexp_func2->SetParameter(2,10);  mexp_func2->SetParameter(3,tctail);  mexp_func2->SetParameter(4,tctail);
					mexp_func2->SetParLimits(3,0,1000);mexp_func2->SetParLimits(4,0,1000);

	multi2_func2 = new TF1("multi2_func2",func2::multifunc2,0,20.e6,10);
		multi2_func2->SetParameter(2,10);  multi2_func2->SetParameter(3,tctail);   multi2_func2->SetParameter(4,tctail);
		multi2_func2->SetParameter(7,10);  multi2_func2->SetParameter(8,tctail);   multi2_func2->SetParameter(9,tctail);
		multi2_func2->SetParLimits(3,0,1000);multi2_func2->SetParLimits(4,0,1000);
		multi2_func2->SetParLimits(8,0,1000);multi2_func2->SetParLimits(9,0,1000);

	multi3_func2 = new TF1("multi3_func2",func2::multifunc3,0,20.e6,15);
		multi3_func2->SetParameter(2,10);  multi3_func2->SetParameter(3,tctail);   multi3_func2->SetParameter(4,tctail);
		multi3_func2->SetParameter(7,10);  multi3_func2->SetParameter(8,tctail);   multi3_func2->SetParameter(9,tctail);
		multi3_func2->SetParameter(12,10);  multi3_func2->SetParameter(13,tctail);   multi3_func2->SetParameter(14,tctail);
		multi3_func2->SetParLimits(3,0,1000);multi3_func2->SetParLimits(4,0,1000);
		multi3_func2->SetParLimits(8,0,1000);multi3_func2->SetParLimits(9,0,1000);
		multi3_func2->SetParLimits(13,0,1000);multi3_func2->SetParLimits(14,0,1000);

	multi2_exp_func2 = new TF1("multi2_exp_func2",func2::multifunc2_exp,0,20.e6,12);
		multi2_exp_func2->SetParameter(2,10);  multi2_exp_func2->SetParameter(3,tctail);  multi2_exp_func2->SetParameter(4,tctail);
		multi2_exp_func2->SetParameter(7,10);  multi2_exp_func2->SetParameter(8,tctail);  multi2_exp_func2->SetParameter(9,tctail);
		multi2_exp_func2->SetParLimits(3,0,1000);multi2_exp_func2->SetParLimits(4,0,1000);
		multi2_exp_func2->SetParLimits(8,0,1000);multi2_exp_func2->SetParLimits(9,0,1000);

	multi2_line_func2 = new TF1("multi2_line_func2",func2::multifunc2_line,0,20.e6,12);
		multi2_line_func2->SetParameter(2,10);  multi2_line_func2->SetParameter(3,tctail);  multi2_line_func2->SetParameter(4,tctail);
		multi2_line_func2->SetParameter(7,10);  multi2_line_func2->SetParameter(8,tctail);  multi2_line_func2->SetParameter(9,tctail);
		multi2_line_func2->SetParLimits(3,0,1000);multi2_line_func2->SetParLimits(4,0,1000);
		multi2_line_func2->SetParLimits(8,0,1000);multi2_line_func2->SetParLimits(9,0,1000);

	funcRb90 = new TF1("funcRb90",func2::Rb90,0,20.e6,6);	
	   funcRb90->SetParameter(2,10);  funcRb90->SetParameter(3,tctail);  funcRb90->SetParameter(4,tctail);
	   funcRb90->SetParLimits(3,0,1000); funcRb90->SetParLimits(4,0,1000);
}


void fitquickly_func2(char whichhisto, Option_t *Fitoption, bool FitAgain, double _RangeL, double _RangeR){

	int histochoice=2;
	if(whichhisto!='x')histochoice=4;
	static int choice = 1;


	if(FitAgain != true){ // first time to fit
	
		//TF1* tem_func=NULL;
		cout<<"fit func choice: 1. gaus_exp; 2. gaus_exp+linearbg; 3. gaus_exp+expbg; 4. double gaus_exp; 5.double gaus_exp+expbg; 6.double gaus_exp+linebg"<<endl;

		cin>>choice;
		if(choice<1||choice>6){cout<<"no fit func match!! "<<endl;return;}

		switch(choice){
			case 1:{NumOfPeaks=1; BackGroundCurve=0;
					char yesorno;
					cout<<"gasu_exp is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = man1_func2;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);
						tem_func->SetParameter(4,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(4,tem_t0);
						}
					}
					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;

			}
			case 2:{NumOfPeaks=1; BackGroundCurve=1;
					char yesorno;
					cout<<"gasu_exp+linearbg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = ml_func2;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);
						tem_func->SetParameter(4,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(4,tem_t0);
						}
					}
					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(5,tem_offset[0]);
					tem_func->SetParameter(6,tem_slop[0]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(5,tem_p0);
						tem_func->SetParameter(6,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 3:{NumOfPeaks=1; BackGroundCurve=2;
					char yesorno;
					cout<<"gasu_exp+expbg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = mexp_func2;  // set func use to fit

					cout<<"draw a line from top of peak to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);
						tem_func->SetParameter(4,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(4,tem_t0);
						}
					}
					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(5,tem_offset[1]);
					tem_func->SetParameter(6,tem_slop[1]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(5,tem_p0);
						tem_func->SetParameter(6,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 4:{NumOfPeaks=2; BackGroundCurve=0;
					char yesorno;
					cout<<"double gasu_exp is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = multi2_func2;  // set func use to fit

					cout<<"draw a line from top of peak1 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);
						tem_func->SetParameter(4,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(4,tem_t0);
						}
					}
					cout<<"draw a second line from top of peak2 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(5,tem_high);
					tem_func->SetParameter(6,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(7,tem_sigma);
						tem_func->SetParameter(8,tctail);
						tem_func->SetParameter(9,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(8,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(9,tem_t0);
						}
					}
					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 5:{NumOfPeaks=2; BackGroundCurve=2;
					char yesorno;
					cout<<"double gasu_exp + expbg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = multi2_exp_func2;  // set func use to fit

					cout<<"draw a line from top of peak1 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);
						tem_func->SetParameter(4,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(4,tem_t0);
						}
					}
					cout<<"draw a second line from top of peak2 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(5,tem_high);
					tem_func->SetParameter(6,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(7,tem_sigma);
						tem_func->SetParameter(8,tctail);
						tem_func->SetParameter(9,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(8,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(9,tem_t0);
						}
					}

					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(10,tem_offset[1]);
					tem_func->SetParameter(11,tem_slop[1]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(10,tem_p0);
						tem_func->SetParameter(11,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}
			case 6:{NumOfPeaks=2; BackGroundCurve=1;
					char yesorno;
					cout<<"double gasu_exp + linebg is chosen continue or no [y/n]?"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='n'){cout<<"abort and return"<<endl;return;}
					tem_func = multi2_line_func2;  // set func use to fit

					cout<<"draw a line from top of peak1 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0,tem_high);
					tem_func->SetParameter(1,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(2,tem_sigma);
						tem_func->SetParameter(3,tctail);
						tem_func->SetParameter(4,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(4,tem_t0);
						}
					}
					cout<<"draw a second line from top of peak2 to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(5,tem_high);
					tem_func->SetParameter(6,tem_cento);
					if(fixPara == false){	
						tem_func->SetParameter(7,tem_sigma);
						tem_func->SetParameter(8,tctail);
						tem_func->SetParameter(9,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(8,tem_t0);
							cout<<"input new t0 exp tail Left:"<<endl;
							cin>>tem_t0;
							tem_func->SetParameter(9,tem_t0);
						}
					}

					cout<<"draw a line along background"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(10,tem_offset[0]);
					tem_func->SetParameter(11,tem_slop[0]);

					cout<<"change offset and slop manually? [y/n]"<<endl;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new offset and slop:"<<endl;
						double tem_p0,tem_p1;
						cin>>tem_p0>>tem_p1;
						tem_func->SetParameter(10,tem_p0);
						tem_func->SetParameter(11,tem_p1);
					}

					cout<<"draw a line to give fiting range: "<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					break;
			}

		}
	}// end of if of FitAgain
	else{

		if(_RangeL != -1 && _RangeR != -1){fitrangeL = _RangeL; fitrangeR = _RangeR;}
	}



	if(whichhisto=='x'){
		//if(!FitAgain)FitHelper(tem_func,h_zoom_x,fitrangeL,fitrangeR);
		h_zoom_x->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
	}
	else{
		//if(!FitAgain)FitHelper(tem_func,h_zoom_ref,fitrangeL,fitrangeR);
		h_zoom_ref->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
	}

	printf("width for fit: %.1f [ns]\n",fitrangeR-fitrangeL);
	printf("centro_1: %.4f(%.4f) \n",tem_func->GetParameter(1),tem_func->GetParError(1));
	if(choice==4 || choice==5 || choice==6){printf("centro_2: %.4f(%.4f) \n",tem_func->GetParameter(6),tem_func->GetParError(6));}

	// output fitting result to global variables
	if(whichhisto!='x'){// fit result for reference peak
			/*double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
			cout<<"reference: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
			cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
				 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;*/

			printChi(h_zoom_ref,whichhisto,fitrangeL,fitrangeR);
			tof_ref_cento = tem_func->GetParameter(1);  tof_ref_cento_err = tem_func->GetParError(1);

	}
	else{// fit result for peak of interest
			/*double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
			cout<<"X ion: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
			cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
				 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;*/

			printChi(h_zoom_x,whichhisto,fitrangeL,fitrangeR);	 

			tof_x_cento[0] = tem_func->GetParameter(1); 
			tof_x_cento_err[0] = tem_func->GetParError(1);
			if(choice==4 || choice==5 || choice==6){
				tof_x_cento[1] = tem_func->GetParameter(6); 
				tof_x_cento_err[1] = tem_func->GetParError(6);
			}

	}  


}


#endif

/*
void fitgeneralcurve(char whichhisto,int func1ORfunc2, int NumPeaks2Fit, int whichbackground, Option_t *Fitoption, bool FitAgain, double _RangeL, double _RangeR){

		SetOneorTwo = func1ORfunc2;   // 0 = gaus;     1 = func1 (gauss + exp Right tail) ;   2 = func2 (gauss + exp Ritht and Left tail)
		NumOfPeaks = NumPeaks2Fit;    // at least 1 peak
		BackGroundCurve = whichbackground;
		if(SetOneorTwo < 0 || SetOneorTwo >2){ cerr<<" error: choice of fit func group 1 or 2 ; or 0 for pure gaus!!"<<endl; return;}
		if(NumOfPeaks <1){ cerr<<" Number of Peaks to fit > = 1   !!"<<endl; return;}
		if(BackGroundCurve<0 || BackGroundCurve >2){cerr<<" error: choice of background 0=> no bg; 1 => linearbg or 2 => expbg !!"<<endl; return;}

		static TF1 * general_curve = NULL;
		if(general_curve != NULL && FitAgain == false){ delete general_curve; general_curve = NULL; }

		// &&&&&&&& define fit function &&&&&&&&&&&
		if(general_curve == NULL){
			if(SetOneorTwo == 0){  // choose pure gaus
				if(BackGroundCurve == 0)	general_curve = new TF1("general_curve",GeneralCurve,0,50e6,3*NumOfPeaks);
				else general_curve = new TF1("general_curve",GeneralCurve,0,50e6,3*NumOfPeaks+2);

			}
			else if(SetOneorTwo == 1){  // choose func1
				if(BackGroundCurve == 0)	general_curve = new TF1("general_curve",GeneralCurve,0,50e6,4*NumOfPeaks);
				else general_curve = new TF1("general_curve",GeneralCurve,0,50e6,4*NumOfPeaks+2);
				for(int i=0;i<NumOfPeaks;i++){general_curve->SetParLimits(i*4+3,0,1000);}

			}
			else{ // choose func2
				if(BackGroundCurve == 0)	general_curve = new TF1("general_curve",GeneralCurve,0,50e6,5*NumOfPeaks);
				else general_curve = new TF1("general_curve",GeneralCurve,0,50e6,5*NumOfPeaks+2);
				for(int i=0;i<NumOfPeaks;i++){
					general_curve->SetParLimits(i*5+3,0,1000);
					general_curve->SetParLimits(i*5+4,0,1000);
				}
			}
		}// end of define fit function

		tem_func = general_curve; // pass pointer of fit function to global handler

		int peak_para_offset=0;  // different group of fit function has different number of paras

		int histochoice=2; // activate corresponding Canvas
		if(whichhisto!='x')histochoice=4;  // activate corresponding Canvas

		if(FitAgain != true){

			// multi peaks section
			if(SetOneorTwo == 0 || SetOneorTwo == 1){ // 0 => pure gaus ; 1 => func1(guas + Right tail)
					
				for(int PeakIndex=0;PeakIndex<NumOfPeaks;peak_para_offset = (4+(SetOneorTwo-1))*(++PeakIndex)){ // already point to next group paras when junp out from loop
					cout<<"draw a line from top of peak "<<PeakIndex+1<<" to FWHM"<<endl;
					cout<<"pending..."<<endl;
					
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0+peak_para_offset,tem_high);
					tem_func->SetParameter(1+peak_para_offset,tem_cento);
					tem_func->SetParameter(2+peak_para_offset,tem_sigma);
					if(SetOneorTwo ==1){
						tem_func->SetParameter(3+peak_para_offset,tctail);
						//tem_func->SetParameter(4,tctail);

						cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
						char yesorno;
						do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
						if(yesorno=='y'){
							cout<<"input new t0 exp tail Right:"<<endl;
							double tem_t0;
							cin>>tem_t0;
							tem_func->SetParameter(3+peak_para_offset,tem_t0);
							//cout<<"input new t0 exp tail Left:"<<endl;
							//cin>>tem_t0;
							//tem_func->SetParameter(4,tem_t0);
						}
					}
				}


			} // end of SetOneorTwo ==0 || 1
			else{ // when SetOneorTwo ==2

				for(int PeakIndex=0;PeakIndex<NumOfPeaks;peak_para_offset = 5*(++PeakIndex)){
					cout<<"draw a line from top of peak "<<PeakIndex+1<<" to FWHM"<<endl;
					cout<<"pending..."<<endl;
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0+peak_para_offset,tem_high);
					tem_func->SetParameter(1+peak_para_offset,tem_cento);
					tem_func->SetParameter(2+peak_para_offset,tem_sigma);
					tem_func->SetParameter(3+peak_para_offset,tctail);
					tem_func->SetParameter(4+peak_para_offset,tctail);

					cout<<"tc = "<<tctail<<" ; change or not [y/n]?"<<endl;
					char yesorno;
					do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
					if(yesorno=='y'){
						cout<<"input new t0 exp tail Right:"<<endl;
						double tem_t0;
						cin>>tem_t0;
						tem_func->SetParameter(3+peak_para_offset,tem_t0);
						cout<<"input new t0 exp tail Left:"<<endl;
						cin>>tem_t0;
						tem_func->SetParameter(4+peak_para_offset,tem_t0);
					}
				}
			}// end of SetOneorTwo ==2

			// background section
			if(BackGroundCurve !=0){
				cout<<"draw a line along background"<<endl;
				cout<<"pending..."<<endl;
				get_para_by_draw(histochoice);
				tem_func->SetParameter(0+peak_para_offset,tem_offset[BackGroundCurve-1]);  //tem_offset[0] => linear; tem_offset[1] => exp
				tem_func->SetParameter(1+peak_para_offset,tem_slop[BackGroundCurve-1]);   // same as above
				cout<<"change offset and slop manually? [y/n]"<<endl;
				char yesorno;
				do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
				if(yesorno=='y'){
					cout<<"input new offset and slop:"<<endl;
					double tem_p0,tem_p1;
					cin>>tem_p0>>tem_p1;
					tem_func->SetParameter(0+peak_para_offset,tem_p0);
					tem_func->SetParameter(1+peak_para_offset,tem_p1);
				}
			}

			cout<<"draw a line to give fiting range: "<<endl;
			cout<<"pending..."<<endl;
			get_para_by_draw(histochoice);

		}// end of fitagain != true
		else{

			if(_RangeL != -1 && _RangeR != -1){fitrangeL = _RangeL; fitrangeR = _RangeR;}
		}


		// start fit procedures
		if(whichhisto=='x'){
			h_zoom_x->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
		}
		else{
			h_zoom_ref->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
		}


		for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
				peak_para_offset = PeakIndex * (4+(SetOneorTwo-1));
				printf("centro_%d: %.4f(%.4f) \n",PeakIndex+1,tem_func->GetParameter(1+peak_para_offset),tem_func->GetParError(1+peak_para_offset));
		}



		// output fitting result to global variables
		if(whichhisto!='x'){// fit result for reference peak
				double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
				cout<<"reference: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
				cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
					 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;

				tof_ref_cento = tem_func->GetParameter(1);  tof_ref_cento_err = tem_func->GetParError(1);

		}
		else{// fit result for peak of interest
				double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
				cout<<"X ion: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
				cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
					 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;

				for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
						peak_para_offset = PeakIndex * (4+(SetOneorTwo-1));
						tof_x_cento[PeakIndex] = tem_func->GetParameter(1+peak_para_offset);
						tof_x_cento_err[PeakIndex] = tem_func->GetParError(1+peak_para_offset);
						
				}
		}  


}*/



void fitgeneralcurve(char whichhisto,int func1ORfunc2, int NumPeaks2Fit, int whichbackground, int Index_MainPeak,bool _ParasLock,Option_t *Fitoption, bool FitAgain, double _RangeL, double _RangeR){


		if(func1ORfunc2 < 0 || func1ORfunc2 >2){ cerr<<" error: choice of fit func group 1 or 2 ; or 0 for pure gaus!!"<<endl; return;}
		if(NumPeaks2Fit <1){ cerr<<" Number of Peaks to fit > = 1   !!"<<endl; return;}
		if(whichbackground<0 || whichbackground >2){cerr<<" error: choice of background 0=> no bg; 1 => linearbg or 2 => expbg !!"<<endl; return;}
		if(Index_MainPeak>NumPeaks2Fit){cerr<<"error: mainpeak index should be <= "<< NumPeaks2Fit<<endl;}

		SetOneorTwo = func1ORfunc2;   // 0 = gaus;     1 = func1 (gauss + exp Right tail) ;   2 = func2 (gauss + exp Ritht and Left tail)
		NumOfPeaks = NumPeaks2Fit;    // at least 1 peak
		BackGroundCurve = whichbackground;
		MainPeakIndex = Index_MainPeak;
		ParasLock = _ParasLock;

		static TF1 * general_curve = NULL;
		if(general_curve != NULL && FitAgain == false){ delete general_curve; general_curve = NULL; }

		int NumOfPars_atTF = NumOfPars_TF_cal(SetOneorTwo , NumOfPeaks,BackGroundCurve, ParasLock);


		// &&&&&&&& define fit function &&&&&&&&&&&
		if(general_curve == NULL){
			general_curve = new TF1("general_curve",GeneralCurve,0,50e6,NumOfPars_atTF);

			if(SetOneorTwo == 1){
				if(ParasLock){ general_curve->SetParLimits((MainPeakIndex-1)*2+3, 0, 1000);}  // tailR limits[0,1000] , para lock case
				else{
					for(int i=0;i<NumOfPeaks;i++){general_curve->SetParLimits(i*(SetOneorTwo+3)+3,0,1000);}  // tailR limit, para unlock case
				}
			}
			else{
				if(ParasLock){
					general_curve->SetParLimits((MainPeakIndex-1)*2+3, 0, 1000); // locked para ; tailR limits
					general_curve->SetParLimits((MainPeakIndex-1)*2+4, 0, 1000); // locked para; tailL limits
				}
				else{
						for(int i=0;i<NumOfPeaks;i++){
							general_curve->SetParLimits(i*(SetOneorTwo+3)+3,0,1000);
							general_curve->SetParLimits(i*(SetOneorTwo+3)+4,0,1000);
						}
				}
			}
		}


		tem_func = general_curve; // pass pointer of fit function to global handler

		int peak_para_offset=0;  // different group of fit function has different number of paras

		int histochoice=2; // activate corresponding Canvas
		if(whichhisto!='x')histochoice=4;  // activate corresponding Canvas

		if(FitAgain != true){

			// multi peaks section
					
				for(int PeakIndex=0;PeakIndex<NumOfPeaks;peak_para_offset = Paras_offset_cal(++PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock) ){ // already point to next group paras when junp out from loop
					cout<<"draw a line from top of peak "<<PeakIndex+1<<" to FWHM"<<endl;
					cout<<"pending..."<<endl;
					
					get_para_by_draw(histochoice);
					tem_func->SetParameter(0+peak_para_offset,tem_high);
					tem_func->SetParameter(1+peak_para_offset,tem_cento);
					tem_func->SetParName(0+peak_para_offset,Form("Amp_Peak%d",PeakIndex+1));
					tem_func->SetParName(1+peak_para_offset,Form("centro_Peak%d",PeakIndex+1));

					if(PeakIndex != (MainPeakIndex-1) && ParasLock) continue;
					else{
							tem_func->SetParameter(2+peak_para_offset,tem_sigma);
							tem_func->SetParName(2+peak_para_offset,Form("sigma_Peak%d",PeakIndex+1));
						if(SetOneorTwo !=0){
							tem_func->SetParameter(3+peak_para_offset,tem_sigma*1.2);//tctail
							tem_func->SetParName(3+peak_para_offset,Form("tctailR_Peak%d",PeakIndex+1));
							//tem_func->SetParameter(4,tctail);
							if(SetOneorTwo ==2){
								tem_func->SetParameter(4+peak_para_offset,tem_sigma*1.2);
								tem_func->SetParName(4+peak_para_offset,Form("tctailL_Peak%d",PeakIndex+1));
							}

							cout<<"tc = "<<tem_sigma*1.2<<" ; change or not [y/n]?"<<endl;
							char yesorno;
							do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
							if(yesorno=='y'){
								cout<<"input new tc exp tail Right:"<<endl;
								double tem_t0;
								cin>>tem_t0;
								tem_func->SetParameter(3+peak_para_offset,tem_t0);

								if(SetOneorTwo ==2){
									cout<<"input new tc exp tail Left:"<<endl;
									cin>>tem_t0;
									tem_func->SetParameter(4+peak_para_offset,tem_t0);
								}					
								//cout<<"input new t0 exp tail Left:"<<endl;
								//cin>>tem_t0;
								//tem_func->SetParameter(4,tem_t0);
							}
						}
					}
				}


			// background section
			if(BackGroundCurve !=0){
				cout<<"draw a line along background"<<endl;
				cout<<"pending..."<<endl;
				get_para_by_draw(histochoice);
				tem_func->SetParameter(0+peak_para_offset,tem_offset[BackGroundCurve-1]);  //tem_offset[0] => linear; tem_offset[1] => exp
				tem_func->SetParameter(1+peak_para_offset,tem_slop[BackGroundCurve-1]);   // same as above
				if(BackGroundCurve ==1){
					tem_func->SetParName(0+peak_para_offset,"linearBG_offset");
					tem_func->SetParName(1+peak_para_offset,"linearBG_slop");
				}
				else{
					tem_func->SetParName(0+peak_para_offset,"expBG_offset");
					tem_func->SetParName(1+peak_para_offset,"expBG_slop");
				}

				cout<<"change offset and slop manually? [y/n]"<<endl;
				char yesorno;
				do{cin>>yesorno;}while(yesorno != 'y' && yesorno != 'n');
				if(yesorno=='y'){
					cout<<"input new offset and slop:"<<endl;
					double tem_p0,tem_p1;
					cin>>tem_p0>>tem_p1;
					tem_func->SetParameter(0+peak_para_offset,tem_p0);
					tem_func->SetParameter(1+peak_para_offset,tem_p1);
				}
			}

			cout<<"draw a line to give fiting range: "<<endl;
			cout<<"pending..."<<endl;
			get_para_by_draw(histochoice);

		}// end of fitagain != true
		else{

			if(_RangeL != -1 && _RangeR != -1){fitrangeL = _RangeL; fitrangeR = _RangeR;}
		}


		// start fit procedures
		printf("width for fit: %.1f [ns]\n",fitrangeR-fitrangeL);
		if(whichhisto=='x'){
			h_zoom_x->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
		}
		else{
			h_zoom_ref->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
		}


		for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
				peak_para_offset = Paras_offset_cal(PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock);
				//peak_para_offset = PeakIndex * (4+(SetOneorTwo-1));
				printf("centro_%d: %.4f(%.4f) \n",PeakIndex+1,tem_func->GetParameter(1+peak_para_offset),tem_func->GetParError(1+peak_para_offset));
		}



		// output fitting result to global variables
		if(whichhisto!='x'){// fit result for reference peak
				/*double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
				cout<<"reference: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
				cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
					 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;*/

				printChi(h_zoom_ref,whichhisto,fitrangeL,fitrangeR);

				tof_ref_cento = tem_func->GetParameter(1);  tof_ref_cento_err = tem_func->GetParError(1);

		}
		else{// fit result for peak of interest
				/*double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
				cout<<"X ion: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
				cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
					 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;*/

				printChi(h_zoom_x,whichhisto,fitrangeL,fitrangeR);

				for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
						peak_para_offset = Paras_offset_cal(PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock);
						//peak_para_offset = PeakIndex * (4+(SetOneorTwo-1));
						tof_x_cento[PeakIndex] = tem_func->GetParameter(1+peak_para_offset);
						tof_x_cento_err[PeakIndex] = tem_func->GetParError(1+peak_para_offset);
						
				}
		}  


}


void initializer(){
	//	ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000);
	//	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
		initializer_func1();
		initializer_func2();
		for(int i=0;i<30;i++){	FHistory[i]=NULL;}

}

void fitquickly(char whichhisto,int func1ORfunc2,Option_t *Fitoption, bool FitAgain, double _RangeL, double _RangeR,int N_loop){
	for(int i =0;i<N_loop;i++){
		if(func1ORfunc2==1){SetOneorTwo=1;ParasLock = false;fitquickly_func1(whichhisto, Fitoption, FitAgain, _RangeL, _RangeR);}
		else if(func1ORfunc2==2){SetOneorTwo=2;ParasLock = false;fitquickly_func2(whichhisto, Fitoption, FitAgain, _RangeL, _RangeR);}
		else{cout<<"choose Set of func1 for 1;  Set of func2 for 2;"<<endl;}

		if(i==N_loop-1) printf("use ?: fitquickly('%c',%d,\"%s\",true,%.2f,%.2f,50)\n",whichhisto,func1ORfunc2,Fitoption,fitrangeL,fitrangeR);

	}
}


void mass_calculator(int Peak1OrPeak2){
	if(Peak1OrPeak2 <1 && Peak1OrPeak2 >10){cout<< "Peak number input err! [1 - 10]"<<endl; return;}

	int lapwidth=6;
	int lapl=laps_ref-lapwidth; int laph=laps_ref+lapwidth;

	// calculate three different charge state
	cout<<"lap_ref"<<laps_ref<<endl;
	printf("laps_xion \t mass_over_q \t Amass_x+ \t Amass_x++ \t Amass_x+++ \n");
  for(int searchlap=lapl;searchlap<=laph;searchlap++){
	for(int _q_x=1;_q_x<4;_q_x++){
		// ion mass of X
		/*double m_x_cal = ((m_ref - m_ele*q_ref)/ q_ref) * TMath::Power((tof_x_cento[Peak1OrPeak2-1] - t0) / (tof_ref_cento+bref*(searchlap-laps_ref) - t0),2)*_q_x;  
		double moq = m_x_cal / _q_x;
		double moq_err = ErrCal_moq(moq,tof_x_cento[Peak1OrPeak2-1],tof_x_cento_err[Peak1OrPeak2-1],tof_ref_cento,tof_ref_cento_err,m_ref,err_ref);
		// atom mass of X
		m_x_cal = m_x_cal + m_ele * _q_x;
		double m_x_cal_err = ErrCal_mass(moq_err,_q_x);*/

		double* mass_extract = mass_calculator(tof_x_cento[Peak1OrPeak2-1],tof_x_cento_err[Peak1OrPeak2-1],searchlap,_q_x,
														tof_ref_cento,tof_ref_cento_err,laps_ref,m_ref,err_ref, bref,bref_err,q_ref,false);
		double moq = mass_extract[2] * mass_extract[2] * (m_ref - m_ele*q_ref)/ q_ref;
		double m_x_cal = mass_extract[0];
		double m_x_cal_err = mass_extract[1];
	
		if(_q_x==1)printf("%d \t %.2f \t %.3f(%.3f) \t",searchlap,moq,m_x_cal,m_x_cal_err);
		if(_q_x==2)printf("%.3f(%.3f) \t",m_x_cal,m_x_cal_err);
		if(_q_x==3)printf("%.3f(%.3f) \n",m_x_cal,m_x_cal_err);
	}
   }
}


void mass_calculator(double tof0,double tof0_err,double tof_ref,double tof_ref_err,int _laps_ref, double mass_ref, double mass_ref_err,
								double b_ref, double b_ref_err, int _q_ref){

	int lapwidth=3;
	int lapl=_laps_ref-lapwidth; int laph=_laps_ref+lapwidth;

	// calculate three different charge state
	cout<<"lap_ref"<<_laps_ref<<endl;
	printf("laps_xion \t  mass_over_q \t Amass_x+ \t Amass_x++ \t Amass_x+++\n");
  for(int searchlap=lapl;searchlap<=laph;searchlap++){
	for(int _q_x=1;_q_x<4;_q_x++){
		// ion mass of X
		/*
		double m_x_cal = ((mass_ref - m_ele*_q_ref)/ _q_ref) * TMath::Power((tof0 - t0) / (tof_ref+b_ref*(searchlap-_laps_ref) - t0),2)*_q_x;  
		double moq = m_x_cal / _q_x;
		double moq_err = ErrCal_moq(moq,tof0,tof0_err,tof_ref,tof_ref_err,mass_ref,mass_ref_err);

		// atom mass of X
		m_x_cal = m_x_cal + m_ele * _q_x;
		double m_x_cal_err = ErrCal_mass(moq_err,_q_x);  */

		double* mass_extract = mass_calculator(tof0,tof0_err,searchlap,_q_x,tof_ref,tof_ref_err,_laps_ref,mass_ref,mass_ref_err,b_ref,b_ref_err,_q_ref);
		double moq = mass_extract[2] * mass_extract[2] * (mass_ref - m_ele*_q_ref)/ _q_ref;
		double m_x_cal = mass_extract[0];
		double m_x_cal_err = mass_extract[1];
	
		if(_q_x==1)printf("%d \t %.2f \t %.3f(%.3f) \t",searchlap,moq,m_x_cal,m_x_cal_err);
		if(_q_x==2)printf("%.3f(%.3f) \t",m_x_cal,m_x_cal_err);
		if(_q_x==3)printf("%.3f(%.3f) \n",m_x_cal,m_x_cal_err);
	}
   }
}


double * mass_calculator(double tof0,double tof0_err,int _laps_x, int _q_x, double tof_ref,double tof_ref_err,int _laps_ref, double mass_ref, double mass_ref_err,
										double b_ref, double b_ref_err,int _q_ref, bool verbal){
	static double mass_return[4];
	mass_return[0]=0;   // [0] => mass value
	mass_return[1]=0;	// [1] => mass_err value;
	mass_return[2]=0;  // [2] =>rho value
	mass_return[3]=0;  //[3] => rho_err value;

	if(verbal){
			cout<<"lap_ref ="<<_laps_ref<<endl;
			printf("laps_xion: %d \t  charge: %d\n\n",_laps_x,_q_x);
	}

		// ion mass of X

		double tof_ref_xlap = tof_ref + b_ref*(_laps_x-_laps_ref);
		double tof_ref_xlap_err = TMath::Sqrt( tof_ref_err * tof_ref_err + (_laps_x-_laps_ref)*(_laps_x-_laps_ref)*b_ref_err*b_ref_err  );

		double rho_cal = (tof0 - t0) / (tof_ref_xlap - t0);
		double rho_cal_err = ErrCal_rho(rho_cal,tof0,tof0_err,tof_ref_xlap,tof_ref_xlap_err); //rho,  tofx, tofx_err, tofref, tofref_err

		/*
		double m_x_cal = ((mass_ref - m_ele*_q_ref)/ _q_ref) * TMath::Power((tof0 - t0) / (tof_ref+b_ref*(_laps_x-_laps_ref) - t0),2)*_q_x;  
		double moq = m_x_cal / _q_x;
		double moq_err = ErrCal_moq(moq,tof0,tof0_err,tof_ref,tof_ref_err,mass_ref,mass_ref_err);*/

		double m_x_cal = ((mass_ref - m_ele*_q_ref)/ _q_ref) *rho_cal * rho_cal * _q_x;
		double moq = m_x_cal / _q_x;
		double moq_err = ErrCal_moq(rho_cal,rho_cal_err,mass_ref,mass_ref_err,_q_ref);// _rho,  _rho_err,  mass_ref,  mass_ref_err,  _q_ref

		// atom mass of X
		m_x_cal = m_x_cal + m_ele * _q_x;
		double m_x_cal_err = ErrCal_mass(moq_err,_q_x);       
	
		if(verbal){
			printf("X_ion laps: %d \t  charge: %d\n",_laps_x,_q_x);
			printf("X_mass: %.3f(%.3f) \n",m_x_cal,m_x_cal_err);
		}

		mass_return[0] = m_x_cal;
		mass_return[1] = m_x_cal_err;

		mass_return[2]=rho_cal;
		mass_return[3]=rho_cal_err;
		
		return mass_return;
	
}

double * mass_calculator(double tof0,double tof0_err,int _q_x, Ion *RefIon, int Num_RefIons){
	// use TGraphErr to fit, find out parameters which can measure t0 using mulitple reference, assuming tof function is:
	// tof = slop * sqrt(m/q) + t0 ==> Y = slop * X + offset;  == t = A *X +B;
	// pol1 => p0 = offset; p1 = slop;
	const int N_refion_allow = 10;
	double X[N_refion_allow] ={0}, X_err[N_refion_allow] ={0} ;
	double Y[N_refion_allow] ={0}, Y_err[N_refion_allow] ={0} ;
	static double mass_return[2];
	mass_return[0] =0;
	mass_return[1] =0;
	if(Num_RefIons<1 || Num_RefIons >10){cerr<<"Number of Ref ions input error!!, [1,10]"; return mass_return;}

	if(Num_RefIons ==1){// noted that bref, and bref_err from global
		return mass_calculator(tof0,tof0_err,laps_ref, _q_x, RefIon[0].tof,RefIon[0].tof_err,laps_ref, RefIon[0].mass, RefIon[0].mass_err,bref,bref_err,RefIon[0].charge);
	}
	else if(Num_RefIons>2){
			for(int i=0;i<Num_RefIons;i++){
				X[i] = TMath::Sqrt((RefIon[i].mass - RefIon[i].charge * m_ele) / RefIon[i].charge);  //    => sqrt(m/q)
				X_err[i] = TMath::Sqrt( RefIon[i].mass_err * RefIon[i].mass_err + (RefIon[i].charge * err_me) * (RefIon[i].charge * err_me) );  // temporary term
				X_err[i] = X_err[i] * TMath::Sqrt( RefIon[i].charge/(RefIon[i].mass-RefIon[i].charge * m_ele) ) / (2*RefIon[i].charge);  
				Y[i] = RefIon[i].tof;
				Y_err[i] = RefIon[i].tof_err;

				printf("X[%d]= %.3f(%.3f) ; Y[%d] = %.3f(%.3f)\n",i,X[i],X_err[i],i,Y[i],Y_err[i]);
			}


			// method 1; use TGraph to fit => Get A and B(t0)  ==>  tof = A * sqrt(m/q) + B ==> m/q = ((tof-B)/A)^2
			TGraphErrors * gr = new TGraphErrors(Num_RefIons,X,Y,X_err,Y_err);
			TF1 * mass_tof_line = new TF1("mass_tof_line","pol1",0,50e9);

			gr->Fit(mass_tof_line,"NQ");			gr->Fit(mass_tof_line,"NQ");
			double A = mass_tof_line->GetParameter(1) , A_err = mass_tof_line->GetParError(1);// slop
			double B = mass_tof_line->GetParameter(0), B_err = mass_tof_line->GetParError(0);// offset	==> t0		

			printf("A = %.3f(%.3f) ; B(t0) = %.3f(%.3f)\n",A,A_err,B,B_err);

			double moq = TMath::Power((tof0 - B)/A,2);

			// calculate Error of moq
			double term1 = TMath::Power(2*(tof0 - B)/(A*A) * tof0_err,2);
			double term2 = TMath::Power(2*(tof0 - B)/(A*A) * B_err,2);
			double term3 = TMath::Power(2*(tof0 - B)*(tof0 - B)/(A*A*A) * A_err,2);

			double moq_err = TMath::Sqrt(term1 + term2 + term3 );

			printf("moq = %.3f(%.3f)\n",moq,moq_err);

			mass_return[0] = moq * _q_x + _q_x * m_ele;  // atomic mass
			mass_return[1] =  _q_x * TMath::Sqrt(moq_err * moq_err + err_me * err_me);

			gr->Delete();
			mass_tof_line->Delete();

			cout<<endl;
			cout<<"Measure result by "<< Num_RefIons << "reference ions:"<<endl;
			printf("mass: %.3f(%.3f)\n",mass_return[0],mass_return[1]);
			printf("t0: %.2f(%.2f)\n",B,B_err);

	
	}
	else{
			// method 2: cancel A and B
			for(int i=0;i<Num_RefIons;i++){
				X[i] = TMath::Sqrt((RefIon[i].mass - RefIon[i].charge * m_ele) / RefIon[i].charge);  //    => sqrt(m/q)
				X_err[i] = TMath::Sqrt( RefIon[i].mass_err * RefIon[i].mass_err + (RefIon[i].charge * err_me) * (RefIon[i].charge * err_me) );  // temporary term
				X_err[i] = X_err[i] * TMath::Sqrt( RefIon[i].charge/(RefIon[i].mass-RefIon[i].charge * m_ele) ) / (2*RefIon[i].charge);  
				Y[i] = RefIon[i].tof;
				Y_err[i] = RefIon[i].tof_err;

				//printf("X[%d]= %.3f(%.3f) ; Y[%d] = %.3f(%.3f)\n",i,X[i],X_err[i],i,Y[i],Y_err[i]);
			}
			// just for easy calculate t0 and t0_err
			TGraphErrors * gr = new TGraphErrors(Num_RefIons,X,Y,X_err,Y_err);
			TF1 * mass_tof_line = new TF1("mass_tof_line","pol1",0,50e9);

			gr->Fit(mass_tof_line,"NQ");			gr->Fit(mass_tof_line,"NQ");
			//double A = mass_tof_line->GetParameter(1) , A_err = mass_tof_line->GetParError(1);// slop
			double B = mass_tof_line->GetParameter(0), B_err = mass_tof_line->GetParError(0);// offset	==> t0	


			double f = ((X[0]-X[1]) * tof0 + X[1]*Y[0] - X[0]*Y[1]) / (Y[0] - Y[1]);
			double moq = f * f;

			// moq_err calculation
			double term1 = TMath::Power(2 * f* (tof0-Y[1])/(Y[0]-Y[1]) * X_err[0],2);
			double term2 = TMath::Power(2 * f* (-tof0+Y[0])/(Y[0]-Y[1]) * X_err[1],2);
			double term3 = TMath::Power( 2* f* (X[1]-X[0])*(tof0-Y[1])/((Y[0]-Y[1])*(Y[0]-Y[1])) * Y_err[0],2);
			double term4 = TMath::Power( 2* f* (X[0]-X[1])*(tof0-Y[0])/((Y[0]-Y[1])*(Y[0]-Y[1])) * Y_err[1],2);
			double term5 = TMath::Power( 2* f* (X[0]-X[1])/(Y[0]-Y[1]) * tof0_err,2);
			double moq_err = TMath::Sqrt(term1 + term2 + term3 + term4 + term5);

			mass_return[0] = moq * _q_x + _q_x * m_ele;  // atomic mass
			mass_return[1] =  _q_x * TMath::Sqrt(moq_err * moq_err + err_me * err_me);

			gr->Delete();
			mass_tof_line->Delete();

			cout<<"Measure result by "<< Num_RefIons << "reference ions:"<<endl;
			printf("mass: %.3f(%.3f)\n",mass_return[0],mass_return[1]);
			printf("t0: %.2f(%.2f)\n",B,B_err);
	}


	return mass_return;


}


double ErrCal_rho(double _rho, double tofx, double tofx_err, double tofref, double tofref_err){

	double term1 = TMath::Power(tofx_err/ (tofx - t0), 2);
	double term2 = TMath::Power(tofref_err/ (tofref - t0), 2);
	double term3 = TMath::Power( (tofx-tofref)* err_t0 / ( (tofx-t0)*(tofref-t0) ), 2);

	return _rho * TMath::Sqrt(term1 + term2 + term3);

}


double ErrCal_moq(double _rho, double _rho_err, double mass_ref, double mass_ref_err, int _q_ref){
	double moq_ref = (mass_ref - _q_ref * m_ele) / _q_ref;
	double moq_ref_err = TMath::Sqrt((mass_ref_err/_q_ref)*(mass_ref_err/_q_ref) + err_me * err_me);

	double term1 = (2*moq_ref*_rho_err) * (2*moq_ref*_rho_err);
	double term2 = (_rho*moq_ref_err) * (_rho*moq_ref_err);

	return _rho * TMath::Sqrt(term1 + term2);

}



/*
double ErrCal_moq(double moq_x, double tofx, double tofx_err, double tofref, double tofref_err, double mass_ref, double mass_ref_err){

		double moq_ref = (mass_ref - q_ref * m_ele) / q_ref;
		double moq_ref_err = TMath::Sqrt((mass_ref_err/q_ref)*(mass_ref_err/q_ref) + err_me * err_me);

		double term1 = TMath::Power(2.0 / (tofx-t0),2) *  (tofx_err*tofx_err);
		double term2 = TMath::Power(2.0 / (tofref-t0),2) *  (tofref_err * tofref_err);
		double term3 = TMath::Power(moq_ref_err / moq_ref ,2);
		double term4 = TMath::Power(2.0 *(tofx-tofref)* err_t0 / ( (tofx-t0)*(tofref-t0) ),2);

		return moq_x * TMath::Sqrt(term1 + term2 + term3 + term4);

}*/

double ErrCal_mass(double moq_x_err,double _q_x){

		return _q_x * TMath::Sqrt( moq_x_err * moq_x_err + err_me * err_me);

}

void PrintInfo(){
	cout<<endl;
	cout<<"reference:"<<endl;
	printf("mass_ref:%.3f(%.3f) \t q_ref: %d \t laps_ref:%d \t b value:%.3f\n",m_ref,err_ref,q_ref,laps_ref,bref);
	printf("tof:%.4f(%.4f)",tof_ref_cento,tof_ref_cento_err);
	cout<<endl;
	cout<<endl;
	cout<<"X Ion:"<<endl;
	printf("laps_x = %d \t q_x = %d\n",laps_x,q_x);
	printf("1: tof_x:%.4f(%.4f) \n",tof_x_cento[0],tof_x_cento_err[0]);
	printf("2: tof_x:%.4f(%.4f) \n",tof_x_cento[1],tof_x_cento_err[1]);

}

void ResetRef(double mass_ref,double mass_ref_err,int _qref,int _laps_ref,double b_ref,double tofref,double tofref_err){
	if(mass_ref != -1) 	m_ref = mass_ref;
	if(mass_ref_err != -1) err_ref = mass_ref_err;
	if(_qref!=-1) q_ref = _qref;
	if(_laps_ref != -1) laps_ref = _laps_ref;
	if(b_ref != -1) bref = b_ref;
	if(tofref != -1) tof_ref_cento = tofref;
	if(tofref_err != -1) tof_ref_cento_err = tofref_err;
	
	PrintInfo();
}

void ResetRefMassTime(int Anum, const char* elename, int _qref,int _laps_ref,double b_ref,double tofref,double tofref_err){

	double * masstoset = SearchAME(Anum,elename);
	if(masstoset[0] == -1){
		cout<<"Fail to set Ref mass!!!"<<endl;
		return;
	}
	else{
		ResetRef(masstoset[0],masstoset[1],_qref,_laps_ref,b_ref,tofref,tofref_err);
	}

	PrintInfo();

	
}

void ResetRefMassTime(string formula, int _qref,int _laps_ref,double b_ref,double tofref,double tofref_err){

	double * masstoset = SearchMolMass(formula);
	if(masstoset[0] == -1){
		cout<<"Fail to set Ref mass!!!"<<endl;
		return;
	}
	else{
		ResetRef(masstoset[0],masstoset[1],_qref,_laps_ref,b_ref,tofref,tofref_err);
	}

	PrintInfo();

}

void ResetX(int Peak1OrPeak2 , double tofx, double tofx_err){
	if(Peak1OrPeak2 !=1 && Peak1OrPeak2 !=2){cout<< "Peak number input err! [1 or 2]"<<endl; return;}

	if(tofx != -1) tof_x_cento[Peak1OrPeak2-1] = tofx;
	if(tofx_err != -1) tof_x_cento_err[Peak1OrPeak2-1] = tofx_err;
	PrintInfo();

}

//******************************************************************************************************
//*********** Set Marker  _lap_x=0 => same laps as laps_ref;  _lap_x=-1 => search in [laps_ref+-50]************
void MarkTof(double mass_xx , int _q_x, const char* IonName, int _lap_x ,int Tag,bool renew_current,bool renew_all){
	Int_t kala[]={633,808,799,417,433,600,617};
	static int Linenum = 0;
	static bool firstrunDone = false;
	//if(Pline == NULL){ Pline = new TLine[20];}
	//if(LabelXion == NULL){LabelXion = new TLatex[20];}
	if(TOFMarker::marker_reset){
		firstrunDone = false;
		renew_all	= true;
		TOFMarker::marker_reset=false;
	}
	if(renew_current== false && firstrunDone ==true ){
		 Linenum++;
		 Linenum = Linenum%40;
	}

	if(renew_all == true ){  // refresh and delete all line already drawen on canvas, 20 is maximum of array Pline and LabelXion
		for(int index=0;index<40;index++){
			marker_tof[index].Clear();
			
		}

		delete[] marker_tof;
		
		c1->Modified();

		marker_tof = new TOFMarker[40]; 

		Linenum=0;

	}

	marker_tof[Linenum].Clear();

	double moq_ref = (m_ref-q_ref*m_ele)/q_ref;
	double moq_x = (mass_xx-_q_x*m_ele)/_q_x;
	double tofx = 0;

	
	if(_lap_x==0){	// 0 => use same lap as ref
		tofx = TMath::Sqrt(moq_x/moq_ref) * (tof_ref_cento - t0) + t0;
		printf("tofx = %.4f @ %d laps\n",tofx,laps_ref);
		marker_tof[Linenum].laps = laps_ref;
	}    
	else if(_lap_x==-1){ // start to search laps for ion x
			for(_lap_x = laps_ref-50;_lap_x<=laps_ref+50;_lap_x++){
	 				tofx= TMath::Sqrt(moq_x/moq_ref) * (tof_ref_cento+ (_lap_x-laps_ref)* bref - t0) + t0;
	 				if(Tag==0){
	 					if(tofx>EJE0){
		 					 printf("tofx = %.4f @ %d laps\n",tofx,_lap_x);
		 					 marker_tof[Linenum].laps = _lap_x;
		 					 break;
	 					}
	 				}
	 				else if(Tag==1){
	 					if(tofx>EJE1){
		 					 printf("tofx = %.4f @ %d laps\n",tofx,_lap_x);
		 					 marker_tof[Linenum].laps = _lap_x;
		 					 break;
	 					}
	 				}
	 		}
	 		if(_lap_x > laps_ref+50){
	 			cout<<"Fail to Mark Tof at lap number ["<<laps_ref-50<<" , "<<laps_ref+50<<"]"<<endl;
	 			return;
	 		}

	}
	else{// designated lap number
		tofx = TMath::Sqrt(moq_x/moq_ref) * (tof_ref_cento+ (_lap_x-laps_ref)* bref - t0) + t0;
		printf("tofx = %.4f @ %d laps\n",tofx,_lap_x);
		marker_tof[Linenum].laps = _lap_x;
	}

	//printf("tofx = %.4f @ %d laps\n",tofx,_lap_x);

	string nametoset = IonName;
	if(_lap_x!=0){nametoset+="@"+ to_string(_lap_x);  /*cout<<nametoset.c_str()<<endl;*/}

	//printf("tofx = %.4f @ %d laps\n",tofx,laps_ref);

	double LineHigh =0;
	if(Tag ==0){
		LineHigh = h_xF->GetBinContent( h_xF->GetMaximumBin() ) * 1.01;
	}
	else{
		LineHigh = h_refF->GetBinContent( h_refF->GetMaximumBin() ) * 1.01;
	}

	marker_tof[Linenum].tag = Tag;
	marker_tof[Linenum].SetLine(tofx,LineHigh);
	marker_tof[Linenum].SetLabel(tofx,LineHigh,nametoset.c_str());
	marker_tof[Linenum].SetColor(kala[Linenum%7]);
	marker_tof[Linenum].SetMassAndCharge(mass_xx,_q_x);


    if(Tag==0){
    	c1->cd(1)->SetEditable(kTRUE);
    	marker_tof[Linenum].Draw();
    	c1->cd(1)->SetEditable(kFALSE);
	}
	else{
		c1->cd(3);
   		marker_tof[Linenum].Draw();
	}

	firstrunDone = true;
	TOFMarker::marker_amount = Linenum;
}

void MarkTof(int Anum,const char* EleName, int _q_x, int _lap_x ,int Tag, bool renew_current, bool renew_all){
		double* mass_xx = SearchAME(Anum,EleName);
		string nametoset = "^{";
		nametoset+= to_string(Anum)+"}";
		string temname= EleName;
		temname.erase(remove(temname.begin(),temname.end(),' '),temname.end());
		nametoset+=temname+"^{"+to_string(_q_x)+"+}"; 
		MarkTof(mass_xx[0],_q_x,nametoset.c_str(),_lap_x ,Tag,renew_current,renew_all);//EleName
}


void MarkTof(string formula, int _q_x, int _lap_x ,int Tag, bool renew_current, bool renew_all){
	Ion ion_tem(formula);
	string nametoset=ion_tem._name;
	nametoset+="^{"+ to_string(_q_x) +"+}";
	MarkTof(ion_tem.mass, _q_x, nametoset.c_str(), _lap_x,Tag, renew_current, renew_all);//ion_tem.name
}


//%%%%%%%%%%% manual add ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// must adjust ROI_index manually.
// "ROI_index" from 1 to 40; which showing the number of ROI available at the moment.
// However the default value of  "ROI_index" is 1;
// Therefore if it is the first ROI: SetROI(central tof, ROI width, [ROI_index or 1 ], Time offset, whether display counts, for which tag);
// Other cases: SetROI(central tof, ROI width, [++ROI_index], Time offset, whether display counts, for which tag);
// Interactive mode can only set ROI for tag0. To set ROI at tag1, using method SetROI() directly and specify at the arguement for tag value
// Other SetROI methods are similar, ROI_index must be increase manually.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void SetROI(double TCento, double ROI_width , int ROI_index , double TOffset, bool showcount,int Tag){

	if(ROI_index<1 || ROI_index>40) cout<<"Maximum ROI number is 40 !!!!"<<endl; 

	Int_t kala[]={633,808,799,417,433,600,617};
	auto LineHighSelection = [Tag]() ->double{
		if(Tag==0){ return h_xF->GetBinContent( h_xF->GetMaximumBin() ) * 1.01;}
		else if(Tag==1){return h_refF->GetBinContent( h_refF->GetMaximumBin() ) * 1.01;}
		else return -10;
	};
	double LineHigh =LineHighSelection();
	if(LineHigh<0){cout<<"Tag input is wrong!! Try again!"<<endl; return;}

	if(LineHigh<1)LineHigh=1; //in case zoom into empty region

	double TextHigh = LineHigh*0.5;

	double LeftEdge =  TCento + TOffset - ROI_width;
	double RightEdge = TCento + TOffset + ROI_width;

	ROI_tag[ROI_index-1] = Tag;
	ROIL[ROI_index-1].Clear();
	ROIL[ROI_index-1].SetX1(LeftEdge);
	ROIL[ROI_index-1].SetX2(LeftEdge);
	ROIL[ROI_index-1].SetY1(0.);
	ROIL[ROI_index-1].SetY2(LineHigh);
	ROIL[ROI_index-1].SetLineColor(kala[(ROI_index-1)%7]);
	ROIL[ROI_index-1].SetLineStyle(7);

	ROIR[ROI_index-1].Clear();
	ROIR[ROI_index-1].SetX1(RightEdge);
	ROIR[ROI_index-1].SetX2(RightEdge);
	ROIR[ROI_index-1].SetY1(0.);
	ROIR[ROI_index-1].SetY2(LineHigh);
	ROIR[ROI_index-1].SetLineColor(kala[(ROI_index-1)%7]);
	ROIR[ROI_index-1].SetLineStyle(7);

	ROI_label[ROI_index-1].Clear();
	ROI_label[ROI_index-1].SetTextAlign(21);
    ROI_label[ROI_index-1].SetTextSize(0.07);
    ROI_label[ROI_index-1].SetText(TCento + TOffset,TextHigh,Form("%d",ROI_index));
    ROI_label[ROI_index-1].SetTextColor(kala[(ROI_index-1)%7]);

	c1->Update();
	c1->Modified();

	if(Tag==0){
		c1->cd(1)->SetEditable(kTRUE);
	}
	else{
		c1->cd(3);
	}

		ROIL[ROI_index-1].Draw();
		ROIR[ROI_index-1].Draw();
		ROI_label[ROI_index-1].Draw();
		c1->Modified();
		c1->cd(1)->SetEditable(kFALSE);

		if(showcount){
			cout<<"Counts in ROI__"<<ROI_index<<" = ";
			if(Tag==0) cout<<h_xF->Integral(h_xF->FindBin(LeftEdge),h_xF->FindBin(RightEdge));
			else cout<<h_refF->Integral(h_refF->FindBin(LeftEdge),h_refF->FindBin(RightEdge));
			cout<<endl;
		}
	
	
}

void SetMassROI(double _massxx,int _q_x, int _laps_x,int ROI_index, double TOffset, double ROI_width, bool showcount,int Tag){

	double mass_xx = _massxx;

	double moq_ref = (m_ref-q_ref*m_ele)/q_ref;
	double moq_x = (mass_xx-_q_x*m_ele)/_q_x;
	double tofx = TMath::Sqrt(moq_x/moq_ref) * (tof_ref_cento +(_laps_x - laps_ref)*bref- t0) + t0;

	if(_laps_x==0){// same _laps_x as laps_ref;
		tofx = TMath::Sqrt(moq_x/moq_ref) * (tof_ref_cento - t0) + t0;
	}
	else if(_laps_x==-1){// auto search _laps_x
		for(_laps_x=laps_ref-50;_laps_x<=laps_ref+50;_laps_x++){
			tofx = TMath::Sqrt(moq_x/moq_ref) * (tof_ref_cento +(_laps_x - laps_ref)*bref- t0) + t0;
			if(tofx>EJE0){break;}
		}

		if(_laps_x>(laps_ref+50)){cout<<"No matching in laps_x = ["<<laps_ref-50<<" , "<<laps_ref+50<<"]"<<endl; return;}

	}
	else{ // designate _laps_x
			tofx = TMath::Sqrt(moq_x/moq_ref) * (tof_ref_cento +(_laps_x - laps_ref)*bref- t0) + t0;
	}



	SetROI(tofx, ROI_width, ROI_index, TOffset, showcount, Tag);


}

void SetEleROI(int Anum,const char *element, int _q_x, int _laps_x, int ROI_index, double TOffset, double ROI_width, bool showcount,int Tag){
	double * mass_xx = SearchAME(Anum,element);

	if(mass_xx[0]==-1){cout<<"fail to get mass from AME!!!"<<endl; return;}
	SetMassROI(mass_xx[0], _q_x, _laps_x, ROI_index,  TOffset,  ROI_width, showcount, Tag);	

}


void SetMoleculeROI(string formula, int _q_x, int _laps_x, int ROI_index, double TOffset, double ROI_width, bool showcount,int Tag){
	double * mass_xx = SearchMolMass(formula);

	if(mass_xx[0]<0){cout<<"fail to get molecule mass from AME!!!"<<endl; return;}
	SetMassROI(mass_xx[0], _q_x, _laps_x, ROI_index,  TOffset,  ROI_width, showcount, Tag);	

}

void CreateComponent(int Compon_index,int Anum,const char* elename,int N_atom, bool addTrue){

	if(Compon_index<1 || Compon_index>20) {cout<<"Compon_index err!!!!!!!! [1-20]"<<endl; return;}

	// clear history
	if(addTrue == false){ 
		component_mass[Compon_index-1] =0; 
		component_mass_err[Compon_index-1]=0; 
		for(int j=0;j<20;j++){component_name[Compon_index-1][j]='\0';}
	}

	double* mass_tem = SearchAME(Anum,elename);

	if(mass_tem[0]==-1){cout<<"fail to get mass from AME!!!"<<endl; return;}

	mass_tem[0] = mass_tem[0] * N_atom; // mean of N atom
	mass_tem[1] = mass_tem[1] * N_atom; // mean err of N atom

	component_mass[Compon_index-1] = component_mass[Compon_index-1] + mass_tem[0];
	component_mass_err[Compon_index-1] = component_mass_err[Compon_index-1] * component_mass_err[Compon_index-1] + mass_tem[1] * mass_tem[1];
	component_mass_err[Compon_index-1] = TMath::Sqrt(component_mass_err[Compon_index-1]);

	char Nnum_c[5];
	if(component_name[Compon_index-1][0] != '\0'){
		sprintf(Nnum_c,":%d",N_atom);
	}
	else sprintf(Nnum_c,"%d",N_atom);
	strcat(&component_name[Compon_index-1][0],Nnum_c);
	strcat(&component_name[Compon_index-1][0],elename);

	char Anum_c[5];
	sprintf(Anum_c,"%d",Anum);
	strcat(&component_name[Compon_index-1][0],Anum_c);
	
	
}

void PrintComponentList(){
	cout<<"component list:"<<endl;
	for(int j=0;j<20;j++){
		if(component_mass[j] == 0) continue;
		printf("%d.\tname  %s  mass:  %.4f(%.4f)\n",j+1,&component_name[j][0],component_mass[j],component_mass_err[j]);

	}

}

void PrintMassTimeList(){

	for(int j=0;j<20;j++){
		if(mass_recorded[j]!= 0) {
			printf("%d.\tname  '%s'  mass:  %.4f(%.4f)\t",j+1,mass_recorded_message[j],mass_recorded[j],mass_recorded_err[j]);
			if(time_recorded[j]!=0) printf("name  %s   time:  %.4f(%.4f)\n",time_recorded_message[j],time_recorded[j],time_recorded_err[j]);
			else printf("name  N/A   time:  N/A(N/A)\n");
		}
		else {
			if(time_recorded[j]!=0){
				printf("%d.\tname  N/A  mass:  N/A(N/A)\t",j+1);
				printf("name  '%s'   time:  %.4f(%.4f)\n",time_recorded_message[j],time_recorded[j],time_recorded_err[j]);
			}
			else continue;
		}

	}

}

void SaveMassResult(int result_index, const char * Note, double mass_result, double mass_result_err){

		if(result_index<1 || result_index>20){cout<<"Result_index err!!!!!!!! [1-20]"<<endl; return;}
		if(!strncmp(Note,"",3)){cout<<"err...Input note!!!!"<<endl; return;}
		if(mass_result != -1){
			mass_recorded[result_index-1] = mass_result;
			mass_recorded_err[result_index-1] = mass_result_err;
			mass_recorded_message[result_index-1]= Note;
		}

		PrintMassTimeList();
		return;

}


void SaveTimeResult(int result_index,const char * Note, double time_result, double time_result_err){

		if(result_index<1 || result_index>20){cout<<"Result_index err!!!!!!!! [1-20]"<<endl; return;}
		if(!strncmp(Note,"",3)){cout<<"err...Input note!!!!"<<endl; return;}
		if(time_result !=-1){
			time_recorded[result_index-1] = time_result;
			time_recorded_err[result_index-1] = time_result_err;
			time_recorded_message[result_index-1]= Note;
		}
		PrintMassTimeList();
		return;

}


double* ExtractMassDeviate(double mass_A, double mass_err_A,double mass_B, double mass_err_B,bool saveMass,int result_index, const char * Note){  // mass A - mass B
		static double mass_dev[2];
		mass_dev[0]=0;
		mass_dev[1]=0;

		double mass_dev_tem = mass_A - mass_B;
		double mass_dev_err_tem = mass_err_A * mass_err_A + mass_err_B * mass_err_B;
		mass_dev_err_tem = TMath::Sqrt(mass_dev_err_tem);
		printf("mass deviate: %.4f(%.4f)\n",mass_dev_tem,mass_dev_err_tem);

		mass_dev[0] = mass_dev_tem;
		mass_dev[1] = mass_dev_err_tem;

		if(saveMass){
			if(!strncmp(Note,"",3)){cout<<"Input message to record mass"<<endl; return mass_dev;}
			SaveMassResult(result_index,Note,mass_dev_tem,mass_dev_err_tem);
		}

		return mass_dev;

}

double* TokeV90(double mass_value,double mass_value_err, bool verbal=true){
		static double massreturn[2];
		massreturn[0]=0; massreturn[1]=0;
		double masskev = mass_value * KEV90;
		double masskev_err = (mass_value*err_KEV90)*(mass_value*err_KEV90) + (KEV90*mass_value_err)*(KEV90*mass_value_err);
		masskev_err = TMath::Sqrt(masskev_err);
		if(verbal)printf("mass in keV90: %.5f(%.5f)keV90\n",masskev,masskev_err);
		massreturn[0] = masskev;
		massreturn[1] = masskev_err;
		return massreturn;
}

double* MassExcess(double mass_value,double mass_value_err, bool verbal=true){
		int A = TMath::Nint(mass_value*1e-6);
		double mass_excess = mass_value - A*1e6;  // unit in micro u
		return TokeV90(mass_excess,mass_value_err,verbal);
		

}


vector<double> PlotMultiMassResult(int NumResult,double* multimass, double* multimass_err, double massAME, double massAME_err,double YMaximum=-1){
		bool keVon=false;
		if(massAME < 931500) keVon = true;

		double* weight_i = new double[NumResult];
		double sum_weight = 0;
		double mass_mean = 0;
		double* N_xaxis = new double[NumResult];

		for(int i=0;i<NumResult;i++){
			weight_i[i]= 1/(multimass_err[i] * multimass_err[i]);
			mass_mean = mass_mean + weight_i[i] * multimass[i];
			sum_weight = sum_weight + weight_i[i];
			N_xaxis[i] = i+1;
		}

		mass_mean = mass_mean / sum_weight;

		double mass_mean_err = TMath::Sqrt(1 / sum_weight);

		//&&&&&&&&&&&&&&  calculate birge ratio  &&&&&&&&&&&&&&&&&&&&&&
		double birge_ratio=0;
		double mass_mean_err_B=-1;
		double mass_mean_err_MB=-1;
		double t_qth=-1;
		if(NumResult>1){
			for(int i=0;i<NumResult;i++){
				birge_ratio = birge_ratio + TMath::Power((multimass[i] - mass_mean)/multimass_err[i],2);
			}

			birge_ratio = TMath::Sqrt(birge_ratio/(NumResult-1));  // result of birge ratio
			// &&&&&&&&&&&&&&& end of calculation of birge ratio  &&&&&&&&&&&&&7

			mass_mean_err_B = mass_mean_err * birge_ratio;  // mean error adjusted by birge ratio: err(mean) * birge ratio

			if(NumResult>3){
				 mass_mean_err_MB = TMath::Sqrt((NumResult-1.)/(NumResult-3.)) * mass_mean_err_B; // by modified birge ratio
				 t_qth = TMath::StudentQuantile(0.975,NumResult-1);
				 
			}

		}

		// output format:
		// mass mean xx(xx)[amu]; mass mean xx(xx)[keV]; mass dva: xx(xx)[amu]; mass dav: xx(xx)[keV]; birge ; modified birge err; 0.975 student distribute region
		static vector<double> outputdata;
		outputdata.clear();
		double* tem_convert;

		if(MassGraph != NULL){  // delete mass graph and pad
			// important: must delete pad and text objects before clear and delete TGraph; 
			// otherwise will cause crash !!!!!!!!!!!!!!!!
			// it seems pad and text object drawn on TGraph also belong to the graph. Clear graph will also delete object on it !!!!!!!!!!!!
			delete padmass;
			delete padbirge;
			delete text_delta_mass;
			delete text_mass_mean;
			delete text_birge;
			MassGraph->Clear();
			MassGraph->Delete();
		}

		MassGraph = new TGraphErrors();

		int i =0; // just for looping purpose

		for(i=0;i<NumResult;i++){
			if(keVon){
				tem_convert = MassExcess(multimass[i],multimass_err[i],false);
				MassGraph->SetPoint(i,N_xaxis[i],tem_convert[0]-massAME);
				MassGraph->SetPointError(i,0,tem_convert[1]);	
			}
			else{
				MassGraph->SetPoint(i,N_xaxis[i],multimass[i]-massAME);
				MassGraph->SetPointError(i,0,multimass_err[i]);	
			}		
		}
		

		double mass_dev = 0;

		if(keVon){  // change to keV unit
			tem_convert = MassExcess(mass_mean,mass_mean_err,false);
			mass_dev = tem_convert[0] - massAME;
		}
		else{
			mass_dev = mass_mean-massAME;

		}


		//MassGraph->SetPoint(i,i+1.0,mass_dev);
		//MassGraph->SetPointError(i,0,mass_mean_err);
		MassGraph->GetXaxis()->SetLimits(0,i+2.0);
		MassGraph->GetYaxis()->SetTitleOffset(1.4);
		MassGraph->GetYaxis()->CenterTitle();
   		//MassGraph->GetYaxis()->SetTitle("#left| #frac{1}{1 - #Delta#alpha}#right|^{2} (1+cos^{2}#theta)");
   		if(keVon) MassGraph->GetYaxis()->SetTitle("#Delta ME = ME - ME_{AME2020} [keV]");
   		else MassGraph->GetYaxis()->SetTitle("#Delta m = mass - mass_{AME2020} [micro amu]");
   		MassGraph->SetMarkerStyle(21);
   		if(YMaximum != -1){MassGraph->SetMaximum(YMaximum); MassGraph->SetMinimum(-YMaximum);}
		
		massLowEdge.Clear();
		massLowEdge.SetX1(0);
		massLowEdge.SetX2(i+2.0);
		if(keVon){
			massLowEdge.SetY1(mass_dev-tem_convert[1]);
			massLowEdge.SetY2(mass_dev-tem_convert[1]);
		}
		else{
			massLowEdge.SetY1(mass_dev-mass_mean_err);
			massLowEdge.SetY2(mass_dev-mass_mean_err);
		}
		massLowEdge.SetLineColor(kRed);
		massLowEdge.SetLineStyle(7);

		massHighEdge.Clear();
		massHighEdge.SetX1(0);
		massHighEdge.SetX2(i+2.0);
		if(keVon){
			massHighEdge.SetY1(mass_dev+tem_convert[1]);
			massHighEdge.SetY2(mass_dev+tem_convert[1]);
		}
		else{
			massHighEdge.SetY1(mass_dev+mass_mean_err);
			massHighEdge.SetY2(mass_dev+mass_mean_err);
		}
		massHighEdge.SetLineColor(kRed);
		massHighEdge.SetLineStyle(7);

		ameLowEdge.Clear();
		ameLowEdge.SetX1(0);
		ameLowEdge.SetX2(i+2.0);
		ameLowEdge.SetY1(-massAME_err);
		ameLowEdge.SetY2(-massAME_err);
		ameLowEdge.SetLineColor(kBlue);
		ameLowEdge.SetLineStyle(7);

		ameHighEdge.Clear();
		ameHighEdge.SetX1(0);
		ameHighEdge.SetX2(i+2.0);
		ameHighEdge.SetY1(massAME_err);
		ameHighEdge.SetY2(massAME_err);
		ameHighEdge.SetLineColor(kBlue);
		ameHighEdge.SetLineStyle(7);

		//c2->Update();
		//c2->Modified();

		if(c2==NULL)c2 = new TCanvas("c2","Mass Result preview",800,600);
		c2->cd();
		c2->Update();
		c2->Modified();
		MassGraph->Draw("AP");
		massLowEdge.Draw();
		massHighEdge.Draw();
		ameLowEdge.Draw();
		ameHighEdge.Draw();

		double mass_mean_amu=0 , mass_mean_err_amu=0;
		double mass_mean_kev=0 , mass_mean_err_kev=0;
		double mass_dev_amu=0 , mass_dev_err_amu=0;
		double mass_dev_kev=0 , mass_dev_err_kev=0;
		double mass_ref_amu=0 , mass_ref_err_amu=0;
		double mass_ref_kev=0 , mass_ref_err_kev=0;
		double* mass_mean_dev;

			mass_mean_amu = mass_mean; mass_mean_err_amu = mass_mean_err;
			tem_convert = MassExcess(mass_mean,mass_mean_err,false);	
			mass_mean_kev = tem_convert[0]; mass_mean_err_kev = tem_convert[1];

		//	tem_convert = MassExcess(mass_mean,mass_mean_err,false);

		if(keVon){ // ref in keV
			mass_ref_kev = massAME ; mass_ref_err_kev = massAME_err;
			mass_ref_amu =-1 ; mass_ref_err_amu =-1; // no information

			mass_mean_dev = ExtractMassDeviate(mass_mean_kev,mass_mean_err_kev,mass_ref_kev, mass_ref_err_kev);
			mass_dev_kev = mass_mean_dev[0]; mass_dev_err_kev = mass_mean_dev[1];
		
			mass_dev_amu =-1 ; mass_dev_err_amu = -1; // no information

		}
		else{ // ref in amu
			mass_ref_amu = massAME ; mass_ref_err_amu =  massAME_err;
			tem_convert = MassExcess(massAME,massAME_err,false);
			mass_ref_kev = tem_convert[0]; 	mass_ref_err_kev = tem_convert[1];
			mass_mean_dev = ExtractMassDeviate(mass_mean_amu,mass_mean_err_amu,mass_ref_amu,mass_ref_err_amu);
			mass_dev_amu = mass_mean_dev[0] ; mass_dev_err_amu = mass_mean_dev[1]; // deviate in amu

			mass_mean_dev = ExtractMassDeviate(mass_mean_kev,mass_mean_err_kev,mass_ref_kev,mass_ref_err_kev);
			mass_dev_kev = mass_mean_dev[0]; mass_dev_err_kev = mass_mean_dev[1];

		}

		printf("Mass Mean= %.4f(%.4f)[micro_amu]; %.4f(%.4f)[keV]\n",mass_mean_amu,mass_mean_err_amu,mass_mean_kev,mass_mean_err_kev);
		printf("Mass_mean - Mass_ame\n");
		printf("Mass deviate: %.4f(%.4f)[micro_amu]; %.4f(%.4f)[keV]\n",mass_dev_amu,mass_dev_err_amu,mass_dev_kev,mass_dev_err_kev);
		printf("Birge ratio: %.4f\n",birge_ratio);
		printf("Birge Adjusted Mass mean: %.4f(%.4f)\n",mass_mean,mass_mean_err_B);
		printf("Modigied Birge Mass mean: %.4f(%.4f); 0.95 coverage interval[%.4f-%.4f,%.4f+%.4f]\n",mass_mean,mass_mean_err_MB,mass_mean,t_qth*mass_mean_err_MB,mass_mean,t_qth*mass_mean_err_MB);

		outputdata.push_back(mass_mean);outputdata.push_back(mass_mean_err);  // mass mean [amu]
		outputdata.push_back(mass_mean_kev);outputdata.push_back(mass_mean_err_kev); // mass mean [kev]
		outputdata.push_back(mass_dev_amu);outputdata.push_back(mass_dev_err_amu); // mass dva [amu]
		outputdata.push_back(mass_dev_kev);outputdata.push_back(mass_dev_err_kev); // mass dva[keV]
		outputdata.push_back(birge_ratio);  // err adjust => birge * mean_err;
		outputdata.push_back(mass_mean_err_MB);  // => modified birge mass mean err = mass_mean_err_MB
		outputdata.push_back(t_qth);                // => 0.95 credible interval +/ t_qth*mass_mean_err_MB

		padmass = new TPad("padmass","padmass",0.1,0.,0.9,0.06,0,0);
		padbirge = new TPad("birge","birge",0.7,0.82,0.9,0.89,0,0);
		text_delta_mass = new TLatex(); // text setting for delta mass
		text_delta_mass->SetTextSize(0.60);
		text_delta_mass->SetTextAlign(32);
		text_mass_mean = new TLatex();  // text setting for mean mass;
		text_mass_mean->SetTextSize(0.60);
		text_mass_mean->SetTextAlign(12);
		text_birge = new TLatex();
		text_birge->SetTextSize(0.60);
		text_birge->SetTextAlign(32);

		if(keVon){
			text_delta_mass->SetText(1.,0.5,Form("#Delta ME=%.4f(%.4f) keV",mass_dev_kev,mass_dev_err_kev));
			text_mass_mean->SetText(0.00,0.5,Form("ME=%.4f(%.4f) keV",mass_mean_kev,mass_mean_err_kev));
		}
		else{
			text_delta_mass->SetText(1.,0.5,Form("#Delta m=%.4f(%.4f) micro_amu",mass_dev_amu,mass_dev_err_amu));
			text_mass_mean->SetText(0.00,0.5,Form("m=%.4f(%.4f) micro_amu",mass_mean,mass_mean_err));
		}
		text_birge->SetText(0.99,0.5,Form("birge:%.4f",birge_ratio));

		padmass->Draw();
		padmass->cd();
		text_delta_mass->Draw();
		text_mass_mean->Draw();

		c2->cd();
		padbirge->Draw();
		padbirge->cd();
		text_birge->Draw();

		c2->cd();
		c2->Modified();
		c2->Update();

		delete[] weight_i;
		delete[] N_xaxis;

		return outputdata;

}


void savefile_current_fitparas(const char* Note, bool NewIndex){
	FILE *fp;
	time_t now = time(0);
	char * CTime = ctime(&now);
	//double mainintensity=0;
	static int TFIndex=1;
	if(NewIndex == true) TFIndex =1;  // Reset TFIndex to 1
    char fname[200];
	sprintf(fname,"%shito_fitpara_%s.txt",FilePath.c_str(),FileName.c_str());

    fp = fopen(fname,"r");
	if(fp==NULL){ // if file doesn't exist. make a new one and give the first line
		cout<< "Create a new file: "<<"hito_fitpara_"<<FileName<<".txt"<<endl;
		fp = fopen(fname,"a+");
		fprintf(fp,"Index \t Note \t FitFuncGroup \t NumOfPeaks \t MainPeakIndex \t ParasLock \t High \t Center \t Sigma \t tctailR \t tctailL \t whatkindofBg \t offset \t slope") ;
		fprintf(fp,"\t fitRangeL \t fitRangeR \t tag \t  histoNBins \t histoRangeL \t histoRangeR \t Date\n");
	}
	else{// when file is existed. add data to the end;
		fclose(fp);
		fp = fopen(fname,"a+");
	}


	double tem_high_file=0, tem_center_file = 0, tem_sigma_file=0, tem_tctailR = 0, tem_tctailL = 0;
	double tem_offset_file=0, tem_slop_file=0;
	int peak_para_offset=0;

	if(BackGroundCurve==0){
		tem_offset_file = 0;
		tem_slop_file = 0;
	}
	else{
		//peak_para_offset = (4+ (SetOneorTwo-1) ) * NumOfPeaks; // SetOneOrTwe = 1 => 4 paras each peak; =2 =>5 paras each peak
		peak_para_offset = Paras_offset_cal(NumOfPeaks,SetOneorTwo,MainPeakIndex,ParasLock);
		tem_offset_file = tem_func->GetParameter(0 + peak_para_offset);
		tem_slop_file = tem_func->GetParameter(1 + peak_para_offset);

	}

	
	int NumOfPars_all = (SetOneorTwo+3)*NumOfPeaks;
	if(BackGroundCurve !=0)NumOfPars_all+=2;

	double * setpar = new double[NumOfPars_all];

	//distribute parameter input from outside; becasue sub peaks share same parameters with main peaks (sigma, taiL_1_Left, tail_1_R, tail_2_R)!!!
	ParaDistributor(SetOneorTwo,NumOfPeaks,BackGroundCurve,MainPeakIndex,ParasLock,tem_func->GetParameters(),setpar);


      for(int index_p=0;index_p<NumOfPeaks;index_p++){
             peak_para_offset = (4+ (SetOneorTwo-1) ) * index_p;
             tem_high_file = setpar[0+peak_para_offset]; 
             tem_center_file = setpar[1+peak_para_offset]; 
             tem_sigma_file = setpar[2+peak_para_offset];
             tem_tctailR = 0;
             tem_tctailL = 0;
             if(SetOneorTwo !=0)tem_tctailR = setpar[3+peak_para_offset];
             if(SetOneorTwo ==2)tem_tctailL = setpar[4+peak_para_offset];
             fprintf(fp,"%d \t %s \t %d \t %d \t %d \t %d \t %.1f \t %.3f \t %.3f \t %.3f \t %.3f \t %d \t %.3f \t %.3f",TFIndex,Note,SetOneorTwo,NumOfPeaks,MainPeakIndex,(int)ParasLock,tem_high_file,tem_center_file,tem_sigma_file,tem_tctailR,tem_tctailL,BackGroundCurve,tem_offset_file,tem_slop_file);
             fprintf(fp," \t %.3f \t %.3f \t %d \t %d \t %.2f \t %.2f \t %s",fitrangeL,fitrangeR,active_tagX ,active_histoX_Nbins, active_histoX_RangeL , active_histoX_RangeR, CTime);
             
       }


 	TFIndex++;
 	delete[] setpar;

	fclose(fp);
}


void savefile_history_fitparas(DefGeneralFunc * history, const char *Note,bool NewIndex){
	// back current setting
	int SetOneorTwo_backup = SetOneorTwo;
	int NumOfPeaks_backup = NumOfPeaks;
	int BackGroundCurve_backup = BackGroundCurve;
	int MainPeakIndex_backup = MainPeakIndex;
	bool ParasLock_backup = ParasLock;
	double fitrangeL_backup = fitrangeL;
	double fitrangeR_backup = fitrangeR;
	double active_tagX_backup = active_tagX;
	double active_histoX_Nbins_backup = active_histoX_Nbins;
	double active_histoX_RangeL_backup = active_histoX_RangeL;
	double active_histoX_RangeR_backup = active_histoX_RangeR;

	SetOneorTwo	= history->SetOneorTwo;   // 0 = pure gaus;  1 = func1 (gauss + exp Right tail) ; 2 = func2 (gauss + exp Ritht and Left tail)
	NumOfPeaks	= history->NumOfPeaks;    // at least 1 peak
	BackGroundCurve	= history->BackGroundCurve; 
	MainPeakIndex = history->MainPeakIndex;
	ParasLock = history->ParasLock;
	active_tagX	= history->TAG;
	fitrangeL	= history->fitRangeL;
	fitrangeR	= history->fitRangeR;
	active_histoX_Nbins	= history->histo_nbins;
	active_histoX_RangeL	= history->histo_L;
	active_histoX_RangeR	= history->histo_R;

	TF1 * tem_func_backup = tem_func;

	tem_func = history->GetTFunc();

	const char* inputNote;
	if(Note != NULL) inputNote = Note;
	else inputNote = history->name;

	savefile_current_fitparas(inputNote,NewIndex);

	// restore current setting
	SetOneorTwo = SetOneorTwo_backup;
	NumOfPeaks = NumOfPeaks_backup;
	BackGroundCurve = BackGroundCurve_backup;
	MainPeakIndex = MainPeakIndex_backup;
	ParasLock = ParasLock_backup;
	fitrangeL = fitrangeL_backup;
	fitrangeR = fitrangeR_backup;
	active_tagX = active_tagX_backup;
	active_histoX_Nbins = active_histoX_Nbins_backup;
	active_histoX_RangeL = active_histoX_RangeL_backup;
	active_histoX_RangeR = active_histoX_RangeR_backup;
	tem_func = tem_func_backup;		

}


void readfile_fitparas_history(int recordIndex_low=1,int recordIndex_high=0){//recordIndex_high=0 => readout all // only readout record between low and high
	// index begin from 1;
	bool shouldreadout=false;
	int recordIndex_current=0;  // default => Index of  the first history is 1 !!!!!!!!!!!!!!!!!!!;
	FILE *fp;
    char fname[200];
	sprintf(fname,"%shito_fitpara_%s.txt",FilePath.c_str(),FileName.c_str());
	fp = fopen(fname,"r");
	if(fp==NULL){

			cout<<"\e[1;33m"<<"Error!! not history file of fit paras existed!!!"<<"\e[0m"<<endl;
			return;
	}
	char filehead[1000]; // just to get file head;
	fgets(filehead,1000,(FILE*)fp);
	static int history_num=0;
	history_num=0;  // only provide container for 30 history records;

	int TFIndex;
	char Note[100];
	int _SetOneorTwo;
	int _NumOfPeaks;
	int _BackGroundCurve;
	int _MainPeakIndex;
	bool _ParasLock;
	int _ParasLock_int;
	double tem_high_file=0, tem_center_file = 0, tem_sigma_file=0, tem_tctailR = 0, tem_tctailL = 0;
	double tem_offset_file=0, tem_slop_file=0;
	double fitrangeL=0,fitrangeR=0;
	int _tagX=0;
	int _histoX_Nbins=0; 
	double _histoX_RangeL=0 , _histoX_RangeR=0; 
	char CTime[100];
	
	bool firstLoadDone=false;
	double* tem_paras=NULL;

	while(!feof(fp)){

		if(firstLoadDone){
			printf("%d \t %s \t %d \t %d \t %d \t %d\n", TFIndex,Note,_SetOneorTwo,_NumOfPeaks,_MainPeakIndex,_ParasLock);
			printf(" \t %.1f \t %.3f \t %.3f \t %.3f \t %.3f \t %d \t %.3f \t %.3f\n",tem_high_file,tem_center_file,tem_sigma_file,tem_tctailR,tem_tctailL,_BackGroundCurve,tem_offset_file,tem_slop_file);
			printf("\t %.3f \t %.3f \t %d \t %d \t %.2f \t %.2f \t %s",fitrangeL,fitrangeR,_tagX ,_histoX_Nbins, _histoX_RangeL ,_histoX_RangeR, CTime);

			int NumOfPars = (_SetOneorTwo+3)*_NumOfPeaks;
			if(_BackGroundCurve !=0)NumOfPars+=2;
			tem_paras = new double[NumOfPars];
			int Par_index=0;
			int loopnumber = _NumOfPeaks;
			do{
				tem_paras[Par_index++]=tem_high_file; tem_paras[Par_index++]=tem_center_file;tem_paras[Par_index++]=tem_sigma_file;
				if(_SetOneorTwo !=0)tem_paras[Par_index++] = tem_tctailR;
				if(_SetOneorTwo==2)tem_paras[Par_index++] = tem_tctailL;
				if((--loopnumber) ==0){ // last subline for multipeak curve
					if(_BackGroundCurve !=0){tem_paras[Par_index++] = tem_offset_file;tem_paras[Par_index++] =tem_slop_file;}
				}
				else{ // cout<<"crash in else"<<endl; // debug
					fscanf(fp,"%d \t %s \t %d \t %d \t %d \t %d", &TFIndex,Note,&_SetOneorTwo,&_NumOfPeaks,&_MainPeakIndex,&_ParasLock_int);
					fscanf(fp," \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \t %lf \t %lf",&tem_high_file,&tem_center_file,&tem_sigma_file,&tem_tctailR,&tem_tctailL,&_BackGroundCurve,&tem_offset_file,&tem_slop_file);
					fscanf(fp," \t %lf \t %lf \t %d \t %d \t %lf \t %lf",&fitrangeL,&fitrangeR,&_tagX ,&_histoX_Nbins, &_histoX_RangeL ,&_histoX_RangeR);
					fgets(CTime,100,(FILE*)fp); //cout<<"crash3"<<endl;  // debug
					_ParasLock = (bool)_ParasLock_int; // 
			// debug purpose only
			/*printf("%d \t %s \t %d \t %d \t %d \t %d\n", TFIndex,Note,_SetOneorTwo,_NumOfPeaks,_MainPeakIndex,_ParasLock);
			printf(" \t %.1f \t %.3f \t %.3f \t %.3f \t %.3f \t %d \t %.3f \t %.3f\n",tem_high_file,tem_center_file,tem_sigma_file,tem_tctailR,tem_tctailL,_BackGroundCurve,tem_offset_file,tem_slop_file);
			printf("\t %.3f \t %.3f \t %d \t %d \t %.2f \t %.2f \t %s",fitrangeL,fitrangeR,_tagX ,_histoX_Nbins, _histoX_RangeL ,_histoX_RangeR, CTime);*/


				}
			}while(loopnumber);
			//cout<<"after while"<<endl; // debug
			if(shouldreadout){
				if(FHistory[history_num%30] !=NULL){cout<<"delete history"<<endl; delete FHistory[history_num%30];}
				int Numberofpars_TF = NumOfPars_TF_cal(_SetOneorTwo ,_NumOfPeaks,_BackGroundCurve,_ParasLock);
				double* setpar = new double[Numberofpars_TF];
				ParaContractor(_SetOneorTwo,_NumOfPeaks,_BackGroundCurve,_MainPeakIndex,_ParasLock,tem_paras,setpar);
				//for(int j=0;j<Numberofpars_TF;j++){cout<<setpar[j]<<endl;} // debug
				FHistory[(history_num++)%30]=new DefGeneralFunc(Note,_tagX,_SetOneorTwo,_NumOfPeaks,_BackGroundCurve,_MainPeakIndex,_ParasLock,setpar,fitrangeL,fitrangeR,_histoX_Nbins,_histoX_RangeL,_histoX_RangeR,CTime);
				delete[] setpar;
				if(tem_paras != NULL){cout<<"delete par"<<endl; delete[] tem_paras;} // clear for next aound
				for(int i=0;i<100;i++){Note[i]='\0'; CTime[i] = '\0';}

				shouldreadout=false;
			}
		}

		fscanf(fp,"%d \t %s \t %d \t %d \t %d \t %d", &TFIndex,Note,&_SetOneorTwo,&_NumOfPeaks,&_MainPeakIndex,&_ParasLock_int);
		fscanf(fp," \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \t %lf \t %lf",&tem_high_file,&tem_center_file,&tem_sigma_file,&tem_tctailR,&tem_tctailL,&_BackGroundCurve,&tem_offset_file,&tem_slop_file);
		fscanf(fp," \t %lf \t %lf \t %d \t %d \t %lf \t %lf",&fitrangeL,&fitrangeR,&_tagX ,&_histoX_Nbins, &_histoX_RangeL ,&_histoX_RangeR);
		fgets(CTime,100,(FILE*)fp);
		_ParasLock = (bool)_ParasLock_int;

		recordIndex_current++;
		if(recordIndex_high==0 || (recordIndex_current>=recordIndex_low && recordIndex_current<=recordIndex_high)) shouldreadout=true;
	
		firstLoadDone=true;
	}// while loop over file

	fclose(fp);

}


void fitlooping(Option_t *Fitoption="MEQ",int loopNum=50, double newRangeL=0,double newRangeR=0, bool fixpars=true){// true => keep default setting; 
	// "fixpars" can not be use to FixParameter;  but it can only set "ReleaseParameters" by set to false !!!!!!!!!
	if(newRangeL>0 && newRangeR>0){fitrangeL = newRangeL; fitrangeR = newRangeR;}

	int peak_para_offset=0;
	if(!fixpars){
		if(ParasLock){// only mainpeak has parameters for sigma and tails
			peak_para_offset = Paras_offset_cal(MainPeakIndex-1,SetOneorTwo,MainPeakIndex,ParasLock);
			for(int i=0;i<=SetOneorTwo;i++){
				tem_func->ReleaseParameter(peak_para_offset+2+i);
			}
		}
		else{
			for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
				peak_para_offset = Paras_offset_cal(PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock);
				for(int i=0;i<=SetOneorTwo;i++){
					tem_func->ReleaseParameter(peak_para_offset+2+i);
				}
			}
		}
	}

	c1->cd(2)->SetEditable(kTRUE);
	printf("width for fit: %.1f [ns]\n",fitrangeR-fitrangeL);
	if(tem_func!=NULL){
		for(int i=0;i<loopNum;i++){
			h_zoom_x->Fit(tem_func,Fitoption,"",fitrangeL,fitrangeR);
		}
	}
	else{
		cout<<"fit function is empty!!! break"<<endl;
		return;
	}

	for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
			peak_para_offset = Paras_offset_cal(PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock);
				//peak_para_offset = PeakIndex * (4+(SetOneorTwo-1));
			printf("centro_%d: %.4f(%.4f) \n",PeakIndex+1,tem_func->GetParameter(1+peak_para_offset),tem_func->GetParError(1+peak_para_offset));
			tof_x_cento[PeakIndex] = tem_func->GetParameter(1+peak_para_offset);
			tof_x_cento_err[PeakIndex] = tem_func->GetParError(1+peak_para_offset);
	}


	double reduce_chisquare = tem_func->GetChisquare() / tem_func->GetNDF();
	cout<<"X ion: chisquare / NDF = "<<tem_func->GetChisquare() <<" / "<< tem_func->GetNDF()<<" = "<< reduce_chisquare <<endl;
	cout<<"Chisquare accept region = ["<<TMath::ChisquareQuantile(0.05,tem_func->GetNDF()) <<","
		 << TMath::ChisquareQuantile(0.95,tem_func->GetNDF())<<"]"<<endl;

	if(fixpars){
		printf("using: fitlooping(\"%s\",%d,%.4f,%.4f,true)\n",Fitoption,loopNum,fitrangeL,fitrangeR);
	}
	else{
		printf("using: fitlooping(\"%s\",%d,%.4f,%.4f,false)\n",Fitoption,loopNum,fitrangeL,fitrangeR);
	}

	c1->cd(2)->SetEditable(kFALSE);

}


void fithistory(DefGeneralFunc* inhistory, bool fixpars=true){ // it can load the info of FHistory[] to global
	char* ionname = inhistory->name;
	SetOneorTwo	= inhistory->SetOneorTwo;   // 0 = pure gaus;  1 = func1 (gauss + exp Right tail) ; 2 = func2 (gauss + exp Ritht and Left tail)
	NumOfPeaks	= inhistory->NumOfPeaks;    // at least 1 peak
	BackGroundCurve	= inhistory->BackGroundCurve; 
	MainPeakIndex = inhistory->MainPeakIndex;
	ParasLock = inhistory->ParasLock;
	active_tagX	= inhistory->TAG;
	fitrangeL	= inhistory->fitRangeL;
	fitrangeR	= inhistory->fitRangeR;
	active_histoX_Nbins	= inhistory->histo_nbins;
	active_histoX_RangeL	= inhistory->histo_L;
	active_histoX_RangeR	= inhistory->histo_R;

	cout<<"Get ion: "<<ionname<<endl;

	tem_func = inhistory->GetTFunc();

	histo_zoom_in_x(active_tagX,active_histoX_Nbins,active_histoX_RangeL,active_histoX_RangeR);

	
	int peak_para_offset=0;

	if(fixpars){
		if(ParasLock){// only mainpeak has parameters for sigma and tails
			peak_para_offset = Paras_offset_cal(MainPeakIndex-1,SetOneorTwo,MainPeakIndex,ParasLock);
			for(int i=0;i<=SetOneorTwo;i++){
				tem_func->FixParameter(peak_para_offset+2+i,tem_func->GetParameter(peak_para_offset+2+i));
			}
		}
		else{
			for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
				peak_para_offset = Paras_offset_cal(PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock);
				for(int i=0;i<=SetOneorTwo;i++){
					tem_func->FixParameter(peak_para_offset+2+i,tem_func->GetParameter(peak_para_offset+2+i));
				}
			}
		}
	}


	fitlooping("ME",1,fitrangeL,fitrangeR,fixpars);

}


void saveTOFMarker(bool AddMode){
    char fname[200];
	sprintf(fname,"%sTOFMarker_%s.txt",FilePath.c_str(),FileName.c_str());
	bool nofile=true;
	FILE* fp=NULL;
	if((fp = fopen(fname,"r")) == NULL){// no file exist
		cout<<"No TOFMarker record exist; create a new file!"<<endl;
	}
	else{// file exist
		nofile = false;
		fclose(fp);
	}

	if(AddMode == false || nofile){
		fp = fopen(fname,"w");
		fprintf(fp,"name \t mass \t charge \t tag \t laps\n");

	}
	else{
		fp = fopen(fname,"a+");
	}

	// write to txt
	int MaxIndex=40;

	for(int index=0;index<MaxIndex;index++){
		if(marker_tof[index].GetX()>1){  // not empty => to save
			string _name = marker_tof[index].name;
			_name.erase(remove(_name.begin(),_name.end(),' '),_name.end());
			double mass = marker_tof[index].mass;
			int charge = marker_tof[index].charge;
			int tag = marker_tof[index].tag;
			int laps = marker_tof[index].laps;
			fprintf(fp,"%s \t %.3lf \t %d \t %d \t %d\n",_name.c_str(),mass,charge,tag,laps);
		}
	}

	fclose(fp);
}

void readfile_TOFMarker(int Index_start=1, int Index_End=-1){ // read from index =Index_start to Index_End; index >=1; Index_End==-1  => until file end
		FILE * fp=NULL;
    	char fname[200];
		sprintf(fname,"%sTOFMarker_%s.txt",FilePath.c_str(),FileName.c_str());
		fp = fopen(fname,"r");
		if(fp == NULL){cout<<"No file of TOFMarker record exist!!!! Abort"<<endl; return;}

		char filehead[1000];
		bool UseLaps=false;
		fgets(filehead,1000,(FILE*)fp);
		string filehead_s = filehead;
		if(filehead_s.find("laps")!=string::npos){UseLaps=true;}    // has laps information in .txt


		char getname[100];
		double getmass;
		int getcharge;
		int gettag;
		int getlaps;
		string getname_s;

		int index=0;
		for(index=0;index<Index_start-1;index++){fgets(filehead,1000,(FILE*)fp);}

		if(UseLaps){
			for(index=0;index<(Index_End-Index_start+1) || Index_End==-1 ;index++){
				fscanf(fp,"%s \t %lf \t %d \t %d \t %d",getname,&getmass,&getcharge,&gettag,&getlaps);
				getname_s = getname;
				if(getname_s.find("@")!=string::npos){getname_s.erase(getname_s.begin()+getname_s.find("@"),getname_s.end());}

				if(!feof(fp)){
					if(index==0){

						MarkTof(getmass,getcharge,getname_s.c_str(),getlaps,gettag,false,true);
					}
					else{
				 			MarkTof(getmass,getcharge,getname_s.c_str(),getlaps,gettag);
				 	}
				}
				else{ //end of file
					if(index>=40) cout<<"Marker record may be larger than the maximum number can show[40], be careful. Only last 40 markers are valid"<<endl;
					 return;
				}

			}
		}
		else{
			for(index=0;index<(Index_End-Index_start+1) || Index_End==-1 ;index++){
				fscanf(fp,"%s \t %lf \t %d \t %d",getname,&getmass,&getcharge,&gettag);
				if(!feof(fp)){
					if(index==0){
						MarkTof(getmass,getcharge,getname,0,gettag,false,true);
					}
					else{
				 			MarkTof(getmass,getcharge,getname,0,gettag);
				 	}
				}
				else{// end of file
					if(index>=40) cout<<"Marker record may be larger than the maximum number can show[40], be careful. Only last 40 markers are valid"<<endl;
					return;
				}

			}
		}

		if(index>=40) cout<<"Marker record may be larger than the maximum number can show[20], be careful.  Only last 40 markers are valid"<<endl;

}


void ShowResHisto(TH1D* inhisto=NULL, TF1 * inFunc =NULL,double checkrangeL=-1,double checkrangeR=-1){// histogram - funcValue
	if(inhisto != NULL && inFunc !=NULL){
		if(histores != NULL){delete histores; histores = NULL;}
			histores = new TH1D("h_res","Residue Histogram",active_histoX_Nbins,active_histoX_RangeL,active_histoX_RangeR);
			if(checkrangeL ==-1 || checkrangeR == -1){
				cout<<"\e[1;33m"<<"using default range as current fit range!!!"<<"\e[0m"<<endl;
				checkrangeL = fitrangeL;
				checkrangeR = fitrangeR;
			}
			for(int ibin=1;ibin<=active_histoX_Nbins;ibin++){
				double xpos = inhisto->GetBinCenter(ibin);
				double funcvalue =0;
				if(xpos>=checkrangeL && xpos<=checkrangeR){funcvalue = inFunc->Eval(xpos);}
				histores->SetBinContent(ibin,inhisto->GetBinContent(ibin) - funcvalue);
			}

			histores->Draw("same");

	}
	else{
		cout<<"invalid input histogram!! abort!!"<<endl;
		return;
	}

}


void ReverseEjeTime(){
	swap(TagBit0,TagBit1);
	
}

void ReleasePar(){

	for(int i=0;i<man1_func1->GetNpar();i++){man1_func1->ReleaseParameter(i);}
	for(int i=0;i<ml_func1->GetNpar();i++){ml_func1->ReleaseParameter(i);}
	for(int i=0;i<mexp_func1->GetNpar();i++){mexp_func1->ReleaseParameter(i);}
	for(int i=0;i<multi2_func1->GetNpar();i++){multi2_func1->ReleaseParameter(i);}
	for(int i=0;i<multi2_line_func1->GetNpar();i++){multi2_line_func1->ReleaseParameter(i);}
	for(int i=0;i<multi2_exp_func1->GetNpar();i++){multi2_exp_func1->ReleaseParameter(i);}


	for(int i=0;i<man1_func2->GetNpar();i++){man1_func2->ReleaseParameter(i);}
	for(int i=0;i<ml_func2->GetNpar();i++){ml_func2->ReleaseParameter(i);}
	for(int i=0;i<mexp_func2->GetNpar();i++){mexp_func2->ReleaseParameter(i);}
	for(int i=0;i<multi2_func2->GetNpar();i++){multi2_func2->ReleaseParameter(i);}
	for(int i=0;i<multi2_line_func2->GetNpar();i++){multi2_line_func2->ReleaseParameter(i);}
	for(int i=0;i<multi2_exp_func2->GetNpar();i++){multi2_exp_func2->ReleaseParameter(i);}

}


///////&&&&&&&&&&&&&&&&&&&&&&& beta-tof coincidence &&&&&&&&&&&&&&&&&&&&&&&&&&&&

#include"Beta_Beta_coin.cc"
#include"BetaEval.h"

TFile * f_beta_tof_tree=NULL;
TTree* tbeta_tof = NULL;// = new TTree("tbeta_tof","tbeta_tof");  // tree for store coincidence event
TH1D* h_beta_decay=NULL;
double* xedgeL=NULL; // variable bin size for decay time; log scale

TF1* dk_curve_ref = NULL;
TF1* dk_curve_m = NULL;
TLegend* dk_legend=NULL;
double Times2Halflive = 4;  // how many times of RI halflife are considered to be coincident with beta
double Beta_bg_rate=0;
double detector_effi = 0.33; // detector efficiency 
int Ntof_counter=0;  // number of TOF cut by cursor for beta-tof coincidence

double Multiple_max =4; // for fitting of efficiency (control time window width in the unit of halflife)
double Multiple_diff=0.2;

//***************** tof and sweep w/ or w/o coin  *************
vector<double> eve_sele_tof;
vector<Long64_t> eve_sele_glo_sweep;

vector<double> eve_woutcoin_tof;
vector<Long64_t> eve_woutcoin_sweep;
TH1D* h_no_coin=NULL;
TH1D* h_sele = NULL;
//*********************************************************

vector<string> MRTOF_hit_Info; // for saving the information (global sweep, global time in [s]) of selected TOF hits for beta-TOF measurement

class RejectionList{
	private:
		vector<Long64_t> tof_gclock;
		vector<Long64_t> time1_gclock;
		vector<Long64_t> time2_gclock;
		int NumofReject;

	public:
		RejectionList():NumofReject(0){;}
		void Clear(){
			tof_gclock.clear();
			time1_gclock.clear();
			time2_gclock.clear();
			NumofReject=0;
		}

		void Add(Long64_t _tof_gclock , Long64_t _time1_gclock, Long64_t _time2_gclock){ // add one reject to list
			tof_gclock.push_back(_tof_gclock);
			time1_gclock.push_back(_time1_gclock);
			time2_gclock.push_back(_time2_gclock);
			NumofReject++;
		}

		void Remove(Long64_t _tof_gclock , Long64_t _time1_gclock, Long64_t _time2_gclock){

			vector<Long64_t>tof_gclock_tem;
			vector<Long64_t>time1_gclock_tem;
			vector<Long64_t>time2_gclock_tem;
			for(unsigned long i=0;i<tof_gclock.size();i++){
				if(_tof_gclock == tof_gclock[i] && _time1_gclock== time1_gclock[i] && _time2_gclock== time2_gclock[i]){ continue;}
				else{
					tof_gclock_tem.push_back(tof_gclock[i]);
					time1_gclock_tem.push_back(time1_gclock[i]);
					time2_gclock_tem.push_back(time2_gclock[i]);
				}
			}

			Clear();

			for(unsigned long i=0;i<tof_gclock_tem.size();i++){
				Add(tof_gclock_tem[i],time1_gclock_tem[i],time2_gclock_tem[i]);
			}

		}

		void Print(){
			printf("tof_gclock \t time1_glock \t time2_gclock\n");
			for(unsigned long i=0;i<tof_gclock.size();i++){
				printf("*%lld \t *%lld \t *%lld*\n",tof_gclock[i],time1_gclock[i],time2_gclock[i]);
			}
		}

		bool IsInside(Long64_t _tof_gclock , Long64_t _time1_gclock, Long64_t _time2_gclock){
			for(unsigned long i=0;i<tof_gclock.size();i++){
				if(_tof_gclock == tof_gclock[i] && _time1_gclock== time1_gclock[i] && _time2_gclock== time2_gclock[i]){ return true;}
			}
			return false;
		}

		string RejectList_s(){
			if(NumofReject==0) return "";
			string output=" && ";
			for(unsigned long i=0;i<tof_gclock.size();i++){
				output += "!((time_gclock_mr == " + to_string(tof_gclock[i]) + ") && ( time1 ==" + to_string(time1_gclock[i]) + ") && ( time2 ==" + to_string(time2_gclock[i]) + ")) && ";
			}
			output.erase(output.begin()+output.find_last_of("&&")-2,output.end());
			return output;
		}

		int GetN(){ return NumofReject;}

};



RejectionList rejectlist;
bool RecordOnly_1st_beta=true;
bool RecordOnly_1st_beta_neg = RecordOnly_1st_beta;


void GetEveNoCoin(vector<double>& sele_tof,vector<Long64_t>& sele_sweep,vector<double>& coin_tof, vector<Long64_t>& coin_sweep){//(imput sele_tof, input sele sweep, input coin_tof, input coin sweep)
	eve_woutcoin_tof.clear();
	eve_woutcoin_sweep.clear();

	if(sele_tof.size() ==0 ){cout<<"\e[1;33m"<<"No selected events"<<"\e[0m"<<endl; return;}
	else if(coin_tof.size() ==0){
		cout<<"\e[1;33m"<<"No coincident events"<<"\e[0m"<<endl;
		for(unsigned long i = 0; i<sele_tof.size(); i++){
			eve_woutcoin_tof.push_back(sele_tof[i]);
			eve_woutcoin_sweep.push_back(sele_sweep[i]);
		}
	}
	else{
			for(unsigned long i = 0; i<sele_tof.size(); i++){
				bool In_coin_list = false;
					for(unsigned long j=0;j<coin_tof.size();j++){
						if(sele_tof[i] == coin_tof[j] && sele_sweep[i] == coin_sweep[j]){
							In_coin_list = true;
							break;
						}
					}
				
				if(!In_coin_list){
					eve_woutcoin_tof.push_back(sele_tof[i]);
					eve_woutcoin_sweep.push_back(sele_sweep[i]);
				}
			}
	}
	
	eve_woutcoin_tof.shrink_to_fit();
	eve_woutcoin_sweep.shrink_to_fit();

}


void Showhisto_except_coin(int which_all1_nocoin2_both3){
	// show histo 0: remove all extra histogram in c1->cd(2)
	// 1: add histogram showing all selected event before coincident measurement
	// 2: add histogram showing events without coincidence
	// 3: add both histogram (1 and 2)
	if(h_zoom_x == NULL){cout<<"\e[1;33m"<<"Break!!! h_zoom_x no exist!!"<<"\e[0m"<<endl; return;}
	if(eve_woutcoin_tof.size() ==0){cout<<"All tof have coincidence"<<endl; return;}
	if(h_zoom_x == NULL){cout<<"h_zoom_x is not exist, create it first!!"<<endl; return;}

	int Nbins = h_zoom_x->GetNbinsX();
	double EdgeL = h_zoom_x->GetBinLowEdge(1);
	double EdgeR = h_zoom_x->GetBinLowEdge(Nbins)+h_zoom_x->GetBinWidth(1);

	if(h_no_coin != NULL){delete h_no_coin; h_no_coin = NULL;}

	if(h_sele != NULL){delete h_sele; h_sele = NULL;}

	if(which_all1_nocoin2_both3 == 0){
		cout<<"Delete h_sele and h_no_coin"<<endl; 
		double Ymax_h_zoom_x = h_zoom_x->GetBinContent(h_zoom_x->GetMaximumBin());
		if(c1->cd(2)->GetLogy()) h_zoom_x->SetMaximum(Ymax_h_zoom_x*3);
        else h_zoom_x->SetMaximum(Ymax_h_zoom_x*1.08);
        c1->cd(2)->SetEditable(kTRUE);
        c1->cd(2)->Modified(); c1->cd(2)->Update();
		return;
	}

	h_no_coin = new TH1D("h_no_coin","h_no_coin",Nbins,EdgeL,EdgeR);
	h_sele = new TH1D("h_sele","h_sele",Nbins,EdgeL,EdgeR);

	long nevents = eve_woutcoin_tof.size();
	if(nevents>0){
		for(long index=0;index< nevents;index++){
			h_no_coin->Fill(eve_woutcoin_tof[index]);
		}
		cout<<"Event w/o coincidence = "<<nevents<<endl;
	}else{
		cout<<"No events for h_no_coin  !!"<<endl;
	}

	nevents = eve_sele_tof.size();

	if(nevents>0){
		for(long index=0;index< nevents;index++){
			h_sele->Fill(eve_sele_tof[index]);
		}
		cout<<"Total selected events = "<<nevents<<endl;
	}
	else{
		cout<<"No events for h_sele !!"<<endl;
	}

	h_no_coin->SetLineColor(kSpring-6);
	h_no_coin->SetLineWidth(2);
	h_sele->SetLineColor(kGray+2);
	h_sele->SetLineWidth(4);
	c1->cd(2)->SetEditable(kTRUE);

	if(which_all1_nocoin2_both3 == 1){
		if(h_sele->GetEntries()>0) h_sele->Draw("same");
	}
	else if(which_all1_nocoin2_both3 == 2){
		if(h_no_coin->GetEntries()>0) h_no_coin->Draw("same");
	}
	else{
		if(h_no_coin->GetEntries()>0) h_no_coin->Draw("same");
		if(h_sele->GetEntries()>0) h_sele->Draw("same");
	}

	double h_max_amp = -1;

 	TH1D* h_get_from_canvas = (TH1D*) c1->cd(2)->GetPrimitive("h_sele");
   
  	if(h_get_from_canvas != NULL){
        int ibin = h_get_from_canvas->GetMaximumBin();
        h_max_amp = h_get_from_canvas->GetBinContent(ibin);
  	}

  	h_get_from_canvas = (TH1D*) c1->cd(2)->GetPrimitive("h_no_coin");
  	if(h_get_from_canvas != NULL){
	      int ibin = h_get_from_canvas->GetMaximumBin();
	    if(h_max_amp<h_get_from_canvas->GetBinContent(ibin) && h_get_from_canvas->GetBinContent(ibin)>0){
	       h_max_amp = h_get_from_canvas->GetBinContent(ibin);
	    }
    }


    if(h_max_amp>h_zoom_x->GetBinContent(h_zoom_x->GetMaximumBin())){
        if(c1->cd(2)->GetLogy()) h_zoom_x->SetMaximum(h_max_amp*3);
        else h_zoom_x->SetMaximum(h_max_amp*1.08);
    } 


	c1->cd(2)->Modified(); c1->cd(2)->Update();

	c1->cd(2)->SetEditable(kFALSE);


}


//bool BetaFastMode=false;

void FindBetaTof_Coin(double tof_L, double tof_R, Long64_t _Halflife, TCutG** usecut=NULL){
		fin->cd();

		T_data_neg.clear(); // defined in chance.cc
		T_data_neg.shrink_to_fit();

		MRTOF_hit_Info.clear(); MRTOF_hit_Info.shrink_to_fit();

		//***** save global sweeps and tof for all events in the selected range************
		eve_sele_tof.clear();  eve_sele_tof.shrink_to_fit(); // clear to store selected event in cd(2)
		eve_sele_glo_sweep.clear();  eve_sele_glo_sweep.shrink_to_fit();


		if(h_beta_accumulate !=NULL){ delete h_beta_accumulate; h_beta_accumulate =NULL;}

		static double tof_L_old=-1;   // for remember setting of last run; skip when setting has no change
		static double tof_R_old=-1;
		static Long64_t _Halflife_old =-1;

		if(tof_L_old==tof_L && tof_R_old==tof_R && _Halflife_old==_Halflife){
			cout<<"\e[1;33m"<<"Beta-Tof coincident condition no change! break"<<"\e[0m"<<endl; 
			return;
		}

		tof_L_old = tof_L;
		tof_R_old = tof_R;
		_Halflife_old = _Halflife;


		vector <Double_t> *mr_time = new vector <Double_t>();
		vector <Long_t> * mr_sweeps_global = new vector <Long_t>;
		vector <Int_t> *mr_tag = new vector <Int_t>();
		vector <Int_t> *mr_channel = new vector <Int_t>();


		intree->ResetBranchAddresses();
		intree->SetBranchAddress("sweeps_global",&mr_sweeps_global);
		intree->SetBranchAddress("tag",&mr_tag);
		intree->SetBranchAddress("channel",&mr_channel);
		string treename = intree->GetName();
		if(treename == "tree") intree->SetBranchAddress("time",&mr_time);
		else{ intree->SetBranchAddress("timec",&mr_time);}



		//**************** Set up tree to store coincident event ***************
	/*	if(h_beta_decay!=NULL)delete h_beta_decay;
		if(h_zoom_x!=NULL){
			int nbins = h_zoom_x->GetNbinsX();
			double tem_left = h_zoom_x->GetBinCenter(1);
			double tem_right = h_zoom_x->GetBinCenter(nbins);
			h_beta_decay = new TH1D("h_beta_decay","beta-tof coincident",nbins,tem_left,tem_right); cout<<"create Bbabab"<<endl; 

			h_beta_decay->GetXaxis()->SetTitle(h_zoom_ref->GetXaxis()->GetTitle());
			h_beta_decay->GetXaxis()->SetTitleSize(h_zoom_ref->GetXaxis()->GetTitleSize());
			h_beta_decay->GetXaxis()->CenterTitle();
			h_beta_decay->GetYaxis()->SetTitle(h_zoom_ref->GetYaxis()->GetTitle());
			h_beta_decay->GetYaxis()->SetTitleSize(h_zoom_ref->GetYaxis()->GetTitleSize());
			h_beta_decay->GetYaxis()->CenterTitle();
			cout<<"oooooooo"<<endl;
		}
		else{
				h_beta_decay = new TH1D();
				h_beta_decay->SetName("h_beta_decay");
				h_beta_decay->SetTitle("beta-tof coincident");
		}
*/

		//h_beta_decay->Reset();

		if(h_beta_decay!=NULL) delete h_beta_decay; // must be deleteed firstly
		if(tbeta_tof!=NULL){ 	delete tbeta_tof;}  // must be deleted first

		if(f_beta_tof_tree!=NULL){
			if(f_beta_tof_tree->IsOpen()){f_beta_tof_tree->Close(); delete f_beta_tof_tree;}
			else{delete f_beta_tof_tree;}
		}

		string f_beta_tof_tem = FilePath + "tbeta_tof_tree_tem.root";
		f_beta_tof_tree = new TFile(f_beta_tof_tem.c_str(),"RECREATE");
		if(f_beta_tof_tree->IsOpen()){
			f_beta_tof_tree->cd();
		}
		else{
			cout<<"\e[1;33m"<<"fail to create tbeta_tof_tree_tem.root!!! Abort!!!"<<"\e[0m"<<endl;
			return;
		}

		
		tbeta_tof = new TTree("tbeta_tof","tbeta_tof");


		Long64_t mr_Bingo_sweeps_global=0; double mr_Bingo_time=0; Long64_t mr_Bingo_time_gclock=0;
		Long64_t time1=0; Long64_t time1_relative=0; Long64_t sweeps_gclock_1=0; int adc1=0; double EkeV_1=0;
		Long64_t time2=0; Long64_t time2_relative=0;  Long64_t sweeps_gclock_2=0; int adc2=0; double EkeV_2=0;
		int delta_t=0;
		Long64_t delta_t_decay=0;
		int tof_index=0;

		tbeta_tof->Branch("Ntof_counts",&Ntof_counter,"Ntof_counts/I");
		tbeta_tof->Branch("tof_index",&tof_index,"tof_index/I");
		tbeta_tof->Branch("sweeps_global_mr",&mr_Bingo_sweeps_global,"sweeps_global_mr/L");
		tbeta_tof->Branch("time_mr",&mr_Bingo_time,"time_mr/D");
		tbeta_tof->Branch("time_gclock_mr",&mr_Bingo_time_gclock,"time_gclock_mr/L");
		tbeta_tof->Branch("sweeps_gclock_1",&sweeps_gclock_1,"sweeps_gclock_1/L");
		tbeta_tof->Branch("time1",&time1,"time1/L");
		tbeta_tof->Branch("time1_relative",&time1_relative,"time1_relative/L");
		tbeta_tof->Branch("adc1",&adc1,"adc1/I");
		tbeta_tof->Branch("EkeV_1",&EkeV_1,"EkeV_1/D");
		tbeta_tof->Branch("sweeps_gclock_2",&sweeps_gclock_2,"sweeps_gclock_2/L");
		tbeta_tof->Branch("time2",&time2,"time2/L");
		tbeta_tof->Branch("time2_relative",&time2_relative,"time2_relative/L");
		tbeta_tof->Branch("adc2",&adc2,"adc2/I");
		tbeta_tof->Branch("EkeV_2",&EkeV_2,"EkeV_2/D");
		tbeta_tof->Branch("delta_t",&delta_t,"delta_t/I");
		tbeta_tof->Branch("delta_t_decay",&delta_t_decay,"delta_t_decay/L");


		//*******************************************************************

		Long64_t N_entries = intree->GetEntriesFast();
		int bingocounter=0;

		//************************
		Long64_t mr_candidate_sweeps_global_old=-1; // for looping , quickly back to same sweep at multihit case;
		Long64_t line_index_mark = -1 ;  // for looping , quickly back to same sweep at multihit case;
			//%%% beta background rate %%%
		bool Cutvalid=false;
		double beta_bg_rate=0; // basing on the first closest beta to the tof
		int counts_beta_bg=0;  // counting the number of first closeest beta to the tof
		vector<Long64_t> mr_sequence; //global tof of all ions in specific the peak
		double Tof_rate=0;  // calculate the counting rate of specific ion by the average time interval
		//***************************

		if(CutModifyCounter>0) Cutvalid=true;


int ion_counts=0;
bool NewFileChecking=false;
bool GetNtofDone = false;

Long64_t sweeps_globalTime_From_betaLST=0;  // for extracting the global time of give sweeps number from "tbeta_sweeps_gtime" tree;

unsigned long First_neg_Beta=0;    // to speed up searching in Match bate list
bool Allow_new_First_neg_Beta=true;

			//%%%%%%%%%%%% get raw beta from tree  %%%%%%%%%%%%%%%%%%%%
				Long64_t sweeps_gclock_AT_line_index =0;
				Long64_t time_AT_line_index=0;
				int channel_AT_line_index=0;
			if(!BetaFastMode && MakeBetaRawTree){
				if(tbeta_Raw!=NULL){
						nevt_bta = tbeta_Raw->GetEntries();
						tbeta_Raw->SetBranchAddress("time",&time_AT_line_index);
						tbeta_Raw->SetBranchAddress("channel",&channel_AT_line_index);
						tbeta_Raw->SetBranchAddress("sweeps_global",&sweeps_gclock_AT_line_index);
				}else{
					cout<<"\e[1;31m Error!! beta raw tree not available, fail to match tof and beta!!! Abort.\e[0m"<<endl;
					return;
				}
			}
			//5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if(!BetaFastMode || MakeBetaRawTree){
			if(!BetaSweepsTimeReady || tbeta_sweeps_gtime ==NULL){
				cout<<endl;
				cout<<"\e[1;31m Error!!! Faile to get sweeps global time information from \"tbeta_sweeps_gtime\" tree!! Abort!!!\e[0m"<<endl;
				return;
			}else{
				tbeta_sweeps_gtime->SetBranchAddress("sweeps_gtime",&sweeps_globalTime_From_betaLST);
			}
		}
		

		//work for both decoders with event at hit and event at sweeps **********************
		for(Long64_t event_i=0;event_i<N_entries;event_i++){
				if(event_i%1000==0) cout<<" Matching Tof and Beta......   "<<event_i<<" / "<<N_entries<<"\r"<<flush;
				mr_time->clear();
				mr_sweeps_global->clear();
				mr_tag->clear();
				mr_channel->clear();
				double tof_last_hit=0;
				intree->GetEntry(event_i);
	
				
				for(unsigned long hit_i=0;hit_i<mr_time->size();hit_i++){
						if(mr_time->at(hit_i) <10 ) continue;  // trigger marker
						if(mr_tag->at(hit_i)==1) break; // reference side no RI
						if(mr_channel->at(hit_i) != Tof_Channel2Read) continue;

						if(mr_time->at(hit_i)-tof_last_hit < Tof_veto) continue;
						else tof_last_hit = mr_time->at(hit_i);
							 
						if(mr_time->at(hit_i)>= tof_L && mr_time->at(hit_i)<=tof_R){  // tof peak candidate;

if(!GetNtofDone){
	ion_counts++;
	continue;
}
Ntof_counter = ion_counts;
tof_index++;
								eve_sele_tof.push_back(mr_time->at(hit_i)); // store selected tof
								Long64_t mr_candidate_sweeps_global = mr_sweeps_global->at(hit_i);
								eve_sele_glo_sweep.push_back(mr_candidate_sweeps_global); // store selected sweep_glo

if(tbeta_sweeps_gtime!=NULL)tbeta_sweeps_gtime->GetEntry(mr_candidate_sweeps_global);


								Long64_t mr_candidate_sweeps_global_gloTime = 0;
								if(mr_candidate_sweeps_global_old!=mr_candidate_sweeps_global){  // new sweeps comes , search position again in beta.lst
									 line_index_mark = mr_candidate_sweeps_global;  // first round to determine corresponding sweeps position at beta .lst
									 mr_candidate_sweeps_global_old=mr_candidate_sweeps_global;
								}

								for(Long64_t line_index = line_index_mark;line_index<=nevt_bta || BetaFastMode;line_index++){ // BetaFastMode load beta file based on tof  // general mode load all beta file in advance
									if(BetaFastMode || BetaSweepsTimeReady ) goto SkipSweepGlobalTimeGetFromBetaListFile;

										if(MakeBetaRawTree && tbeta_Raw!=NULL){
												tbeta_Raw->GetEntry(line_index-1);
										}else{
											sweeps_gclock_AT_line_index = sweeps_gclock_v->at(line_index-1);
											time_AT_line_index = time_v->at(line_index-1);
											channel_AT_line_index = channel_v->at(line_index-1);
										}


									if(sweeps_gclock_AT_line_index > mr_candidate_sweeps_global){ // line_index-1 because sweeps start from 1
										cout<<"line in beta list = "<<line_index<<endl;
										cout<<"find out sweeps in beta list = "<<sweeps_gclock_AT_line_index<<endl;
										cout<<"mr sweeps_global candidate = "<<mr_candidate_sweeps_global<<endl;
											cout<<"sweeps error in beta or mrtof .lst file. Abort!!!!"<<endl;
											return;
									}
									else if(sweeps_gclock_AT_line_index < mr_candidate_sweeps_global){										
										continue;
									}
									else{  // sweeps_gclock_v = mr_candidate_sweeps_global
										// get the global time of that sweep
										if(channel_AT_line_index ==syn_CH){
											line_index_mark = line_index;

SkipSweepGlobalTimeGetFromBetaListFile:
if(!BetaFastMode){//precise mode, get the global time of sweep from beta.lst data

									if(BetaSweepsTimeReady){
										mr_candidate_sweeps_global_gloTime = sweeps_globalTime_From_betaLST;
									}
									else{
										mr_candidate_sweeps_global_gloTime = time_AT_line_index;   // get mrtof sweep at global clock
									}
											
											mr_candidate_sweeps_global_gloTime = mr_candidate_sweeps_global_gloTime - ModuleDelay; // nim  form synchronizer to TTL show some delay
				//	cout<<mr_candidate_sweeps_global_gloTime<<endl;


}
else{// base on 25 ms per sweep
		Long64_t timepersweep = 25000000;

	//	if(BetaSweepsTimeReady){
	//		mr_candidate_sweeps_global_gloTime = sweeps_globalTime_From_betaLST;
	//	}
	//	else{
			mr_candidate_sweeps_global_gloTime = (mr_candidate_sweeps_global-1)*timepersweep;
	//	}
		
//		cout<<"global sweep time = "<<mr_candidate_sweeps_global_gloTime/1e9<<"global sweep = "<<mr_candidate_sweeps_global<<endl;
		//break;



}

											Long64_t mr_candidate_gloTime = (Long64_t)mr_time->at(hit_i) + mr_candidate_sweeps_global_gloTime;
											//cout<<endl;
											//cout<<endl;
					MRTOF_hit_Info.push_back(Form("Found MRTOF hit #%d at global sweep =%lld \t  at global clock = %.3f [s]\n",ion_counts,mr_candidate_sweeps_global,mr_candidate_gloTime/1e9));
//break;
											Long64_t time1_temp = 0;
											Long64_t time_diff_mr_beta=0;

											Long64_t time_diff_mr_background_beta=0;
											bool Isnoise=false,usetimecut=false;

											mr_sequence.push_back(mr_candidate_gloTime);

if(BetaFastMode){
	bool IsNewBetaFile=false;
	        if(Current_First_beta_file_FastMode != ToLoad_beta_file ){  // beta file read to vector The name of Current first beta file in the series update by Read_Beta_lst_batch()
	        																// once Read_Beta_lst_batch() has been executed under BetaFastMode for a beta file series
	        																// Then Current_First_beta_file_FastMode always be == ToLoad_beta_file, in this case, the list of all file name and start time of each file at global clock have been created already
	        																// set IsNewBetaFile = false to skip the reading for creating the list of all file name and start time at global clock
	            if(ToLoad_beta_file=="+++"){
	              cout<<"warning: no beta file to be load; please load file first"<<endl;
	              cout<<"\e[1;33m"<<"beta filename with .lst"<<"\e[0m"<<endl;
	              cout<<"using:  LoadNewBetaFile("<<endl;
	              return;
	            }
	            else{
	//              Read_Beta_lst(FilePath+"../","LST/"+ToLoad_beta_file);
	 //               Read_Beta_lst_batch(FilePath+"../","LST/"+ToLoad_beta_file);
	            	IsNewBetaFile = true;
	            }

            }

            if(IsNewBetaFile && !NewFileChecking){
            	 cout<<"\e[1;33m"<<"Loading key Info of file series start at: "<<ToLoad_beta_file<<endl;
            }
            else if(!IsNewBetaFile && !NewFileChecking){
            	cout<<"\e[1;33m"<<"Key Info of file series start at: "<<ToLoad_beta_file<<" had been loaded. Skip loading....."<<endl;
            }

            NewFileChecking = true;  // finish the checking to see whether going to read a new series of file, avoid to print repetitive information

    bool ReadBetaStatus =false;  // to know the reading of beta list is successful or not

     //%%%%%%%%%%% offline  %%%%%%%%%%%%%5
   if(!ReadBetaOnline) ReadBetaStatus = Read_Beta_lst_batch(FilePath+"../","LST/"+ToLoad_beta_file,PreSetADC1_low_Thres_ch ,PreSetADC2_low_Thres_ch ,mr_candidate_gloTime-Times2Halflive*_Halflife, mr_candidate_gloTime+Times2Halflive*_Halflife, IsNewBetaFile);
	//%%%%%%%%%%% online %%%%%%%%%%%%%5
	else ReadBetaStatus = Read_Beta_lst_batch(OnlineBetaFilePath,ToLoad_beta_file,PreSetADC1_low_Thres_ch ,PreSetADC2_low_Thres_ch ,mr_candidate_gloTime-Times2Halflive*_Halflife, mr_candidate_gloTime+Times2Halflive*_Halflife, IsNewBetaFile);
	
	if(!ReadBetaStatus){cout<<"\e[1;31m Fail to read beta list file. Abort!!!! Please check. \e[0m"<<endl; return;	}

	Compare_Beta_Beta();
}

												bool Tof_with_coincidence = false; // True: either negative or positive beta, this tof will be save in tree; False ==> No any negative or positive beta found

												Allow_new_First_neg_Beta =true;

												for(unsigned long index=First_neg_Beta;index< Match_hit.size();index++){  // search in beta-beta coinci list
														// histo_beta_deltaT->Fill(Match_hit[index].delta_t);	
													time1_temp = Match_hit[index].B1_time;
													time_diff_mr_beta = time1_temp - mr_candidate_gloTime;
													usetimecut=false;	Isnoise=false;

													//********** checking background rate at 8 half live before ***************
													if( time_diff_mr_beta<=0 && true){// background beta candidate
														//time_diff_mr_beta>=-8*_Halflife && time_diff_mr_beta<=-1000 (original condition)

														/*
														if(usecut!=NULL){ // use cut

															if(usecut[0]!=NULL){ // cut available
																Cutvalid=true;
																for(int ic=1;ic<5;ic++){
																	if(usecut[ic]!=NULL){
																		usetimecut = true;
																		Isnoise = Isnoise || usecut[ic]->IsInside((double)Match_hit[index].B1_time_relative,(double)Match_hit[index].B2_time_relative);  // in any of the time noise cut
																	} 
																}
															

																if(usetimecut)	Isnoise = Isnoise && usecut[0]->IsInside(ADC2keV_1(Match_hit[index].B1_adc),ADC2keV_2(Match_hit[index].B2_adc)); // use timecut:  noise should be in time cut and energy cut at the same time
																else Isnoise = usecut[0]->IsInside(ADC2keV_1(Match_hit[index].B1_adc),ADC2keV_2(Match_hit[index].B2_adc));   // only energy cut;

																if(!Isnoise){
																	if(time_diff_mr_background_beta == 0 || time_diff_mr_beta>time_diff_mr_background_beta){// close to 0(Tof moment) get the beta closet to Tof
																		time_diff_mr_background_beta = time_diff_mr_beta;
																	}

																	if(time_diff_mr_beta > -Times2Halflive * _Halflife && !RecordOnly_1st_beta_neg){
																		 T_data_neg.push_back( (double)time_diff_mr_beta*1e-9 );
																	}

																}
															}

														}*/

																	
									if(time_diff_mr_background_beta == 0 || time_diff_mr_beta>time_diff_mr_background_beta){// close to 0(Tof moment) get the beta closet to Tof
										time_diff_mr_background_beta = time_diff_mr_beta;
									}


									if(time_diff_mr_beta > -Times2Halflive * _Halflife && !RecordOnly_1st_beta_neg){ // record all beta in [-halflife*times2halflife,0]
										 T_data_neg.push_back( (double)time_diff_mr_beta*1e-9 );
									}




													}
													//************* matchiing Tof-Beta pair ******************
													if(time_diff_mr_beta>-Times2Halflive*_Halflife && time_diff_mr_beta< Times2Halflive*_Halflife){ // decay signal from 1 us before to 5* half life

														if(Allow_new_First_neg_Beta){First_neg_Beta = index; Allow_new_First_neg_Beta = false;}

														Tof_with_coincidence = true; // negative and positive are both regarded as coincidence
															// bingo => save to tree
														if(time_diff_mr_beta>0) bingocounter++;  // get counts of beta candidate at time position side
														mr_Bingo_sweeps_global = mr_candidate_sweeps_global;
														mr_Bingo_time = mr_time->at(hit_i);
														mr_Bingo_time_gclock = mr_candidate_gloTime;
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
														delta_t_decay=(time1+time2)/2 - mr_Bingo_time_gclock;
														tbeta_tof->Fill();
														//if(BetaFastMode) tbeta_tof->Write(); // IMPORTANT!!! write to file avoid crash when variable delete
														if(time_diff_mr_beta>0 && RecordOnly_1st_beta)break; // --> to get next tof

													}

													else if(time_diff_mr_beta>=Times2Halflive*_Halflife) break;


												}// for loop search in beta-beta list



												if(!Tof_with_coincidence){ // just to save tof without any coincidence in tree as well 
														mr_Bingo_sweeps_global = mr_candidate_sweeps_global;
														mr_Bingo_time = mr_time->at(hit_i);
														mr_Bingo_time_gclock = mr_candidate_gloTime;
														sweeps_gclock_1 = 0;
														time1 = 0;
														time1_relative = 0;
														adc1 = -1000;
														EkeV_1 = ADC2keV_1(adc1);
														sweeps_gclock_2 = 0;
														time2 = 0;
														time2_relative = 0;
														adc2 = -1000;
														EkeV_2 = ADC2keV_2(adc2);
														delta_t = 0;
														delta_t_decay=0;
														tbeta_tof->Fill();
												}


					if(BetaFastMode){f_beta_tof_tree->cd(); tbeta_tof->Write(); }// IMPORTANT!!! write to file avoid crash when variable delete


												if(time_diff_mr_background_beta!=0){// record time of background beta
														beta_bg_rate+=time_diff_mr_background_beta;
														if(RecordOnly_1st_beta_neg)T_data_neg.push_back((double)time_diff_mr_background_beta *1e-9); // record the first beta at negative time side from Tof
														time_diff_mr_background_beta=0;
														counts_beta_bg++;
														
												}

											break;// finish the searching for the global time of the mrtof sweeps in beta.lst raw file

										}// get sweeps of mrtof in global time
										
									}// sweeps_gclock_v = mr_candidate_sweeps_global  (beta raw file : mrtof tree)

								}// for loop searching matching sweep in beta raw file.

//break; // skip multihit in same sweep
						}// the end of if tof candidate is in [tof_L, tof_R] (selected peak region)

				}// loop for mrtof multihit in a single sweep
//cout<<"ion_real_counts: "<<ion_counts<<endl;

				if(!GetNtofDone && event_i==N_entries-1){ // first loop to get the Ntof; second loop to get tof- beta coincident events
					event_i=-1;  // so that event_i++ makes event_i back to 0; restart the loop for getting tree one by one again
					GetNtofDone=true;
				}


		}//loop all event in mrtof
cout<<endl;
cout<<"ion_real_counts: "<<ion_counts<<endl;
cout<<endl;

if(!BetaFastMode && MakeBetaRawTree && tbeta_Raw != NULL) tbeta_Raw->ResetBranchAddresses(); 

if((!BetaFastMode || MakeBetaRawTree) && tbeta_sweeps_gtime != NULL ) tbeta_sweeps_gtime->ResetBranchAddresses();


//return;
if(BetaFastMode)ShowBeta_Beta_Coin();		

if(!BetaFastMode)	tbeta_tof->Write();	

		//delete tbeta_tof;	
		f_beta_tof_tree->Close(); delete f_beta_tof_tree;

		f_beta_tof_tree = new TFile(f_beta_tof_tem.c_str(),"UPDATE");		
		
		tbeta_tof = (TTree*)f_beta_tof_tree->Get("tbeta_tof");

		//tbeta_tof = (TTree*)f_beta_tof_tree->Get("tbeta_tof");

		c1->cd(2)->SetEditable(kTRUE);
		fin->cd();
		
		if(tbeta_tof->GetEntriesFast()>0){
				 tbeta_tof->Draw("time_mr>>h_zoom_x");   // this kind of draw must switch back to the file where the histogram locates
				 active_tree_name = "tbeta_tof";
		}
	//	h_zoom_x->Draw();
		c1->cd(2)->Update();

		cout<<" Create beta-tof tree:  tbeta_tof->Print() ;  tbeta_tof->Scan()"<<endl;
		cout<<"Found Number of positive beta-tof coincident events: "<< bingocounter <<endl;

		//************** to show time spectrum of deta decay **************
		f_beta_tof_tree->cd();
		double h_uplimit = Times2Halflive * RIHalflive * 1e-9;  // from [ns] to [s]
		double h_lowlimit=-2;

		h_uplimit = (int)TMath::Log10(h_uplimit)+2;

		//h_uplimit = TMath::Power(10,h_uplimit);

		if(xedgeL!=NULL){ delete xedgeL; xedgeL=NULL;}

		int numbins = (int)((h_uplimit-h_lowlimit)*9);
		
		xedgeL = new double[numbins+1];

		int xindex=0;
		for(int order = h_lowlimit;order<h_uplimit; order++){
			int digit=1;
			while(digit<10){
				xedgeL[xindex++] = (digit++)*TMath::Power(10,order);
			}
		}

		xedgeL[xindex++] = TMath::Power(10,h_uplimit); // right edge of last bin


		//h_beta_decay = new TH1D("h_beta_decay","beta decay time distribution;decay time [s];counts",(int)((h_uplimit-0.01)/0.02),0.01,h_uplimit);
				//h_beta_decay = new TH1D("h_beta_decay","beta decay time distribution;decay time [s];counts",numbins,xedgeL); // various bin width for log scale
		h_beta_decay = new TH1D("h_beta_decay","beta decay time distribution;Decay time [s];Counts",50,-3,2);
		h_beta_decay->GetXaxis()->SetTitleSize(0.05);
		h_beta_decay->GetXaxis()->CenterTitle();
		h_beta_decay->GetYaxis()->SetTitleSize(0.05);
		h_beta_decay->GetYaxis()->CenterTitle();

		c_beta_beta->cd(4);
		   //c_beta_beta->cd(4)->SetLogx();
	//	tbeta_tof->Draw("delta_t_decay/1e9>>h_beta_decay")	;

		tbeta_tof->ResetBranchAddresses();
		Long64_t decay_time=0;
		tbeta_tof->SetBranchAddress("delta_t_decay",&decay_time);
		for(Long64_t index=0;index<tbeta_tof->GetEntriesFast();index++){
			tbeta_tof->GetEntry(index);
			h_beta_decay->Fill(TMath::Log10((double)decay_time*1e-9));
		}
		h_beta_decay->Draw();

		c1->cd(2)->Update();
		c1->cd(2)->SetEditable(kFALSE);



	Ntof_counter = (int)eve_sele_tof.size(); // store number of tof cut by cursor for beta-tof coincidence to GLOBAL variable

	if(mr_sequence.size()>0){

		if(mr_sequence.size()==1){// only one count of tof in spectrum
			Tof_rate = 1./(Maxsweeps/40.); // rate per second;
		}
		else{
			sort(mr_sequence.begin(),mr_sequence.end());
			for(long i=0;i<(long)mr_sequence.size()-1;i++){
					Tof_rate+=mr_sequence[i+1]-mr_sequence[i];
			}
			Tof_rate = Tof_rate/(mr_sequence.size()-1)*1e-9; // average interval ns ==> s
			Tof_rate = 1/Tof_rate; // counting rate per second;
		}
	}
	else{
		if(Ntof_counter==0) cout<<"\e[1;31m"<<"No tof event found in the region specificed by cursor. Abort!"<<"\e[0m"<<endl;
		else cout<<"\e[1;31m"<<"Please check whether ALL beta files are loaded when using normal beta-tof coincident mode. Abort!"<<"\e[0m"<<endl;
				
		f_beta_tof_tree->cd();
		h_beta_decay->Write();
	
		intree->ResetBranchAddresses();            // Very important, variables address linked to tree CAN NOT be delete or destroy untill tree is deleted or 
		tbeta_tof->ResetBranchAddresses();		// reset !!!!!!!!!!!!!!!!!!!!!!!!
	
		fin->cd();
		delete mr_time;
		delete mr_sweeps_global;
		delete mr_tag;
		delete mr_channel;

		Cutvalid = false; // reset
		return;

	}

if(Cutvalid){
	//%%%%%%%%%%%%% method 1: evaluate beta background by the average time difference of the closest negative beta to the tof %%%%%%%%%%%%%%%%%%%%
	cout<<"Evaluate background beta by average interval of the closest negative beta to the TOF event"<<endl;
	if(counts_beta_bg>0){
			cout<<"couts bg beta = "<<counts_beta_bg<<endl;
			cout<<"sum time beta = "<<-beta_bg_rate<<endl;
				beta_bg_rate = 1./( -beta_bg_rate/counts_beta_bg * 1e-9);  // counts rate per second!!!!!!!!!!
				Beta_bg_rate = beta_bg_rate; // stored for globally use ; in [s]
				cout<<"Beta background rate = "<<beta_bg_rate<<" /s"<<endl;
				beta_bg_rate = beta_bg_rate * Halflife2Lifetime(RIHalflive/1e6)/1000; //Halflife2Lifetime(ms in)=> ms out
				cout<<"Beta background rate = "<<beta_bg_rate<<" /tao (lifetime)"<<endl;
	}else{
		Beta_bg_rate=0;			
		cout<<"Beta background rate = 0 / s"<<endl;			
		cout<<"Beta background rate = 0 /tao (lifetime)"<<endl;
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//%%%%%%%%%%%%% method 2: evaluate background beta by total number of negative in a given time region [halflife x time2halflife, 0] during observation time (halflife x time2halflife x Ntof) %%%%%%%%%%%%%%
	if(!RecordOnly_1st_beta_neg){
		counts_beta_bg = T_data_neg.size();
		if(counts_beta_bg ==0){
			Beta_bg_rate=0;
			cout<<endl;
			cout<<endl;
			cout<<"\e[1;31m"<<"No beta found in [- halflife x Times2Halflive , 0] region. Set Beta_bg_rate=0"<<endl;
		} 
		else{
			beta_bg_rate = counts_beta_bg / (Ntof_counter * Times2Halflive*_Halflife * 1e-9);
			Beta_bg_rate = beta_bg_rate;
			cout<<endl;
			cout<<endl;
			cout<<"*************** evaluate beta background by average beta counts in [- halflife x Times2Halflive , 0] x Ntof"<<endl;
			cout<<Form("replace Beta background rate by: Beta_bg_rate = %.6f /s \n",Beta_bg_rate)<<endl;
			beta_bg_rate = beta_bg_rate * Halflife2Lifetime(RIHalflive/1e6)/1000; //Halflife2Lifetime(ms in)=> ms out
			cout<<"Beta background rate = "<<beta_bg_rate<<" /tao (lifetime)"<<endl;

		}

	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double tem_par[2] ={(double)Ntof_counter, detector_effi}; //par0 = Ntof; par1 = efficiency   ; for FOM
cout<<"Ntof = "<<tem_par[0]<<endl;

double tem_lamda_ref = 1/(Halflife2Lifetime(RIHalflive/1e6)/1000);
	double tem_par_2[2]={tem_lamda_ref,detector_effi}; //for FOM_2;

	//double vb[1]= {beta_bg_rate};
double vb_2[1] = {Beta_bg_rate};
//	double bestwindow = BestDeltaT(vb,tem_par); // best window for coincidnece in the unit of life time!!!!
double bestwindow = BestDeltaT_2(vb_2,tem_par_2);// best window for coincidnece in the unit of second!!!!
	//bestwindow = bestwindow / TMath::Log(2);  // convert to the unit of half life
bestwindow = bestwindow /(RIHalflive*1e-9); //convert to the unit of half life
	cout<<"Best window for coincidence(in the unit of half life): "<<"\e[5m\e[1;37m"<<bestwindow<<"\e[0m"<<endl;
	
}
else{
	cout<<"\e[5m\e[1;31m"<<"Cut is not available yet OR NO background counts found!!!! No beta background rate display!!"<<"\e[0m"<<endl;
}


		//delete tempp;
		f_beta_tof_tree->cd();
	//	tbeta_tof->Write();
		h_beta_decay->Write();
	//	f_beta_tof_tree->Write();

		intree->ResetBranchAddresses();            // Very important, variables address linked to tree CAN NOT be delete or destroy untill tree is deleted or 
		tbeta_tof->ResetBranchAddresses();		// reset !!!!!!!!!!!!!!!!!!!!!!!!
	
		delete mr_time;
		delete mr_sweeps_global;
		delete mr_tag;
		fin->cd();

	Cutvalid = false; // reset


	cout<<"\e[1;33m Option: using PrintSelectedTofInfo() to print Informaton of selected Tof. Allowing to print All and set offset by setting arguements \e[0m"<<endl;

}

void PrintSelectedTofInfo(bool PrintAll=false, unsigned long offsetNum=0){
	cout<<MRTOF_hit_Info.size()<<" records of selected TOF"<<endl<<endl;

	for(unsigned long index=offsetNum;index<MRTOF_hit_Info.size();index++){
			if( (index - offsetNum)%30==0 && !PrintAll){
				cout<<"Show next 30 records or not? ('y' for continue; 'n' for stop)"<<endl;
				char inputc[3]={'\0'};
				while(1){
					cin>>inputc;
					if(inputc[0]!='y' && inputc[0]!='n') cout<<"input 'y' or 'n' "<<endl;
					else break;
				}
				if(inputc[0] == 'n') return;
			}

			cout<<MRTOF_hit_Info[index]<<endl;
	}

}


void Beta_Beta_counter(TCutG** usecut=NULL){

	bool usetimecut = false, Isnoise =false;
	int counter=0;

	for(unsigned long index=0;index< Match_hit.size();index++){

		usetimecut = false;
		Isnoise = false;

			if(usecut!=NULL){ // use cut

				if(usecut[0]!=NULL){ // cut available
					
					for(int ic=1;ic<5;ic++){
						if(usecut[ic]!=NULL){
							usetimecut = true;
							Isnoise = Isnoise || usecut[ic]->IsInside((double)Match_hit[index].B1_time_relative,(double)Match_hit[index].B2_time_relative);  // in any of the time noise cut
						} 
					}
				

					if(usetimecut)	Isnoise = Isnoise && usecut[0]->IsInside(ADC2keV_1(Match_hit[index].B1_adc),ADC2keV_2(Match_hit[index].B2_adc)); // use timecut:  noise should be in time cut and energy cut at the same time
					else Isnoise = usecut[0]->IsInside(ADC2keV_1(Match_hit[index].B1_adc),ADC2keV_2(Match_hit[index].B2_adc));   // only energy cut;

					if(!Isnoise){ counter++;	}

				}
				else{ counter++; }

			}
			else{ counter++; }

	}

	cout<<"Total beta-beta counts = "<<counter<<endl;
	cout<<"Background beta rate = "<<counter/(Maxsweeps/40.0)<<"/s"<<"  Total running time= "<<Maxsweeps/40.0<<"/s"<<endl;

}



string CutCondition(TCutG** _cut=NULL){
	if(_cut==NULL) return ""; // no cut input
	string condition_e_s="";
	string condition_t_s=" ";
	if(_cut[0]!=NULL) condition_e_s = string(_cut[0]->GetName());
	else return condition_e_s;  // empty cut

	for(int ic=1;ic<5;ic++){
		if(_cut[ic]!=NULL) condition_t_s += string(_cut[ic]->GetName()) + " || " ;
	}
	
	if(condition_t_s != " "){ // has E cut and time cut
		 condition_t_s += "false";
		 condition_t_s.erase(condition_t_s.begin()+condition_t_s.find(" || false"),condition_t_s.end());
		 condition_e_s += " && (";  // when using " && " event in energy cut and time cut at the same time will be considered as noise
		 condition_e_s = condition_e_s + condition_t_s +")";
		 
	}
	
	return condition_e_s = "!( "  + condition_e_s + ")";

}



double Decay_curve(double* miu, double* para){// miu = lg(decay time)
	// para[0] => N0  ; para[1] => miu0 = lg (livetime)

	double lamda_dk = 1/TMath::Power(10,para[1]);
	double ten_power_miu = TMath::Power(10,miu[0]);

	return para[0]*lamda_dk*ten_power_miu*TMath::Exp(-lamda_dk*ten_power_miu);
}


/*double Decay_curve(double* miu, double* para){
				// para[0] => N0  ; para[1] => miu0 = log (livetime)
	miu[0] = TMath::Log10(miu[0]);
	double result = (miu[0] - para[1])*TMath::Ln10() - TMath::Exp((miu[0] - para[1])*TMath::Ln10());
	return para[0]*TMath::Ln10()*TMath::Exp(result);

}*/

/*double Decay_curve(double* miu, double* para){
				// para[0] => N0  ; para[1] => miu0 = log (livetime)
	//miu[0] = TMath::Log10(miu[0]);
	double result = (miu[0] - para[1])*2.3 - TMath::Exp((miu[0] - para[1])*2.3);
	return para[0]*2.3*TMath::Exp(result);

}*/

TFile* f_beta_tof_history=NULL;
TTree* tbeta_tof_history=NULL;
string SeletIonCommand(string command);
TTree* CopyTreeFromFile(string _selection, TFile** F2keeptree);  // select and copytree from beta-tof_tree_raw.root with multi ion records
TTree* RefineAllTree(TTree* _treein, TFile** _fileout);
vector<double> tof_unbinned;  // store tof of ion selected for unbinned fitting
vector<double> dk_time_s_cp;  // store dk time for 2D histogram
vector<Long64_t> sweep_glo_coin; // sweeps of coin event;



vector<double> dk_time_test;
long Nbeta_test = 10;

void GenerateDKtime(vector<double> &_dk_time_test, long _ncounts){
	TRandom3* ran = new TRandom3((unsigned int)time(0));
	double mylifetime = Halflife2Lifetime(RIHalflive/1000000.0)/1000; //==>> lifetime in second

	_dk_time_test.clear();
	_dk_time_test.shrink_to_fit();

	for(long i=0;i<_ncounts;i++){
		_dk_time_test.push_back(ran->Exp(mylifetime));
	}

}



bool IsNewTOF(double _intof, Long64_t _insweep_global, bool Reset = false){ // Judge the input TOF is a new event or not; Reset = true to clear the memory
	static vector<double> tof_history;
	static vector<Long64_t> sweep_history;
	if(Reset){
		tof_history.clear(); tof_history.shrink_to_fit();
		sweep_history.clear(); sweep_history.shrink_to_fit();
		return false;
	}

	if(tof_history.size()==0){
			tof_history.push_back(_intof);  
			sweep_history.push_back(_insweep_global);
	 		return true;
	}

	bool IsNew = true;
	for(unsigned long index=0;index<tof_history.size();index++){
		if(tof_history[index] == _intof && sweep_history[index] == _insweep_global){ IsNew = false; break;} // found repeat event
	}

	if(IsNew){tof_history.push_back(_intof);  sweep_history.push_back(_insweep_global); } 

	return IsNew;

}



void ShowDecayHisto(TCutG** usecut=NULL, RejectionList* _rejectlist=NULL,string SelectIonList="IonA;index IonB;index", 
		string GeneralSelection="", string DrawIonName="",bool refinalltree=false){


		//*********** test with fake beta decay (simulation only)******************

		//	GenerateDKtime(dk_time_test,Nbeta_test);

		//******************************************************




	if(tbeta_tof == NULL){
		cout<<"\e[1;33m"<<"tbeta_tof is empty. Abort!!!!"<<"\e[0m"<<endl;
		return;
	}

	if(f_beta_tof_tree->IsOpen()){f_beta_tof_tree->cd();}


	TTree* tbeta_tof_backup_tem = tbeta_tof; // this pointer is temporarily borrowed by tbeta_tof_history

	if(SelectIonList!="IonA;index IonB;index"){ // load tree form history record
		string select_s = SeletIonCommand(SelectIonList);
		
		if(f_beta_tof_history!=NULL && f_beta_tof_history->IsOpen()) f_beta_tof_history->Close();
		TFile** fptr= &f_beta_tof_history;
		tbeta_tof = CopyTreeFromFile(select_s,fptr);
		if(tbeta_tof==NULL){
			 cout<<"\e[1;33m"<<"fail to get tree form beta-tof_tree_raw.root"<<"\e[0m"<<endl;
			 tbeta_tof = tbeta_tof_backup_tem;
			 return;
		}
		else{
			tbeta_tof_history=tbeta_tof;
			cout<<"get tbeta_tof to tbeta_tof_history from beta-tof_tree_raw.root"<<endl;
		}

		if(DrawIonName==""){
			 cout<<"\e[1;32m"<<"Specify the IonName to Draw"<<"\e[0m"<<endl;
			 tbeta_tof = tbeta_tof_backup_tem;
			 return;
		}

		TFile* f_refine_all_final;
		TFile** fout_temp = &f_refine_all_final;
		
		if(refinalltree){
			tbeta_tof=RefineAllTree(tbeta_tof, fout_temp);
			tbeta_tof_history=tbeta_tof;
			f_beta_tof_history->Close();
			f_beta_tof_history = f_refine_all_final;
		}

	}


	if(f_beta_tof_tree->IsOpen()){f_beta_tof_tree->cd();}

		tbeta_tof = tbeta_tof->CopyTree(GeneralSelection.c_str());
	//TTree* tbeta_tof_selected_history = tbeta_tof_history->CloneTree(); return;
	//tbeta_tof_selected_history->ResetBranchAddresses();
		//tbeta_tof = tbeta_tof_selected_history->CopyTree("");


	double h_uplimit = Times2Halflive * RIHalflive * 1e-9;  // from [ns] to [s]
	
	h_uplimit = (int)TMath::Log10(h_uplimit)+2;

	h_uplimit = TMath::Power(10,h_uplimit);

	if(h_beta_decay!=NULL)h_beta_decay->Reset(); // created in findcoincidence function
	else{cout<<"\e[1;33m"<<"beta decay plot not Exist!!! Press 'c' to activate it firstly"<<"\e[0m"<<endl;return;}

	if(h_zoom_x!=NULL)h_zoom_x->Reset();
	else{cout<<"\e[1;31m"<<"Please recreate h_zoom_x first by using \"space\" key!!!"<<"\e[0m"<<endl; return;}

	if(histo_beta_EkeV_CutAndRej!=NULL)histo_beta_EkeV_CutAndRej->Reset();


	T_data.clear();// T_data global at chance.cc ==>store decay time
	T_data.shrink_to_fit();
//T_data_neg.clear();
//T_data_neg.shrink_to_fit();
	vector<double> Tof_without_coin;
	IsNewTOF(0,0,true); // reset the tof history inside the function

	//************** fill time spectrum of decay with cut **************
	bool Isnoise=false;
	bool usetimecut=false;
	bool IsReject=false;
	bool IsThisIon=false;
	double Ncounts=0;
	double dk_mean_s = 0;
	Long64_t dk_time_gclock_ns=0;
	double dk_time_s=0;
	int adc1=0; // identify tof without coincident beta in the tree ==> adc1==-1000 case
	double E1=0;
	double E2=0;
	double tof=0;
	Long64_t time_gclock_mr=0;
	Long64_t time1=0;
	Long64_t time1_relative=0;
	Long64_t time2=0;
	Long64_t time2_relative=0;
	Long64_t mr_Bingo_sweeps_global=0;
	string* IonName_ptr=NULL;

	tbeta_tof->SetBranchAddress("Ntof_counts",&Ntof_counter);
	tbeta_tof->SetBranchAddress("delta_t_decay",&dk_time_gclock_ns);
	tbeta_tof->SetBranchAddress("adc1",&adc1);
	tbeta_tof->SetBranchAddress("EkeV_1",&E1);
	tbeta_tof->SetBranchAddress("EkeV_2",&E2);
	tbeta_tof->SetBranchAddress("time_gclock_mr",&time_gclock_mr);
	tbeta_tof->SetBranchAddress("time_mr",&tof);
	tbeta_tof->SetBranchAddress("time1",&time1);
	tbeta_tof->SetBranchAddress("time1_relative",&time1_relative);
	tbeta_tof->SetBranchAddress("time2",&time2);
	tbeta_tof->SetBranchAddress("time2_relative",&time2_relative);
	tbeta_tof->SetBranchAddress("sweeps_global_mr",&mr_Bingo_sweeps_global);

	string treename = tbeta_tof->GetName();
	if(treename.find("refine")!=string::npos) tbeta_tof->SetBranchAddress("IonName",&IonName_ptr);

	tof_unbinned.clear(); tof_unbinned.shrink_to_fit();  // store tof of ion selected
	dk_time_s_cp.clear(); dk_time_s_cp.shrink_to_fit(); // store dk time of ion selected
	sweep_glo_coin.clear(); sweep_glo_coin.shrink_to_fit(); // used together with tof_unbinned

	Long64_t entries = tbeta_tof->GetEntriesFast();



//	for(long i=0;i<Nbeta_test;i++){// for simulation


	for(Long64_t index=0;index<entries;index++){
		tbeta_tof->GetEntry(index);

		if(_rejectlist!=NULL){ // use rejectlist
			IsReject = _rejectlist->IsInside(time_gclock_mr,time1,time2);
		}

		if(IonName_ptr==NULL){	IsThisIon=true;	}
		else{
			if(DrawIonName==(*IonName_ptr))IsThisIon=true;
			else IsThisIon=false;
		}


		if(usecut!=NULL){ // use cut
				if(usecut[0]!=NULL){ // cut available
					for(int ic=1;ic<5;ic++){
						if(usecut[ic]!=NULL){
							usetimecut = true;
							Isnoise = Isnoise || usecut[ic]->IsInside((double)time1_relative,(double)time2_relative);  // in any of the time noise cut
						} 
					}

					if(usetimecut)	Isnoise = Isnoise && usecut[0]->IsInside(E1,E2); // use timecut:  noise should be in time cut and energy cut at the same time
					else Isnoise = usecut[0]->IsInside(E1,E2);   // only energy cut;

					if(!Isnoise && !IsReject && IsThisIon){
						if(adc1 == -1000) Tof_without_coin.push_back(tof);
						else{
							dk_time_s = (double)(dk_time_gclock_ns * 1e-9);
							h_beta_decay->Fill(TMath::Log10(dk_time_s));
							if( IsNewTOF(tof,mr_Bingo_sweeps_global) ){  // block repeated tof
								h_zoom_x->Fill(tof);
								tof_unbinned.push_back(tof);
								sweep_glo_coin.push_back(mr_Bingo_sweeps_global);
							}
							dk_time_s_cp.push_back(dk_time_s);
							if(histo_beta_EkeV_CutAndRej!=NULL)histo_beta_EkeV_CutAndRej->Fill(E1,E2);
							if(dk_time_s>0) dk_mean_s+= dk_time_s;
							Ncounts++;
							T_data.push_back(dk_time_s);
						}
						
					}

					Isnoise = false; // clear
				}
				else{ // empty cut
					if(!IsReject && IsThisIon){
						if(adc1 == -1000) Tof_without_coin.push_back(tof);
						else{
							dk_time_s = (double)(dk_time_gclock_ns * 1e-9);
							h_beta_decay->Fill(TMath::Log10(dk_time_s));
							if( IsNewTOF(tof,mr_Bingo_sweeps_global) ){
								h_zoom_x->Fill(tof);
								tof_unbinned.push_back(tof);
								sweep_glo_coin.push_back(mr_Bingo_sweeps_global);
							}
							dk_time_s_cp.push_back(dk_time_s);
							if(histo_beta_EkeV_CutAndRej!=NULL)histo_beta_EkeV_CutAndRej->Fill(E1,E2);
							if(dk_time_s>0) dk_mean_s+= dk_time_s; 
							Ncounts++;
							T_data.push_back(dk_time_s);
						}
					}
				}  
		}
		else{ // no cut input
			if(!IsReject && IsThisIon){
					if(adc1 == -1000) Tof_without_coin.push_back(tof);
					else{
						dk_time_s = (double)(dk_time_gclock_ns * 1e-9);
		//dk_time_s = dk_time_test[i];// for simulation
						h_beta_decay->Fill(TMath::Log10(dk_time_s));
						if( IsNewTOF(tof,mr_Bingo_sweeps_global) ){
							h_zoom_x->Fill(tof);
							tof_unbinned.push_back(tof);
							sweep_glo_coin.push_back(mr_Bingo_sweeps_global);
						}
						dk_time_s_cp.push_back(dk_time_s);
		//sweep_glo_Beta_Beta.()
						if(histo_beta_EkeV_CutAndRej!=NULL)histo_beta_EkeV_CutAndRej->Fill(E1,E2);
						if(dk_time_s>0) dk_mean_s+= dk_time_s; 
						Ncounts++;
						T_data.push_back(dk_time_s);
					}
			}
		}  
	}

	tbeta_tof->ResetBranchAddresses();

	c_beta_beta->cd(4);//->SetLogx();
	h_beta_decay->Draw();
	h_beta_decay->GetXaxis()->SetNdivisions(005);
	for(int i=1; i<=6;i++){
		h_beta_decay->GetXaxis()->ChangeLabel(i,-1,-1,-1,-1,-1,Form("10^{%d}",i-4));
	}

	h_beta_decay->GetXaxis()->SetLabelOffset(0.015);


	// get the count of tof event (total tof counts)
//Ntof_counter = (int)tof_unbinned.size() + (int)Tof_without_coin.size();

	// evaluate the beta background, using negative decay counts
	if(RecordOnly_1st_beta_neg){// single negative decay counts for each tof
		double dk_sum_s=0;
		int tem_count=0;
		for(unsigned long index=0; index< T_data.size();index++){
			if(T_data[index]<0){
				 dk_sum_s += -T_data[index];
				 tem_count++;
			}
		}

		if(tem_count>0) Beta_bg_rate = tem_count/dk_sum_s;
		else Beta_bg_rate = 0;

		cout<<"Beta background rate = "<<Beta_bg_rate<<" /s"<<endl;
		cout<<"Beta background rate = "<<Beta_bg_rate * Halflife2Lifetime(RIHalflive/1e6)/1000<<" /tao (lifetime)"<<endl;//Halflife2Lifetime(ms in)=> ms out

	}else{// multiple negative decay counts for each tof
			int tem_count=0;
			for(unsigned long index=0; index< T_data.size();index++){
				if(T_data[index]<0) tem_count++;
			}

			if(tem_count==0) Beta_bg_rate = 0;
			else{
				Beta_bg_rate= tem_count / (Ntof_counter * Times2Halflive * RIHalflive* 1e-9);
			}

			cout<<"Beta background rate = "<<Beta_bg_rate<<" /s"<<endl;
			cout<<"Beta background rate = "<<Beta_bg_rate * Halflife2Lifetime(RIHalflive/1e6)/1000<<" /tao (lifetime)"<<endl;//Halflife2Lifetime(ms in)=> ms out
	}



double setpar[4]={Beta_bg_rate,Times2Halflive*RIHalflive*1e-9,detector_effi,(double)Ntof_counter};// vb in [s]; window in [s] ; efficiency; Ntof

ExtractLifetime(setpar, RIHalflive, RecordOnly_1st_beta);

	c_beta_beta->cd(4);

	//double Ncounts = h_beta_decay->Integral(1,h_beta_decay->GetNbinsX());
cout<<"Total Ncounts_b_tof(negative + positive) = "<<Ncounts<<endl;

	Ncounts=0;
	for(unsigned long index=0; index< T_data.size();index++){ if(T_data[index]>0) Ncounts++;	}

	if(RecordOnly_1st_beta) cout<<"Ncounts_b_tof (positive) at \"RecordOnly_1st_beta\" Mode = "<<Ncounts<<endl;		//old method for calculating lifetime using only the first positive beta
	else cout<<"Ncounts_b_tof (positive) at \"All beta in window\" Mode = "<<Ncounts<<endl;


	if(dk_curve_ref!=NULL) delete dk_curve_ref;
	dk_curve_ref = new TF1("dk_curve_ref",Decay_curve,-3,2,2);
	//dk_curve_ref = new TF1("dk_curve_ref",Decay_curve,0.01,h_uplimit,2);
	//dk_curve_ref->SetParameter(0,Ncounts);
			//*******************
			//dk_curve_ref->SetParameter(0,h_beta_decay->GetBinContent(h_beta_decay->GetMaximumBin())); // method 1
			//*******************

			//scale based on the sum of the y value of the curve
			double scale_factor=0;
			dk_curve_ref->SetParameter(0,1);
			for(int i=1;i<=h_beta_decay->GetNbinsX();i++){
				scale_factor+= dk_curve_ref->Eval(h_beta_decay->GetBinCenter(i));
			}

			scale_factor = Ncounts/scale_factor;
			dk_curve_ref->SetParameter(0,scale_factor);

	double miu0 = TMath::Log10((RIHalflive * 1e-9) / TMath::Log(2));
	dk_curve_ref->SetParameter(1,miu0);
	/*double Ncounts_expected = Ncounts/(1-TMath::Power(0.5,Times2Halflive)); // method 2
	dk_curve_ref->SetParameter(0,Ncounts_expected);

	double miu0 = TMath::Log10((RIHalflive * 1e-9) / TMath::Log(2));
	dk_curve_ref->SetParameter(1,miu0);

	double scalefactor=0;
	double N_curve=0;
	for(int ibin=1;ibin<=h_beta_decay->GetNbinsX();ibin++){
			N_curve+=dk_curve_ref->Eval(h_beta_decay->GetBinCenter(ibin));
	}

	scalefactor = N_curve/Ncounts_expected;
	dk_curve_ref->SetParameter(0,Ncounts_expected/scalefactor);*/

	dk_curve_ref->SetLineColor(kGreen);
	dk_curve_ref->Draw("same");


	//****** life time measurement*************
	dk_mean_s = dk_mean_s / Ncounts;  // life time

	if(dk_curve_m!=NULL) delete dk_curve_m;
	dk_curve_m = new TF1("dk_curve_m",Decay_curve,-3,2,2);
	//dk_curve_m = new TF1("dk_curve_m",Decay_curve,0.01,h_uplimit,2);
	//dk_curve_ref->SetParameter(0,Ncounts);
			//dk_curve_m->SetParameter(0,h_beta_decay->GetBinContent(h_beta_decay->GetMaximumBin())); // method 1
	//dk_curve_m->SetParameter(0,Ncounts_expected/scalefactor); // method 2
			dk_curve_m->SetParameter(0,scale_factor);

	double miu0_m = TMath::Log10(dk_mean_s);
	dk_curve_m->SetParameter(1,miu0_m);
	dk_curve_m->SetLineColor(kRed);
	dk_curve_m->Draw("same");

	dk_mean_s *= TMath::Log(2); // to half lifetime
	char ref_out[30]={'\0'};
	char m_out[30]={'\0'};

	if(Ncounts>=20){
		cout<<"mean T1/2 = "<<dk_mean_s<<" ("<<dk_mean_s /TMath::Sqrt(Ncounts)<<")"<<endl;
		sprintf(m_out,"T1/2_m = %.2f(%.2f)",dk_mean_s,dk_mean_s /TMath::Sqrt(Ncounts));
	}
	else if(Ncounts>2){
		cout<<"mean T1/2 = "<<dk_mean_s<<" (-"<<dk_mean_s /(1+TMath::Sqrt(Ncounts))<<" , +"<<dk_mean_s /(TMath::Sqrt(Ncounts)-1)<<")"<<endl;
		sprintf(m_out,"T1/2_m = %.2f(-%.2f,+%.2f)",dk_mean_s,dk_mean_s /(1+TMath::Sqrt(Ncounts)),dk_mean_s /(TMath::Sqrt(Ncounts)-1));
	}
	else if(Ncounts==2){
		cout<<"mean T1/2 = "<<dk_mean_s<<" (-"<<dk_mean_s *(1-0.606)<<" , +"<<dk_mean_s *(2.82-1)<<")"<<endl;
		sprintf(m_out,"T1/2_m = %.2f(-%.2f,+%.2f)",dk_mean_s,dk_mean_s *(1-0.606),dk_mean_s *(2.82-1));
	}
	else{
		cout<<"mean T1/2 = "<<dk_mean_s<<" (-"<<dk_mean_s *(1-0.543)<<" , +"<<dk_mean_s *(5.79-1)<<")"<<endl;
		sprintf(m_out,"T1/2_m = %.2f(-%.2f,+%.2f)",dk_mean_s,dk_mean_s *(1-0.543),dk_mean_s *(5.79-1));
	}


	double lifetime = (RIHalflive * 1e-9) / TMath::Log(2);  // lifetime in second

	double Expect_halflife = ExpectHalflife(lifetime, Times2Halflive * RIHalflive * 1e-9);
	
	printf("Simulate: Calculated half life by average decay time(without background) within limits(%.4f x halflife) if infinite statistic: %.4f\n",Times2Halflive,Expect_halflife);

	sprintf(ref_out,"T1/2_ref = %.2f",RIHalflive * 1e-9);

	if(dk_legend!=NULL) delete dk_legend;

	dk_legend = new TLegend(0.12,0.7,0.4,0.89);
	dk_legend->AddEntry(dk_curve_ref,ref_out,"l");
	dk_legend->AddEntry(dk_curve_m,m_out,"l");
	dk_legend->Draw();

	if(c_eject!=NULL && c_eject->GetCanvasImp()!=NULL) c_eject->cd(5); // window is valid
	else{
			ShowEjectionBeta(EJE0,EJE1); // define in Beta_Beta_coin.cc
			c_eject->cd(5);
	}

	gStyle->SetPalette(kBird);

	f_beta_tree->cd();
	tbeta_tof->Draw("EkeV_2:EkeV_1>>histo_beta_EkeV_withcut",CutCondition(usecut).c_str(),"colz");
	if(histo_beta_EkeV_CutAndRej!=NULL){
		c_eject->cd(6);
		histo_beta_EkeV_CutAndRej->Draw("colz");
	}

	c1->cd(2)->SetEditable(kTRUE);
	fin->cd();
	h_zoom_x->Draw();
	tof_unbinned.shrink_to_fit();
	c1->cd(2)->SetEditable(kFALSE);

	string this_tree_name = tbeta_tof->GetName();

	if(this_tree_name != "tbeta_tof") active_tree_name = "tbeta_tof_history";
	else{active_tree_name = "tbeta_tof";}

	tbeta_tof = tbeta_tof_backup_tem;

	GetEveNoCoin(eve_sele_tof,eve_sele_glo_sweep,tof_unbinned,sweep_glo_coin); // generate list of event without coin


}


double* xedgeL_2D=NULL;
double* yedgeL=NULL;

void ShowBeta_Tof_2D(bool _makecut=false){

	double h_uplimit = Times2Halflive * RIHalflive * 1e-9;  // from [ns] to [s]
	h_uplimit = (int)TMath::Log10(h_uplimit)+2;
	double h_lowlimit=-2;

	if(h_beta_tof_2D!=NULL) {delete h_beta_tof_2D; h_beta_tof_2D=NULL;}
	
	if(xedgeL_2D!=NULL){ delete xedgeL_2D;xedgeL_2D=NULL;}

				// various bin width for log scale X
				//**************************
				/*int numbins = (int)((h_uplimit-h_lowlimit)*9);
					
				xedgeL_2D = new double[numbins+1];
				int xindex=0;
				for(int order = h_lowlimit;order<h_uplimit; order++){
					int digit=1;
					while(digit<10){
						xedgeL_2D[xindex++] = (digit++)*TMath::Power(10,order);
					}
				}

				h_uplimit = TMath::Power(10,h_uplimit); // right edge of last bin
				xedgeL_2D[xindex] = h_uplimit;*/
				//****************************

	xedgeL_2D = new double[51];
	for(int i=0;i<51;i++){
		xedgeL_2D[i] = -3+i*0.1;
	}


	if(yedgeL!=NULL) delete yedgeL;
	yedgeL = new double[h_zoom_x->GetNbinsX()+1];
	for(int ibin=1;ibin<=h_zoom_x->GetNbinsX();ibin++){
		yedgeL[ibin-1] = h_zoom_x->GetBinLowEdge(ibin); // left edge of each bin
	}
	yedgeL[h_zoom_x->GetNbinsX()] = yedgeL[h_zoom_x->GetNbinsX()-1] + h_zoom_x->GetBinWidth(1);


	fin->cd();
	h_beta_tof_2D = new TH2D("h_beta_tof_2D","",h_zoom_x->GetNbinsX(),yedgeL,50,xedgeL_2D); // X: tof; Y: decay time
	//h_beta_tof_2D = new TH2D("h_beta_tof_2D","",h_zoom_x->GetNbinsX(),yedgeL,numbins,xedgeL_2D); // X: tof; Y: decay time
	h_beta_tof_2D->SetTitle("delta_t_decay:time_mr;Tof [ns]; Decay time [s]");
	h_beta_tof_2D->GetXaxis()->CenterTitle();
	h_beta_tof_2D->GetYaxis()->CenterTitle();

	for(long index=0;index<(long)tof_unbinned.size();index++){
		h_beta_tof_2D->Fill(tof_unbinned[index],TMath::Log10(dk_time_s_cp[index]));
	}


	h_beta_tof_2D->GetYaxis()->SetNdivisions(005);
	for(int i=1; i<=6;i++){
		h_beta_tof_2D->GetYaxis()->ChangeLabel(i,-1,-1,-1,-1,-1,Form("10^{%d}",i-4));
	}

	h_beta_tof_2D->GetYaxis()->SetLabelOffset(0.015);



		TCutG* my2d=NULL;
		static TCanvas* C_beta_tof_2D=NULL;
		static TCutG* my2d_cp=NULL;
		if(C_beta_tof_2D!=NULL && C_beta_tof_2D->GetCanvasImp()!=NULL) delete C_beta_tof_2D; // window is closed
		C_beta_tof_2D = new TCanvas("C_beta_tof_2D","C_beta_tof_2D",700,500); 
		h_beta_tof_2D->Draw("colz");
		//C_beta_tof_2D->SetLogy();
		C_beta_tof_2D->Modified();
		C_beta_tof_2D->Update();
		if(_makecut){
			my2d = (TCutG*)C_beta_tof_2D->WaitPrimitive("CUTG");
			int Npoints = my2d->GetN();
			my2d->SetPoint(0,my2d->GetPointX(Npoints-1),my2d->GetPointY(Npoints-1));
			if(gROOT->GetListOfSpecials()->FindObject("my2d") != NULL){gROOT->GetListOfSpecials()->Remove((TObject *)my2d_cp);}
			if(my2d_cp!=NULL) delete my2d_cp;
			my2d_cp = new TCutG(*my2d);
			for(int i=0;i<Npoints;i++){
				//my2d_cp->SetPointY(i,my2d_cp->GetPointY(i)*1E9); // convert the Y value of cut to the unit of delta_t_decay [ns]
				my2d_cp->SetPointY(i,TMath::Power(10,my2d_cp->GetPointY(i))*1E9);

			}
			my2d_cp->SetName("my2d");
			gROOT->GetListOfSpecials()->Add((TObject *)my2d_cp);
		}
}

double* gtof=NULL; // test for tof veto
double* gsweesps_globel=NULL; // test for tof veto


void UnbinnedFit(TF1* func_handle,string _treename, TH1D* outhisto, double _fitL, double _fitR){

	if(outhisto==NULL){cout<<"\e[1;31m"<<"h_zoom_x is not exist; abort!!!"<<"\e[0m"<<endl; return;}

	long Numdata=0;
	double* _intof=NULL;	

	//*********** Get data for fit ***************************
	if(_treename == "intree"){ // tree without beta correlation
		/*intree->SetEstimate(-1);
		string selected_s="";
		string intree_name = intree->GetName();
		if(intree_name == "tree"){ // no drift corrected
			selected_s = Form("time >=%.4f && time<=%.4f",_fitL,_fitR);
			intree->Draw("time:sweeps_global",selected_s.c_str(),"goff");
			_intof = intree->GetV1(); gtof=_intof; gsweesps_globel = intree->GetV2();
			Numdata = intree->GetSelectedRows(); cout<<Numdata<<endl;
		}
		else{
			selected_s = Form("timec >=%.4f && timec<=%.4f",_fitL,_fitR);
			intree->Draw("timec:sweeps_global",selected_s.c_str(),"goff");
			_intof = intree->GetV1(); gtof=_intof; gsweesps_globel = intree->GetV2();
			Numdata = intree->GetSelectedRows(); cout<<Numdata<<endl;

		}*/

		if(Tof_array_hzoomx_entries==0){
				cout<<"\e[1;31m"<<"No TOF found!! Unbinned fitting abort!!"<<"\e[0m"<<endl;
				return;			
		}

		for(long long index=0;index<Tof_array_hzoomx_entries;index++){
			if(Tof_array_hzoomx[index]>=_fitL && Tof_array_hzoomx[index]<=_fitR) Numdata++;
		}

		if(Numdata==0){
				cout<<"\e[1;31m"<<"No TOF selected!! Unbinned fitting abort!!"<<"\e[0m"<<endl;
				return;			
		}


		_intof = new double[Numdata];
		long long index_tem=0;

		for(long long index=0;index<Tof_array_hzoomx_entries;index++){
			if(Tof_array_hzoomx[index]>=_fitL && Tof_array_hzoomx[index]<=_fitR) _intof[index_tem++] = Tof_array_hzoomx[index];
		}



	}
	else{// tof with beta correlation
			Numdata = tof_unbinned.size();
			if(Numdata == 0){
				cout<<"\e[1;31m"<<"No beta-correlated TOF found!! Unbinned fitting abort!!"<<"\e[0m"<<endl;
				return;
			}
			_intof = new double[Numdata];
			for(long index=0;index<Numdata; index++){_intof[index] = tof_unbinned[index];}

	}

	//if(Numdata==0){cout<<"Error!!!! No data input!!!"<<endl; return;}

	ROOT::Fit::DataRange range(_fitL,_fitR);
    ROOT::Fit::UnBinData data((unsigned int)Numdata, _intof, range);
    cout<<endl;
    cout <<"\e[1;32m"<<"ncounts in datasample for fit = " << Numdata <<"\e[0m"<< endl;


    ROOT::Math::WrappedMultiTF1 fitFunction( *func_handle, 1 );
    ROOT::Fit::Fitter fitter;
    fitter.SetFunction( fitFunction, false);

    int Npars = func_handle->GetNpar();
    double* Pars = func_handle->GetParameters();
    fitter.Config().SetParamsSettings(Npars,Pars);
  //  printf("original par0=%.4f, par1=%.4f\n",Pars[0],Pars[1]);

  //  cout<<"npar="<<fitter.Config().NPar()<<endl;
 //   printf("par0 =%.4f  par1=%.4f\n ",fitter.Config().ParSettings(0).Value(),fitter.Config().ParSettings(1).Value());
    //cout<<"ct= "<<fitter.Config().ParSettings(1).Value()<<endl;


    for(int ipar=0;ipar<Npars; ipar++){
    	if(func_handle->GetParError(ipar)==0){
    		fitter.Config().ParSettings(ipar).Fix();
    	}
    	else{
    		double parvalue = func_handle->GetParameter(ipar);
    		if(parvalue>_fitL && parvalue<_fitR){
    			fitter.Config().ParSettings(ipar).SetLimits(_fitL,_fitR);
    		}
    		else{fitter.Config().ParSettings(ipar).SetLimits(0,parvalue*10); } //!!!!! very important!!!! otherwise fitting fail
    		
    	}
    }

	fitter.Config().SetUpdateAfterFit();
	fitter.LikelihoodFit(data,true);
    TFitResult r=fitter.Result();
    r.Print();
    cout<<endl;
    func_handle->SetParameters(r.Parameters().data());
    func_handle->SetParErrors(r.Errors ().data());


    double sum=0;
    for(int nbin=outhisto->FindBin(_fitL);nbin<=outhisto->FindBin(_fitR);nbin++){
    	sum+=func_handle->Eval(outhisto->GetBinCenter(nbin));
    }

    //&********** Adjust a proper amplitude  of fititng curve for display ***************
    string funcname = func_handle->GetName();
    if(funcname == "fsample"){
    	for(int ip=0;ip<fs->NumOfPeaks;ip++){
    		double ampnow = func_handle->GetParameter(2*ip);
    		func_handle->SetParameter(2*ip,ampnow/(sum/Numdata));
    		tof_x_cento[ip] = func_handle->GetParameter(2*ip+1);
    		tof_x_cento_err[ip] = func_handle->GetParError(2*ip+1);
    		printf("cento_%d: %.4f(%.4f)\n",ip+1,tof_x_cento[ip],tof_x_cento_err[ip]);
    	}
    }
    else if(funcname == "fextend_N"){
    	fext->UpdateParsFromTF();
    	int para_offset=0;
    	double ampnow=0;
    	for(int ip=0;ip<fext->GetNumberOfPeaks();ip++){
    		para_offset = fext->Paras_offset_cal(ip,fext->GetMainPeakIndex());
    		ampnow = fext->GetAmp(ip+1);
    		func_handle->SetParameter(para_offset,ampnow/(sum/Numdata));
    		tof_x_cento[ip] = fext->GetTofCenter(ip+1);
    		tof_x_cento_err[ip] = fext->GetTofCenterErr(ip+1);
    		
    	}
    }
    else{
    	for(int ip=0;ip<NumOfPeaks;ip++){
    		int ipar = Paras_offset_cal(ip,SetOneorTwo,MainPeakIndex,ParasLock);
    		double ampnow = func_handle->GetParameter(ipar);
    		func_handle->SetParameter(ipar,ampnow/(sum/Numdata));
    		tof_x_cento[ip] = func_handle->GetParameter(ipar+1);
    		tof_x_cento_err[ip] = func_handle->GetParError(ipar+1);
    		printf("cento_%d: %.4f(%.4f)\n",ip+1,tof_x_cento[ip],tof_x_cento_err[ip]);

    	}
    	
    }
    /*
	cout<<r.Chi2()<<endl;
	cout<<r.Ndf()<<endl;
    cout<<"Show Probability > Chi2 ==> "<<TMath::Prob(r.Chi2(),r.Ndf())<<endl;*/

    func_handle->SetRange(_fitL,_fitR);
    c1->cd(2)->SetEditable(kTRUE);
    func_handle->SetLineColor(kBlue);
    func_handle->Draw("same");
    c1->cd(2)->Modified();
    c1->cd(2)->Update();
    c1->cd(2)->SetEditable(kFALSE);
  //  func_handle->SetLineColor(kRed);

	//if(_treename != "intree") delete[] _intof;

	delete[] _intof;

}





string SeletIonCommand(string command){
	   string iname[100];  // iostop name 12C
	  int ifillindex[100];  // number of this iostop in a molecule
	  int compon = 0;

	  istringstream is(command);
	  char buffer[100];
	  char *p;
	  while(is>>buffer){
		    p = strtok(buffer,";");
		    iname[compon] = p;
		    p = strtok(NULL,";");
		    if(p) ifillindex[compon] = atoi(p);
		    else  ifillindex[compon] = 1;
		    compon++;
  	   }

  	   string command_s="";
  	   for(int inum=0;inum<compon;inum++){
  	   		if(inum<compon-1){
  	   			command_s+= Form("(IonName == \"%s\" && IonFillIndex==%d) || ",iname[inum].c_str(),ifillindex[inum]);
  	   		}
  	   		else{
				command_s+= Form("(IonName == \"%s\" && IonFillIndex==%d)",iname[inum].c_str(),ifillindex[inum]);
  	   		}

  	   }

  	   cout<<command_s.c_str()<<endl;
  	   return command_s;
}



TTree* CopyTreeFromFile(string _selection, TFile** F2keeptree){

		static TTree* t_return;
		string outpath = FilePath + "beta-tof_tree_raw.root";
		if(gSystem->AccessPathName(outpath.c_str())){
			cout<<"\e[1;32m"<<"Can not open beta-tof_tree_raw.root (raw refine)!! Path or File not exist!"<<"\e[0m"<<endl;
			return NULL;
		}
		TFile* fin_beta_tof_refine = new TFile(outpath.c_str(),"UPDATE");
		if(!fin_beta_tof_refine->IsOpen()){
			cout<<"\e[1;32m"<<"Can not open beta-tof_tree_raw.root (raw refine)!! "<<"\e[0m"<<endl;
			return NULL;
		}
		TTree* gettree = (TTree*)fin_beta_tof_refine->Get("tbeta_tof_refine");
		fin_beta_tof_refine->cd();
		TTree* treecopy = gettree->CopyTree(_selection.c_str());
	cout<<"treecopy entry = "<<treecopy->GetEntriesFast()<<endl;
		vector<string> ionlist;
		vector<vector<int>> ion_list_fillindex;
		vector<int> new_species_fillindex;

		string* getIonName = NULL;  // very important! must be initialized otherwise it causes crash when getentry from tree
		int getFillIndex;
		treecopy->SetBranchAddress("IonName",&getIonName);
		treecopy->SetBranchAddress("IonFillIndex",&getFillIndex);


		for(int index=0;index<(int)treecopy->GetEntriesFast();index++){
			treecopy->GetEntry(index);

			if(ionlist.size()==0){ // initialize

				ionlist.push_back(*getIonName);
				new_species_fillindex.clear();
				new_species_fillindex.push_back(getFillIndex);
				ion_list_fillindex.push_back(new_species_fillindex);
			}

			bool isnewion=true;
			bool isnewindex=true;

			for(int row=0;row<(int)ionlist.size();row++){
				if(ionlist[row] == *getIonName){
					isnewion=false;					
					for(int column=0;column<(int)ion_list_fillindex[row].size();column++){
						if(ion_list_fillindex[row][column] == getFillIndex){
							isnewindex=false; break;
						}
					}
				}

				if(isnewindex==false)break; // old
				else if(isnewion==false && isnewindex==true){ // old ion new index
					ion_list_fillindex[row].push_back(getFillIndex);
					break;
				}

			}

			if(isnewion){
				ionlist.push_back(*getIonName);
				new_species_fillindex.clear();
				new_species_fillindex.push_back(getFillIndex);
				ion_list_fillindex.push_back(new_species_fillindex);
			}

		}


		bool valid_tree=true;
		for(int row=0;row<(int)ionlist.size();row++){
			printf("%s\t",ionlist[row].c_str());

			if(ion_list_fillindex[row].size()>1) valid_tree=false;

			for(int column=0;column<(int)ion_list_fillindex[row].size();column++){
				printf("%d\t",ion_list_fillindex[row][column]);
			}
			printf("\n");

		}

		if(!valid_tree){
				 cout<<"\e[1;33m"<<"Please apply selection option to specify unique \"FillIndex\" for each species!!! Abort!!"<<"\e[0m"<<endl;
				 return NULL;
		}
		else{
			//if((*F2keeptree) !=NULL && (*F2keeptree)->IsOpen()){(*F2keeptree)->Close(); delete (*F2keeptree);}

			outpath = FilePath + "beta-tof_tree_all_temp.root";

			TFile* fout_temp = new TFile(outpath.c_str(),"RECREATE");

			*F2keeptree = fout_temp;

			
			t_return = treecopy->CloneTree();
			t_return->Write();

			
			fin_beta_tof_refine->Close();  // must be closed;


			return t_return;

		}


}


TTree* RefineAllTree(TTree* _treein, TFile** _fileout){
	vector<match_event> repeatlist;
	match_event repeatA;
	match_event repeatB;
	vector<match_event> repeatBgroup;


		Long64_t mr_Bingo_sweeps_global=0; double mr_Bingo_time=0; Long64_t mr_Bingo_time_gclock=0;
		Long64_t time1=0; Long64_t time1_relative=0; Long64_t sweeps_gclock_1=0; int adc1=0; double EkeV_1=0;
		Long64_t time2=0; Long64_t time2_relative=0;  Long64_t sweeps_gclock_2=0; int adc2=0; double EkeV_2=0;
		int delta_t=0;
		Long64_t delta_t_decay=0;
		string* IonName_ptr = NULL;
		Long64_t B_B_gate_in_ns;
		Long64_t Halflive;
		double TimesOfHalflive;
		int IonFillIndex=1;

		_treein->SetBranchAddress("sweeps_global_mr",&mr_Bingo_sweeps_global);
		_treein->SetBranchAddress("time_mr",&mr_Bingo_time);
		_treein->SetBranchAddress("time_gclock_mr",&mr_Bingo_time_gclock);
		_treein->SetBranchAddress("sweeps_gclock_1",&sweeps_gclock_1);
		_treein->SetBranchAddress("time1",&time1);
		_treein->SetBranchAddress("time1_relative",&time1_relative);
		_treein->SetBranchAddress("adc1",&adc1);
		_treein->SetBranchAddress("EkeV_1",&EkeV_1);
		_treein->SetBranchAddress("sweeps_gclock_2",&sweeps_gclock_2);
		_treein->SetBranchAddress("time2",&time2);
		_treein->SetBranchAddress("time2_relative",&time2_relative);
		_treein->SetBranchAddress("adc2",&adc2);
		_treein->SetBranchAddress("EkeV_2",&EkeV_2);
		_treein->SetBranchAddress("delta_t",&delta_t);
		_treein->SetBranchAddress("delta_t_decay",&delta_t_decay);
		_treein->SetBranchAddress("IonName",&IonName_ptr);
		_treein->SetBranchAddress("IonFillIndex",&IonFillIndex);
		_treein->SetBranchAddress("B_B_gate_in_ns",&B_B_gate_in_ns);
		_treein->SetBranchAddress("Halflive",&Halflive);
		_treein->SetBranchAddress("TimesOfHalflive",&TimesOfHalflive);

/*
repeatA.B1_time=1;
repeatA.B2_time=2;
repeatA.B2_time_relative=3;
repeatlist.push_back(repeatA);*/

	auto Isrepeat=[&repeatlist](Long64_t mrtime,Long64_t _time1, Long64_t _time2)->bool{
		bool skipreturn=false;
		for(int list=0; list< (int)repeatlist.size(); list++){
			if(	mrtime==repeatlist[list].B1_time && _time1 == repeatlist[list].B2_time && _time2==repeatlist[list].B2_time_relative){
				skipreturn=true; break;
			}
		}
		return skipreturn;
	};

//cout<<Isrepeat(0,2,3)<<endl;
//return NULL;

	for(int indexi=0; indexi < (int)_treein->GetEntriesFast(); indexi++){
		_treein->GetEntry(indexi);
		
		if(Isrepeat(mr_Bingo_time_gclock,time1,time2))continue;
			
		repeatA.B1_time = mr_Bingo_time_gclock; repeatA.B2_time=time1; repeatA.B2_time_relative=time2;

			for(int indexj = indexi+1; indexj<(int)_treein->GetEntriesFast(); indexj++){
				_treein->GetEntry(indexj);
				if(Isrepeat(mr_Bingo_time_gclock,time1,time2)) continue;
				repeatB.B1_time = mr_Bingo_time_gclock; repeatB.B2_time=time1; repeatB.B2_time_relative=time2;
				if(repeatA.B2_time == repeatB.B2_time){
					repeatBgroup.push_back(repeatB);
				}
			}

		if(repeatBgroup.size()>0){
			repeatlist.push_back(repeatA);
			for(int indexk=0; indexk<(int)repeatBgroup.size(); indexk++){
				repeatlist.push_back(repeatBgroup[indexk]);
			}
		}

		repeatBgroup.clear();
	}

	vector<match_event> acceptlist;
	
	RefineMatchHit(repeatlist,acceptlist);

	RejectionList temlist;

	for(int indexi=0;indexi<(int)repeatlist.size();indexi++){
		temlist.Add(repeatlist[indexi].B1_time,repeatlist[indexi].B2_time,repeatlist[indexi].B2_time_relative);
	}

	for(int indexi=0;indexi<(int)acceptlist.size();indexi++){// have final rejection list
		temlist.Remove(acceptlist[indexi].B1_time, acceptlist[indexi].B2_time, acceptlist[indexi].B2_time_relative);
	}


		string outpath = FilePath + "beta-tof_tree_all_final.root";

		TFile* fout_temp = new TFile(outpath.c_str(),"RECREATE");

		*_fileout = fout_temp;

		TTree* t_return = new TTree("tbeta_tof_refine_all","tbeta_tof_refine_all");

		string IonName;

		t_return->Branch("sweeps_global_mr",&mr_Bingo_sweeps_global,"sweeps_global_mr/L");
		t_return->Branch("time_mr",&mr_Bingo_time,"time_mr/D");
		t_return->Branch("time_gclock_mr",&mr_Bingo_time_gclock,"time_gclock_mr/L");
		t_return->Branch("sweeps_gclock_1",&sweeps_gclock_1,"sweeps_gclock_1/L");
		t_return->Branch("time1",&time1,"time1/L");
		t_return->Branch("time1_relative",&time1_relative,"time1_relative/L");
		t_return->Branch("adc1",&adc1,"adc1/I");
		t_return->Branch("EkeV_1",&EkeV_1,"EkeV_1/D");
		t_return->Branch("sweeps_gclock_2",&sweeps_gclock_2,"sweeps_gclock_2/L");
		t_return->Branch("time2",&time2,"time2/L");
		t_return->Branch("time2_relative",&time2_relative,"time2_relative/L");
		t_return->Branch("adc2",&adc2,"adc2/I");
		t_return->Branch("EkeV_2",&EkeV_2,"EkeV_2/D");
		t_return->Branch("delta_t",&delta_t,"delta_t/I");
		t_return->Branch("delta_t_decay",&delta_t_decay,"delta_t_decay/L");
		t_return->Branch("IonName",&IonName);
		t_return->Branch("IonFillIndex",&IonFillIndex,"IonFillIndex/I");
		t_return->Branch("B_B_gate_in_ns",&B_B_gate_in_ns,"B_B_gate_in_ns/L");
		t_return->Branch("Halflive",&Halflive,"Halflive/L");
		t_return->Branch("TimesOfHalflive",&TimesOfHalflive,"TimesOfHalflive/D");

		for(int nevt=0;nevt<(int)_treein->GetEntriesFast();nevt++){
			_treein->GetEntry(nevt);
			if(temlist.IsInside(mr_Bingo_time_gclock,time1,time2)) continue;
			else{
				IonName = *IonName_ptr;
				t_return->Fill();
			}
		}

		t_return->Write();
		t_return->ResetBranchAddresses();
		return t_return;

}





void ScanBetaTof_tree(TCutG** _usecut=NULL, string speciesname="", RejectionList* _rejectlist=NULL, string GeneralSelection="",bool IsSave=false, bool IsAdd=true){
	if(!gSystem->AccessPathName(FilePath.c_str()) && IsSave){
		speciesname = "start: " + speciesname;
		string command_s ;
		if(IsAdd){command_s = ".!echo " + speciesname + " >> " + FilePath + "beta-tof.log";}
		else{command_s = ".!echo " + speciesname + " > " + FilePath + "beta-tof.log";}
		gROOT->ProcessLine(command_s.c_str());

	//	string command_s = ".>> " + FilePath + "beta-tof.log";
	//	gROOT->ProcessLine(command_s);
		//gROOT->ProcessLine(".>> tree.log");
		string command_coindtion = CutCondition(_usecut);

		if(_rejectlist!=NULL){
			if(_rejectlist->GetN()) command_coindtion += _rejectlist->RejectList_s();
		}

		if(GeneralSelection!=""){
			if(command_coindtion!=""){
				command_coindtion = "(" + command_coindtion + ") && (" +GeneralSelection + ")";
			}
			else command_coindtion = GeneralSelection;
		}

fin->cd();
		//string 
		//gROOT->ProcessLine("tbeta_tof->Scan(\"sweeps_global_mr:time_mr:sweeps_gclock_1:time1:adc1:sweeps_gclock_2:time2:adc2:delta_t:delta_t_decay\",CutCondition(_usecut).c_str(),\"col=lld:12.3f:lld:12lld::lld:12lld:::.3f\")");
		//gROOT->ProcessLine(".>");
		command_s = ".>> " + FilePath +"beta-tof.log";
		gROOT->ProcessLine(command_s.c_str());
		command_s = "tbeta_tof->Scan(\"sweeps_global_mr:time_mr:time_gclock_mr:sweeps_gclock_1:time1:adc1:time2:adc2:delta_t:delta_t_decay\",\"" 
						+command_coindtion + "\",\"col=lld:12.3f:15lld:lld:15lld::15lld:::13lld\");" ;
		tbeta_tof->SetScanField(0);
		gROOT->ProcessLine(command_s.c_str());
		gROOT->ProcessLine(".>");
		command_s = ".!echo end >> " + FilePath + "beta-tof.log";
		gROOT->ProcessLine(command_s.c_str());
		//cout<<command_s.c_str()<<endl;
		cout<<"save beta-tof file to "<<FilePath.c_str()<<"beta-tof.log"<<endl;
	}
	else{	cout<<"\e[1;33m"<<"file without save"<<"\e[0m"<<endl;}

	string condition_s;
	tbeta_tof->SetScanField(0);
	if(_rejectlist!=NULL){
		if(_rejectlist->GetN()){
			condition_s = CutCondition(_usecut) + _rejectlist->RejectList_s();
			//cout<<"condition total: "<<condition_s.c_str()<<endl;
		}
		else{condition_s = CutCondition(_usecut);}
	}
	else{ condition_s = CutCondition(_usecut);}

	if(GeneralSelection!=""){
		if(condition_s!=""){
			condition_s = "(" + condition_s + ") && (" +GeneralSelection + ")";
		}
		else condition_s = GeneralSelection;
	}

	//f_beta_tree->cd();
	//f_beta_tof_tree->cd();
	fin->cd();
	tbeta_tof->Scan("sweeps_global_mr:time_mr:time_gclock_mr:sweeps_gclock_1:time1:adc1:time2:adc2:delta_t:delta_t_decay",
	condition_s.c_str(),"col=lld:12.3f:15lld:lld:15lld::15lld:::13lld");

	cout<<endl;
	cout<<condition_s.c_str()<<endl;

	fin->cd();
	
}


void GenerateRejectList(TTree* inRawTree = NULL, TCutG** _usecut=NULL, RejectionList* _rejectlist=NULL){
	if(inRawTree ==NULL){cout<<"Err! nullptr of input tree!!! break!!"<<endl; return;}
	if(_rejectlist == NULL){cout<<"Err! invalid input of _rejectlist"<<endl; return;}

	string command_coindtion = CutCondition(_usecut);

	f_beta_tof_tree->cd();
	TTree* inRawTree_copy = inRawTree->CopyTree(command_coindtion.c_str());

	//cout<<"netry = "<<inRawTree_copy->GetEntriesFast()<<endl;

	match_event OneEvent;
	vector<match_event> match_group;

	inRawTree_copy->SetBranchAddress("time_gclock_mr",&OneEvent.B1_time);  // just borrow the sturcture of match_event// mr_gclock
	inRawTree_copy->SetBranchAddress("time1",&OneEvent.B2_time);				//  time1
	inRawTree_copy->SetBranchAddress("time2",&OneEvent.B2_time_relative);  //time2

	for(int index=0;index<(int)inRawTree_copy->GetEntriesFast();index++){
		OneEvent.clear();
		inRawTree_copy->GetEntry(index);
		match_group.push_back(OneEvent);
	}
/*
	for(int index=0;index<(int)match_group.size();index++){
		printf("%15lld \t %15lld \t %15lld\n",match_group[index].B1_time,match_group[index].B2_time,match_group[index].B2_time_relative);
	}

	cout<<endl;

	inRawTree_copy->Scan("time_gclock_mr:time1:time2","","col=15lld:15lld:15lld",5);*/

	match_group.shrink_to_fit();

	vector<match_event> match_group_refine;

	//cout<<"get events = "<<match_group.size()<<endl;

	RefineMatchHit(match_group,match_group_refine); 

	cout<<endl;
	
	_rejectlist->Clear();

	for(int index=0;index<(int)match_group.size();index++){
		_rejectlist->Add(match_group[index].B1_time,match_group[index].B2_time,match_group[index].B2_time_relative);
	}

	for(int index=0;index<(int)match_group_refine.size();index++){
		_rejectlist->Remove(match_group_refine[index].B1_time, match_group_refine[index].B2_time, match_group_refine[index].B2_time_relative);
	}

	cout<<"nevts rejected = "<<_rejectlist->GetN()<<endl;
	cout<<endl;

	inRawTree_copy->ResetBranchAddresses();
	delete inRawTree_copy;

	fin->cd();

}


void SaveBeta_tof_tree(TCutG** _usecut=NULL, string speciesname="", RejectionList* _rejectlist=NULL, string GeneralSelection="",bool IsAdd=true){
		if(_usecut==NULL || _usecut[0]==NULL){
			cout<<"No cut is used??? continue??		[y/n]"<<endl;
			char key=' ';
			while(key!='y' && key!='n') cin>>key;
			if(key=='n'){cout<<"Abort!"<<endl; return;}
		}

		if( speciesname==""){
			cout<<"No IonName is specified ??? continue??		[y/n]"<<endl;
			char key=' ';
			while(key!='y' && key!='n') cin>>key;
			if(key=='n'){cout<<"Abort"<<endl; return;}
		}

		string command_coindtion = CutCondition(_usecut);

if(GeneralSelection!=""){
	if(command_coindtion!="") command_coindtion = "(" + command_coindtion + ") && (" + GeneralSelection + ")";
	else command_coindtion = GeneralSelection;
}

		TTree* tbeta_tof_tem = tbeta_tof->CopyTree(command_coindtion.c_str()); // with Cut
		if(tbeta_tof_tem==NULL){
			cout<<"\e[1;31m Fail to save the tree!! Check the selection condition!! Abort!! \e[0m"<<endl;
			return;
		}

		tbeta_tof_tem->SetName("tbeta_tof_tem");

		Long64_t mr_Bingo_sweeps_global=0; double mr_Bingo_time=0; Long64_t mr_Bingo_time_gclock=0;
		Long64_t time1=0; Long64_t time1_relative=0; Long64_t sweeps_gclock_1=0; int adc1=0; double EkeV_1=0;
		Long64_t time2=0; Long64_t time2_relative=0;  Long64_t sweeps_gclock_2=0; int adc2=0; double EkeV_2=0;
		int delta_t=0;
		Long64_t delta_t_decay=0; int tof_index=0;

	tbeta_tof_tem->SetBranchAddress("Ntof_counts",&Ntof_counter);
	tbeta_tof_tem->SetBranchAddress("tof_index",&tof_index);
		tbeta_tof_tem->SetBranchAddress("sweeps_global_mr",&mr_Bingo_sweeps_global);
		tbeta_tof_tem->SetBranchAddress("time_mr",&mr_Bingo_time);
		tbeta_tof_tem->SetBranchAddress("time_gclock_mr",&mr_Bingo_time_gclock);
		tbeta_tof_tem->SetBranchAddress("sweeps_gclock_1",&sweeps_gclock_1);
		tbeta_tof_tem->SetBranchAddress("time1",&time1);
		tbeta_tof_tem->SetBranchAddress("time1_relative",&time1_relative);
		tbeta_tof_tem->SetBranchAddress("adc1",&adc1);
		tbeta_tof_tem->SetBranchAddress("EkeV_1",&EkeV_1);
		tbeta_tof_tem->SetBranchAddress("sweeps_gclock_2",&sweeps_gclock_2);
		tbeta_tof_tem->SetBranchAddress("time2",&time2);
		tbeta_tof_tem->SetBranchAddress("time2_relative",&time2_relative);
		tbeta_tof_tem->SetBranchAddress("adc2",&adc2);
		tbeta_tof_tem->SetBranchAddress("EkeV_2",&EkeV_2);
		tbeta_tof_tem->SetBranchAddress("delta_t",&delta_t);
		tbeta_tof_tem->SetBranchAddress("delta_t_decay",&delta_t_decay);

		TFile* fout=NULL;
		string outpath = FilePath + "beta-tof_tree_raw.root";

		TTree* tbeta_tof_refine;
		string IonName = speciesname;
		string* IonName_ptr = &IonName;
		Long64_t B_B_gate_in_ns = gate_time;
		Long64_t Halflive = RIHalflive/1000000; // in [ms]  originally this is in ns in the program need to convert
		double TimesOfHalflive = Times2Halflive; // window for coincidence with width equal to "TimesOfHalflive X Halflive";
		int IonFillIndex=1;   // record how many time tht same species of ion fill

		if(IsAdd==false || gSystem->AccessPathName(outpath.c_str())){
			 cout<<"Create new file for saving...."<<endl;
			 fout= new TFile(outpath.c_str(),"RECREATE");
			 tbeta_tof_refine = new TTree("tbeta_tof_refine","tbeta_tof_refine");


		tbeta_tof_refine->Branch("Ntof_counts",&Ntof_counter,"Ntof_counts/I");
		tbeta_tof_refine->Branch("tof_index",&tof_index,"tof_index/I");
			tbeta_tof_refine->Branch("sweeps_global_mr",&mr_Bingo_sweeps_global,"sweeps_global_mr/L");
			tbeta_tof_refine->Branch("time_mr",&mr_Bingo_time,"time_mr/D");
			tbeta_tof_refine->Branch("time_gclock_mr",&mr_Bingo_time_gclock,"time_gclock_mr/L");
			tbeta_tof_refine->Branch("sweeps_gclock_1",&sweeps_gclock_1,"sweeps_gclock_1/L");
			tbeta_tof_refine->Branch("time1",&time1,"time1/L");
			tbeta_tof_refine->Branch("time1_relative",&time1_relative,"time1_relative/L");
			tbeta_tof_refine->Branch("adc1",&adc1,"adc1/I");
			tbeta_tof_refine->Branch("EkeV_1",&EkeV_1,"EkeV_1/D");
			tbeta_tof_refine->Branch("sweeps_gclock_2",&sweeps_gclock_2,"sweeps_gclock_2/L");
			tbeta_tof_refine->Branch("time2",&time2,"time2/L");
			tbeta_tof_refine->Branch("time2_relative",&time2_relative,"time2_relative/L");
			tbeta_tof_refine->Branch("adc2",&adc2,"adc2/I");
			tbeta_tof_refine->Branch("EkeV_2",&EkeV_2,"EkeV_2/D");
			tbeta_tof_refine->Branch("delta_t",&delta_t,"delta_t/I");
			tbeta_tof_refine->Branch("delta_t_decay",&delta_t_decay,"delta_t_decay/L");
			tbeta_tof_refine->Branch("IonName",&IonName);
			tbeta_tof_refine->Branch("IonFillIndex",&IonFillIndex,"IonFillIndex/I");
			tbeta_tof_refine->Branch("B_B_gate_in_ns",&B_B_gate_in_ns,"B_B_gate_in_ns/L");
			tbeta_tof_refine->Branch("Halflive",&Halflive,"Halflive/L");
			tbeta_tof_refine->Branch("TimesOfHalflive",&TimesOfHalflive,"TimesOfHalflive/D");

		} 
		else{

			cout<<"Update file"<<endl;
			fout= new TFile(outpath.c_str(),"UPDATE");
			if(!fout->IsOpen()){cout<<"fail to open file!!! Abort!!!"<<endl; return;}

			tbeta_tof_refine= (TTree*)fout->Get("tbeta_tof_refine");

tbeta_tof_refine->SetBranchAddress("Ntof_counts",&Ntof_counter);
tbeta_tof_refine->SetBranchAddress("tof_index",&tof_index);
			tbeta_tof_refine->SetBranchAddress("sweeps_global_mr",&mr_Bingo_sweeps_global);
			tbeta_tof_refine->SetBranchAddress("time_mr",&mr_Bingo_time);
			tbeta_tof_refine->SetBranchAddress("time_gclock_mr",&mr_Bingo_time_gclock);
			tbeta_tof_refine->SetBranchAddress("sweeps_gclock_1",&sweeps_gclock_1);
			tbeta_tof_refine->SetBranchAddress("time1",&time1);
			tbeta_tof_refine->SetBranchAddress("time1_relative",&time1_relative);
			tbeta_tof_refine->SetBranchAddress("adc1",&adc1);
			tbeta_tof_refine->SetBranchAddress("EkeV_1",&EkeV_1);
			tbeta_tof_refine->SetBranchAddress("sweeps_gclock_2",&sweeps_gclock_2);
			tbeta_tof_refine->SetBranchAddress("time2",&time2);
			tbeta_tof_refine->SetBranchAddress("time2_relative",&time2_relative);
			tbeta_tof_refine->SetBranchAddress("adc2",&adc2);
			tbeta_tof_refine->SetBranchAddress("EkeV_2",&EkeV_2);
			tbeta_tof_refine->SetBranchAddress("delta_t",&delta_t);
			tbeta_tof_refine->SetBranchAddress("delta_t_decay",&delta_t_decay);
			tbeta_tof_refine->SetBranchAddress("IonName",&IonName_ptr);
			tbeta_tof_refine->SetBranchAddress("IonFillIndex",&IonFillIndex);
			tbeta_tof_refine->SetBranchAddress("B_B_gate_in_ns",&B_B_gate_in_ns);
			tbeta_tof_refine->SetBranchAddress("Halflive",&Halflive);
			tbeta_tof_refine->SetBranchAddress("TimesOfHalflive",&TimesOfHalflive);

			int MaxFillIndex=0;
			for(int index=0;index<(int)tbeta_tof_refine->GetEntriesFast();index++){
				tbeta_tof_refine->GetEntry(index);
				if(*IonName_ptr == speciesname && IonFillIndex>MaxFillIndex ) MaxFillIndex = IonFillIndex;
			}

			*IonName_ptr= speciesname; // update name
			IonFillIndex = MaxFillIndex+1; // index for next entry
		}


		bool shouldkeep=true;

		for(int index=0;index<(int)tbeta_tof_tem->GetEntriesFast();index++){
			tbeta_tof_tem->GetEntry(index);
			shouldkeep=true;
			if(_rejectlist!=NULL){
				shouldkeep = !(_rejectlist->IsInside(mr_Bingo_time_gclock ,time1,time2));
			}

			if(shouldkeep){
				tbeta_tof_refine->Fill();
			}
		}

	
		tbeta_tof_refine->Write();
		fout->Write();
		
		tbeta_tof_refine->Scan("sweeps_global_mr:time_mr:time_gclock_mr:sweeps_gclock_1:time1:adc1:time2:adc2:delta_t:delta_t_decay","",
			"col=lld:12.3f:15lld:lld:15lld::15lld:::13lld");

		cout<<"Num of Event:   "<<tbeta_tof_refine->GetEntriesFast()<<endl;
		cout<<"tbeta_tof_refine is saved"<<endl;

		fout->Close();
		

}


//%%%%%%%%%%%%%%%%%%%% show how the number of beta-tof counts varies with energy **********************
//for scanning the energy threshold of Si
TCanvas* c_keV_decaytime=NULL;
TGraph* gKev_Decay_Si1=NULL;
TGraph* gKev_Decay_Si2=NULL;
void showEnergyVsDecayTime(string GeneralSelection="", int Draw1Si1_2Si2_3All=1){
		
		long long NEntriesGet =0;
		if(!alphaMode){
			NEntriesGet = tbeta_tof->Draw("delta_t_decay/1e9:EkeV_1:EkeV_2",GeneralSelection.c_str(),"goff");
		}
		else{
			NEntriesGet = tbeta_tof->Draw("delta_t_decay/1e9:EkeV_1",GeneralSelection.c_str(),"goff");
		}


		if(NEntriesGet == 0){
			cout<<"\e[1;31m"<<"Find no responding data, Please check!!"<<"\e[0m"<<endl;
			return;
		} 

		if(Draw1Si1_2Si2_3All<1 || Draw1Si1_2Si2_3All>3){
			cout<<"\e[1;31m Please select the correct draw option: 1 for Si1 only; 2 for Si2 only; 3 for Si1 and Si2. Abort!!! \e[0m"<<endl;
			return;
		}

		tbeta_tof->SetEstimate(NEntriesGet);
		if(c_keV_decaytime!=NULL && c_keV_decaytime->GetCanvasImp()!=NULL) delete c_keV_decaytime;

		c_keV_decaytime = new TCanvas("c_keV_decaytime","Beta Energy Vs. Decay time");	
		
		if(!alphaMode){
			NEntriesGet = tbeta_tof->Draw("delta_t_decay/1e9:EkeV_1:EkeV_2",GeneralSelection.c_str(),"goff");
			if(gKev_Decay_Si1!=NULL) delete gKev_Decay_Si1;
				gKev_Decay_Si1 = new TGraph(tbeta_tof->GetSelectedRows(),tbeta_tof->GetV1(),tbeta_tof->GetV2());
			if(gKev_Decay_Si2!=NULL) delete gKev_Decay_Si2;
				gKev_Decay_Si2= new TGraph(tbeta_tof->GetSelectedRows(),tbeta_tof->GetV1(),tbeta_tof->GetV3());
				gKev_Decay_Si2->SetMarkerColor(kRed);

			if(Draw1Si1_2Si2_3All ==1 || Draw1Si1_2Si2_3All==3){
					gKev_Decay_Si1->GetXaxis()->SetTitle("Decay time [s]");
					gKev_Decay_Si1->GetXaxis()->CenterTitle();
				if(Draw1Si1_2Si2_3All ==1) gKev_Decay_Si1->GetYaxis()->SetTitle("Deposit Energy of Beta1 [keV]");
				else gKev_Decay_Si1->GetYaxis()->SetTitle("Deposit Energy of Beta1 and Beta2 [keV]");
					gKev_Decay_Si1->GetYaxis()->CenterTitle();
			}else{

				gKev_Decay_Si2->GetXaxis()->SetTitle("Decay time [s]");
				gKev_Decay_Si2->GetXaxis()->CenterTitle();
				gKev_Decay_Si2->GetYaxis()->SetTitle("Deposit Energy of Beta2 [keV]");
				gKev_Decay_Si2->GetYaxis()->CenterTitle();
				gKev_Decay_Si2->Draw("AP*");
				return;
			}


				gKev_Decay_Si1->Draw("AP*");
				if(Draw1Si1_2Si2_3All ==3){
				 	gKev_Decay_Si2->Draw("P* same");
				 	double MaximumSi2Y = gKev_Decay_Si2->GetHistogram()->GetMaximum();
				 	double MinimumSi2Y = gKev_Decay_Si2->GetHistogram()->GetMinimum();
				 	if( MaximumSi2Y> gKev_Decay_Si1->GetHistogram()->GetMaximum()) gKev_Decay_Si1->GetHistogram()->SetMaximum(MaximumSi2Y);
				 	if(MinimumSi2Y < gKev_Decay_Si1->GetHistogram()->GetMinimum()) gKev_Decay_Si1->GetHistogram()->SetMinimum(MinimumSi2Y);
				 }
		}
		else{
			NEntriesGet = tbeta_tof->Draw("delta_t_decay/1e9:EkeV_1",GeneralSelection.c_str(),"goff");
			if(gKev_Decay_Si1!=NULL) delete gKev_Decay_Si1;
			gKev_Decay_Si1 = new TGraph(tbeta_tof->GetSelectedRows(),tbeta_tof->GetV1(),tbeta_tof->GetV2());
			gKev_Decay_Si1->GetXaxis()->SetTitle("Decay time [s]");
			gKev_Decay_Si1->GetXaxis()->CenterTitle();
			gKev_Decay_Si1->GetYaxis()->SetTitle("Deposit Energy of Alpha [keV]");
			gKev_Decay_Si1->GetYaxis()->CenterTitle();
			gKev_Decay_Si1->Draw("AP*");
			if(Draw1Si1_2Si2_3All>1){
				cout<<"\e[1;33m Only energy at Si1 is shown under alpha mode \e[0m"<<endl;
			}
		}

		return;
}





//************** Setcuts ******************
TCutG *mycut[5];

bool MakeCut(){  
  //TCutG *cut[5];
  for(int i=0;i<5;i++){
  	if(mycut[i]!=NULL){
  	 	gROOT->GetListOfSpecials()->Remove((TObject*)mycut[i]);
  	 	delete mycut[i];
  	}
  	mycut[i] = NULL;
  }

  TObject* obj;
  TIter iter(gROOT->GetListOfSpecials());
  while((obj=(TObject*)iter())){
  		string cutname = obj->GetName();
  		if(cutname=="no_0" || cutname=="no_1" || cutname=="no_2" || cutname=="no_3" || cutname=="no_4"){gROOT->GetListOfSpecials()->Remove(obj);}
  }


  if(c_beta_beta == NULL){cout<<"canvas is not existed, ShowBeta_Beta_Coin() first!!"<<endl; return false;}
  if(c_eject == NULL){cout<<"canvas is not existed, run  ShowBeta_Beta_Coin() then ShowEjectionBeta(EJE0, EJE1)!!"<<endl; return false;}

  cout << "You can draw a cut by clicking View->Toolbar->Graphical Cut" << endl;
  cout << "Double-click to close cut" << endl<<endl;
  cout<< "\e[1;33m"<<"No.0 Set Energy noise cut at c_beta_beta->cd(2)....."<<"\e[0m"<<endl;
  mycut[0] = (TCutG*)c_beta_beta->cd(2)->WaitPrimitive("CUTG");
  mycut[0]->SetName("no_0"); 
  gROOT->GetListOfSpecials()->Add((TObject *)mycut[0]);
 /* cout << "Cut " << name << " with "<< cut->GetN() << " points "<< endl;
  vector<double> x;
  vector<double> y;
  x.resize(cut->GetN());
  y.resize(cut->GetN());
  for(int n=0;n<cut->GetN();n++){
    cut->GetPoint(n,x[n],y[n]);
    cout << x[n] << "\t" << y[n] << endl;
  }*/

  int cutindex=1;  // in case skip some cut on TGraph
  for(int i=1;i<5;i++){
  	cout<< "\e[1;33m"<<"No."<<i<<" Set Time noise cut at c_eject->cd("<<i<<")? ('y' or 'n')"<<"\e[0m"<<endl;
  	char yesorno='\0';
  	while(1){
  		cin>>yesorno;
  		if(yesorno=='n') break;
  		if(yesorno=='y'){
			  cout<< "\e[1;33m"<<"Draw Time noise cut at c_eject->cd("<<i<<")........"<<"\e[0m"<<endl;
			  char cutname[5]={'\0'};
			  sprintf(cutname,"no_%d",i);
			  mycut[cutindex] = (TCutG*)c_eject->cd(i)->WaitPrimitive("CUTG");
			  mycut[cutindex]->SetName(cutname); 
			  gROOT->GetListOfSpecials()->Add((TObject *)mycut[cutindex++]);
			  break;
  		}
  	}  	

  }


  TFile *cutf;
  
  	string cutf_name = FilePath + "cutfile.root";
    cutf = new TFile(cutf_name.c_str(),"RECREATE");
    if(!cutf->IsOpen()){ cout<<"Can not recreate file.... Please check the path.... Abort!!!!"<<endl; return false;}
	  cutf->cd();
	  for(int i=0;i<5;i++){
	  		if(mycut[i]!=NULL)mycut[i]->Write(mycut[i]->GetName());
	  }
	  cutf->ls();
	  cutf->Close();
	  fin->cd();

	  return true;
}


bool LoadCut(){
fin->cd();
 for(int i=0;i<5;i++){
  	if(mycut[i]!=NULL){
  	 	gROOT->GetListOfSpecials()->Remove((TObject*)mycut[i]);
  	 	delete mycut[i];
  	}
  	mycut[i] = NULL;
  }

  TObject* obj;
  TIter iter(gROOT->GetListOfSpecials());
  while((obj=(TObject*)iter())){
  		string cutname = obj->GetName();
  		if(cutname=="no_0" || cutname=="no_1" || cutname=="no_2" || cutname=="no_3" || cutname=="no_4"){gROOT->GetListOfSpecials()->Remove(obj);}
  }

  int ncut=0;
  string cutf_name = FilePath + "cutfile.root";
  TFile *fcut_in = new TFile(cutf_name.c_str(),"READ");
  if(!fcut_in->IsOpen()){ cout<<"Can not cut file is not existed!!!! Please check .... Abort!!!!"<<endl; return false;}
 // TObject *obj;
  TIter nextobj(fcut_in->GetListOfKeys());
  while( (obj = (TObject *)nextobj()) ){
    if(fcut_in->Get(obj->GetName())->InheritsFrom("TCutG")){
      mycut[ncut] = (TCutG *)fcut_in->Get(obj->GetName());
      ncut++;
    }
  }
/*
  for(int cutindex=0;cutindex<5;cutindex++){
  	if(mycut[cutindex]!=NULL){
  		gROOT->GetListOfSpecials()->Add( (TObject*)mycut[cutindex] );
  	}
  }*/

  fcut_in->ls();
  fcut_in->Close();
  fin->cd();

  return true;
}

//*******************************************


vector<int> countrate_i_cp[40];
vector<double> tof_ion[40][5]; //[roi_index][1 st count,2nd count,3rd,4rd]; store tof info each hit in ROI

vector<Long64_t> sweep_glo_Beta_Beta;

//************ rate meter ************************************
void RefreshRateMeter(int rate_in_second=1, bool showExtra=false, bool ShowB2Bcoin=false){ // set showB2Bcoin to show beta-beta coin rate
		static TMultiGraph* mtg=NULL;
		static TLegend* lg=NULL;

		if(mtg!=NULL){delete lg; mtg->Clear(); delete mtg;}
		mtg = new TMultiGraph("mtg","mtg");
		mtg->SetTitle(Form("Rate meter;Time [s];Counts per %ds",rate_in_second));
	mtg->SetMinimum(0.01);
		lg = new TLegend(0.75,0.7,0.9,0.9);



		static TCanvas* C_ratermeter=NULL;
		if(C_ratermeter!=NULL && C_ratermeter->GetCanvasImp()!=NULL) delete C_ratermeter; // window is closed
		C_ratermeter = new TCanvas("ratermeter","RaterMeter",900,700); 

		Int_t kala[]={633,808,799,417,433,600,617,1};
		Int_t MarkerStyle[]={20,25,22,27,28};

		vector<Long64_t> *ptr_extra_vector[2];
		ptr_extra_vector[0]=NULL;
		ptr_extra_vector[1]=NULL;

		vector<string> extra_name;
		int extra_index=0;
		if(showExtra){
				if(sweep_glo_coin.size()>0){extra_index++; extra_name.push_back("w/_coin");  ptr_extra_vector[0]= &sweep_glo_coin; }
				if(eve_woutcoin_sweep.size()>0){extra_index++; extra_name.push_back("w/o_coin"); ptr_extra_vector[1]= &eve_woutcoin_sweep;}
		}

		const int sweeps_per_second=40;
		long update_interval = rate_in_second * sweeps_per_second;  // convert update frequency in the unit of second to sweep number;
		long sweeps_max=0;
		long counts_max=0;

		vector <Double_t> *ion_time = new vector <Double_t>();
		vector <Long_t> * ion_sweeps_global = new vector <Long_t>;
		vector <Int_t> *ion_tag = new vector <Int_t>();
		vector <Int_t> *ion_channel = new vector <Int_t>();

		
		intree->SetBranchAddress("sweeps_global",&ion_sweeps_global);
		intree->SetBranchAddress("tag",&ion_tag);
		intree->SetBranchAddress("channel",&ion_channel);
		string treename = intree->GetName();
		if(treename == "tree") intree->SetBranchAddress("time",&ion_time);
		else{ intree->SetBranchAddress("timec",&ion_time);}

		vector<Long_t> sweeps_global_ROI_i[40];
		vector<int> countrate_i[40];

for(int roi_index=0;roi_index<40;roi_index++){
	tof_ion[roi_index][0].clear();
	tof_ion[roi_index][1].clear();
	tof_ion[roi_index][2].clear();
	tof_ion[roi_index][3].clear();
	tof_ion[roi_index][4].clear();

}
		
		double* countrate_i_array[40];
		double* countrate_i_err_array[40];
		double* time_in_second[40];
		double* time_in_second_err[40];

		static TGraphErrors* ratemeter[40]={NULL};
/*		static TMultiGraph* mtg=NULL;
		static TLegend* lg=NULL;

		if(mtg!=NULL){delete lg; delete mtg;}
		mtg = new TMultiGraph("mtg","mtg");
		mtg->SetTitle(Form("Rater meter;Time [s];Counts per %ds",rate_in_second));
		mtg->GetXaxis()->CenterTitle();
		mtg->GetYaxis()->CenterTitle();
		lg = new TLegend(0.75,0.7,0.9,0.9);*/

		for(long nevt=0;nevt<intree->GetEntriesFast();nevt++){ // generate vector of sweeps corresponding ROI
				intree->GetEntry(nevt);
int multi_hit_index[40]={0};
				double tof_last_hit=0;  // for reject multi-hit events
				for(unsigned long nevt_sub=0;nevt_sub<ion_time->size();nevt_sub++){
					if(ion_channel->at(nevt_sub) != Tof_Channel2Read) continue; // reject hit event from other channel
					if(ion_time->at(nevt_sub)-tof_last_hit < Tof_veto) continue; // reject multiple hits

					tof_last_hit = ion_time->at(nevt_sub); // update tof of last hit

					for(int i=0;i<ROI_INDEX && i<40;i++){
						if(ion_time->at(nevt_sub) >= ROIL[i].GetX1() && ion_time->at(nevt_sub)<=ROIR[i].GetX1() && ion_tag->at(nevt_sub)==ROI_tag[i]){// match with ROI[i]
if(multi_hit_index[i]==0)		sweeps_global_ROI_i[i].push_back(ion_sweeps_global->at(nevt_sub)); //delete  if(multi_hit_index==0) to accept all hit
// in the same sweep
if(multi_hit_index[i]<5)tof_ion[i][multi_hit_index[i]].push_back(ion_time->at(nevt_sub));
multi_hit_index[i]++;

							break;
						}
					}
				}
				if(ion_sweeps_global->at(0)>sweeps_max) sweeps_max = ion_sweeps_global->at(0);  // find out the maximum of sweeps of the whole file
				ion_sweeps_global->clear();	ion_tag->clear();	ion_time->clear();
		}

int extragraph_num=0;

		if(showExtra){
			int j= ROI_INDEX;
			for(int i=0; i<2 && j<40;i++){
				if(ptr_extra_vector[i] == NULL) continue;
				for(unsigned long hit_index=0; hit_index<ptr_extra_vector[i]->size();hit_index++){
					sweeps_global_ROI_i[j].push_back(ptr_extra_vector[i]->at(hit_index));
				}
				j++;
				extragraph_num++;
			}
		}

int betabetaGraph_num=0;

if(ShowB2Bcoin){
		int j= ROI_INDEX;
		if(j ==0){
			for(unsigned long betabetaIndex=0;betabetaIndex<sweep_glo_Beta_Beta.size();betabetaIndex++){
				sweeps_global_ROI_i[j].push_back(sweep_glo_Beta_Beta[betabetaIndex]);
			}
			betabetaGraph_num=1;
		}
		else if(j+extra_index<40){
			for(unsigned long betabetaIndex=0;betabetaIndex<sweep_glo_Beta_Beta.size();betabetaIndex++){
				sweeps_global_ROI_i[j+extra_index].push_back(sweep_glo_Beta_Beta[betabetaIndex]);
			}
			betabetaGraph_num=1;
		}

}



		for(int i=0;i<ROI_INDEX+extra_index+betabetaGraph_num && i<40;i++){
			sort(sweeps_global_ROI_i[i].begin(),sweeps_global_ROI_i[i].end());  // make sure global sweeps keep increase
			int counter=0;
			long marker_num=1;
//cout<<"ROI_"<<i<<" :"<<sweeps_global_ROI_i[i].size()<<endl;  // test purpose
			for(unsigned long nevt=0;nevt<sweeps_global_ROI_i[i].size();nevt++){
				long num_should_fill = sweeps_global_ROI_i[i][nevt]/update_interval;
				long num_record_lack = num_should_fill-(long)countrate_i[i].size(); // set 0 count at some interval
				for(long Nloop=0; Nloop<num_record_lack;Nloop++){
					countrate_i[i].push_back(counter);
countrate_i_cp[i].push_back(counter);					
					if(counts_max<counter)counts_max=counter; // for displaying Y range
					counter=0;
				}
				counter++;
			}
			countrate_i[i].push_back(counter);
countrate_i_cp[i].push_back(counter);			


			long num_record_lack = sweeps_max/update_interval-(long)countrate_i[i].size();
			if(sweeps_max%update_interval>0)num_record_lack++; // last point
			for(long fill0_i=0; fill0_i < num_record_lack;fill0_i++){countrate_i[i].push_back(0);
countrate_i_cp[i].push_back(0);
			}


			countrate_i_array[i] = new double[countrate_i[i].size()];
			countrate_i_err_array[i] = new double[countrate_i[i].size()];
			time_in_second[i] = new double[countrate_i[i].size()];
			time_in_second_err[i] = new double[countrate_i[i].size()];

			//if(ratemeter[i]!=NULL) delete ratemeter[i]; // Important!!!!! When TMultiGraph deleted, graph added to it will be deleted as well


			for(unsigned long index=0;index<countrate_i[i].size(); index++){
				if(countrate_i[i][index]>0) countrate_i_array[i][index]=countrate_i[i][index];
				else countrate_i_array[i][index] = 0.01;

				countrate_i_err_array[i][index] = TMath::Sqrt(countrate_i[i][index]);
				time_in_second[i][index]=(index+1)*rate_in_second;
				time_in_second_err[i][index]=0;
			}

			ratemeter[i] = new TGraphErrors(countrate_i[i].size(),time_in_second[i],countrate_i_array[i],time_in_second_err[i],countrate_i_err_array[i]);
			
			if(i<ROI_INDEX){

				ratemeter[i]->SetName(Form("ROI%d",i+1));
				ratemeter[i]->SetTitle( Form("ROI%d:%d",i+1,(int)tof_ion[i][0].size()) );
				ratemeter[i]->SetLineColor(kala[i%8]);
				ratemeter[i]->SetMarkerStyle(MarkerStyle[i%5]);
				ratemeter[i]->SetMarkerColor(kala[i%8]);
				for(int j=0;j<= TOFMarker::marker_amount && j<40;j++){
					if(marker_tof[j].GetX()>= ROIL[i].GetX1() && marker_tof[j].GetX()<= ROIR[i].GetX1()){
						ratemeter[i]->SetName(Form("%s:%d",marker_tof[j].name,(int)tof_ion[i][0].size()));
						ratemeter[i]->SetTitle(Form("%s:%d",marker_tof[j].name,(int)tof_ion[i][0].size()));//only this set the name of roi in legend
					}
				}
			}
			else{
				if(showExtra && extra_index>0 && i<ROI_INDEX+extra_index){				

						ratemeter[i]->SetName(extra_name[i-ROI_INDEX].c_str());
						ratemeter[i]->SetTitle(extra_name[i-ROI_INDEX].c_str());
						ratemeter[i]->SetLineColor(kala[i%8]);
						ratemeter[i]->SetMarkerStyle(MarkerStyle[i%5]);
						ratemeter[i]->SetMarkerColor(kala[i%8]);					
				}
				else if(ShowB2Bcoin && i>=ROI_INDEX+extra_index){
						ratemeter[i]->SetTitle(Form("B1B2coin:%lu",sweep_glo_Beta_Beta.size()));
						ratemeter[i]->SetLineColor(kala[i%8]);
						ratemeter[i]->SetMarkerStyle(MarkerStyle[i%5]);
						ratemeter[i]->SetMarkerColor(kala[i%8]);		
				}
			}



			//C_ratermeter->cd();
			/*
			if(i>0)ratemeter[i]->Draw("lp same");
			else{ 
				ratemeter[i]->GetXaxis()->SetTitle("Time [s]");
				ratemeter[i]->GetXaxis()->CenterTitle();
				ratemeter[i]->GetYaxis()->SetTitle(Form("Counts per %ds",rate_in_second));
				ratemeter[i]->GetYaxis()->CenterTitle();
				ratemeter[i]->Draw("APL");
			}*/

			mtg->Add(ratemeter[i],"lep");
			lg->AddEntry(ratemeter[i],"","lep");

			delete[] countrate_i_array[i];
			delete[] countrate_i_err_array[i];
			delete[] time_in_second[i];
			delete[] time_in_second_err[i];


		}

		C_ratermeter->cd();
		mtg->GetHistogram()->SetMaximum(counts_max*1.12);
		mtg->GetHistogram()->SetMinimum(0.001);
		mtg->GetXaxis()->SetLimits(0,((sweeps_max/sweeps_per_second)/rate_in_second+1)*rate_in_second);
		mtg->GetXaxis()->CenterTitle();
		mtg->GetYaxis()->CenterTitle();
		mtg->Draw("A");
		lg->Draw();
		//C_ratermeter->BuildLegend();

		intree->ResetBranchAddresses();

		c1->cd(2);

}
//********************************************************************************


//////&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



///////&&&&&&&&&&&&&&&&&&&&&&&& Interaction Mode &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#ifndef _EXEC_H_
#include "Exec.h"
#endif
void StartInteractionMode(int canvasID =0){

	if(canvasID>2 || canvasID<0){cout<<"0=>activiate both 1 and 2 pad;   1=> pad1;   2=> pad2"<<endl; return;}
	if(canvasID==0){
		c1->cd(1)->AddExec("exec1","Exec(1)");
		c1->cd(2)->AddExec("exec2","Exec(2)");
	}
	else if(canvasID == 1){
		c1->cd(1)->AddExec("exec1","Exec(1)");
	}
	else c1->cd(2)->AddExec("exec2","Exec(2)");

}

void StopInteractionMode(int canvasID =0){

	if(canvasID>2 || canvasID<0){cout<<"0=>activiate both 1 and 2 pad;   1=> pad1;   2=> pad2"<<endl; return;}
	if(canvasID==0){
		c1->cd(1)->DeleteExec("exec1");
		c1->cd(2)->DeleteExec("exec2");
	}
	else if(canvasID == 1){
		c1->cd(1)->DeleteExec("exec1");
	}
	else c1->cd(2)->DeleteExec("exec2");

}


#endif // end of #ifdef _PREVIEWER_