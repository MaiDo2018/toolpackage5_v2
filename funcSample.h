// Sample function class
#ifndef _FUNCSAMPLE_H_
#define _FUNCSAMPLE_H_

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TH1.h"
#include "TFitResultPtr.h"
#include "TSpline.h"
#include <iostream>
#include <stdlib.h>
#include <vector>



using namespace std;
//TCanvas ctest;  // for showing the MC sim. fitting result
class funcS{
	private:
		vector<double> X_tof;
		vector<double> Y_count;
		double sAmp;   // Amp of sampled peak
		double sPeakCenter; // peak center of sampled peak;
		double sPeakCenter_err;
		double sFWHM;
		double sSigma;  // using 0.6826 define +/- sigma region
		double Amp[10];
		double tof_center[10];
		double tof_center_err[10];
		double sample_range_L;
		double sample_range_R;
		int MaxNPeaks;
		int OldNPeaks;
		TSpline3* spl;
		TF1* fsample;
		TF1* fsample_1p;
		TF1* fresult[10];
		int MC_sim_counts; // number of spectra to be simulated
		double* MC_err_ptr;
//char ctemp;		
	public:
		static int smoothlevel;
		static int NumOfPeaks;
		bool FreeRange;  // FreeRange=false, fix left and right weight ratio
		double range_L;  // left width for fit
		double range_R; // right width for fit
		double bins_width;
		int Nbins;
		TH1D* h_sample;
		bool useMC;


		double* GetMC_err_x(){return MC_err_ptr;}

		void SetMC_sim_counts(int Ncounts){
			MC_sim_counts = Ncounts;
		}

		int GetMC_sim_counts(bool IsPrint=true){ if(IsPrint){cout<<"MC simulation Nloops = "<<MC_sim_counts<<endl;} return MC_sim_counts;	}


		void SetPars(int Peakindex, double _Amp, double _tof_center){// Peakindex >=1
			if(Peakindex>10 || Peakindex<1){
				cout<<"Error!!! Peakindex should be [1,10], abort!!!!!"<<endl;
				return;
			}
			else{
				Amp[Peakindex-1]=_Amp;
				tof_center[Peakindex-1] = _tof_center;
			}
		}

		double GetAmp(int Peakindex=1){
			if(Peakindex>10 || Peakindex<1){
				cout<<"Error!!! Peakindex should be [1,10], abort!!!!!"<<endl;
				return -1;
			}
			else{ return Amp[Peakindex-1];}// Peakindex >=1
		}

		double GetTofCenter(int Peakindex=1){
			if(Peakindex>10 || Peakindex<1){
				cout<<"Error!!! Peakindex should be [1,10], abort!!!!!"<<endl;
				return -1;
			}
			else{ return tof_center[Peakindex-1];} //Peakindex >=1
		}

		double GetTofCenterErr(int Peakindex=1){
			if(Peakindex>10 || Peakindex<1){
				cout<<"Error!!! Peakindex should be [1,10], abort!!!!!"<<endl;
				return -1;
			}
			else{ return tof_center_err[Peakindex-1];} //Peakindex >=1

		}

		double GetsPeakCenter(){return sPeakCenter;}

		double GetsPeakCenter_err(){return sPeakCenter_err;}

		int GetOldNPeaks(){return OldNPeaks;}

		double GetFWHM(){return sFWHM;}


		void UpdateSigma(){
			if(fsample_1p != NULL){  //fsample_1p is valid after Draw() the reference function
				double Area_x = fsample_1p->Integral(sPeakCenter-20*sFWHM, sPeakCenter+20*sFWHM);
				double LowLimitX = sPeakCenter-2*sFWHM;
				double UpLimitX = sPeakCenter;
				double MiddleX = 0;
				double temSigma = 1;
				while(1){
				  MiddleX = (LowLimitX+UpLimitX)*0.5;
				  temSigma = sPeakCenter - MiddleX;
				  double Area_ratio_now = fsample_1p->Integral(MiddleX,sPeakCenter+temSigma)/Area_x;
				  double Area_diff = Area_ratio_now - 0.6826;
				  if(  TMath::Abs(Area_diff) < 0.00001) break;
				  else{
				  	if(Area_diff >0 ) LowLimitX = MiddleX;
				  	else UpLimitX =MiddleX;
				  }
				}

				sSigma = temSigma;

				return;

			}else{
				cout<<"\e[1;32m Fitting function of reference is not available, return sigma = 10\e[0m"<<endl;
				sSigma = 10;
				return ;
			}
		}


		double GetSigma(){
			if(sSigma==10)cout<<"\e[1;32m Fitting function of reference is not available, return sigma = 10\e[0m"<<endl;
			return sSigma;
		}



		double MC_err_ref(TH1D* _hin_sample = NULL, double _half_sigma = 10, int NLoops=10){
			double* Px_fit_result = new double[NLoops];
			int simNcounts = (int)_hin_sample->Integral(1,_hin_sample->GetNbinsX());
			TH1D* _hin_sample_cp = dynamic_cast<TH1D*>(_hin_sample->Clone("_hin_sample_cp"));
			TF1* gfit = new TF1("gfit","gaus",0,25e6);
			gfit->	SetParError(0,0.1);
			gfit->SetParError(1,0.01);
			gfit->SetParError(2,0.1);
		//	gfit->SetParameter(0,sAmp);
			gfit->SetParameter(0,_hin_sample_cp->GetMaximum());
			gfit->SetParLimits(0,0,_hin_sample_cp->GetMaximum()*3);
		//	gfit->SetParameter(1,sPeakCenter);
			double temcenter = _hin_sample_cp->GetBinCenter(_hin_sample_cp->GetMaximumBin());
			gfit->SetParameter(1,temcenter);
			gfit->SetParameter(2,sFWHM/2);

			int ProcessRatio=0;
			for(int index=0;index<NLoops;index++){

				if((int)index*100.0/NLoops>=ProcessRatio){
					cout<<"\e[1;37m MC_sim progress = "<<ProcessRatio<<" %\e[0m"<<"\r"<<flush;
					ProcessRatio+=5;
				}else if(index== NLoops-1){
					cout<<"\e[1;37m MC_sim progress = 100 %\e[0m"<<"\r"<<endl;
					cout<<endl;
				}

				_hin_sample_cp->Reset();
				_hin_sample_cp->FillRandom(_hin_sample,simNcounts);
				_hin_sample_cp->Smooth(smoothlevel);
/*
::ctest.cd();
_hin_sample_cp->Draw();
*/
				for(int i=0;i<20;i++){
						if(i<8){
						//	gfit->SetParLimits(0,0.,10*sAmp);
														
							gfit->SetParLimits(1,temcenter-sFWHM,temcenter+sFWHM);
						if(i%2==0){ _hin_sample_cp->Fit(gfit,"MLBQN","QN",temcenter-sFWHM*0.5,temcenter+sFWHM*0.5);}
						else{_hin_sample_cp->Fit(gfit,"MBQN","QN",temcenter-sFWHM*0.5,temcenter+sFWHM*0.5);}
						}
						else{
							if(_half_sigma<bins_width*15){ //half_sigma<bins_width*15  //low statistics
								double temPy = gfit->GetParameter(0);
								double temPx = gfit->GetParameter(1);
								gfit->SetParLimits(0,0,temPy*3);
								gfit->SetParLimits(1,temPx-_half_sigma*1.5,temPx+_half_sigma*1.5);
								gfit->SetParLimits(2,0,sFWHM*2);
								_hin_sample_cp->Fit(gfit,"MLBQN","NQ",temPx-_half_sigma*1.5,temPx+_half_sigma*1.5);  // option B --> use user-defined parameter limits
	//	if(i==19)						_hin_sample_cp->Fit(gfit,"MLB","",temPx-_half_sigma*1.5,temPx+_half_sigma*1.5);  

							}
							else{
								double temPy = gfit->GetParameter(0);
								double temPx = gfit->GetParameter(1);
								gfit->SetParLimits(0,0,temPy*10 );
								gfit->SetParLimits(1,temPx-bins_width*15,temPx+bins_width*15);
								gfit->SetParLimits(2,0,sFWHM*2);
								_hin_sample_cp->Fit(gfit,"MLBQN","NQ",temPx-bins_width*15,temPx+bins_width*15);
	//if(i==19)							_hin_sample_cp->Fit(gfit,"MLB","",temPx-bins_width*15,temPx+bins_width*15); 
							}
						}
				}

				Px_fit_result[index] = gfit->GetParameter(1);
/*
::ctest.cd();
_hin_sample_cp->Draw();
::ctest.Modified();
::ctest.Update();
cin>>ctemp;
*/
			}

			

			double Err_result = TMath::RMS(NLoops,Px_fit_result);

			delete[] Px_fit_result;
			delete gfit;
			delete _hin_sample_cp;

			return Err_result;

		}


		void Sampling(TH1D* h_in, double _range_L, double _range_R){
			int bin_start = h_in->GetXaxis()->FindBin(_range_L);
			int bin_end = h_in->GetXaxis()->FindBin(_range_R);
			X_tof.clear(); X_tof.shrink_to_fit();
			Y_count.clear(); Y_count.shrink_to_fit();

			for(int i=bin_start;i<=bin_end;i++){
				X_tof.push_back(h_in->GetBinCenter(i));
				Y_count.push_back(h_in->GetBinContent(i));
			}

			Nbins = (int)X_tof.size();
			bins_width = h_in->GetBinWidth(1);
			sample_range_L = h_in->GetBinLowEdge(bin_start);
			sample_range_R = h_in->GetBinLowEdge(bin_end)+bins_width;

		}

		void RecreateHisto(){
			if(h_sample != NULL){delete h_sample;}
			h_sample = new TH1D("h_sample","h_sample",Nbins,sample_range_L,sample_range_R);
			for(int index=0;index<Nbins;index++){
				h_sample->SetBinContent(index+1,Y_count[index]);
			}
			
			FreeRange=true;
		}

		void MakeSplinefunc(){
			const int Npoints = Nbins;
			double X[Npoints], Y[Npoints];
			for(int i=0;i<Npoints;i++){
				X[i] = X_tof[i];
				Y[i] = h_sample->GetBinContent(i+1);
			}

			if(spl !=NULL) delete spl;
			spl = new TSpline3("spl",X,Y,Npoints);

		}


		bool SmoothHisto(){
			if(smoothlevel<0){cout<<"\e[1;33m"<<"error!! smoothlevel must be >=0"<<"\e[0m"<<endl; return false;}

			double FWHM_L=0, FWHM_R =0;
			h_sample->Smooth(smoothlevel);

/*
::ctest.cd();
h_sample->Draw();
::ctest.Modified();
::ctest.Update();
cin>>ctemp;
*/
			
			int Bin_max = h_sample->GetMaximumBin();
			sAmp = h_sample->GetBinContent(Bin_max);
			sPeakCenter = h_sample->GetBinCenter(Bin_max);
			double height_half = 0.5*sAmp;

			double bin_i_y1;
			double bin_i_y2;
			double candidate_x;

			for(int i=1;i<Nbins;i++){
				bin_i_y1 = h_sample->GetBinContent(i);
				bin_i_y2 = h_sample->GetBinContent(i+1);
				candidate_x = (height_half - bin_i_y1) / (bin_i_y2 - bin_i_y1) * bins_width + h_sample->GetBinCenter(i);
				if(i<Bin_max){
					if(bin_i_y1<height_half && bin_i_y2 >height_half){ 	FWHM_L = candidate_x;	}
				}
				else{
					if(bin_i_y1>height_half && bin_i_y2 <height_half){		FWHM_R = candidate_x; break;	}
				}
			}

			sFWHM = FWHM_R - FWHM_L;

			TF1* gfit = new TF1("gfit","gaus",0,25e6);
			gfit->SetParameter(0,sAmp);
			gfit->SetParameter(1,sPeakCenter);
			h_sample->GetXaxis()->SetRangeUser(sPeakCenter-sFWHM,sPeakCenter+sFWHM);
			double half_sigma = h_sample->GetStdDev();
			gfit->SetParameter(2,half_sigma);
			h_sample->GetXaxis()->UnZoom();
			//gfit->SetParameter(2,sFWHM/2.36);
			half_sigma*=0.8; // =>80% of sigma
			double couts_for_fit=0;

			for(int i=0;i<20;i++){
				if(i<8){
		//			gfit->SetParLimits(0,TMath::Max(0., sAmp- 3*TMath::Sqrt(sAmp)), sAmp+3*TMath::Sqrt(sAmp) );
					gfit->SetParLimits(1,sPeakCenter-sFWHM,sPeakCenter+sFWHM);
					if(i%2==0){h_sample->Fit(gfit,"MELQN","NQ",sPeakCenter-sFWHM*0.5,sPeakCenter+sFWHM*0.5);}
					else{h_sample->Fit(gfit,"MEQN","NQ",sPeakCenter-sFWHM*0.5,sPeakCenter+sFWHM*0.5);}

				}
				else{
					if(half_sigma<bins_width*15){ //half_sigma<bins_width*15
						double temPy = gfit->GetParameter(0);
						double temPx = gfit->GetParameter(1);
	//					gfit->SetParLimits(0,TMath::Max(0., temPy- 3*TMath::Sqrt(temPy)), temPy+3*TMath::Sqrt(temPy) );
						gfit->SetParLimits(0,0,temPy*3 );
						gfit->SetParLimits(1,temPx-half_sigma*1.5,temPx+half_sigma*1.5);
						gfit->SetParLimits(2,0,sFWHM*1.5);
						h_sample->Fit(gfit,"MELBQN","NQ",temPx-half_sigma*1.5,temPx+half_sigma*1.5);  // option B --> use user-defined parameter limits
	//	if(i==19) h_sample->Fit(gfit,"MELB","",temPx-half_sigma*1.5,temPx+half_sigma*1.5);
					//h_sample->Fit(gfit,"MELQN","NQ",temPx-half_sigma*0.7,temPx+half_sigma*0.7);
						//couts_for_fit = h_sample->Integral(h_sample->GetXaxis()->FindBin(sPeakCenter-half_sigma),h_sample->GetXaxis()->FindBin(sPeakCenter+half_sigma));
						couts_for_fit = h_sample->Integral(h_sample->GetXaxis()->FindBin(temPx-half_sigma*1.5),h_sample->GetXaxis()->FindBin(temPx+half_sigma*1.5));
					}
					else{
						double temPy = gfit->GetParameter(0);
						double temPx = gfit->GetParameter(1);
					//	gfit->SetParLimits(0,TMath::Max(0., temPy- 3*TMath::Sqrt(temPy)), temPy+3*TMath::Sqrt(temPy) );
					//	gfit->SetParLimits(1,temPx-bins_width*15,temPx+bins_width*15);						
						gfit->SetParLimits(0,0,temPy*3 );
						gfit->SetParLimits(1,temPx-bins_width*15,temPx+bins_width*15);
						gfit->SetParLimits(2,0,sFWHM);
						h_sample->Fit(gfit,"MELBQN","NQ",temPx-bins_width*15,temPx+bins_width*15); 
//if(i==19)h_sample->Fit(gfit,"MEL","",temPx-bins_width*15,temPx+bins_width*15);
						//couts_for_fit = h_sample->Integral(h_sample->GetXaxis()->FindBin(sPeakCenter-bins_width*15),h_sample->GetXaxis()->FindBin(sPeakCenter+bins_width*15));
						couts_for_fit = h_sample->Integral(h_sample->GetXaxis()->FindBin(temPx-bins_width*15),h_sample->GetXaxis()->FindBin(temPx+bins_width*15));
					}
				}

			}


			double error_estimate = sFWHM/(2.36*TMath::Sqrt(couts_for_fit));

			sAmp = gfit->GetParameter(0);
			sPeakCenter = gfit->GetParameter(1);
			sPeakCenter_err =gfit->GetParError(1);
			cout<<endl;
	//		cout<<"\e[1;37m MCtoy error estimated.......\e[0m"<<"\r"<<flush;
			double MC_err=0;
			if(useMC) MC_err = MC_err_ref(h_sample,half_sigma,MC_sim_counts);   // fit 500 random histograms

			h_sample->Scale(1/sAmp);
			cout<<"\e[1;33m"<<"Sampling result:"<<"\e[0m"<<endl;
			printf("Amp = %.2f \t Tof= %.4f(%.4f), err_estimated by counts = %.4f\n",sAmp,sPeakCenter,sPeakCenter_err,error_estimate);
			
			if(useMC)	{
				sPeakCenter_err = MC_err;
				printf("MC estimated err =\033[1;37m %.4f \033[;37m <--Adopt this error with MC_sim activated, N_loop = %d \n",MC_err,MC_sim_counts);
				printf("result: Tof= %.4f(%.4f)\033[0m\n",sPeakCenter,sPeakCenter_err);
			}			
			
			printf("Rm = %f\n",sPeakCenter/(2*sFWHM));
			printf("FWHM = %f\n",sFWHM);
			delete gfit;
			cout<<endl;
			cout<<"smoothlevel= "<<smoothlevel<<endl;
			MakeSplinefunc();
			return true;
		}

		double  fitfunc(double* x, double* par){ // multiple peak depends on global varibale: NumOfPeak
			double value =0;
			for(int i=0;i<OldNPeaks && i<MaxNPeaks;i++){
				double X_convert = x[0]-par[1+i*2]+sPeakCenter;

				//value +=	par[i*2]*h_sample->GetBinContent(h_sample->GetXaxis()->FindBin(x[0]-par[1+i*2]+sPeakCenter));
				if(X_tof[0]<=X_convert && X_convert<=X_tof[Nbins-1]){ value +=  par[i*2]* TMath::Max(0.,spl->Eval(X_convert)); }
				else{value+=0.;}
				//else{printf("x = %.4f, par[%d] = %.4f, par[%d]= %.4f , sPeakCenter=%.4f\n",x[0],1+i*2,par[1+i*2],i*2,par[i*2],sPeakCenter);}
			}
			return value;
			
		}

		double fitfunc_cpy(double* x, double* par){ // just for draw at reference histo
			 //return par[0]*h_sample->GetBinContent(h_sample->GetXaxis()->FindBin(x[0]-par[1]+sPeakCenter));
			double X_convert = x[0]-par[1]+sPeakCenter;
			if(X_convert<X_tof[0] || X_convert>X_tof[Nbins-1]){return 0;}
			else{return par[0]*spl->Eval(X_convert);}
		}

		TF1* Getfitfunc(){return fsample;}

		void Makefitfunc(){
			if(fsample != NULL){delete fsample; fsample=NULL;}
			if(NumOfPeaks>MaxNPeaks){
				cout<<endl;
				cout<<"\e[1;31m NumOfPeaks to be set exceed the limit of MaxNPeaks=10;"<<endl;
				cout<<"Set NumOfPeaks = MaxNPeaks by default\e[0m"<<endl;
				NumOfPeaks = MaxNPeaks;
			}else if(NumOfPeaks<1){
				cout<<endl;
				cout<<"\e[1;31m NumOfPeaks must be >=1; Set to NumOfPeaks=1 by default !!!\e[0m"<<endl;
				NumOfPeaks =1;
			}

			fsample = new TF1("fsample",this,&funcS::fitfunc,0,25e6,NumOfPeaks*2,"1funcS","1fitfunc");
			OldNPeaks = NumOfPeaks;
			//sleep(1);
		}

		bool Draw(TCanvas* c_todraw=NULL,int Padindex=4){
			if(Padindex>4 || Padindex<1){cout<<"Peak index at canvas is wrong [1,4]"<<endl; return false;}

			if(c_todraw!=NULL){
				c_todraw->cd(Padindex)->SetEditable(kTRUE);
			}
			else{cout<<"canvas is not available!!!!"<<endl; return false;}

			if(fsample_1p!=NULL){delete fsample_1p; fsample_1p = NULL;}
			fsample_1p = new TF1("fsample_1p",this,&funcS::fitfunc_cpy,0,25e6,2,"1funcS","1fitfunc_cpy");
			fsample_1p->SetParameter(0,sAmp);
			fsample_1p->SetParameter(1,sPeakCenter);
			fsample_1p->Draw("same");
			c_todraw->cd(Padindex)->Modified();
			c_todraw->cd(Padindex)->Update();
			if(Padindex==2)c_todraw->cd(Padindex)->SetEditable(kFALSE);
			return true;
		}

		bool Draw_subline(TCanvas* c_todraw=NULL,int Padindex=4){
			if(Padindex>4 || Padindex<1){cout<<"Peak index at canvas is wrong [1,4]"<<endl; return false;}

			if(c_todraw!=NULL){
				c_todraw->cd(Padindex)->SetEditable(kTRUE);
			}
			else{cout<<"canvas is not available!!!!"<<endl; return false;}

			for(int i=0;i<10;i++){ if(fresult[i]!=NULL)delete fresult[i]; fresult[i] = NULL; }

			for(int i=0;i<OldNPeaks && i<MaxNPeaks;i++){
				int kala[]={633,808,799,417,433,600,617};
				//if(fresult[i]!=NULL) delete fresult[i];
				fresult[i] = new TF1(Form("fresult_%d",i+1),this,&funcS::fitfunc_cpy,0,25e6,2,"1funcS","1fitfunc_cpy");
				fresult[i]->SetParameter(0,Amp[i]);
				fresult[i]->SetParameter(1,tof_center[i]);
				fresult[i]->SetLineColor(kala[i%7]);
				fresult[i]->SetLineStyle(2);
				fresult[i]->Draw("same");
			}

			c_todraw->cd(Padindex)->Modified();
			c_todraw->cd(Padindex)->Update();
			if(Padindex==2)c_todraw->cd(Padindex)->SetEditable(kFALSE);
			return true;
		}


		double* MC_err_x(TH1D* _hx_in = NULL, TF1* _fsample_in = NULL, int NLoops=10,string fitopt="LMEQ"){  // NOTE: this will not update the errors of fitting result automatically
			TF1* fsim_x = new TF1("fsim_x",this,&funcS::fitfunc,0,25e6,OldNPeaks*2,"1funcS","1fitfunc");
			fsim_x->SetParameters(_fsample_in->GetParameters());
			fsim_x->SetParErrors(fsim_x->GetParErrors());

			static double* return_err = NULL;
			if(return_err != NULL) delete[] return_err;
			return_err = new double[OldNPeaks];


			double** temfitresult = new double*[OldNPeaks];


			for(int ip=0;ip<OldNPeaks;ip++){
				temfitresult[ip] = new double[NLoops];
			}


			TH1D* _hx_in_cp = dynamic_cast<TH1D*>(_hx_in->Clone("_hx_in_cp"));
			int Nstatistics = _hx_in->Integral(1,_hx_in->GetNbinsX());

			for(int index=0;index<NLoops;index++){
				static int ProcessRatio=0;
				if(index==0) ProcessRatio=0;
				if((int)index*100.0/NLoops>=ProcessRatio){
					cout<<"\e[1;37m MC_sim progress = "<<ProcessRatio<<" %\e[0m"<<"\r"<<flush;
					ProcessRatio+=5;
				}else if(index== NLoops-1){
					cout<<"\e[1;37m MC_sim progress = 100 %\e[0m"<<"\r"<<endl;
					cout<<endl;
				}

				for(int ip=0;ip<OldNPeaks;ip++){
						fsim_x->SetParLimits(ip*2,0,fsim_x->GetParameter(ip*2) * 3 );
						double temPx = fsim_x->GetParameter(ip*2+1);
						fsim_x->SetParLimits(ip*2+1,temPx-2*sFWHM, temPx+2*sFWHM);
				}

				_hx_in_cp->Reset();
				_hx_in_cp->FillRandom(_hx_in,Nstatistics);

/*::ctest.cd();
_hx_in_cp->Draw();*/

				double fit_L=0;
				double fit_R=0;
				int h_in_Nbins=_hx_in_cp->GetNbinsX();

				fit_L = TMath::Max(tof_center[0]-range_L,_hx_in_cp->GetBinCenter(1));
				fit_R = TMath::Min(tof_center[OldNPeaks-1]+range_R,_hx_in_cp->GetBinCenter(h_in_Nbins));

				if(fitopt.find("Q")==string::npos || fitopt.find("q")==string::npos) fitopt+="Q";

				for(int i=0;i<(50+10*(OldNPeaks-1));i++){ // star 50 fitting for each histogram
				  //  _hx_in_cp->Fit(fsim_x,fitopt.c_str(),"",fit_L,fit_R);
//if(i==(50+10*(OldNPeaks-1))-1) _hx_in_cp->Fit(fsim_x,"LMEQ","",fit_L,fit_R);
//else _hx_in_cp->Fit(fsim_x,fitopt.c_str(),"",fit_L,fit_R);
						if(i<20){
							if(i<5){
								_hx_in_cp->Fit(fsim_x,(fitopt+"N").c_str(),"N",fit_L,fit_R);
							}
							else if(i%2==0){
								_hx_in_cp->Fit(fsim_x,"MEQN","N",fit_L,fit_R); // chi-square
							}
							else{
								_hx_in_cp->Fit(fsim_x,"LMEQN","N",fit_L,fit_R); // likelyhood
							}

						}
						else{
							if(i==(50+10*(OldNPeaks-1))-1) _hx_in_cp->Fit(fsim_x,(fitopt+"N").c_str(),"N",fit_L,fit_R); // for draw fitting on histogram in debug mode
							else _hx_in_cp->Fit(fsim_x,(fitopt+"N").c_str(),"N",fit_L,fit_R);
						}

				}// end of 50 fittings 

				for(int ip=0;ip<OldNPeaks;ip++){
					temfitresult[ip][index]=fsim_x->GetParameter(ip*2+1);
				}
/*
::ctest.Modified();
::ctest.Update();
cin>>ctemp;*/


			}// end of Nloop for fitting Nloop histograms


			for(int ip=0;ip<OldNPeaks;ip++){
				return_err[ip] = TMath::RMS(NLoops,temfitresult[ip]);
			}



			for(int ip=0;ip<OldNPeaks;ip++){ delete[] temfitresult[ip]; }
			delete[] temfitresult;
			delete fsim_x;
			delete _hx_in_cp;
			return return_err;
		}




		void Fit(TH1D* h_in, double _range_L, double _range_R, string fitopt="LMEQ"){ //=> fit 50 times
			for(int i=0;i<OldNPeaks  && i<MaxNPeaks;i++){
				fsample->SetParameter(i*2,Amp[i]);
				fsample->SetParameter(1+i*2,tof_center[i]);
			}

			double fit_L=0;
			double fit_R=0;
			int h_in_Nbins=h_in->GetNbinsX();

			for(int i=0;i<(50+10*(OldNPeaks-1));i++){ 
				if(range_L==-1 || range_R==-1 || FreeRange){// new fitting
					if(i==0)cout<<"\e[1;33m"<<"free range"<<"\e[0m"<<endl;
					bool atMinimum=false; // _range_L < sPeakCenter-sample_range_L
					bool atMaximum=false; // _range_R > sample_range_R - sPeakCenter
					double left_width =tof_center[0]-_range_L; 
					double right_width= -(tof_center[OldNPeaks-1]-_range_R);


					if(left_width > (sPeakCenter-sample_range_L)){
						range_L=sPeakCenter-sample_range_L;  // left side fitting width
						atMinimum = true;
					}

					if(right_width > (sample_range_R-sPeakCenter)){
						range_R=sample_range_R-sPeakCenter;   // right side fitting width
						atMaximum = true;
					}

					if(atMinimum && atMaximum){;}
					else if(atMinimum){
						range_R=sample_range_R-sPeakCenter;
					}
					else if(atMaximum){
						range_L=sPeakCenter-sample_range_L;
					}
					else{ // keep the ratio as sample_range_L : sample_range_R;
						if(left_width>right_width){ 
							right_width=left_width /( (sPeakCenter-sample_range_L)/(sample_range_R-sPeakCenter) ); 
						}
						else{
							left_width =right_width *( (sPeakCenter-sample_range_L)/(sample_range_R-sPeakCenter) ); 
						}

						range_L = left_width;
						range_R = right_width;
					}

					
				}


				int NpeaksNow = TMath::Min(OldNPeaks,MaxNPeaks); // make sure the peak sequence is correct
                double * temPx = new double[NpeaksNow];
                double * temPy = new double[NpeaksNow];
                int * sequencelist = new int[NpeaksNow];

				TMath::Sort(NpeaksNow,tof_center,sequencelist,kFALSE);

				for(int ip = 0; ip<NpeaksNow; ip++){
					temPx[ip] = tof_center[sequencelist[ip]];
					temPy[ip] = Amp[sequencelist[ip]];
				}


				for(int ip = 0; ip<NpeaksNow; ip++){
					tof_center[ip] = temPx[ip];
					Amp[ip] = temPy[ip];
					fsample->SetParameter(ip*2,Amp[ip]);
					fsample->SetParameter(ip*2+1,tof_center[ip]);
				}

				delete[] temPx;
				delete[] temPy;
				delete[] sequencelist;

			// fit 50 times
				fit_L = TMath::Max(tof_center[0]-range_L,h_in->GetBinCenter(1));
				fit_R = TMath::Min(tof_center[OldNPeaks-1]+range_R,h_in->GetBinCenter(h_in_Nbins));

				if(i<20){
					if(OldNPeaks>1){
						for(int index=0;index<OldNPeaks  && index<MaxNPeaks;index++){ // set the limit of position X for each peak
							if(index==0) fsample->SetParLimits(1+index*2,fit_L,tof_center[1]);
							else if(index==OldNPeaks-1) fsample->SetParLimits(1+index*2, tof_center[OldNPeaks-2], fit_R);
							else fsample->SetParLimits(1+index*2, tof_center[index-1], tof_center[index+1]);

							double countUpper=0;
							double countLower=0;
							double tem_high = fsample->GetParameter(index*2);
							countLower =TMath::Max( tem_high-TMath::Sqrt(tem_high)*2 , 0.);
							countUpper = tem_high+TMath::Sqrt(tem_high)*2;
							fsample->SetParLimits(index*2,countLower,countUpper);
		fsample->SetParLimits(index*2+1,tof_center[index]-sFWHM*2.5,tof_center[index]+sFWHM*2.5);
						}
					}else{ // single peak case
							double countUpper=0;
							double countLower=0;
							double tem_high = fsample->GetParameter(0);
							countLower =TMath::Max( tem_high-TMath::Sqrt(tem_high)*2 , 0.);
							countUpper = tem_high+TMath::Sqrt(tem_high)*2;
							fsample->SetParLimits(0,countLower,countUpper);
							fsample->SetParLimits(1,fit_L,fit_R);
		fsample->SetParLimits(1,tof_center[0]-sFWHM*2.5,tof_center[0]+sFWHM*2.5);
					}
				}else{ // remove limits of paras
					for(int index=0;index<OldNPeaks  && index<MaxNPeaks;index++){
						fsample->SetParLimits(index*2,0,fsample->GetParameter(index*2)*10);  // Amp limit for each peak [0,10 x current value]
						fsample->SetParLimits(index*2+1,0.,0.); // free position of each peak
		fsample->SetParLimits(index*2+1,tof_center[index]-sFWHM*2.5,tof_center[index]+sFWHM*2.5);
					}

				}

			
				if(i==0){
					if(fit_L==h_in->GetBinCenter(1)){
						cout<<"\e[1;2m"<<"Warning: fitting range is limited to the \"Left\" edge of histogram!!!!"<<"\e[0m"<<endl;
					}
					if(fit_R==h_in->GetBinCenter(h_in_Nbins)){
						cout<<"\e[1;2m"<<"Warning: fitting range is limited to the \"Right\" edge of histogram!!!!"<<"\e[0m"<<endl;
					}
				}


				if(i<20){
					if(i<5){
						h_in->Fit(fsample,(fitopt+"N").c_str(),"N",fit_L,fit_R);
					}
					else if(i%2==0){
						h_in->Fit(fsample,"MEQN","N",fit_L,fit_R); // chi-square
					}
					else{
						h_in->Fit(fsample,"LMEQN","N",fit_L,fit_R); // likelyhood
					}

				}
				else{
					h_in->Fit(fsample,(fitopt+"N").c_str(),"N",fit_L,fit_R);
				}


				for(int index=0;index<OldNPeaks  && index<MaxNPeaks;index++){
					tof_center[index]=fsample->GetParameter(index*2+1);
					Amp[index] = fsample->GetParameter(index*2);
				}

			}// end of fit 50 times



			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'Q'),fitopt.end());
			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'q'),fitopt.end());

				int NpeaksNow = TMath::Min(OldNPeaks,MaxNPeaks);
                double * temPx = new double[NpeaksNow];
                double * temPy = new double[NpeaksNow];
                double * temPx_err = new double[NpeaksNow];
                int * sequencelist = new int[NpeaksNow];

				TMath::Sort(NpeaksNow,tof_center,sequencelist,kFALSE);

				for(int ip = 0; ip<NpeaksNow; ip++){
					temPx[ip] = tof_center[sequencelist[ip]];
					temPy[ip] = Amp[sequencelist[ip]];
				}


				for(int ip = 0; ip<NpeaksNow; ip++){
					tof_center[ip] = temPx[ip];
					Amp[ip] = temPy[ip];
				}



			fit_L = TMath::Max(tof_center[0]-range_L,h_in->GetBinCenter(1));
			fit_R = TMath::Min(tof_center[OldNPeaks-1]+range_R,h_in->GetBinCenter(h_in_Nbins));

			h_in->Fit(fsample,fitopt.c_str(),"",fit_L,fit_R);

			for(int index=0;index<OldNPeaks  && index<MaxNPeaks;index++){ // save fitting result
				tof_center[index]=fsample->GetParameter(index*2+1);
				Amp[index] = fsample->GetParameter(index*2);
				tof_center_err[index] = fsample->GetParError(index*2+1);
				tof_center_err[index] = TMath::Sqrt(tof_center_err[index]*tof_center_err[index]+sPeakCenter_err*sPeakCenter_err);

			}


			TMath::Sort(NpeaksNow,tof_center,sequencelist,kFALSE); // sort fitting result

			for(int ip = 0; ip<NpeaksNow; ip++){
				temPx[ip] = tof_center[sequencelist[ip]];
				temPy[ip] = Amp[sequencelist[ip]];
				temPx_err[ip] = tof_center_err[sequencelist[ip]];
			}


			for(int ip = 0; ip<NpeaksNow; ip++){
				tof_center[ip] = temPx[ip];
				Amp[ip] = temPy[ip];
				tof_center_err[ip] = temPx_err[ip];
			}


			range_L = tof_center[0] - fit_L;
			range_R = fit_R-tof_center[OldNPeaks-1];
			
			printf("width for fit: %.1f [ns] ==>[%.2f,%.2f]\n",fit_R-fit_L,fit_L-tof_center[0],fit_R-tof_center[OldNPeaks-1]);

if(fit_R-fit_L<0){
	cout<<endl;
	cout<<endl;
	cout<<"\e[1;31m Fitting range Error!!!!!!!!!!!!!!!!!!!!\e[0m"<<endl;
	cout<<"\e[1;33m First and last peak\e[0m"<<endl;
	printf("%.4f \t %.4f\n",tof_center[0],tof_center[OldNPeaks-1]);
	printf("fit from %.4f   to  %.4f\n",fit_L, fit_R);
	cout<<endl;
}

			// update the params of fsample as sequence
			double * tempars= new double[fsample->GetNpar()];
			double * tempars_err = new double[fsample->GetNpar()];
			fsample->GetParameters(tempars);
			const double * pars_err_ptr = fsample->GetParErrors();
			for(int ipar=0;ipar<fsample->GetNpar();ipar++) tempars_err[ipar] = pars_err_ptr[ipar];

			for(int i=0;i<OldNPeaks && i<MaxNPeaks;i++){
				fsample->SetParameter( 2*i, tempars[sequencelist[i]*2] );
				fsample->SetParError( 2*i, tempars_err[sequencelist[i]*2] );
				fsample->SetParameter( 2*i+1, tempars[sequencelist[i]*2+1] );
				fsample->SetParError( 2*i+1, tempars_err[sequencelist[i]*2+1] );
			}



			if(useMC)	 MC_err_ptr = MC_err_x(h_in,fsample,MC_sim_counts,fitopt);


            delete[] temPx ;
            delete[] temPy;
            delete[] temPx_err;
            delete[] sequencelist;
			delete[] tempars;
			delete[] tempars_err;

			FreeRange=false;

		}


		funcS(){
			sAmp = 0;   // Amp of sampled peak
			sPeakCenter=0; // peak center of sampled peak;
			sPeakCenter_err=0;
			sFWHM=0;
			sSigma=10;
			for(int i=0;i<10;i++){Amp[i]=0;tof_center[i]=0;tof_center_err[i]=0;fresult[i]=NULL;}
			sample_range_L=0;
			sample_range_R=0;
			MaxNPeaks=10;

			spl=NULL;
			fsample=NULL;
			fsample_1p=NULL;
			FreeRange=true;  // FreeRange=false, fix left and right weight ratio
			range_L=-1;  //real raange for fit
			range_R=-1;
			bins_width=0;
			Nbins=0;
			h_sample=NULL;
			MC_sim_counts=50;
			useMC=false;
			MC_err_ptr=NULL;
		}

		~funcS(){
			if(spl!=NULL) delete spl;
			if(fsample!=NULL) delete fsample;
			if(fsample_1p!=NULL) delete fsample_1p;
			if(h_sample!=NULL) delete h_sample;
			for(int i=0;i<10;i++){if(fresult[i]!=NULL)delete fresult[i];}
		}

};

int funcS::smoothlevel=2;
int funcS::NumOfPeaks=1;

#endif