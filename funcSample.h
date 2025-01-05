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

class funcS{
	private:
		vector<double> X_tof;
		vector<double> Y_count;
		double sAmp;   // Amp of sampled peak
		double sPeakCenter; // peak center of sampled peak;
		double sPeakCenter_err;
		double sFWHM;
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
		
	public:
		static int smoothlevel;
		static int NumOfPeaks;
		bool FreeRange;  // FreeRange=false, fix left and right weight ratio
		double range_L;  // left width for fit
		double range_R; // right width for fit
		double bins_width;
		int Nbins;
		TH1D* h_sample;


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
					h_sample->Fit(gfit,"MELQN","NQ",sPeakCenter-sFWHM*0.5,sPeakCenter+sFWHM*0.5);
				}
				else{
					if(half_sigma<bins_width*15){ //half_sigma<bins_width*15
						h_sample->Fit(gfit,"MELQN","NQ",sPeakCenter-half_sigma,sPeakCenter+half_sigma); 
						couts_for_fit = h_sample->Integral(h_sample->GetXaxis()->FindBin(sPeakCenter-half_sigma),h_sample->GetXaxis()->FindBin(sPeakCenter+half_sigma));
					}
					else{ 
						h_sample->Fit(gfit,"MELQN","NQ",sPeakCenter-bins_width*15,sPeakCenter+bins_width*15); 
						couts_for_fit = h_sample->Integral(h_sample->GetXaxis()->FindBin(sPeakCenter-bins_width*15),h_sample->GetXaxis()->FindBin(sPeakCenter+bins_width*15));
					}
				}
			}

			double error_estimate = sFWHM/(2.36*TMath::Sqrt(couts_for_fit));

			sAmp = gfit->GetParameter(0);
			sPeakCenter = gfit->GetParameter(1);
			sPeakCenter_err =gfit->GetParError(1);
			h_sample->Scale(1/sAmp);
			cout<<"\e[1;33m"<<"Sampling result:"<<"\e[0m"<<endl;
			printf("Amp = %.2f \t Tof= %.4f(%.4f), err_estimated by counts = %.4f\n",sAmp,sPeakCenter,sPeakCenter_err,error_estimate);
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
				if(X_tof[0]<=X_convert && X_convert<=X_tof[Nbins-1]){value +=  par[i*2]* spl->Eval(X_convert);}
				else{value+=0;}
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
			if(fsample != NULL){delete fsample;}
			fsample = new TF1("fsample",this,&funcS::fitfunc,0,25e6,NumOfPeaks*2,"1funcS","1fitfunc");
			OldNPeaks = NumOfPeaks;
			
		}

		bool Draw(TCanvas* c_todraw=NULL,int Padindex=4){
			if(Padindex>4 || Padindex<1){cout<<"Peak index at canvas is wrong [1,4]"<<endl; return false;}

			if(c_todraw!=NULL){
				c_todraw->cd(Padindex)->SetEditable(kTRUE);
			}
			else{cout<<"canvas is not available!!!!"<<endl; return false;}

			if(fsample_1p!=NULL) delete fsample_1p;
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

			for(int i=0;i<OldNPeaks && i<MaxNPeaks;i++){
				int kala[]={633,808,799,417,433,600,617};
				if(fresult[i]!=NULL) delete fresult[i];
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

		void Fit(TH1D* h_in, double _range_L, double _range_R, string fitopt="LMEQ"){ //=> fit 50 times
			for(int i=0;i<OldNPeaks  && i<MaxNPeaks;i++){
				fsample->SetParameter(i*2,Amp[i]);
				fsample->SetParameter(1+i*2,tof_center[i]);
			}

			double fit_L=0;
			double fit_R=0;
			int h_in_Nbins=h_in->GetNbinsX();

			for(int i=0;i<50;i++){ 
				if(range_L==-1 || range_R==-1 || FreeRange){// new fitting
					if(i==0)cout<<"\e[1;33m"<<"free range"<<"\e[0m"<<endl;
					bool atMinimum=false; // _range_L < sPeakCenter-sample_range_L
					bool atMaximum=false; // _range_R > sample_range_R - sPeakCenter
					double left_width =tof_center[0]-_range_L; 
					double right_width= -(tof_center[OldNPeaks-1]-_range_R);


					if(left_width > (sPeakCenter-sample_range_L)){
						range_L=sPeakCenter-sample_range_L;
						atMinimum = true;
					}

					if(right_width > (sample_range_R-sPeakCenter)){
						range_R=sample_range_R-sPeakCenter;
						atMaximum = true;
					}

					if(atMinimum && atMaximum){;}
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


			// fit 50 times
				fit_L = TMath::Max(tof_center[0]-range_L,h_in->GetBinCenter(1));
				fit_R = TMath::Min(tof_center[OldNPeaks-1]+range_R,h_in->GetBinCenter(h_in_Nbins));
			
				if(i==0){
					if(fit_L==h_in->GetBinCenter(1)){
						cout<<"\e[1;2m"<<"Warning: fitting range is limited to the \"Left\" edge of histogram!!!!"<<"\e[0m"<<endl;
					}
					if(fit_R==h_in->GetBinCenter(h_in_Nbins)){
						cout<<"\e[1;2m"<<"Warning: fitting range is limited to the \"Right\" edge of histogram!!!!"<<"\e[0m"<<endl;
					}
				}

				h_in->Fit(fsample,fitopt.c_str(),"N",fit_L,fit_R);

				for(int i=0;i<OldNPeaks  && i<MaxNPeaks;i++){
					tof_center[i]=fsample->GetParameter(i*2+1);
				}

			}

			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'Q'),fitopt.end());
			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'q'),fitopt.end());

			h_in->Fit(fsample,fitopt.c_str(),"",fit_L,fit_R);
			printf("width for fit: %.1f [ns] ==>[%.2f,%.2f]\n",fit_R-fit_L,fit_L-tof_center[0],fit_R-tof_center[OldNPeaks-1]);

			for(int i=0;i<OldNPeaks && i<MaxNPeaks;i++){
				Amp[i]=fsample->GetParameter(i*2);
				tof_center[i]=fsample->GetParameter(i*2+1);
				tof_center_err[i] = fsample->GetParError(i*2+1);
				tof_center_err[i] = TMath::Sqrt(tof_center_err[i]*tof_center_err[i]+sPeakCenter_err*sPeakCenter_err);
			}

			FreeRange=false;

		}


		funcS(){
			sAmp = 0;   // Amp of sampled peak
			sPeakCenter=0; // peak center of sampled peak;
			sPeakCenter_err=0;
			sFWHM=0;
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