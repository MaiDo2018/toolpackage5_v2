#ifndef _FUNCN_H_
#define _FUNCN_H_
#include "FitHelper.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>



using namespace std;

//TCanvas ctest2;
class funcN{
private:
		double sAmp;   // Amp of sampled peak
		double sPeakCenter; // peak center of sampled peak;
		double sPeakCenter_err;
		double sFWHM;
		double sSigma;  // using 0.6826 define +/- sigma region
		double Amp[10];
		double tof_center[10];
		double tof_center_err[10];

		TF1* func_ref;
		TF1* func_x;
		TF1* fresult[10];
		TF1* func_ref_cp;
		TF1* func_x_cp;

		int Norder;
		int NumOfPeaks;
		int MainpeakIndex;
		int pars_per_peak;
		double* par_ref;
		int MC_sim_counts;

public:
		static int NPs2Set;
		static int Norder2Set;
		bool LockPars; // lock shape pars; true ==> set to lock; false ==> release
		bool FreeRange;
		bool AutoUpdateMainPeakIndex;
		bool useMC; // MC simulation for getting fitting error
		double sample_range_L; // half range for fit
		double sample_range_R;
		bool EnableFunc_refInitialize;  // set to false to skip the parameter initialize for reference function


		void SetMC_sim_counts(int Ncounts){
			MC_sim_counts = Ncounts;
		}

		int GetMC_sim_counts(bool IsPrint=true){ if(IsPrint){cout<<"MC simulation Nloops = "<<MC_sim_counts<<endl;} return MC_sim_counts;	}


		int GetNumberOfPeaks(){return NumOfPeaks;}
		bool SetNumberOfPeaks(int NPeaks){
			if(NPeaks>10 || NPeaks<1){cout<<"Npeaks allowed = 1-10"<<endl; return false;}
			else{ NumOfPeaks = NPeaks; return true;}
		}

		int GetNorder(){return Norder;}
		void SetNorder(int _Norder){Norder = _Norder;}

		int GetMainPeakIndex(){return MainpeakIndex;}


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

		double GetFWHM(){
			return sFWHM;
		}

		double UpdateSigma(){
			if(func_ref != NULL){
				double Area_x = func_ref->Integral(sPeakCenter-20*sFWHM, sPeakCenter+20*sFWHM);
				double LowLimitX = sPeakCenter-2*sFWHM;
				double UpLimitX = sPeakCenter;
				double MiddleX = 0;
				double temSigma = 1;
				while(1){
				  MiddleX = (LowLimitX+UpLimitX)*0.5;
				  temSigma = sPeakCenter - MiddleX;
				  double Area_ratio_now = func_ref->Integral(MiddleX,sPeakCenter+temSigma)/Area_x;
				  double Area_diff = Area_ratio_now - 0.6826;
				  if(  TMath::Abs(Area_diff) < 0.00001) break;
				  else{
				  	if(Area_diff >0 ) LowLimitX = MiddleX;
				  	else UpLimitX =MiddleX;
				  }
				}

				return temSigma;

			}else{
				cout<<"\e[1;32m Fitting function of reference is not available, return sigma = 10\e[0m"<<endl;
				return 10;
			}
		}

		double GetSigma(){
			if(sSigma==10)cout<<"\e[1;32m Fitting function of reference is not available, return sigma = 10\e[0m"<<endl;
			return sSigma;
		}

		void Initial_par_ref(TH1D* _hin, double left_range, double right_range){
			if(!EnableFunc_refInitialize && func_ref == NULL){
				EnableFunc_refInitialize=true;
				cout<<endl;
				cout<<"\e[1;31m Warnning: Reference function not exist; Automatically set EnableFunc_refInitialize-> true  by default \e[0m"<<endl;
			}

			if(!EnableFunc_refInitialize) return;
			if(func_ref!=NULL){delete func_ref; func_ref=NULL;}
			if(par_ref!=NULL) delete[] par_ref;
			pars_per_peak = 5+Norder*2;
			par_ref = new double[pars_per_peak];
			int binl = _hin->FindBin(left_range);
			int binr = _hin->FindBin(right_range);
			int binmax=0;
			int high = 0;
			for(int ibin=binl;ibin<=binr;ibin++){
				if(_hin->GetBinContent(ibin)>high){high = _hin->GetBinContent(ibin); binmax=ibin;}
			}

			par_ref[0]=high;
			par_ref[1]=_hin->GetBinCenter(binmax);

			//******* get FWHM *******************

			double FWHM_L=0, FWHM_R =0;
			double bin_i_y1;
			double bin_i_y2;
			double candidate_x;
			double height_half=0.5*high;
			double bins_width = _hin->GetBinWidth(1);
			TH1D* _hin_cp = (TH1D*)_hin->Clone();
			_hin_cp->Smooth(4);

			for(int i=1;i<_hin_cp->GetNbinsX();i++){
				bin_i_y1 = _hin_cp->GetBinContent(i);
				bin_i_y2 = _hin_cp->GetBinContent(i+1);
				candidate_x = (height_half - bin_i_y1) / (bin_i_y2 - bin_i_y1) * bins_width + _hin_cp->GetBinCenter(i);
				if(i<binmax){
					if(bin_i_y1<height_half && bin_i_y2 >height_half){ 	FWHM_L = candidate_x;	}
				}
				else{
					if(bin_i_y1>height_half && bin_i_y2 <height_half){		FWHM_R = candidate_x; break;	}
				}
			}

			//**********************************
			sFWHM = FWHM_R-FWHM_L;
			_hin_cp->GetXaxis()->SetRangeUser(par_ref[1]-sFWHM,par_ref[1]+sFWHM);


			double sigma_ref = _hin_cp->GetStdDev();
			par_ref[2] = sigma_ref; par_ref[3] = sigma_ref; par_ref[4] = -sigma_ref;
			double ka = 1/sigma_ref;
			for(int iorder=1;iorder<=Norder;iorder++){
					par_ref[5+2*(iorder-1)]= TMath::Power(-1,iorder-1)*ka; //odd order==> position; even order=>negative
					par_ref[6+2*(iorder-1)]= TMath::Power(-1,iorder-1)*((iorder-1)/2+2)*sigma_ref;//odd order==> position; even order=>negative
			}

			delete _hin_cp;


			for(int i=0;i<pars_per_peak;i++){ // print the initial pars
				if(i<5){
					if(i==0) printf("High = %.2f\n",par_ref[i]);
					if(i==1) printf("center = %.4f\n",par_ref[i]);
					if(i==2) printf("sigma = %.4f\n",par_ref[i]);
					if(i==3) printf("tR = %.4f\n",par_ref[i]);
					if(i==4) printf("tL = %.4f\n",par_ref[i]);
				}
				else{
					if(i%2 ==1) printf("k%d = %.4f\n",(i-5)/2+1,par_ref[i]); // odd index ==> k
					else printf("t%d = %.4f\n",(i-5)/2+1,par_ref[i]); // even index ==> t
				}
			}

		}

		double GetsPeakCenter(){return sPeakCenter;}

		double GetsPeakCenter_err(){return sPeakCenter_err;}

		TF1* GetFuncRef(){return func_ref;}

		TF1* GetFuncX(){return func_x;}


		double A_n(double t,double tm,int iorder, double* par){// iorder>=1 k1t1(++),k2t2(--),k3t3(++) ==> R, L, R
			double tR = par[3], tL=par[4], sigma = par[2], k1=par[5],t1=par[6], k2=par[7],t2=par[8];
			double A0,k0;
			if(t>tm){// iorder==1,3,5,7

					A0 = par[0]* TMath::Exp(tR*(tR-t1)/(2*sigma*sigma));
					k0 = tR/(2*sigma*sigma);
				if(iorder==1){
					return A0* TMath::Exp((-k0+par[5+2*(iorder-1)])*par[6+2*(iorder-1)]);
				}
				else{
					return A_n(t,tm,iorder-2,par) * TMath::Exp((-par[5+2*(iorder-2-1)]+par[5+2*(iorder-1)]) * par[6+2*(iorder-1)]);
				}
			}
			else{//iorder==2,4,6,8
					A0 = par[0] * TMath::Exp(tL*(tL-t2)/(2*sigma*sigma));
					k0 = tL/(2*sigma*sigma);
				if(iorder==2){
					return A0* TMath::Exp((-k0 + par[5+2*(iorder-1)])*par[6+2*(iorder-1)]);
				}
				else{
					return A_n(t,tm,iorder-2,par) * TMath::Exp((-par[5+2*(iorder-2-1)]+par[5+2*(iorder-1)]) * par[6+2*(iorder-1)]);

				}
			}


		}

		int WhichOrder(double t,double* par){
			double tR = par[3], tL=par[4];
			int iorder=0;
			if(Norder==0) return iorder;
			if(t>par[1]+tR){
				for(iorder=1;iorder<=Norder;iorder+=2){
					if(t<=par[1]+par[5+2*(iorder-1)+1]){
						if(iorder==1)return 0;
						else return iorder-2;
					}
				}
				return iorder-2;
			}
			else{//t<tL
				for(iorder=2;iorder<=Norder;iorder+=2){
					if(t>=par[1]+par[5+2*(iorder-1)+1]){
						if(iorder==2)return 0;
						else return iorder-2;
					}
				}
				return iorder-2;
			}
		}

		Double_t gas(double*x,double*par){
		   return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2]));  

		}

		Double_t exR(double*x,double*par){
			return par[0]*TMath::Exp(par[3]*(par[3]-2*(x[0]-par[1]))/(2*par[2]*par[2]));   // right long tail
		 
		}

		Double_t exL(double*x,double*par){
			return par[0]*TMath::Exp(par[4]*(par[4]-2*(x[0]-par[1]))/(2*par[2]*par[2]));// left long tail
		 
		}

		Double_t fitfunc(double*x,double*par){  // par sequence=> high, cento, sigma, right tail tc, left tail tc
		   /* for right long tail, all parameters>0 */
		    if(x[0]<=par[1]+par[3]&&x[0]>=par[1]+par[4]) return gas(x,par); //+ linearbg(x,&par[4]);

		    int getorder = WhichOrder(x[0],par);

		    if(getorder==0){//getorder==0
			/* for left long tail , right tail para>0*/
		    	if(x[0]>par[1]+par[3]) return exR(x,par); //+ linearbg(x,&par[4]);
		    
		 	/* for left long tail , left tail para<0*/
		    	if(x[0]<par[1]+par[4]) return exL(x,par);
		    }
		    //if(x[0]>par[1]-par[4]) return gas(x,par);
		    else{// get order>0

			return A_n(x[0],par[1],getorder, par)*TMath::Exp(-par[5+2*(getorder-1)]*(x[0]-par[1]));
		    }

		    return 0;
		   
		}

		void UpdateMainPeakIndexAuto(){
			double amp_tem=-1;
			int MP_index=0;
			for(int i=0;i<NumOfPeaks;i++){
				if(Amp[i]>amp_tem){
					amp_tem = Amp[i];
					MP_index = i;
					//cout<<"peak "<<i<<"  Amp "<<Amp[i]<<"  tof "<<tof_center[i]<<endl;
				}
			}
			MainpeakIndex = MP_index+1;  // MainpeakIndex [1,10]
			
		}

		bool IsChange_MainpeakIndex(){
			static int MainpeakIndex_record=-10;
			if(MainpeakIndex_record==-10){MainpeakIndex_record = MainpeakIndex;}
			else if(MainpeakIndex_record != MainpeakIndex){
					MainpeakIndex_record = MainpeakIndex;
					return true;
			}
			return false;
		}


		int Paras_offset_cal(int ip,int _MainpeakIndex){// ip from 0; _MainpeakIndex from1
			
			if(ip<_MainpeakIndex) return ip*2;
			else{
				return (_MainpeakIndex-1)*2+pars_per_peak+(ip-_MainpeakIndex)*2;
			}

		}

		void ParaDistributor(int ip, int _MainpeakIndex, double* inParas_all, double* outParas){
			int ip_par_offset = Paras_offset_cal(ip,_MainpeakIndex);
			int mainp_par_offset = Paras_offset_cal(_MainpeakIndex-1,_MainpeakIndex);
			for(int ipar=0;ipar<pars_per_peak;ipar++){
				if(ipar<2) outParas[ipar] = inParas_all[ip_par_offset+ipar];
				else outParas[ipar] = inParas_all[mainp_par_offset+ipar];
			}

		}

		void UpdateParsFromTF(){
			int ip_par_offset=0;
			for(int ip=0;ip<NumOfPeaks;ip++){
				ip_par_offset = Paras_offset_cal(ip,MainpeakIndex); // ip from 0; MainpeakIndex from 1
				Amp[ip] = func_x->GetParameter(ip_par_offset);
				tof_center[ip] = func_x->GetParameter(ip_par_offset+1);
				tof_center_err[ip] = func_x->GetParError(ip_par_offset+1);
				printf("cento_%d:  %.4f (%.4f)\n",ip+1,tof_center[ip],tof_center_err[ip]);
			}
		}


		void LockPars_Ref(){
			for(int i=2;i<pars_per_peak;i++){
				func_ref->FixParameter(i,func_ref->GetParameter(i));
				func_ref_cp->FixParameter(i,func_ref_cp->GetParameter(i));
			}
		}


		void LockPars_X(){ // lock shape pars of func_x
			if(func_x==NULL){cout<<"function for fitting is not existed! MakeFitFunc_N() first!!!"<<endl; return;}
			double para_offset = Paras_offset_cal(MainpeakIndex-1,MainpeakIndex);
			for(int i=2;i<pars_per_peak;i++){
				func_x->FixParameter(para_offset+i,func_x->GetParameter(para_offset+i));
				func_x_cp->FixParameter(para_offset+i,func_x_cp->GetParameter(para_offset+i));
			}
		}

		void Reles_X(){
			if(func_x==NULL){cout<<"function for fitting is not existed! MakeFitFunc_N() first!!!"<<endl; return;}
			for(int i=0;i<(NumOfPeaks-1)*2+pars_per_peak;i++){
				func_x->ReleaseParameter(i);
				func_x_cp->ReleaseParameter(i);

			}
			cout<<"\e[1;31m"<<"Warning: shaping parameters are RELEASED!!!! "<<"\e[0m"<<endl;
		}

		double fitfunc_N(double* x, double* par){
			double* paras_ip = new double[pars_per_peak];
			double value =0;
			for(int ip=0;ip<NumOfPeaks;ip++){
				ParaDistributor(ip,MainpeakIndex,par,paras_ip);
				value+= fitfunc(x,paras_ip);
			}
			delete[] paras_ip;
			return value;
		}

		void MakeFitFunc_ref(){
			if(!EnableFunc_refInitialize){return;}			

			if(func_ref!=NULL){delete func_ref; delete func_ref_cp;}
			func_ref = new TF1("fextend_ref",this,&funcN::fitfunc,0,25e6,pars_per_peak,"1func_ref","1fitfunc_ref");
			func_ref_cp = new TF1("fextend_ref_cp",this,&funcN::fitfunc,0,25e6,pars_per_peak,"1func_ref","1fitfunc_ref");

			func_ref->SetParameters(par_ref);
			for(int i=0;i<=Norder;i++){
				if(i==0){
					func_ref->SetParName(3,"tR");
					func_ref->SetParName(4,"tL");
				}
				else{
					if( i%2 !=0){//odd order
						func_ref->SetParName(5+2*(i-1),Form("k%d(+)",i));
						func_ref->SetParName(6+2*(i-1),Form("t%d(+)",i));
					}
					else{
						func_ref->SetParName(5+2*(i-1),Form("k%d(-)",i));
						func_ref->SetParName(6+2*(i-1),Form("t%d(-)",i));
					}

				}
			}
		}

		void MakeFitFunc_N(){//to use auto get MainPeakIndex, initial peaks height in Amp[10] is necessary!!!!
			if(func_ref==NULL){cout<<"Reference is not existed yet. Abort!!!!"<<endl; return;}
			if(func_x!=NULL){ delete func_x; func_x=NULL; delete func_x_cp; func_x_cp = NULL;}
			
			if(!AutoUpdateMainPeakIndex){
				int MP_index=0;
				cout<<endl;
				while(MP_index<1 || MP_index>10){
					cout<<"Please specify the Main peak index[1,10]:"<<endl;
					cin>>MP_index;
				}

				MainpeakIndex = MP_index;
				
			}
			else{UpdateMainPeakIndexAuto();}


			func_x = new TF1("fextend_N",this,&funcN::fitfunc_N,0,25e6,pars_per_peak + 2*(NumOfPeaks-1),"1func_ref","1fitfunc_ref");
			func_x->SetNpx(500);
			func_x_cp = new TF1("fextend_N_cp",this,&funcN::fitfunc_N,0,25e6,pars_per_peak + 2*(NumOfPeaks-1),"1func_ref","1fitfunc_ref");
			func_x_cp->SetNpx(500);

			int paroffset = Paras_offset_cal(MainpeakIndex-1,MainpeakIndex);
			for(int i=0;i<pars_per_peak;i++){
				func_x->SetParameter(paroffset+i,func_ref->GetParameter(i));
				func_x_cp->SetParameter(paroffset+i,func_ref->GetParameter(i));
			}

			for(int i=0;i<=Norder;i++){
				if(i==0){
					func_x->SetParName(paroffset+3,"tR");
					func_x->SetParName(paroffset+4,"tL");
				}
				else{
					if( i%2 !=0){//odd order
						func_x->SetParName(paroffset+5+2*(i-1),Form("k%d(+)",i));
						func_x->SetParName(paroffset+6+2*(i-1),Form("t%d(+)",i));
					}
					else{
						func_x->SetParName(paroffset+5+2*(i-1),Form("k%d(-)",i));
						func_x->SetParName(paroffset+6+2*(i-1),Form("t%d(-)",i));
					}

				}
			}
			
		}


		void MC_ErrEval(TH1D* h_in, char whichfunc, double _fitL, double _fitR, string fitopt="LMEQ"){

			 TF1* fhandle_cp=NULL, *fhandle=NULL;

			 if(whichfunc == 'r'){ fhandle_cp = func_ref_cp; fhandle = func_ref; }
			 else{ fhandle_cp = func_x_cp; fhandle = func_x; }


			 int Ncouts_fill = h_in->Integral(1,h_in->GetNbinsX() );   // h_in->Integral(h_in->FindBin(_fitL),h_in->FindBin(_fitR));

			/* for(int ibin =h_in->FindBin(_fitL); ibin<=h_in->FindBin(_fitR); ibin++){
			 	Ncouts_fill+= (int)fhandle_cp->Eval( h_in->GetBinCenter(ibin) );
			 }*/


			if(Ncouts_fill == 0){cout<<"No counts input for fit!!! Abort"<<endl; return;}
			if(NumOfPeaks==1){	 
				if( sFWHM/2.36/TMath::Sqrt( Ncouts_fill)< 0.1){
			 		cout<<"\e[1;37m Skip MC simulation for intense peak\e[0m"<<endl;
			 		return; // skip strong pick
			 	}
			}

			TH1D* h_in_cp =NULL, *h_in_makeRandom =NULL;
			h_in_cp = dynamic_cast<TH1D*>(h_in->Clone("h_in_cp")) ;
			h_in_makeRandom = dynamic_cast<TH1D*>(h_in_cp->Clone("h_in_makeRandom") );

		//	h_in_makeRandom->Smooth(3);  // seem better to not smooth. Bin with counts less than one will generate nothing !!!!!!!!
			 
			 h_in_cp->Reset();
			 double temL=0,temR=0;
			 fhandle->GetRange(temL,temR);
			 fhandle->SetRange(h_in->GetBinCenter(1),h_in->GetBinCenter(h_in->GetNbinsX()));
			 fhandle->SetNpx(500);

			 vector<double> *fitcenter = new vector<double>[NumOfPeaks];

			 if(fitopt.find("q") == string::npos && fitopt.find("Q") == string::npos) fitopt+="Q";
			 if(fitopt.find("n") == string::npos && fitopt.find("N") == string::npos) fitopt+="N";

			 time_t start = time(0);
			 time_t now;
			 
//char temc;			 
			 for(int nloop=0;nloop<MC_sim_counts;nloop++){
			 	static int ProcessRatio=0;
				if(nloop==0) ProcessRatio=0;
				if((int)nloop*100.0/MC_sim_counts>=ProcessRatio){
					cout<<"\e[1;37m MC_sim progress = "<<ProcessRatio<<" %\e[0m"<<"\r"<<flush;
					ProcessRatio+=5;
				}else if(nloop== MC_sim_counts-1){
					cout<<"\e[1;37m MC_sim progress = 100 %\e[0m"<<endl;
					cout<<endl;
				}

			 	now = time(0);

			 	if(whichfunc == 'r'){
					if(now-start>(int)TMath::Max(60*MC_sim_counts*1./100 , 60.) ){cout<<"\e[1;32m Time out or because too much fitting failure of simulation!! Stop!!\e[0m"<<endl; break;} // >1 min during fitting(fail)
			 	}else{
			 		if(now-start>(int)TMath::Max(60*NumOfPeaks*MC_sim_counts*1./100 , 60.) ){cout<<"Time out or because too much fitting failure of simulation!! Stop!!"<<endl; break;} // >1 min during fitting(fail)
			 	}


				h_in_cp->FillRandom(h_in_makeRandom,Ncouts_fill);
				if(whichfunc == 'r'){
						h_in_cp->GetXaxis()->SetRangeUser(_fitL,_fitR);
				}else{
					int Paroffset_MainPeak = Paras_offset_cal(MainpeakIndex-1,MainpeakIndex);
					double Px_main = fhandle_cp->GetParameter(Paroffset_MainPeak+1);
					h_in_cp->GetXaxis()->SetRangeUser(TMath::Max(_fitL,Px_main-sFWHM*1.5), TMath::Min(_fitR,Px_main+sFWHM*1.5) );
				}


				/* for(int i=0;i<Ncouts_fill;i++){
				 	h_in_cp->Fill(fhandle->GetRandom(h_in->GetBinCenter(1),h_in->GetBinCenter( h_in->GetNbinsX() )) );

				 	//cout<<_fhandle->GetRandom(_fitL,_fitR)<<endl;
				 }*/

				 double t_max = h_in_cp->GetBinCenter(h_in_cp->GetMaximumBin());

				 double t_shift=0;
				 h_in_cp->GetXaxis()->UnZoom();
/*
::ctest2.cd();
h_in_cp->Draw();*/

				 if(whichfunc == 'r'){
				 	fhandle_cp->SetParameter(1,t_max);
				 }
				 else{
				 	int para_offset = Paras_offset_cal(MainpeakIndex-1,MainpeakIndex);
				 	t_shift = t_max - fhandle_cp->GetParameter(para_offset+1);
				 	for(int ip=0;ip<NumOfPeaks;ip++){
				 		para_offset = Paras_offset_cal(ip,MainpeakIndex);
				 		fhandle_cp->SetParameter(para_offset+1, fhandle_cp->GetParameter(para_offset+1)+t_shift);
				 	}
				 }


				 int stableNcounts =0;
				 double* Px_record = new double[NumOfPeaks];

				 // fit 35 times for each random spectrum
				 for(int i=0;i<35;i++){	 //h_in_cp->Fit(fhandle_cp,(fitopt+"N").c_str(),"N", h_in->GetBinCenter(1),h_in->GetBinCenter(h_in->GetNbinsX())); 
					int paroffset=0;
				 	if(whichfunc == 'r'){
				 		fhandle_cp->SetParLimits(0,0,fhandle_cp->GetParameter(0)*5);
				 		double Px_now = fhandle_cp->GetParameter(1);
				 		fhandle_cp->SetParLimits(1, TMath::Max(_fitL,Px_now-sFWHM*1.5), TMath::Min(Px_now+sFWHM*1.5, _fitR) );
				 	}
					else{// func of 'x'

						// start setting limit for fit func
						if(NumOfPeaks>1){
									for(int ip=0;ip<NumOfPeaks;ip++){

										paroffset = Paras_offset_cal(ip,MainpeakIndex);
										fhandle_cp->SetParLimits(paroffset, 0, fhandle_cp->GetParameter(paroffset)*5); // Py limit
										double Px_ip = fhandle_cp->GetParameter(paroffset+1);

										if(ip==0){
											int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
											double Px_next = fhandle_cp->GetParameter(paroffset_nextPeak+1);
											fhandle_cp->SetParLimits(paroffset+1, TMath::Max(_fitL,Px_ip-sFWHM*1.5), TMath::Min(Px_ip+sFWHM*1.5, Px_next) );
										}else if(ip == NumOfPeaks-1){
											int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
											double Px_last = fhandle_cp->GetParameter(paroffset_lastPeak+1);
											fhandle_cp->SetParLimits(paroffset+1, TMath::Max(Px_last,Px_ip-sFWHM*1.5), TMath::Min(Px_ip+sFWHM*1.5, _fitR) );
										}else{
											int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
											double Px_last = fhandle_cp->GetParameter(paroffset_lastPeak+1);
											int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
											double Px_next = fhandle_cp->GetParameter(paroffset_nextPeak+1);
											fhandle_cp->SetParLimits(paroffset+1, TMath::Max(Px_last,Px_ip-sFWHM*1.5),  TMath::Min(Px_ip+sFWHM*1.5, Px_next));
										}
									}
						}else{// single peak case
								fhandle_cp->SetParLimits(0,0,fhandle_cp->GetParameter(0)*5);
								double temPx = fhandle_cp->GetParameter(1);
								fhandle_cp->SetParLimits(1, TMath::Max(_fitL, temPx-sFWHM*1.5), TMath::Min(temPx+sFWHM*1.5, _fitR)  );
						}

					}// end of setting limit for fit func


					if(i<20){
						if(i<3){h_in_cp->Fit(fhandle_cp,"EQN","N",_fitL,_fitR);}
						else if(i<5){h_in_cp->Fit(fhandle_cp,"LEQN","N",_fitL,_fitR);}
						else if(i%2==0){
								h_in_cp->Fit(fhandle_cp,"EQN","N",_fitL,_fitR); 
						}
						else{
							h_in_cp->Fit(fhandle_cp,"LEQN","N",_fitL,_fitR);
						}
					}
					else{
							
							h_in_cp->Fit(fhandle_cp,fitopt.c_str(),"",_fitL,_fitR);	

					}


					// update fitting range
					_fitL  = TMath::Max(fhandle_cp->GetParameter(1)-sample_range_L,h_in_cp->GetBinCenter(1));
					if(whichfunc == 'r'){  _fitR = TMath::Min( fhandle_cp->GetParameter(1)+sample_range_R, h_in_cp->GetBinCenter(h_in_cp->GetNbinsX()) );  }
					else{
						paroffset = Paras_offset_cal(NumOfPeaks-1,MainpeakIndex);
						_fitR = TMath::Min(fhandle_cp->GetParameter(paroffset+1)+sample_range_R,h_in_cp->GetBinCenter(h_in_cp->GetNbinsX()) );
					}


					bool SameResult = true;

					for(int ip=0;ip<NumOfPeaks;ip++){
						paroffset = Paras_offset_cal(ip,MainpeakIndex);
						if(Px_record[ip] !=  fhandle_cp->GetParameter(paroffset+1)){
							Px_record[ip] = fhandle_cp->GetParameter(paroffset+1);
							SameResult = false;
						}

						if(whichfunc == 'r'){ break; }
					}

					if(SameResult){
						stableNcounts++;						
					}
					else{
						stableNcounts=0;
					}


					if(stableNcounts == 6){ 
						//cout<<"\e[1;33mFit stable --> Stop\e[0m"<<endl;
						break;
					}


				}// end of 35 fitting for same histogram



/*
								fitopt.erase(remove(fitopt.begin(),fitopt.end(),'Q'),fitopt.end());
								fitopt.erase(remove(fitopt.begin(),fitopt.end(),'q'),fitopt.end());
								fitopt.erase(remove(fitopt.begin(),fitopt.end(),'N'),fitopt.end());
								fitopt.erase(remove(fitopt.begin(),fitopt.end(),'n'),fitopt.end());
								h_in_cp->Fit(fhandle_cp,fitopt.c_str(),"",_fitL,_fitR);		


::ctest2.Modified();
::ctest2.Update();
//cin>>temc;*/
				



				 bool ip_success[10];

				 for(int ip=0;ip<NumOfPeaks;ip++){
				 	int paroffset = Paras_offset_cal(ip,MainpeakIndex);
				 	double ip_result = fhandle_cp->GetParameter(paroffset+1);
				 	double ip_high = fhandle_cp->GetParameter(paroffset);

				 	if(ip_result > fhandle->GetParameter(paroffset+1)-2*sFWHM && ip_result < fhandle->GetParameter(paroffset+1)+2*sFWHM && ip_high>0){
				 		// fit successful
				 		ip_success[ip] = true;
				 	}
				 	else{
				 		ip_success[ip] = false;
				 		//cout<<"Peak "<<ip+1<<" fail"<<endl;
//cin>>temc;
				 	}

				 	if(whichfunc == 'r') break;

				 }

				 bool total_success = true;

				 for(int ip=0;ip<NumOfPeaks;ip++){total_success = total_success && ip_success[ip];  if(whichfunc == 'r') break;  }

				if(total_success){
					 for(int ip=0;ip<NumOfPeaks;ip++){
					 	int paroffset = Paras_offset_cal(ip,MainpeakIndex);
					 	fitcenter[ip].push_back( fhandle_cp->GetParameter(paroffset+1) );
					 	//cout<<fhandle_cp->GetParameter(paroffset+1)<<endl;
					 	if(whichfunc == 'r') break;
					 }
				}
				else{fhandle_cp->SetParameters(fhandle->GetParameters()); } // recover after bad fit
				 
				 h_in_cp->Reset();

				 delete[] Px_record;
//cin>>temc;
			}//end of loop of fitting random histogram



			cout<<"\e[1;32m"<<"MC simulate error (loop= "<<fitcenter[0].size()<<"):"<<"\e[0m"<<endl;
			for(int ip=0;ip<NumOfPeaks;ip++){
				printf("\033[1;32mPeak_%d: \033[5m\033[1;37m%.4f \033[0m\n",ip+1,TMath::RMS(fitcenter[ip].begin(),fitcenter[ip].end()));
				if(whichfunc == 'r') break;
			}

			cout<<"\e[1;37m MC simulate results are not be saved by default. Please update tof_ref_cento_err  or tof_x_cento_err[ipeak] manually!!!\e[0m"<<endl;

			fhandle->SetRange(temL,temR);

			delete h_in_cp;
			delete h_in_makeRandom;
			delete[] fitcenter;
		}



		void Fit_ref(TH1D* h_in, double _rangeL, double _rangeR, string fitopt="LMEQ"){

			if(EnableFunc_refInitialize) FitHelper(func_ref,h_in,_rangeL,_rangeR,Norder,fitopt);//return;
			else{
				for(int index=0;index< pars_per_peak;index++){
					func_ref->ReleaseParameter(index);
				}
			}

			for(int i=0;i<50;i++){
				h_in->Fit(func_ref,(fitopt+"N").c_str(),"",_rangeL,_rangeR);
				if(i<40){SetParLimitAuto(func_ref,Norder);}
				else{ 
					SetParLimitAuto(func_ref,Norder,0,false);
					for(int index=0;index< pars_per_peak;index++){
						func_ref->ReleaseParameter(index);
					}
				}
			}


			func_ref_cp->SetParameters(func_ref->GetParameters());

			LockPars_Ref();

			for(int i=0;i<5;i++){
				double Py_tem = func_ref->GetParameter(0);
				double Px_tem = func_ref->GetParameter(1);
				func_ref->SetParLimits(0, Py_tem*0.5,Py_tem*3);
				func_ref->SetParLimits(1, _rangeL,_rangeR  );
				h_in->Fit(func_ref,(fitopt+"N").c_str(),"",_rangeL,_rangeR);
			}

			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'Q'),fitopt.end());
			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'q'),fitopt.end());
			
			h_in->Fit(func_ref,fitopt.c_str(),"",_rangeL,_rangeR);

			sAmp = func_ref->GetParameter(0);
			sPeakCenter = func_ref->GetParameter(1);
			sPeakCenter_err = func_ref->GetParError(1);

			double tem_FWHM_L = func_ref->GetX(sAmp*0.5, sPeakCenter-1.5*sFWHM,sPeakCenter);
			double tem_FWHM_R = func_ref->GetX(sAmp*0.5, sPeakCenter,sPeakCenter+1.5*sFWHM);

			sFWHM = tem_FWHM_R-tem_FWHM_L;
			sSigma = UpdateSigma();

			cout<<"\e[1;33m"<<"Reference fitting result:"<<"\e[0m"<<endl;
			printf("Amp = %.2f \t Tof= %.4f(%.4f)\n",sAmp,sPeakCenter,sPeakCenter_err);
			printf("Rm = %f\n",sPeakCenter/(2*sFWHM));
			printf("FWHM = %f\n",sFWHM);

			func_ref_cp->SetParameters(func_ref->GetParameters());


			MC_ErrEval(h_in,'r',_rangeL,_rangeR,fitopt);

			sample_range_L = sPeakCenter - _rangeL; // get half range of fitting
			sample_range_R=  _rangeR - sPeakCenter;

			FreeRange = false;

		}

		void Fit_N(TH1D* h_in, double _rangeL, double _rangeR, string fitopt="LMEQ"){
			if(func_x==NULL){cout<<"function for fitting is not existed! MakeFitFunc_N() first!!!"<<endl; return;}
			if(LockPars){LockPars_X();}
			else{ Reles_X();}

			int paroffset=0;

			for(int ip=0;ip<NumOfPeaks;ip++){
				paroffset = Paras_offset_cal(ip,MainpeakIndex); // ip from 0; MainpeakIndex from 1
				func_x->SetParameter(paroffset, Amp[ip]);
				func_x->SetParameter(paroffset+1,tof_center[ip]);
			}

			double fit_L=0;
			double fit_R=0;
			int h_in_Nbins=h_in->GetNbinsX();

			int stableNcounts=0;  // if fitting results are the same for sequential 6 loops, stop fitting

			for(int i=0;i<50+10*(NumOfPeaks-1);i++){ 
				if(sample_range_L==-1 || sample_range_R==-1 || FreeRange){// new fitting
					if(i==0)cout<<"\e[1;33m"<<"free range"<<"\e[0m"<<endl;
					sample_range_L =tof_center[0]-_rangeL; // get half range of fitting
					sample_range_R= -(tof_center[NumOfPeaks-1]-_rangeR); // get right half range of fitting
				}
				
				fit_L = TMath::Max(tof_center[0]-sample_range_L,h_in->GetBinCenter(1));
				fit_R = TMath::Min(tof_center[NumOfPeaks-1]+sample_range_R,h_in->GetBinCenter(h_in_Nbins));

				if(i==0){
					if(fit_L==h_in->GetBinCenter(1)){
						cout<<"\e[1;2m"<<"Warning: fitting range is limited to the \"Left\" edge of histogram!!!!"<<"\e[0m"<<endl;
					}
					if(fit_R==h_in->GetBinCenter(h_in_Nbins)){
						cout<<"\e[1;2m"<<"Warning: fitting range is limited to the \"Right\" edge of histogram!!!!"<<"\e[0m"<<endl;
					}
				}

				if(!LockPars){// set limit for shape params of main peak and Amp of each peak
					for(int ip=0;ip<NumOfPeaks;ip++){
						if(ip == MainpeakIndex-1)	 SetParLimitAuto(func_x,Norder, Paras_offset_cal(MainpeakIndex-1,MainpeakIndex) ); // meain peak
						else{
								double countUpper=0;
								double countLower=0;
								paroffset = Paras_offset_cal(ip,MainpeakIndex);
								double tem_high = func_x->GetParameter(paroffset);
								countLower =TMath::Max( tem_high-TMath::Sqrt(tem_high)*2 , 0.);
								countUpper = tem_high+TMath::Sqrt(tem_high)*2;
								func_x->SetParLimits(paroffset,countLower,countUpper);
						}
					}
				}else{ // if shape parameters are locked, set the lower limit of Amp for each peak only

					for(int ip=0;ip<NumOfPeaks;ip++){
						paroffset = Paras_offset_cal(ip,MainpeakIndex);
						func_x->SetParLimits(paroffset,0,func_x->GetParameter(paroffset)*5);
					}
				}

				// set limit of position for each peak
				if(NumOfPeaks>1){
						for(int ip=0;ip<NumOfPeaks;ip++){
							paroffset = Paras_offset_cal(ip,MainpeakIndex);
							if(ip==0){
								int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
								func_x->SetParLimits(paroffset+1,fit_L,func_x->GetParameter(paroffset_nextPeak+1));
							}else if(ip == NumOfPeaks-1){
								int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
								func_x->SetParLimits(paroffset+1,func_x->GetParameter(paroffset_lastPeak+1),fit_R);
							}else{
								int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
								int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
								func_x->SetParLimits(paroffset+1, func_x->GetParameter(paroffset_lastPeak+1), func_x->GetParameter(paroffset_nextPeak+1) );
							}
						}
				}else{// single peak case
					double temPx = func_x->GetParameter(1);
					func_x->SetParLimits(1, TMath::Max(fit_L,temPx-sFWHM*1.5), TMath::Min(temPx+sFWHM*1.5,fit_R)  );
//					func_x->SetParLimits(0, func_x->GetParameter(0)*0.2, func_x->GetParameter(0)*5 );
					//func_x->SetParLimits(1, fit_L,fit_R);

				}


				if(i<20){
					if(i<2){h_in->Fit(func_x,"LEQN","N",fit_L,fit_R);}
					else if(i<5){h_in->Fit(func_x,"MEQN","N",fit_L,fit_R);}
					else if(i%2==0){
							h_in->Fit(func_x,"MEQN","N",fit_L,fit_R); 
					}
					else{
						h_in->Fit(func_x,"LEQN","N",fit_L,fit_R);
					}
				}
				else{
						h_in->Fit(func_x,(fitopt+"N").c_str(),"",fit_L,fit_R);						
				}

/*
double temPx = func_x->GetParameter(1);
double limitL, limitR;
func_x->GetParLimits(1,limitL, limitR);
printf("%.4f\t%.4f\n",temPx-limitL,limitR-temPx);

char ctest;
cin>>ctest;*/
			
				bool PeakChanged = false;
					for(int ip=0;ip<NumOfPeaks;ip++){
						paroffset = Paras_offset_cal(ip,MainpeakIndex);
						if(tof_center[ip] != func_x->GetParameter(paroffset+1) ){ tof_center[ip] = func_x->GetParameter(paroffset+1); PeakChanged = true;}
					}

				if(!PeakChanged){
					stableNcounts++;
					if(stableNcounts>=6){ cout<<"\e[1;33m Fitting reach stable --> stop\e[0m"<<endl; break;  }
				}else{
					stableNcounts =0;
				}


			}//end of for 50 fitting


			if(!LockPars){
				for(int i=0;i<20;i++){// free high of each peak to fit again

					for(int ip=0;ip<NumOfPeaks;ip++){
						if(ip == MainpeakIndex-1)	 SetParLimitAuto(func_x,Norder, Paras_offset_cal(MainpeakIndex-1,MainpeakIndex), false ); // meain peak
						else{
								paroffset = Paras_offset_cal(ip,MainpeakIndex);
								//func_x->SetParLimits(paroffset,0.,0.); //// completely free
								func_x->SetParLimits(paroffset,0.,func_x->GetParameter(paroffset)*10); 								
						}
					}

					// set limit of position for each peak
					if(NumOfPeaks>1){
							for(int ip=0;ip<NumOfPeaks;ip++){
								paroffset = Paras_offset_cal(ip,MainpeakIndex);
								if(ip==0){
									int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
									func_x->SetParLimits(paroffset+1,fit_L,func_x->GetParameter(paroffset_nextPeak+1));
								}else if(ip == NumOfPeaks-1){
									int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
									func_x->SetParLimits(paroffset+1,func_x->GetParameter(paroffset_lastPeak+1),fit_R);
								}else{
									int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
									int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
									func_x->SetParLimits(paroffset+1, func_x->GetParameter(paroffset_lastPeak+1), func_x->GetParameter(paroffset_nextPeak+1) );
								}
							}
					}else{
						func_x->SetParLimits(1,fit_L,fit_R);
					}

					h_in->Fit(func_x,(fitopt+"N").c_str(),"",fit_L,fit_R);
				}
			}


			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'Q'),fitopt.end());
			fitopt.erase(remove(fitopt.begin(),fitopt.end(),'q'),fitopt.end());

			// set limit of position for each peak
			if(NumOfPeaks>1){
					for(int ip=0;ip<NumOfPeaks;ip++){
						paroffset = Paras_offset_cal(ip,MainpeakIndex);
						if(ip==0){
							int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
							func_x->SetParLimits(paroffset+1,fit_L,func_x->GetParameter(paroffset_nextPeak+1));
						}else if(ip == NumOfPeaks-1){
							int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
							func_x->SetParLimits(paroffset+1,func_x->GetParameter(paroffset_lastPeak+1),fit_R);
						}else{
							int paroffset_lastPeak = Paras_offset_cal(ip-1,MainpeakIndex);
							int paroffset_nextPeak = Paras_offset_cal(ip+1,MainpeakIndex);
							func_x->SetParLimits(paroffset+1, func_x->GetParameter(paroffset_lastPeak+1), func_x->GetParameter(paroffset_nextPeak+1) );
						}
					}
			}


			// set limit of Amp for each peak
			for(int ip=0;ip<NumOfPeaks;ip++){
					paroffset = Paras_offset_cal(ip,MainpeakIndex);
					func_x->SetParLimits(paroffset,0,func_x->GetParameter(paroffset)*10);
			}


			h_in->Fit(func_x,fitopt.c_str(),"",fit_L,fit_R);


			func_x_cp->SetParameters(func_x->GetParameters());

			if(useMC)MC_ErrEval(h_in,'x',fit_L,fit_R,fitopt);

			printf("Main Peak Index = %d ; Norder = %d\n",MainpeakIndex,Norder);

			/*for(int ip=0; ip<NumOfPeaks;ip++){
				paroffset = Paras_offset_cal(ip,MainpeakIndex); // ip from 0; MainpeakIndex from 1
				Amp[ip] = func_x->GetParameter(paroffset);
				tof_center[ip] = func_x->GetParameter(paroffset+1);
				tof_center_err[ip] = func_x->GetParError(paroffset+1);
				printf("cento_%d:  %.4f (%.4f)\n",ip+1,tof_center[ip],tof_center_err[ip]);
			}*/

			UpdateParsFromTF();

			sample_range_L = tof_center[0]-fit_L;
			sample_range_R = -(tof_center[NumOfPeaks-1]-fit_R);


			if(!LockPars)cout<<"shape parameters are free"<<endl;

			FreeRange = false;

		}


		bool Draw_subline(TCanvas* c_todraw=NULL,int Padindex=4){
			if(Padindex>4 || Padindex<1){cout<<"Peak index at canvas is wrong [1,4]"<<endl; return false;}

			if(c_todraw!=NULL){
				c_todraw->cd(Padindex)->SetEditable(kTRUE);
			}
			else{cout<<"canvas is not available!!!!"<<endl; return false;}

			for(int i=0;i<10;i++){ if(fresult[i]!=NULL){ delete fresult[i]; fresult[i]=NULL;} }

			for(int i=0;i<NumOfPeaks;i++){
				int kala[]={633,808,799,417,433,600,617};				
				fresult[i] = new TF1(Form("fext_result_%d",i+1),this,&funcN::fitfunc,0,25e6,pars_per_peak,"1func_sub","1fitfunc_sub");
				fresult[i]->SetParameter(0,Amp[i]);
				fresult[i]->SetParameter(1,tof_center[i]);
				int paroffest = Paras_offset_cal(MainpeakIndex-1,MainpeakIndex);
				for(int j=2;j<pars_per_peak;j++){
					fresult[i]->SetParameter(j,func_x->GetParameter(paroffest+j)); 
				}
				fresult[i]->SetLineColor(kala[i%7]);
				fresult[i]->SetLineStyle(2);
				fresult[i]->SetNpx(500);
				fresult[i]->Draw("same");
			}

			c_todraw->cd(Padindex)->Modified();
			c_todraw->cd(Padindex)->Update();
			if(Padindex==2)c_todraw->cd(Padindex)->SetEditable(kFALSE);
			return true;
		}


		funcN(){
			sAmp=0;   // Amp of sampled peak
			sPeakCenter=0; // peak center of sampled peak;
			sPeakCenter_err=0;
			sFWHM=0;
			sSigma=10;
			for(int i=0;i<10;i++){
				Amp[i]=0;
				tof_center[i]=0;
				tof_center_err[i]=0;
				fresult[i]=NULL;
			}
			sample_range_L=-1; // half range for fit
			sample_range_R=-1;
			func_ref=NULL;
			func_x=NULL;
			func_ref_cp=NULL;
			func_x_cp=NULL;
			MC_sim_counts=100;


			Norder=1;
			NumOfPeaks=1;
			MainpeakIndex=1;
			pars_per_peak = 5+Norder*2;
			par_ref=NULL;

			LockPars=true; //lock shape parameter
			FreeRange =true;
			AutoUpdateMainPeakIndex=true;
			useMC=false;

			EnableFunc_refInitialize=true;
		}

		~funcN(){
			if(func_ref!=NULL){delete func_ref; func_ref = NULL; delete func_ref_cp; func_ref_cp = NULL;}
			if(func_x !=NULL){delete func_x; func_x = NULL; delete func_x_cp; func_x_cp=NULL;}
			if(par_ref != NULL){delete[] par_ref; par_ref=NULL;}
		}


};


int funcN::NPs2Set=1;
int funcN::Norder2Set=0;


#endif