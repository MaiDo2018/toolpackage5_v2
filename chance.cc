#include<iostream>
#include<stdlib.h>
#include<vector>
#include"TMath.h"
#include"TRandom3.h"
#include<ctime>
#include"TF1.h"
#include"TH1D.h"


#include "Math/WrappedMultiTF1.h"
#include "Fit/DataRange.h"
#include "Fit/UnBinData.h"
#include "Fit/Fitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

using namespace std;


void FakeCoincidence(int numtof=3,int numbeta=2,double Time_measure = 20,double coingate=1, int Maxloop=5){
	
		TRandom3* ranTOF = new TRandom3((unsigned int)time(0));
		TRandom3* ranBeta= new TRandom3(5192+(unsigned int)time(0)/209);

		/*	double Time_measure = 20; // second;
			double coingate= 1;  // second;

			int numtof=3;  // num of tof events
			int numbeta=2;  // num of beta events*/

		double* tof_time = new double[numtof];
		int* sequence_tof_list = new int[numtof];
		double* tof_time_s = new double[numtof];

		double* beta_time = new double[numbeta];
		int* sequence_beta_list = new int[numbeta];
		double* beta_time_s  = new double[numbeta];

		vector<double> tof_list;
		vector<double> beta_list;

		int numcase = TMath::Min(numtof,numbeta)+1;

		int* result_counts = new int[numcase];

		for(int index = 0; index<numcase;index++){
			result_counts[index]=0;
		}


		for(int nloop=0;nloop< Maxloop; nloop++){
				tof_list.clear();
				beta_list.clear();

				if(nloop%1000==0) cout<<"processing......"<<(double)nloop/Maxloop*100<<"%"<<"\r"<<flush;

				for(int index_tof=0; index_tof<numtof; index_tof++){

					tof_time[index_tof] =ranTOF->Uniform(Time_measure);

				}

				for(int index_beta=0;index_beta<numbeta;index_beta++){
					beta_time[index_beta] = ranBeta->Uniform(Time_measure);

				}


				TMath::Sort(numtof,tof_time,sequence_tof_list,kFALSE);
				TMath::Sort(numbeta,beta_time,sequence_beta_list,kFALSE);


				for(int index_tof=0; index_tof<numtof; index_tof++){

					tof_time_s[index_tof] = tof_time[sequence_tof_list[index_tof]];
					//cout<<tof_time_s[index_tof]<<"\t";

				}
				//cout<<endl;

				for(int index_beta=0; index_beta<numbeta; index_beta++){

					beta_time_s[index_beta] = beta_time[sequence_beta_list[index_beta]];
					//cout<<beta_time_s[index_beta]<<"\t";

				}
				//cout<<endl;


					for(int index_beta=0;index_beta<numbeta;index_beta++){

						for(int index_tof=0;index_tof<numtof;index_tof++){
							if(tof_time_s[index_tof]<beta_time_s[index_beta] && index_tof<numtof-1) continue;
							else if(tof_time_s[index_tof]>=beta_time_s[index_beta]){
										if( index_tof>=1){index_tof--;}
										else{break;}
							}
														
								
							
							if(beta_time_s[index_beta]-tof_time_s[index_tof] < coingate){ //candidate
								bool newtof = true, newbeta =true;
								for(int ilist=0;ilist<(int)tof_list.size();ilist++){
									if(tof_time_s[index_tof] == tof_list[ilist]){ newtof=false;break;} // repeat tof event
								}
								for(int jlist=0;jlist<(int)beta_list.size();jlist++){
									if(beta_time_s[index_beta] == beta_list[jlist]){newbeta=false;break;}
								}
								if(newtof && newbeta){
									tof_list.push_back(tof_time_s[index_tof]);
									beta_list.push_back(beta_time_s[index_beta]);
								}

							}
								
								break; // to next beta
							
						}// tof loop

					}// beta loop


					result_counts[tof_list.size()]++;
		}// end of nloop

			cout<<endl;
			cout<<"Num of coincidence \t Happens in loop \t Probability"<<endl;
			double mean=0;
			for(int index=0;index<numcase;index++){
				printf("%d \t %d \t %.2f%%\n",index,result_counts[index],(double)result_counts[index]/Maxloop*100);
				mean+=index*(double)result_counts[index]/Maxloop;
			}
			cout<<"mean = "<<mean<<endl;

	delete ranTOF;
	delete ranBeta;
	delete[] tof_time;
	delete[] sequence_tof_list;
	delete[] tof_time_s;
	delete[] beta_time;
	delete[] sequence_beta_list;
	delete[] beta_time_s;
	delete[] result_counts;




}




void ToyBg(int numtof=3,int numbeta=2,double Time_measure = 20, double coin_window=10, int Maxloop=5){ //coin_window==> window for coincidence in second
	
		TRandom3* ranTOF = new TRandom3((unsigned int)time(0));
		TRandom3* ranBeta= new TRandom3(7192+(unsigned int)time(0)/209);

		/*	double Time_measure = 20; // second;
			double coingate= 1;  // second;

			int numtof=3;  // num of tof events
			int numbeta=2;  // num of beta events*/

		double* tof_time = new double[numtof];
		int* sequence_tof_list = new int[numtof];
		double* tof_time_s = new double[numtof];

		double* beta_time = new double[numbeta];
		int* sequence_beta_list = new int[numbeta];
		double* beta_time_s  = new double[numbeta];

		vector<double> tof_list;
		vector<double> beta_list;

		int numcase = TMath::Min(numtof,numbeta)+1;

		vector<double> mean_store;// negative side
		vector<double> mean_store_posi; // positive side

		vector<double>store_rate; // = Ncounts_in_window / (Ntof*window) 


		for(int nloop=0;nloop< Maxloop; nloop++){
				tof_list.clear();
				beta_list.clear();

				if(nloop%1000==0) cout<<"processing......"<<(double)nloop/Maxloop*100<<"%"<<"\r"<<flush;

				for(int index_tof=0; index_tof<numtof; index_tof++){

					tof_time[index_tof] =ranTOF->Uniform(Time_measure);

				}

				for(int index_beta=0;index_beta<numbeta;index_beta++){
					beta_time[index_beta] = ranBeta->Uniform(Time_measure);

				}


				TMath::Sort(numtof,tof_time,sequence_tof_list,kFALSE);
				TMath::Sort(numbeta,beta_time,sequence_beta_list,kFALSE);


				for(int index_tof=0; index_tof<numtof; index_tof++){

					tof_time_s[index_tof] = tof_time[sequence_tof_list[index_tof]];
					//cout<<tof_time_s[index_tof]<<"\t";

				}
				//cout<<endl;

				for(int index_beta=0; index_beta<numbeta; index_beta++){

					beta_time_s[index_beta] = beta_time[sequence_beta_list[index_beta]];
					//cout<<beta_time_s[index_beta]<<"\t";

				}
				//cout<<endl;

				double time_diff=0;
				int index_beta_near=-1;
				vector<double>store_diff;
				vector<double>store_diff_posi; // demonstrate the measurement by average beta interval is symmetric

				int Ncounts_in_gate=0;


				for(int index_tof=0;index_tof<numtof;index_tof++){
						index_beta_near=-1;

						
					for(int index_beta=0;index_beta<numbeta;index_beta++){

						if(beta_time_s[index_beta] < tof_time_s[index_tof])index_beta_near = index_beta;
						else if(beta_time_s[index_beta] >=tof_time_s[index_tof] && beta_time_s[index_beta] <=tof_time_s[index_tof]+coin_window){
							Ncounts_in_gate++;
						}
						else break;

					}// beta loop

					if(index_beta_near !=-1){
						store_diff.push_back(tof_time_s[index_tof]-beta_time_s[index_beta_near]);
						if(index_beta_near+1<numbeta) store_diff_posi.push_back(beta_time_s[index_beta_near+1]-tof_time_s[index_tof]);

					}

				}

				double mean_diff=0;

				for(int index=0;index<(int)store_diff.size();index++){ // negative side

					mean_diff+=store_diff[index];

				}

				mean_diff/=(int)store_diff.size();

				mean_store.push_back(mean_diff);

				mean_diff=0;

				for(int index=0;index<(int)store_diff_posi.size();index++){
						mean_diff+=store_diff_posi[index];
				}

				mean_diff/=(int)store_diff_posi.size();

				mean_store_posi.push_back(mean_diff);

				store_rate.push_back(Ncounts_in_gate/(coin_window*numtof));
					
		}// end of nloop


			double mean=0;
			for(int index=0;index<(int)mean_store.size();index++){
				mean+=mean_store[index];
			}
mean/=(int)mean_store.size();

			double mean2=0;// by constant window
			for(int index=0;index<(int)store_rate.size();index++){
				mean2+=store_rate[index];
			}
mean2/=(int)store_rate.size();

			double mean_posi=0;
			for(int index=0;index<(int)mean_store_posi.size();index++){
				mean_posi+=mean_store_posi[index];
			}
mean_posi/=(int)mean_store_posi.size();

			cout<<"mean bg interval_negative= "<<mean<<" ; by positive interval= "<<mean_posi<<endl;
			cout<<"==> mean bg rate_negative = "<<1/mean<<" ; rate_positive = "<<1/mean_posi<<endl;
			cout<<"mean bg rate by const window "<<coin_window<<" s = "<<mean2<<endl;

	delete ranTOF;
	delete ranBeta;
	delete[] tof_time;
	delete[] sequence_tof_list;
	delete[] tof_time_s;
	delete[] beta_time;
	delete[] sequence_beta_list;
	delete[] beta_time_s;
	//delete[] result_counts;




}





double ExpectHalflife(double tao=1, double tUpperLimit=10, int nloop=100000){
		double sum=0;
		int N_counts=0;
		TRandom3 * r = new TRandom3((unsigned int)time(0));
		for(int i=0; i<nloop;i++){
			double getvalue = r->Exp(tao);
			if(getvalue<=tUpperLimit){
				sum+=getvalue;
				N_counts++;
			}			
		}

		return (sum/N_counts) * TMath::Log(2);
		
}





// find the zero point of a function (double x) , with single variable to solve.
double GetZeroX (double(*fcn)(double), double LeftX, double RightX, double pricision=0.01){ // relative_err=1%

	if(LeftX>RightX){
		double tem_x = LeftX;
		LeftX = RightX;
		RightX = tem_x;
	}
	
	double Y_Left= fcn(LeftX);
	double Y_Right= fcn(RightX);

	if(Y_Left*Y_Right>0){// no zero point

		cout<<"no zero point found in the region ["<<LeftX<<" , "<<RightX<<"], return LeftX-1"<<endl;
		return LeftX-1; // return a value out of the region
		
	}
	if(Y_Left==0) return LeftX;
	if(Y_Right==0) return RightX;

	double Solution_last=0;
	double X_tem,Y_tem;

	do{

		X_tem=(LeftX+RightX)/2;
		Y_tem = fcn(X_tem);
		if(Y_tem==0) return X_tem;
		else if(Y_Left * Y_tem<0){
			RightX = X_tem;
			Y_Right = Y_tem;
		}
		else{
			LeftX = X_tem;
			Y_Left = Y_tem;
		}
//cout<<X_tem<<","<<LeftX<<","<<RightX<<endl;
	}while((RightX-LeftX)/X_tem>pricision);

	return X_tem;

}



double GetZeroX (double(*fcn)(double*,double*), double LeftX, double RightX, double* par,double pricision=0.01){ // relative_err=1%

	if(LeftX>RightX){
		double tem_x = LeftX;
		LeftX = RightX;
		RightX = tem_x;
	}
	
	double tao[1];

	tao[0]=LeftX;
	double Y_Left= fcn(tao,par);
	tao[0]=RightX;
	double Y_Right= fcn(tao, par);

	if(Y_Left*Y_Right>0){// no zero point

		cout<<"no zero point found in the region ["<<LeftX<<" , "<<RightX<<"], return LeftX-1"<<endl;
		return LeftX-1; // return a value out of the region
		
	}
	if(Y_Left==0) return LeftX;
	if(Y_Right==0) return RightX;

	
	double X_tem,Y_tem;

	do{

		X_tem=(LeftX+RightX)/2;
		tao[0] = X_tem;
		Y_tem = fcn(tao,par);
		if(Y_tem==0) return X_tem;
		else if(Y_Left * Y_tem<0){
			RightX = X_tem;
			Y_Right = Y_tem;
		}
		else{
			LeftX = X_tem;
			Y_Left = Y_tem;
		}
//cout<<X_tem<<","<<LeftX<<","<<RightX<<endl;
	}while((RightX-LeftX)/X_tem>pricision);

	return X_tem;

}

/*
double Beta_time_AVG=0;
double DK_gate_window=1; // window to get coincidence unit in [s]
double effici = 0.33; // detector efficiency

double Fcn_no_bg(double tao){// lifetime correction function without considering background
	double lamda_T = DK_gate_window/tao;
	return (1-TMath::Exp(-lamda_T)*lamda_T / (1-TMath::Exp(-lamda_T)))*tao-Beta_time_AVG;
}


double P_density_bg(double * t,double* par){ // probability of the first beta signal at time t with background rate b counts /[s]
// par[0] = lamda; par[1]= ratio of b/lamda; par[2]= coincident window in [s]

	double C_normalized = effici*(1-TMath::Exp(-(1+par[1])*par[0]*par[2])) + (1-effici) * (1- TMath::Exp(-par[1]*par[0]*par[2]));

	return ( effici*(1+par[1])*par[0] *TMath::Exp(-(1+par[1])*par[0]* t[0]) + (1-effici) *par[1]*par[0]*TMath::Exp(-par[1]*par[0]*t[0]) )/C_normalized;



}


double P_density_no_bg(double * t,double* par){//par[0]= lamda; par[1] = time window in [s]

		double C_normalized = 1- TMath::Exp(-par[0]*par[1]);

		return par[0]* TMath::Exp(-par[0]*t[0]) / C_normalized;

}


double P_for_fit(double* t, double* par){// par[0] = TAO !!!!; par[1]= ratio of b/lamda; par[2]= coincident window in [s]; par[3] = effici

	if(par[3]!=-1) effici = par[3];

	double _par[3];
	_par[0] = 1/par[0]; // tao => lamda
	_par[1]= par[1];
	_par[2]= par[2];

	return P_density_bg(t,_par);



}


double Expect_DK_time_bg(double *tao, double* par){// par[0] = ration b/lamda; par[1]  = time window in [s]; par[2] = effici;

	double term1 = par[2]*(1-TMath::Exp(-(1+par[0])*par[1]/tao[0])) + (1-par[2]) * (1- TMath::Exp(-par[0]*par[1]/tao[0]));
	double term2 = par[2]/(1+par[0])*tao[0] * (1-TMath::Exp(-(1+par[0])*par[1]/tao[0]) - (1+par[0])*par[1]/tao[0] * TMath::Exp(-(1+par[0])*par[1]/tao[0]));
	double term3 = (1-par[2])/par[0]*tao[0] *( 1- TMath::Exp(-par[0]*par[1]/tao[0])- par[0]*par[1]/tao[0]* TMath::Exp(-par[0]*par[1]/tao[0]) );

	cout<<term1<<" , "<<term2<<" , "<<term3<<endl;

	return (term2+term3)/term1;

}


double Expect_DK_time_no_bg(double *tao, double* par){// par[0]= time window in [s]
	double term1 = 1 - TMath::Exp(-par[0]/tao[0]);
	double term2 = (1 - TMath::Exp(-par[0]/tao[0]) - par[0]/tao[0]*TMath::Exp(- par[0]/tao[0]))* tao[0];

	cout<<term1<<" , "<<term2<<endl;

	return term2/term1;
}




double Fcn_w_bg(double* tao, double* par){
	return Expect_DK_time_bg(tao,par) - Beta_time_AVG;
}


TH1D* h = new TH1D("h","h",100,0,10);*/

vector<double> T_data; // all beta decay time (positive and negative) in unite of [s]
vector<double> T_data_neg; // only negative decay time. (the closest negative decay to TOF or all negative decay in [times2halflife * halflife, 0 ]) in unite of [s]
string MinmizerName= "Minuit";
string AlgorithmName = "";

TH1D* h_beta_accumulate = NULL;
double binwidth_fraction_halflife=10; // binwidth = halflife/10;


/*double LikelihoodFunction(double * tao, double* par){ // variable = tao lifetime; par[0] = ratio b/lamda, par[1]=coincidence window; par[2]=effici efficiency

	double t[1];
	double Log_P_t=0; // possibilit at t

	for(unsigned long i=0;i<T_data.size();i++){
			t[0] = T_data[i];
			double _par[4];
			_par[0] = tao[0];
			_par[1] = par[0];
			_par[2] = par[1];
			_par[3] = par[2];
			Log_P_t+= TMath::Log(P_for_fit(t,_par));
	}

	return Log_P_t;

}


TF1* flog=NULL;


void try_fit(int N_simulate=500){
h->Reset();
	T_data.clear();
Beta_time_AVG=0;


	TF1* f1 = new TF1("f1",P_density_bg,0,100,3);

	f1->SetParameters(1.25,0.001,50);

	const int N_dk_counts=N_simulate;

	double t_data[N_dk_counts];

	for(int i=0;i<N_dk_counts;i++){
		t_data[i] = f1->GetRandom(0,50);
		T_data.push_back(t_data[i]);
		h->Fill(t_data[i]);
		Beta_time_AVG+=t_data[i];
	}

	Beta_time_AVG /= N_dk_counts;

	double Par_set[3]={0.1,50,0.33};

	double tao_result = GetZeroX(Fcn_w_bg,0.001,10,Par_set);

	cout<<"\e[1;33m"<<"tao_result = "<<tao_result<<"\e[0m"<<endl;
	cout<<endl;


	if(flog!=NULL) delete flog;
	flog = new TF1("flog",LikelihoodFunction,0,1000,3);

	flog->SetParameters(0.001,50,-1);
	flog->Draw();

double gettao = flog->GetMaximumX(0.001,1000);
cout<<"tao = "<<gettao<<endl;
cout<<"tao up="<<flog->GetX(flog->Eval(gettao)-0.5,gettao,1000)<<endl;
cout<<"tao low="<<flog->GetX(flog->Eval(gettao)-0.5,0.001,gettao)<<endl;



	ROOT::Fit::DataRange range(0,50);
    ROOT::Fit::UnBinData data((unsigned int)N_dk_counts, t_data, range);

    TF1* f2 = new TF1("f2",P_for_fit,0,100,4);
    f2->SetParameters(4,0.001,50,-1);

    ROOT::Math::WrappedMultiTF1 fitFunction( *f2, f2->GetNdim());

    ROOT::Fit::Fitter fitter;
    fitter.SetFunction( fitFunction, false);// false ==> not use customer defined gradient

    int Npars = f2->GetNpar();
    double* Pars = f2->GetParameters();
    fitter.Config().SetParamsSettings(Npars,Pars);
	fitter.Config().ParSettings(0).SetLimits(0,100);
    fitter.Config().ParSettings(1).Fix();
    fitter.Config().ParSettings(2).Fix();
    fitter.Config().ParSettings(3).Fix();


  //  fitter.Config().SetMinimizer("Minuit2","Migrad");
	fitter.Config().SetUpdateAfterFit();
	fitter.LikelihoodFit(data);

    TFitResult r=fitter.Result();
    r.Print();

    cout<<r.Parameter(0)<<endl;
    cout<<r.LowerError(0)<<endl;
    cout<<r.UpperError(0)<<endl;

    delete f2;
    delete f1;




}
*/

double Ratio_true(double* t, double* par){// par[0]=lamda; par[1]=ratio b/lamda; par[2]=efficiency
	double P_true = par[2]*par[0]*TMath::Exp(-(1+par[1])*par[0]*t[0]);
	double P_eve_1st = par[2] *(1+par[1])* par[0]* TMath::Exp(-(1+par[1])* par[0]*t[0]) + (1-par[2])*par[1]*par[0]*TMath::Exp(-par[1]*par[0]*t[0]);
	return P_true/P_eve_1st;

}






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rewrite functioins with "independent background rate" %%%%%%%%%%%%%%%%%%%%%%%%%%%%
// above function use the ratio background to lamda to describe background
double P_density_dk(double* t, double* par){// par0= lamda; par1= background /[s]; par2= time window in [s]; par3= efficienncy

	double lamda = par[0];
	double vb = par[1];
	double T_window = par[2];
	double effi = par[3];
	double C_normalized=0;
	double result=0;

	if(vb==0){// no background
		C_normalized = 1- TMath::Exp(-lamda*T_window);
		result = lamda* TMath::Exp(-lamda*t[0]) / C_normalized;

	}
	else{
		C_normalized = effi*(1-TMath::Exp(-(vb+lamda)*T_window)) + (1-effi) * (1- TMath::Exp(-vb*T_window));
		result = ( effi*(vb+lamda)* TMath::Exp(-(vb+lamda)*t[0]) + (1-effi)*vb*TMath::Exp(-vb*t[0])) / C_normalized;
	}

	return result;

}


double Expect_DK_time(double * tao, double* par){// par0=background rate /s ; par1= time window in [s]; par2 = efficiency;

	double lamda = 1/ tao[0];
	double vb = par[0];
	double T_window = par[1];
	double effi = par[2];
	double result=0;

	if(vb==0){// no background

		result = (1- lamda*T_window*TMath::Exp(-lamda*T_window)/(1- TMath::Exp(-lamda*T_window)))/lamda;

	}
	else{

		double term1,term2,term3;
		term1 = effi*(1-TMath::Exp(-(vb+lamda)*T_window)) + (1-effi) * (1- TMath::Exp(-vb*T_window));
		term2 = (1-TMath::Exp(-(vb+lamda)*T_window) -(vb+lamda)*T_window * TMath::Exp(-(vb+lamda)*T_window)) * effi / (vb+lamda);
		term3 = ( 1-TMath::Exp(-vb*T_window) -vb*T_window*TMath::Exp(-vb*T_window) ) *(1-effi)/vb;
		result = (term2+term3)/term1;

	}

	return result;


}


double ExtractTaobyDKtime(double* tao, double*par){// par0=background rate /s ; par1= time window in [s]; par2 = efficiency; par3 = decay Time_avg

	return Expect_DK_time(tao,par)- par[3]; // To get the zero point, which is the solution of tao

}


double P_density_dk_tao(double* t, double* par){// par[0] = TAO !!!!; par[1]= backgrounds / s; par[2]= coincident window in [s]; par[3] = effici
	// for unbinned fitting as Pdf

	double _par[4];
	_par[0] = 1/par[0]; // tao => lamda
	_par[1]= par[1];
	_par[2]= par[2];
	_par[3]= par[3];

	return P_density_dk(t,_par);

}


double LikelyhoodFcn(double* tao,double* par){// likelyhood vs tao ==> maximized
	// par0= backgrounds in [s]; par1 = coincident window in [s]; par2 = efficiency
	//par3=Ntof;

	double t[1];
	double Log_P_t=0; // possibilit at t
	double ntof = par[3];

	double vb = par[0];
	double T_window = par[1];
	double effi = par[2];
	double lamda = 1/tao[0];

	for(unsigned long i=0;i<T_data.size();i++){
			t[0] = T_data[i];
			double _par[4];
			_par[0] = tao[0];
			_par[1] = par[0];
			_par[2] = par[1];
			_par[3] = par[2];
			Log_P_t+= TMath::Log(P_density_dk_tao(t,_par));
	}

	double n_beta_expect = ntof* par[2]*(1- TMath::Exp(-par[1]/tao[0]));
	double n_bg = ntof*par[0]*par[1];

    double mean_n = ntof * (effi*(1-TMath::Exp(-(vb+lamda)*T_window)) + (1-effi) * (1- TMath::Exp(-vb*T_window)) );


	//double mean_n = n_beta_expect + n_bg;

	int ncounts = (int)T_data.size();

	Log_P_t += TMath::Log(TMath::PoissonI(ncounts,mean_n));//TMath::Exp(-mean_n)
	//Log_P_t += -mean_n;


	return -Log_P_t;

}

bool fit_vb=false;
bool fit_efficiency=false;
bool fit_tao=true;

double Multi_LikelyhoodFcn(const double* inpar){// for use in minimizer to fit tao and efficiency
// inpar0 = tao ; inpar1= backgrounds in [s]; 
	//inpar2 = coincident window in [s]; inpar3 = efficiency; inpar4=Ntof;

	double* _inpar = (double *)inpar;
	return LikelyhoodFcn(_inpar,&_inpar[1]);
}



int NumericalMinimization(const double* inpar, const char * minName = "Minuit2",
                           const char *algoName = "" ,
                           int randomSeed = -1)
{
    // create minimizer giving a name and a name (optionally) for the specific
    // algorithm
    // possible choices are:
    //     minName                  algoName
    // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
    //  Minuit2                     Fumili2
    //  Fumili
    //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
    //                              BFGS2, SteepestDescent
    //  GSLMultiFit
    //   GSLSimAn
    //   Genetic
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    if(minimum==NULL){
    	cout<<"\e[1;31m"<<"Fail to define minimizer, try another one!!"<<"\e[0m"<<endl;
    	return -1;
    }

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.0001);
  //  minimum->SetPrintLevel(1);
    minimum->SetErrorDef(0.5); // 0.5 for likely hood; 1 for chi2
    minimum->SetStrategy(2);//2==> more reliabe
    minimum->SetValidError(true); // accurate error analysis

    // create function wrapper for minimizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&Multi_LikelyhoodFcn,5);
    double step[5] = {0.00001,0.00001,0.01,0.00001,0.01};
    // starting point

    double variable[5] = {inpar[0],inpar[1],inpar[2],inpar[3],inpar[4]};// inpar0 = tao ; inpar1= backgrounds in [s]; 
	//inpar2 = coincident window in [s]; inpar3 = efficiency; inpar4=Ntof;
    if (randomSeed >= 0) {
       TRandom2 r(randomSeed);
       variable[0] = r.Uniform(0.001,20);
       variable[3] = r.Uniform(0.001,0.34);
    }

    minimum->SetFunction(f);

    // Set the free variables to be minimized !
    minimum->SetVariable(0,"tao",variable[0], step[0]);
    minimum->SetVariable(1,"background/[s]",variable[1], step[1]);
    minimum->SetVariable(2,"T_window [s]",variable[2], step[2]);
    minimum->SetVariable(3,"efficiency",variable[3], step[3]);
    minimum->SetVariable(4,"Ntof",variable[4], step[4]);

	if(fit_tao){
		//minimum->ReleaseVarible(0);
		minimum->SetLimitedVariable(0,"tao",variable[0], step[0],0.001,50);
	}
	else{
		minimum->FixVariable(0);
	}

	if(fit_vb){
		//minimum->ReleaseVarible(1);
		minimum->SetLimitedVariable(1,"background/[s]",variable[1], step[1],1e-5,5);
	}
	else{
		minimum->FixVariable(1);
	}

	if(fit_efficiency){
		//minimum->ReleaseVarible(3);
		minimum->SetLimitedVariable(3,"efficiency",variable[3], step[3],0.001,0.34);
	}
	else{
		minimum->FixVariable(3);
	}

	minimum->FixVariable(2);
	minimum->FixVariable(4);


    // do the minimization
    minimum->Minimize();

    const double *xs = minimum->X();
    const double *xs_err = minimum->Errors();

    double errL,errU;

    cout<<endl;
	cout<<endl;
    cout<<"*******************************************************"<<endl;
    if(fit_tao){
    	minimum->GetMinosError(0,errL,errU);
    	errL*=TMath::Log(2);
    	errU*=TMath::Log(2);
    	printf("T1/2 [s]= %.4f (%.4f) or (%.4f,+%.4f)\n",xs[0]*TMath::Log(2),xs_err[0]*TMath::Log(2),errL,errU);
    }
    else{
    	printf("T1/2 [s]= %.4f fixed\n",xs[0]*TMath::Log(2));
    }

    if(fit_vb){
    	minimum->GetMinosError(1,errL,errU);
 		printf("vb background /s = %.4f (%.4f) or (%.4f,+%.4f)\n",xs[1],xs_err[1],errL,errU);
    }
    else{
    	printf("vb background /s =%.4f fixed\n",xs[1]);
    }
   

   	if(fit_efficiency){
		minimum->GetMinosError(3,errL,errU);
		printf("efficiency = %.4f (%.4f) or (%.4f,+%.4f)\n",xs[3],xs_err[3],errL,errU);
   	}
   	else{
   		printf("efficiency = %.4f fixed\n",xs[3]);
   	}

	cout<<"*******************************************************"<<endl;
	cout<<endl;
 

	int status = minimum->Status();
	cout<<"status= "<<status<<endl;
	switch(status){
		case 0:
				cout<<"fit has problem"<<endl;break;
		case 1:
				cout<<"covariance is made pos defined"<<endl;break;
		case 2:
				cout<<"Hesse is invalid"<<endl;break;
		case 3:
				cout<<"Edm is about max"<<endl;break;
		case 4:
				cout<<"Reached call limit"<<endl;break;
		case 5:
				cout<<"Covariance is not positive defined"<<endl;break;
	}

	delete minimum;

    return 0;
}



void FitBetaOption(bool _fit_vb, bool _fit_efficiency, bool _fit_tao=true){
	fit_vb = _fit_vb;
	fit_efficiency = _fit_efficiency;
	fit_tao = _fit_tao;
}

double Calculated_efficiency(double _Ntof_beta, double _Ntof, double vb_s, double lamda_ref, double _timewindow_s){
	double ratio = _Ntof_beta/_Ntof;
	double exp_vb_T = TMath::Exp(-vb_s*_timewindow_s);
	double exp_lamda_T = TMath::Exp(-lamda_ref*_timewindow_s);
	return (ratio-(1-exp_vb_T))/(exp_vb_T*(1-exp_lamda_T)) ;
}


double Expect_b_tof_counts(double* T_window, double* par){//par0= lamda_ref; par1= vb /s background; par2 = efficiency; par3= Ntof; 
	double lamda_ref = par[0];
	double vb = par[1];
	double effi = par[2];
	double ntof = par[3];
	
	return ntof * (effi*(1-TMath::Exp(-(vb+lamda_ref)*T_window[0])) + (1-effi) * (1- TMath::Exp(-vb*T_window[0])) );

}


double Sensitive_T(double vb_s, double lamda_ref){

		return - TMath::Log(vb_s/(vb_s+lamda_ref))/lamda_ref;
}

//&&&&&&&&&&&&&&&&&& accumulated spectrum &&&&&&&&&&&&&&&&&&&&&&&&

double P_accu_curve(double* t, double* par){
//par0 = const background in [s]; par1 = efficiency; par2 = T1/2 [s]; par3=Ntof; par4= binwidth;

	double C_bg = par[0];
	double effi = par[1];
	double lamda = TMath::Log(2) / par[2];
	double ntof = par[3];
	double binwidth = par[4];

	return (C_bg + effi*lamda* TMath::Exp(-lamda*t[0]))*binwidth*ntof;  // return the expected beta count in the bin of t[0];

}



void Fit_accumulate_histo(TH1D** hin, double* par){
// par0= backgrounds in [s]; par1 = coincident window in [s]; par2 = efficiency;  par3=Ntof; par4 = halflife in ns;


	static TF1* f_dk_curve = NULL;
	static TCanvas* c_decay = NULL;

	if(f_dk_curve != NULL) delete f_dk_curve;
	
	if(c_decay!=NULL && c_decay->GetCanvasImp()!=NULL) delete c_decay; // window is closed
    c_decay = new TCanvas("c_decay","c_decay",600,400); 

    double halflife_s = par[4] * 1e-9; //==> from ns to s
    double lamda_ref = TMath::Log(2)/halflife_s ;

    double binwidth = halflife_s / binwidth_fraction_halflife;
    int nbins = TMath::Nint(par[1]/binwidth);

    double timeL = -nbins*binwidth;
    double timeR = nbins*binwidth;

    nbins*=2;
cout<<"delete hin"<<endl;
    if(*hin != NULL) delete *hin;
    *hin = new TH1D("hin","hin",nbins,timeL,timeR);

//    for(unsigned long i =0;i<T_data_neg.size();i++){(*hin)->Fill(T_data_neg[i]);}
    for(unsigned long i =0;i<T_data.size();i++){ (*hin)->Fill(T_data[i]);}

	f_dk_curve = new TF1("f_dk_curve",P_accu_curve,0,20*halflife_s,5);

	f_dk_curve->SetParameters(par[0],par[2],halflife_s,par[3],binwidth);

	f_dk_curve->SetParLimits(0,0,1);
	f_dk_curve->SetParLimits(1,0,0.80);
	f_dk_curve->SetParLimits(2,0.01,20);

	if(!fit_tao) f_dk_curve->FixParameter(2,halflife_s);
	if(!fit_vb) f_dk_curve->FixParameter(0,par[0]);
	if(!fit_efficiency) f_dk_curve->FixParameter(1,par[2]);
	f_dk_curve->FixParameter(3,par[3]);
	f_dk_curve->FixParameter(4,binwidth);

	f_dk_curve->SetParNames("C_background","efficiency","T1/2","Ntof","binwidth");

	TFitResultPtr result_fit = (*hin)->Fit(f_dk_curve,"LMES","",0,timeR);

	if(fit_tao){
		cout<<"T1/2 = "<<f_dk_curve->GetParameter(2)<<" ("<<f_dk_curve->GetParError(2)<<")  or (-"<<result_fit->LowerError(2)<<", +"
				<<result_fit->UpperError(2)<<")"<<endl;
	}
	else{
		cout<<"T1/2 = "<<f_dk_curve->GetParameter(2)<<" (fixed)"<<endl;
	}

	if(fit_vb){
		cout<<"C_background = "<<f_dk_curve->GetParameter(0)<<" ("<<f_dk_curve->GetParError(0)<<")  or (-"<<result_fit->LowerError(0)<<", +"
				<<result_fit->UpperError(0)<<")"<<endl;

	}
	else{
		cout<<"C_background = "<<f_dk_curve->GetParameter(0)<<" (fixed)"<<endl;
	}


	if(fit_efficiency){
		cout<<"efficiency = "<<f_dk_curve->GetParameter(1)<<" ("<<f_dk_curve->GetParError(1)<<")  or (-"<<result_fit->LowerError(1)<<", +"
				<<result_fit->UpperError(1)<<")"<<endl;		
	}
	else{
		cout<<"efficiency = "<<f_dk_curve->GetParameter(1)<<" (fixed)"<<endl;
	}

	c_decay->cd();
	(*hin)->SetTitle(Form("Beta decay;Decay time [s]; Counts / %.2f",binwidth));
	(*hin)->GetXaxis()->CenterTitle();
	(*hin)->GetYaxis()->CenterTitle();
	(*hin) ->Draw();
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




//%%%%%%%%%%%%%% method 1 ~ 4: only the first beta %%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%% method 5: all beta appear in window %%%%%%%%%%%%%%%

double ExtractLifetime(double* par, double Halflife_in, bool FirstPositiveBetaMode){// Make sure data is loaded in T_data vector in global;
// par0= backgrounds in [s]; par1 = coincident window in [s]; par2 = efficiency;  par3=Ntof
// halflife in ns

	vector<double> T_data_clone = T_data; // method 1-4 accept positive decay only; method 5 CAN be modified to accept both negative and positive, but the fitting range now still for dk_time >0 [s] !!!!

	T_data.clear(); T_data.shrink_to_fit();

	for(unsigned long i=0;i<T_data_clone.size();i++){
		if(T_data_clone[i] >0) T_data.push_back(T_data_clone[i]); // get positive only
	}

	double dk_time_avg=0;

	for(unsigned long i=0;i<T_data.size();i++){
			dk_time_avg+=T_data[i];
	}

	dk_time_avg /= T_data.size();


	double tao_ref[1]={Halflife_in*1e-9/TMath::Log(2)};

if(FirstPositiveBetaMode){

// method One: get only the expected life time
	double Par_set[4]={par[0],par[1],par[2],dk_time_avg};

	double tao_result = GetZeroX(ExtractTaobyDKtime,0.00001,100,Par_set,1e-6);

	cout<<"********************************"<<endl;
	cout<<"\e[1;33m"<<"Ref: T1/2 = "<<Halflife_in*1e-9<<" [s]; tao = "<<tao_ref[0]<<" [s]; lamda = "<<1/tao_ref[0]<<" /s"<<endl;
	cout<<"Expected Avg decay time = "<<Expect_DK_time(tao_ref,Par_set)<<endl;
	cout<<endl;
	cout<<"Data: vb = "<<Par_set[0]<<" /s; time window = "<<Par_set[1]<<" [s]; efficiency = "<<Par_set[2]<<endl;
	cout<<"Avg decay time = "<<Par_set[3]<<endl;
	cout<<"Ntof = "<<par[3]<<endl;

	double N_bg = par[3]*Par_set[0]*Par_set[1];
	double N_decay = (1-TMath::Exp(-Par_set[1]/tao_ref[0]))*par[3];
	double N_b_tof = (double)T_data.size();
	double Effi_err = TMath::Sqrt( TMath::Power(TMath::Sqrt(N_b_tof)/N_decay,2) + TMath::Power(TMath::Sqrt(N_bg)/N_decay,2) 
						  + TMath::Power((N_b_tof-N_bg)*TMath::Sqrt(N_decay)/(N_decay*N_decay),2) );


	cout<<"Expected background counts N_bg = "<<N_bg<<endl;
	cout<<"Expected N_decay = "<< N_decay<<endl;
	cout<<"efficiency = (N_b_tof - N_bg)/N_decay_expected = "<<(T_data.size() - par[3]*Par_set[0]*Par_set[1])/((1-TMath::Exp(-Par_set[1]/tao_ref[0]))*par[3])<<
			" ("<<Effi_err<<")"<<"\e[0m"<<endl;
cout<<"Sensitive T for efficiency measurement: T ="<<Sensitive_T(par[0],1./tao_ref[0])<<endl;
cout<<"Efficiency by calculate = "<<Calculated_efficiency(N_b_tof,par[3],par[0],1./tao_ref[0],par[1])<<endl;
	cout<<"********************************"<<endl;
	cout<<endl;


	cout<<"********************************"<<endl;

	cout<<"\e[1;33m"<<"By expected by decay time: T1/2 = "<<tao_result*TMath::Log(2)<<"\e[0m"<<endl;

	cout<<"********************************"<<endl;


//method Two: Likelyhood curve


	TF1* flog = new TF1("flog",LikelyhoodFcn,0.001,100,4);

	flog->SetParameters(par[0],par[1],par[2],par[3]);
	//flog->Draw();

	double gettao = flog->GetMinimumX(0.001,100);
	double gettao_up = flog->GetX(flog->Eval(gettao)+0.5,gettao,100);
	double gettao_low = flog->GetX(flog->Eval(gettao)+0.5,0.001,gettao);

	if(gettao>99.99){
		cout<<"**********************************************"<<endl;
		cout<<"\e[1;32m"<<"By Likelyhood curve maximum; Error==> Lmax - 0.5"<<endl;
		cout<<"Tao is larger than 100 s ????"<<"\e[0m"<<endl;
	}
	else if(gettao<0.0011){
		cout<<"**********************************************"<<endl;
		cout<<"\e[1;32m"<<"By Likelyhood curve maximum; Error==> Lmax - 0.5"<<endl;
		cout<<"Tao is larger than 0.001 s ????"<<"\e[0m"<<endl;
	}
	else{
		gettao *= TMath::Log(2); // to T1/2
		gettao_up *= TMath::Log(2);
		gettao_low *= TMath::Log(2);

		cout<<endl;
		cout<<"**********************************************"<<endl;
		cout<<"\e[1;32m"<<"By Likelyhood curve maximum; Error==> Lmax - 0.5"<<endl;
		cout<<"T1/2 = "<<gettao<<endl;
		cout<<"+sigma ="<< gettao_up - gettao<<endl;
		cout<<"-sigma ="<<gettao - gettao_low<<"\e[0m"<<endl;
	}


//method three: likelyhood unbinned fitting

	const int N_dk_counts=(int)T_data.size();

	double t_data[N_dk_counts];

	for(int i=0;i<N_dk_counts;i++){
		t_data[i] = T_data[i];
	}

	//ROOT::Fit::DataRange range(0,50);
    ROOT::Fit::UnBinData data((unsigned int)N_dk_counts, t_data);


    TF1* fdk_density = new TF1("fdk_density",P_density_dk_tao,0,100,4);
    fdk_density->SetParameters(Halflife_in*1e-9/TMath::Log(2),par[0],par[1],par[2]);

    ROOT::Math::WrappedMultiTF1 fitFunction( *fdk_density, fdk_density->GetNdim());

    ROOT::Fit::Fitter fitter;
    fitter.SetFunction( fitFunction, false);// false ==> not use customer defined gradient

    int Npars = fdk_density->GetNpar();
    double* Pars = fdk_density->GetParameters();
    fitter.Config().SetParamsSettings(Npars,Pars);
	fitter.Config().ParSettings(0).SetLimits(0,1000);
	if(!fit_tao)fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).SetLimits(0.00001,10);;
    if(!fit_vb)fitter.Config().ParSettings(1).Fix();
    fitter.Config().ParSettings(2).Fix();
    fitter.Config().ParSettings(3).SetLimits(0,0.34);
    if(!fit_efficiency)fitter.Config().ParSettings(3).Fix();

	fitter.Config().SetUpdateAfterFit();
	for(int iii=0;iii<10;iii++){
		fitter.LikelihoodFit(data);
	}

    TFitResult r=fitter.Result();
    r.Print();

    cout<<endl;
    cout<<"******************************** By unbinned fitting***********************"<<endl;
    cout<<"\e[1;37m"<<"T1/2 = "<<r.Parameter(0)* TMath::Log(2)<<endl;
    cout<<"-sigma = "<<r.LowerError(0)*TMath::Log(2)<<endl;
    cout<<"+sigma = "<<r.UpperError(0)*TMath::Log(2)<<"\e[0m"<<endl;
    cout<<"**************************************************************************"<<endl;
    cout<<endl;

    delete flog;
    delete fdk_density;



//method four:
    double par_all[5] ={tao_ref[0],par[0],par[1],par[2],par[3]};

    int getstatus = NumericalMinimization(par_all, MinmizerName.c_str(),AlgorithmName.c_str(),-1);// inpar0 = tao ; inpar1= backgrounds in [s]; 
	//inpar2 = coincident window in [s]; inpar3 = efficiency; inpar4=Ntof;

}// the end of if(FirstPositiveBetaMode)

//method five: fit the accumulate histogram (all beta between negative and positive window, not just the first positive beta)
    T_data.clear();  T_data.shrink_to_fit();
    T_data = T_data_clone;
	T_data_clone.clear();  T_data_clone.shrink_to_fit();

    cout<<"******************************************"<<endl;
    double par_fit_accumu[] = {par[0],par[1],par[2],par[3],Halflife_in};
    Fit_accumulate_histo(&h_beta_accumulate,par_fit_accumu);
    // par0= backgrounds in [s]; par1 = coincident window in [s]; par2 = efficiency;  par3=Ntof; par4 = halflife in ns;
    cout<<"******************************************"<<endl;
    cout<<endl;

    return 0;

}