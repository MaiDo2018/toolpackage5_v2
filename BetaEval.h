/* //////////////////////////////////////////////////////////////////////////////////////////////

Calculate suitable window for beta coincidence basing on beta background rate and Tof rate 

/////////////////////////////////////////////////////////////////////////////////////////////////
*/

#ifndef _BETAEVAL_H_
#define _BETAEVAL_H_ 1
#include "TF1.h"



double FOM(double* x, double* par){
	//x[0]=>deltaT in tao ; 
	//par0=>Ntof par1=>vb in tao; par2=> efficiency
	double num_tof = par[0];
	double vb = par[1];
	double effici= par[2];
	double Ndecay = effici * num_tof *(1-TMath::Exp(-x[0]));
	double Nback = num_tof * vb*x[0];// original without par[0]
	return Ndecay/(TMath::Sqrt(Ndecay+Nback));
}


double FOM_2(double* x, double* par){
	//x[0]=> time window in [s]
	//par0=vb in [s]; par1= lamda in [s]; par2= efficiency

	double T_w = x[0];
	double vb = par[0];
	double lamda = par[1];
	double effici = par[2];

	double P_true_beta = (1- TMath::Exp(-(vb+lamda)*T_w)) * effici*lamda / (vb+lamda);
	double P_eve = effici *(1- TMath::Exp(-(vb+lamda)*T_w)) + (1-effici)*(1-TMath::Exp(-vb*T_w));

	return P_true_beta/TMath::Sqrt(P_eve);


}



double BestDeltaT(double *x,double* par){// best delta T in unit of tao
	// par0=>NTof; par1 = efficiency
	//x[0]=>vb in tao
	TF1* f_FOM = new TF1("f_FOM",FOM,1e-2,10,3);
	f_FOM->SetParameters(par[0],x[0],par[1]);
	f_FOM->SetNpx(10000);
	double value_return = f_FOM->GetMaximumX(1e-2,10);
	delete f_FOM;
	return value_return;

}


double BestDeltaT_2(double *x,double* par){// best delta T in unit of second
	// par0=>lamda in second; par1 = efficiency
	//x[0]=>vb in second
	TF1* f_FOM = new TF1("f_FOM",FOM_2,1e-2,30,3);
	f_FOM->SetParameters(x[0],par[0],par[1]);
	f_FOM->SetNpx(10000);
	double value_return = f_FOM->GetMaximumX(1e-2,30);
	delete f_FOM;
	return value_return;

}

double BestDeltaT_const_background(double* x, double* par){
	//x[0] =>NTof
	//par[0] => vb
	return BestDeltaT(par,x);
}

double NTof_LowLimit(double* x,double* par){
	//x[0]=>NTof
	//par[0] =>vb; par[1]=> coefficience of deviation of Nback
	double deltaT = BestDeltaT_const_background(x,par);
	double NTof_limit = par[1] * TMath::Sqrt(par[0] * deltaT); //LD or Lc
	double epsilon=0.33;
	return NTof_limit / (epsilon*(1-TMath::Exp(-deltaT))); // low limit of NTof
}

double NTof_Fullfilled(double* x, double* par){
	//x[0]=>NTof
	//par[0] =>vb; par[1]=> coefficience of deviation of Nback
	return TMath::Abs(x[0] - NTof_LowLimit(x,par));

}

double NTof_OK_Vs_vb(double* x, double* par){
	//x[0]=>vb
	//par[0] => coefficience of deviation of Nback

	TF1* f_ntof_fullfilled = new TF1("f_ntof_fullfilled",NTof_Fullfilled,1e-2,1e5,2);
	f_ntof_fullfilled->SetParameters(x[0],par[0]);
	double value_return = f_ntof_fullfilled->GetMinimumX(1e-2,1e5);
	delete f_ntof_fullfilled;
	return value_return;

}


#endif