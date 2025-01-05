#ifndef _FUNCSJOHNSON_H_
#define _FUNCSJOHNSON_H_

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

class funcJohnson{

	private:
		TF1* funcJohn;

	public:
		static double Johnson_func(double* x, double* par);
		TF1* GetTF1(){return funcJohn;}

		funcJohnson(){
			//funcJohn = new TF1("funcJohn",this,&funcJohnson::Johnson_func,0,25e6,5,"1funcJohn","1fitfuncJohn");
			funcJohn = new TF1("funcJohn",funcJohnson::Johnson_func,0,25e6,5);
			funcJohn->SetParameters(0,0,0,1,-1);
		}

		~funcJohnson(){
			if(funcJohn!=NULL) delete funcJohn;
		}




};


double funcJohnson::Johnson_func(double* x, double* par){
	double Amp = par[0];
	double mu = par[1];
	double lambda=par[2];
	double gamma = par[3];
	double delta = par[4];

	double Z = gamma + delta*TMath::ASinH((x[0]-mu)/lambda);

	return Amp*TMath::Exp(-0.5*Z*Z)/TMath::Sqrt(1+TMath::Power(((x[0]-mu)/lambda),2));
}


#endif

