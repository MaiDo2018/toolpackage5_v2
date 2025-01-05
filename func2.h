#ifndef _FUNC2_h_
#define _FUNC2_h_
#include "TMath.h"


class func2{
 	public:
		func2();
		~func2();
		static Double_t linearbg(double *x, double*par);
		static Double_t expbg(double *x, double* par);
		static Double_t gas(double*x,double*par);
		static Double_t exR(double*x,double*par);
		static Double_t exL(double*x,double*par);
		static Double_t fitfunc(double*x,double*par);
		static Double_t multifunc2(double *x,double *par);
		static Double_t multifunc3(double *x,double *par);
		static Double_t fitfunc_line(double*x,double*par);
		static Double_t fitfunc_exp(double*x,double*par);
		static Double_t multifunc2_exp(double*x,double*par);
		static Double_t multifunc2_line(double*x,double*par);
		static double intensity_cal(TF1* temline,double rangeL,double rangeR, double binwidth);
		static Double_t find_x_to_y(TF1* inputline,double y_value,double xpoint1,double xpoint2);
		static Double_t Rm_cal(TF1* inputline,int peaknum_index=1);

		static Double_t Rb90(double*x,double*par);

};




func2::func2(){}
func2::~func2(){}


Double_t func2::linearbg(double *x, double*par){

    return par[0] + par[1] * x[0];

}

Double_t func2::expbg(double *x, double* par){

	return TMath::Exp(par[0]+par[1]*x[0]);
}

Double_t func2::gas(double*x,double*par){
   return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2]));  

}

Double_t func2::exR(double*x,double*par){

	return par[0]*TMath::Exp(par[3]*(2*par[1]-2*x[0]+par[3])/(2*par[2]*par[2]));   // right long tail
	//return par[0]*TMath::Exp(par[3]*(-2*par[1]+2*x[0]+par[3])/(2*par[2]*par[2]));   // left long tail
 
}

Double_t func2::exL(double*x,double*par){

	//return par[0]*TMath::Exp(par[3]*(2*par[1]-2*x[0]+par[3])/(2*par[2]*par[2]));   // right long tail
	return par[0]*TMath::Exp(par[4]*(-2*par[1]+2*x[0]+par[4])/(2*par[2]*par[2]));   // left long tail
 
}

Double_t func2::fitfunc(double*x,double*par){  // par sequence=> high, cento, sigma, right tail tc, left tail tc
   /* for right long tail, all parameters>0 */
    if(x[0]<par[1]+par[3]&&x[0]>par[1]-par[4]) return gas(x,par); //+ linearbg(x,&par[4]);
    if(x[0]>=par[1]+par[3]) return exR(x,par); //+ linearbg(x,&par[4]);
    
 /* for left long tail , all parameters > 0*/
    if(x[0]<=par[1]-par[4]) return exL(x,par);
    //if(x[0]>par[1]-par[4]) return gas(x,par);
   return 0;
}



Double_t func2::multifunc2(double *x,double *par){
   
    return fitfunc(x,par) + fitfunc(x,&par[5]);// + fitfunc(x,&par[8]);// + fitfunc(x,&par[12]);

}

Double_t func2::multifunc3(double *x,double *par){
   
    return fitfunc(x,par) + fitfunc(x,&par[5]) + fitfunc(x,&par[10]);// + fitfunc(x,&par[12]);

}


Double_t func2::fitfunc_line(double*x,double*par){ // linear background + 1 gaus_exp
   
    return fitfunc(x,par) + linearbg(x,&par[5]);
    
}


Double_t func2::fitfunc_exp(double*x,double*par){  // exp background + 1 gaus_exp

     return fitfunc(x,par) + expbg(x,&par[5]);

}

Double_t func2::multifunc2_exp(double*x,double*par){

	return fitfunc(x,par) + fitfunc(x,&par[5]) + expbg(x,&par[10]);

}

Double_t func2::multifunc2_line(double*x,double*par){

	return fitfunc(x,par) + fitfunc(x,&par[5]) + linearbg(x,&par[10]);

}


Double_t func2::Rb90(double*x,double*par){
	double par_p1[5],par_p2[5];
	for(int i=0;i<5;i++){ par_p1[i] = par[i];}
	par_p2[0]=par[5];
	par_p2[1] = 1.000000638174*(par_p1[1]-130)+130; // 130 is assumed to be t0
	for(int i=2;i<5;i++){par_p2[i] = par_p1[i];}

	return fitfunc(x,par_p1) + fitfunc(x,par_p2);
}




double func2::intensity_cal(TF1* temline,double rangeL,double rangeR, double binwidth){
	double sum = 0;
	for(double bin_pos=rangeL;bin_pos<rangeR+binwidth;bin_pos=bin_pos+binwidth){
		sum = sum+temline->Eval(bin_pos);
	}
	return sum;
}

// enable if use in preview2.C ; disable when use in preview1.C because it is defined in preview1.h already
Double_t func2::find_x_to_y(TF1* inputline,double y_value,double xpoint1,double xpoint2){  // default xpoint1 close to cento
	if(y_value > inputline->Eval(xpoint1) && y_value > inputline->Eval(xpoint2) ) return -1;
	if(y_value < inputline->Eval(xpoint1) && y_value < inputline->Eval(xpoint2) ) return -1;

	if(y_value > inputline->Eval(xpoint1)){//default xpoint1 = 1 sigma position, close to center
		double tem_point = xpoint1;
		xpoint1 = xpoint2;
		xpoint2=tem_point;
	}

	while(TMath::Abs(xpoint1-xpoint2)>0.02){  // solution accuracy = 0.02
		double tem_middle = (xpoint1+xpoint2)*0.5;
		if(inputline->Eval(tem_middle) > y_value) xpoint1 = tem_middle;
		else xpoint2 = tem_middle;
	}
	
	return (xpoint1+xpoint2)*0.5;
}

Double_t func2::Rm_cal(TF1* inputline,int peaknum_index){
	TF1* temline = new TF1("temline",fitfunc,0,1e9,5);
	int para_offset=5*(peaknum_index-1);
	temline->SetParameter(0,inputline->GetParameter(0+para_offset));
	temline->SetParameter(1,inputline->GetParameter(1+para_offset));
	temline->SetParameter(2,inputline->GetParameter(2+para_offset));
	temline->SetParameter(3,inputline->GetParameter(3+para_offset));
	temline->SetParameter(4,inputline->GetParameter(4+para_offset));

	double peak_high = temline->GetParameter(0);
	double peak_center = temline->GetParameter(1);
	double peak_sigma = temline->GetParameter(2);
	double peak_tcR = temline->GetParameter(3);
	double peak_tcL = temline->GetParameter(4);

	if(peak_high<0 || peak_center <0 || peak_sigma<0 || peak_tcR<0 || peak_tcL<0){ temline->Delete(); return -1;}  // fitcurve fail

	double left_points_high[2] = { temline->Eval(peak_center-peak_sigma) , temline->Eval(peak_center-3.5*peak_sigma) };
	double right_points_high[2] = { temline->Eval(peak_center+peak_sigma) , temline->Eval(peak_center+3.5*peak_sigma) };

	bool fit_suc1 = peak_high > left_points_high[0] && left_points_high[0]>left_points_high[1];   // left half feature, judge fit success or not
	bool fit_suc2 = peak_high > right_points_high[0] && right_points_high[0]>right_points_high[1]; // right half feature,judge fit success or not

	if(fit_suc1 && fit_suc2){
		double left_edge = find_x_to_y(temline,peak_high*0.5,peak_center-peak_sigma,peak_center-3.5*peak_sigma);
		double right_edge = find_x_to_y(temline,peak_high*0.5,peak_center+peak_sigma,peak_center+3.5*peak_sigma);
		cout<<"FWHM= "<<right_edge-left_edge<<endl;
		if(left_edge>0 && right_edge >0){ return peak_center/(2*(right_edge-left_edge));} // success output, -1=> fail
		else{ temline->Delete(); return -1; }
	}
	else{  // fit curve unreasonable, break
		temline->Delete(); 
		return -1;
	}
}

#endif
