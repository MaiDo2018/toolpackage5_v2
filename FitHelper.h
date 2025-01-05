#ifndef _FITHELPER_H
#define _FITHELPER_H 1
#include <iostream>
#include "stdlib.h"
#include "TCanvas.h"

using namespace std;

//TCanvas ccc_fit;
//TH1D* h_random=NULL;


void SetParLimitAuto(TF1* ftemp, int _Norder,int _par_offset=0, bool fixh0=true){
	double tR = 0;
	double tL = 0;
	double countUpper=0;
	double countLower=0;

	ftemp->GetParLimits(_par_offset,countUpper,countLower);
	if(fixh0){
		if(countUpper==0 && countLower==0){
			double tem_high = ftemp->GetParameter(_par_offset);
			countLower =TMath::Max( tem_high-TMath::Sqrt(tem_high)*2 , 0.);
			countUpper = tem_high+TMath::Sqrt(tem_high)*2;
			ftemp->SetParLimits(_par_offset,countLower,countUpper);
		}
	}
	else{
		ftemp->SetParLimits(_par_offset,0,0);
	}


	for(int i=0;i<=_Norder;i++){
		if(i==0){
				ftemp->SetParLimits(_par_offset+2,0.0001,10000); //sigma>0
				tR = ftemp->GetParameter(_par_offset+3);
				tL = ftemp->GetParameter(_par_offset+4);
				ftemp->SetParLimits(_par_offset+3,0.01,10000); //tR>0
				ftemp->SetParLimits(_par_offset+4,-10000,-0.01);
		}
		else if(i==1){
				ftemp->SetParLimits(_par_offset+5,0.01,10); //k>0
				//ftemp->SetParLimits(_par_offset+6,tR,1000); // old method
				double temdiff = ftemp->GetParameter(_par_offset+6) - tR;
				if(temdiff>0) ftemp->SetParLimits(_par_offset+6,tR, ftemp->GetParameter(_par_offset+6)+2*tR);
				else ftemp->SetParLimits(_par_offset+6,tR, 3*tR);

		}
		else if(i==2){
				ftemp->SetParLimits(_par_offset+7,-10,-0.01); //k<0
				//ftemp->SetParLimits(_par_offset+8,-1000,tL); //old method
				double temdiff = ftemp->GetParameter(_par_offset+8) - tL;
				if(temdiff<0) ftemp->SetParLimits(_par_offset+8, ftemp->GetParameter(_par_offset+8)+2*tL ,tL);
				else ftemp->SetParLimits(_par_offset+8, 3*tL,tL);
		}
		else{
			if( i%2 !=0){//odd order
				ftemp->SetParLimits(_par_offset+5+2*(i-1),0.01,10); //k>0
				if(i>3)ftemp->SetParLimits(_par_offset+5+2*(i-1),0.0001,2);
				//ftemp->SetParLimits(_par_offset+6+2*(i-1), ftemp->GetParameter(_par_offset+6+2*(i-1)-4) , 1000);	//t3_5_7.....
				ftemp->SetParLimits(_par_offset+6+2*(i-1), ftemp->GetParameter(_par_offset+6+2*(i-1)-4) , ftemp->GetParameter(_par_offset+6+2*(i-1)-4)*3);
			}
			else{
				ftemp->SetParLimits(_par_offset+5+2*(i-1),-10,-0.01);
				if(i>4)ftemp->SetParLimits(_par_offset+5+2*(i-1),-2,-0.0001);
				//ftemp->SetParLimits(_par_offset+6+2*(i-1), -1000,ftemp->GetParameter(_par_offset+6+2*(i-1)-4));//t4_6_8
				ftemp->SetParLimits(_par_offset+6+2*(i-1), ftemp->GetParameter(_par_offset+6+2*(i-1)-4)*3,ftemp->GetParameter(_par_offset+6+2*(i-1)-4));
			}

		}
	}
}

void FitHelper(TF1* ftemp, TH1D* htemp, double fitlowedge, double fithighedge, int _Norder=0, string fitopt="LMEQ"){//ccc_fit.cd();
	double Maxcount = htemp->GetBinContent(htemp->GetMaximumBin());

	int Nbins = htemp->GetNbinsX();
	double hlowedge = htemp->GetBinLowEdge(1);
	double hupedge = htemp->GetBinLowEdge(Nbins)+ htemp->GetBinWidth(1);
	TH1D*  h_random= new TH1D("h_random","h_random",Nbins,hlowedge,hupedge);

	double ftemp_range_L=0,ftemp_range_R=0;
	ftemp->GetRange(ftemp_range_L,ftemp_range_R);
	ftemp->SetRange(hlowedge,hupedge);
	ftemp->SetNpx(500);

	int Nrandoms =0; double usescale=1;
	if(Maxcount<100){Nrandoms=100;  usescale = Maxcount*1./Nrandoms;}
	else Nrandoms = Maxcount;
	

//long counter=0;
	while(h_random->GetBinContent(h_random->GetMaximumBin()) < Nrandoms){
		h_random->Fill(ftemp->GetRandom(hlowedge,hupedge));
		//counter++;
		//if(counter%1000)cout<<"Fill = "<<counter<<"\r"<<flush;
	}

/*h_random->Draw();
ccc_fit.Modified();ccc_fit.Update();*/
char ctem;

	h_random->Scale(usescale);

	SetParLimitAuto(ftemp,_Norder);
	fitopt= fitopt + "N";
	h_random->Fit(ftemp,fitopt.c_str(),"",fitlowedge,fithighedge);

	SetParLimitAuto(ftemp,_Norder);

	for(int i=0;i<10;i++){
		if(i<2 || i>9) h_random->Fit(ftemp,"MNQ","",fitlowedge,fithighedge);
		else h_random->Fit(ftemp,"LMNQ","",fitlowedge,fithighedge);
		SetParLimitAuto(ftemp,_Norder);

	}


//ccc_fit.Modified();ccc_fit.Update();
//cout<<"Press any key to continue......"<<endl;
//cin>>ctem;
sleep(1); // Seem this pause for the system can help to prevent the clash of the fit

	TH1D* h_sub = (TH1D*)h_random->Clone();

	h_sub->Add((TH1*)htemp,-1);

	int scalefactor=70;

	h_sub->Scale(1./scalefactor);

	for(int i=0;i<scalefactor;i++){//char pc;
		h_random->Add((TH1*)h_sub,-1); h_random->Draw();
		cout<<"i ="<<i<<"\r"<<flush;
		int j=0;
		while(j++<10){
			if(j<3 ||j >9) h_random->Fit(ftemp,"MEQN","",fitlowedge,fithighedge);//fitopt.c_str()
			else h_random->Fit(ftemp,"LMEQN","",fitlowedge,fithighedge);
			SetParLimitAuto(ftemp,_Norder);

		}
//cin>>pc;
//ccc_fit.Update();ccc_fit.Modified();
		SetParLimitAuto(ftemp,_Norder);

	}


	ftemp->SetRange(ftemp_range_L,ftemp_range_R);

sleep(1);

	delete h_sub;
	delete h_random;


}

#endif

