#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <time.h>
#include "TROOT.h"


#ifndef _BATCHMASSCAL_
#include "BatchMassCal.C"
#endif


using namespace std;

double* WeightMean(double* value, double* err, int Ninputs);
double* MOQ2mass(double moq,double moq_err, int _charge, double _me, double _me_err);
double* Mass2moq(double _mass, double _mass_err, int _charge, double _me, double _me_err);
double* MOQ_x_cal(double moq_r, double moq_r_err, double _rho, double _rho_err);


//this version now can accept both amu or Mass excess in keV as reference value

// %%%%%%%%%%%% Using input range ; not DrawingBit %%%%%%%%%%%%%%%%
void batchplot(string path="../",string filename="not include .txt",double mass_amu=1,double mass_err_amu=1,int indexStart=0,int indexStop=-1,double YMaximum=-1){ // index begins from 0 here

	struct tm *ptr_time;
	time_t time_now= time(0);
	ptr_time = localtime(&time_now);
	char time_str[20];
	strftime(time_str,20,"%Y%m%d%H%M%S",ptr_time);

	string filein = path + filename + ".txt";

	FILE *fp;
    fp = fopen(filein.c_str(),"r");
	if(fp==NULL){
			cout<<"\e[1;33m"<<"Error!! No .txt file import for plotting!!!"<<"\e[0m"<<endl;
			return;
	}

	char buffer[1000];
	char runname[20];
	char ionname[100];
	vector<double> mass_cal;
	vector<double> mass_cal_err;
	vector<double> delta_mass;
	vector<double> delta_mass_err;


	int charge_x;
	int ncounts;
	double tem_mass_cal;
	double tem_mass_cal_err;
	double tem_delta_mass;
	double tem_delta_mass_err;

	char refname[100];
	int charge_r;
	double massref_ame;
	double massref_err_ame;

	double getrho;
	double getrho_err;

	bool firstloopDone=false;

	fgets(buffer,1000,fp);

	//%%%%% just for read this bit from file; not really use %%%%
	int DrawingBit = 0;
	string buffer_s = buffer;
	bool haveDrawBit=false;
	if(buffer_s.find("Drawing") != string::npos) haveDrawBit=true;


	IonRecorder_perIon ion_each_run;
	Ionrecorder_MTIon ion_all_run;


	while(!feof(fp)){
		if(firstloopDone){
			mass_cal.push_back( tem_mass_cal);
			mass_cal_err.push_back(tem_mass_cal_err);
			delta_mass.push_back(tem_delta_mass);
			delta_mass_err.push_back(tem_delta_mass_err);

			Ref_Ion ref_info_run(massref_ame,massref_err_ame,charge_r,refname);
			ion_each_run.LoadData(runname,ionname,tem_mass_cal,tem_mass_cal_err,tem_delta_mass,tem_delta_mass_err,
										getrho,getrho_err,charge_x,ncounts,ref_info_run);

			// rename an ion with Xion+qx_Rion+qr
			ion_each_run.ionname = Form("%s+%d_%s+%d",ion_each_run.ionname.c_str(),ion_each_run.charge,ion_each_run.ref_info.name.c_str(),ion_each_run.ref_info.charge);
			ion_all_run.LoadPerIon(ion_each_run);

			//clear
			for(int i=0;i<100;i++){
				ionname[i]='\0';  refname[i]='\0';
			}

			ion_each_run.clear();

		}

		if(haveDrawBit){
			fscanf(fp,"%s\t%s\t%d\t%d\t%lf\t%lf\t",runname,ionname,&charge_x,&ncounts,&tem_mass_cal,&tem_mass_cal_err);
			fscanf(fp,"%lf\t%lf\t%s\t%d\t%lf\t%lf\t",&tem_delta_mass,&tem_delta_mass_err,refname,&charge_r,&massref_ame,&massref_err_ame);
			fscanf(fp,"%lf\t%lf%d",&getrho,&getrho_err,&DrawingBit);
		}
		else{
			fscanf(fp,"%s\t%s\t%d\t%d\t%lf\t%lf\t",runname,ionname,&charge_x,&ncounts,&tem_mass_cal,&tem_mass_cal_err);
			fscanf(fp,"%lf\t%lf\t%s\t%d\t%lf\t%lf\t",&tem_delta_mass,&tem_delta_mass_err,refname,&charge_r,&massref_ame,&massref_err_ame);
			fscanf(fp,"%lf\t%lf",&getrho,&getrho_err);
		}
		firstloopDone=true;
/*
printf("%s\t%s\t%d\t%d\t%.2f\t%.2f\n",runname,ionname,charge_x,ncounts,tem_mass_cal,tem_mass_cal_err);
printf("%f\t%f\t%s\t%d\t%.2f\t%.2f\n",tem_delta_mass,tem_delta_mass_err,refname,charge_r,massref_ame,massref_err_ame);
printf("%.10f\t%.10f\t%d\n",getrho,getrho_err,DrawingBit);
kkkk++;
		if(kkkk==2)return;*/
	}

	fclose(fp);
	printf("%s \t %s\n",runname,ionname);

	double * mass_cal2plot=NULL, *mass_cal_err2plot=NULL;

	int counter=0;    // record how many groups of data will be plotted actually;
	if(indexStop==-1){
		indexStop = mass_cal.size()-1;  // index begins from 0 here;
	}

	mass_cal2plot = new double[indexStop-indexStart+1];
	mass_cal_err2plot = new double[indexStop-indexStart+1];
	int index=0;
	for(unsigned int i=0;i<mass_cal.size();i++){
		if(i>=indexStart && i<=indexStop){
				mass_cal2plot[index]=mass_cal[i];
				mass_cal_err2plot[index++]=mass_cal_err[i];
				counter++;
		}
	}
	

	for(int i=0;i<counter;i++){
		printf("%.4f(%.4f)\n",mass_cal2plot[i],mass_cal_err2plot[i]);
	}


	string fileout = Form("%splotnote_%s.txt",path.c_str(),time_str);
	fp = fopen(fileout.c_str(),"w");
	fprintf(fp,"input file: %s \n",filename.c_str());

	if(indexStop== mass_cal.size()-1){
		fprintf(fp,"plot index range(begins from 0): %d - end \n",indexStart);
	}
	else{
		fprintf(fp,"plot index range(begins from 0): %d - %d \n",indexStart,indexStop);
	}

		
	vector<double>getprintdata = PlotMultiMassResult(counter,mass_cal2plot,mass_cal_err2plot,mass_amu,mass_err_amu,YMaximum);
	fprintf(fp,"Mass Mean[micro_amu]= %.4f(%.4f);\t\n",getprintdata[0],getprintdata[1]);
	fprintf(fp,"Mass Mean[ME in keV]= %.4f(%.4f);\t\n",getprintdata[2],getprintdata[3]);
	fprintf(fp,"Mass_mean - Mass_ame\t\n");
	fprintf(fp,"mass deviate[micro_amu]:  %.4f(%.4f);\t\n",getprintdata[4],getprintdata[5]);
	fprintf(fp,"mass deviate[ME in keV]: %.4f(%.4f);\t\n",getprintdata[6],getprintdata[7]);
	fprintf(fp,"Birge ratio: %.4f\t\n",getprintdata[8]);
	fprintf(fp,"Modified Birge Mass mean[micro_amu]: %.4f(%.4f)\t\n",getprintdata[0],getprintdata[9]);
	fprintf(fp,"0.95 confidential interval coefficient q= %.4f ==> [+/- q * MB mass mean err]\t\n",getprintdata[10]);

	fclose(fp);

	MassGraph->SetTitle(ionname);
	gStyle->SetTitleY(0.95);
	c2->Modified();
	c2->Update();
	c2->SaveAs(Form("%splot_%s.png",path.c_str(),time_str));


	delete[] mass_cal2plot;
	delete[] mass_cal_err2plot;


	fileout = Form("%srho_%s.txt",path.c_str(),time_str);
	fp = fopen(fileout.c_str(),"w");
	fprintf(fp,"input file: %s \n",filename.c_str());
	fprintf(fp,"casename\tchargeX\tRefmass\tRefmass_err\tchargeR\trho\trho_err\tmassX_uamu\t massX_err_uamu\tmassX_kev\tmassX_err_kev\n");


	vector<double> mass_cal_item;  // for calculate mass mean for all
	vector<double> mass_cal_err_item;
	double mass_combine;
	double mass_err_combine;


	for(int i=0;i<ion_all_run.NumIonSpecies;i++){
		string name_case = ion_all_run.ionlist[i];
		double ref_mass;
		double ref_mass_err;
		int ref_charge;
		int x_charge;

		double* rho_result;
		double mean_rho;
		double mean_rho_err;

		double mean_m_amu;
		double mean_m_err_amu;

		double mean_m_kev;
		double mean_m_err_kev;

		int N_same_case = (int)ion_all_run.MTIon[i].size();
		double* get_rho_i = new double[N_same_case]; // store rhos for each case of Xion+qx_Rion+qr
		double* get_rho_i_err = new double[N_same_case];

		for(int j=0;j<N_same_case;j++){
			ref_mass = ion_all_run.MTIon[i][j].ref_info.mass;
			ref_mass_err = ion_all_run.MTIon[i][j].ref_info.mass_err;
			ref_charge = ion_all_run.MTIon[i][j].ref_info.charge;
			
			x_charge = ion_all_run.MTIon[i][j].charge;
			get_rho_i[j] = ion_all_run.MTIon[i][j].rho;
			get_rho_i_err[j] = ion_all_run.MTIon[i][j].rho_err;
		}

		rho_result = WeightMean(get_rho_i,get_rho_i_err,N_same_case);

		mean_rho = rho_result[0];
		mean_rho_err = rho_result[1];


		double* get_moq_r = Mass2moq(ref_mass,ref_mass_err,ref_charge,m_ele,err_me); // m_ele,err_me from preview.C

		double* get_moq_x = MOQ_x_cal(get_moq_r[0],get_moq_r[1],mean_rho,mean_rho_err);
		double* get_case_mass = MOQ2mass(get_moq_x[0],get_moq_x[1],x_charge,m_ele,err_me);

		mean_m_amu = get_case_mass[0];
		mean_m_err_amu = get_case_mass[1];

		mass_cal_item.push_back(mean_m_amu);
		mass_cal_err_item.push_back(mean_m_err_amu);

		double* get_mass_kev = MassExcess(mean_m_amu,mean_m_err_amu,false);

		mean_m_kev = get_mass_kev[0];
		mean_m_err_kev = get_mass_kev[1];

		fprintf(fp,"%s\t%d\t%.4f\t%.4f\t%d\t%.10f\t%.10f\t",name_case.c_str(),x_charge,ref_mass,ref_mass_err,ref_charge,mean_rho,mean_rho_err);
		fprintf(fp,"%.4f\t%.4f\t%.4f\t%.4f\n",mean_m_amu,mean_m_err_amu,mean_m_kev,mean_m_err_kev);

		delete[] get_rho_i;
		delete[] get_rho_i_err;

	}// end of loop for all case


	int N_mass_item = (int)mass_cal_item.size();
	double* mass_item_i=NULL;
	double* mass_item_err = NULL;



	if(N_mass_item>1){
		mass_item_i = new double[N_mass_item];
		mass_item_err = new double[N_mass_item];
		for(int i=0;i<N_mass_item;i++){
			mass_item_i[i] = mass_cal_item[i];
			mass_item_err[i] = mass_cal_err_item[i];
		}

		double* get_m_final = WeightMean(mass_item_i,mass_item_err,N_mass_item);
		mass_combine = get_m_final[0];
		mass_err_combine = get_m_final[1];
	}
	else{
		mass_combine = mass_cal_item[0];
		mass_err_combine = mass_cal_err_item[0];
	}

	delete[] mass_item_i;
	delete[] mass_item_err;

	double* get_ME_final = MassExcess(mass_combine,mass_err_combine,false);

	double mass_combine_kev = get_ME_final[0];
	double mass_err_combine_kev = get_ME_final[1];

	string summary = "combine";

	fprintf(fp,"%s\t%d\t%.4f\t%.4f\t%d\t%.10f\t%.10f\t",summary.c_str(),0,0.,0.,0,0.,0.);
	fprintf(fp,"%.4f\t%.4f\t%.4f\t%.4f\n",mass_combine,mass_err_combine,mass_combine_kev,mass_err_combine_kev);

	fclose(fp);

}




//%%%%%%%%%%%%%%%% Using Drawing Bit %%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%% can accept both amu or Mass excess in keV as reference value %%%%%%%
void batchplot_Bit(string path="../",string filename="not include .txt",double mass_amu=1,double mass_err_amu=1,double YMaximum=-1){ // index begins from 0 here

	struct tm *ptr_time;
	time_t time_now= time(0);
	ptr_time = localtime(&time_now);
	char time_str[20];
	strftime(time_str,20,"%Y%m%d%H%M%S",ptr_time);

	string filein = path + filename + ".txt";

	FILE *fp;
    fp = fopen(filein.c_str(),"r");
	if(fp==NULL){
			cout<<"\e[1;33m"<<"Error!! No .txt file import for plotting!!!"<<"\e[0m"<<endl;
			return;
	}

	char buffer[1000];
	char runname[20];
	char ionname[100];
	vector<double> mass_cal;
	vector<double> mass_cal_err;
	vector<double> delta_mass;
	vector<double> delta_mass_err;
	vector<string> runname_v;
	vector<int> drawingbit_v;

	double tem_mass_cal;
	double tem_mass_cal_err;
	double tem_delta_mass;
	double tem_delta_mass_err;

	bool firstloopDone=false;

	fgets(buffer,1000,fp);

	//%%%%% just for read this bit from file; not really use %%%%
	int DrawingBit = 0;
	string buffer_s = buffer;
	bool haveDrawBit=false;
	if(buffer_s.find("Drawing") != string::npos){ haveDrawBit=true;}
	else{cout<<"\e[1;33m<<"<<"This file has no drawingbit; Break!!"<<"\e[0m"<<endl; return;}


	while(!feof(fp)){
		if(firstloopDone){
			mass_cal.push_back( tem_mass_cal);
			mass_cal_err.push_back(tem_mass_cal_err);
			delta_mass.push_back(tem_delta_mass);
			delta_mass_err.push_back(tem_delta_mass_err);
			runname_v.push_back(runname);
			drawingbit_v.push_back(DrawingBit);
		}

		if(haveDrawBit){
			fscanf(fp,"%s\t%s\t%lf\t%lf\t%lf\t%lf%d",runname,ionname,&tem_mass_cal,&tem_mass_cal_err,&tem_delta_mass,&tem_delta_mass_err,&DrawingBit);
		}
		else{
			fscanf(fp,"%s\t%s\t%lf\t%lf\t%lf\t%lf",runname,ionname,&tem_mass_cal,&tem_mass_cal_err,&tem_delta_mass,&tem_delta_mass_err);
		}
		firstloopDone=true;

	}

	fclose(fp);
	printf("%s \t %s\n",runname,ionname);

	double * mass_cal2plot=NULL, *mass_cal_err2plot=NULL;

	int counter=0;    // record how many groups of data will be plotted actually;
	for(unsigned int i=0;i<drawingbit_v.size();i++){
		if(drawingbit_v[i]==1) counter++;
	}

	printf("Get %d run files to draw.\n",counter);
	cout<<endl;

	mass_cal2plot = new double[counter];
	mass_cal_err2plot = new double[counter];
	int index=0;
	for(unsigned int i=0;i<mass_cal.size();i++){
		if(drawingbit_v[i] == 1){
				mass_cal2plot[index]=mass_cal[i];
				mass_cal_err2plot[index++]=mass_cal_err[i];
				printf("Get index %d: %.4f(%.4f)\n",i,mass_cal[i],mass_cal_err[i]);
		}
	}
	

	string fileout = Form("%splotnote_%s.txt",path.c_str(),time_str);
	fp = fopen(fileout.c_str(),"w");
	fprintf(fp,"input file: %s \n\n",filename.c_str());

	for(unsigned int i=0;i<drawingbit_v.size();i++){
		if(drawingbit_v[i] == 1){
			fprintf(fp,"%s \n",runname_v[i].c_str());
		}
	}


	
	vector<double>getprintdata = PlotMultiMassResult(counter,mass_cal2plot,mass_cal_err2plot,mass_amu,mass_err_amu,YMaximum);
	fprintf(fp,"Mass Mean[micro_amu]= %.4f(%.4f);\t\n",getprintdata[0],getprintdata[1]);
	fprintf(fp,"Mass Mean[ME in keV]= %.4f(%.4f);\t\n",getprintdata[2],getprintdata[3]);
	fprintf(fp,"Mass_mean - Mass_ame\t\n");
	fprintf(fp,"mass deviate[micro_amu]:  %.4f(%.4f);\t\n",getprintdata[4],getprintdata[5]);
	fprintf(fp,"mass deviate[ME in keV]: %.4f(%.4f);\t\n",getprintdata[6],getprintdata[7]);
	fprintf(fp,"Birge ratio: %.4f\t\n",getprintdata[8]);
	fprintf(fp,"Modified Birge Mass mean[micro_amu]: %.4f(%.4f)\t\n",getprintdata[0],getprintdata[9]);
	fprintf(fp,"0.95 confidential interval coefficient q= %.4f ==> [+/- q * MB mass mean err]\t\n",getprintdata[10]);

	fclose(fp);

	MassGraph->SetTitle(ionname);
	gStyle->SetTitleY(0.95);
	c2->Modified();
	c2->Update();
	c2->SaveAs(Form("%splot_%s.png",path.c_str(),time_str));


	delete[] mass_cal2plot;
	delete[] mass_cal_err2plot;

}



// common formula


double* WeightMean(double* value, double* err, int Ninputs){
	static double result_cal[2];
	result_cal[0]=0; // mean
	result_cal[1]=0;  // mean_err;

	double sum_weight=0;
	double sum_weighted_value=0;

	for(int i=0;i<Ninputs;i++){
		double weight_i = 1/(err[i]*err[i]);
		sum_weight+=weight_i;
		sum_weighted_value += value[i] * weight_i;
	}
	
	result_cal[0] = sum_weighted_value / sum_weight;
	result_cal[1] = TMath::Sqrt(1/sum_weight);

	return result_cal;

}


double* MOQ2mass(double moq,double moq_err, int _charge, double _me, double _me_err){
	static double mass_result[2];  // value; err

	mass_result[0] = moq *_charge + _charge * _me;

	mass_result[1] = _charge * TMath::Sqrt(moq_err*moq_err + _me_err*_me_err);

	return mass_result;

}


double* Mass2moq(double _mass, double _mass_err, int _charge, double _me, double _me_err){
	static double moq_result[2]; // value; err
	moq_result[0] = (_mass - _charge*_me) / _charge;
	moq_result[1] = TMath::Sqrt((_mass_err/_charge) * (_mass_err/_charge) + _me_err*_me_err);

	return moq_result;

}

double* MOQ_x_cal(double moq_r, double moq_r_err, double _rho, double _rho_err){
	static double moqx[2];
	moqx[0] = _rho * _rho * moq_r;
	moqx[1] = _rho *TMath::Sqrt(TMath::Power(2*moq_r*_rho_err,2) +TMath::Power(_rho*moq_r_err,2));

	return moqx;
}