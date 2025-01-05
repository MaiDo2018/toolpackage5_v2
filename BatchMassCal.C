#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <time.h>
#include "TROOT.h"

#define _BATCHMASSCAL_

using namespace std;

struct Ref_Ion{
	double mass;   // atomic mass, including all electron
	double mass_err;
	int charge;
	string name;

	void Load(double _mass,double _mass_err, int _charge, string _name){
		mass = _mass;
		mass_err = _mass_err;
		charge = _charge;
		name = _name;
	}


	Ref_Ion(double _mass,double _mass_err, int _charge, string _name){
		Load(_mass,_mass_err,_charge,_name);
	}

	Ref_Ion(){
		Load(-1,-1,-1," ");
	}

	void Clear(){
		Load(-1,-1,-1," ");
	}

};


class IonRecorder_perIon{  // store info of a single ion
	public:
		string runname;
		string ionname;
		int charge;
		int ncounts;
		double massresult;
		double massresult_err;
		double mass_delta;
		double mass_delta_err;
		double rho;
		double rho_err;

		Ref_Ion ref_info;
		

		void clear(){
			runname=" ";
			ionname=" ";
			charge=-1;
			ncounts=0;
			massresult=-1;
			massresult_err=-1;
			mass_delta=-1;
			mass_delta_err=-1;
			rho = -1;
			rho_err = -1;
			ref_info.Clear();
			
		}

		void LoadData(string _runname,string _ionname, double _massresult,double _massresult_err, double _mass_delta,double _mass_delta_err, 
							double _rho, double _rho_err, int _charge, int _ncounts, Ref_Ion _refion){
			runname= _runname;
			ionname= _ionname;
			massresult= _massresult;
			massresult_err=	 _massresult_err;
			mass_delta= _mass_delta;
			mass_delta_err= _mass_delta_err;
			rho = _rho;
			rho_err = _rho_err;
			charge = _charge;
			ncounts = _ncounts;
			ref_info = _refion;
		}

};

typedef vector<IonRecorder_perIon> VT_sameIon;  // store info of ions of the same species

class Ionrecorder_MTIon{ //classify ions as species MTIon[j][k]: ion species => jth in ionlist; k => the kth ion
	public:
		vector<string> ionlist;
		int NumIonSpecies;
		IonRecorder_perIon perIon;
		VT_sameIon sameIon;   // vector of IonRecorder_perIon
		vector<VT_sameIon> MTIon;
		
		void ClassifySave(){
			int IonSpeciesIndex=-1;
			for(int i=0;i<ionlist.size();i++){
				if(strcmp(perIon.ionname.c_str(),ionlist[i].c_str()) == 0){ IonSpeciesIndex=i; break;	}
			}

			if(IonSpeciesIndex==-1){ // new species ion input;
				ionlist.push_back(perIon.ionname);
				NumIonSpecies=ionlist.size();
				sameIon.push_back(perIon);
				MTIon.push_back(sameIon);
				sameIon.clear();
				perIon.clear();
			}
			else{// old ion species
				MTIon[IonSpeciesIndex].push_back(perIon);
				perIon.clear();
			}
		}


		void LoadPerIon(const char* _runname,const char* _ionname,double _massresult,double _massresult_err,double _mass_delta,double _mass_delta_err,
									double _rho, double _rho_err, int _charge, int _ncounts, Ref_Ion _refion){
			perIon.runname = _runname;
			perIon.ionname = _ionname;
			perIon.massresult = _massresult;
			perIon.massresult_err = _massresult_err;
			perIon.mass_delta = _mass_delta;
			perIon.mass_delta_err = _mass_delta_err;
			perIon.rho = _rho;
			perIon.rho_err = _rho_err;
			perIon.charge = _charge;
			perIon.ncounts = _ncounts;
			perIon.ref_info = _refion;
			ClassifySave();
		}

		void LoadPerIon(IonRecorder_perIon &_PerIon){
			LoadPerIon(_PerIon.runname.c_str(),_PerIon.ionname.c_str(),_PerIon.massresult,_PerIon.massresult_err,_PerIon.mass_delta,_PerIon.mass_delta_err,
								_PerIon.rho, _PerIon.rho_err, _PerIon.charge, _PerIon.ncounts, _PerIon.ref_info);
		}


		void clear(){
			ionlist.clear();
			NumIonSpecies=0;
			perIon.clear();
			sameIon.clear();
			MTIon.clear();
		}

		Ionrecorder_MTIon(){
			clear();
		}

};


class RunRecorder{ // store info of each run in .txt
	public:
		string runname;
		int NumofIon;
		vector<string> ionname;
		vector<int> charge;
		vector<double> tof;
		vector<double> tof_err;
		vector<int> counts;
		vector<double> mass_ame;
		vector<double> mass_ame_err;
		vector<int> refbit;

		vector<double> massresult;
		vector<double> massresult_err;
		vector<double> delta_mass;
		vector<double> delta_mass_err;
		vector<double> rho;
		vector<double> rho_err;

		void clear(){
			runname=" ";
			NumofIon=0;
			ionname.clear();
			charge.clear();
			tof.clear();
			tof_err.clear();
			counts.clear();
			mass_ame.clear();
			mass_ame_err.clear();
			refbit.clear();
			massresult.clear();
			massresult_err.clear();
			delta_mass.clear();
			delta_mass_err.clear();
			rho.clear();
			rho_err.clear();
		}

		void LoadData(IonRecorder_perIon &_perIon,int _NumofIon,int _charge,double _tof,double _tof_err,int _counts,double _mass_ame,double _mass_ame_err,int _refbit){
			runname= _perIon.runname;
			NumofIon = _NumofIon;
			ionname.push_back( _perIon.ionname);
			charge.push_back(_charge);
			tof.push_back(_tof);
			tof_err.push_back(_tof_err);
			counts.push_back(_counts);
			mass_ame.push_back(_mass_ame);
			mass_ame_err.push_back(_mass_ame_err);
			refbit.push_back(_refbit);
			massresult.push_back(_perIon.massresult);
			massresult_err.push_back(_perIon.massresult_err );
			delta_mass.push_back(_perIon.mass_delta);
			delta_mass_err.push_back(_perIon.mass_delta_err);
			rho.push_back(_perIon.rho);
			rho_err.push_back(_perIon.rho_err);
		}

		RunRecorder(){
			clear();
		}

};


IonRecorder_perIon EachIon;  // contain mass result of a single ion
Ionrecorder_MTIon MTIon_recorder; // contain run and mass info of multi species
RunRecorder Run_recorder; // contain info for each run
vector<RunRecorder> MTRun_recorder;  // contain info for all run



void BatchMassCal(string path="../to the folder of .txt/",int massnum=0){

	struct tm *ptr_time;
	time_t time_now= time(0);
	ptr_time = localtime(&time_now);
	char time_str[20];
	strftime(time_str,20,"%Y%m%d%H%M%S",ptr_time);

	string filename = Form("%sA%d.txt",path.c_str(),massnum);
	FILE *fp;
    fp = fopen(filename.c_str(),"r");
	if(fp==NULL){
			cout<<"\e[1;33m"<<"Error!! .txt file for mass calculate!!!"<<"\e[0m"<<endl;
			return;
	}

	char buffer[1000];
	fgets(buffer,1000,fp);
	if(strncmp(buffer,"filename",8) !=0){rewind(fp);}  // have no file head


	char runname[100];
	int NumInput=0;    // number of ion to input in this file
	char ionname[100];
	int charge=0;
	double ion_tof=0;
	double ion_tof_err=0;
	int counts=0;
	double ion_ame_mass=0;
	double ion_ame_mass_err=0;
	int refbit=0;  // 0 for calculate; 1 for reference;

	int NumOfRef=0;  // how many ions are reference

	bool FirstLoadDone=false;

	Ion *inputIon = NULL;
	Ion *rion = NULL;


	while(!feof(fp)){
		if(FirstLoadDone){
			int IndexRef=0;
			NumOfRef =0;

			for(int i=0;i<NumInput;i++){
				fscanf(fp,"%s\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%d",ionname,&charge,&ion_tof,&ion_tof_err,&counts,&ion_ame_mass,&ion_ame_mass_err,&refbit);
				//printf("%s;%d;%.4f;%.4f;%d;%.4f;%.4f;%d\n",ionname,charge,ion_tof,ion_tof_err,counts,ion_ame_mass,ion_ame_mass_err,refbit);

				inputIon[i].mass=ion_ame_mass;
				inputIon[i].mass_err=ion_ame_mass_err;
				inputIon[i].tof=ion_tof;
				inputIon[i].tof_err = ion_tof_err;
				inputIon[i].charge = charge;
				inputIon[i]._name = ionname;
				inputIon[i].counts = counts;
				//printf("%s;%d;%.4f;%.4f;%d;%.4f;%.4f;%d\n",inputIon[i]._name.c_str(),inputIon[i].charge,inputIon[i].tof,inputIon[i].tof_err,counts,inputIon[i].mass,inputIon[i].mass_err,refbit);

				if(refbit==1){  // is ref ion
					if(IndexRef==5){cout<<"warning!! too much ref input; omit inputs beyond 5"<<endl;}
					else{
						NumOfRef++;
						rion[IndexRef] = inputIon[i];
						IndexRef++;
					}
				}

				for(int j=0;j<100;j++){ // clear name container
					ionname[j]='\0';
				}

			}// end of for loop over all ion in a Run

			// handle data in a run
			printf("&&&&&&&&&&&&&&&&& %s &&&&&&&&&&&&&&&&&\n",runname);
			//FILE* fout;
			//string outfilename=Form("%sA%d_result.txt",path.c_str(),massnum);
			//fout=fopen()
			for(int i=0;i<NumInput;i++){ // calcalute mass
				cout<<i+1<<":\t"<<inputIon[i]._name<<endl;
				double* massresult = mass_calculator(inputIon[i].tof,inputIon[i].tof_err,inputIon[i].charge,rion,NumOfRef);
				//massresult0 =>mass_cal; massresult1 => mass_cal_err; massresult2 => rho_cal; massreslt3 => rho_cal_err;
				double mass_delta=massresult[0] - inputIon[i].mass;
				double mass_delta_err = TMath::Sqrt(massresult[1]*massresult[1] + inputIon[i].mass_err*inputIon[i].mass_err);
				printf("Delta mass = %.4f( %.4f ) \n ",mass_delta,mass_delta_err);
				cout<<endl;
				cout<<endl;
				cout<<endl;

				Ref_Ion rion_info(rion[0].mass, rion[0].mass_err, rion[0].charge, rion[0]._name);

				EachIon.LoadData(runname,inputIon[i]._name,massresult[0],massresult[1], mass_delta,mass_delta_err, massresult[2], massresult[3],
										inputIon[i].charge, inputIon[i].counts,rion_info);
				MTIon_recorder.LoadPerIon(EachIon);
				int bitref=0;
				for(int j=0;j<NumOfRef;j++){
					if( strcmp(inputIon[i]._name.c_str(),rion[j]._name.c_str()) ==0 ){bitref=1;break;}
				}

				// save current ion in the same run
				Run_recorder.LoadData(EachIon,NumInput,inputIon[i].charge,inputIon[i].tof,inputIon[i].tof_err,inputIon[i].counts,inputIon[i].mass,inputIon[i].mass_err,bitref);
			}


			MTRun_recorder.push_back(Run_recorder);		// save per run	

			Run_recorder.clear();	// clear to accept next run!!!!!!!!!!!!! important

			for(int j=0;j<100;j++){ // clear name container
					runname[j]='\0';
			}


			delete[] inputIon;
			delete[] rion;
		}


		fscanf(fp,"%s\t%d",runname,&NumInput);
		if(!feof(fp)){
				inputIon = new Ion[NumInput];  // Ion class from preview3.C
				rion = new Ion[5];   // maximum number of Ref ion can be accepted
		}
		else break;

		FirstLoadDone=true;
	}

cout<<"mark1"<<endl;

	//&&&&&&&&&&& save result to file on hard disk  &&&&&&&&&&&&&&&&&&&&&&
	string outfile = Form("%sA%d_%s.txt",path.c_str(),massnum,time_str);
	FILE* fout;
	fout=fopen(outfile.c_str(),"w");
	if(fout==NULL){printf("faile to save A%d_%s.txt\n",massnum,time_str); goto saveion;}
	fprintf(fout,"ionname\tcharge\ttof\ttof_err\tcounts\tmassAME\tmassAME_err\trefbit\tmassCal\tmassCal_err\tdelta_mass\tdelta_mass_err\trho\trho_err\n");

	for(int i=0;i<MTRun_recorder.size();i++){
		//printf("%s \t %d\n",MTRun_recorder[i].runname.c_str(),MTRun_recorder[i].NumofIon);
		fprintf(fout,"%s\t%d\n",MTRun_recorder[i].runname.c_str(),MTRun_recorder[i].NumofIon);
		double mass_cal=0;
		double mass_cal_err=0;
		double delta_mass=0;
		double delta_mass_err=0;
		double rho = 0;
		double rho_err=0;

		for(int j=0;j<MTRun_recorder[i].NumofIon;j++){
			strcpy(ionname,MTRun_recorder[i].ionname[j].c_str());
			charge = MTRun_recorder[i].charge[j];
			ion_tof=MTRun_recorder[i].tof[j];
			ion_tof_err=MTRun_recorder[i].tof_err[j];
			counts=MTRun_recorder[i].counts[j];
			ion_ame_mass=MTRun_recorder[i].mass_ame[j];
			ion_ame_mass_err=MTRun_recorder[i].mass_ame_err[j];
			refbit=MTRun_recorder[i].refbit[j];

			mass_cal=MTRun_recorder[i].massresult[j];
			mass_cal_err=MTRun_recorder[i].massresult_err[j];
			delta_mass=MTRun_recorder[i].delta_mass[j];
			delta_mass_err=MTRun_recorder[i].delta_mass_err[j];

			rho=MTRun_recorder[i].rho[j];
			rho_err=MTRun_recorder[i].rho_err[j];


			//printf("%s;%d;%.4f;%.4f;%d;%.4f;%.4f;%d;%.4f;%.4f;%.4f;%.4f\n",ionname,charge,ion_tof,ion_tof_err,counts,ion_ame_mass,ion_ame_mass_err,refbit,mass_cal,mass_cal_err,delta_mass,delta_mass_err);
			fprintf(fout,"%s\t%d\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t",ionname,charge,ion_tof,ion_tof_err,counts,ion_ame_mass,ion_ame_mass_err);
			fprintf(fout,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.10f\t%.10f\n",refbit,mass_cal,mass_cal_err,delta_mass,delta_mass_err,rho,rho_err);

		}
	}
	fclose(fout);

	cout<<"mark2"<<endl;

saveion:
	for(int i=0;i<MTIon_recorder.NumIonSpecies;i++){
		string nameDir = Form("A%d_%s",massnum,MTIon_recorder.ionlist[i].c_str());
		string openDir = ".!ls " + path + nameDir;
		string makeDir = ".!mkdir " + path + nameDir;
		int DirOK=0;
		gROOT->ProcessLine(openDir.c_str(),&DirOK);
		if(DirOK){gROOT->ProcessLine(makeDir.c_str());}
		outfile=Form("%s%s/A%d_%s_%s.txt",path.c_str(),nameDir.c_str(),massnum,MTIon_recorder.ionlist[i].c_str(),time_str);
		fout=fopen(outfile.c_str(),"w");
		fprintf(fout,"filename\tIon\tcharge\tncounts\tmassCal\tmassCal_err\tdelta_mass\tdelta_mass_err\tIon_r\tCharge_r\tmass_r_ame\tmass_r_err_ame\trho\trho_err\tDrawingBit(1:on;0:off)\n");
		if(fout==NULL){printf("faile to save A%d_%s_%s.txt\n",massnum,MTIon_recorder.ionlist[i].c_str(),time_str); goto theend;}
		string runname;
		string ionname;
		int charge_x;
		int count_x;
		double massresult;
		double massresult_err;
		double mass_delta;
		double mass_delta_err;
		string ion_r_name;
		int charge_r;
		double mass_r_ame;
		double mass_r_err_ame;
		double rho;
		double rho_err;

		for(int j=0;j<MTIon_recorder.MTIon[i].size();j++){
			runname=MTIon_recorder.MTIon[i][j].runname;
			ionname=MTIon_recorder.MTIon[i][j].ionname;
			charge_x = MTIon_recorder.MTIon[i][j].charge;
			count_x =MTIon_recorder.MTIon[i][j].ncounts;
			massresult=MTIon_recorder.MTIon[i][j].massresult;
			massresult_err=MTIon_recorder.MTIon[i][j].massresult_err;
			mass_delta=MTIon_recorder.MTIon[i][j].mass_delta;
			mass_delta_err=MTIon_recorder.MTIon[i][j].mass_delta_err;

			ion_r_name = MTIon_recorder.MTIon[i][j].ref_info.name;
			charge_r = MTIon_recorder.MTIon[i][j].ref_info.charge;
			mass_r_ame = MTIon_recorder.MTIon[i][j].ref_info.mass;
			mass_r_err_ame = MTIon_recorder.MTIon[i][j].ref_info.mass_err;

			rho = MTIon_recorder.MTIon[i][j].rho;
			rho_err = MTIon_recorder.MTIon[i][j].rho_err;

			//printf("%s;%s;%.4f(%.4f);%.4f(%.4f)\n",runname.c_str(),ionname.c_str(),massresult,massresult_err,mass_delta,mass_delta_err);
			fprintf(fout,"%s\t%s\t%d\t%d\t%.4f\t%.4f\t",runname.c_str(),ionname.c_str(),charge_x,count_x,massresult,massresult_err);
			fprintf(fout,"%.4f\t%.4f\t%s\t%d\t%.4f\t%.4f\t",mass_delta,mass_delta_err,ion_r_name.c_str(),charge_r,mass_r_ame,mass_r_err_ame);
			fprintf(fout,"%.10f\t%.10f\t%d\n",rho,rho_err,1);
		}
		//cout<<endl;
		//cout<<endl;
		fclose(fout);
	}


theend:
	fclose(fp);
	MTRun_recorder.clear();
	MTIon_recorder.clear();
}

//mass_calculator(double tof0, double tof0_err, int q_x, Ion* RefIon, int Num_RefIons)