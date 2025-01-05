// define global function for preview1.C
#ifndef _PREVIEW1_H_
#define _PREVIEW1_H_

//&&&&&&&&&&&&&& selection of AME database to be used &&&&&&&&&&&&&
//#define AME2016
#define AME2020
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#ifdef AME2016
#define TXT "AME2016.txt"
#endif
#ifdef AME2020
#define TXT "AME2020.txt"
#endif

#include "TF1.h"
#include "TMath.h"
#include <stdlib.h>
#include <sstream>
#include <string>
#include <algorithm>


Double_t* SearchEleByMass(double mass_input=0,const char* filename=TXT);
Double_t* SearchAME(int Anum =0,const char* element="" ,bool verbose =false,const char* filename= TXT);

/*
// calculate the x value corresponding to a given Y of a curve
Double_t find_x_to_y(TF1* inputline,double y_value,double xpoint1,double xpoint2){  // default xpoint1 close to cento
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

// resolution power calculator according to different fit func setting
#ifdef _SETFUNC1_
Double_t Rm_cal(TF1* inputline,int peaknum_index=1){
	TF1* temline = new TF1("temline",func1::fitfunc,0,1e9,4);
	int para_offset=4*(peaknum_index-1);
	temline->SetParameter(0,inputline->GetParameter(0+para_offset));
	temline->SetParameter(1,inputline->GetParameter(1+para_offset));
	temline->SetParameter(2,inputline->GetParameter(2+para_offset));
	temline->SetParameter(3,inputline->GetParameter(3+para_offset));

	double peak_high = temline->GetParameter(0);
	double peak_center = temline->GetParameter(1);
	double peak_sigma = temline->GetParameter(2);
	double peak_tc = temline->GetParameter(3);

	if(peak_high<0 || peak_center <0 || peak_sigma<0 || peak_tc<0){ temline->Delete(); return -1;}  // fitcurve fail

	double left_points_high[2] = { temline->Eval(peak_center-peak_sigma) , temline->Eval(peak_center-3.5*peak_sigma) };
	double right_points_high[2] = { temline->Eval(peak_center+peak_sigma) , temline->Eval(peak_center+3.5*peak_sigma) };

	bool fit_suc1 = peak_high > left_points_high[0] && left_points_high[0]>left_points_high[1];   // left half feature, judge fit success or not
	bool fit_suc2 = peak_high > right_points_high[0] && right_points_high[0]>right_points_high[1]; // right half feature,judge fit success or not

	if(fit_suc1 && fit_suc2){
		double left_edge = find_x_to_y(temline,peak_high*0.5,peak_center-peak_sigma,peak_center-3.5*peak_sigma);
		double right_edge = find_x_to_y(temline,peak_high*0.5,peak_center+peak_sigma,peak_center+3.5*peak_sigma);
		if(left_edge>0 && right_edge >0){ return peak_center/(2*(right_edge-left_edge));} // success output, -1=> fail
		else{ temline->Delete(); return -1; }
	}
	else{  // fit curve unreasonable, break
		temline->Delete(); 
		return -1;
	}

}

#endif



#ifdef _SETFUNC2_
Double_t Rm_cal(TF1* inputline,int peaknum_index=1){
	TF1* temline = new TF1("temline",func2::fitfunc,0,1e9,5);
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
		if(left_edge>0 && right_edge >0){ return peak_center/(2*(right_edge-left_edge));} // success output, -1=> fail
		else{ temline->Delete(); return -1; }
	}
	else{  // fit curve unreasonable, break
		temline->Delete(); 
		return -1;
	}
}
#endif
*/

#ifdef AME2016
Double_t* SearchAME(int Anum ,const char* element ,bool verbose,const char* filename){  // return an array finally, massreturn[0]=> mass; massreturn[1]=>err

	FILE *fp = NULL;
	static Double_t massreturn[2];
	massreturn[0]=-1;
	massreturn[1]=-1;

	int fscanstatus=0;

	fp = fopen(filename,"r");
	if(fp==NULL){ cout<< "fail to open file"<<endl;return massreturn; }

	//a1,i3, i5,i5,i5,1x,a3,a4,1x,f13.5,f11.5,f11.3,     f9.3,1x,a2,f11.3,f9.3,1x,i3,1x,f12.5,f11.5
	//cc NZ  N  Z  A    el  o     mass  unc binding unc     B  beta  unc    atomic_mass   unc

	bool initial_state = true;
	char cc[3]={0};
	int NZ=0;
	int N=0;
	int Z=0;
	int A=0;
	char space1=0;
	char el[3]={0};
	char orig[5]={0};
	char space2=0;
	double mass=0;
	char space3=0;
	double mass_unc=0;
	char space4=0;
	double binding=0;
	char space5=0;
	double binding_unc=0;
	char space6=0;
	char B[3]={0};
	//double beta=0;
	//char space7=0;
	//double beta_unc=0;
	//char space8=0;
	char betainfo[30];
	int headdig=0;
	double atomic_mass=0;
	char space9=0;
	double atomic_mass_unc=0;
	char space10=0;


while(!feof(fp)){

	if(!initial_state){
 		/*printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
 		printf("%2s%20s ",B,betainfo);
 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
 		if(space9=='#'){printf("%c\n",space10);}
 		else printf("\n");*/
 		
 	
 	
		cc[0]=0;cc[1]=0;cc[2]=0;
		NZ=0;
		N=0;
		Z=0;
		A=0;
		space1=0;
		el[0]=0;el[1]=0;el[2]=0;
		orig[0]=0;orig[1]=0;orig[2]=0;orig[3]=0;orig[4]=0;
		space2=0;
		mass=0;
		space3=0;
		mass_unc=0;
		space4=0;
		binding=0;
		space5=0;
		binding_unc=0;
		space6=0;
		B[0]=0;B[1]=0;B[2]=0;
		//beta=0;
		//space7=0;
		//beta_unc=0;
		//space8=0;
		for(int j=0;j<30;j++){betainfo[j] = 0;}
		headdig=0;
		atomic_mass=0;
		space9=0;
		atomic_mass_unc=0;
		space10=0;

	}

	fscanstatus=fscanf(fp,"%2c%2d %4d %4d %4d%1c%2c%1c%4c ",cc,&NZ,&N,&Z,&A,&space1,el,&space2,orig);  // %2c%2d => input 2 grid char and 2 grid int number without interval
	fscanstatus=fscanf(fp,"%13lf%1c %11lf%1c %11lf%1c %9lf%1c ",&mass,&space3,&mass_unc,&space4,&binding,&space5,&binding_unc,&space6);
	//fscanf(fp,"%2c %lf%1c %lf%1c ",B,&beta,&space7,&beta_unc,&space8);
	fscanstatus=fscanf(fp,"%2c%20c ",B,betainfo);
	fscanstatus=fscanf(fp,"%3d %12lf%1c %11lf",&headdig,&atomic_mass,&space9,&atomic_mass_unc);
	if(space9=='#')fscanstatus=fscanf(fp,"%1c",&space10);   // !!!! actually there is a '\n' at the end of each line that is why using %2c%2d to fix something wrong

	initial_state=false;

	if(A==Anum && strncmp(el,element,2)==0){
		if(verbose){
			printf("1N-Z \t N \t Z \t A \t EL \t Orig \t MASS EXCESS(keV) \t BINDING ENERGY/A (keV) \t  BETA-DECAY ENERGY(keV) \t ATOMIC MASS(micro-u)\n");
	 		printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
	 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
	 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
	 		printf("%2s%20s ",B,betainfo);
	 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
	 		if(space9=='#'){printf("%c\n",space10);}
 			else printf("\n");
 		}
 		massreturn[0]=headdig*1e6+atomic_mass;
 		massreturn[1]=atomic_mass_unc;
 		return massreturn;

	}

}	



	cout<<"no match Isotope, make sure 1st letter is CAPITAL	like: He "<<endl;
 	fclose(fp);
 	return massreturn;

}

#endif


#ifdef AME2020
Double_t* SearchAME(int Anum ,const char* element ,bool verbose,const char* filename){  // return an array finally, massreturn[0]=> mass; massreturn[1]=>err

	FILE *fp = NULL;
	static Double_t massreturn[2];
	massreturn[0]=-1;
	massreturn[1]=-1;

	int fscanstatus=0;

	fp = fopen(filename,"r");
	if(fp==NULL){ cout<< "fail to open file"<<endl;return massreturn; }

	//a1,i3, i5,i5,i5,1x,a3,a4,1x,f13.5,f11.5,f11.3,     f9.3,1x,a2,f11.3,f9.3,1x,i3,1x,f12.5,f11.5
	//cc NZ  N  Z  A    el  o     mass  unc binding unc     B  beta  unc    atomic_mass   unc

	bool initial_state = true;
	char cc[3]={0};
	int NZ=0;
	int N=0;
	int Z=0;
	int A=0;
	char space1=0;
	char el[3]={0};
	char orig[5]={0};
	char space2=0;
	double mass=0;
	char space3=0;
	double mass_unc=0;
	char space4=0;
	double binding=0;
	char space5=0;
	double binding_unc=0;
	char space6=0;
	char B[3]={0};
	//double beta=0;
	//char space7=0;
	//double beta_unc=0;
	//char space8=0;
	char betainfo[30];
	int headdig=0;
	double atomic_mass=0;
	char space9=0;
	double atomic_mass_unc=0;
	char space10=0;


while(!feof(fp)){

	if(!initial_state){
 		/*printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
 		printf("%2s%20s ",B,betainfo);
 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
 		if(space9=='#'){printf("%c\n",space10);}
 		else printf("\n");*/
 		
 	
 	
		cc[0]=0;cc[1]=0;cc[2]=0;
		NZ=0;
		N=0;
		Z=0;
		A=0;
		space1=0;
		el[0]=0;el[1]=0;el[2]=0;
		orig[0]=0;orig[1]=0;orig[2]=0;orig[3]=0;orig[4]=0;
		space2=0;
		mass=0;
		space3=0;
		mass_unc=0;
		space4=0;
		binding=0;
		space5=0;
		binding_unc=0;
		space6=0;
		B[0]=0;B[1]=0;B[2]=0;
		//beta=0;
		//space7=0;
		//beta_unc=0;
		//space8=0;
		for(int j=0;j<30;j++){betainfo[j] = 0;}
		headdig=0;
		atomic_mass=0;
		space9=0;
		atomic_mass_unc=0;
		space10=0;

	}

	fscanstatus=fscanf(fp,"%2c%2d %4d %4d %4d%1c%2c%1c%4c ",cc,&NZ,&N,&Z,&A,&space1,el,&space2,orig);  // %2c%2d => input 2 grid char and 2 grid int number without interval
	fscanstatus=fscanf(fp,"%14lf%1c %11lf%1c %11lf%1c %10lf%1c ",&mass,&space3,&mass_unc,&space4,&binding,&space5,&binding_unc,&space6);
	//fscanf(fp,"%2c %lf%1c %lf%1c ",B,&beta,&space7,&beta_unc,&space8);
	fscanstatus=fscanf(fp,"%2c%24c ",B,betainfo);
	fscanstatus=fscanf(fp,"%3d %13lf%1c %12lf",&headdig,&atomic_mass,&space9,&atomic_mass_unc);
	if(space9=='#')fscanstatus=fscanf(fp,"%1c",&space10);   // !!!! actually there is a '\n' at the end of each line that is why using %2c%2d to fix something wrong

	initial_state=false;

	if(A==Anum && strncmp(el,element,2)==0){
		if(verbose){
			printf("1N-Z \t N \t Z \t A \t EL \t Orig \t MASS EXCESS(keV) \t BINDING ENERGY/A (keV) \t  BETA-DECAY ENERGY(keV) \t ATOMIC MASS(micro-u)\n");
	 		printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
	 		printf("%14.6lf%c %11.6lf%c %13.5lf%c %11.5lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
	 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
	 		printf("%2s%20s ",B,betainfo);
	 		printf("%3d %013.6f%c %11.6f",headdig,atomic_mass,space9,atomic_mass_unc);
	 		if(space9=='#'){printf("%c\n",space10);}
 			else printf("\n");
 		}
 		massreturn[0]=headdig*1e6+atomic_mass;
 		massreturn[1]=atomic_mass_unc;
 		return massreturn;

	}

}	



	cout<<"no match Isotope, make sure 1st letter is CAPITAL	like: He "<<endl;
 	fclose(fp);
 	return massreturn;

}
#endif


Double_t* SearchMolMass(string formula="12C;1 16Cl;4", bool verbose=false){
	   string iso[100];  // iostop name 12C
	  int weight[100];  // number of this iostop in a molecule
	  int compon = 0;

	  istringstream is(formula);
	  char buffer[100];
	  char *p;
	  while(is>>buffer){
		    p = strtok(buffer,";");
		    iso[compon] = p;
		    p = strtok(NULL,";");
		    if(p) weight[compon] = atoi(p);
		    else  weight[compon] = 1;
		    compon++;
  	   }

  	   static double massreturn[2]; // [0] for mass; [1] for err
  	   massreturn[0] = 0;
  	   massreturn[1] =0;

  	   for(int icompon=0; icompon<compon; icompon++){
    		string sA, El;
    		for(unsigned long int i=0; i<iso[icompon].size(); i++){
	     		 char c = iso[icompon][i];
	      		if(isdigit(c)) sA+=c;
	      		else El+=c;
	    	}
	    		if(El.size() == 1) El+=" ";

	    		int A=atoi(sA.c_str());

	    		double* getmass = SearchAME(A,El.c_str());

	    		getmass[0]*=weight[icompon];
	    		getmass[1]*=weight[icompon];

	    		massreturn[0]+=getmass[0];
	    		massreturn[1] = massreturn[1]*massreturn[1] + getmass[1]*getmass[1];
	    		massreturn[1] = TMath::Sqrt(massreturn[1]);
	    }

	    if(verbose)printf("mass= %.4f(%.4f)\n",massreturn[0],massreturn[1]);

	    return massreturn;

}


#ifdef AME2016
Double_t* SearchEleByMass(double mass_input,const char* filename){  // return an array finally, massreturn[0]=> mass; massreturn[1]=>err

	FILE *fp = NULL;
	static Double_t massreturn[2];
	massreturn[0]=-1;
	massreturn[1]=-1;
	double mass_dev=0;
	double mass_abs_dev=0;

	int fscanstatus=0;

	struct isotope{

		char name[5];
		int Anum;
		double mass_dev;
		double mass_abs_dev;
		double mass;
		double mass_err;
		isotope(){
			for(int i=0;i<5;i++){
				name[i]='\0';	
			}
			Anum=0;
			mass_dev=1e10;
			mass_abs_dev=1e10;
			mass = 0;
			mass_err=0;
		}
	};

	isotope iso[4];  	// iso[3] is used to temporary store new candidate form AME2016.txt; then sort => only output 3 closest ones !!!!!!
	isotope isotem;    // temporary container for iso exchange

	fp = fopen(filename,"r");
	if(fp==NULL){ string txfile = TXT;cout<< "fail to open file AME2016.txt to search isotope by mass"<<endl;return massreturn; }

	//a1,i3, i5,i5,i5,1x,a3,a4,1x,f13.5,f11.5,f11.3,     f9.3,1x,a2,f11.3,f9.3,1x,i3,1x,f12.5,f11.5
	//cc NZ  N  Z  A    el  o     mass  unc binding unc     B  beta  unc    atomic_mass   unc

	bool initial_state = true;
	char cc[3]={0};
	int NZ=0;
	int N=0;
	int Z=0;
	int A=0;
	char space1=0;
	char el[3]={0};
	char orig[5]={0};
	char space2=0;
	double mass=0;
	char space3=0;
	double mass_unc=0;
	char space4=0;
	double binding=0;
	char space5=0;
	double binding_unc=0;
	char space6=0;
	char B[3]={0};
	//double beta=0;
	//char space7=0;
	//double beta_unc=0;
	//char space8=0;
	char betainfo[30];
	int headdig=0;
	double atomic_mass=0;
	char space9=0;
	double atomic_mass_unc=0;
	char space10=0;



	while(!feof(fp)){

	if(!initial_state){  // clear container
 		/*printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
 		printf("%2s%20s ",B,betainfo);
 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
 		if(space9=='#'){printf("%c\n",space10);}
 		else printf("\n");*/
 		
 		for(int i=0;i<5;i++){
 			iso[3].name[i]='\0';
 		}
		iso[3].Anum=0;
		iso[3].mass_dev=1e10;
		iso[3].mass_abs_dev=1e10;
		iso[3].mass = 0;
		iso[3].mass_err=0;
 	
 	
		cc[0]=0;cc[1]=0;cc[2]=0;
		NZ=0;
		N=0;
		Z=0;
		A=0;
		space1=0;
		el[0]=0;el[1]=0;el[2]=0;
		orig[0]=0;orig[1]=0;orig[2]=0;orig[3]=0;orig[4]=0;
		space2=0;
		mass=0;
		space3=0;
		mass_unc=0;
		space4=0;
		binding=0;
		space5=0;
		binding_unc=0;
		space6=0;
		B[0]=0;B[1]=0;B[2]=0;
		//beta=0;
		//space7=0;
		//beta_unc=0;
		//space8=0;
		for(int j=0;j<30;j++){betainfo[j] = 0;}
		headdig=0;
		atomic_mass=0;
		space9=0;
		atomic_mass_unc=0;
		space10=0;

	}

	fscanstatus=fscanf(fp,"%2c%2d %4d %4d %4d%1c%2c%1c%4c ",cc,&NZ,&N,&Z,&A,&space1,el,&space2,orig);  // %2c%2d => input 2 grid char and 2 grid int number without interval
	fscanstatus=fscanf(fp,"%13lf%1c %11lf%1c %11lf%1c %9lf%1c ",&mass,&space3,&mass_unc,&space4,&binding,&space5,&binding_unc,&space6);
	//fscanf(fp,"%2c %lf%1c %lf%1c ",B,&beta,&space7,&beta_unc,&space8);
	fscanstatus=fscanf(fp,"%2c%20c ",B,betainfo);
	fscanstatus=fscanf(fp,"%3d %12lf%1c %11lf",&headdig,&atomic_mass,&space9,&atomic_mass_unc);
	if(space9=='#')fscanstatus=fscanf(fp,"%1c",&space10);   // !!!! actually there is a '\n' at the end of each line that is why using %2c%2d to fix something wrong

	initial_state=false;

		/*	if(A==Anum && strncmp(el,element,2)==0){
		printf("1N-Z \t N \t Z \t A \t EL \t Orig \t MASS EXCESS(keV) \t BINDING ENERGY/A (keV) \t  BETA-DECAY ENERGY(keV) \t ATOMIC MASS(micro-u)\n");
 		printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
 		printf("%2s%20s ",B,betainfo);
 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
 		if(space9=='#'){printf("%c\n",space10);}
 		else printf("\n");
		*/
 		massreturn[0]=headdig*1e6+atomic_mass;
 		massreturn[1]=atomic_mass_unc;
 		mass_abs_dev = TMath::Abs(massreturn[0]-mass_input);
 		mass_dev = mass_input - massreturn[0];

		iso[3].name[0]=el[0];
		iso[3].name[1]=el[1];
		iso[3].Anum=A;
		iso[3].mass_dev=mass_dev;
		iso[3].mass_abs_dev = mass_abs_dev;
		iso[3].mass = massreturn[0];
		iso[3].mass_err=massreturn[1];

		for(int i=0;i<4;i++){
			for(int j=i+1;j<4;j++){
				if(iso[j].mass_abs_dev<iso[i].mass_abs_dev){
						memcpy(&isotem, &iso[i], sizeof(struct isotope));
						memcpy(&iso[i], &iso[j], sizeof(struct isotope));
						memcpy(&iso[j], &isotem, sizeof(struct isotope));
				}
			}
		}



	}	// end of while loop


	cout<<"Possible candidates are: "<<endl;
	printf("%d%s\t%.5f(%.5f), [X-candidate] Dev=>  %.5f [micro-u]\n",iso[0].Anum,iso[0].name,iso[0].mass,iso[0].mass_err,iso[0].mass_dev);
	printf("%d%s\t%.5f(%.5f), [X-candidate] Dev=>  %.5f [micro-u]\n",iso[1].Anum,iso[1].name,iso[1].mass,iso[1].mass_err,iso[1].mass_dev);
	printf("%d%s\t%.5f(%.5f), [X-candidate] Dev=>  %.5f [micro-u]\n",iso[2].Anum,iso[2].name,iso[2].mass,iso[2].mass_err,iso[2].mass_dev);

	massreturn[0] = iso[0].mass;
	massreturn[1] = iso[0].mass_err;

 	fclose(fp);
 	return massreturn;

}
#endif


#ifdef AME2020
Double_t* SearchEleByMass(double mass_input,const char* filename){  // return an array finally, massreturn[0]=> mass; massreturn[1]=>err

	FILE *fp = NULL;
	static Double_t massreturn[2];
	massreturn[0]=-1;
	massreturn[1]=-1;
	double mass_dev=0;
	double mass_abs_dev=0;

	int fscanstatus=0;

	struct isotope{

		char name[5];
		int Anum;
		double mass_dev;
		double mass_abs_dev;
		double mass;
		double mass_err;
		isotope(){
			for(int i=0;i<5;i++){
				name[i]='\0';	
			}
			Anum=0;
			mass_dev=1e10;
			mass_abs_dev=1e10;
			mass = 0;
			mass_err=0;
		}
	};

	isotope iso[4];  	// iso[3] is used to temporary store new candidate form AME2016.txt; then sort => only output 3 closest ones !!!!!!
	isotope isotem;    // temporary container for iso exchange

	fp = fopen(filename,"r");
	if(fp==NULL){ cout<< "fail to open file AME2020.txt to search isotope by mass"<<endl;return massreturn; }

	//a1,i3, i5,i5,i5,1x,a3,a4,1x,f13.5,f11.5,f11.3,     f9.3,1x,a2,f11.3,f9.3,1x,i3,1x,f12.5,f11.5
	//cc NZ  N  Z  A    el  o     mass  unc binding unc     B  beta  unc    atomic_mass   unc

	bool initial_state = true;
	char cc[3]={0};
	int NZ=0;
	int N=0;
	int Z=0;
	int A=0;
	char space1=0;
	char el[3]={0};
	char orig[5]={0};
	char space2=0;
	double mass=0;
	char space3=0;
	double mass_unc=0;
	char space4=0;
	double binding=0;
	char space5=0;
	double binding_unc=0;
	char space6=0;
	char B[3]={0};
	//double beta=0;
	//char space7=0;
	//double beta_unc=0;
	//char space8=0;
	char betainfo[30];
	int headdig=0;
	double atomic_mass=0;
	char space9=0;
	double atomic_mass_unc=0;
	char space10=0;



	while(!feof(fp)){

	if(!initial_state){  // clear container
 		/*printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
 		printf("%2s%20s ",B,betainfo);
 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
 		if(space9=='#'){printf("%c\n",space10);}
 		else printf("\n");*/
 		
 		for(int i=0;i<5;i++){
 			iso[3].name[i]='\0';
 		}
		iso[3].Anum=0;
		iso[3].mass_dev=1e10;
		iso[3].mass_abs_dev=1e10;
		iso[3].mass = 0;
		iso[3].mass_err=0;
 	
 	
		cc[0]=0;cc[1]=0;cc[2]=0;
		NZ=0;
		N=0;
		Z=0;
		A=0;
		space1=0;
		el[0]=0;el[1]=0;el[2]=0;
		orig[0]=0;orig[1]=0;orig[2]=0;orig[3]=0;orig[4]=0;
		space2=0;
		mass=0;
		space3=0;
		mass_unc=0;
		space4=0;
		binding=0;
		space5=0;
		binding_unc=0;
		space6=0;
		B[0]=0;B[1]=0;B[2]=0;
		//beta=0;
		//space7=0;
		//beta_unc=0;
		//space8=0;
		for(int j=0;j<30;j++){betainfo[j] = 0;}
		headdig=0;
		atomic_mass=0;
		space9=0;
		atomic_mass_unc=0;
		space10=0;

	}

	fscanstatus=fscanf(fp,"%2c%2d %4d %4d %4d%1c%2c%1c%4c ",cc,&NZ,&N,&Z,&A,&space1,el,&space2,orig);  // %2c%2d => input 2 grid char and 2 grid int number without interval
	fscanstatus=fscanf(fp,"%14lf%1c %11lf%1c %11lf%1c %10lf%1c ",&mass,&space3,&mass_unc,&space4,&binding,&space5,&binding_unc,&space6);
	//fscanf(fp,"%2c %lf%1c %lf%1c ",B,&beta,&space7,&beta_unc,&space8);
	fscanstatus=fscanf(fp,"%2c%24c ",B,betainfo);
	fscanstatus=fscanf(fp,"%3d %13lf%1c %12lf",&headdig,&atomic_mass,&space9,&atomic_mass_unc);
	if(space9=='#')fscanstatus=fscanf(fp,"%1c",&space10);   // !!!! actually there is a '\n' at the end of each line that is why using %2c%2d to fix something wrong

	initial_state=false;

		/*	if(A==Anum && strncmp(el,element,2)==0){
		printf("1N-Z \t N \t Z \t A \t EL \t Orig \t MASS EXCESS(keV) \t BINDING ENERGY/A (keV) \t  BETA-DECAY ENERGY(keV) \t ATOMIC MASS(micro-u)\n");
 		printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
 		printf("%2s%20s ",B,betainfo);
 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
 		if(space9=='#'){printf("%c\n",space10);}
 		else printf("\n");
		*/
 		massreturn[0]=headdig*1e6+atomic_mass;
 		massreturn[1]=atomic_mass_unc;
 		mass_abs_dev = TMath::Abs(massreturn[0]-mass_input);
 		mass_dev = mass_input - massreturn[0];

		iso[3].name[0]=el[0];
		iso[3].name[1]=el[1];
		iso[3].Anum=A;
		iso[3].mass_dev=mass_dev;
		iso[3].mass_abs_dev = mass_abs_dev;
		iso[3].mass = massreturn[0];
		iso[3].mass_err=massreturn[1];

		for(int i=0;i<4;i++){
			for(int j=i+1;j<4;j++){
				if(iso[j].mass_abs_dev<iso[i].mass_abs_dev){
						memcpy(&isotem, &iso[i], sizeof(struct isotope));
						memcpy(&iso[i], &iso[j], sizeof(struct isotope));
						memcpy(&iso[j], &isotem, sizeof(struct isotope));
				}
			}
		}



	}	// end of while loop


	cout<<"Possible candidates are: "<<endl;
	printf("%d%s\t%.5f(%.5f), [X-candidate] Dev=>  %.5f [micro-u]\n",iso[0].Anum,iso[0].name,iso[0].mass,iso[0].mass_err,iso[0].mass_dev);
	printf("%d%s\t%.5f(%.5f), [X-candidate] Dev=>  %.5f [micro-u]\n",iso[1].Anum,iso[1].name,iso[1].mass,iso[1].mass_err,iso[1].mass_dev);
	printf("%d%s\t%.5f(%.5f), [X-candidate] Dev=>  %.5f [micro-u]\n",iso[2].Anum,iso[2].name,iso[2].mass,iso[2].mass_err,iso[2].mass_dev);

	massreturn[0] = iso[0].mass;
	massreturn[1] = iso[0].mass_err;

 	fclose(fp);
 	return massreturn;

}
#endif



int ReadParaLst(string PATH,string filename,double* paracontain){// return 1 for seccess; -1 for fail;  paracontain: eje0, eje1, timeRef, massRef, bref, qx
	
  const int kMaxBufLen = 0x800000;
  int fscanstatus = 0;

  string inputfile = PATH +"../LST/" + filename + ".lst";//"LST/"
  //ifstream fin_lst;
  FILE *fin_lst = NULL;
  //fin_lst.open(inputfile.c_str());
  fin_lst = fopen(inputfile.c_str(),"r");
	if(fin_lst==NULL){ cerr<<"\e[1;33m"<< "fail to open file .LST for parameters loading !!!"<<"\e[0m"<<endl;return -1; }

  char *buffer = new char[kMaxBufLen];


  while(!feof(fin_lst)){
		char* fget_status = fgets(buffer,kMaxBufLen,fin_lst);
    	if(feof(fin_lst)){
      		cerr<<"cannot find data!!!"<<endl;
	      	return -1;
		}
		if(strncmp(buffer,"cmline1",7)==0){  // strncmp: compare buffer with "cmline1" for 7 character, if equal will be 0; buffer> cmline1 -> value>0
                                                    // buffer< cmline1  -> value<0; based on ACSII sequence
      		cout<<"\e[1;33m reading Parameters from "<<inputfile<<"\e[0m"<<endl;
	      	break;
    	}
   }


	//cout<<buffer<<endl;
	double timeA=0;
	double timeB=0;
	char cmline[50]; for(int i=0;i<50;i++){cmline[i] = '\0';}
	char comma;
	double timeRef=0;
	double massRef=0;
	double timeX=0;
	double massX=0;


	fscanstatus=fscanf(fin_lst,"%12c:%lf,%lf,%lf%2c",cmline,&timeA,&timeB,&timeA,cmline);
//cout<<cmline<<"+"<<timeA<<"+"<<timeB<<"+"<<cmline<<endl;
//return 0;
	if(cmline[1]=='R'){
		fscanstatus=fscanf(fin_lst,":%lf,%lf",&timeRef,&massRef);
		cout<<"************** Set Eje time; Ref.mass and Ref.tof for spectrum preview automatically****************"<<endl;
		printf("Ejection time candidates: %.3f ;  %.3f\n",timeA,timeB);
		printf("Reference ion TOF: %.3f ; Mass: %.6f\n",timeRef,massRef);
		cout<<"***************************************************************************************"<<endl;

		paracontain[0] = timeA;
		paracontain[1] = timeB;
		paracontain[2]= timeRef;
		paracontain[3]=massRef;

		return 1;
	}

	fscanstatus=fscanf(fin_lst,"%12c:%lf,%lf,T:%lf,%lf",cmline,&timeRef,&massRef,&timeX,&massX); // use ",T:" in .lst as interval symbol

	if(timeRef ==0 || massRef ==0 ){
		cout<<"************** Only Eje time available for preview  ****************"<<endl;
		printf("Ejection time candidates: %.3f ;  %.3f\n",timeA,timeB);
		cout<<"Reference ion TOF: N/A ; Mass: N/A"<<endl;
		cout<<"***************************************************************************************"<<endl;
		paracontain[0] = timeA;
		paracontain[1] = timeB;

		return 1;

	}

	// calculate the possible charge of X ion
	auto ChargeX = [](double time_r,double mass_r,double time_x, double mass_x)->int {
		double moq_x = TMath::Power(time_x/time_r,2) * mass_r;
		double different=10000;
		int qx =0;
		for(int i=1;i<6;i++){
			if(TMath::Abs(mass_x/i-moq_x) < different){
				different = TMath::Abs(mass_x/i-moq_x);
				qx=i;
			}
		}
		return qx;
	};

	double ar=0,br=0,ax=0,bx=0;
	fscanstatus=fscanf(fin_lst,"%32c: %lf, %lf, %lf, %lf",cmline,&ar,&br,&ax,&bx);

	cout<<"*********** Successfully set parameter for spectrum preview automatically**************"<<endl;
	printf("Ejection time candidates: %.3f ;  %.3f\n",timeA,timeB);
	printf("Reference ion TOF: %.3f ; Mass: %.6f\n",timeRef,massRef);
	printf("Target ion TOF: %.3f ; Mass: %.6f\n",timeX,massX);
	printf("Charge: q_x = %d\n",ChargeX(timeRef,massRef,timeX,massX));
	if(ar!=0)	printf("ar: %.5f ; br: %.5f\nax: %.5f ; bx: %.5f\n",ar,br,ax,bx);
	cout<<"***************************************************************************************"<<endl;

	// paracontain=> time eje0,eje1; time TOF ref ; mass Ref
	if(timeX<timeRef){
		paracontain[0] =TMath::Min(timeA,timeB); 
		paracontain[1]=TMath::Max(timeA,timeB); 
	}
	else{
		paracontain[0] =TMath::Max(timeA,timeB); 
		paracontain[1]=TMath::Min(timeA,timeB);
	}

	
	paracontain[2]= timeRef;
	paracontain[3]=massRef;
	paracontain[4]=br;   // br in [us]
	paracontain[5] = ChargeX(timeRef,massRef,timeX,massX);

	fclose(fin_lst);

	return 1;

}



#endif
