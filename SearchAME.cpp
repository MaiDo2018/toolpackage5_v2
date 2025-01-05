#include "SearchAME.h"

Double_t* SearchAME(const char* filename, int Anum,const char* element ){  // return an array finally, massreturn[0]=> mass; massreturn[1]=>err

	FILE *fp = NULL;
	Double_t *massreturn = new Double_t[2];
	massreturn[0]=-1;
	massreturn[1]=-1;

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
 		else printf("\n");
 		*/
 	
 	
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

	fscanf(fp,"%2c%2d %4d %4d %4d%1c%2c%1c%4c ",cc,&NZ,&N,&Z,&A,&space1,el,&space2,orig);  // %2c%2d => input 2 grid char and 2 grid int number without interval
	fscanf(fp,"%13lf%1c %11lf%1c %11lf%1c %9lf%1c ",&mass,&space3,&mass_unc,&space4,&binding,&space5,&binding_unc,&space6);
	//fscanf(fp,"%2c %lf%1c %lf%1c ",B,&beta,&space7,&beta_unc,&space8);
	fscanf(fp,"%2c%20c ",B,betainfo);
	fscanf(fp,"%3d %12lf%1c %11lf",&headdig,&atomic_mass,&space9,&atomic_mass_unc);
	if(space9=='#')fscanf(fp,"%1c",&space10);   // !!!! actually there is a '\n' at the end of each line that is why using %2c%2d to fix something wrong

	initial_state=false;

	if(A==Anum && strncmp(el,element,2)==0){
		printf("1N-Z \t N \t Z \t A \t EL \t Orig \t MASS EXCESS(keV) \t BINDING ENERGY/A (keV) \t  BETA-DECAY ENERGY(keV) \t ATOMIC MASS(micro-u)\n");
 		printf("%2s %2d %4d %4d %4d %2s %4s ",cc,NZ,N,Z,A,el,orig);  // using %ns output n number of char in array
 		printf("%13.5lf%c %11.5lf%c %11.3lf%c %9.3lf%c ",mass,space3,mass_unc,space4,binding,space5,binding_unc,space6);
 		//printf("%2s %11.3lf%c %9.3lf%c ",B,beta,space7,beta_unc,space8);
 		printf("%2s%20s ",B,betainfo);
 		printf("%3d %012.5f%c %11.5f",headdig,atomic_mass,space9,atomic_mass_unc);
 		if(space9=='#'){printf("%c\n",space10);}
 		else printf("\n");

 		massreturn[0]=headdig*1e6+atomic_mass;
 		massreturn[1]=atomic_mass_unc;
 		return massreturn;

	}

}	



	cout<<"no match,make sure 1st letter is CAPITAL	like: He "<<endl;
 	fclose(fp);
 	return massreturn;

}