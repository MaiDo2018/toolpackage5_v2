

void batchfitagain(int NumofIon,const char* RefName){ // how many ions record in history file; and which ion peak to determine shape of fit function
	int IndexOfRef=-1;
	for(int i=0;i<NumofIon;i++){
		if(strcmp(FHistory[i]->name,RefName)==0){ IndexOfRef=i; break;}
	}

	if(IndexOfRef==-1){cout<<"fail to get "<<RefName<<" !! check again!!"<<endl;}

	fithistory(FHistory[IndexOfRef]);
	fitlooping("LME",1,fitrangeL,fitrangeR,false);
	fitlooping("LMEQ",50,fitrangeL,fitrangeR,false);

	int numtailpar = SetOneorTwo +1;  // number of pars of the fit function for Ref
	double *shappar = new double[numtailpar];
	int peak_para_offset = Paras_offset_cal(MainPeakIndex-1,SetOneorTwo,MainPeakIndex,ParasLock);

	for(int i=0;i<=SetOneorTwo;i++){
		shappar[i] = tem_func->GetParameter(peak_para_offset+2+i);
		cout<<"par"<<i<<" = "<<shappar[i]<<endl;
	}

	cout<<"Ref fitting finish. continue?? [y/n]"<<endl;
	char keepgoing;
	do{
		cin>>keepgoing;
	}while(keepgoing != 'y' && keepgoing!= 'n');

	if(keepgoing !='y') {cout<<"break!!"<<endl; return;}

	for(int i=0;i<NumofIon;i++){

		fithistory(FHistory[i]);

		// fix parameter shared by Ref fit function
		if(ParasLock){// only mainpeak has parameters for sigma and tails
			peak_para_offset = Paras_offset_cal(MainPeakIndex-1,SetOneorTwo,MainPeakIndex,ParasLock);
			for(int j=0;j<=SetOneorTwo && j<numtailpar ;j++){
				tem_func->FixParameter(peak_para_offset+2+j,shappar[j]);
			}
		}
		else{
			for(int PeakIndex=0;PeakIndex<NumOfPeaks;PeakIndex++){
				peak_para_offset = Paras_offset_cal(PeakIndex,SetOneorTwo,MainPeakIndex,ParasLock);
				for(int j=0;j<=SetOneorTwo && j<numtailpar ;j++){
					tem_func->FixParameter(peak_para_offset+2+j,shappar[j]);
				}
			}
		}

		fitlooping("LME",1,fitrangeL,fitrangeR,true);
		fitlooping("LMEQ",50,fitrangeL,fitrangeR,true);
		cout<<"this is "<<FHistory[i]->name<<endl;
		cout<<endl;
		while(c1->WaitPrimitive()){;}

	}



}