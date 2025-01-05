



void LoadToCurrentBetaTree(TTree* tread3, string _ionname, int _ionIndex, double tof_main, double tof_i){
		//_ionname ; _ionIndex ==> name and index of ion will be extracted from history tbeta_tof file
		//

		Long64_t mr_Bingo_sweeps_global=0; double mr_Bingo_time=0; Long64_t mr_Bingo_time_gclock=0;
		Long64_t time1=0; Long64_t time1_relative=0; Long64_t sweeps_gclock_1=0; int adc1=0; double EkeV_1=0;
		Long64_t time2=0; Long64_t time2_relative=0;  Long64_t sweeps_gclock_2=0; int adc2=0; double EkeV_2=0;
		int delta_t=0;
		Long64_t delta_t_decay=0;
		string* IonName_ptr=NULL;
		int Ionindex =0;

		tbeta_tof->SetBranchAddress("sweeps_global_mr",&mr_Bingo_sweeps_global); //tbeta_tof from global
		tbeta_tof->SetBranchAddress("time_mr",&mr_Bingo_time);
		tbeta_tof->SetBranchAddress("time_gclock_mr",&mr_Bingo_time_gclock);
		tbeta_tof->SetBranchAddress("sweeps_gclock_1",&sweeps_gclock_1);
		tbeta_tof->SetBranchAddress("time1",&time1);
		tbeta_tof->SetBranchAddress("time1_relative",&time1_relative);
		tbeta_tof->SetBranchAddress("adc1",&adc1);
		tbeta_tof->SetBranchAddress("EkeV_1",&EkeV_1);
		tbeta_tof->SetBranchAddress("sweeps_gclock_2",&sweeps_gclock_2);
		tbeta_tof->SetBranchAddress("time2",&time2);
		tbeta_tof->SetBranchAddress("time2_relative",&time2_relative);
		tbeta_tof->SetBranchAddress("adc2",&adc2);
		tbeta_tof->SetBranchAddress("EkeV_2",&EkeV_2);
		tbeta_tof->SetBranchAddress("delta_t",&delta_t);
		tbeta_tof->SetBranchAddress("delta_t_decay",&delta_t_decay);

		tread3->SetBranchAddress("sweeps_global_mr",&mr_Bingo_sweeps_global);
		tread3->SetBranchAddress("time_mr",&mr_Bingo_time);
		tread3->SetBranchAddress("time_gclock_mr",&mr_Bingo_time_gclock);
		tread3->SetBranchAddress("sweeps_gclock_1",&sweeps_gclock_1);
		tread3->SetBranchAddress("time1",&time1);
		tread3->SetBranchAddress("time1_relative",&time1_relative);
		tread3->SetBranchAddress("adc1",&adc1);
		tread3->SetBranchAddress("EkeV_1",&EkeV_1);
		tread3->SetBranchAddress("sweeps_gclock_2",&sweeps_gclock_2);
		tread3->SetBranchAddress("time2",&time2);
		tread3->SetBranchAddress("time2_relative",&time2_relative);
		tread3->SetBranchAddress("adc2",&adc2);
		tread3->SetBranchAddress("EkeV_2",&EkeV_2);
		tread3->SetBranchAddress("delta_t",&delta_t);
		tread3->SetBranchAddress("delta_t_decay",&delta_t_decay);
		tread3->SetBranchAddress("IonName",&IonName_ptr);
		tread3->SetBranchAddress("IonFillIndex",&Ionindex);

		f_beta_tof_tree->cd(); // from global 

		for(int i=0; i<(int)tread3->GetEntries();i++){
			tread3->GetEntry(i);

			if(_ionname != (*IonName_ptr) || _ionIndex != Ionindex) continue;

			mr_Bingo_time/=(tof_i/tof_main);
			tbeta_tof->Fill();
		}

		tbeta_tof->ResetBranchAddresses();

}


void CombineBeta(string PATH = "./to the control file.txt", string RunPath = "../parent path of folder with Runxx"){
	///RunPath =home/xian/winfoder/MT2021/ 

	if(gSystem->AccessPathName(PATH.c_str()) == true){
				cout<<"\e[1;31m"<<"Error!!!! file of beta loading list not exist!!!!"<<"\e[0m"<<endl;
				return;
	}

	FILE* f_b_list = NULL;

	f_b_list = fopen(PATH.c_str(),"r");
	char buffer[1000];
	fgets(buffer,1000,f_b_list);

	vector<string> folderlist;
	vector<int> runlist;
	char _ionname[100];
	vector<double> toflist;
	vector<int> indexlist;
	vector<double> tof_reflist;

	double maintof=0;
	double maintof_ref=0;

	char tem_folder[100];
	int tem_run;
	double tem_tof;
	double tem_tof_err;
	int mainbit;
	int tem_index;
	int tem_count;
	double tem_tof_ref;

	while(!feof(f_b_list)){
		fscanf(f_b_list,"%s\t%d\t%s\t%lf\t%lf\t%d\t%d\t%d\t%lf",tem_folder,&tem_run,_ionname,&tem_tof,&tem_tof_err,&mainbit,&tem_index,&tem_count,tem_tof_ref);
// there is not tem_tof_ref information in some old  control file.txt file
		if(!feof(f_b_list)){
			if(mainbit == 1){maintof = tem_tof;maintof_ref = tem_tof_ref;continue;}
			folderlist.push_back(tem_folder);
			runlist.push_back(tem_run);
			toflist.push_back(tem_tof);
			indexlist.push_back(tem_index);
tof_reflist.push_back(tem_tof_ref);
		}
	}

	fclose(f_b_list);

	for(int i=0;i<(int)folderlist.size();i++){
		//printf("%s\t%d\t%s\t%f\t%d\n",folderlist[i].c_str(),runlist[i],_ionname,toflist[i],indexlist[i]);

		string rootfilepath = RunPath + folderlist[i] + "/Run" + to_string(runlist[i]) +"/rootfiles/beta-tof_tree_raw.root" ;


		if(gSystem->AccessPathName(rootfilepath.c_str()) == true){
				cout<<"\e[1;31m"<<"Error!!!! file of \"beta-tof_tree_raw.root\" not exist!!!!"<<"\e[0m"<<endl;
				continue;
		}

		TFile* f_get_root = new TFile(rootfilepath.c_str(),"READ");

		TTree* gettree = (TTree*)f_get_root->Get("tbeta_tof_refine");

		//cout<<gettree->GetEntries()<<endl;

		LoadToCurrentBetaTree(gettree, _ionname, indexlist[i], maintof_ref, tof_reflist[i]);
	}

}