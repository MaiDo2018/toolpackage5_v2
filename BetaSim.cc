double TaoCorrEff(double tao, double timewindow){

	double lamdaT = timewindow/tao;

	return 1- lamdaT*TMath::Exp(-lamdaT)/(1-TMath::Exp(-lamdaT));
}


double ExpectHalflife(double tao=1, double tUpperLimit=10, int nloop=100000){
		double sum=0;
		int N_counts=0;
		TRandom3 * r = new TRandom3((unsigned int)time(0));
		for(int i=0; i<nloop;i++){
			double getvalue = r->Exp(tao);
			if(getvalue<=tUpperLimit){
				sum+=getvalue;
				N_counts++;
			}			
		}

		double lifetime_m = sum/N_counts;

		double lifetime_up = lifetime_m/(1 - 1/TMath::Sqrt(N_counts));

		double lifetime_down = lifetime_m/(1 + 1/TMath::Sqrt(N_counts));

		double lamdaT = tUpperLimit/tao;

		double expect_tao = tao*(1- lamdaT*TMath::Exp(-lamdaT)/(1-TMath::Exp(-lamdaT)));

		printf("expect tao with window = %f\n",expect_tao);

		printf("down = %f; up= %f\n",lifetime_down,lifetime_up);

		return lifetime_m;//* TMath::Log(2);
		
}




void fillhisto(TH1D* hin, long npoints, long nloops, double window){

		TRandom3* ran = new TRandom3((unsigned int)time(0));
		double  sum=0;

		double * value = new double[npoints];
		for(int i=0;i<npoints;i++){
			value[i]=0;
		}

		for(long index_loop=0;index_loop<nloops;index_loop++){
			sum=0;
			for(long index_point=0;index_point<npoints;index_point++){
					sum+=ran->Uniform(window);
			}

			hin->Fill(sum/npoints);

		}

}



RIHalflive
Halflife2Lifetime(double Inhalflife_ms=150)



vector<double> dk_time_test;
long Nbeta_test = 10;
void GenerateDKtime(vector<double> &_dk_time_test, long _ncounts){
	TRandom3* ran = new TRandom3((unsigned int)time(0));
	double mylifetime = Halflife2Lifetime(RIHalflive/1000000.0)/1000; //==>> lifetime in second

	_dk_time_test.clear();
	_dk_time_test.shrink_to_fit();

	for(long i=0;i<_ncounts;i++){
		_dk_time_test.push_back(ran->Exp(mylifetime));
	}


}