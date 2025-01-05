#ifndef _REFINEBETA_H
#define _REFINEBETA_H
#include <iostream>
#include <vector>
#include "TMath.h"
#include <algorithm>


using namespace std;

struct match_event{
	Long64_t B1_sweeps_gclock;
	Long64_t B1_time;
	Long64_t B1_time_relative;
	int B1_adc;
	Long64_t B2_sweeps_gclock;
	Long64_t B2_time;
	Long64_t B2_time_relative;
	int B2_adc;
	int delta_t;
	void clear(){
		B1_sweeps_gclock=-1;
		B1_time=-1;
		B1_time_relative=-1;
		B1_adc=-1;
		B2_sweeps_gclock=-1;
		B2_time=-1;
		B2_time_relative=-1;
		B2_adc=-1;
		delta_t=-1;
	}

	match_event(){
		clear();
	}

};

struct hit_info{  // for getting one event in beta.lst file
	bool empty;   // container empty or not empty=> true; full=>false
	Long64_t sweeps_gclock;
	Long64_t time;  // beta event time. time=0 ns at the start of USB ADC acqusistion
	Long64_t time_relative;
	int channel;
	int adc;
};

void Cij(int i,int j, vector<int> &r,int num,vector<vector<int>>& result);
double Evaluate(Long64_t time1_in , Long64_t time2_in);
void GetOutputList(vector<match_event>& in_subgroup, vector<vector<int>>& in_result_Cij, vector<int>& outlist_index);
void RefineBetaBeta(vector<match_event>& in_subgroup, vector<match_event>& Add_to_final);
void SelectByBeta(vector<match_event>& in_subgroup, vector<Long64_t>& b2_time_list,vector<match_event>& Add_to_final,int num2keep=3);// if Num of b2<= b1
void SelectByTof(vector<match_event>& in_subgroup, vector<Long64_t>& b1_time_list,vector<match_event>& Add_to_final,int num2keep=3);// if Num of b1<b2

// (total num, seleted num, temp result with length==num, num=j, vector with final set of combination for return)
// for example (6,3,resulttemp,3,result)
void Cij(int i,int j, vector<int> &r,int num,vector<vector<int>>& result){
	if(j==1){
		for(int k=0;k<i;k++){
			vector<int>temp(num);
			r[num-1] = k;
			for(int index=0;index<num;index++){
				temp[index] = r[index];
				//cout<<r[index]<<' ';
			}
			result.push_back(temp);
			//cout<<endl;
		}
	}
	else if(j==0){

	}
	else{
			for(int k=i;k>=j;k--){
				r[j-2] = k-1;
				Cij(k-1,j-1,r,num,result);
			}
	}
}

int printindex;
void RefineMatchHit(vector<match_event>& in_Match_hit, vector<match_event>& Match_hit_refine){// selected pair in in_Match_hit will be added to Match_hit_refine
		vector<match_event> match_temp_group; // pick up match event close related
		bool have_repeat=false;
cout<<endl;

		for(Long64_t index_match=0;index_match<(Long64_t)in_Match_hit.size();index_match++){
				have_repeat=false;
//printindex = index_match	;			
//cout<<"print index = "<<printindex<<endl;
				if(index_match%20==0) cout<<"Refine precessing...  "<<(double)index_match/in_Match_hit.size()*100<<"%"<<'\r'<<flush;
				if(match_temp_group.size()==0){match_temp_group.push_back(in_Match_hit[index_match]);continue;} // the first pair beta-beta
				else{
					for(Long64_t index_temp_group=(Long64_t)match_temp_group.size()-1; index_temp_group>=0; index_temp_group--){
						bool B1_equal = ( in_Match_hit[index_match].B1_time == match_temp_group[index_temp_group].B1_time );
						bool B2_equal = ( in_Match_hit[index_match].B2_time == match_temp_group[index_temp_group].B2_time );
						if(B1_equal || B2_equal){
							match_temp_group.push_back(in_Match_hit[index_match]);
							have_repeat=true;
							break;
						}
					}

					if(!have_repeat){ // new unlated event, belong to next group and precess current group
/*
cout<<"do moore"<<endl;
								for(int ijk=0;ijk<(int)match_temp_group.size();ijk++){
			printf("%lld,%lld\n",match_temp_group[ijk].B1_time,match_temp_group[ijk].B2_time);
		}
cout<<endl;*/ 
				
						RefineBetaBeta(match_temp_group,Match_hit_refine);


						match_temp_group.clear(); match_temp_group.shrink_to_fit();

						match_temp_group.push_back(in_Match_hit[index_match]);
					}
					
				}

		}// run precess again
/*
		for(int ijk=0;ijk<(int)match_temp_group.size();ijk++){
			printf("%lld,%lld\n",match_temp_group[ijk].B1_time,match_temp_group[ijk].B2_time);
		}
*/


		RefineBetaBeta(match_temp_group,Match_hit_refine);


}



void RefineBetaBeta(vector<match_event>& in_subgroup, vector<match_event>& Add_to_final){ //(subgroup with repeat data; out to final list)

	if(in_subgroup.size()==0) return;
	if(in_subgroup.size()==1){
		Add_to_final.push_back(in_subgroup[0]);
		return;
	}

	vector<Long64_t> b1_time_list; bool B1_equal=false;
	vector<Long64_t> b2_time_list; bool B2_equal=false;
	vector<Long64_t> b1_time_list_clone;
	vector<Long64_t> b2_time_list_clone;
	bool HaveReduced=false;

	int num2keep_for_each_beta=1;  // how many beta1 caandidate to keep for each beta2 ; or beta2 candidates keep for each beta1

	vector<match_event> in_subgroup_clone;

	for(unsigned long index=0; index<in_subgroup.size();index++){
		in_subgroup_clone.push_back(in_subgroup[index]);
	}



After_reduce:	b1_time_list.clear();b2_time_list.clear();



	b1_time_list.push_back(in_subgroup[0].B1_time);
	b2_time_list.push_back(in_subgroup[0].B2_time);
	


	for(long index_subgroup=1; index_subgroup< (long)in_subgroup.size(); index_subgroup++){
		B1_equal=false;
		B2_equal=false;

		for(long index_b1_list=0; index_b1_list< (long)b1_time_list.size();index_b1_list++){
				if(b1_time_list[index_b1_list] == in_subgroup[index_subgroup].B1_time){
					B1_equal=true;
					break;
				}
		}

		if(!B1_equal)b1_time_list.push_back(in_subgroup[index_subgroup].B1_time); // new B1_time

		for(long index_b2_list=0; index_b2_list< (long)b2_time_list.size();index_b2_list++){
				if(b2_time_list[index_b2_list] == in_subgroup[index_subgroup].B2_time){
					B2_equal=true;
					break;
				}
		}

		if(!B2_equal)b2_time_list.push_back(in_subgroup[index_subgroup].B2_time); // new B2_time

	}


	if(b1_time_list_clone.size()==0){
		for(unsigned long index=0;index<b1_time_list.size();index++){
			b1_time_list_clone.push_back(b1_time_list[index]);
		}
		for(unsigned long index=0;index<b2_time_list.size();index++){
			b2_time_list_clone.push_back(b2_time_list[index]);
		}
	}



	int Max_num_keep = TMath::Min( (int)b1_time_list.size() , (int)b2_time_list.size() );

	int num_subgroup = (int)in_subgroup.size();


//if(printindex>139)cout<<"number in this group= "<<num_subgroup<<"  Max to keep= "<<Max_num_keep<<endl;
/*for(int ij=0; ij<(int)in_subgroup.size();ij++){
	printf("%lld \t %lld\n",in_subgroup[ij].B1_time,in_subgroup[ij].B2_time);
}*/



	if(TMath::Binomial(num_subgroup,Max_num_keep)*Max_num_keep*8*2>50000000){ // 50Mb too many combinaiton, hanle with simple method
		cout<<endl;
cout<<"too much"<<endl;
		if(num2keep_for_each_beta>1){
			if(b2_time_list.size()<=b1_time_list.size()){
					vector<match_event> Bybeta;
					SelectByBeta(in_subgroup,b2_time_list,Bybeta,num2keep_for_each_beta--);
					in_subgroup.clear(); in_subgroup.shrink_to_fit();
					for(unsigned long index=0;index<Bybeta.size();index++){
						in_subgroup.push_back(Bybeta[index]);
					}
				goto After_reduce;
			}
			else{
					vector<match_event> ByTof;
					SelectByTof(in_subgroup,b1_time_list,ByTof,num2keep_for_each_beta--);
					in_subgroup.clear(); in_subgroup.shrink_to_fit();
					for(unsigned long index=0;index<ByTof.size();index++){
						in_subgroup.push_back(ByTof[index]);
					}
				goto After_reduce;
			}

			//cout<<"reduce"<<endl;
			//if(b2_time_list.size()<=b1_time_list.size())

		}
		else{
			if(b2_time_list_clone.size()<=b1_time_list_clone.size()){ cout<<"by b2"<<endl;
				SelectByBeta(in_subgroup_clone,b2_time_list_clone,Add_to_final,1);
				return;
			}
			else{cout<<"by b1"<<endl;
				SelectByTof(in_subgroup_clone,b1_time_list_clone,Add_to_final,1);
				return;
			}
		}
	}




vector<vector<int>> result_Cij_refine; // delete combinations with repeated information


	if(Max_num_keep>1){
		for(;Max_num_keep>=1;Max_num_keep--){
				vector<vector<int>> result_Cij;
				vector<int> resulttemp_Cij(Max_num_keep);
				Cij(num_subgroup,Max_num_keep,resulttemp_Cij,Max_num_keep,result_Cij);  
				for(int index_result=0;index_result<(int)result_Cij.size(); index_result++){
						bool found_repeat=false;
					for(int row_i=0; row_i<Max_num_keep; row_i++){
						for(int row_j=row_i+1; row_j<Max_num_keep; row_j++){
							if(in_subgroup[ result_Cij[index_result][row_i] ].B1_time  == in_subgroup[ result_Cij[index_result][row_j] ].B1_time){found_repeat=true; break;}
							if(in_subgroup[ result_Cij[index_result][row_i] ].B2_time  == in_subgroup[ result_Cij[index_result][row_j] ].B2_time){found_repeat=true; break;}
						}
						if(found_repeat)break;
					}

					if(!found_repeat)result_Cij_refine.push_back( result_Cij[index_result] ); // possible combination

				}// finish get result_Cij_refine

				if(result_Cij_refine.size()>0) break;
		}

	}
	else{// keep only one event in this subgroup

		vector<int> single;
		for(int index=0;index<num_subgroup;index++){
			single.clear();
			single.push_back(index);
			result_Cij_refine.push_back(single);
		}

	} 

/*
for(int ij=0; ij<(int)in_subgroup.size();ij++){
	printf("%lld \t %lld\n",in_subgroup[ij].B1_time,in_subgroup[ij].B2_time);
}*/

/*
	for(int irow=0;irow<(int)result_Cij_refine.size();irow++){
		cout<<"[";
		for(int icolumn=0;icolumn<Max_num_keep;icolumn++){
			printf("%d,",result_Cij_refine[irow][icolumn]);
		}
		cout<<"]"<<endl;
	}*/

	if(result_Cij_refine.size()==0) return;

	
	vector<int> outlist_index; // store the index of event in subgroup which should be kept



	GetOutputList(in_subgroup,result_Cij_refine,outlist_index);

	sort(outlist_index.begin(),outlist_index.end());  // sort to have a sequence "" early "" to ""later""

	for(int index=0;index< Max_num_keep; index++){
		Add_to_final.push_back(in_subgroup[ outlist_index[index] ]);
	}

}


Long64_t b2tob1_standard = -200;  // standard value of b2_time - b1_time depend on shaping time if shape t1 > shape t2 ==> t2-t1 <0
	

double Evaluate(Long64_t time1_in , Long64_t time2_in){ // for beta-beta: prefer short time difference
	//return (double) TMath::Abs(time1_in-time2_in);   // select the closest b1 and b2;
	return (double) TMath::Abs(time2_in - time1_in - b2tob1_standard);
}


void GetOutputList(vector<match_event>& in_subgroup, vector<vector<int>>& in_result_Cij, vector<int>& outlist_index){
			vector<double> list_evaluate_value; 
			for(Long64_t row_index=0; row_index<(Long64_t)in_result_Cij.size();row_index++){
				
				double temp_row_evalue_value=0;

				for(int column_index=0; column_index<(int)(in_result_Cij[0].size()); column_index++){
				
						int index_look_at_subgroup = in_result_Cij[row_index][column_index];
						temp_row_evalue_value+= Evaluate(in_subgroup[index_look_at_subgroup].B1_time, in_subgroup[index_look_at_subgroup].B2_time);
				}// get final evaluate value for each row (ONE possible combination)

				list_evaluate_value.push_back(temp_row_evalue_value);
			}// fill whole list of evaluate value of all combination

			//select to choose the most possible combination and return their index

			double* list_copy = new double[ list_evaluate_value.size() ];
			for(Long64_t index=0; index<(Long64_t)list_evaluate_value.size(); index++){list_copy[index]=list_evaluate_value[index]; }

			Long64_t* index_of_sort_list = new Long64_t[ list_evaluate_value.size() ];

			TMath::Sort((Long64_t)list_evaluate_value.size(), list_copy, index_of_sort_list,kFALSE); // first one with the index of "smallest"!! value in evaluate list

			outlist_index=in_result_Cij[ index_of_sort_list[0] ];

			delete[] list_copy;
			delete[] index_of_sort_list;

}



void SelectByBeta(vector<match_event>& in_subgroup, vector<Long64_t>& b2_time_list,vector<match_event>& Add_to_final, int num2keep){// b2_time_list sorted to "later" to "earlier"

	int Reverse_direct=0;
	
	vector<vector<Long64_t>> B1_candidate_group;
	vector<Long64_t> B1_candidate_each_B2;

	vector<match_event> selected_event; // earlier to later
	vector<match_event> selected_event_reverse; // later to earlier

Rever:
	B1_candidate_group.clear(); B1_candidate_group.shrink_to_fit();

	if(Reverse_direct==1){
		for(unsigned long index_sele=0; index_sele<selected_event.size();index_sele++){
			selected_event_reverse.push_back(selected_event[index_sele]);
		}
			selected_event.clear(); selected_event.shrink_to_fit();
	}

	if(Reverse_direct==0){
		sort(b2_time_list.rbegin(),b2_time_list.rend()); // the first round is later to earlier
	}
	else{
		sort(b2_time_list.begin(),b2_time_list.end());// the second round is earlier to later
	}

	for(int b2_list_index=0;b2_list_index<(int)(b2_time_list.size());b2_list_index++){// get all B1 related to each B2
		B1_candidate_each_B2.clear();
		for(long index_subgroup=0; index_subgroup< (long)in_subgroup.size(); index_subgroup++){
			if(in_subgroup[index_subgroup].B2_time == b2_time_list[b2_list_index]) B1_candidate_each_B2.push_back(in_subgroup[index_subgroup].B1_time);
		}

		if(Reverse_direct==0)sort(B1_candidate_each_B2.rbegin(),B1_candidate_each_B2.rend()); // "later" to "earlier"
		else{sort(B1_candidate_each_B2.begin(),B1_candidate_each_B2.end());}

		B1_candidate_group.push_back(B1_candidate_each_B2);


	}

	match_event tem_event;	
	vector<double> time_difference;

	for(int b2_list_index=0;b2_list_index<(int)(b2_time_list.size());b2_list_index++){
		time_difference.clear(); 
		time_difference.shrink_to_fit();
		for(int index_column_each_B2=0;index_column_each_B2<(int)B1_candidate_group[b2_list_index].size();index_column_each_B2++){

			time_difference.push_back(Evaluate(B1_candidate_group[b2_list_index][index_column_each_B2],b2_time_list[b2_list_index]));
		}

		double* time_difference_clone = new double[ time_difference.size() ];
		Long64_t* index_of_sort_list = new Long64_t[ time_difference.size() ];
		for(unsigned long in_cp=0;in_cp<time_difference.size();in_cp++){time_difference_clone[in_cp] = time_difference[in_cp];}

		TMath::Sort((Long64_t)time_difference.size(), time_difference_clone, index_of_sort_list,kFALSE); // first one with the index of "smallest"!! value in evaluate list

		if(num2keep>1){
			for(unsigned long index_keep=0;index_keep<time_difference.size() && index_keep< (unsigned int)num2keep; index_keep++){// smallest B1-B2
					tem_event.B1_time =	B1_candidate_group[b2_list_index][index_of_sort_list[index_keep]];
					tem_event.B2_time = b2_time_list[b2_list_index];
					// get corresponding ADC
					for(long index_subgroup=0; index_subgroup< (long)in_subgroup.size(); index_subgroup++){
						if(in_subgroup[index_subgroup].B1_time == tem_event.B1_time && in_subgroup[index_subgroup].B2_time == tem_event.B2_time){
							tem_event.B1_sweeps_gclock = in_subgroup[index_subgroup].B1_sweeps_gclock;
							tem_event.B1_time_relative = in_subgroup[index_subgroup].B1_time_relative;
							tem_event.B1_adc = in_subgroup[index_subgroup].B1_adc;

							tem_event.B2_sweeps_gclock = in_subgroup[index_subgroup].B2_sweeps_gclock;
							tem_event.B2_time_relative = in_subgroup[index_subgroup].B2_time_relative;
							tem_event.B2_adc = in_subgroup[index_subgroup].B2_adc;
							tem_event.delta_t = in_subgroup[index_subgroup].delta_t;
							break;
						}
					}
					
					selected_event.push_back(tem_event);
			}
		}
		else{
			for(unsigned long index_keep=0;index_keep<time_difference.size(); index_keep++){// both B1 and B2 has no repeat
					bool Have_selected=false;
					tem_event.B1_time =	B1_candidate_group[b2_list_index][index_of_sort_list[index_keep]];
					tem_event.B2_time = b2_time_list[b2_list_index];
					for(unsigned long index_select=0; index_select < selected_event.size();index_select++){
						if(tem_event.B1_time == selected_event[index_select].B1_time){
							Have_selected=true;
							break;
						}

					}
					if(!Have_selected){
						// get corresponding ADC, relative time .......
						for(long index_subgroup=0; index_subgroup< (long)in_subgroup.size(); index_subgroup++){
							if(in_subgroup[index_subgroup].B1_time == tem_event.B1_time && in_subgroup[index_subgroup].B2_time == tem_event.B2_time){
								tem_event.B1_sweeps_gclock = in_subgroup[index_subgroup].B1_sweeps_gclock;
								tem_event.B1_time_relative = in_subgroup[index_subgroup].B1_time_relative;
								tem_event.B1_adc = in_subgroup[index_subgroup].B1_adc;

								tem_event.B2_sweeps_gclock = in_subgroup[index_subgroup].B2_sweeps_gclock;
								tem_event.B2_time_relative = in_subgroup[index_subgroup].B2_time_relative;
								tem_event.B2_adc = in_subgroup[index_subgroup].B2_adc;
								tem_event.delta_t = in_subgroup[index_subgroup].delta_t;
								break;
							}
						}
						
								cout<<"B1: "<<tem_event.B1_time<<" ; B2: "<<tem_event.B2_time<<" deltat:"<<tem_event.delta_t<<endl;
						selected_event.push_back(tem_event); break;
					}
			}
		}

	}

	if(num2keep == 1 && Reverse_direct==0){Reverse_direct=1;  goto Rever;} // to search from B2 earlier to later

	if(num2keep>1){
		for(int index_stored=(int)selected_event.size()-1;index_stored>=0;index_stored--){
			selected_event[index_stored].delta_t = selected_event[index_stored].B2_time - selected_event[index_stored].B1_time;
			Add_to_final.push_back(selected_event[index_stored]);
		}
		return;
	}

/*for(int ij=0; ij<(int)Add_to_final.size();ij++){
	printf("%lld \t %lld\n",Add_to_final[ij].B1_time,Add_to_final[ij].B2_time);
}*/

	if(num2keep == 1 && Reverse_direct==1){ // compare and select to better data set

		double avg_time_diff=0, avg_time_diff_rev=0;
		for(unsigned long index_sele=0;index_sele<selected_event.size();index_sele++){
			avg_time_diff+=Evaluate(selected_event[index_sele].B1_time,selected_event[index_sele].B2_time);
		}
		avg_time_diff/=selected_event.size()*1.0;

		for(unsigned long index_sele=0;index_sele<selected_event_reverse.size();index_sele++){
			avg_time_diff_rev+=Evaluate(selected_event_reverse[index_sele].B1_time,selected_event_reverse[index_sele].B2_time);
		}
		avg_time_diff_rev/=selected_event_reverse.size()*1.0;

		if(avg_time_diff<avg_time_diff_rev){// B2 in sequent better
cout<<"B2 in sequence better"<<endl;
			for(unsigned long index_stored=0;index_stored<selected_event.size();index_stored++){
					selected_event[index_stored].delta_t = selected_event[index_stored].B2_time - selected_event[index_stored].B1_time;
					Add_to_final.push_back(selected_event[index_stored]);
			}

		}
		else{// B2 reversed seqence better
cout<<"B2 reversed sequence better"<<endl;
			for(long index_stored=(long)selected_event_reverse.size()-1;index_stored>=0;index_stored--){
					selected_event_reverse[index_stored].delta_t = selected_event_reverse[index_stored].B2_time - selected_event_reverse[index_stored].B1_time;
					Add_to_final.push_back(selected_event_reverse[index_stored]);
			}
		}

		return;

	}



}//selectbybeta



void SelectByTof(vector<match_event>& in_subgroup, vector<Long64_t>& b1_time_list,vector<match_event>& Add_to_final,int num2keep){// b1_time_list sorted to "earlier" to "later"
	int Reverse_direct=0;
	
	vector<vector<Long64_t>> B2_candidate_group;
	vector<Long64_t> B2_candidate_each_B1;

	vector<match_event> selected_event; // earlier to later
	vector<match_event> selected_event_reverse; // later to earlier

Rever:
	B2_candidate_group.clear(); B2_candidate_group.shrink_to_fit();

	if(Reverse_direct==1){
		for(unsigned long index_sele=0; index_sele<selected_event.size();index_sele++){
			selected_event_reverse.push_back(selected_event[index_sele]);
		}
		selected_event.clear(); selected_event.shrink_to_fit();
	}

	if(Reverse_direct==0){
		sort(b1_time_list.rbegin(),b1_time_list.rend()); // first round B1 later to earlier
	}
	else{
		sort(b1_time_list.begin(),b1_time_list.end());// second round B1 earlier to later
	}


	for(int b1_list_index=0;b1_list_index<(int)(b1_time_list.size());b1_list_index++){// find out all B2 correspond to each B1
		B2_candidate_each_B1.clear();
		for(long index_subgroup=0; index_subgroup< (long)in_subgroup.size(); index_subgroup++){
			if(in_subgroup[index_subgroup].B1_time == b1_time_list[b1_list_index]) B2_candidate_each_B1.push_back(in_subgroup[index_subgroup].B2_time);
		}

		if(Reverse_direct==0){sort(B2_candidate_each_B1.rbegin(),B2_candidate_each_B1.rend());} // "later" to "earlier"
		else{sort(B2_candidate_each_B1.begin(),B2_candidate_each_B1.end());} // "earlier" 2 later

		B2_candidate_group.push_back(B2_candidate_each_B1);

	}

	match_event tem_event;
	vector<double> time_difference;

	for(int b1_list_index=0;b1_list_index<(int)(b1_time_list.size());b1_list_index++){
		time_difference.clear(); 
		time_difference.shrink_to_fit();
		for(int index_column_each_B1=0;index_column_each_B1<(int)B2_candidate_group[b1_list_index].size();index_column_each_B1++){

			time_difference.push_back(Evaluate(b1_time_list[b1_list_index],B2_candidate_group[b1_list_index][index_column_each_B1]));
		}

		double* time_difference_clone = new double[ time_difference.size() ];
		Long64_t* index_of_sort_list = new Long64_t[ time_difference.size() ];
		for(unsigned long in_cp=0;in_cp<time_difference.size();in_cp++){time_difference_clone[in_cp] = time_difference[in_cp];}

		TMath::Sort((Long64_t)time_difference.size(), time_difference_clone, index_of_sort_list,kFALSE); // first one with the index of "smallest"!! value in evaluate list


		if(num2keep>1){
			for(unsigned long index_keep=0;index_keep<time_difference.size() && index_keep< (unsigned int)num2keep; index_keep++){// smallest B1-B2
					tem_event.B1_time =	b1_time_list[b1_list_index];
					tem_event.B2_time = B2_candidate_group[b1_list_index][index_of_sort_list[index_keep]];
					 // get corresponding ADC, relative time .......
					for(long index_subgroup=0; index_subgroup< (long)in_subgroup.size(); index_subgroup++){
						if(in_subgroup[index_subgroup].B1_time == tem_event.B1_time && in_subgroup[index_subgroup].B2_time == tem_event.B2_time){
							tem_event.B1_sweeps_gclock = in_subgroup[index_subgroup].B1_sweeps_gclock;
							tem_event.B1_time_relative = in_subgroup[index_subgroup].B1_time_relative;
							tem_event.B1_adc = in_subgroup[index_subgroup].B1_adc;

							tem_event.B2_sweeps_gclock = in_subgroup[index_subgroup].B2_sweeps_gclock;
							tem_event.B2_time_relative = in_subgroup[index_subgroup].B2_time_relative;
							tem_event.B2_adc = in_subgroup[index_subgroup].B2_adc;
							tem_event.delta_t = in_subgroup[index_subgroup].delta_t;
							break;
						}
					}
					
					selected_event.push_back(tem_event);
			}
		}
		else{
			for(unsigned long index_keep=0;index_keep<time_difference.size(); index_keep++){// both B1 and B2 has no repeat
					bool Have_selected=false;
					tem_event.B1_time =	b1_time_list[b1_list_index];
					tem_event.B2_time = B2_candidate_group[b1_list_index][index_of_sort_list[index_keep]];
					for(unsigned long index_select=0; index_select < selected_event.size();index_select++){
						if(tem_event.B2_time == selected_event[index_select].B2_time){
							Have_selected=true;
							break;
						}

					}
					if(!Have_selected){
						// get corresponding ADC, relative time .......
						for(long index_subgroup=0; index_subgroup< (long)in_subgroup.size(); index_subgroup++){
							if(in_subgroup[index_subgroup].B1_time == tem_event.B1_time && in_subgroup[index_subgroup].B2_time == tem_event.B2_time){
								tem_event.B1_sweeps_gclock = in_subgroup[index_subgroup].B1_sweeps_gclock;
								tem_event.B1_time_relative = in_subgroup[index_subgroup].B1_time_relative;
								tem_event.B1_adc = in_subgroup[index_subgroup].B1_adc;

								tem_event.B2_sweeps_gclock = in_subgroup[index_subgroup].B2_sweeps_gclock;
								tem_event.B2_time_relative = in_subgroup[index_subgroup].B2_time_relative;
								tem_event.B2_adc = in_subgroup[index_subgroup].B2_adc;
								tem_event.delta_t = in_subgroup[index_subgroup].delta_t;
								break;
							}
						}
						
						cout<<"B1: "<<tem_event.B1_time<<" ; B2: "<<tem_event.B2_time<<" deltat:"<<tem_event.delta_t<<endl;
						selected_event.push_back(tem_event); break;

					}
			}
		}

	}

	if(num2keep == 1 && Reverse_direct==0){Reverse_direct=1;  goto Rever;} // to search from B2 earlier to later


	if(num2keep>1){
		for(long index_stored=(long)selected_event.size()-1;index_stored>=0;index_stored--){
			selected_event[index_stored].delta_t = selected_event[index_stored].B2_time - selected_event[index_stored].B1_time;
			Add_to_final.push_back(selected_event[index_stored]);
		}
		return;
	}



	if(num2keep == 1 && Reverse_direct==1){ // compare and select to better data set

			double avg_time_diff=0, avg_time_diff_rev=0;
			for(unsigned long index_sele=0;index_sele<selected_event.size();index_sele++){
				avg_time_diff+=Evaluate(selected_event[index_sele].B1_time,selected_event[index_sele].B2_time);
			}
			avg_time_diff/=selected_event.size()*1.0;

			for(unsigned long index_sele=0;index_sele<selected_event_reverse.size();index_sele++){
				avg_time_diff_rev+=Evaluate(selected_event_reverse[index_sele].B1_time,selected_event_reverse[index_sele].B2_time);
			}
			avg_time_diff_rev/=selected_event_reverse.size()*1.0;


			if(avg_time_diff<avg_time_diff_rev){// B1 in sequent better
cout<<"B1 in sequence better"<<endl;
				for(unsigned long index_stored=0;index_stored<selected_event.size();index_stored++){
						selected_event[index_stored].delta_t = selected_event[index_stored].B2_time - selected_event[index_stored].B1_time;
						Add_to_final.push_back(selected_event[index_stored]);
				}

			}
			else{// B1 reversed seqence better
cout<<"B1 reversed sequence better"<<endl;
				for(long index_stored=(long)selected_event_reverse.size()-1;index_stored>=0;index_stored--){
						selected_event_reverse[index_stored].delta_t = selected_event_reverse[index_stored].B2_time - selected_event_reverse[index_stored].B1_time;
						Add_to_final.push_back(selected_event_reverse[index_stored]);
				}
			}

			return;

	}
/*
for(int ij=0; ij<(int)Add_to_final.size();ij++){
	printf("%lld \t %lld\n",Add_to_final[ij].B1_time,Add_to_final[ij].B2_time);
}*/


}


#endif