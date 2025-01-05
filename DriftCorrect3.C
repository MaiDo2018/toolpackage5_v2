//%%%%%%%%%% for debut, enable showgraph and showfitting %%%%%%%%%%%%%%%

#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TSystem.h"

//#include "TProfile.h"
//#include "TAxis.h"

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "TGProgressBar.h"
#include "TGClient.h"
#include "funcJohnson.h"

//#define _LASER_


#if defined (MAKECINT)
#pragma link C++ class vector<Long_t>+;   // in order to identify vector<Long_t> 
#endif


using namespace std;

bool abort_corr_loop=false;

Double_t gas(double*x,double*par){
   return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2]));  

}

Double_t exR(double*x,double*par){
  return par[0]*TMath::Exp(par[3]*(par[3]-2*(x[0]-par[1]))/(2*par[2]*par[2]));   // right long tail
 
}

Double_t exL(double*x,double*par){
  return par[0]*TMath::Exp(par[4]*(par[4]-2*(x[0]-par[1]))/(2*par[2]*par[2]));// left long tail
 
}

Double_t fitfunc(double*x,double*par){  // par sequence=> high, cento, sigma, right tail tc, left tail tc
   /* for right long tail, all parameters>0 */
    if(x[0]<=par[1]+par[3]&&x[0]>=par[1]+par[4]) return gas(x,par); //+ linearbg(x,&par[4]);

    /* for Right long tail , right tail para>0*/
    else if(x[0]>par[1]+par[3]) return exR(x,par); //+ linearbg(x,&par[4]);
    
  /* for left long tail , left tail para<0*/
    else  if(x[0]<par[1]+par[4]) return exL(x,par);
    
    return 0;
   
}


double* FeaturePar(TH1D* _inhisto,Long64_t _min_evt, Long64_t _max_evt);
double fit_result(TTree *t1, Long64_t min_evt, Long64_t max_evt, int nbins, double lowx, double highx,Double_t Sigma, int selectTag=1, int selectCH=1,double Time_veto=0);
TH1D *h;
TF1 *func=NULL, *fgaus_exp=NULL;
funcJohnson* john=NULL;

string fitname ="john";

//%%%%%%%%%%% display for debut %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TCanvas *c1= NULL; // = new TCanvas("c1","c1",1000,800); 
bool showgraph=false; // alow to show histo of empty slice
bool showfitting=false; // showgraph && showfitting to show each fit
bool showbadfit= false; // set to true --> show graph when it is bad fit
//TCanvas *c_p1_p2 = new TCanvas("c_p1_p2","t0 correction",1200,1000);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TGraph *g_mean_p1 = new TGraph();
TGraph *g_mean_p2 = new TGraph();
TGraph *g_sigma_p1 = new TGraph();
TGraph *g_sigma_p2 = new TGraph();

double event_akill_global=0;
bool isnewfile=true;

Long64_t stop_point = -1;  // stop_point marker set to -1
int N_byfit=0;      // Number of slices by fit
int N_byappro=0;    // Number of slices by approximating

Long_t empty_slice_gsweeps_L=0; // gsweeps of empty slice
Long_t empty_slice_gsweeps_R=0; // gsweeps of empty slice

double fit_half_rangeL=0;
double fit_half_rangeR=0;
bool fit_by_free_pars=true; // true for rough fit;

int DriftCorrect(string PATH="../",string filename = "mcs_39Kvs143X2plus@500_174927", double centro=17e6, int nbins = 200, int event_akill =600, double ref_sigma=10,
                          double Half_hiswidth = 150,double _inT0=130,TGProgressBar *bar=NULL, int selectTag=1, int selectCH=1,double Time_veto=0,
                          double * outpars=NULL,vector<double> *timesANDsweep=NULL){

    //%%%%%%%%%%%%%%%%%%%%%% display purpose for debut %%%%%%%%%%%%%%%%%%%%%55

          if(c1==NULL && showgraph) c1 = new TCanvas("c1","c1",1000,800);
          if(!showgraph){ 
            if(c1!=NULL && c1->GetCanvasImp()!=NULL){ delete c1; c1=NULL; }
            if(c1!=NULL && c1->GetCanvasImp()==NULL) c1=NULL;
          } 

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      string inputfile = PATH + "rootfiles/" + filename + ".root";
      TFile *fin = new TFile(inputfile.c_str(),"READ");
      if(!fin->IsOpen()) {cout<<"Can not open file: "<<inputfile<<endl; cout<<"Break!!"<<endl; return 0;}
      TTree *intree = (TTree *)fin->Get("tree");
	if(intree==NULL)intree = (TTree *)fin->Get("tree0");
	if(intree==NULL){cout<<"Can not find tree or tree0; Break!!!"<<endl; return 0;}
      TTree* tree_copy = intree;

      //%%%%%% initialize %%%%%%%%%%%%%%%
      isnewfile=true;
 
	event_akill_global = event_akill;

    if(fit_by_free_pars){ 
        if(fitname == "john"){
          if(john != NULL) delete john;
          john = new funcJohnson();
        }
        else if(fitname == "fgaus_exp"){
          if(fgaus_exp != NULL) delete fgaus_exp;
          fgaus_exp = new TF1("fgaus_exp",fitfunc,0,26e6,5);
          fgaus_exp->SetParLimits(2,0,3*ref_sigma);
          fgaus_exp->SetParLimits(3,ref_sigma*0.2,ref_sigma*5);
          fgaus_exp->SetParLimits(4,-ref_sigma*5,-ref_sigma*0.2);
          fgaus_exp->SetParameters(0,0,ref_sigma,ref_sigma,-ref_sigma);
        }

    }

    if(fitname == "john" && john != NULL){
        func = john->GetTF1();
        func->SetName("john");
    }
    else if(fitname == "fgaus_exp" && fgaus_exp!= NULL){
        func = fgaus_exp;
        func->SetName("fgaus_exp");
    } 

    if(func == NULL){
      cout<<"\e[1;31m Fail to create fitting function!!! Check the \"fitname\" should be either \"john\" or \"fgaus_exp\" . Abort!!\e[0m"<<endl;
      return 0;
    }

      N_byfit=0;
      N_byappro=0;
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      cout<<"\e[1;33m reading data from "<<inputfile<<"\e[0m"<<endl;


  Long64_t nentries = intree->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  if(outpars!=NULL)outpars[0] = nentries;

  //Long64_t stop_point = -1;        // stop_point marker set to -1
  stop_point=-1;
  int counter = 1;    // display purpose 
  int SlideWidth = event_akill;  // number of events to load in each time
  int num_of_SlideWidth = (int) (nentries / SlideWidth);
  cout<<"num_of_SlideWidth = "<<num_of_SlideWidth<<endl;
  double binlowedge1 = centro - Half_hiswidth;    double binlowedge2 = centro - Half_hiswidth;// total histowidth should be 300 usually
  double binupedge1 = centro + Half_hiswidth;     double binupedge2 = centro + Half_hiswidth;
  int    nbins1    = nbins;           int    nbins2    = nbins;
  double sigma_ref_measured1 = ref_sigma;   double sigma_ref_measured2 = ref_sigma; //;6
  int accumu_num = 0;
  double reference1 =0;              double reference2 =0;           // fit first slidewidth data as reference
  double n_slide_central1 = 0;       double n_slide_central2 = 0;
  int    node = 0 ;
                                           //(1000,0,794000,1000,9.006e6,9.007e6) ;  1600,0,793001,100,9.0064e6,9.0067e6                             
                                            //(1000,0,1771000,1000,9.019e6,9.021e6)
                                             //(1000,0,181000,1000,9.0423e6,9.0426e6);

  if(outpars!=NULL){
    outpars[0] = nentries;
    outpars[1] = num_of_SlideWidth;
    outpars[2] = nbins;
    outpars[3] = binlowedge1;
    outpars[4] = binupedge1;
  }


  bool fitTrue = false;
  bool fitTrue3 = false;
//  double offset_stop =0;
//  vector <Double_t> *offset = new vector <Double_t>();

  const int t0_num = 1;         // check effect of 11 kinds of daq delay t0 to drift correction
  // !!!!!!!!!!!! modified to allow setting a single t0
  // in the following:  
  //double t0 = index*t0_delta + t0_delta; index=[0,t0_num) => index=0 only; t0 = t0_delta finally; 
  // by setting t0_delta = setting t0
  ///////////////////////////////////
  double t0_delta = _inT0; //130;    // step size of t0 = 20 ns for multi t0 test situation, t0 delay is about 130 ns usually
  double offset_stop[t0_num] ={0};
  vector <Double_t> **offset = new vector <Double_t> *[t0_num];

	for(int index =0;index<t0_num;index++){offset[index] = new vector <Double_t>[1];}

  if(bar!=NULL) bar->Reset();
 
  for(Long64_t jentry=0; jentry<nentries; jentry++){       //nentries
    if(abort_corr_loop) break;

    if(jentry == (int)counter*nentries/20){  // divide into 20 piese , each piese is equivalent to 5%
      cout<< '\r' << "finish running : "<< counter*5 << "%" <<"\t ratio by fit= "<<N_byfit/(N_byfit+2.0*N_byappro)*100.<<" %";
      cout<<"\t ratio by approximate= "<<2*N_byappro/(N_byfit+2.0*N_byappro)*100.<<" %"<<flush; 
      if(bar != NULL){bar->SetPosition((Float_t)counter*5.0); bar->RaiseWindow();gSystem->ProcessEvents();} // gSystem->ProcessEvents() keep window alive
      counter++;
     }

//     cout<<endl;
 //    cout<<"1_nevt = "<<num<<endl; //return 1;
       
     if(jentry == 1){
           cout<<"continue or not: 'n' for break!!!!"<<endl;
             char go;
             cin>>go;
            if( go =='n'){cerr<<"\e[1;33m"<<"program break!!!!!"<<"\e[0m"<<endl; return 0;}
            cout<<"Pending......"<<endl;
      }

    // correct drift time
    if((jentry < SlideWidth) && fitTrue == false){//%%%%%%%%%%% first piece
    cout<<"do 1"<<endl;
    
   // cout<<"do 1 num1 , nevt =  "<<num<< " , "<<nevt<<endl;     
  
      reference1 = fit_result(tree_copy,0,((Long64_t) SlideWidth),nbins1,binlowedge1,binupedge1,sigma_ref_measured1,selectTag,selectCH,Time_veto);  //(TTree *t1, Long64_t min_evt, Long64_t max_evt, 
              // cout<<"do 1 num2, nevt = "<<num << " , "<<nevt<<endl;
               cout<<"referen 1 = "<<reference1<<endl;                                         //int nbins, double lowx, double highx);
      reference2 = fit_result(tree_copy,0,((Long64_t) SlideWidth),nbins2,binlowedge2,binupedge2,sigma_ref_measured2,selectTag,selectCH,Time_veto);
               cout<<"referen 2 = "<<reference2<<endl; 
            //offset->push_back(1);  // shift==0; ration ==1
         for(int index =0;index<t0_num;index++){offset[index]->push_back(1);}
       fitTrue = true;
                 
    }
    else if( (node<(int)(jentry/SlideWidth)) && (int)(jentry/SlideWidth)< num_of_SlideWidth){ //%%%%%%%%% middle piece
   // cout<<"do 2"<<endl;

          if(stop_point < 0){
              n_slide_central1 = fit_result(tree_copy,jentry,(jentry + (Long64_t) SlideWidth),nbins1,binlowedge1,binupedge1,sigma_ref_measured1,selectTag,selectCH,Time_veto);// node here delays by 1;
              n_slide_central2 = fit_result(tree_copy,jentry,(jentry + (Long64_t) SlideWidth),nbins2,binlowedge2,binupedge2,sigma_ref_measured2,selectTag,selectCH,Time_veto);
          }else{
              n_slide_central1 = fit_result(tree_copy,stop_point,(jentry + (Long64_t) SlideWidth),nbins1,binlowedge1,binupedge1,sigma_ref_measured1,selectTag,selectCH,Time_veto); 
              n_slide_central2 = fit_result(tree_copy,stop_point,(jentry + (Long64_t) SlideWidth),nbins2,binlowedge2,binupedge2,sigma_ref_measured2,selectTag,selectCH,Time_veto);
          }


          if(n_slide_central1>0 && n_slide_central2>0){   // good fit
                           accumu_num++;
                 for(int akill=0;akill<accumu_num;akill++){
                    // offset_stop = ((n_slide_central1 - reference1)+(n_slide_central2 - reference2))/2;     //update offset record
				//	offset_stop = ((n_slide_central1 / reference1)+(n_slide_central2 / reference2))/2;
                       	//	offset->push_back( offset_stop ) ;  // update offset
              				for(int index=0;index<t0_num;index++){
              					double t0 = index*t0_delta + t0_delta;   // + t0_delta at the end just want t0 to start from 130 ns 
              					offset_stop[index] = ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
              					//offset_stop[index] = ((n_slide_central1 - reference1)+(n_slide_central2 - reference2))/2; // by difference method
              					//offset_stop[index] = ((reference1 / n_slide_central1)+(reference2 / n_slide_central2))/2;   // by ratio method
                                     		offset[index]->push_back( offset_stop[index] ) ;  // update offset
              				}

                  }

                        accumu_num =0;          // reset to 0;
                        stop_point = -1;       // reset stop_point

          }
          else if(n_slide_central1==-20 || n_slide_central2==-20){ // empty slice!!!!!
                accumu_num++;
                for(int akill=0;akill<accumu_num;akill++){
                          for(int index=0;index<t0_num;index++){
                               // double t0 = index*t0_delta + t0_delta;   // + t0_delta at the end just want t0 to start from 130 ns 
                               // offset_stop[index] = ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
                              //  offset_stop[index] = offset[index]->at(offset[index]->size()-1);
                                offset[index]->push_back( offset_stop[index] ) ;  
                          }

                }
                cout<<"\e[1;32m"<<"warning: empty slice founded!!!  "<<empty_slice_gsweeps_L<<" ~ "<<empty_slice_gsweeps_R<<"\e[0m"<<endl;

                accumu_num =0;          // reset to 0;
                stop_point = -1;       // reset stop_point
          }
          else{ // bad fit to accumulate more slideWidth
                       if(accumu_num ==0){ stop_point = jentry;}  // set stop_point
                  accumu_num ++;
                /* if(accumu_num>20) {
                      cerr<<"warning !!!!!! find no peak!!! check point, from nevt= "<<stop_point<<"accumulate num = "<<accumu_num<<endl; 
                           return 0;
                     }*/
          }

       //     cout<<"offset = , node = "<<offset_stop[0]<<"  ,  "<<node+1<<".  nevt range = "<<(node+1)*SlideWidth<<" ~ "<<(node+2)*SlideWidth<<endl; 
       //        cout<<endl;
            //cout<<"jentry= "<<jentry<<endl;

            node++;  //code above only execute once at the beginning of each piece; point to next piece       
              

    }
    else if((int)(jentry/SlideWidth) == num_of_SlideWidth && fitTrue3 == false){ //%%%%%%%%%%% last piese
    //   cout<<"do final"<<endl;

          if(stop_point < 0){
             n_slide_central1 = fit_result(tree_copy,jentry,nentries,nbins1,binlowedge1,binupedge1,sigma_ref_measured1,selectTag,selectCH,Time_veto);
             n_slide_central2 = fit_result(tree_copy,jentry,nentries,nbins2,binlowedge2,binupedge2,sigma_ref_measured2,selectTag,selectCH,Time_veto);
          }else{
              n_slide_central1 = fit_result(tree_copy,stop_point,nentries,nbins1,binlowedge1,binupedge1,sigma_ref_measured1,selectTag,selectCH,Time_veto); 
              n_slide_central2 = fit_result(tree_copy,stop_point,nentries,nbins2,binlowedge2,binupedge2,sigma_ref_measured2,selectTag,selectCH,Time_veto);
          }

             // %%%%%%%%%%%%   good fit or bad fit %%%%%%%%%%%%%%
          if(n_slide_central1>0 && n_slide_central2>0){   // good fit
                           accumu_num++;
                 for(int akill=0;akill<accumu_num;akill++){
                   //  offset_stop = ((n_slide_central1 - reference1) + (n_slide_central2 - reference2))/2;
			//	offset_stop = ((n_slide_central1 / reference1)+(n_slide_central2 / reference2))/2;
                   //    offset->push_back( offset_stop ) ;  // update offset

              				for(int index=0;index<t0_num;index++){
              					double t0 = index*t0_delta + t0_delta;  // + t0_delta at the end just want t0 to start from 130 ns
              					offset_stop[index] =  ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
              					//offset_stop[index] = ((n_slide_central1 - reference1)+(n_slide_central2 - reference2))/2;  // by difference method
              					//offset_stop[index] = ((reference1 / n_slide_central1)+(reference2 / n_slide_central2 ))/2;   // by ratio method
                                     		offset[index]->push_back( offset_stop[index] ) ;  // update offset
              				}

                  }
                        accumu_num =0;          // reset to 0;
                        stop_point = -1;       // reset stop_point

          }
         else if(n_slide_central1==-20 || n_slide_central2==-20){ // empty slice!!!!!
                accumu_num++;
                for(int akill=0;akill<accumu_num;akill++){
                          for(int index=0;index<t0_num;index++){
                               // double t0 = index*t0_delta + t0_delta;   // + t0_delta at the end just want t0 to start from 130 ns 
                               // offset_stop[index] = ((n_slide_central1-t0)/(reference1-t0) + (n_slide_central2-t0)/(reference2-t0) )/2;
                               // offset_stop[index] = offset[index]->at(offset[index]->size()-1);
                                offset[index]->push_back( offset_stop[index] ) ;  // update offset with offset from last slice
                          }

                }
                cout<<"\e[1;32m"<<"warning: empty slice founded!!!  "<<empty_slice_gsweeps_L<<" ~ "<<empty_slice_gsweeps_R<<"\e[0m"<<endl;

                accumu_num =0;          // reset to 0;
                stop_point = -1;       // reset stop_point
          }
          else{ // bad fit 

                           accumu_num++;
                 for(int akill=0;akill<accumu_num;akill++){

                     //  offset->push_back( offset_stop ) ;  // use old offset to approximate
                     //  offset->push_back( offset_stop[index] ) ; 
            				for(int index=0;index<t0_num;index++){
            					offset[index]->push_back( offset_stop[index] ) ;  // update offset with offset from last slice
            				}
                  }
          }


           // cout<<"end of fit : offset = "<<offset_stop[0]<<endl;    
          // cout<<"jentry= "<<jentry<<endl; 

              fitTrue3 = true;  // make sure last piece fit only once
        
    }else{;}
  }// end of for offset extraction


  ////%%%%%%%%%%%%%%%%%%%%%%% end of offset extraction %%%%%%%%%%%%%%%%%%%%%%%%

    if(abort_corr_loop){
       cout<<"external interrupt"<<endl;
       return 0;
    }

     cout<<endl;
     cout<<"nentries = "<<nentries<<endl;
     cout<<"offset num = "<< offset[0]->size()<<endl;
     cout<<endl;
     cout<<endl;
   
      if(bar != NULL) bar->SetPosition((Float_t)counter*5.0); 


      Long64_t nevt = 0;
 //     Long64_t num = 0;
      int nhits[10];
//      int glo = 0;
      vector <Int_t> *channel = new vector <Int_t>();
      vector <Int_t> *edge = new vector <Int_t>();
      vector <Int_t> *value = new vector <Int_t>();
      vector <Double_t> *time = new vector <Double_t>();
      vector <Double_t> *timec = new vector <Double_t>();
      vector <Int_t> *sweeps = new vector <Int_t>();
      vector <Long_t> * sweeps_global = new vector <Long_t>;
      vector <Int_t> *glo = new vector <Int_t>();
      vector <Int_t> *tag = new vector <Int_t>();
      vector <Int_t> *lost = new vector <Int_t>();
      #ifdef _LASER_
      vector <Double_t> *time_gcDrift = new vector<Double_t>();
      #endif
      double offset_newtree = 0;
      
      intree->SetBranchAddress("nevt", &nevt); 
//      intree->SetBranchAddress("nevt", &num); 
      intree->SetBranchAddress("nhits", nhits);
      intree->SetBranchAddress("channel", &channel);
      intree->SetBranchAddress("edge", &edge);
      intree->SetBranchAddress("value", &value);
      intree->SetBranchAddress("time", &time);
      intree->SetBranchAddress("sweeps", &sweeps);
      intree->SetBranchAddress("sweeps_global", &sweeps_global);
      intree->SetBranchAddress("glo", &glo);
      intree->SetBranchAddress("tag", &tag);
      intree->SetBranchAddress("lost", &lost);
      #ifdef _LASER_
      intree->SetBranchAddress("time_gcDrift", &time_gcDrift);
      #endif
  // correct drift time

  string outputfile = PATH + "rootfiles/" + filename +"_dcorrect.root";
  TFile *fout = new TFile(outputfile.c_str(),"RECREATE");
//  TTree *tree = new TTree("tree","analyzed tree");
  TTree **tree = new TTree *[t0_num];
  TH2D *h2p1 = new TH2D("h2p1","h2p1",1000,0,794000,200,9.0063e6,9.0068e6);
  TH2D *h2p2 = new TH2D("h2p2","h2p2",1000,0,794000,100,9.0134e6,9.0139e6);
  TH1D *h1p1 = new TH1D("h1p1","h1p1",200,9.0063e6,9.0068e6);
  TH1D *h1p2 = new TH1D("h1p2","h1p2",100,9.0134e6,9.0139e6);


for(int index_t0=0;index_t0<t0_num;index_t0++){
	double t0 = index_t0 * t0_delta + t0_delta;
  tree[index_t0]= new TTree(Form("tree%d",index_t0),"analyzed tree");
  tree[index_t0]->Branch("nevt",&nevt);
//  tree->Branch("nevt1",&num);
  tree[index_t0]->Branch("nhits",nhits,"nhits[10]/I");
  tree[index_t0]->Branch("channel",&channel);
  tree[index_t0]->Branch("edge",&edge);
  tree[index_t0]->Branch("value",&value); // in unit of 100ps
  tree[index_t0]->Branch("time",&time); // in unit of ns
  tree[index_t0]->Branch("timec",&timec); // in unit of ns
  tree[index_t0]->Branch("sweeps",&sweeps);
  tree[index_t0]->Branch("sweeps_global",&sweeps_global);
  tree[index_t0]->Branch("glo",&glo);
  tree[index_t0]->Branch("tag",&tag);
  tree[index_t0]->Branch("lost",&lost);
  tree[index_t0]->Branch("offset",&offset_newtree);
  #ifdef _LASER_
  tree[index_t0]->Branch("time_gcDrift",&time_gcDrift);
  #endif

 
       counter =0;

  for(Long64_t jentry=0; jentry<nentries; jentry++){
    if(jentry == (int)counter*nentries/20){  // divide into 20 piese , each piese is equivalent to 5%
      cout<< '\r' << "finish running : "<< counter*5 << "%" <<flush;
      counter++;
     }
   
       node = (int)(jentry/SlideWidth);

        intree->GetEntry(jentry);
         
     //   offset_newtree = offset->at(node);


       // offset_newtree = offset->at(node);
        offset_newtree = offset[index_t0]->at(node);
     // calculate timec
    double itime=0;
    double itime_last=0;

    for(int j=0; j<((int)time->size()); j++){
       itime = time->at(j);
      //timec->push_back(itime-offset_newtree);  // shift method correct
      if(itime==0){timec->push_back(0);}
      else timec->push_back((itime-t0) / offset_newtree + t0);    // ratio method correct

        if(timec->size()>=1){
          if(tag->at(j)==selectTag && channel->at(j)==selectCH){
            if(itime - itime_last > Time_veto){
                timesANDsweep[0].push_back(  timec->at(timec->size()-1)  );
                timesANDsweep[1].push_back(jentry);
                itime_last = itime;
            }
          }
        }

      // cout<<"offset = "<<offset<<", timec = "<<itime-offset<<endl;
    }// calculate timec



    tree[index_t0]->Fill();
    
    channel->clear();
    edge->clear();

    value->clear();
    time->clear();
    timec->clear();
    sweeps->clear();
    sweeps_global->clear();
    glo->clear();
    tag ->clear();
    lost->clear(); 
    #ifdef _LASER_
    time_gcDrift->clear();
    #endif
    for(int i=0; i<10; i++) nhits[i] = 0;


  }// end of for getEvent loop



  cout<<endl;

  cout<<"\e[1;33m writing tree to "<<outputfile<<"\e[0m"<<endl;
  fout->cd();
  tree[index_t0]->Write();
  

  fout->Close();
  fin->Close();

  
   /*for(int index_t0=0;index_t0<t0_num;index_t0++){
     
       delete tree[index_t0];
    }*/

    delete[] tree;

    cout<<endl;
    cout<<"nentries = "<<nentries<<endl;
    cout<<"num_of_SlideWidth = "<<num_of_SlideWidth<<endl;
    cout<<"ratio by fit= "<<N_byfit/(N_byfit+2.0*N_byappro)*100.<<"%";
    cout<<"\t ratio by approximate= "<<N_byappro*2./(N_byfit+2.0*N_byappro)*100.<<"%"<<endl;
    cout<<"Maximum sweeps = "<<empty_slice_gsweeps_R<<endl;
    if(bar != NULL)bar->SetPosition(100);

/*
   if(index_t0==0){ c_p1_p2->Divide(2,2);}
   //if(index_t0 !=0) {delete fun_g1;delete fun_g2;}
   c_p1_p2->cd(1);
   tree[index_t0]->Draw("timec:nevt>>h2p1","","colz");
   c_p1_p2->cd(2);
   tree[index_t0]->Draw("timec:nevt>>h2p2","","colz");
   c_p1_p2->cd(3);
   tree[index_t0]->Draw("timec>>h1p1","","colz");
   fun_g1 = new TF1("fun_g1","gaus",0,1e9);
   h1p1->Fit(fun_g1,"","",9006.4e3,9006.6e3);
   c_p1_p2->cd(4);
   tree[index_t0]->Draw("timec>>h1p2","","colz");
   fun_g2 = new TF1("fun_g2","gaus",0,1e9);
   h1p2->Fit(fun_g2,"","",9013.4e3,9013.8e3);

   g_mean_p1->SetPoint(index_t0,index_t0*t0_delta,fun_g1->GetParameter(1));
   g_sigma_p1->SetPoint(index_t0,index_t0*t0_delta,fun_g1->GetParameter(2));   
   g_mean_p2->SetPoint(index_t0,index_t0*t0_delta,fun_g2->GetParameter(1));
   g_sigma_p2->SetPoint(index_t0,index_t0*t0_delta,fun_g2->GetParameter(2));

   c_p1_p2->Update();
   if(index_t0==0)  c_p1_p2->WaitPrimitive();
   delete fun_g1;
   delete fun_g2;
*/

} // end of for loop of t0
/*
   c1->Divide(2,2);
   c1->cd(1);
   g_mean_p1->SetTitle("p1 mean;t0 [ns];mean");
   g_mean_p1->Draw("APL*");
   c1->cd(3);
   g_sigma_p1->SetTitle("p1 sigma;t0 [ns];sigma");
   g_sigma_p1->Draw("APL"); 
   c1->cd(2);  
   g_mean_p2->SetTitle("p2 mean;t0 [ns];mean");
   g_mean_p2->Draw("APL*");
   c1->cd(4);
   g_sigma_p2->SetTitle("p2 sigma;t0 [ns];sigma");
   g_sigma_p2->Draw("APL");
*/

  return 1;
}


// dynamic hiso range and dynamic sigam 

double fit_result(TTree *t1, Long64_t min_evt, Long64_t max_evt, int nbins, double lowx, double highx,Double_t Sigma, int selectTag, int selectCH,double Time_veto){
    //gROOT->SetBatch(true);
       double n_slide_central = 0;
       static double active_lowx = lowx;
       static double active_highx = highx;
       static bool appro_lasttime = false; // use approximate last time
  
       if(isnewfile){
        active_lowx = lowx;
        active_highx = highx;
        appro_lasttime = false;
         isnewfile=false;
       }
        double Half_hiswidth = (highx-lowx)*0.5;  

    // for search peak and temp result of fitting;//&&&&&&&&&&&  not use now; cover by sigma_StdDev

       // active_lowx = h->GetXaxis()->GetBinCenter(h->GetMaximumBin())-Half_hiswidth;
       // active_highx = n_slide_central+Half_hiswidth;
             
       h = new TH1D("h","h",nbins,active_lowx,active_highx);
        
 
         
       vector <Double_t> *tmp_time = new vector <Double_t>();
       vector <Long_t> * tmp_sweeps_global = new vector <Long_t>;
       vector <Int_t> *tmp_channel = new vector <Int_t>();
       vector <Int_t> *tmp_tag = new vector <Int_t>();
       t1->SetBranchAddress("time",&tmp_time);
       t1->SetBranchAddress("sweeps_global",&tmp_sweeps_global);
       t1->SetBranchAddress("tag",&tmp_tag);
       t1->SetBranchAddress("channel",&tmp_channel);       

       tmp_time->clear();
       tmp_sweeps_global->clear();
       tmp_tag->clear();
       tmp_channel->clear();

       vector<Long_t> temprecord_sweeps;  // record global sweeps of all hits 


        for(Long64_t index = min_evt;index< max_evt;index++){

               t1->GetEvent(index);
               int nhits = (int) (tmp_time->size());

               double last_tof=0;

              for(int n_hit=0;n_hit < nhits;n_hit++){
                  temprecord_sweeps.push_back( tmp_sweeps_global->at(n_hit) );
                  if(tmp_tag->at(n_hit)!=selectTag) continue;
                  if(tmp_channel->at(n_hit)!=selectCH && tmp_channel->at(n_hit)!=6) continue;
                  if(tmp_time->at(n_hit)-last_tof > Time_veto){
                      h->Fill( (Double_t) tmp_time->at(n_hit) );
                      last_tof=tmp_time->at(n_hit);
                  }
              }
                  
             tmp_time->clear();
             tmp_sweeps_global->clear();
             tmp_tag->clear();
             tmp_channel->clear();


        }     

       // t1->ResetBranchAddresses();

      if(temprecord_sweeps.size()>0){
        empty_slice_gsweeps_L = *min_element(temprecord_sweeps.begin(),temprecord_sweeps.end());
        empty_slice_gsweeps_R = *max_element(temprecord_sweeps.begin(),temprecord_sweeps.end());
      }
        temprecord_sweeps.clear();

         // h->Draw();
//c1->Modified();



        double* getfeature = FeaturePar(h,min_evt,max_evt);  // [0] Amp; [1] Peak center; [2] sigma; [3] bin width

        if(getfeature[6]<(max_evt-min_evt)*0.5*0.1){// empty slice; x0.5 ==> only half of them from tag1

if(showgraph){
c1->cd();
h->Draw();
c1->Modified();
c1->Update();
char tempc;
while(1){
  cout<<"empty slice!!!!!!!!! 'y' to continue, 'n' to abort! "<<endl;// " <<empty_slice_gsweeps_L<<" ~ "<<  empty_slice_gsweeps_R <<"continue ?"<<endl;
  cin>>tempc;
  if(tempc!='y' && tempc!='n') continue;
  if(tempc =='y') break;
  else{
    abort_corr_loop = true;
    break;
  }
}
}

            delete h;
            t1->ResetBranchAddresses();
            delete tmp_time;
            delete tmp_sweeps_global;
            delete tmp_tag;
            delete tmp_channel;
            return -20;  // value represents empty slice
        }

        //%%%%%%%%%%%%%%%%%%%%   renew histogram with new feature paras  %%%%%%%%%%%%%%%%%%
        /*
                    int bin_N = TMath::Nint((active_highx-active_lowx)/getfeature[3]);
                    if(bin_N<10){
                      cout<<"binN="<<bin_N<<endl;
                      cout<<"exp binwidth = "<<getfeature[3]<<endl;
                      cout<<"exp std ="<<getfeature[2]<<endl;
                      c1->cd();
                      h->Draw();
                    c1->Update();
                    c1->WaitPrimitive();
                    return -10;
                    }
          */

        int Nbins_cal = TMath::Nint((active_highx - active_lowx) / (getfeature[3]*0.5) ); // just to add a 0.8 factor to use finner bin size

        if (Nbins_cal > 10) {  // use new Nbins define histogram

              delete h;
              h = new TH1D("h","h",Nbins_cal,active_lowx,active_highx);

             for(Long64_t index = min_evt;index< max_evt;index++){

                     t1->GetEvent(index);
                     int nhits = (int) (tmp_time->size());

                     double last_tof=0;

                 for(int n_hit=0;n_hit < nhits;n_hit++){
                      if(tmp_tag->at(n_hit)!=selectTag) continue;
                      if(tmp_channel->at(n_hit)!=selectCH && tmp_channel->at(n_hit)!=6) continue;
                      if(tmp_time->at(n_hit)-last_tof > Time_veto){
                        h->Fill( (Double_t) tmp_time->at(n_hit) );
                        last_tof=tmp_time->at(n_hit);
                      }

                  }
                        
                   tmp_time->clear();
                   tmp_sweeps_global->clear();
                   tmp_tag->clear();
                   tmp_channel->clear();

              }  

              getfeature = FeaturePar(h,min_evt,max_evt);
        }
        else{
                      cout<<"binN="<<Nbins_cal<<endl;
                      cout<<"exp binwidth = "<<getfeature[3]<<endl;
                      cout<<"exp std ="<<getfeature[2]<<endl;
                      printf("exp fwhm_l = %.4f exp fwhm_r = %.4f\n",getfeature[4],getfeature[5]);
                 /*     c1->cd();
                      h->Draw();
                      c1->Update();
                      c1->WaitPrimitive(); */
        }



        double lambda_old = func->GetParameter(2); // corresponding to sigma in fgaus_exp
        double gamma_old = func->GetParameter(3);  // corresponding to tR in fgaus_exp
        double delta_old = func->GetParameter(4);  // corresponding to tL in fgaus_exp


        if(fit_by_free_pars){ 
            if(fitname == "john"){ func->SetParameters(getfeature[0],getfeature[1],getfeature[2]);
             // cout<<"set john: "<<getfeature[0]<<" , "<<getfeature[1]<<" , "<<getfeature[2]<<endl;
            }
            else if(fitname == "fgaus_exp"){
             func->SetParameters(getfeature[0],getfeature[1],getfeature[2],getfeature[2],-getfeature[2]);
            // cout<<"set fgaus_exp: "<<func->GetParameter(0)<<" , "<<func->GetParameter(1)<<" , "<<func->GetParameter(2)<<endl;
             //return 0;
           }
        }  
        else{
              func->SetParameter(0,getfeature[0]);
              func->SetParameter(1,getfeature[1]);
                for(int ipar=2;ipar<func->GetNpar();ipar++){
                    func->FixParameter(ipar,func->GetParameter(ipar));
                }
        }



        double fit_rangeL, fit_rangeR;

        if(fit_by_free_pars){
             fit_rangeL = (getfeature[1] - (getfeature[2]*2.36)*0.5*2.5 );   //0.5* FWHM * 1.5
             fit_rangeR = (getfeature[1] + (getfeature[2]*2.36)*0.5*3 );
        }
        else{
             fit_rangeL = getfeature[1] - fit_half_rangeL;
             fit_rangeR = getfeature[1] + fit_half_rangeR;
        }


        fit_rangeL = (fit_rangeL>active_lowx)? fit_rangeL : active_lowx;
        fit_rangeR = (fit_rangeR < active_highx)? fit_rangeR : active_highx;




        for(int i=0;i<20;i++){
           if(fit_by_free_pars){
                   h->Fit(func,"LMQN","",fit_rangeL,fit_rangeR);
                   continue;
           }
           else{
              if(i<5){  h->Fit(func,"LMQN","",fit_rangeL,fit_rangeR);}
              else if(i>=5 && i<=10){
                if(i%2==0)   h->Fit(func,"LMQN","",fit_rangeL,fit_rangeR);
                else h->Fit(func,"LMQN","",fit_rangeL,fit_rangeR);
              }
              else h->Fit(func,"LMQN","",fit_rangeL,fit_rangeR);
           }


             fit_rangeL = func->GetParameter(1) - fit_half_rangeL;
             fit_rangeR = func->GetParameter(1) + fit_half_rangeR;
             fit_rangeL = (fit_rangeL>active_lowx)? fit_rangeL : active_lowx;
             fit_rangeR = (fit_rangeR < active_highx)? fit_rangeR : active_highx;
            

        }




       // n_slide_central = func->GetMaximumX(fit_rangeL,fit_rangeR,1e-13,500);
        double MaxY=0;

        if(fitname == "fgaus_exp"){
           MaxY = func->GetParameter(0);
           n_slide_central = func->GetParameter(1);
        }
        else{
           n_slide_central = func->GetMaximumX(fit_rangeL,fit_rangeR,1e-13,500);
           MaxY = func->GetMaximum(fit_rangeL,fit_rangeR,1e-13,500);
        }
        
        double func_FWHM_L = func->GetX(MaxY*0.5,n_slide_central-10000,n_slide_central);
        double func_FWHM_R = func->GetX(MaxY*0.5,n_slide_central,n_slide_central+10000);


        bool condition1 = (lambda_old==0)? true : (TMath::Abs((func->GetParameter(2)-lambda_old)/lambda_old) < 4);  // first fit: lambda_old==0
        bool condition2 = (gamma_old==1)? true : (TMath::Abs((func->GetParameter(3)-gamma_old)/gamma_old) < 9);
        bool condition3 = (delta_old==-1)? true : (TMath::Abs((func->GetParameter(4)-delta_old)/delta_old) < 4);
        bool condition4 = TMath::Abs(n_slide_central - getfeature[1]) <getfeature[2] * 2.35;  // check peak position
        bool condition5 = TMath::Abs(  (MaxY -getfeature[0])/TMath::Sqrt(getfeature[0])  ) < 3;  // check the amplitude of peak of func and histogram
        bool condition6 =  (MaxY > 1); // peak heigh should be >1
        bool condition7 = TMath::Abs(n_slide_central - getfeature[1]) < (func_FWHM_R-func_FWHM_L);

        if(fitname == "fgaus_exp"){
          condition2 = (TMath::Abs((func->GetParameter(3)-gamma_old)/gamma_old) < 15);
          condition3 = (TMath::Abs((func->GetParameter(4)-delta_old)/delta_old) < 100);
          if(lambda_old == Sigma) condition1 = true;
          if(gamma_old == Sigma) condition2 = true;
          if(delta_old == -Sigma) condition3 = true;
          //cout<<"here"<<endl;
        }



if(showgraph && showfitting){//show fit for each slice

cout<<endl;
cout<<endl;
cout<<"fit func MaxY= "<<MaxY<<" , FWHM= "<<(func_FWHM_R-func_FWHM_L)<<endl;

cout<<"histo heigh "<<getfeature[0]<<" Amp difference = "<<TMath::Abs(  (MaxY -getfeature[0]))<<endl;
  
  if(condition7) cout<<"pos in"<<endl;
  else cout<<"\e[1;33m"<<"out"<<"\e[0m"<<endl;
  if(condition5) cout<<"heigh match"<<endl;
  else cout<<"\e[1;33m"<<"heigh out"<<"\e[0m"<<endl;

cout<<endl;

c1->cd();
h->Draw();
h->Fit(func,"LM","",fit_rangeL,fit_rangeR);
h->SetTitle(Form("Event range %lld ~ %lld",min_evt,max_evt));
c1->Modified();
c1->Update();
char tempc;
while(1){
  cout<<" 'y' to next fit; 'n' to abort "<<endl;
  cin>>tempc;
  if(tempc!='n' && tempc!='y') continue;
  else{
    if(tempc=='n'){
      abort_corr_loop = true; 
      delete h;
      t1->ResetBranchAddresses();
      delete tmp_time;
      delete tmp_sweeps_global;
      delete tmp_tag;
      delete tmp_channel;
      return -10;
    } 
    if(tempc=='y')break;
  }
}
}


        if((condition1 && condition2 && condition3 && condition4 && condition5 && condition6 && condition7) || !fit_by_free_pars){ // fit successful
                active_lowx = n_slide_central-Half_hiswidth;
                active_highx = n_slide_central+Half_hiswidth;
                appro_lasttime = false;
                if(stop_point>0)N_byfit+=2;
                else N_byfit++;

              delete h;
              t1->ResetBranchAddresses();
              delete tmp_time;
      delete tmp_sweeps_global;
      delete tmp_tag;
      delete tmp_channel;
              return n_slide_central;
        }
        else{ // fit fail


if(showgraph && showfitting) cout<<"\e[1;33m"<<"bad fit"<<"\e[0m"<<endl; // print when showing each fit
if(showbadfit){
  printf("%.4f,%.4f,%.4f\n",lambda_old,gamma_old,delta_old);
  printf("%d , %d, %d, %d, %d, %d, %d\n",condition1 , condition2 , condition3 , condition4 , condition5 , condition6 , condition7);
  //if(c1==NULL){  c1 = new TCanvas("c1","c1",1000,800); }
 // else if(c1->GetCanvasImp()==NULL){ c1 = new TCanvas("c1","c1",1000,800); }
  c1->cd();
  h->Draw();
  h->Fit(func,"LM","",fit_rangeL,fit_rangeR);
  h->SetTitle(Form("Event range %lld ~ %lld",min_evt,max_evt));
  c1->Modified();
  c1->Update();
  char ctem;
  cin>>ctem;
}

            if(stop_point>0 || t1->GetEntriesFast()==max_evt){ // already fail last slice, ==> approximate by mean or last slice

                  if(fit_by_free_pars){ h->GetXaxis()->SetRangeUser(fit_rangeL+(getfeature[2]*2.36)*0.5*1.5,fit_rangeR-(getfeature[2]*2.36)*0.5*1.5);}
                  else{  h->GetXaxis()->SetRangeUser(fit_rangeL,fit_rangeR);}
                  n_slide_central = h->GetMean();

                  active_lowx = n_slide_central-Half_hiswidth;
                  active_highx = n_slide_central+Half_hiswidth;
                  if(fit_by_free_pars){
                      if(appro_lasttime){
                         if(fitname == "john") func->SetParameters(0,0,0,1,-1); //(0,0,0,0,-1)
                         else if(fitname == "fgaus_exp") {func->SetParameters(0,0,Sigma,Sigma,-Sigma); cout<<"reset"<<endl;}
                      }  // reset to initial pars if fitting of continuous slices are both failed
                      else{ func->SetParameters(0,0,lambda_old,gamma_old,delta_old);}
                  }
                  appro_lasttime = true;
                  N_byappro++;

                  delete h;
                  t1->ResetBranchAddresses();
                  delete tmp_time;
                  delete tmp_sweeps_global;
                  delete tmp_tag;
                  delete tmp_channel;
                  return n_slide_central;

            }
            else{ // fail current lice ==> accumulate to next slice

                if(fit_by_free_pars){
                  if(appro_lasttime){
                     if(fitname == "john") func->SetParameters(0,0,0,1,-1);
                     else if(fitname == "fgaus_exp") {func->SetParameters(0,0,Sigma,Sigma,-Sigma); cout<<"reset"<<endl;}
                  }  // reset to initial pars
                  else{ func->SetParameters(0,0,lambda_old,gamma_old,delta_old);}
                }
                  delete h;
                  t1->ResetBranchAddresses();
                  delete tmp_time;
                  delete tmp_sweeps_global;
                  delete tmp_tag;
                  delete tmp_channel;
                  return -10;
            }

        }
}




double* FeaturePar(TH1D* _inhisto,Long64_t _min_evt, Long64_t _max_evt){
    static double par_return[7];  // Peak height , peak center, sigma, good bin size, FWHM_L, FWHM_R, integral counts
    TH1D* hc = (TH1D*) _inhisto->Clone();
    par_return[6] = hc->Integral(1,hc->GetNbinsX());
    if(par_return[6]<(_max_evt-_min_evt)*0.5*0.1){ // empty slice
       delete hc;
       return par_return;
    }

    hc->Smooth(4);

    int Bin_max = hc->GetMaximumBin();
    par_return[0] = hc->GetBinContent(Bin_max);  // Amp
    par_return[1] = hc->GetBinCenter(Bin_max);   // peak center
    double height_half = par_return[0]*0.5;

    double bin_i_y1;
    double bin_i_y2;
    double candidate_x;
    double FWHM_L=0,FWHM_R=0;

    for(int i=1;i<hc->GetNbinsX();i++){
      bin_i_y1 = hc->GetBinContent(i);
      bin_i_y2 = hc->GetBinContent(i+1);
      candidate_x = (height_half - bin_i_y1) / (bin_i_y2 - bin_i_y1) * hc->GetBinWidth(1) + hc->GetBinCenter(i);
      if(i<Bin_max){
          if(bin_i_y1<height_half && bin_i_y2 >height_half){  FWHM_L = candidate_x; }
      }
      else{
          if(bin_i_y1>height_half && bin_i_y2 <height_half){    FWHM_R = candidate_x; break;  }
      }
    }

    double FWHM = FWHM_R - FWHM_L;

    hc->GetXaxis()->SetRangeUser(FWHM_L-1.5*FWHM,FWHM_R+1.5*FWHM);

    par_return[2] = hc->GetStdDev();

    _inhisto->GetXaxis()->SetRangeUser(FWHM_L-1.5*FWHM,FWHM_R+1.5*FWHM);//FWHM_L-0.25*FWHM,FWHM_R+0.25*FWHM

    double maincounts = _inhisto->Integral();

    //scott rule for bin size:  bin size= 3.49*sigma*(counts)^(-1/3)

    par_return[3] = 3.49 * par_return[2] *TMath::Power(maincounts,-1./3.);

    _inhisto->GetXaxis()->UnZoom();

par_return[4]=FWHM_L;
par_return[5]=FWHM_R;

    delete hc;

    return par_return;

}




#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc>1){
    DriftCorrect( string(argv[1]) , string(argv[2]),atof(argv[3]),atoi(argv[4]),atoi(argv[5]),atof(argv[6]),atof(argv[7]),atof(argv[8]));
  }
  return 0;
}
#endif







