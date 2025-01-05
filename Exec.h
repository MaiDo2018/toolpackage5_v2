#ifdef _PREVIEWER_
#ifndef _EXEC_H_
#define _EXEC_H_


#include "TSystem.h"
#include "TGClient.h"
#include "Buttons.h"
#include "TObject.h"
#include "TFrame.h"
#include "TVirtualX.h"
#include <stack>

TCanvas* c_handle;
bool stack2clear[2] = {1,1};
// fit setting
string AFitOption= "LMEQ";
int AFitSetFunc=2;
char AFit_x_r = 'x';
bool FixFitRange=false;
double fitL_width=0;  // fix fit width for normal fit
double fitR_width=0;

double unbinned_fitL=0;  // fit left edge tof value
double unbinned_fitR=0; // fit right edge tof value

// detector efficiency measurement
TH1D* effigraph = NULL; 
TCanvas* c_effi =NULL;
TF1* f_N_b_tof = NULL;
bool UseRejectionList=false;
double* timewidth = NULL;
double* N_b_tof = NULL;


void SetRange(int index,TCanvas* c_get,TH1D* h,double xmin, double xmax);
void SetLogy(TPad* p,bool flag);
void MakerAdjY(TH1D* h_in, TOFMarker* marker,double h_left,double h_right);

//********** for sampling fitting ********************
bool Sampling(int padindex, double _rangeL, double _rangeR);

void Exec(int index=1){
    if (!gPad) {
      Error("exec", "gPad is null, you are not supposed to run this macro");
      return;
    }  

    static double xmouse[2] = {0};
    static int imouse = 0;
    static TH1D* h_history[2] = {};
    static stack<double> lowlimit[2];
    static stack<double> highlimit[2];

    c_handle = (TCanvas*)gPad->GetCanvas();
    c_handle->FeedbackMode(kTRUE);


    //TObject *Sele = gPad->GetSelected();
   //if(!select) return;
   //if (!select->InheritsFrom(TH1::Class())) {gPad->SetUniqueID(0); return;}
/*
   TH1D * h_handle = (TH1D*)Sele;
   int index=0;
   const char* getname = h_handle->GetName();

   if(strcmp(getname,"h_xF") ==0){index=0;}
   else if(strcmp(getname,"h_zoom_x") ==0){index=1;}
   else if(strcmp(getname,"h_refF") ==0){index=2;}
   else return;*/

    
  TPad* Pd = (TPad*)c_handle->cd(index);

  //cout<<"\r"<<"index = "<<index<<flush;



    // new histogram setting, clear stack;
    if(index==1){
        if(h_history[index-1] != h_xF){
          h_history[index-1] = h_xF;
          stack2clear[index-1] = true;
        }
        else stack2clear[index-1] = false;

    }
    else if(index ==2){
        if(Pd->GetPrimitive("h_zoom_x")!=NULL){   // case of h_zoom_x
            if(h_history[index-1] != h_zoom_x){
              h_history[index-1] = h_zoom_x;
              stack2clear[index-1] = true;
            }
            else stack2clear[index-1] = false;
            
        }// case of h_refS
        else if(Pd->GetPrimitive("h_refS")!=NULL){
            if(h_history[index-1] != h_refS){
              h_history[index-1] = h_refS;
              stack2clear[index-1] = true;
            }
            else stack2clear[index-1] = false;
        }
        else{
            Error("exec", "No histogram found in c1->cd(2)");
            return;
        }
    }
    else{
        Error("exec", "index =1 or 2; exec program only limit to h_xF and h_zoom_x!!");
        return;
    }


    if(stack2clear[index-1]){
      for(unsigned int i=0;i<lowlimit[index-1].size();i++){
        lowlimit[index-1].pop();
        highlimit[index-1].pop();
      }
      stack2clear[index-1]=false;
    }





  //TPad* Pd = (TPad*) gPad->GetPad(index);

  int pxold = Pd->GetUniqueID();
  int px =  Pd->GetEventX();
  int py =  Pd->GetEventY();

  float uxmin = Pd->GetUxmin();
  float uxmax = Pd->GetUxmax();
  int pxmin = Pd->XtoAbsPixel(uxmin);
  int pxmax = Pd->XtoAbsPixel(uxmax);
  float uymin = Pd->GetUymin();
  float uymax = Pd->GetUymax();
  int pymin = Pd->YtoAbsPixel(uymin);
  int pymax = Pd->YtoAbsPixel(uymax);

  if(px<pxmin) px = pxmin+1;
  if(px>pxmax) px = pxmax-1;
  
  if(pxold) gVirtualX->DrawLine(pxold,pymin,pxold,pymax);
  gVirtualX->DrawLine(px,pymin,px,pymax);
  Pd->SetUniqueID(px);
  Double_t upx =  Pd->AbsPixeltoX(px);
  Double_t x =  Pd->PadtoX(upx);

  EEventType event = static_cast<EEventType> ( Pd->GetEvent());
  if(event == kButton1Down){ // select range
    //cout<<"\r x: "<<Form("%.2f",upx)<<flush;
    xmouse[imouse] = upx;  imouse = (imouse+1)%2;
   /* TObject *select = gPad->GetSelected();
   TH1D * h_handle = (TH1D*)select;
   int index=0;
   const char* getname = h_handle->GetName();

   if(strcmp(getname,"h_xF") ==0){index=0;}
   else if(strcmp(getname,"h_zoom_x") ==0){index=1;}
   else if(strcmp(getname,"h_refF") ==0){index=2;}
   else return;*/
    
    gVirtualX->DrawLine(px,uymin,px,uymax);

   // printf("px=%d , py=%d, upx=%f, x=%f\n",px,py,upx,x);
  }
  else if(event == kButton1Double){
    printf("tof position = %.3f\n\n",x);
    gXposition = x;   //SetROI use

    if(m_ref==0){cout<<"No ref mass setting, no mass calculation"<<endl;}
    else if(tof_ref_cento <1){cout<<"No ref tof setting, no mass calculation"<<endl;}
    else{
      double tem_tof = tof_x_cento[0];
      double tem_tof_err = tof_x_cento_err[0];
      tof_x_cento[0] = x;
      tof_x_cento_err[0]=0.7;
      mass_calculator(1);
      tof_x_cento[0] =tem_tof;
      tof_x_cento_err[0] = tem_tof_err;
      
    }
  }
  else if(event == kKeyPress){ // Key action
      int press = Pd->GetEventX(); //cout<<"press: "<<press<<endl;
      if(press==' '){ // set new range
        double Tmin = TMath::Min(xmouse[0],xmouse[1]);
        double Tmax = TMath::Max(xmouse[0],xmouse[1]);
        histo_zoom_in_x(0,100,Tmin,Tmax);
        printf("used: histo_zoom_in_x(%d,%d,%.1f,%.1f)\n",0,100,Tmin,Tmax);
        c_handle->cd(2)->Modified();
        c_handle->cd(2)->Update();      
        c_handle->Modified();  c_handle->Update();
        
      }else if(press=='z'){ // zoom
        if(h_history[index-1] ==NULL) return;
        lowlimit[index-1].push(TMath::Min(xmouse[0],xmouse[1]));
        highlimit[index-1].push(TMath::Max(xmouse[0],xmouse[1]));
        SetRange(index,c_handle,h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top()); 
        if(index==1){ //"h_xF"
            MakerAdjY(h_history[index-1],marker_tof,lowlimit[index-1].top(),highlimit[index-1].top());
            ROIadjY(h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top());
            c_handle->cd(1)->Modified();
            c_handle->cd(1)->Update();
            c_handle->Modified();  c_handle->Update();
        }     
      }else if(press=='x'){ // unzoom
        if(h_history[index-1] ==NULL) return;
        if(lowlimit[index-1].size() !=0){
              lowlimit[index-1].pop();
              highlimit[index-1].pop();
        }
        if(lowlimit[index-1].size() ==0){ 
              h_history[index-1]->GetXaxis()->UnZoom(); Pd->Modified();Pd->Update();
              double histo_L = EJE0-300;
              double histo_H = histo_L + sptrFW;
              MakerAdjY(h_history[index-1],marker_tof,histo_L,histo_H);
              ROIadjY(h_history[index-1],histo_L,histo_H);
              c_handle->cd(1)->Modified();
              c_handle->cd(1)->Modified();
              c_handle->cd(1)->Update();
              c_handle->Modified();  c_handle->Update();
        }
        else{
             SetRange(index,c_handle,h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top()); 
              if(index==1){ //"h_xF"
                  MakerAdjY(h_history[index-1],marker_tof,lowlimit[index-1].top(),highlimit[index-1].top());
                  ROIadjY(h_history[index-1],lowlimit[index-1].top(),highlimit[index-1].top());
                  c_handle->cd(1)->Modified();
                  c_handle->cd(1)->Modified();
                  c_handle->cd(1)->Update();
                  c_handle->Modified();  c_handle->Update();
              } 
        }
      }else if(press=='l'){ // logy
        SetLogy(Pd,!Pd->GetLogy());
        c_handle->Modified(); c_handle->Update();
      }else if(press=='r'){// add new ROI
        if(ROI_initial){
                  ROI_INDEX++;
                  SetROI(gXposition,ROI_WIDTH,(ROI_INDEX-1)%40+1,0,true,0);
                  c_handle->Modified(); c_handle->Update();
                  //printf("%.4f,%f,%d,%d",gXposition,ROI_WIDTH,ROI_INDEX,(ROI_INDEX-1)%20+1); // test purpose
        }
        else{
                  SetROI(gXposition,ROI_WIDTH,ROI_INDEX,0,true,0);
                  c_handle->Modified(); c_handle->Update();
                  //printf("%.4f,%f,%d,%d",gXposition,ROI_WIDTH,ROI_INDEX,(ROI_INDEX-1)%20+1); // test purpose
                  ROI_initial=true;
        }
      }else if(press=='R'){// correct current ROI
                  SetROI(gXposition,ROI_WIDTH,(ROI_INDEX-1)%40+1,0,true,0);
                  c_handle->Modified(); c_handle->Update();
      }else if(press=='s'){
                if(h_zoom_ref==NULL){cout<<"Warning!!! h_zoom_ref is not exist!!!! Ahort!!!"<<endl; return;}
                cout<<"Sampling at pad(4)............"<<endl;
                 if(Sampling(4,-1,-1))cout<<"Successful sampling"<<endl; // sample at tag1 , have to draw arrow to define sampling range

      }else if(press=='S'){
                if(h_zoom_x==NULL){cout<<"Warning!!! h_zoom_x is not exist!!!! Ahort!!!"<<endl; return;}
                cout<<"Sampling at pad(2)............"<<endl;
                if(Sampling(2,TMath::Min(xmouse[0],xmouse[1]),TMath::Max(xmouse[0],xmouse[1])))cout<<"Successful sampling"<<endl; // sample at tag1 , have to draw arrow to define sampling range

      }else if(press=='e'){
                if(h_zoom_ref==NULL){cout<<"Warning!!! h_zooom_ref is not exist!!!! Ahort!!!"<<endl; return;}
                cout<<"Draw arrow to set the range for fitting of reference......"<<endl;
                //c_handle->cd(4);
                get_para_by_draw(4);
                fext->SetNorder(funcN::Norder2Set); // update norder before "Initial_par_ref()"
                fext->Initial_par_ref(h_zoom_ref, fitrangeL , fitrangeR);                
                fext->MakeFitFunc_ref();
                fext->Fit_ref(h_zoom_ref, fitrangeL , fitrangeR,AFitOption);
                tem_func = fext->GetFuncRef(); //c_handle->cd(4);tem_func->Draw("same");// printChi() uses global handle of func
                printChi(h_zoom_ref,'r',fitrangeL , fitrangeR);
                tof_ref_cento = fext->GetsPeakCenter();
                tof_ref_cento_err = fext->GetsPeakCenter_err();
                fext->MakeFitFunc_N();
                tem_func = fext->GetFuncX(); //make the current fitting recognised
                cout<<endl;
                c_handle->Update(); c_handle->Modified();
                c_handle->cd(2);
                


      }else if(press=='E'){
                if(h_zoom_x==NULL){cout<<"Warning!!! h_zoom_x is not exist!!!! Ahort!!!"<<endl; return;}
                c_handle->cd(2)->SetEditable(kTRUE);

                double Tmin = TMath::Min(xmouse[0],xmouse[1]);
                double Tmax = TMath::Max(xmouse[0],xmouse[1]);
                if(Tmax-Tmin<0.5) return; // fitting range should >0.5 ns for normal fitting

                fext->SetNorder(funcN::Norder2Set);
                fext->Initial_par_ref(h_zoom_x, Tmin ,Tmax );
                fext->MakeFitFunc_ref();
                fext->Fit_ref(h_zoom_x, Tmin ,Tmax, AFitOption);
                tem_func = fext->GetFuncRef(); //c_handle->cd(2);tem_func->Draw("same");
                printChi(h_zoom_x,'x',Tmin ,Tmax);
                tof_ref_cento = fext->GetsPeakCenter();
                tof_ref_cento_err = fext->GetsPeakCenter_err();
                fext->MakeFitFunc_N();
                tem_func = fext->GetFuncX();
                c_handle->Update(); c_handle->Modified();
                c_handle->cd(2)->SetEditable(kFALSE);
        
      }else if(press=='f'){// fitquickly
              if(tem_func == NULL) cout<<"No fitting curve Available!!!!!!!!! Abort!!"<<endl;
              else tem_func->SetLineColor(kRed); //unbinned fitting change color to blue;

              c_handle->cd(2);
              if(index==2){
                  double Max_high =0;  // maximum count in selected range
                  double Max_tof=0;    // tof of maximum bin
                  double Tmin = TMath::Min(xmouse[0],xmouse[1]);
                  double Tmax = TMath::Max(xmouse[0],xmouse[1]);
                  int Bin_L = h_zoom_x->FindBin(Tmin);
                  int Bin_R = h_zoom_x->FindBin(Tmax);

                  for(int i=Bin_L;i<=Bin_R;i++){ // find the position and count of maximum bin in the range
                    if( h_zoom_x->GetBinContent(i) > Max_high ){Max_high = h_zoom_x->GetBinContent(i); Max_tof=h_zoom_x->GetBinCenter(i);}
                  }

                  tem_func->SetParameter(0,Max_high);
                  tem_func->SetParameter(1,Max_tof);

                  if(FixFitRange && fitL_width!=0 && fitR_width!=0){ // fix fit range
                      Tmin = Max_tof-fitL_width;
                      Tmax = Max_tof+fitR_width;
                  }
                  else{ // free fix fit range
                      fitL_width = Max_tof - Tmin;  
                      fitR_width = Tmax - Max_tof;
                  }


                  string tem_func_name = tem_func->GetName();
                  cout<<"tem_func_anme = "<<tem_func_name<<endl;


                  //*************** fit with sampling function ******************
                  if(tem_func_name == "fsample"){
                      if(fs->GetOldNPeaks() != fs->NumOfPeaks){fs->Makefitfunc();}

                       if(NumOfPeaks==1){ // one peak
                          fs->SetPars(1,Max_high,Max_tof);
                       }
                       else{// more than one peak

                          for(int i=0;i<NumOfPeaks;i++){
                            cout<<"draw an arrow from top of "<<"\e[1;33m"<<"Peak_"<<i+1<<"\e[0m"<<" to FWHM"<<endl;
                            cout<<"pending......"<<endl;
                            get_para_by_draw(2);
                            fs->SetPars(i+1,tem_high,tem_cento);
                            
                          }

                       }

                      fs->Fit(h_zoom_x,Tmin,Tmax,AFitOption);
                     // tem_func = fs->Getfitfunc();
                      c_handle->Modified(); c_handle->Update();
                      

                      for(int i=0;i<NumOfPeaks;i++){
                           tof_x_cento[i] = fs->GetTofCenter(i+1); // peakindex from 1
                           tof_x_cento_err[i] = fs->GetTofCenterErr(i+1);
                           printf("cento_%d: %.4f(%.4f)\n",i+1,tof_x_cento[i],tof_x_cento_err[i]);
                           unbinned_fitR = tof_x_cento[i] + fs->range_R; // update fit right edge for unbinned fit
                      }

                          unbinned_fitL = tof_x_cento[0] - fs->range_L;
                          unbinned_fitR = tof_x_cento[NumOfPeaks-1] + fs->range_R; // update fit right edge for unbinned fit

                      printChi(h_zoom_x,'x',unbinned_fitL,unbinned_fitR);

                      if(fs->NumOfPeaks>1){  fs->Draw_subline(c1,2);  }

                  }// end of sampling fit
                  else if(tem_func_name == "fextend_N"){ //extend gaus-exp fitting
                      c_handle->cd(2)->SetEditable(kTRUE);

                      if(NumOfPeaks==1){ // one peak
                          fext->SetPars(1,Max_high,Max_tof);
                       }
                       else{// more than one peak

                          for(int i=0;i<NumOfPeaks;i++){
                            cout<<"draw an arrow from top of "<<"\e[1;33m"<<"Peak_"<<i+1<<"\e[0m"<<" to FWHM"<<endl;
                            cout<<"pending......"<<endl;
                            get_para_by_draw(2);
                            fext->SetPars(i+1,tem_high,tem_cento);
                            
                          }

                       }

                       if(fext->GetNumberOfPeaks() != funcN::NPs2Set){
                            if(fext->SetNumberOfPeaks(funcN::NPs2Set)) fext->MakeFitFunc_N();
                            else return;
                        }

                       fext->Fit_N(h_zoom_x,Tmin,Tmax,AFitOption);                      
                       c_handle->Modified(); c_handle->Update();
                       tem_func = fext->GetFuncX();
                       

                      for(int i=0;i<NumOfPeaks;i++){
                           tof_x_cento[i] = fext->GetTofCenter(i+1); // peakindex from 1
                           tof_x_cento_err[i] = fext->GetTofCenterErr(i+1);
                           printf("cento_%d: %.4f(%.4f)\n",i+1,tof_x_cento[i],tof_x_cento_err[i]);
                           //unbinned_fitR = tof_x_cento[i] + fext->range_R; // update fit right edge for unbinned fit
                      }

                          unbinned_fitL = tof_x_cento[0] - fext->sample_range_L;
                          unbinned_fitR = tof_x_cento[NumOfPeaks-1] + fext->sample_range_R; // update fit right edge for unbinned fit

                      printChi(h_zoom_x,'x',unbinned_fitL,unbinned_fitR);

                      if(fext->GetNumberOfPeaks()>1){  fext->Draw_subline(c1,2);  }
                  }
                  else{//******************** fit with gaus_exp function *********************
                      fitquickly(AFit_x_r,AFitSetFunc,AFitOption.c_str(),true,Tmin,Tmax,50);
                      unbinned_fitL = Tmin;
                      unbinned_fitR = Tmax;
                      c_handle->Modified(); c_handle->Update();
                      printf("using: fitquickly('%c',%d,\"%s\",true,%.4f,%.4f,50)\n",AFit_x_r,AFitSetFunc,AFitOption.c_str(),Tmin,Tmax);  
                  }// end of gaus_exp fit                
              }
              else cout<<"Only canvas 2 has fitting function"<<endl;

      }
     else if(press=='-'){
             double Tmin = TMath::Min(xmouse[0],xmouse[1]);
             double Tmax = TMath::Max(xmouse[0],xmouse[1]);
          	 printf("x1=%.4f;\t x2=%.4f\n",Tmin,Tmax);
          	 printf("distance:%.4f\n",Tmax-Tmin);
      }
      else if(press =='c'){

          gROOT->ProcessLine(".x BetaCoin.confg");

          if(MakeBetaRawTree && BetaFastMode){
                  cout<<endl;
                  cout<<endl;
                  cout<<"\e[1;33m MakeBetaRawTree option is invalid under \"BetaFastMode\". By default, BetaFastMode is set to false to continue.\e[0m"<<endl;
                  BetaFastMode = false;
                  sleep(2);
          }

          static int PreSetADC1_low_Thres_ch_old=-1;
          static int PreSetADC2_low_Thres_ch_old=-1;
          static bool MakeBetaRawTree_old=false;

          static bool BetaFastMode_last_time=false;
          bool SwitchingFastMode2Normal=false;
          if(BetaFastMode_last_time && !BetaFastMode) SwitchingFastMode2Normal = true;
          BetaFastMode_last_time = BetaFastMode;

          printf("Searching for beta - TOF coincidence .......\n\n");
          ShowCoinCondition();
          ShowRIHalflife();
          if(ReadBetaOnline){
              Path_beta_file=FilePath+"../LST/";
              cout<<"\e[1;33m"<<"Loading beta file online"<<"\e[0m"<<endl;
          }
          else{
              //%%%%%%%%%%%%%%%%%%  the Path_beta_file will be set to FilePath+"../LST/" in Read_beta_lst()
              cout<<"\e[1;33m"<<"Loading beta file from local."<<"\e[0m"<<endl;
          }

          

          //%%%%%%%%%%%%%%%%              only handle Normal reading mode, Fast mode will skip this part,                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          //%%%%%%%%%%%%%%%% reading beta.lst and do beta-beta and beta-tof analysis in FindBetaTof_Coin() function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if(Current_beta_file != ToLoad_beta_file || SwitchingFastMode2Normal || PreSetADC1_low_Thres_ch_old != PreSetADC1_low_Thres_ch || PreSetADC2_low_Thres_ch_old != PreSetADC2_low_Thres_ch
                || MakeBetaRawTree_old != MakeBetaRawTree){  
                                                                                                      // beta file read to vector
                                                                                                      // Current_beta_file depends on the last time reading a file with "loadmore = false" in Read_Beta_lst(), All raw data
                                                                                                      // in vector is cleared. Therefore:
                                                                                                      // After runing Normal mode:  Current_beta_file = the first file in the batch of beta.lst
                                                                                                      // After BetaFastMode:  Current_beta_file = the first file in the batch of beta.lst that covering the global of the last Tof event
                                                                                                      // SwitchingFastMode2Normal --> make sure to read all raw data when reading mode is switched from BetaFastMode to normal mode
                                                                                                      // To get rid of a case: only one Tof event happen during the first beta.lst. Leading Current_beta_file == ToLoad_beta_file after
                                                                                                      // executing Reading beta.lst process by BetaFastMode last time.
            if(ToLoad_beta_file=="+++"){
              cout<<"warning: no beta file to be load; please load file first"<<endl;
              cout<<"\e[1;33m"<<"beta filename with .lst"<<"\e[0m"<<endl;
              cout<<"using:  LoadNewBetaFile("<<endl;
              return;
            }
            else{// really going to read new beta file or reload beta.lst because of changing threshold
              //  Path_beta_file = FilePath;
                if(!BetaFastMode){
                  bool RunStatus=false;  
                  if(ReadBetaOnline) RunStatus = Read_Beta_lst_batch(OnlineBetaFilePath,ToLoad_beta_file,PreSetADC1_low_Thres_ch,PreSetADC2_low_Thres_ch); //load all data as default//"/mnt/beta_driver/F11_MRTOF/240604MT/usb/"
                  else RunStatus = Read_Beta_lst_batch(FilePath+"../","LST/"+ToLoad_beta_file,PreSetADC1_low_Thres_ch,PreSetADC2_low_Thres_ch);
  //              Read_Beta_lst(FilePath+"../","LST/"+ToLoad_beta_file);
   //               Read_Beta_lst_batch(FilePath+"../","LST/"+ToLoad_beta_file);

                  if(!RunStatus){
                    cout<<endl;
                    cout<<"\e[1;31m Fail to load beta list file!!!!!!!!! Abort!!! \e[0m"<<endl;
                    return;
                  }

                    PreSetADC1_low_Thres_ch_old = PreSetADC1_low_Thres_ch;
                    PreSetADC2_low_Thres_ch_old = PreSetADC2_low_Thres_ch;
                    MakeBetaRawTree_old = MakeBetaRawTree;
                }
            }
          }
          if(!BetaFastMode){ //return;
             if(Compare_Beta_Beta()){
                  cout<<"show beta beta"<<endl;
                ShowBeta_Beta_Coin();
            }
          }  // which will create tbeta_beta coincident tree, can be an option

          if(!BetaFastMode){if(Match_hit.size()==0){cout<<"\e[1;23m"<<"No coincidence in store"<<"\e[0m"<<endl; return;} }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          IsBetaFromNewFile(0,0,false,true); // clear the history of beta file have been read and compare for beta_beta coincidence; allow to have an accummulated histogram in "BeatFastMode"
          if(BetaFastMode){ // Beta_Graph_Initialize() is executed by ShowBeta_Beta_Coin() in nnrmal mode
            if(!Beta_Graph_Initialize()){
                  cout<<"\e[1;33m"<<"Fail to initialize the Beta_Beta_coin Graph beacaue of failure of creating file to save tbeta_beta tree!!"<<"\e[0m"<<endl;
                  return;
            }
          }

          double Tmin = TMath::Min(xmouse[0],xmouse[1]);
          double Tmax = TMath::Max(xmouse[0],xmouse[1]);
          FindBetaTof_Coin(Tmin,Tmax, RIHalflive,mycut); 
         // if(BetaFastMode)  ShowBeta_Beta_Coin();      
      }
      else if(press =='m'){ // mass calculate
        for(int index_m=0;index_m<NumOfPeaks && index_m<10;index_m++){
          cout<<"Unknow mass "<<index_m+1<<":"<<endl;
          mass_calculator(tof_x_cento[index_m],tof_x_cento_err[index_m],laps_x,q_x,tof_ref_cento,tof_ref_cento_err,laps_ref,m_ref,err_ref,bref,bref_err,q_ref);
          cout<<endl;
        }

      }
      else if(press==','){
          if(tem_func == NULL)return;
          string curvename = tem_func->GetName();

          if(curvename == "fsample"){
              fs->FreeRange = !(fs->FreeRange); // flag of FreeRange initialized to false
              
              if(fs->FreeRange){ cout<<"Set to FreeRange== true for sampling fitting, fitting range FREE now!!!"<<endl; FixFitRange =false;}
              else{ cout<<"Set to FreeRange== false"<<endl;
                      cout<<"fitting range fix with  same left and right tail ratio as sampling peak"<<endl;
                      FixFitRange=true;
              }
          }
          else if(curvename == "fextend_N"){
              fext->FreeRange = !(fext->FreeRange);
              if(fext->FreeRange){ cout<<"Set to FreeRange== true for sampling fitting, fitting range FREE now!!!"<<endl; FixFitRange =false;}
              else{ cout<<"Set to FreeRange== false"<<endl;
                      cout<<"fitting range fix with  same left and right tail ratio as sampling peak"<<endl;
                      FixFitRange=true;
              }
          }
          
          
      }
      else if(press=='u'){
        if(tem_func ==NULL) cout<<"No fitting curve available!! Abort!!!"<<endl;
        else tem_func->SetLineColor(kBlue);
        UnbinnedFit(tem_func,active_tree_name, h_zoom_x, unbinned_fitL,unbinned_fitR);
      }
      else if(press=='k'){//beta tof SSD efficiency scan VS time window
          
         
          

          if(effigraph!=NULL && Efficiency_histo_reset){delete effigraph;  effigraph=NULL;}

          Efficiency_histo_reset = true;
//cout<<"after delete"<<endl;


          double Tmin = TMath::Min(xmouse[0],xmouse[1]);
          double Tmax = TMath::Max(xmouse[0],xmouse[1]);
          
          TRandom3 ranran((unsigned int)time(0));
          
          double Time2halflife_old = Times2Halflive;

          const int Npoints = TMath::Nint(Multiple_max / Multiple_diff);


          if(timewidth!=NULL){delete timewidth; delete N_b_tof;}
          timewidth = new double[Npoints];
          N_b_tof = new double [Npoints];

          for(int i=0;i<Npoints;i++){
              Times2Halflive = Multiple_diff*(i+1);
              FindBetaTof_Coin(ranran.Uniform(Tmin-1,Tmin+1),ranran.Uniform(Tmax-1,Tmax+1), RIHalflive,mycut);
            if(UseRejectionList)  GenerateRejectList(tbeta_tof,mycut,&rejectlist);
              ShowDecayHisto(mycut,&rejectlist);
              timewidth[i] = Times2Halflive*RIHalflive*1e-9;
              N_b_tof[i] = (double)T_data.size();

          }

//cout<<"delet 1canvas"<<endl;
          if(c_effi!=NULL && c_effi->GetCanvasImp()!=NULL) delete c_effi; // window is closed
          c_effi = new TCanvas("c_effi","c_effi",600,400); 
//cout<<"delet canvas"<<endl;


//cout<<"effigraph = N= ; l= r="<<effigraph<<";"<<Npoints<<";"<<timewidth[0]*0.5<<";"<<timewidth[Npoints-1]+timewidth[0]*0.5<<endl;
          effigraph = new TH1D("effigraph","effigraph",Npoints,timewidth[0]*0.5,timewidth[Npoints-1]+timewidth[0]*0.5);
//cout<<"effigraph = N= ; l= r="<<effigraph<<";"<<Npoints<<";"<<timewidth[0]*0.5<<";"<<timewidth[Npoints-1]+timewidth[0]*0.5<<endl;
if(effigraph==NULL) return;


          effigraph->SetTitle("get efficiency;T window width [s];N_beta_tof counts");
          effigraph->GetXaxis()->CenterTitle();
          effigraph->GetYaxis()->CenterTitle();
          
          for(int i=1;i<=Npoints;i++){effigraph->SetBinContent(i,N_b_tof[i-1]);}


          c_effi->cd();
          effigraph->Draw("e");

          if(f_N_b_tof!=NULL) delete f_N_b_tof;

          f_N_b_tof = new TF1("f_N_b_tof",Expect_b_tof_counts,0,timewidth[Npoints-1]+Multiple_diff*RIHalflive*1e-9,4);
          f_N_b_tof->SetParameter(2,detector_effi);f_N_b_tof->SetParLimits(2,0,0.34);
          if(!fit_efficiency) f_N_b_tof->FixParameter(2,detector_effi);

          f_N_b_tof->SetParameter(0,TMath::Log(2)/(RIHalflive*1e-9));  f_N_b_tof->SetParLimits(0,0.01,70);
          if(!fit_tao)f_N_b_tof->FixParameter(0,TMath::Log(2)/(RIHalflive*1e-9));

          f_N_b_tof->SetParameter(1,Beta_bg_rate); f_N_b_tof->SetParLimits(1,0.00001,10);
          if(!fit_vb)f_N_b_tof->FixParameter(1,Beta_bg_rate);

          f_N_b_tof->FixParameter(3,(double)Ntof_counter);

          f_N_b_tof->SetParName(0,"lamda_ref");
          f_N_b_tof->SetParName(1,"vb /s");
          f_N_b_tof->SetParName(2,"efficiency");
          f_N_b_tof->SetParName(3,"Ntof");

          effigraph->Fit(f_N_b_tof,"LE");
          Times2Halflive =Time2halflife_old;

          c1->cd(2);
          cout<<endl;
          cout<<endl;

      }
            //fit_tao;fit_vb;fit_efficiency
  }
 
  

  funcS::NumOfPeaks = NumOfPeaks; // keep update to global NumOfPeaks
  funcN::NPs2Set = NumOfPeaks;
  MainPeakIndex = fext->GetMainPeakIndex();

  gClient->GetDefaultRoot();
}//end of Exec


void SetRange(int index,TCanvas* c_get,TH1D* h,double xmin, double xmax){
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  c_get->cd(index)->Modified();
  c_get->cd(index)->Update();      
  c_get->Modified();  c_get->Update();
  return;
}

void SetLogy(TPad* p,bool flag){
  p->SetLogy(flag);
  p->Modified();
  p->Update();      
  return;
}

void MakerAdjY(TH1D* h_in, TOFMarker* marker, double h_left,double h_right){
  string histo_name = h_in->GetName();
  if(histo_name == "h_xF"){
      double Y2Set = h_in->GetBinContent( h_in->GetMaximumBin() );
      if(Y2Set ==0) Y2Set=1.01;
      else Y2Set *=1.01;
      double lineX=0;
      for(int index=0;index<40;index++){
          lineX = marker[index].GetX();
          if(lineX>= h_left && lineX<=h_right){
            marker[index].SetY(Y2Set);
          }
      }
  }
}

bool Sampling(int padindex, double _rangeL, double _rangeR){
    if(fs==NULL){cout<<"sample function is no exist, Abort!!!"<<endl; return false;}
  /*  int padindex=2;
    while(1){
      cout<<"Which Pad for sampling: 2 or 4"<<endl;
      cin>>padindex;
      if(padindex==2 || padindex==4) break;
    }*/

    TH1D* getHisto=NULL;
    if(padindex==2){getHisto=h_zoom_x;}
    else if(padindex==4){
          getHisto=h_zoom_ref;
          cout<<"draw a line for sampling range"<<endl;
          cout<<"pending..."<<endl;
          get_para_by_draw(padindex);
          _rangeL = fitrangeL;
          _rangeR = fitrangeR;
    }
    else{cout<<"Error: padindex!!!"<<endl; return false;}

    fs->Sampling(getHisto,_rangeL,_rangeR);
    fs->RecreateHisto();

    if(fs->SmoothHisto()){
        fs->Makefitfunc();
        tem_func=fs->Getfitfunc(); // get handle of fs
        tof_ref_cento = fs->GetsPeakCenter();
        tof_ref_cento_err = fs->GetsPeakCenter_err();
    }
    else{cout<<"Error in Smooth histogram!!! Abort!!!"<<endl; return false;}

    if(fs->Draw(c1,padindex)){c1->cd(2);}
    else{cout<<"Error, faile to draw sampling line on histogram!!!"<<endl; return false;}

    return true;

}


#endif  // ifndef _Exec_h_

#endif// ifdef _PREVIEWER_
