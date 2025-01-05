//*************************************************************************
//* Offline decoder
//*
//*************************************************************************

//*****************************************************
// note: channel value = channel number directly
//    0 = ch 0 (not exist)
//    1 =  ch 1
//    2 = ch 2
//    3 = ch 3
//    ...... (ch 5 is maximum)
//    6 =  ch 6 (start channel)
//******************************************************
#include <vector>
#include <TRint.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TSystem.h>

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

//define _LASER_   // switch to old verson(event by sweeps) by disable this line
#define LaserFequence 0.774

#if defined (MAKECINT)
//#pragma extra include "vector";
#pragma link C++ class std::vector<long>+;
#endif 

using namespace std;
/*
int sweepsMax=0; // for using in Online.C to known the maximum sweeps number allowed;
long gsweeps_tof=0;
int gsweeps_old_tof=-1;
bool HasData=false;*/

class EventContainer{
  public:
      vector <Int_t> *channel;    
      vector <Int_t> *edge;    
      vector <long> *value;   
      vector <Double_t> *time;  
      vector <Int_t> *sweeps;    
      vector <long> * sweeps_global;  
      vector <Int_t> *tag;   
      vector <Int_t> *lost;
      vector <Int_t> *glo;   

      Long64_t evt_index = 0; // record the index of current event in the whole file, determine by outside
      int nhits[10]; // hits at each channel

      unsigned long num_entry;

      EventContainer(){ // initialize
        for(int i=0; i<10;i++){ nhits[i]=0; }
        
        num_entry=0;

        channel = new vector <Int_t>();    
        edge = new vector <Int_t>();    
        value = new vector <long>();   
        time = new vector <Double_t>();  
        sweeps = new vector <Int_t>();    
        sweeps_global = new vector <long>;  
        tag = new vector <Int_t>();   
        lost = new vector <Int_t>();
        glo = new vector <Int_t>();   

      }

      void GetStats(){
        num_entry=sweeps->size();

        for(int i=0; i<10;i++){ nhits[i]=0; }

        for(unsigned long ihit=0;ihit<num_entry;ihit++){
          nhits[channel->at(ihit)]++;
        }
      }

      unsigned long GetNumEntry(){
        GetStats();
        return num_entry;
      }


      void Sort(){

            if(GetNumEntry()==0) return;

            vector <Int_t> channel_tem = *channel;    
            vector <Int_t> edge_tem = *edge;    
            vector <long> value_tem = *value;   
            vector <Double_t> time_tem = *time;  
            vector <Int_t> sweeps_tem = *sweeps;    
            vector <long> sweeps_global_tem = *sweeps_global;  
            vector <Int_t> tag_tem = *tag;   
            vector <Int_t> lost_tem = *lost;
            vector <Int_t> glo_tem = *glo;

            long long Numhit=(long long)time_tem.size();

            double * time_cp = new double[Numhit];
            long long* tof_sequence_list = new long long[Numhit];

            for(long long index=0;index<Numhit;index++){
              time_cp[index] = time_tem[index];
            }

            TMath::Sort(Numhit,time_cp,tof_sequence_list,kFALSE); // tof from less to greater

            for(long long index=0;index<Numhit;index++){
                channel->at(index) = channel_tem[ tof_sequence_list[index] ];
                edge->at(index) = edge_tem[ tof_sequence_list[index] ];
                value->at(index) = value_tem[ tof_sequence_list[index] ];
                time->at(index) = time_tem[ tof_sequence_list[index] ];
                sweeps->at(index) = sweeps_tem[ tof_sequence_list[index] ];    
                sweeps_global->at(index) = sweeps_global_tem[ tof_sequence_list[index] ];
                tag->at(index) = tag_tem[ tof_sequence_list[index] ];
                lost->at(index) = lost_tem[ tof_sequence_list[index] ];
                glo->at(index) = glo_tem[ tof_sequence_list[index] ];
            }

            delete[] time_cp;
            delete[] tof_sequence_list;

      }


      void AddData(Int_t _channel, Int_t _edge, long _value, Double_t _time, Int_t _sweeps, long _sweeps_global, Int_t _tag, Int_t _lost, Int_t _glo){
                if(_value==0) _tag = !_tag;  //reverse the tag value of the trigger of sweep to have a correct tag value
                channel->push_back(_channel);
                edge->push_back(_edge);
                value->push_back(_value);
                time->push_back(_time);
                sweeps->push_back(_sweeps);    
                sweeps_global->push_back(_sweeps_global);
                tag->push_back(_tag);
                lost->push_back(_lost);
                glo->push_back(_glo);
      }

      void AddData_vector(vector <Int_t>* _channel, vector <Int_t>* _edge,vector <long>* _value, vector <Double_t>* _time, vector <Int_t>* _sweeps, vector <long>* _sweeps_global,
                                    vector <Int_t>* _tag, vector <Int_t>* _lost, vector <Int_t>* _glo){
                unsigned long Numhit = _channel->size();
                if(Numhit==0) return;

                for(unsigned long index=0;index<Numhit;index++){
                  AddData(_channel->at(index), _edge->at(index), _value->at(index), _time->at(index), _sweeps->at(index), 
                                  _sweeps_global->at(index), _tag->at(index), _lost->at(index), _glo->at(index));
                }

      }

      void Print(){
                if(GetNumEntry()==0) return;
                printf("time \t channel \t sweeps \t sweeps_global \t tag\n");
                for(unsigned long index=0; index<num_entry;index++){
                  printf("%.4f \t %d \t %d \t %ld \t %d\n",time->at(index),channel->at(index),sweeps->at(index),sweeps_global->at(index),tag->at(index));
                }
      }


      void FillTree(TTree* _tree=nullptr){
            if(_tree == nullptr){
              cout<<"\e[1;33m"<<"Error!!!!!! tree for data storing is not exist!!! Stop FillTree!"<<"\e[0m"<<endl;
              return;
            }

              _tree->SetBranchAddress("nevt", &evt_index); 
              _tree->SetBranchAddress("nhits", nhits);
              _tree->SetBranchAddress("channel", &channel);
              _tree->SetBranchAddress("edge", &edge);
              _tree->SetBranchAddress("value", &value);
              _tree->SetBranchAddress("time", &time);
              _tree->SetBranchAddress("sweeps", &sweeps);
              _tree->SetBranchAddress("sweeps_global", &sweeps_global);
              _tree->SetBranchAddress("glo", &glo);
              _tree->SetBranchAddress("tag", &tag);
              _tree->SetBranchAddress("lost", &lost);
              _tree->Fill();
              _tree->ResetBranchAddresses();
      }



      void DeleteContainer(){
                if(channel!=nullptr){ delete channel; channel = nullptr;}
                if(edge!=nullptr){ delete edge; edge= nullptr;}
                if(value!=nullptr){ delete value; value = nullptr;}
                if(time!=nullptr){ delete time; time = nullptr;}
                if(sweeps!=nullptr){ delete sweeps; sweeps = nullptr;}
                if(sweeps_global!=nullptr){ delete sweeps_global; sweeps_global = nullptr;}
                if(tag!=nullptr){ delete tag; tag = nullptr;}
                if(lost!=nullptr){ delete lost; lost = nullptr;}
                if(glo!=nullptr){ delete glo; glo = nullptr;}
      }

      ~EventContainer(){
                DeleteContainer();
               
      }

};


//************* Can combine different files *******************
// Set AddMode = true to use record; set AddMode = false to Reset
// Set IsTheLast = false to get ready for combining with next file
// first file-->  AddMode =false, IsTheLast=false
// middle files--> AddMode = true,IsTheLast=false
// the last file --> AddMode = true, IsTheLast=true;
// Decode a single file -->  AddMode =false, IsTheLast=true;
//***********************************************************
unsigned long int ParserMCS6A_3(string PATH = "../Run1/",string filename = "mcs_39Kvs143X2plus@502_160142",bool SelectChannel=true,int ChannelSelected=1,bool AddMode=false,bool IsTheLast=true){
  const int kMaxBufLen = 0x800000; // 8MB , get line of char in file
  static const int kstart = 6;     // the maximum of channel number, used in a "if()" condition
  unsigned long long int kChannel=0;    // = 0x00000007;
  int kShiftChannel         = 0;
  unsigned long long int kEdge=0;        //= 0x00000008;
  int kShiftEdge            =0;              //3;
  unsigned long long int kTimeData=0;    //= 0xfffffffff0;  //0xfffffff0;
  int kShiftData            = 0;              //4;
  unsigned long long int kSweeps=0;      //= 0x7f0000000000;  //0xffff00000000;
  int kShiftSweeps          = 0;            //40; //32;
  unsigned long long int kTagBit=0;      //= 0xffff000000000000;  //0x7fff000000000000;
  int kShiftTagBit          = 0;             //48;
  unsigned long long int kLost=0;        //= 0x800000000000; //0x8000000000000000;
  int kShiftLost            = 0;             //47;  //63;

  unsigned long long int sweepsclearbit=0;
  unsigned long long int sweepsfull=0;   // value = value of sweeps full + 1;  => 0b 1111 + 1 = 0b 10000;

  EventContainer * eve_container[2];
  eve_container[0] = nullptr;
  eve_container[1] = nullptr;

  bool EveContinerStat[2]={false,false}; // false --> no data in event container array

  Long64_t sweeps_delta_s=0;
  Long64_t sweeps_delta_l=0;


  string inputfile = PATH + "LST/" + filename + ".lst";//"LST/"
  ifstream fin;
  fin.open(inputfile.data(),ios::in);
  if(!fin){
    cerr<<"cannot open LST input file: "<<inputfile<<" from decorder !!! "<<endl;
    return 0;
  }
  

  char *buffer = new char[kMaxBufLen];

  int Nbitshift=0;


  while(!fin.eof()){
      fin.getline(buffer,kMaxBufLen);
      if(fin.eof()){
        cerr<<"cannot find data!!!"<<endl;
	  fin.close();
        return 0;
      }

      //&&&&&&&&&&&&&&&7 read data format description  &&&&&&&&&&&&&&&&&&&&&&
      if(strncmp(buffer,";bit",4)==0){
            unsigned int bitL=0, bitH=0;

            string sbitL;
            sbitL=buffer[4]; 
            if(isdigit(buffer[5])) sbitL = sbitL + buffer[5];
            bitL = atoi(sbitL.c_str());

            string _buffer(buffer);

            string::size_type position =  _buffer.find("..");

            if(position<10){ // using more than one bit
                string sbitH;
                for(int i=(int)position +2;i<11;i++){
                    if(isdigit(buffer[i])) sbitH +=buffer[i];
                    if(buffer[i]==':') break;
                }

                bitH = atoi(sbitH.c_str());
                if(bitH==0){cout<<"error at reading bits format for decorder, break!!"<<endl; fin.close(); return 0;}


            }
            else{ bitH = bitL;} // using only one bit


            unsigned long long int GetBitRange=0;
            unsigned long long int one=1;

            for(unsigned int i=bitL;i<=bitH;i++){GetBitRange+= (one<<i);}  // generate bits to use & for saving corresponding bits



            if(_buffer.find("channel") != _buffer.npos){
                //printf("channel:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);
                kChannel = GetBitRange;
                kShiftChannel = Nbitshift;
            }
            else if(_buffer.find("edge") != _buffer.npos){
                //printf("kEdge:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);  
                kEdge = GetBitRange;
                kShiftEdge = Nbitshift;         
            }
            else if(_buffer.find("timedata") != _buffer.npos){
                 //printf("timedata:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);
                 kTimeData = GetBitRange;
                 kShiftData = Nbitshift;        
            }
            else if(_buffer.find("sweeps") != _buffer.npos){
                
                for(unsigned int i=0;i<(bitH-bitL)+1;i++){sweepsclearbit+= (one<<i);}

                sweepsfull = sweepsclearbit + 1;   // value = value of sweeps full + 1;  => 0b 1111 + 1 = 0b 10000;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		 //   sweepsMax = sweepsfull;  // only use in Online.C ==> to know the maximum number allowed sweeps number 
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sweepsclearbit=~sweepsclearbit;
                //printf("sweeps:%d .. %d; %llx ;kshift=%d ;sweepclear=%llx ;sweepsfull value= %llx \n ",bitL,bitH,GetBitRange,Nbitshift,sweepsclearbit,sweepsfull);  
                kSweeps= GetBitRange;
                kShiftSweeps = Nbitshift;

            }
            else if(_buffer.find("tag") != _buffer.npos){
                 //printf("tag:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift);    
                 kTagBit = GetBitRange;
                 kShiftTagBit = Nbitshift; 
            }
            else if(_buffer.find("data_lost") != _buffer.npos){
                //printf("data_lost:%d .. %d; %llx ;kshift=%d \n ",bitL,bitH,GetBitRange,Nbitshift); 
                kLost= GetBitRange;
                kShiftLost = Nbitshift;         
            }
            else{cout<<"Error; find unknow property of bits. Break!!!!"<<endl; fin.close(); return 0;}

            Nbitshift+= bitH-bitL +1;

      }


      if(strncmp(buffer,"[DATA]",6)==0){  // strncmp: compare buffer with "[DATA]" for 6 character, if equal will be 0; buffer> data -> value>0
                                                      // buffer< data  -> value<0; based on ACSII sequence
        cout<<"\e[1;33m reading data from "<<inputfile<<"\e[0m"<<endl;
        break;
      }
  }//end of file head reading for data format

  
  Long64_t nevt = 0;
 static Long64_t evtIndex=0;
  int nhits[10];

  int isweep_old=-1;
  int isweep_old_cp = isweep_old;
  int iglo =0;
  long isweep_global =0;
  long isweep_global_old=-1;
  long isweep_global_max=0;

  static long isweeps_global_offset=0;
  static int isweep_offset=0;

  if(!AddMode){
    isweeps_global_offset=0;
    isweep_offset=0;
    evtIndex=0;
  }


  vector <Int_t> *channel = new vector <Int_t>();		vector< vector <Int_t> > eve_channel;	
  vector <Int_t> *edge = new vector <Int_t>();		vector< vector <Int_t> > eve_edge;
  vector <long> *value = new vector <long>();		vector< vector <long> > eve_value;
  vector <Double_t> *time = new vector <Double_t>();	vector< vector <Double_t> > eve_time;
  vector <Int_t> *sweeps = new vector <Int_t>();		vector< vector <Int_t> > eve_sweeps;
  vector <long> * sweeps_global = new vector <long>;	vector< vector <long> > eve_sweeps_global;
  vector <Int_t> *tag = new vector <Int_t>();		vector< vector <Int_t> > eve_tag;
  vector <Int_t> *lost = new vector <Int_t>();		vector< vector <Int_t> > eve_lost;
vector <Int_t> *glo = new vector <Int_t>();		vector< vector <Int_t> > eve_glo;

  string outputfile = PATH + "rootfiles/" + filename + ".root";
  static TFile *fout =NULL;
  if(!AddMode)fout = new TFile(outputfile.c_str(),"RECREATE");
  if(!fout->IsOpen()) {cout<<"Can not Recreate file: "<<outputfile<<" from decorder."<<endl; cout<<"Break!!"<<endl; return 0;}

  static TTree *tree=NULL;

  if(!AddMode){
      tree = new TTree("tree","rawdata tree");
      tree->Branch("nevt",&nevt,"nevt/L");
      tree->Branch("nhits",nhits,"nhits[10]/I");
      tree->Branch("channel",&channel);
      tree->Branch("edge",&edge);
      tree->Branch("value",&value); // in unit of 100ps
      tree->Branch("time",&time); // in unit of ns
      tree->Branch("sweeps",&sweeps);
      tree->Branch("sweeps_global",&sweeps_global);
    tree->Branch("glo",&glo);
      tree->Branch("tag",&tag);
      tree->Branch("lost",&lost);
  }

  unsigned long long int evtdata;
  string datastring;
  char datachar[100];
  istringstream is;
  while(!fin.eof()){

    if(nevt%1000==0) cout<<'\r'<<"read "<<nevt<<" events"<<flush;
    fin>>datachar; 
    datastring = string(datachar);
    is.clear();
    is.str(datastring);
    if(datastring.size()!=16){
	if(nevt==0 && strncmp(datachar,"[DATA]",6)==0) return 0; // in case there is no data after [DATA]
	else{	cout<<"read unknown data"<<endl; continue;}
     }
    is>>hex>>evtdata; 
   // fin >> hex >> evtdata;     // value in evtdata will be ox.... but anyway cout<< evtdata is in decimal; data is always binary in computer
    int ich = (evtdata&kChannel) >> kShiftChannel;  // extract 0...2 bits: channel:1~6 channel
    int iedge = (evtdata&kEdge) >> kShiftEdge;   // extract 4th bit; right shift 3 bit
    long ivalue = (evtdata&kTimeData) >> kShiftData;
    double itime = ivalue*0.1;   // resolution is 0.1 nanosecond
    int isweeps = (evtdata&kSweeps) >> kShiftSweeps;
    int itag = (evtdata&kTagBit) >> kShiftTagBit;
    int ilost = (evtdata&kLost) >> kShiftLost;

    // using tag0
    itag = itag&0x1;   // when itag==1 -> itag=1; when itag==0; -> itag =0; extract the first bit of tag value
    
    if(evtdata==0){
      continue;
    }else if(ich<1 || ich>kstart ){
      cerr<<"read unknow data: "<<evtdata<<endl;
      continue;
    }else if(SelectChannel && ich!=ChannelSelected && ich!=6){ // ch6 => sweep start trigger
	continue;
    }else if(isweeps != isweep_old){       //(ich == kstart && ivalue == 0){
      // new event
/***********************************************************************
        if(isweeps > isweep_old || isweeps>(isweep_old-10)){  
//cout<<"global sweep = "<<isweep_global<<" , isweeps "<<isweeps<<endl;
          //  isweep_global = (isweep_global&0xFFFFFFFFFFFF0000) + isweeps; // clear lower 4 digit   // 0xFFFFFFFFFFFF0000
          //      iglo = (iglo&0xFFFF0000) + isweeps;

          //  isweep_global = (isweep_global&0xFFFFFFFFFFFFFF80) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
          //      iglo = (iglo&0xFFFFFF80) + isweeps;

            isweep_global = (isweep_global&sweepsclearbit) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
                iglo = (iglo&sweepsclearbit) + isweeps;

        }else{

             // isweep_global = (isweep_global&0xFFFFFFFFFFFF0000) + isweeps + 0x10000;
//iglo = (iglo&0xFFFF0000) + isweeps + 0x10000;

              //isweep_global = (isweep_global&0xFFFFFFFFFFFFFF80) + isweeps + 0x80;
//iglo = (iglo&0xFFFFFF80) + isweeps + 0x80;

              isweep_global = (isweep_global&sweepsclearbit) + isweeps + (int)sweepsfull;
iglo = (iglo&sweepsclearbit) + isweeps + (int)sweepsfull;
          }
*********************************************************************************/

	if(isweeps > isweep_old){
		sweeps_delta_l = isweeps - isweep_old;   // if same around as isweep_old
		sweeps_delta_s = sweepsfull - isweeps + isweep_old; //if one around before isweep_old; round -1
		if(sweeps_delta_l > sweeps_delta_s){//last round;
			isweep_global = (isweep_global&sweepsclearbit) + isweeps - sweepsfull;
			iglo = (iglo&sweepsclearbit) + isweeps - sweepsfull;
		}
		else{ //same around as isweep_old;
			    isweep_global = (isweep_global&sweepsclearbit) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
                iglo = (iglo&sweepsclearbit) + isweeps;
		}
	}
	else{
		sweeps_delta_s = isweep_old - isweeps;  // if same around as isweep_old
		sweeps_delta_l = sweepsfull - isweep_old + isweeps;  //if one around after isweep_old; round+1
		if(sweeps_delta_s>sweeps_delta_l){ //next round
				isweep_global = (isweep_global&sweepsclearbit) + isweeps + sweepsfull;
				iglo = (iglo&sweepsclearbit) + isweeps + sweepsfull;
		}
		else{//same around as isweep_old
				isweep_global = (isweep_global&sweepsclearbit) + isweeps; // clear lower 7 digit (binary)  // 0xFFFFFFFFFFFF0000
		        iglo = (iglo&sweepsclearbit) + isweeps;
		}
	}


      isweep_old = isweeps;  //update isweep_old to new one,  new event: nevt++
//cout<<"global sweep = "<<isweep_global<<" , isweeps "<<isweeps<<endl;
//return ;
   if(isweep_global_max<isweep_global) isweep_global_max=isweep_global;

  if(nevt>0){ 
       //   tree->Fill(); // do not fill empty data at the beginning

/*//%%%%%%%%%%%%% modified for making sequence
		eve_channel.push_back(*channel);	
		eve_edge.push_back(*edge);
		eve_value.push_back(*value);
		eve_time.push_back(*time);
		eve_sweeps.push_back(*sweeps);
		eve_sweeps_global.push_back(*sweeps_global);
		eve_tag.push_back(*tag);
		eve_lost.push_back(*lost);
		eve_glo.push_back(*glo);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

          int SectionIndex = (isweep_global_old/sweepsfull)%2;
          if(eve_container[SectionIndex]==nullptr){ eve_container[SectionIndex] = new EventContainer[sweepsfull];  EveContinerStat[SectionIndex]=false;}
          eve_container[SectionIndex][isweep_old_cp].AddData_vector(channel, edge, value, time, sweeps, sweeps_global, tag, lost, glo);
          EveContinerStat[SectionIndex]=true;

          int SectionIndex_old = (SectionIndex+1)%2;
          if((unsigned int)isweep_old_cp > sweepsfull/2 && (unsigned long)isweep_global_old>sweepsfull && EveContinerStat[SectionIndex_old]){//fill data in last container to tree,  > sweepsfull/2 guarantee no more update need to that array
              
              for(unsigned long long int index=0; index<sweepsfull;index++){
                    if( eve_container[SectionIndex_old][index].GetNumEntry() > 0 ){
                        evtIndex++;
                        eve_container[SectionIndex_old][index].evt_index = evtIndex;
                        eve_container[SectionIndex_old][index].Sort();
                        eve_container[SectionIndex_old][index].FillTree(tree);
                    }
              }

              delete[] eve_container[SectionIndex_old];
              eve_container[SectionIndex_old] = nullptr;
              EveContinerStat[SectionIndex_old]=false;
          }

	}

      // clean container
      channel->clear();
      edge->clear();
      value->clear();
      time->clear();
      sweeps->clear();
      sweeps_global->clear();
glo->clear();
      tag->clear();
      lost->clear();
      for(int i=0; i<10; i++) nhits[i] = 0;
      // load data
      nhits[ich]++;
      channel->push_back(ich);
      edge->push_back(iedge);
      value->push_back(ivalue);
      time->push_back(itime);

      sweeps->push_back( (isweeps + isweep_offset)%((int)sweepsfull) );
      sweeps_global->push_back(isweep_global + isweeps_global_offset);
glo->push_back( (int)(iglo + isweeps_global_offset) );
      tag->push_back(itag);
      lost->push_back(ilost);
      nevt++;
      isweep_old_cp = isweeps;
      isweep_global_old = isweep_global;

    }else{
      if(!fin.eof()){  // prevent to fill repeated data  // multi-hits
      // real time data
      nhits[ich]++;
      channel->push_back(ich);
      edge->push_back(iedge);
      value->push_back(ivalue);
      time->push_back(itime);
      sweeps->push_back( (isweeps + isweep_offset)%((int)sweepsfull) );
      sweeps_global->push_back( isweep_global + isweeps_global_offset );
glo->push_back( (int)(iglo + isweeps_global_offset) );
      tag->push_back(itag);
      lost->push_back(ilost);
       }
    }
  }
  cout<<'\r'<<flush<<"read "<<nevt<<" events"<<endl;
  cout<<"\e[1;33m writing tree to "<<outputfile<<"\e[0m"<<endl; // show color
  fout->cd();

  //tree->Fill(); // fill last data

  //fill last piece of data from vector to container
    int SectionIndex = (isweep_global_old/sweepsfull)%2;
    if(eve_container[SectionIndex]==nullptr){ eve_container[SectionIndex] = new EventContainer[sweepsfull];  EveContinerStat[SectionIndex]=false;}
    eve_container[SectionIndex][isweep_old_cp].AddData_vector(channel, edge, value, time, sweeps, sweeps_global, tag, lost, glo);
    EveContinerStat[SectionIndex]=true;

    int SectionIndex_old = (SectionIndex+1)%2;
    if(EveContinerStat[SectionIndex_old]){//fill data in last container to tree
        
        for(unsigned long long int index=0; index<sweepsfull;index++){
              if( eve_container[SectionIndex_old][index].GetNumEntry() > 0 ){
                  evtIndex++;
                  eve_container[SectionIndex_old][index].evt_index = evtIndex;
                  eve_container[SectionIndex_old][index].Sort();
                  eve_container[SectionIndex_old][index].FillTree(tree);
              }
        }

        delete[] eve_container[SectionIndex_old];
        eve_container[SectionIndex_old] = nullptr;
        EveContinerStat[SectionIndex_old]=false;
    }

    unsigned long long int index_tmp=0;

    for(unsigned long long int index=0; index<sweepsfull;index++){
      if( eve_container[SectionIndex][index].GetNumEntry() > 0 ){
          evtIndex++;
          eve_container[SectionIndex][index].evt_index = evtIndex;
          eve_container[SectionIndex][index].Sort();
          eve_container[SectionIndex][index].FillTree(tree);
          index_tmp = index;
      }
    }

long last_isweep_global = (eve_container[SectionIndex][index_tmp].sweeps_global)->at(0);
int last_sweeps = (eve_container[SectionIndex][index_tmp].sweeps)->at(0);

cout<<"last sweep global = "<<last_isweep_global<<endl;
cout<<"last sweep ="<<last_sweeps<<endl;


    delete[] eve_container[SectionIndex];
    eve_container[SectionIndex] = nullptr;
    EveContinerStat[SectionIndex]=false;


/*
	//%%%%%%%modified for sequence
		//save the last record of data
		eve_channel.push_back(*channel);	
		eve_edge.push_back(*edge);
		eve_value.push_back(*value);
		eve_time.push_back(*time);
		eve_sweeps.push_back(*sweeps);
		eve_sweeps_global.push_back(*sweeps_global);
		eve_tag.push_back(*tag);
		eve_lost.push_back(*lost);
		eve_glo.push_back(*glo);

		const unsigned long N_data_eve = eve_sweeps_global.size();
		// combine data with same global sweeps
		for(unsigned long event_i=0;event_i<N_data_eve;event_i++){
			long glo_sweeps_i = eve_sweeps_global[event_i][0];
			for(unsigned long event_i_before=0; event_i_before< event_i; event_i_before++){
				if(glo_sweeps_i == eve_sweeps_global[event_i_before][0] ){ // find data with same global sweeps ==> combine
					for(unsigned long ihit=0;ihit<(unsigned long)eve_sweeps_global[event_i].size();ihit++){

						eve_channel[event_i_before].push_back( eve_channel[event_i][ihit] );	
						eve_edge[event_i_before].push_back( eve_edge[event_i][ihit] );
						eve_value[event_i_before].push_back( eve_value[event_i][ihit] );
						eve_time[event_i_before].push_back( eve_time[event_i][ihit] );
						eve_sweeps[event_i_before].push_back( eve_sweeps[event_i][ihit] );
						eve_sweeps_global[event_i_before].push_back( eve_sweeps_global[event_i][ihit] );
						eve_tag[event_i_before].push_back( eve_tag[event_i][ihit] );
						eve_lost[event_i_before].push_back( eve_lost[event_i][ihit] );
						eve_glo[event_i_before].push_back( eve_glo[event_i][ihit] );


					}
					eve_sweeps_global[event_i][0]=-99; // has been combine (marker)
					break; // next current event;
				}
			}

		}

		// make sequence for tof in each global sweep
			// first: sort the tof in a global sweep

		long long tem_gsweeps[N_data_eve];
		long long gsweeps_sequence_list[N_data_eve];

		for(unsigned long event_i=0;event_i<N_data_eve;event_i++){
			if(eve_sweeps_global[event_i][0]>0){// first global sweep ==1; reject glo sweeps ==-99
				if(eve_sweeps_global[event_i].size() >1 ){ // more than one hit
					const unsigned long hits_num = eve_time[event_i].size();
					double time_tem[hits_num];
					long long tof_sequence_list[hits_num];
					int channel_seq[hits_num];
					int edge_seq[hits_num];		
					long value_seq[hits_num];		
					int tag_seq[hits_num]; 		
					int lost_seq[hits_num]; 		
					int glo_seq[hits_num]; 						

					for(unsigned long ihit=0;ihit<(unsigned long)eve_time[event_i].size();ihit++){
						time_tem[ihit]=eve_time[event_i][ihit];
					}
					sort(eve_time[event_i].begin(),eve_time[event_i].end());

					TMath::Sort((long long)hits_num,time_tem,tof_sequence_list,kFALSE); // tof from less to greater

					for(unsigned long indexi=0;indexi<(unsigned long)eve_time[event_i].size();indexi++){
						long long get_ch_index =tof_sequence_list[indexi];
						 channel_seq[indexi] = eve_channel[event_i][get_ch_index];
						 edge_seq[indexi]=eve_edge[event_i][get_ch_index];		
						 value_seq[indexi]=eve_value[event_i][get_ch_index];		
						 tag_seq[indexi]=eve_tag[event_i][get_ch_index]; 		
						 lost_seq[indexi]=eve_lost[event_i][get_ch_index]; 		
						 glo_seq[indexi]=eve_glo[event_i][get_ch_index]; 
					}

					for(unsigned long indexi=0;indexi<(unsigned long)eve_time[event_i].size();indexi++){
						
						 eve_channel[event_i][indexi]=channel_seq[indexi];
						 eve_edge[event_i][indexi]=edge_seq[indexi];		
						 eve_value[event_i][indexi]=value_seq[indexi];		
						 eve_tag[event_i][indexi]=tag_seq[indexi]; 		
						 eve_lost[event_i][indexi]=lost_seq[indexi]; 		
						 eve_glo[event_i][indexi]=glo_seq[indexi]; 
					}
			
					
				}
			}

			tem_gsweeps[event_i]=eve_sweeps_global[event_i][0];
		}

			// second: sort the global sweeps
		TMath::Sort((long long)N_data_eve,tem_gsweeps,gsweeps_sequence_list,kFALSE); // sort by less => large

			// third: save data to tree
		long long getindex_content=0;
		nevt=0; //clear;
		for(unsigned long event_i=0;event_i<N_data_eve;event_i++){
			getindex_content = gsweeps_sequence_list[event_i];
			if(eve_sweeps_global[getindex_content][0]!=-99){// save to tree

			      channel->clear();
			      edge->clear();
			      value->clear();
			      time->clear();
			      sweeps->clear();
			      sweeps_global->clear();
			glo->clear();
			      tag->clear();
			      lost->clear();
			      for(int i=0; i<10; i++) nhits[i] = 0;


				for(unsigned long ihit=0;ihit<(unsigned long)eve_sweeps_global[getindex_content].size();ihit++){
						//load data of one sweep
						channel->push_back( eve_channel[getindex_content][ihit] );	
						edge->push_back( eve_edge[getindex_content][ihit] );
						value->push_back( eve_value[getindex_content][ihit] );
						time->push_back( eve_time[getindex_content][ihit] );
						sweeps->push_back( eve_sweeps[getindex_content][ihit] );
						sweeps_global->push_back( eve_sweeps_global[getindex_content][ihit] );
						tag->push_back( eve_tag[getindex_content][ihit] );
						lost->push_back( eve_lost[getindex_content][ihit] );
						glo->push_back( eve_glo[getindex_content][ihit] );


				}
				
				for(unsigned long ihit=0;ihit<(unsigned long)channel->size();ihit++){// hit statics of channel
					nhits[channel->at(ihit)]++;
				}

				tree->Fill();
				
				nevt++;
			}
		}


	//%%%%%%%%%%%%%%%% end of modified for sequence

*/


if(!IsTheLast){
        isweep_offset =  (last_sweeps + isweep_offset)%((int)sweepsfull);
        isweeps_global_offset = last_isweep_global;
        return 1 ;
    
}


  tree->Write();
  
  delete [] buffer;

   fout->Write();
  //fout->Close();
delete fout;
tree=NULL;
fout=NULL;
isweep_offset=0;
isweeps_global_offset=0;

cout<<"total events= "<<evtIndex<<endl;

evtIndex=0;

  fin.clear();
  fin.seekg(0,ios::end);
  unsigned long int FileEnd = fin.tellg();

  fin.close();
   is.clear();
cout<<"check global sweeps: "<<isweep_global<<" Maxsweep= "<<isweep_global_max<<endl;
  return FileEnd;
}





void DecodeAndCombine(string PATH = "../",string filename = "mcs_39Kvs143X2plus@502_160142",bool SelectChannel=true,int ChannelSelected=1){

    if( gSystem->AccessPathName(PATH.c_str()) ){ cout<<"Path is not exist! Abort!!"<<endl; return;}

    vector<string>filenamelist;

    string linux_command = "ls -tr " + PATH + "LST/" + filename +"*.lst";

    char getfilename[1000]={'\0'};

    FILE* pipe = popen(linux_command.c_str(), "r");

    if (!pipe) {
      std::cerr << "popen() failed!" << std::endl;
      return;
    }

    try{
      while(!feof(pipe))
      {fgets(getfilename,1000,pipe);
        if(!feof(pipe))filenamelist.push_back(getfilename);
        else break;
      }
    }
    catch(...){
      pclose(pipe);
      throw;
    }


    pclose(pipe);

    int FileNumber = (int)filenamelist.size();

    if(FileNumber==0){
      cout<<"\e[1;31m"<<filename+"*.lst"<<" is not exist!!! Please check again!!\e[0m"<<endl;
      return;
    }else if(FileNumber==1){
        ParserMCS6A_3(PATH , filename, SelectChannel, ChannelSelected);
    }else{
        string FilePrimary="";
        vector<string> FileScan_tmp;
        for(int i=0;i<FileNumber;i++){
            while(filenamelist[i].find("/") != string::npos){filenamelist[i].erase(filenamelist[i].begin(),filenamelist[i].begin()+filenamelist[i].find("/")+1); }
            if(filenamelist[i].find(".lst") != string::npos){filenamelist[i].erase(filenamelist[i].begin()+filenamelist[i].find(".lst"),filenamelist[i].end()); }
            if(filenamelist[i] == filename) FilePrimary =filename;
            else FileScan_tmp.push_back(filenamelist[i]);
        }

        if(FileScan_tmp.size()>1){// sort
            for(int i=0;i<(int)FileScan_tmp.size();i++){
                for(int j=i+1;j<(int)FileScan_tmp.size();j++){
                  if(FileScan_tmp[j]<FileScan_tmp[i]){
                    string stem = FileScan_tmp[i];
                    FileScan_tmp[i] = FileScan_tmp[j];
                    FileScan_tmp[j]= stem;
                  }
                }
            }
          }

        filenamelist.clear();
        if(FilePrimary!="") filenamelist.push_back(FilePrimary);

        for(int i=0;i<(int)FileScan_tmp.size();i++){
          filenamelist.push_back(FileScan_tmp[i]);
        }

        cout<<"\e[1;33m Get the list of files to decode:\e[0m"<<endl;
        for(int i=0;i<FileNumber;i++){
          cout<<filenamelist[i]<<".lst"<<endl;
        }

        for(int i=0;i<FileNumber;i++){
              if(i==0){// first file
                  if(  ParserMCS6A_3(PATH, filenamelist[i], SelectChannel, ChannelSelected, false, false) == 0) return;
              }else if(i==FileNumber-1){//the last file
                  if(  ParserMCS6A_3(PATH, filenamelist[i], SelectChannel, ChannelSelected, true, true) == 0) return;
              }
              else{// middle files
                  if(  ParserMCS6A_3(PATH, filenamelist[i], SelectChannel, ChannelSelected, true, false) == 0) return;
              }
        }

    }

}

/*
#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc>1){
    ParserMCS6A_2( string(argv[1]) );
  }
  return 0;
}
#endif
*/



