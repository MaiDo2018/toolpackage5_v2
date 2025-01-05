// Mainframe macro generated from application: /home/xian/rootbuild/root_v6_22_06/bin/root.exe
// By ROOT version 6.22/06 on 2020-12-31 09:54:37

#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TGMdiDecorFrame
#include "TGMdiDecorFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGMdiFrame
#include "TGMdiFrame.h"
#endif
#ifndef ROOT_TGMdiMainFrame
#include "TGMdiMainFrame.h"
#endif
#ifndef ROOT_TGMdiMenu
#include "TGMdiMenu.h"
#endif
#ifndef ROOT_TGColorDialog
#include "TGColorDialog.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGNumberEntry
#include "TGNumberEntry.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGuiBldHintsEditor
#include "TGuiBldHintsEditor.h"
#endif
#ifndef ROOT_TRootBrowser
#include "TRootBrowser.h"
#endif
#ifndef ROOT_TGuiBldNameFrame
#include "TGuiBldNameFrame.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGShutter
#include "TGShutter.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCommandPlugin
#include "TGCommandPlugin.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGuiBldEditor
#include "TGuiBldEditor.h"
#endif
#ifndef ROOT_TGColorSelect
#include "TGColorSelect.h"
#endif
#ifndef ROOT_TGTextEdit
#include "TGTextEdit.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TRootContextMenu
#include "TRootContextMenu.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGView
#include "TGView.h"
#endif
#ifndef ROOT_TGProgressBar
#include "TGProgressBar.h"
#endif
#ifndef ROOT_TGMsgBox
#include "TGMsgBox.h"
#endif
#ifndef ROOT_TRootGuiBuilder
#include "TRootGuiBuilder.h"
#endif
#ifndef ROOT_TGFileBrowser
#include "TGFileBrowser.h"
#endif
#ifndef ROOT_TGTab
#include "TGTab.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGTextEditor
#include "TGTextEditor.h"
#endif
#ifndef ROOT_TRootCanvas
#include "TRootCanvas.h"
#endif
#ifndef ROOT_TGStatusBar
#include "TGStatusBar.h"
#endif
#ifndef ROOT_TGListTree
#include "TGListTree.h"
#endif
#ifndef ROOT_TGuiBldGeometryFrame
#include "TGuiBldGeometryFrame.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif
#ifndef ROOT_TGToolBar
#include "TGToolBar.h"
#endif
#ifndef ROOT_TRootEmbeddedCanvas
#include "TRootEmbeddedCanvas.h"
#endif
#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif
#ifndef ROOT_TGuiBldDragManager
#include "TGuiBldDragManager.h"
#endif
#ifndef ROOT_TGHtmlBrowser
#include "TGHtmlBrowser.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGTab
#include "TGTab.h"
#endif


#include "Riostream.h"


#include "TROOT.h"
#include <vector>
#include <thread>
#include "DriftCorrect3.C"
#include "TCanvas.h"
#include "TGClient.h"
#include <fstream>
#include "TVirtualX.h"
#include <stack>
#include <TTimer.h>
#include "TLine.h"

using namespace std;


class MassAnaGUI : public TGMainFrame {

private:
   TGMainFrame *fMainFrame2609;
   TGCompositeFrame *fMainFrame1933, *fMainFrame1817;
   TGHorizontalFrame *MassAnalysis_v1;
   TGFont *ufont;
   TGGC   *uGC;
   TGTextEntry *fTextEntryPATH ,*fTextEntryFilename;
   TGLabel *fLabelPATH , *fLabelFilename, *fLabelHotkey;
   TGTextButton     *fTextButtonDecoder, *fTextButtonLoadTree, *fTextButtonExit;
   TGRadioButton *fTextButtonOriginal , *fTextButtonCorrected;
   TGCheckButton *fTextButtonCHSelect, *fTextButtonCombineLST, *fTextButtonFixPars;
   TGComboBox *fComboBoxCHlist, *fComboBoxTaglist, *fComboBoxCHlist_drif, *fComboBoxFitFunc;
   
   //TGHProgressBar *fHProgressBar1546;

   // drift correct part
   TGGroupFrame *fGroupFrameDC;
   TGNumberEntry *fNumberEntryRefCento;
   TGNumberEntry *fNumberEntryNbins;
   TGNumberEntry *fNumberEntryEveAkill;
   TGNumberEntry *fNumberEntryRefSigma;
   TGNumberEntry *fNumberEntryHisto_halfwidth;
   TGNumberEntry *fNumberEntryTveto;
   
   TGLabel *fLabelRefCenter;  // share
   TGLabel *fLabelnbins;		
   TGLabel *fLabelEveAkill;
   TGLabel *fLabelRefSigma;
   TGLabel *fLabelHalf_histoWidth; // share
   TGTextButton *fTextButtonOK;
   TGLabel *fLabelTveto;
   TGLabel  *fLabelSeleTag;
   TGLabel *fLabelSeleCh;
   TGLabel *fLabelFitFun;

TGTab *fTab662;
TRootEmbeddedCanvas *fRootEmbeddedCanvas676, *fRootEmbeddedCanvas2;
TCanvas *cref; // for 2D histo
//TCanvas *cref2; // for 1D histo
TH1D* href1D;
TH2D* href2D;
TLine* lcursor;


     //T0 entry;
   TGLabel *fLabelT0;
   TGNumberEntry *fNumberEntryT0;
//thread* t;
   TCanvas* cc;
   bool runcorrect;
   void runthread(string _PAHT,string _filename,double _centro, int _nbins, int _event_akill, double _ref_sigma,double _Half_hiswidth,double __intTo,TGProgressBar* _bar,
      int _SeleTag, int _SeleCH, double _Timeveto, double * _outpars,vector<double> *_timesANDsweep);
   void exc(string _PAHT="",string _filename="",double _centro=0, int _nbins=0, int _event_akill=0, double _ref_sigma=0,double _Half_hiswidth=0,double __intTo=0,TGProgressBar* _bar=NULL);

	//decode channel selection
	bool SelectChannel;
	int CHselected;

	// saving history
	char buffer_Path[1000];
	char buffer_filename[1000];
	char buffer_num[1000];
	double tof_history;
	double t0_history;
	bool kSuccess;
	void ReadLog();
	void WriteLog();
   void SetRange(int index,TCanvas* c_get,TH1D* h,double xmin, double xmax);
   void SetLinePos(double x1,double y1,double x2,double y2);


public:
   MassAnaGUI();
   virtual ~MassAnaGUI();
   // slots
	void ActionSetDriftON();
	void ActionSetDriftOFF();
	void ActivateDecoder();
	void ActivateDrifCorrect();
	void CloseWindow();
void cursor();
void subthread();
void loop();
   TTimer* timer;
   TGHProgressBar *fHProgressBar1546;
   TCanvas *cref2;
   thread* t2;

   vector<double> tof_and_enteries[2]; // [0]==> tof; [1] ==> entery index
   double get_drift_corr_pars[20];
   /*
   0=> tree enteries
   1 => num_of_SlideWidth; //bin number of sweeps axis
   2 => nbins;  // bins of tof histogram
   3 => binlowedge1; // low limit of tof histogram
   4 => binupedge1; // upper limit of tof histogram
   */

   ClassDef(MassAnaGUI, 0)
};


MassAnaGUI :: MassAnaGUI():cc(NULL),runcorrect(false),tof_history(-1),t0_history(-1),kSuccess(false),SelectChannel(true),CHselected(1),href1D(NULL),href2D(NULL),lcursor(NULL)
{

   // Read history log
   ReadLog();

   // main frame
	
   fMainFrame2609 = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
   fMainFrame2609->SetCleanup(kDeepCleanup);
   fMainFrame2609->SetName("fMainFrame2609");

  // composite frame
   fMainFrame1933 = new TGCompositeFrame(fMainFrame2609,781,452,kVerticalFrame);
   fMainFrame1933->SetName("fMainFrame1933");
   fMainFrame1933->SetLayoutBroken(kTRUE);


// composite frame
   fMainFrame1817 = new TGCompositeFrame(fMainFrame1933,876,539,kVerticalFrame);
   fMainFrame1817->SetName("fMainFrame1817");

   // horizontal frame
   MassAnalysis_v1 = new TGHorizontalFrame(fMainFrame1817,874,737,kVerticalFrame);
   MassAnalysis_v1->SetName("MassAnalysis_v1");
   MassAnalysis_v1->SetLayoutBroken(kTRUE);

   //TGFont *ufont;         // will reflect user font changes
   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   //TGGC   *uGC;           // will reflect user GC changes
   // PATH input bar
   GCValues_t valEntryPATH;
   valEntryPATH.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",valEntryPATH.fForeground);
   gClient->GetColorByName("#e8e8e8",valEntryPATH.fBackground);
   valEntryPATH.fFillStyle = kFillSolid;
   valEntryPATH.fFont = ufont->GetFontHandle();
   valEntryPATH.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valEntryPATH, kTRUE);
   fTextEntryPATH = new TGTextEntry(MassAnalysis_v1, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
   fTextEntryPATH->SetMaxLength(4096);
   fTextEntryPATH->SetAlignment(kTextLeft);
   if(kSuccess)fTextEntryPATH->SetText(buffer_Path);
   else fTextEntryPATH->SetText("../");
   fTextEntryPATH->Resize(296,fTextEntryPATH->GetDefaultHeight());
   MassAnalysis_v1->AddFrame(fTextEntryPATH, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextEntryPATH->MoveResize(120,8,296,20);


   // Text label of "PATH" & "filename"
   fLabelPATH = new TGLabel(MassAnalysis_v1,"PATH");
   fLabelPATH->SetTextJustify(36);
   fLabelPATH->SetMargins(0,0,0,0);
   fLabelPATH->SetWrapLength(-1);
   MassAnalysis_v1->AddFrame(fLabelPATH, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelPATH->MoveResize(0,0,112,32);
   fLabelFilename = new TGLabel(MassAnalysis_v1,"File Name");
   fLabelFilename->SetTextJustify(36);
   fLabelFilename->SetMargins(0,0,0,0);
   fLabelFilename->SetWrapLength(-1);
   MassAnalysis_v1->AddFrame(fLabelFilename, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelFilename->MoveResize(8,32,96,32);


   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   // Filename input bar
   GCValues_t valEntryFilename;
   valEntryFilename.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",valEntryFilename.fForeground);
   gClient->GetColorByName("#e8e8e8",valEntryFilename.fBackground);
   valEntryFilename.fFillStyle = kFillSolid;
   valEntryFilename.fFont = ufont->GetFontHandle();
   valEntryFilename.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valEntryFilename, kTRUE);
   fTextEntryFilename = new TGTextEntry(MassAnalysis_v1, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
   fTextEntryFilename->SetMaxLength(4096);
   fTextEntryFilename->SetAlignment(kTextLeft);

   if(kSuccess) fTextEntryFilename->SetText(buffer_filename);
   else fTextEntryFilename->SetText("Input filename without .lst or .root");

   fTextEntryFilename->Resize(296,fTextEntryFilename->GetDefaultHeight());
   MassAnalysis_v1->AddFrame(fTextEntryFilename, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextEntryFilename->MoveResize(120,32,296,20);

//CheckBox for channel

   fComboBoxCHlist = new TGComboBox(MassAnalysis_v1,"CHlist",-1,kHorizontalFrame | kSunkenFrame | kOwnBackground);
   fComboBoxCHlist->SetName("CHlist");
   fComboBoxCHlist->AddEntry("CH 1 ",0);
   fComboBoxCHlist->AddEntry("CH 2",1);
   fComboBoxCHlist->AddEntry("CH 3 ",2);
   fComboBoxCHlist->AddEntry("CH 4 ",3);
   fComboBoxCHlist->AddEntry("CH 5 ",4);
   fComboBoxCHlist->AddEntry("CH 6 ",5);
   //fComboBoxCHlist->AddEntry("CH 7 ",6);
 //  fComboBoxCHlist->Resize(80,20);
   fComboBoxCHlist->Select(0);
   MassAnalysis_v1->AddFrame(fComboBoxCHlist, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
fComboBoxCHlist->MoveResize(616,65,50,20);

   fTextButtonCHSelect = new TGCheckButton(MassAnalysis_v1,"CHselect");
   fTextButtonCHSelect->SetState(kButtonUp);
   fTextButtonCHSelect->SetTextJustify(36);
   fTextButtonCHSelect->SetMargins(0,0,0,0);
   fTextButtonCHSelect->SetWrapLength(-1);
  // fTextButtonCHSelect->Connect("Clicked()", "MassAnaGUI", this, "ActionSetDriftON()");
   MassAnalysis_v1->AddFrame(fTextButtonCHSelect, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonCHSelect->MoveResize(540,68,75,17);
   //gClient->GetColorByName("#ffffff",ucolor);


   fTextButtonCombineLST = new TGCheckButton(MassAnalysis_v1,"CombineLST");
   fTextButtonCombineLST->SetState(kButtonUp);
   fTextButtonCombineLST->SetTextJustify(36);
   fTextButtonCombineLST->SetMargins(0,0,0,0);
   fTextButtonCombineLST->SetWrapLength(-1);
  // fTextButtonCHSelect->Connect("Clicked()", "MassAnaGUI", this, "ActionSetDriftON()");
   MassAnalysis_v1->AddFrame(fTextButtonCombineLST, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonCombineLST->MoveResize(440,65,96,24);


   // "decorder" and "load tree" button
   fTextButtonDecoder = new TGTextButton(MassAnalysis_v1,"Decoder",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fTextButtonDecoder->SetTextJustify(36);
   fTextButtonDecoder->SetMargins(0,0,0,0);
   fTextButtonDecoder->SetWrapLength(-1);
   fTextButtonDecoder->Resize(120,24);
   fTextButtonDecoder->Connect("Clicked()", "MassAnaGUI", this, "ActivateDecoder()");
   MassAnalysis_v1->AddFrame(fTextButtonDecoder, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonDecoder->MoveResize(440,8,120,24);
   fTextButtonLoadTree = new TGTextButton(MassAnalysis_v1,"Load Tree",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fTextButtonLoadTree->SetTextJustify(36);
   fTextButtonLoadTree->SetMargins(0,0,0,0);
   fTextButtonLoadTree->SetWrapLength(-1);
   fTextButtonLoadTree->Resize(120,24);
   MassAnalysis_v1->AddFrame(fTextButtonLoadTree, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonLoadTree->MoveResize(440,32,120,24);

   // "original" & "correct" button
   fTextButtonOriginal = new TGRadioButton(MassAnalysis_v1,"Run ON");
   fTextButtonOriginal->SetTextJustify(36);
   fTextButtonOriginal->SetMargins(0,0,0,0);
   fTextButtonOriginal->SetWrapLength(-1);
   fTextButtonOriginal->Connect("Clicked()", "MassAnaGUI", this, "ActionSetDriftON()");
   MassAnalysis_v1->AddFrame(fTextButtonOriginal, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonOriginal->MoveResize(576,32,113,24);

   fTextButtonCorrected = new TGRadioButton(MassAnalysis_v1,"Run OFF");
   fTextButtonCorrected->SetState(kButtonDown);
   fTextButtonCorrected->SetTextJustify(36);
   fTextButtonCorrected->SetMargins(0,0,0,0);
   fTextButtonCorrected->SetWrapLength(-1);
   fTextButtonCorrected->Connect("Clicked()", "MassAnaGUI", this, "ActionSetDriftOFF()");
   MassAnalysis_v1->AddFrame(fTextButtonCorrected, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonCorrected->MoveResize(696,32,113,24);

   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   //T0 input
   if(kSuccess)
   fNumberEntryT0 = new TGNumberEntry(MassAnalysis_v1, (Double_t) t0_history,3,-1,(TGNumberFormat::EStyle) 5);
   else fNumberEntryT0 = new TGNumberEntry(MassAnalysis_v1, (Double_t) 0,3,-1,(TGNumberFormat::EStyle) 5);
   fNumberEntryT0->SetName("fNumberEntryT0");
   MassAnalysis_v1->AddFrame(fNumberEntryT0, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNumberEntryT0->MoveResize(735,66,40,20);
   fLabelT0= new TGLabel(MassAnalysis_v1,"T0 [ns]");
   fLabelT0->SetTextJustify(36);
   fLabelT0->SetMargins(0,0,0,0);
   fLabelT0->SetWrapLength(-1);
   MassAnalysis_v1->AddFrame(fLabelT0, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelT0->MoveResize(688,64,40,24);

   // Exit button
   GCValues_t valButtonExit;
   valButtonExit.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#ff0000",valButtonExit.fForeground);
   gClient->GetColorByName("#e8e8e8",valButtonExit.fBackground);
   valButtonExit.fFillStyle = kFillSolid;
   valButtonExit.fFont = ufont->GetFontHandle();
   valButtonExit.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valButtonExit, kTRUE);
   fTextButtonExit = new TGTextButton(MassAnalysis_v1,"Exit",-1,uGC->GetGC(),ufont->GetFontStruct(),kRaisedFrame);
   fTextButtonExit->SetTextJustify(36);
   fTextButtonExit->SetMargins(0,0,0,0);
   fTextButtonExit->SetWrapLength(-1);
   fTextButtonExit->Resize(99,24);
   fTextButtonExit->Connect("Clicked()", "MassAnaGUI", this, "CloseWindow()");
   MassAnalysis_v1->AddFrame(fTextButtonExit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonExit->MoveResize(680,6,99,24);

   // ProgressBar
   TGProgressBar::EBarType ProgressBarType = TGProgressBar::EBarType::kStandard;
   fHProgressBar1546 = new TGHProgressBar(MassAnalysis_v1,ProgressBarType,296);  
   fHProgressBar1546->SetName("fHProgressBar1546");
   fHProgressBar1546->SetFillType(TGProgressBar::kBlockFill);

   ULong_t ucolor;        // will reflect user color changes
   gClient->GetColorByName("#ffffff",ucolor);
   fHProgressBar1546->SetBackgroundColor(ucolor);
   fHProgressBar1546->SetPosition(5.);  // Set Bar "Percentage"!!!!!!!!!!
   fHProgressBar1546->SetBarColor("#0000ff");
   MassAnalysis_v1->AddFrame(fHProgressBar1546, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fHProgressBar1546->MoveResize(120,56,296,17);


ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   // graphics context changes
   GCValues_t valpFrame1634;
   valpFrame1634.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",valpFrame1634.fForeground);
   gClient->GetColorByName("#e8e8e8",valpFrame1634.fBackground);
   valpFrame1634.fFillStyle = kFillSolid;
   valpFrame1634.fFont = ufont->GetFontHandle();
   valpFrame1634.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&valpFrame1634, kTRUE);

   // "Drift Correct" group frame and number for input
   fGroupFrameDC = new TGGroupFrame(MassAnalysis_v1,"Drift Correct",kVerticalFrame,uGC->GetGC());
   fGroupFrameDC->SetLayoutBroken(kTRUE);
   if(kSuccess)
   fNumberEntryRefCento = new TGNumberEntry(fGroupFrameDC, (Double_t) tof_history,6,-1,(TGNumberFormat::EStyle) 0);
   else  fNumberEntryRefCento = new TGNumberEntry(fGroupFrameDC, (Double_t) 1e+07,6,-1,(TGNumberFormat::EStyle) 0);
   fGroupFrameDC->AddFrame(fNumberEntryRefCento, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNumberEntryRefCento->MoveResize(93,16,79,20);
   fNumberEntryNbins = new TGNumberEntry(fGroupFrameDC, (Double_t) 200,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryNbins, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNumberEntryNbins->MoveResize(232,16,58,20);
   fNumberEntryEveAkill = new TGNumberEntry(fGroupFrameDC, (Double_t) 250,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryEveAkill, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fNumberEntryEveAkill->MoveResize(368,16,56,20);
   fNumberEntryRefSigma = new TGNumberEntry(fGroupFrameDC, (Double_t) 10,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryRefSigma, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fNumberEntryRefSigma->MoveResize(504,16,56,20);
   fNumberEntryHisto_halfwidth = new TGNumberEntry(fGroupFrameDC, (Double_t) 150,6,-1,(TGNumberFormat::EStyle) 5);
   fGroupFrameDC->AddFrame(fNumberEntryHisto_halfwidth, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fNumberEntryHisto_halfwidth->MoveResize(658,15,56,20);
   fNumberEntryTveto = new TGNumberEntry(fGroupFrameDC, (Double_t) 0,6,-1,(TGNumberFormat::EStyle) 0);
   fGroupFrameDC->AddFrame(fNumberEntryTveto, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNumberEntryTveto->MoveResize(93,42,79,20);

   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   // graphics context changes
   GCValues_t vall2224;
   vall2224.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000",vall2224.fForeground);
   gClient->GetColorByName("#e8e8e8",vall2224.fBackground);
   vall2224.fFillStyle = kFillSolid;
   vall2224.fFont = ufont->GetFontHandle();
   vall2224.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&vall2224, kTRUE);

   //Drift correct group label and "OK" buttorn
   fLabelRefCenter = new TGLabel(fGroupFrameDC,"Ref. center[ns]",uGC->GetGC());
   fLabelRefCenter->SetTextJustify(36);
   fLabelRefCenter->SetMargins(0,0,0,0);
   fLabelRefCenter->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelRefCenter, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelRefCenter->MoveResize(8,16,80,24);

   fLabelnbins = new TGLabel(fGroupFrameDC,"nbins");
   fLabelnbins->SetTextJustify(36);
   fLabelnbins->SetMargins(0,0,0,0);
   fLabelnbins->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelnbins, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelnbins->MoveResize(174,20,48,16);

   fLabelEveAkill = new TGLabel(fGroupFrameDC,"EveAkill");
   fLabelEveAkill->SetTextJustify(36);
   fLabelEveAkill->SetMargins(0,0,0,0);
   fLabelEveAkill->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelEveAkill, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelEveAkill->MoveResize(312,16,48,24);

   fLabelRefSigma = new TGLabel(fGroupFrameDC,"Sigma");
   fLabelRefSigma->SetTextJustify(36);
   fLabelRefSigma->SetMargins(0,0,0,0);
   fLabelRefSigma->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelRefSigma, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelRefSigma->MoveResize(447,15,56,24);

   fLabelHalf_histoWidth = new TGLabel(fGroupFrameDC,"Half_histoWidth");
   fLabelHalf_histoWidth->SetTextJustify(36);
   fLabelHalf_histoWidth->SetMargins(0,0,0,0);
   fLabelHalf_histoWidth->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelHalf_histoWidth, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelHalf_histoWidth->MoveResize(561,16,96,22);

   fTextButtonOK = new TGTextButton(fGroupFrameDC,"OK",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fTextButtonOK->SetTextJustify(36);
   fTextButtonOK->SetMargins(0,0,0,0);
   fTextButtonOK->SetWrapLength(-1);
   fTextButtonOK->Resize(38,22);
   fTextButtonOK->Connect("Clicked()", "MassAnaGUI", this, "ActivateDrifCorrect()");
   fGroupFrameDC->AddFrame(fTextButtonOK, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButtonOK->MoveResize(724,16,38,22);

   fLabelTveto = new TGLabel(fGroupFrameDC,"Time_veto [ns]",uGC->GetGC());
   fLabelTveto->SetTextJustify(36);
   fLabelTveto->SetMargins(0,0,0,0);
   fLabelTveto->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelTveto, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelTveto->MoveResize(8,42,80,24);

   fLabelSeleTag = new TGLabel(fGroupFrameDC,"Tag",uGC->GetGC());
   fLabelSeleTag->SetTextJustify(36);
   fLabelSeleTag->SetMargins(0,0,0,0);
   fLabelSeleTag->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelSeleTag, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelSeleTag->MoveResize(174,42,48,24);

   fLabelSeleCh = new TGLabel(fGroupFrameDC,"Channel",uGC->GetGC());
   fLabelSeleCh->SetTextJustify(36);
   fLabelSeleCh->SetMargins(0,0,0,0);
   fLabelSeleCh->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelSeleCh, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelSeleCh->MoveResize(312,42,50,24);

   fLabelFitFun = new TGLabel(fGroupFrameDC,"Fit Func.",uGC->GetGC());
   fLabelFitFun->SetTextJustify(36);
   fLabelFitFun->SetMargins(0,0,0,0);
   fLabelFitFun->SetWrapLength(-1);
   fGroupFrameDC->AddFrame(fLabelFitFun, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelFitFun->MoveResize(454,42,50,24);



   fComboBoxTaglist = new TGComboBox(fGroupFrameDC,"Taglist",-1,kHorizontalFrame | kSunkenFrame | kOwnBackground);
   fComboBoxTaglist->SetName("Taglist");
   fComboBoxTaglist->AddEntry("Tag 0 ",0);
   fComboBoxTaglist->AddEntry("Tag 1",1);
   fComboBoxTaglist->Select(1);
   fGroupFrameDC->AddFrame(fComboBoxTaglist, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fComboBoxTaglist->MoveResize(234,42,54,20);


   fComboBoxCHlist_drif = new TGComboBox(fGroupFrameDC,"CHlist_drif",-1,kHorizontalFrame | kSunkenFrame | kOwnBackground);
   fComboBoxCHlist_drif->SetName("CHlist_drif");
   fComboBoxCHlist_drif->AddEntry("CH 1",0);
   fComboBoxCHlist_drif->AddEntry("CH 2",1);
   fComboBoxCHlist_drif->AddEntry("CH 3",2);
   fComboBoxCHlist_drif->AddEntry("CH 4",3);
   fComboBoxCHlist_drif->AddEntry("CH 5",4);
   fComboBoxCHlist_drif->Select(0);
   fGroupFrameDC->AddFrame(fComboBoxCHlist_drif, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fComboBoxCHlist_drif->MoveResize(372,42,50,20);


   fComboBoxFitFunc = new TGComboBox(fGroupFrameDC,"FitFuncList",-1,kHorizontalFrame | kSunkenFrame | kOwnBackground);
   fComboBoxFitFunc->SetName("FitFuncList");
   fComboBoxFitFunc->AddEntry("Johnson",0);
   fComboBoxFitFunc->AddEntry("GausExp",1);
   fComboBoxFitFunc->Select(0);
   fGroupFrameDC->AddFrame(fComboBoxFitFunc, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fComboBoxFitFunc->MoveResize(503,42,74,20);


   fTextButtonFixPars = new TGCheckButton(fGroupFrameDC,"Fix Parameters of Func.");
   fTextButtonFixPars->SetState(kButtonUp);
   fTextButtonFixPars->SetTextJustify(36);
   fTextButtonFixPars->SetMargins(0,0,0,0);
   fTextButtonFixPars->SetWrapLength(-1);
  // fTextButtonCHSelect->Connect("Clicked()", "MassAnaGUI", this, "ActionSetDriftON()");
   fGroupFrameDC->AddFrame(fTextButtonFixPars, new TGLayoutHints(kLHintsLeft |  kLHintsExpandY,2,2,2,2));
   fTextButtonFixPars->MoveResize(590,42,150,24);


 // tab widget
   fTab662 = new TGTab(MassAnalysis_v1,770,348);//fMainFrame1817



   // container of "Tab1"
   TGCompositeFrame *fCompositeFrame665;
   fCompositeFrame665 = fTab662->AddTab("Tab1");
   fCompositeFrame665->SetLayoutManager(new TGVerticalLayout(fCompositeFrame665));

   // embedded canvas
   fRootEmbeddedCanvas676 = new TRootEmbeddedCanvas(0,fCompositeFrame665,762,255,kSunkenFrame);
   Int_t wfRootEmbeddedCanvas676 = fRootEmbeddedCanvas676->GetCanvasWindowId();
   cref = new TCanvas("cref", 10, 40, wfRootEmbeddedCanvas676);
   fRootEmbeddedCanvas676->AdoptCanvas(cref);
   fCompositeFrame665->AddFrame(fRootEmbeddedCanvas676,new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fRootEmbeddedCanvas676->Draw();
 cref->Draw();
   cref->Modified();
   cref->Update();

   fRootEmbeddedCanvas2 = new TRootEmbeddedCanvas(0,fCompositeFrame665,762,176,kSunkenFrame);
  // Int_t wfRootEmbeddedCanvas2 = fRootEmbeddedCanvas2->GetCanvasWindowId();
  // cref2 = new TCanvas("cref2", 10, 40, wfRootEmbeddedCanvas2);
   //fRootEmbeddedCanvas2->AdoptCanvas(cref2);

cref2 =fRootEmbeddedCanvas2->GetCanvas();

   //fRootEmbeddedCanvas2->MoveResize(100,200,720,176);
   fCompositeFrame665->AddFrame(fRootEmbeddedCanvas2,new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,1,2)); // the number and change the positon of canvas (x,unknow,y,unknow)
fRootEmbeddedCanvas2->MoveResize(100,400,762,255); // can not change the position of canvas, only size?
fRootEmbeddedCanvas2->Draw();


 cref2->Draw();
   cref2->Modified();
   cref2->Update();

   // container of "Tab2"
   TGCompositeFrame *fCompositeFrame667;
   fCompositeFrame667 = fTab662->AddTab("Tab2");
   fCompositeFrame667->SetLayoutManager(new TGVerticalLayout(fCompositeFrame667));


   fTab662->SetTab(0);

   fTab662->Resize(fTab662->GetDefaultSize());
   MassAnalysis_v1->AddFrame(fTab662, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTab662->MoveResize(8,166,770,540);


//cref2->AddExec("exe", "((MassAnaGUI*)gROOT->ProcessLine(\"gROOT->GetListOfSpecials()->FindObject(\\\"myui\\\");\"))->cursor();"  );
//gROOT->GetListOfSpecials()->Add(this);


//loop();

   subthread();



   //%%%%%%%%% label of hot key %%%%%%%%%%%%%%%%%%%%%%
 //  TGFont *ufont;         // will reflect user font changes
  // ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

   TGGC   *uGC2;           // will reflect user GC changes
   // graphics context changes
   GCValues_t hotkeyColor;
   hotkeyColor.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#ff0000",hotkeyColor.fForeground);
   gClient->GetColorByName("#e8e8e8",hotkeyColor.fBackground);
   hotkeyColor.fFillStyle = kFillSolid;
   hotkeyColor.fFont = ufont->GetFontHandle();
   hotkeyColor.fGraphicsExposures = kFALSE;
   uGC2 = gClient->GetGC(&hotkeyColor, kTRUE);

   fLabelHotkey = new TGLabel(MassAnalysis_v1,"Hot key: ' z ' -> zoom in , ' x ' upon moving cursor -> zoom out, ' f ' -> fit; cursor click for specifying region",uGC2->GetGC());
   fLabelHotkey->SetTextJustify(kTextLeft | kTextCenterY);
   fLabelHotkey->SetMargins(0,0,0,0);
   //fLabelHotkey->SetTextColor(kRed);
   fLabelHotkey->SetWrapLength(-1);
   MassAnalysis_v1->AddFrame(fLabelHotkey, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabelHotkey->MoveResize(8,710,770,32);
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



   // main frame setting....

   fGroupFrameDC->SetLayoutManager(new TGVerticalLayout(fGroupFrameDC));
   fGroupFrameDC->Resize(770,56);
   MassAnalysis_v1->AddFrame(fGroupFrameDC, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fGroupFrameDC->MoveResize(8,88,770,76);

   fMainFrame1817->AddFrame(MassAnalysis_v1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,1,1,1,1));
   fMainFrame1933->AddFrame(fMainFrame1817, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fMainFrame1817->MoveResize(0,0,876,739);

   fMainFrame2609->AddFrame(fMainFrame1933, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fMainFrame1933->MoveResize(0,0,781,752);


   fMainFrame2609->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
   fMainFrame2609->SetWindowName("MassAnalyser");
   fMainFrame2609->MapSubwindows();

   fMainFrame2609->Resize(fMainFrame2609->GetDefaultSize());
   fMainFrame2609->MapWindow();
   fMainFrame2609->Resize(876,739);


    cref->Draw();
   cref->Modified();
   cref->Update();
 cref2->Draw();
   cref2->Modified();
   cref2->Update();
}  


MassAnaGUI :: ~MassAnaGUI(){

	MassAnalysis_v1->Cleanup();
	fMainFrame1817->Cleanup();
	fMainFrame1933->Cleanup();
	fMainFrame2609->Cleanup();
}


void MassAnaGUI ::ActionSetDriftON(){

	fTextButtonOriginal->SetState(kButtonDown);
	fTextButtonCorrected->SetState(kButtonUp);
   ::abort_corr_loop=false;
	std::cout<<"\e[1;33m"<<"Enable drift correct"<<"\e[0m"<<std::endl;

}


void MassAnaGUI ::ActionSetDriftOFF(){

	fTextButtonOriginal->SetState(kButtonUp);
	fTextButtonCorrected->SetState(kButtonDown);
   ::abort_corr_loop=true;
	std::cout<<"\e[1;33m"<<"Stop dirft correct"<<"\e[0m"<<std::endl;

}

void MassAnaGUI:: CloseWindow(){ // close the main window return to root CINT
   ActionSetDriftOFF();
   fMainFrame2609->DeleteWindow();
   gROOT->ProcessLine(".q");
}

void MassAnaGUI:: ActivateDecoder(){

	bool SelectStat = fTextButtonCHSelect->IsDown();
   bool IsCombine = fTextButtonCombineLST->IsDown();
	int  CHselected = fComboBoxCHlist->GetSelected()+1;

//cout<<"selected stat="<<SelectStat<<" ch sle="<<CHselected<<endl;
	string LstFilePATH = fTextEntryPATH->GetText();
	string LstFilename = fTextEntryFilename->GetText();

   if(LstFilename.find(".lst") != string::npos) LstFilename.erase(LstFilename.begin()+LstFilename.find(".lst"),LstFilename.end());

   if(IsCombine) gROOT->ProcessLine( Form("DecodeAndCombine(\"%s\",\"%s\",%d,%d)",LstFilePATH.c_str(),LstFilename.c_str(),SelectStat,CHselected) );
   else gROOT->ProcessLine( Form("ParserMCS6A_3(\"%s\",\"%s\",%d,%d)",LstFilePATH.c_str(),LstFilename.c_str(),SelectStat,CHselected) );
   
   if(SelectStat)  cout<<"\e[1;33m"<<"selected channel -> "<<CHselected<<"\e[0m"<<endl;
   else cout<<"\e[1;33m"<<"All channels are decoded"<<"\e[0m"<<endl;

}

void MassAnaGUI::ActivateDrifCorrect(){
   ActionSetDriftON();
	string LstFilePATH = fTextEntryPATH->GetText();
	string LstFilename = fTextEntryFilename->GetText();

   if(LstFilename.find(".lst") != string::npos) LstFilename.erase(LstFilename.begin()+LstFilename.find(".lst"),LstFilename.end());

   double T0 = fNumberEntryT0->GetNumber();
   double RefCentoTof = fNumberEntryRefCento->GetNumber();
   int Nbins = fNumberEntryNbins->GetNumber();
   int EveAkill = fNumberEntryEveAkill->GetNumber();
   double RefSigma = fNumberEntryRefSigma->GetNumber();
   double Histo_halfwidth = fNumberEntryHisto_halfwidth->GetNumber();

   double Timeveto = fNumberEntryTveto->GetNumber();
   int SeleTag = fComboBoxTaglist->GetSelected();
   int SeleChannel = fComboBoxCHlist_drif->GetSelected()+1;
   bool fixpar = fTextButtonFixPars->IsDown();

   int SeleFunc = fComboBoxFitFunc->GetSelected();

   if(SeleFunc==0) ::fitname ="john";
   else ::fitname ="fgaus_exp";

cout<<"time veto = "<<Timeveto<<endl;
cout<<"sele tag = "<<SeleTag<<endl;
cout<<"sele channel = "<<SeleChannel<<endl;
cout<<"Fitting function = "<<fitname<<endl;
cout<<"whether fix par = "<<fixpar<<endl;

   ::fit_by_free_pars =!fixpar; //  fit_by_free_pars is defined in DriffCorrect3.C

   fHProgressBar1546->Reset();
   std::cout<<"\e[1;33m"<<"T0 is set to "<<T0<<" [ns]"<<"\e[0m"<<std::endl;
   //const char * command = Form("DriftCorrect(\"%s\",\"%s\",%f,%i,%i,%f,%f,%f,ui->fHProgressBar1546)",LstFilePATH.c_str(),LstFilename.c_str(),RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0); 
   //gROOT->ProcessLine(command);
   //fHProgressBar1546->Reset();
   //if(t!=NULL)delete t;
   tof_and_enteries[0].clear();
   tof_and_enteries[1].clear();
  thread t(&MassAnaGUI ::runthread,this,LstFilePATH,LstFilename,RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0,fHProgressBar1546,SeleTag,SeleChannel,Timeveto,get_drift_corr_pars,tof_and_enteries);
  // thread t(::DriftCorrect,LstFilePATH,LstFilename,RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0,fHProgressBar1546);
  t.join();
  // gROOT->ProcessLine("\n");


      if(href1D!=NULL){delete href1D; href1D=NULL;}
      href1D = new TH1D("href1D","Tof [ns]",(int)get_drift_corr_pars[2],get_drift_corr_pars[3],get_drift_corr_pars[4]);
      href1D->GetXaxis()->SetLabelSize(0.06);
      href1D->GetYaxis()->SetLabelSize(0.06);
      if(href2D!=NULL){delete href2D; href2D=NULL;} 
      href2D = new TH2D("href2D","Tof vs. sweeps",(int)get_drift_corr_pars[1],0,get_drift_corr_pars[0],(int)get_drift_corr_pars[2],get_drift_corr_pars[3],get_drift_corr_pars[4]);
      href2D->GetXaxis()->SetLabelSize(0.06);
      href2D->GetYaxis()->SetLabelSize(0.06);
      href2D->GetZaxis()->SetLabelSize(0.06);

      //cout<<"vctor = "<<tof_and_enteries[0].size()<<endl;

      for(unsigned long index=0;index<tof_and_enteries[0].size();index++){
         href2D->Fill(tof_and_enteries[1][index],tof_and_enteries[0][index]);
         href1D->Fill(tof_and_enteries[0][index]);
      }

      cref->cd();
      href2D->Draw("colz");
      href2D->SetStats(0);
      cref->Modified();
      cref->Update();
      cref2->cd();
      href1D->Draw();
      href1D->SetStats(0);
      cref2->Modified();
      cref2->Update();
      gSystem->ProcessEvents();

   WriteLog();
   gClient->GetDefaultRoot();
   return;

/*
   if(cc!=NULL) delete cc;
   gROOT->SetBatch();
   cc = new TCanvas();gROOT->SetBatch(kFALSE);
   exc(LstFilePATH,LstFilename,RefCentoTof,Nbins,EveAkill,RefSigma,Histo_halfwidth,T0,fHProgressBar1546);
   runcorrect=true;
   cc->AddExec("exc","exc()");
   //gROOT->SetBatch(kFALSE);*/


}
//string PATH="../",string filename = "mcs_39Kvs143X2plus@500_174927", double centro=17e6, int nbins = 200, int event_akill =600, double ref_sigma=10,double Half_hiswidth = 150,double _inT0=130,TGProgressBar *bar


void MassAnaGUI ::runthread(string _PAHT,string _filename,double _centro, int _nbins, int _event_akill, double _ref_sigma,double _Half_hiswidth,double __intTo,TGProgressBar* _bar,
   int _SeleTag, int _SeleCH, double _Timeveto, double * _outpars, vector<double> *_timesANDsweep){

   ::DriftCorrect(_PAHT, _filename, _centro, _nbins, _event_akill, _ref_sigma, _Half_hiswidth, __intTo, _bar, _SeleTag, _SeleCH, _Timeveto, _outpars, _timesANDsweep);
   fHProgressBar1546->RaiseWindow();
   //sleep(1);
  // gROOT->ProcessLine("cout<<\"done\"<<endl;");
     gSystem->ProcessEvents();
    gClient->GetDefaultRoot();
	cout<<"done"<<endl;
   //fHProgressBar1546->Reset();
}

void MassAnaGUI ::exc(string _PAHT,string _filename,double _centro, int _nbins, int _event_akill, double _ref_sigma,double _Half_hiswidth,double __intTo,TGProgressBar* _bar){
   static string _PAHT_;
   static string _filename_;
   static double _centro_;
   static int _nbins_;
   static int _event_akill_;
   static double _ref_sigma_;
   static double _Half_hiswidth_;
   static double __intTo_;
   static TGProgressBar* _bar_;

   if(_PAHT!=""){
      _PAHT_ = _PAHT;
      _filename_ = _filename;
      _centro_ = _centro;
      _nbins_ = _nbins;
      _event_akill_ = _event_akill;
      _ref_sigma_ = _ref_sigma;
      _Half_hiswidth_ = _Half_hiswidth;
      __intTo_ = __intTo;
      _bar_ = _bar;
   }


   if(runcorrect){
      ::DriftCorrect(_PAHT_, _filename_, _centro_, _nbins_, _event_akill_, _ref_sigma_, _Half_hiswidth_, __intTo_, _bar_);
      runcorrect=false;
   }
   else{;}
}


void MassAnaGUI::ReadLog(){
	ifstream fin;
	fin.open("DriftCorrect.log",ios::in);
	if(fin.is_open()){
		fin>>buffer_Path;
		fin>>buffer_filename;
		fin>>buffer_num;
		tof_history = atof(buffer_num);
		fin>>buffer_num;
		t0_history = atof(buffer_num);
		fin.close();
	}
	else{
		gClient->GetDefaultRoot();
		cout<<"No drift history Log"<<endl;
	}

	if(tof_history>0 && tof_history<25000000)kSuccess = true;
	
}



void MassAnaGUI::WriteLog(){

	string LstFilePATH = fTextEntryPATH->GetText();
	string LstFilename = fTextEntryFilename->GetText();

   	t0_history = fNumberEntryT0->GetNumber();
   	tof_history = fNumberEntryRefCento->GetNumber();
      long tof_int = (long) tof_history;

	ofstream fout;
	fout.open("DriftCorrect.log",ios::out | ios::trunc);
	fout<<LstFilePATH<<endl;
	fout<<LstFilename<<endl;
	fout<<tof_int<<endl;
	fout<<t0_history<<endl;
	fout.close();

}


void MassAnaGUI::SetRange(int index,TCanvas* c_get,TH1D* h,double xmin, double xmax){
     h->GetXaxis()->SetRangeUser(xmin,xmax);
     c_get->cd(index)->Modified();
     c_get->cd(index)->Update();      
     c_get->Modified();  c_get->Update();
     return;

}


void MassAnaGUI::SetLinePos(double x1,double y1,double x2,double y2){
   static double xold = 0;
   if(xold==x1) return;
   else{xold = x1;}


   if(lcursor==NULL){ lcursor = new TLine(x1,y1,x2,y2);  cref2->cd();
   lcursor->Draw();}
   else{
      //delete lcursor;
      //lcursor = new TLine(x1,y1,x2,y2);
      lcursor->SetX1(x1);
      lcursor->SetY1(y1);
      lcursor->SetX2(x2);
      lcursor->SetY2(y2);

   }
   cref2->cd();
   lcursor->Draw();
   cref2->Modified();
   cref2->Update();
   return;
}

void MassAnaGUI::cursor(){
   //gSystem->ProcessEvents();

int indexPad=1;
   static double xmouse[2] = {0};
    static int imouse = 0;
    static TH1D* h_history[2] = {};
    static stack<double> lowlimit[2];
    static stack<double> highlimit[2];
    static bool stack2clear[2] = {1,1};

    cref2->FeedbackMode(kTRUE);
    TPad* Pd = (TPad*)cref2->cd();

    if(indexPad==1){
        if(h_history[indexPad-1] != href1D){
          h_history[indexPad-1] = href1D;
          stack2clear[indexPad-1] = true;
        }
        else stack2clear[indexPad-1] = false;

    }

    if(stack2clear[indexPad-1]){
      for(unsigned int i=0;i<lowlimit[indexPad-1].size();i++){
        lowlimit[indexPad-1].pop();
        highlimit[indexPad-1].pop();
      }
      stack2clear[indexPad-1]=false;
    }

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
  static EEventType last_event=kNoEvent;
  static int last_press=0;

  if(px<pxmin) px = pxmin+1;
  if(px>pxmax) px = pxmax-1;
  


  if(pxold) {
  //  gVirtualX->DrawLine(pxold,pymin,pxold,pymax);
   double x_old = Pd->AbsPixeltoX(pxold);
   double ymax = Pd->AbsPixeltoY(pymax);
   double ymin = Pd->AbsPixeltoY(pymin);
   SetLinePos(x_old,ymin,x_old,ymax);
  }

 // gVirtualX->DrawLine(px,pymin,px,pymax);
  //SetLinePos(px,pymin,px,pymax);

  Pd->SetUniqueID(px);
  Double_t upx =  Pd->AbsPixeltoX(px);
  Double_t x =  Pd->PadtoX(upx);
  double Maxy = Pd->AbsPixeltoY(pymax);
  double Miny = Pd->AbsPixeltoY(pymin);

SetLinePos(x,Miny,x,Maxy);

  EEventType event = static_cast<EEventType> ( Pd->GetEvent() );
  

//if(1)SetLinePos(x,Miny,x,Maxy);

/*
  if(last_event == event && event==kKeyPress){
      if(last_press == Pd->GetEventX()) return;       
  }*/

//%%%%%%%%%%%%%% important %%%%%%%%%%%%%%%%
 if(last_event == event) return; //without this; this function is always activated and so the whole program get stuck since the
                               //event of canvas does clear in time; especially for AddExe() method; timer->Connect() method
                              // can survival without this but cause multiple press unitl new event comes like moving mouse.

  last_event = event;
  last_press = Pd->GetEventX();

//cout<<event<<endl;
//cout<<last_press<<endl;

  if(event == kButton1Down){ // select range
    //cout<<"\r x: "<<Form("%.2f",upx)<<flush;
    xmouse[imouse] = upx;  imouse = (imouse+1)%2;
   // cout<<"upx = "<<upx<<endl;
//SetLinePos(x,Miny,x,Maxy);
    //gVirtualX->DrawLine(px,uymin,px,uymax);
  }
  else if(event == kKeyPress){ // Key action
        // int press = Pd->AbsPixeltoX(Pd->GetEventX());
         int press = Pd->GetEventX();    //  cout<<press<<endl;
         if(press=='z'){ // zoom  
//cout<<"zomm in"<<endl;
           if(h_history[indexPad-1] ==NULL) return;
           lowlimit[indexPad-1].push(TMath::Min(xmouse[0],xmouse[1]));
           highlimit[indexPad-1].push(TMath::Max(xmouse[0],xmouse[1]));
           SetRange(indexPad,cref2,h_history[indexPad-1],lowlimit[indexPad-1].top(),highlimit[indexPad-1].top()); 

         }else if(press=='x'){ // unzoom
        //cout<<"zomm out"<<endl;
              if(h_history[indexPad-1] ==NULL) return;
              if(lowlimit[indexPad-1].size() !=0){
                    lowlimit[indexPad-1].pop();
                    highlimit[indexPad-1].pop();
              }
              if(lowlimit[indexPad-1].size() ==0){ 
                    h_history[indexPad-1]->GetXaxis()->UnZoom(); Pd->Modified();Pd->Update();
                    cref2->Modified();
                    cref2->Update();
                    
              }
              else{
                   SetRange(indexPad,cref2,h_history[indexPad-1],lowlimit[indexPad-1].top(),highlimit[indexPad-1].top()); 
              }
         }else if(press=='f'){
               double fit_left_edege= TMath::Min(xmouse[0],xmouse[1]);
               double fit_right_edege= TMath::Max(xmouse[0],xmouse[1]);
               double RefSigma = fNumberEntryRefSigma->GetNumber();

            if(fitname=="john"){
               if(john != NULL) func = john->GetTF1(); // fit function has already been created
               else{
                  john = new funcJohnson();
                  func = john->GetTF1();
                  func->SetName("john");
               }
            }else if(fitname=="fgaus_exp"){
               if(fgaus_exp != NULL) func = fgaus_exp; // fit function has already been created
               else{
                  fgaus_exp = new TF1("fgaus_exp",fitfunc,0,26e6,5);
                  fgaus_exp->SetParLimits(2,0,3*RefSigma);
                  fgaus_exp->SetParLimits(3,0,30000);
                  fgaus_exp->SetParLimits(4,-30000,0);
                  fgaus_exp->SetParameters(0,0,RefSigma,RefSigma,-RefSigma);
                  func = fgaus_exp;
                  func->SetName("fgaus_exp");
               }
            }  


       for(int ipar=0;ipar<func->GetNpar();ipar++){
            func->ReleaseParameter(ipar);
      }
double amp_tem = h_history[indexPad-1]->GetMaximum();
double peak_cen = (fit_left_edege+fit_right_edege)/2;
               if(fitname=="john") func->SetParameters(amp_tem,peak_cen,0,0,-1);
               else if(fitname=="fgaus_exp") func->SetParameters(amp_tem,peak_cen,RefSigma,RefSigma,-RefSigma);

               for(int i=0;i<10;i++){
                  if(i<5){
                     if(i%2==0){
                         h_history[indexPad-1]->Fit(func,"LMQ","",fit_left_edege,fit_right_edege);
                         continue;
                     }
                     else{
                        h_history[indexPad-1]->Fit(func,"MQ","",fit_left_edege,fit_right_edege);
                        continue;
                     }
                  }
                  h_history[indexPad-1]->Fit(func,"LMQ","",fit_left_edege,fit_right_edege);
               }
      
            double peakcenter = func->GetMaximumX(fit_left_edege,fit_right_edege,1e-13,500); 

            ::fit_half_rangeL=peakcenter-fit_left_edege;
            ::fit_half_rangeR=fit_right_edege-peakcenter;

            double MaxY = func->GetMaximum(fit_left_edege,fit_right_edege,1e-13,500);
            double FWHM_L = func->GetX(MaxY*0.5,peakcenter-10000,peakcenter);
            double FWHM_R = func->GetX(MaxY*0.5,peakcenter,peakcenter+10000);
            cout<<"Rm = "<<peakcenter/( (FWHM_R-FWHM_L)*2)<<endl;

            cref2->Modified();
            cref2->Update();

         }

    }

    //Pd->GetCanvas()->HandleInput(event,eventx,eventy);
    cref2->Modified();
   cref2->Update();
  //  cref->cd();


}

void MassAnaGUI::subthread(){

timer = new TTimer(50);
timer->Connect("Timeout()","MassAnaGUI",this,"cursor()");
timer->Start();

}

void MassAnaGUI::loop(){
   t2 = new thread(&MassAnaGUI::subthread,this);
   t2->detach();
}



MassAnaGUI * ui = NULL;

MassAnaGUI * MassAnaGUI4(){
	gROOT->ProcessLine(".L ParserMCS6A_3.C+");
//	gROOT->ProcessLine(".L DriftCorrect2.C+");
	ui = new MassAnaGUI();
	return ui;

}
