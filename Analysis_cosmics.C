
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <TCanvas.h>
#include <iostream>

#include "TVector3.h"

using namespace std;

void Analysis_cosmics()
{

  //Flags
  Bool_t TrackVisuFlag=0;
  Bool_t MultiFlag=1;
  Bool_t ChargeFlag=1;
    
  if(TrackVisuFlag){
    Can_track=new TCanvas("Can_track","Can_track",900,900);
    //Can_track->Divide(1,2);
    //Can_track->cd(2)->Divide(2,1);
  }

  gSystem->Load("libactar.sl");

  Bool_t StopTrack,CutPadE,CutTrackE;

  Bool_t simuFlag;
  Bool_t gasflag;
  Char_t *simname;
  Char_t *digname;
  Char_t *gasname;

  //simname="./root_files/simFile.root";
  //simname="./root_files/simFile_HeiC4H10.root";
  //simname="./root_files/simFile_10M.root";
  //simname="./root_files/simFile_10M_GasTrig.root";
  //simname="./root_files/simFile_GasTrig.root";
  //simname="./root_files/simFile_GasTrig_NoDiam.root";
  //simname="./root_files/simFile_world3m.root";
  //simname="./root_files/simFile_cryZinv.root";
  //simname="./root_files/simFile_crybox1.root";
  //simname="./root_files/simFile_gasboxpad.root";
  //simname="./root_files/simFile_gasboxpad_cryallinv.root";
  //simname="./root_files/simFile_nogaschamber_noentrancewin.root";
  //simname="./root_files/simFile_HeIso5_GasTrig.root";
  //simname="./root_files/simFile_HeIso5_GasTrig_CRY500x268.root";
  //simname="./root_files/simFile_HeIso5_GasTrig_CRY500x500.root";
  //simname="./root_files/simFile_HeIso5_GasSciTrig_GoodLat.root";
  simname="./root_files/simFile_HeIso5_GasTrig_GoodLat.root";

  //Event info;
  TFile *simFile=new TFile(simname);

  //getting the trees
  TTree *simTree=(TTree*)simFile->Get("The_ACTAR_Event_Tree");
  Int_t nentries=simTree->GetEntries();
  cout<<"Number of sim event : "<<nentries<<endl;
  //TTree *T=(TTree*)simFile->Get("digiTree");
  /*
  //ClonesArray to the silicon Hits
  TClonesArray *silHitsCA=new TClonesArray("ActarSimSilHit",200);
  TBranch *branchSilHits=simTree->GetBranch("silHits");
  branchSilHits->SetAddress(&silHitsCA);
  branchSilHits->SetAutoDelete(kTRUE);
  ActarSimSilHit *silHit=new ActarSimSilHit;
  */     
  //ClonesArray to the scintillator Hits
  TClonesArray *sciHitsCA=new TClonesArray("ActarSimSciHit",200);
  TBranch *branchSciHits=simTree->GetBranch("sciHits");
  branchSciHits->SetAddress(&sciHitsCA);
  branchSciHits->SetAutoDelete(kTRUE);
  ActarSimSciHit *sciHit=new ActarSimSciHit;
         
  //cout<<"digFile to use: ";
  //cin >> digname;

  //digname="./dig_files/digFile.root";
  //digname="./dig_files/digFile_HeiC4H10.root";
  //digname="./dig_files/digFile_10M.root";
  //digname="./dig_files/digFile_10M_2.root";
  //digname="./dig_files/digFile_10M_GasTrig.root";
  //digname="./dig_files/digFile_GasTrig.root";
  //digname="./dig_files/digFile_GasTrig_NoDiam.root";
  //digname="./dig_files/digFile_world3m.root";
  //digname="./dig_files/digFile_world3m_2.root";
  //digname="./dig_files/digFile_cryZinv.root";
  //digname="./dig_files/digFile_crybox1.root";
  //digname="./dig_files/digFile_gasboxpad.root";
  //digname="./dig_files/digFile_gasboxpad_cryallinv.root";
  //digname="./dig_files/digFile_nogaschamber_noentrancewin.root";
  //digname="./dig_files/digFile_HeIso5_GasTrig.root";
  //digname="./dig_files/digFile_HeIso5_GasTrig_CRY500x268.root";
  //digname="./dig_files/digFile_HeIso5_GasTrig_CRY500x500.root";
  //digname="./dig_files/digFile_HeIso5_GasSciTrig_GoodLat.root";
  digname="./dig_files/digFile_HeIso5_GasTrig_GoodLat.root";

  // cout<<"Gas: isobutane (1) or deuterium gas (0)? ";
  // cin >> gasflag;
  // if(gasflag==1)gasname="isobutane";
  // if(gasflag==0)gasname="deuterium";

  gROOT->ProcessLine(".L digit_piotr.h+");
  //gROOT->ProcessLine(".L digit.h+");
  
  padsGeometry thePadsGeometry;
  thePadsGeometry.SetGeometryValues(37.,85.,69.,2.,5.,5.);//digit+
  //thePadsGeometry.SetGeometryValues(32.,85.,64.,2.,0.,0.);//digit+
  //thePadsGeometry.SetGeometryValues(0,0,0,32.,85.,64.,100.,2);
  driftManager theDriftManager;

  Char_t dummy[256];
  Char_t gas[256];
  Double_t v_drift,sigma_trans,sigma_long;
  ifstream *gasfile=new ifstream(gasname);
  gasfile->getline(gas,256);
  cout<<gas<<endl;
  gasfile->getline(dummy,256);
  cout<<dummy<<endl;
  *gasfile>>dummy>>dummy>>v_drift;
  *gasfile>>dummy>>dummy>>sigma_trans;
  *gasfile>>dummy>>dummy>>sigma_long;

  //theDriftManager.SetDriftParameters(2015.,170.,147.5,gasname);
  //theDriftManager.SetDriftParameters(2015.,170.,980.66,gasname);

  //Magboltz Drift paramters for Deuterium Gas
  //theDriftManager.SetDriftVelocity(4.7e-3);
  //theDriftManager.SetDiffusionParameters(1.146e-5,2.342e-5);

  //Magboltz Drift paramters for HeiC4H10 (9:1) Gas
  //theDriftManager.SetDriftVelocity(9.084e-3);
  //theDriftManager.SetDiffusionParameters(2.356e-5,3.105e-5);

  //#######HeiC4H10 (9:1) Gas MagBoltz Values E=2090/17#########
  theDriftManager.SetDriftVelocity(6.865e-3);
  theDriftManager.SetDiffusionParameters(2.369e-5,2.892e-5);

  cout<<"Drift Parameters are:"<<endl;  
  cout<<"v_drift---------> "<<theDriftManager.GetDriftVelocity()<<"mm/ns"<<endl;  
  cout<<"D_long----------> "<<theDriftManager.GetLongitudinalDiffusion()<<"mm^2/ns"<<endl;  
  cout<<"D_trans---------> "<<theDriftManager.GetTransversalDiffusion()<<"mm^2/ns"<<endl;  

  Double_t padSize  = thePadsGeometry.GetPadSize();
  Double_t xLength  = thePadsGeometry.GetXLength();
  Double_t yLength  = thePadsGeometry.GetYLength();
  Double_t zLength  = thePadsGeometry.GetZLength();
 
  cout<<"X length===> "<<xLength<<endl;
  cout<<"Y length===> "<<yLength<<endl;
  cout<<"Z length===> "<<zLength<<endl;

  const Int_t numberOfRows   = thePadsGeometry.GetNumberOfRows();
  const Int_t numberOfColumns= thePadsGeometry.GetNumberOfColumns();

  cout<<"Number of Rows: "<<numberOfRows<<", number of Columns: "<<numberOfColumns<<endl;

  Double_t driftVelocity = theDriftManager.GetDriftVelocity();

  TFile *digFile=new TFile(digname);
  cout<<"Opening digitization file: "<<digname<<endl;

  TTree *digiTree=(TTree*)digFile->Get("digiTree");
  Int_t dentries=digiTree->GetEntries();
  cout<<"Number of digit event : "<<dentries<<endl;
  //Int_t digentries=digiTree->GetEntries();
  //cout<<"Number of digit event : "<<digentries<<endl;
 
  //ClonesArray to the signal
  TClonesArray *padSignalCA=new TClonesArray("ActarPadSignal",4000);
  digiTree->SetBranchAddress("padSignals",&padSignalCA);
  ActarPadSignal *padSignal=new ActarPadSignal;

  Double_t PI=3.1415926535897932384626433;
  Double_t deg=180./PI;
  Double_t rad=PI/180.;
  //==================================================================================//

  //read the Tree generated by tree1w and fill two histograms
   
  //note that we use "new" to create the TFile and TTree objects !
  //because we want to keep these objects alive when we leave this function.

  ///////////////////////////////////////////////
  ///////// HISTOGRAMS //////////////////////////
  ///////////////////////////////////////////////

  TH2F *visu_charge=new TH2F("visu_charge","visu_charge",64,0,64,32,0,32);
  TH2F *visu_time=new TH2F("visu_time","visu_time",64,0,64,32,0,32);

  TH2F *h_PadChargeCumul_Gas=new TH2F("h_PadChargeCumul_Gas","h_PadChargeCumul_Gas",64,0,64,32,0,32);
  TH2F *h_PadMultiplicity_Gas=new TH2F("h_PadMultiplicity_Gas","h_PadMultiplicity_Gas",64,0,64,32,0,32);
  TH2F *h_PadQAv_Gas=new TH2F("h_PadQAv_Gas","h_PadQAv_Gas",64,0,64,32,0,32);

  TH2F *h_PadChargeCumul_SciGas=new TH2F("h_PadChargeCumul_SciGas","h_PadChargeCumul_SciGas",64,0,64,32,0,32);
  TH2F *h_PadMultiplicity_SciGas=new TH2F("h_PadMultiplicity_SciGas","h_PadMultiplicity_SciGas",64,0,64,32,0,32);
  TH2F *h_PadQAv_SciGas=new TH2F("h_PadQAv_SciGas","h_PadQAv_SciGas",64,0,64,32,0,32);

  TH1F *h_PadQTot_Gas=new TH1F("h_PadQTot_Gas","Pad Sum Spectrum. Mesh trigger.",2000,0,2000);
  TH1F *h_AllPadQ_Gas=new TH1F("h_AllPadQ_Gas","All Pads Charge Spectrum. Mesh trigger.",500,0,500);
  TH1F *h_CentPadQTot_Gas=new TH1F("h_CentPadQTot_Gas","Pad Row 16 Col 32 Charge Spectrum. Mesh trigger.",1000,0,100);
  TH1F *h_ExcPadQTot_Gas=new TH1F("h_ExcPadQTot_Gas","Pad Row 16 Col 0 Charge Spectrum. Mesh trigger.",1000,0,100);

  TH1F *h_PadQTot_SciGas=new TH1F("h_PadQTot_SciGas","Pad Sum Spectrum. Sci+Mesh trigger.",2000,0,2000);
  TH1F *h_AllPadQ_SciGas=new TH1F("h_AllPadQ_SciGas","All Pads Charge Spectrum. Sci+Mesh trigger.",500,0,500);
  TH1F *h_CentPadQTot_SciGas=new TH1F("h_CentPadQTot_SciGas","Pad Row 16 Col 32 Charge Spectrum. Sci+Mesh trigger.",1000,0,100);
  TH1F *h_ExcPadQTot_SciGas=new TH1F("h_ExcPadQTot_SciGas","Pad Row 16 Col 0 Charge Spectrum. Sci+Mesh trigger.",1000,0,100);

  TH2F *h_PadChargeCumulCutPad_Gas=new TH2F("h_PadChargeCumulCutPad_Gas","h_PadChargeCumulCutPad_Gas",64,0,64,32,0,32);
  TH2F *h_PadMultiplicityCutPad_Gas=new TH2F("h_PadMultiplicityCutPad_Gas","h_PadMultiplicityCutPad_Gas",64,0,64,32,0,32);
  TH2F *h_PadQAvCutPad_Gas=new TH2F("h_PadQAvCutPad_Gas","h_PadQAvCutPad_Gas",64,0,64,32,0,32);

  TH2F *h_PadChargeCumulCutTrack_Gas=new TH2F("h_PadChargeCumulCutTrack_Gas","h_PadChargeCumulCutTrack_Gas",64,0,64,32,0,32);
  TH2F *h_PadMultiplicityCutTrack_Gas=new TH2F("h_PadMultiplicityCutTrack_Gas","h_PadMultiplicityCutTrack_Gas",64,0,64,32,0,32);
  TH2F *h_PadQAvCutTrack_Gas=new TH2F("h_PadQAvCutTrack_Gas","h_PadQAvCutTrack_Gas",64,0,64,32,0,32);

  TH2F *h_PadChargeCumulCutPadTrack_Gas=new TH2F("h_PadChargeCumulCutPadTrack_Gas","h_PadChargeCumulCutPadTrack_Gas",64,0,64,32,0,32);
  TH2F *h_PadMultiplicityCutPadTrack_Gas=new TH2F("h_PadMultiplicityCutPadTrack_Gas","h_PadMultiplicityCutPadTrack_Gas",64,0,64,32,0,32);
  TH2F *h_PadQAvCutPadTrack_Gas=new TH2F("h_PadQAvCutPadTrack_Gas","h_PadQAvCutPadTrack_Gas",64,0,64,32,0,32);

  ///////////////////////////////////////////////

  Int_t nbsci,numberofpads;
  Double_t Qtot;

  Int_t Rmin=0, Rmax=31;
  Int_t Cmin=0, Cmax=63;

  Double_t threshold = 0;
  Double_t Tthreshold = 1.;

  Double_t charge, Ccharge;
  Double_t Rm,Cm,Zcor,Zm;

  //Matrix for the charge map
  Double_t **padCharge=new Double_t*[numberOfRows];
  Double_t **padTime=new Double_t*[numberOfRows];
  Double_t **padHeight=new Double_t*[numberOfRows];
  for(Int_t j=0;j<numberOfRows;j++){
    padCharge[j]=new Double_t[numberOfColumns];
    padTime[j]=new Double_t[numberOfColumns];
  }
   
  //*********************************************************************************************************//
  //*********************************************************************************************************//
  //**                                                                                                     **//
  //**                                        Event Loop                                                   **//
  //**                                                                                                     **//
  //*********************************************************************************************************//
  //*********************************************************************************************************//

  //for (Long64_t jentry=0;jentry<nentries;jentry++) {
  for (Long64_t jentry=0;jentry<dentries;jentry++) {
  //for (Long64_t jentry=0;jentry<10000;jentry++) {
    //for (Long64_t jentry=5000;jentry<nentries;jentry++) {
    //for (Long64_t jentry=15000;jentry<nentries;jentry++) {
    if(jentry%1000==0)cout<<jentry<<endl;
    //if(jentry%2==0)cout<<"¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡   NEW (2nd) EVENT : "<<jentry<<"   !!!!!!!!!!!!!!!!!!!!"<<endl;

    for(Int_t r=0;r<numberOfRows;r++){
      for(Int_t c=0;c<numberOfColumns;c++){
	padCharge[r][c]=0.;
	padTime[r][c]=0.;
      }
    }
    nbsci = numberofpads = 0;
    StopTrack=0;
    CutPadE=0;
    CutTrackE=0;
    Qtot=0;

    visu_charge->Reset();
    visu_time->Reset();

    sciHitsCA->Clear();
    simTree->GetEvent(jentry);

    nbsci = sciHitsCA->GetEntries();
    //cout<<" SCINT "<<nbsci<<endl;

    padSignalCA->Clear();

    digiTree->GetEvent(jentry);
    
    numberofpads = padSignalCA->GetEntries();
    //cout<<"# Pads fired "<< numberofpads<<endl;

    //if(numberofpads>0 && nbsci>0){
    if(numberofpads>0){

      for(Int_t k=0;k<numberofpads;k++){
	padSignal=(ActarPadSignal*)padSignalCA->At(k);
	Double_t thisCharge = padSignal->GetChargeDeposited();
	Double_t thisTime   = (padSignal->GetInitTime()+padSignal->GetFinalTime())/2.;
	Double_t thisSigmaTime=padSignal->GetSigmaTime();
	Int_t PadRow=padSignal->GetPadRow();
	Int_t PadColumn=padSignal->GetPadColumn();
	//if(k%20==0)cout<<"thisSigmaTime: "<<thisSigmaTime<<endl;
	if (thisCharge!=0){
	    padCharge[padSignal->GetPadRow()-1][padSignal->GetPadColumn()-1]
	      += thisCharge;
	    padTime[padSignal->GetPadRow()-1][padSignal->GetPadColumn()-1]
	      += thisCharge*thisTime;
	    //if(thisCharge>1000000){ StopTrack=1; }//cout<<"BIG CHARGE "<<thisCharge<<" event "<<jentry<<endl;}
	    //if(padCharge[padSignal->GetPadRow()-1][padSignal->GetPadColumn()-1]*30/1e6>3){CutPadE=1; }//cout<<"BIG CHARGE "<<thisCharge<<" event "<<jentry<<endl;}
	    Qtot+=thisCharge*30/1e6;
	}
      }//Loop on ActarPadSignals

      if(Qtot>10)CutTrackE=1;

      /*
	      if(padChargeTest(r,c)*30/1e6<20)PadESig=0.085;
	      else if(padChargeTest(r,c)*30/1e6>=20 && padChargeTest(r,c)*30/1e6<=5000)PadESig=-1.6215*1e-5*(padChargeTest(r,c)*30/1e6-20)+0.085;
	      else if(padChargeTest(r,c)*30/1e6>5000)PadESig=0.00425;

	      padCharge[r][c]=gRandom->Gaus(padChargeTest(r,c),PadESig*padChargeTest(r,c));
	      padCharge[r][c]*=30/1e9;
      */

      //cout<<"Qtot "<<Qtot2<<" keV"<<endl;

      //if(!StopTrack){
	// Loop on rows & columns
	for(Int_t c=0;c<numberOfColumns;c++){ //LOOP on Columns
	  for(Int_t r=0;r<numberOfRows;r++){ //LOOP on Rows
	    padTime[r][c]=padTime[r][c]/padCharge[r][c];
	    //padTime[r][c]=gRandom->Gaus(padTime[r][c],sigma_time);
	    padCharge[r][c]*=30/1e6;//energy in keV

	    CutPadE=0;
	    if(padCharge[r][c]>0.1)CutPadE=1;

	    if(padCharge[r][c]>0){
	      //Qtot+=padCharge[r][c];
	      h_AllPadQ_Gas->Fill(padCharge[r][c]);
	      if(r==16 && c==32)h_CentPadQTot_Gas->Fill(padCharge[r][c]);
	      if(r==16 && c==0)h_ExcPadQTot_Gas->Fill(padCharge[r][c]);

	      visu_charge->Fill(c+0.5,r+0.5,padCharge[r][c]);
	      visu_time->Fill(c+0.5,r+0.5,padTime[r][c]);

	      h_PadChargeCumul_Gas->Fill(c+0.5,r+0.5,padCharge[r][c]);
	      h_PadMultiplicity_Gas->Fill(c+0.5,r+0.5,1);

	      //if(!CutPadE){
	      if(CutPadE){
		h_PadChargeCumulCutPad_Gas->Fill(c+0.5,r+0.5,padCharge[r][c]);
		h_PadMultiplicityCutPad_Gas->Fill(c+0.5,r+0.5,1);
	      }

	      if(CutTrackE){
		h_PadChargeCumulCutTrack_Gas->Fill(c+0.5,r+0.5,padCharge[r][c]);
		h_PadMultiplicityCutTrack_Gas->Fill(c+0.5,r+0.5,1);
	      }

	      if(CutPadE && CutTrackE){
		h_PadChargeCumulCutPadTrack_Gas->Fill(c+0.5,r+0.5,padCharge[r][c]);
		h_PadMultiplicityCutPadTrack_Gas->Fill(c+0.5,r+0.5,1);
	      }

	      if(nbsci>0){
		h_PadChargeCumul_SciGas->Fill(c+0.5,r+0.5,padCharge[r][c]);
		h_PadMultiplicity_SciGas->Fill(c+0.5,r+0.5,1);
		h_AllPadQ_SciGas->Fill(padCharge[r][c]);
		if(r==16 && c==32)h_CentPadQTot_SciGas->Fill(padCharge[r][c]);
		if(r==16 && c==0)h_ExcPadQTot_SciGas->Fill(padCharge[r][c]);
	      }
	    }
	  }
	}//End of Loop on rows & columns

	//cout<<"Total charge on pads: "<<Qtot<<" keV"<<endl;

	h_PadQTot_Gas->Fill(Qtot);
	if(nbsci>0)h_PadQTot_SciGas->Fill(Qtot);

	//}//end of cuts (StopTrack)
    }
   
    if(TrackVisuFlag){
      //if(TrackVisuFlag && StopTrack){

	Can_track->cd();
	visu_charge->Draw("colz");

	/*
	Can_track->cd(2);
	visu_time->Draw("colz");
	*/
	Can_track->Update();
	Can_track->WaitPrimitive();
      }

  }
    
  //*********************************************************************************************************//
  //*********************************************************************************************************//
  //**                                                                                                     **//
  //**                                       End of Event Loop                                             **//
  //**                                                                                                     **//
  //*********************************************************************************************************//
  //*********************************************************************************************************//



    if(MultiFlag){

      Double_t QAv_Gas,QAvCutPad_Gas,QAvCutTrack_Gas,QAvCutPadTrack_Gas;
      Double_t QAv_SciGas;

      for(Int_t c=0;c<numberOfColumns;c++){ //LOOP on Columns
	for(Int_t r=0;r<numberOfRows;r++){ //LOOP on Rows

	  QAv_Gas=QAvCutPad_Gas=QAvCutTrack_Gas=QAvCutPadTrack_Gas=0;
	  QAv_SciGas=0;

	  /*
	  if(h_PadMultiplicity_Gas->GetBinContent(c+0.5,r+0.5)!=0)QAv_Gas=h_PadChargeCumul_Gas->GetBinContent(c+0.5,r+0.5)/h_PadMultiplicity_Gas->GetBinContent(c+0.5,r+0.5);
	  else QAv_Gas=0;

	  if(h_PadMultiplicity_SciGas->GetBinContent(c+0.5,r+0.5)!=0)QAv_SciGas=h_PadChargeCumul_SciGas->GetBinContent(c+0.5,r+0.5)/h_PadMultiplicity_SciGas->GetBinContent(c+0.5,r+0.5);
	  else QAv_SciGas=0;
	  */
	  
	  if(h_PadMultiplicity_Gas->GetBinContent(c+1,r+1)!=0)QAv_Gas=h_PadChargeCumul_Gas->GetBinContent(c+1,r+1)/h_PadMultiplicity_Gas->GetBinContent(c+1,r+1);
	  else QAv_Gas=0;

	  if(h_PadMultiplicity_SciGas->GetBinContent(c+1,r+1)!=0)QAv_SciGas=h_PadChargeCumul_SciGas->GetBinContent(c+1,r+1)/h_PadMultiplicity_SciGas->GetBinContent(c+1,r+1);
	  else QAv_SciGas=0;
		  
	  if(h_PadMultiplicityCutPad_Gas->GetBinContent(c+1,r+1)!=0)QAvCutPad_Gas=h_PadChargeCumulCutPad_Gas->GetBinContent(c+1,r+1)/h_PadMultiplicityCutPad_Gas->GetBinContent(c+1,r+1);
	  else QAvCutPad_Gas=0;
  		  
	  if(h_PadMultiplicityCutTrack_Gas->GetBinContent(c+1,r+1)!=0)QAvCutTrack_Gas=h_PadChargeCumulCutTrack_Gas->GetBinContent(c+1,r+1)/h_PadMultiplicityCutTrack_Gas->GetBinContent(c+1,r+1);
	  else QAvCutTrack_Gas=0;
  		  
	  if(h_PadMultiplicityCutPadTrack_Gas->GetBinContent(c+1,r+1)!=0)QAvCutPadTrack_Gas=h_PadChargeCumulCutPadTrack_Gas->GetBinContent(c+1,r+1)/h_PadMultiplicityCutPadTrack_Gas->GetBinContent(c+1,r+1);
	  else QAvCutPadTrack_Gas=0;
  		  
	  //if(c<10 && r>15)cout<<"Row "<<r<<" Col "<<c<<" PadChargeCumul "<<h_PadChargeCumul_SciGas->GetBinContent(c+0.5,r+0.5)<<" PadMultiplicity "<<h_PadMultiplicity_SciGas->GetBinContent(c+0.5,r+0.5)<< " QAv "<<QAv_SciGas<<endl;

	  h_PadQAv_Gas->Fill(c,r,QAv_Gas);
	  h_PadQAv_SciGas->Fill(c,r,QAv_SciGas);
	  h_PadQAvCutPad_Gas->Fill(c,r,QAvCutPad_Gas);
	  h_PadQAvCutTrack_Gas->Fill(c,r,QAvCutTrack_Gas);
	  h_PadQAvCutPadTrack_Gas->Fill(c,r,QAvCutPadTrack_Gas);
	  //h_PadQAv_Gas->Fill(c+0.5,r+0.5,QAv_Gas);
	  //h_PadQAv_SciGas->Fill(c+0.5,r+0.5,QAv_SciGas);
	  //h_PadQAv_Gas->Fill(c-0.5,r-0.5,QAv_Gas);
	  //h_PadQAv_SciGas->Fill(c-0.5,r-0.5,QAv_SciGas);
	}
      }//End of Loop on rows & columns

      
      TCanvas* Can_multi_Gas=new TCanvas("Can_multi_Gas","Can_multi_Gas",900,900);
      Can_multi_Gas->Divide(1,2);
      Can_multi_Gas->cd(1)->Divide(2,1);

      //Can_multi_Gas->Divide(2,1);
      Can_multi_Gas->cd(1)->cd(1);
      h_PadChargeCumul_Gas->Draw("colz");
      Can_multi_Gas->cd(1)->cd(2);
      h_PadMultiplicity_Gas->Draw("colz");

      Can_multi_Gas->cd(2);
      h_PadQAv_Gas->Draw("colz");
          
      TCanvas* Can_multiCutPad_Gas=new TCanvas("Can_multiCutPad_Gas","Can_multiCutPad_Gas",900,900);
      Can_multiCutPad_Gas->Divide(1,2);
      Can_multiCutPad_Gas->cd(1)->Divide(2,1);

      Can_multiCutPad_Gas->cd(1)->cd(1);
      h_PadChargeCumulCutPad_Gas->Draw("colz");
      Can_multiCutPad_Gas->cd(1)->cd(2);
      h_PadMultiplicityCutPad_Gas->Draw("colz");

      Can_multiCutPad_Gas->cd(2);
      h_PadQAvCutPad_Gas->Draw("colz");
      
      TCanvas* Can_multiCutTrack_Gas=new TCanvas("Can_multiCutTrack_Gas","Can_multiCutTrack_Gas",900,900);
      Can_multiCutTrack_Gas->Divide(1,2);
      Can_multiCutTrack_Gas->cd(1)->Divide(2,1);

      Can_multiCutTrack_Gas->cd(1)->cd(1);
      h_PadChargeCumulCutTrack_Gas->Draw("colz");
      Can_multiCutTrack_Gas->cd(1)->cd(2);
      h_PadMultiplicityCutTrack_Gas->Draw("colz");

      Can_multiCutTrack_Gas->cd(2);
      h_PadQAvCutTrack_Gas->Draw("colz");
          
      TCanvas* Can_multiCutPadTrack_Gas=new TCanvas("Can_multiCutPadTrack_Gas","Can_multiCutPadTrack_Gas",900,900);
      Can_multiCutPadTrack_Gas->Divide(1,2);
      Can_multiCutPadTrack_Gas->cd(1)->Divide(2,1);

      Can_multiCutPadTrack_Gas->cd(1)->cd(1);
      h_PadChargeCumulCutPadTrack_Gas->Draw("colz");
      Can_multiCutPadTrack_Gas->cd(1)->cd(2);
      h_PadMultiplicityCutPadTrack_Gas->Draw("colz");

      Can_multiCutPadTrack_Gas->cd(2);
      h_PadQAvCutPadTrack_Gas->Draw("colz");
      
      /*
      TCanvas* Can_multi_SciGas=new TCanvas("Can_multi_SciGas","Can_multi_SciGas",900,900);
      Can_multi_SciGas->Divide(1,2);
      Can_multi_SciGas->cd(1)->Divide(2,1);

      //Can_multi_SciGas->Divide(2,1);
      Can_multi_SciGas->cd(1)->cd(1);
      h_PadChargeCumul_SciGas->Draw("colz");
      Can_multi_SciGas->cd(1)->cd(2);
      h_PadMultiplicity_SciGas->Draw("colz");

      Can_multi_SciGas->cd(2);
      h_PadQAv_SciGas->Draw("colz");
      */
    }

    if(ChargeFlag){
    
      TCanvas* Can_charge_Gas=new TCanvas("Can_charge_Gas","Can_charge_Gas",900,900);
      Can_charge_Gas->Divide(2,2);
      //Can_charge_Gas->cd(2)->Divide(2,1);

      Can_charge_Gas->cd(1);
      h_PadQTot_Gas->Draw();
      Can_charge_Gas->cd(2);
      h_AllPadQ_Gas->Draw();
      Can_charge_Gas->cd(3);
      h_CentPadQTot_Gas->Draw();
      Can_charge_Gas->cd(4);
      h_ExcPadQTot_Gas->Draw();
      /*
      TCanvas* Can_charge_SciGas=new TCanvas("Can_charge_SciGas","Can_charge_SciGas",900,900);
      Can_charge_SciGas->Divide(2,2);
      //Can_charge_SciGas->cd(2)->Divide(2,1);

      Can_charge_SciGas->cd(1);
      h_PadQTot_SciGas->Draw();
      Can_charge_SciGas->cd(2);
      h_AllPadQ_SciGas->Draw();
      Can_charge_SciGas->cd(3);
      h_CentPadQTot_SciGas->Draw();
      Can_charge_SciGas->cd(4);
      h_ExcPadQTot_SciGas->Draw();
      */
    }

} 
