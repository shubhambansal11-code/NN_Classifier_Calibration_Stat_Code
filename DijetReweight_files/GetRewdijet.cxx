#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TF1.h"
using namespace RooFit ;

void GetRewdijet()
{ 
  TFile *f_Rew_dij = new TFile("RewMC_SB_pT610_3rdPoly.root", "RECREATE");
  //TFile *f_trial = new TFile("DataOvDij.root", "RECREATE");  
  
  TString dij = "QCD_MergeHist.root";


  TFile File_dij(base+dij);
  TH1D *h_dij_M=(TH1D*)File_dij.Get("Zcand_Xbb60_mass_pT610");

  //TF1 *Rew = new TF1("Rew","(0.00101563*x+0.817068)",50,150);
  //TF1 *Rew = new TF1("Rew","(3.10479e-06*x*x+0.000420232*x+0.841603 )",50,150);
  //TF1 *Rew = new TF1("Rew","(-5.65524e-05*x*x+0.00612228*x+0.67357+1.96711e-07*x*x*x )",50,150);
  
  //TF1 *Rew = new TF1("Rew","(6.27411e-06*x*x+(-0.000891689)*x+0.9341+2.35103e-08*x*x*x )",50,150);
  //TF1 *Rew = new TF1("Rew","(1.33974e-05*x*x+(-0.00157201)*x+0.954137)",50,150);
  //TF1 *Rew = new TF1("Rew","((0.000994608)*x+0.848401 )",50,150);
   
  //TF1 *Rew = new TF1("Rew","(( 0.00113736)*x+0.811508)",50,150);
  //TF1 *Rew = new TF1("Rew","((0.00223196)*x+0.766296+(-5.70096e-06)*x*x)",50,150);
  //TF1 *Rew = new TF1("Rew","((0.0121637)*x+0.47254+(-0.000109135)*x*x+(3.39503e-07)*x*x*x)",50,150);
  
   TF1 *Rew = new TF1("Rew","((0.00731618)*x+0.591077+(-6.85876e-05)*x*x+(2.20523e-07)*x*x*x)",50,150);
   //TF1 *Rew = new TF1("Rew","((0.000937837)*x+0.778931+(-1.7833e-06)*x*x)",50,150);
   //TF1 *Rew = new TF1("Rew","((0.000595762)*x+0.793029)",50,150);
  
  Double_t scale_bin=0;
  for(int i=0;i<h_dij_M->GetNbinsX();i++)
  {
   scale_bin= Rew->Eval(h_dij_M->GetBinCenter(i+1));
   h_dij_M->SetBinContent(i+1,scale_bin*h_dij_M->GetBinContent(i+1));
  }
  
  f_Rew_dij->cd();
  h_dij_M->Write();
  
  
  
}
