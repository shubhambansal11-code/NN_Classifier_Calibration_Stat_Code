#include "TFile.h"
#include "TTree.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TLegend.h"
#include "TColor.h"
#include "TDatime.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLine.h"
#include "THStack.h"
#include "TStyle.h"
#include "math.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TVector.h"
#include "TMatrix.h"
#include "TMatrixFUtils.h"
#include <TMatrixDSymEigen.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "string.h"
   
   Double_t fitFunction(Double_t *x, Double_t *par) {
      //return background(x,par) + lorentzianPeak(x,&par[3]);
       //return TMath::Exp(par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0]);
       //return TMath::Exp(par[0] + par[1]*x[0] + par[2]*x[0]*x[0]);
       //return TMath::Exp(par[0]*x[0]/100);
       return par[0]+par[1]*x[0];
       //return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
       //return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
       //return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
       //return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
   }

   void fittodatavMC() {
   TString base = "/eos/user/s/shbansal/scripts/";
   TString input_fit[1]= {"DatavMC_SB_pT610.root"};
   for(int i=0; i<1; i++ )
   {
   TFile File(base+input_fit[i]);
   //TFile File(base_hist);
   TH1D *hbkg=(TH1D*)File.Get("Zcand_Xbb60_mass_pT610_SB");
   TFile *file_fit = new TFile("./OptFitOutput/1OPolyFitto"+input_fit[i], "RECREATE");
   
   // create a TF1 with the range from 0 to 3 and 6 parameters
   //fitFunction->SetParameter(0,1);
   TF1 *fitFcn = new TF1("fitFcn",fitFunction,50,150,2); //keep changing the parameters hereof 
   hbkg->Fit("fitFcn");
   //int fitStatus = hbkg->Fit("fitFcn", Range("SL,SU"));
   TF1 *fit = hbkg->GetFunction("fitFcn");
   TFitResultPtr r = hbkg->Fit("fitFcn","S"); // TFitResultPtr contains the TFitResult
   TMatrixDSym cov = r->GetCovarianceMatrix();
   Double_t chi2 = r->Chi2();
   Double_t ndf = 0;
   
   ndf= 58;  //100 bins -2 degree of freedom
   TString text_Chi2=Form("#chi2/ndf = %0.3f",chi2/ndf);
   //Double_t redchi2 = chi2/NDf;
   //string chi2 = to_string(chi2);
   Double_t prob = r->Prob();
   Double_t par0 = r->Parameter(0);
   Double_t err0 = r->ParError(0);
   Double_t par1 = r->Parameter(1);
   Double_t err1 = r->ParError(1);
   std::cout << "Probability" << prob << std::endl;
   TString text_Chi2prob=Form("Percen. Prob = %0.3f",prob*100);
   //std::cout<<"All plots have been read. Beginning to draw plots."<<std::endl;
   r->Print("V");
   file_fit->cd();
   fit->Write();
   hbkg->Write();
   r->Write();
   

   std::shared_ptr<TCanvas> canv(new TCanvas("canv","",800,800));
    canv.get()->cd(); 
    std::shared_ptr<TPad> pad1(new TPad("pad1","pad1",0.007,0.,1,1));
    pad1.get()->SetLeftMargin(0.17);
    pad1.get()->SetBottomMargin(0.15);
    pad1.get()->Draw();
   

   TFile File_r("./OptFitOutput/1OPolyFitto"+input_fit[i]);
   TH1D *hbkg_fit=(TH1D*)File_r.Get("Zcand_Xbb60_mass_pT610_SB");
   TF1 *bkg_fit = (TF1*)File_r.Get("fitFcn");
   pad1->cd();
   hbkg_fit->SetStats(0);
   hbkg_fit->SetTitle("");
   hbkg_fit->GetXaxis()->SetTitle("Zcand mass [GeV]");
  
   hbkg_fit->GetYaxis()->SetTitle("Data/MC");
   hbkg_fit->GetXaxis()->SetLabelSize(.04);
   hbkg_fit->GetYaxis()->SetLabelSize(.04);
   bkg_fit->SetLineWidth(3);
   bkg_fit->SetLineColor(kBlue);
   hbkg_fit->Draw();
   bkg_fit->Draw("SAME");
   
   TLatex *tex00 = new TLatex();
  float lx ; float ly;
  lx=0.41;
  ly=0.63;
  tex00= new TLatex(lx,ly,"#bf{#it{ATLAS} Internal}");
  tex00->SetNDC();
  tex00->SetTextSize(0.025);
  tex00->SetTextColor(1);
  tex00->SetTextFont(42);

   TLatex *tex0 = new TLatex();
  //float lx ; float ly;
  lx=0.41;
  ly=0.59;
  tex0= new TLatex(lx,ly,"Xbb 60%, 600 <= p_{T} < 1000 GeV");
  tex0->SetNDC();
  tex0->SetTextSize(0.025);
  tex0->SetTextColor(1);
  tex0->SetTextFont(42);

   TLatex *tex1 = new TLatex();
  //float lx ; float ly;
  lx=0.41;
  ly=0.55;
  //tex1= new TLatex(lx,ly,"dijet FF: N_{bkg}#Sigma_{i=0 to 3} a_{i}*X^{i}, X=(x-50)/100");
  tex1= new TLatex(lx,ly,"Function: #Sigma_{i=0 to 1} a_{i}*x^{i}");
  tex1->SetNDC();
  tex1->SetTextSize(0.025);
  tex1->SetTextColor(1);
  tex1->SetTextFont(42);

  TLatex *tex2= new TLatex();
  lx=0.41;
  ly=0.51;
  tex2= new TLatex(lx,ly,text_Chi2.Data()); 
  tex2->SetNDC();
  tex2->SetTextSize(0.025);
  tex2->SetTextColor(1);
  tex2->SetTextFont(42);

  TLatex *tex3= new TLatex();
  lx=0.41;
  ly=0.49;
  tex3= new TLatex(lx,ly,text_Chi2prob.Data()); 
  //tex3= new TLatex(lx,ly,""); 
  tex3->SetNDC();
  tex3->SetTextSize(0.025);
  tex3->SetTextColor(1);
  tex3->SetTextFont(42);

  /*TLatex *tex2= new TLatex();
  lx=0.35;
  ly=0.75;
  //tex2= new TLatex(lx,ly,"");
  tex2= new TLatex(lx,ly,"");
  tex2->SetNDC();
  tex2->SetTextSize(0.025*1.5);
  tex2->SetTextColor(1);
  tex2->SetTextFont(42);*/

  //pad1->cd();
  tex00->Draw("same");
  tex0->Draw("same");
  tex1->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  //tex4->Draw("same");

  //leg->Draw();
  TString name;
  name= TString("./OptFitPlots/1OrderPolyFittoDvMC"+input_fit[i]+"_trial.pdf");
  canv->SaveAs(name.Data());

   }
   }
