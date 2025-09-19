#include "RooDataSet.h"
#include "RooHistFunc.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooChebychev.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooAbsArg.h"
#include "RooWorkspace.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
#include <string>

using namespace RooFit;
using namespace RooStats;

void makePlot()

{
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit/quickFit/Vbbfull_Untag_ZTag.root");
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit/quickFit/workspace/hbbj/Vbbfull_Untag_ZTag/Vbbfull_Untag_ZTag_model_test.root");

  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit/quickFit/workspace/hbbj/Vbbfull_2b60_ZTag/Vbbfull_2b60_ZTag_model_test.root");
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit/quickFit/workspace/hbbj/Vbbfull/Vbbfull_model_test.root");
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit/quickFit/workspace/hbbj/Vbbfull_4thexp/Vbbfull_4thexp_model_test.root");
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit_Xbb/quickFit/output_Vbbfull_Xbb50_139.root");
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit_Xbb/quickFit/output_Vbbfull_Xbb60_ttcat2_ttconstr.root");
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit_Xbb/quickFit/output_Vbbfull_Xbb60_TrigCut.root");
  //TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit_Xbb/quickFit/output_Vbbfull_Xbb60_trigcut_dscb.root");
  TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit_Xbb/quickFit/output_Vbbfull_Xbb60_48_wotrigCut.root");
//TFile *file = new TFile("/afs/cern.ch/work/s/shbansal/QT_2020/qt-2019-2020/qt-2019-2020/fit_Xbb/quickFit/output_Vbbfull_Xbb60_3gaus.root");

RooWorkspace *w = (RooWorkspace*)file->Get("combWS");
RooDataSet* data = (RooDataSet*)w->data("combData");
RooAddPdf* model = (RooAddPdf*)w->pdf("_modelSB_Vbbfullchannel");
RooRealVar *invm = (RooRealVar*)w->var("obs_x_channel");

  // prepare frame for plotting
  RooPlot* invmframe = invm->frame();
  TCanvas * can = new TCanvas("can ","can ",800,600);
  // plot data on frame
  data->plotOn(invmframe, Name("data"));
  model->plotOn(invmframe,Name("Combined model"),LineStyle(ELineStyle::kSolid),LineColor(kBlue),LineWidth(2));
  model->plotOn(invmframe,Components("pdf__Zboson_Vbbfullchannel"),Name("Z"),LineStyle(ELineStyle::kSolid),LineColor(kBlue-2),LineWidth(2));
  //model->plotOn(invmframe,Components("pdf__Wboson_Vbbfullchannel"),Name("W"),LineStyle(ELineStyle::kSolid),LineColor(kYellow),LineWidth(2));
  model->plotOn(invmframe,Components("pdf__ttbar1_Vbbfullchannel"),Name("tt Cat1"),LineStyle(ELineStyle::kSolid),LineColor(kYellow),LineWidth(2));
  //model->plotOn(invmframe,Components("pdf__ttbar_Vbbfullchannel"),Name("tt"),LineStyle(ELineStyle::kSolid),LineColor(kYellow),LineWidth(2));
  model->plotOn(invmframe,Components("pdf__ttbar2_Vbbfullchannel"),Name("tt Cat2"),LineStyle(ELineStyle::kSolid),LineColor(kGreen),LineWidth(2));
  model->plotOn(invmframe,Components("pdf__dijets_Vbbfullchannel"),Name("dijet"),LineStyle(ELineStyle::kSolid),LineColor(kRed),LineWidth(2));
  
  // draw frame
  invmframe->SetTitle("");
  invmframe->SetXTitle("Mass/ GeV");
  invmframe->SetYTitle("Events/ 2 GeV");
  invmframe->Draw();
  
  TLegend *leg1 = new TLegend(0.50,0.65,0.65,0.85);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.03);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.3);
  //leg1->SetFillColor(kWhite);
  //leg1->SetLineColor(kWhite);
  leg1->AddEntry("data","Data", "P");
  leg1->AddEntry("Combined model","Combined (S+B)", "L");
  leg1->AddEntry("Z","Z", "L");
  //leg1->AddEntry("W","W","L");
  //leg1->AddEntry("tt","tt", "L");
  //leg1->AddEntry("tt Cat2","tt Cat2", "L");
  leg1->AddEntry("tt Cat1","tt Cat1", "L");
  leg1->AddEntry("tt Cat2","tt Cat2", "L");
  leg1->AddEntry("dijet","dijets", "L");
  leg1->Draw();
  
  TLatex * AI = new TLatex(0.17,0.83, "#bf{#it{ATLAS} Internal}"); 
  AI->SetNDC();
  AI->SetTextFont(42);
  AI->SetTextColor(1);
  AI->SetTextSize(0.04);
  AI->Draw();

  TLatex * ra = new TLatex(0.15,0.50, "#sqrt{s} = 13 TeV, 80.4 fb^{-1}"); 
  ra->SetNDC();
  ra->SetTextFont(42);
  ra->SetTextColor(1);
  ra->SetTextSize(0.04);
  ra->Draw();

  TLatex * Tag = new TLatex(0.15,0.45,"Xbb 60% Eff, pT = [450,1000] GeV");
  Tag->SetNDC();
  Tag->SetTextFont(42);
  Tag->SetTextColor(1);
  Tag->SetTextSize(0.04);
  Tag->Draw();

  TLatex * mu = new TLatex(0.15,0.41,"#mu_{Z}=1.112 #pm 0.102");
  mu->SetNDC();
  mu->SetTextFont(42);
  mu->SetTextColor(1);
  mu->SetTextSize(0.04);
  mu->Draw();

  TLatex * chi2 = new TLatex(0.15,0.36,"#chi^{2}/ndf =1.816 , #chi^{2} Prob: 12.26% ");
  chi2->SetNDC();
  chi2->SetTextFont(42);
  chi2->SetTextColor(1);
  chi2->SetTextSize(0.04);
  chi2->Draw();
  //*save frame for further investigation*//
  gPad->SetTicks();
  gPad->SaveAs("XbbPlots_1611/Plot_Postfit_Xbb60_pTinc450_tt2split_2GeV_48bin.pdf");
  w->Print();
}
