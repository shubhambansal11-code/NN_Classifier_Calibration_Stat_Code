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
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit ;

TH1* getTH1() ;
TTree* getTTree() ;
Double_t fitFunction(Double_t *x, Double_t *par);

void fitsForFTest_pT610_QCD()
{
  ////////////////////////////////////////////////
  // I m p o r t i n g   R O O T  T T r e e s   //
  ////////////////////////////////////////////////


  // I m p o r t   T T r e e   i n t o   a   R o o D a t a S e t
  // -----------------------------------------------------------

  //TTree* tree = getTTree() ;
  TString data_1517 = "../QCD_610.root";
  //TString data_1517 = "./Z+jet_YHE/800600-pTfull-nominal_MC16a.root";
  TFile File_Ttree(data_1517);
  //TTreeReader readtree("Ztag", &File);
  TTree *tree = (TTree*)File_Ttree.Get("tree_mass");
  float lj1_Xbb60_mass; 
  float weight_Xbb60;
  TFile *f_trial = new TFile("trial.root", "RECREATE");
  tree->SetBranchAddress("lj1_Xbb60_mass", &lj1_Xbb60_mass);
  tree->SetBranchAddress("weight_Xbb60", &weight_Xbb60); 
  TTree* tree_mass = new TTree("tree_mass","tree_mass") ;
  //Double_t* lj1_Xbb60 = new Double_t ;
  float lj1_Xbb60;
  float weight_60;
  tree_mass->Branch("lj1_Xbb60_M",&lj1_Xbb60) ;
  tree_mass->Branch("weight_Xbb60_W",&weight_60);
  
  int N = tree->GetEntries();
  for (int i = 0; i < N; ++i) {
    tree->GetEntry(i);
    
      lj1_Xbb60 = lj1_Xbb60_mass;
      weight_60 = weight_Xbb60;
      //RooDataSet dh_mass("dh_mass","dh_mass",RooArgSet(lj1_mass),Import(*tree)) ;
      tree_mass->Fill();
    
                              }

  // Define the desired observable that needs to be fitted
  RooRealVar lj1_Xbb60_M("lj1_Xbb60_M","lj1_Xbb60_M",50,150) ;
  //RooRealVar weight_Xbb60("weight_Xbb60","weight_Xbb60",-150,150) ;
  
  //RooRealVar  mass = new RooRealVar(caliConfig::_Category+"_mass",caliConfig::_Category+"_mass",mmin,mmax) ;
  auto * weight_Xbb60_W = new RooRealVar("weight_Xbb60_W","weight_Xbb60_W",-1.0e+08,1.0e+08) ;
  auto * dh_mass = new RooDataSet("dh_mass","dh_mass",tree_mass,RooArgSet(lj1_Xbb60_M,*weight_Xbb60_W),0,"weight_Xbb60_W") ;

  //RooRealVar lj1_mass("lj1_mass","lj1_mass",50,150) ;
  // Construct unbinned dataset importing tree branches x and y matching between branches and RooRealVars 
  // is done by name of the branch/RRV 
  // 
  // Note that ONLY entries for which x,y have values within their allowed ranges as defined in 
  // RooRealVar x and y are imported. Since the y values in the import tree are in the range [-15,15]
  // and RRV y defines a range [-10,10] this means that the RooDataSet below will have less entries than the TTree 'tree'
  
 
  //double xmin=0.55;
  //double ymin=0.65;
  //double ymax=0.85;
  //dh_mass.plotOn(frame3,Binning(100)) ;
  
  //RooDataSet dh_mass("dh_mass","dh_mass",Import(*tree_mass),RooArgSet(lj1_Xbb60_mass, weight_Xbb60),0,"weight_Xbb60") ;
    
  //lj1_Xbb60_mass.setRange(50,150);
  //lj1_Xbb60_mass.setRange("SU", 110, 150);
    
  //TString FuncType[10] = { "2nd_Poly", "3rd_Poly", "4th_Poly","5th_Poly","1st_exp", "2nd_exp", "3rd_exp", "4th_exp", "5th_exp","6th_exp"};
  //TString FuncType[4] = { "2nd_Poly", "3rd_Poly","1st_exp", "2nd_exp"};
  //TString FuncType[2] = { "3rd_Poly","5th_Poly"};
   //TString FuncType[6] = { "1st_exp", "2nd_exp", "3rd_exp", "4th_exp", "5th_exp","6th_exp"};
   TString FuncType[1] = { "6th_Poly"};
   //TString FuncType[9] = { "2nd_Poly", "3rd_Poly", "4th_Poly","5th_Poly","1st_exp", "2nd_exp", "3rd_exp", "4th_exp", "5th_exp"};
  //TString FuncType[6] = { "1st_exp", "2nd_exp", "3rd_exp", "4th_exp", "5th_exp","6th_exp"};
  //TString FuncType[4] = { "2nd_Poly", "3rd_Poly", "4th_Poly","5th_Poly"};
  //TString FuncType[2] = { "3rd_Poly","5th_Poly"};
  //TString FuncType[1] = {  "4th_exp"};
  //TString FuncType[1] = {"5th_Poly" };
  double npars = 0;
  RooGenericPdf gp;
  for(int i=0; i<1; i++)
  { 
    //TString frame = FuncType[i]
    RooPlot* frame = lj1_Xbb60_M.frame(Title("Unbinned data shown in default frame binning")) ;
    dh_mass->plotOn(frame,Binning(100)) ;
    dh_mass->Print();
    
    if(FuncType[i]=="2nd_Poly")
    {   
        npars = 3;
        
        /*RooRealVar c0("c0","c0",10,0,1e5);
        RooRealVar c1("c1","c1",10,-1e3,1e5);
        RooRealVar c2("c2","c2",1,-100,1e3);*/

        RooRealVar c0("c0","c0",100,0,1e5);
        RooRealVar c1("c1","c1",10,-100,1e3);
        RooRealVar c2("c2","c2",1,-100,1e3);

        /*RooRealVar c0("c0","c0",100,0,1e8);
        RooRealVar c1("c1","c1",-10,-1e5,1e7);
        RooRealVar c2("c2","c2",-10,-1e4,1e7);*/
        /*RooRealVar c0("c0","c0",10,0,1e5);
        RooRealVar c1("c1","c1",10,-1e3,1e5);
        RooRealVar c2("c2","c2",1,-100,1e3);*/
        RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_M) + c2*(lj1_Xbb60_M)*(lj1_Xbb60_M))", RooArgSet(c0,c1,c2,lj1_Xbb60_M)) ;
        gp.fitTo(*dh_mass);
        gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="3rd_Poly")
    {  
       npars = 4; 
      
       /*RooRealVar c0("c0","c0",100,0,1e8);
       RooRealVar c1("c1","c1",-10,-1e5,1e7);
       RooRealVar c2("c2","c2",-10,-1e4,1e7);
       RooRealVar c3("c3","c3",-1,-100,1e7);*/
       /*RooRealVar c0("c0","c0",100,-1e5,1e5);
       RooRealVar c1("c1","c1",-10,-1e5,1e5);
       RooRealVar c2("c2","c2",1,-100,1e3);
       RooRealVar c3("c3","c3",10,-100,1e3);*/

       /*RooRealVar c0("c0","c0",10,0,1e5);
       RooRealVar c1("c1","c1",10,-1e3,1e5);
       RooRealVar c2("c2","c2",1,-100,1e3);
       RooRealVar c3("c3","c3",10,-100,1e3);*/

       RooRealVar c0("c0","c0",100,0,1e5);
       RooRealVar c1("c1","c1",10,-10,1e3);
       RooRealVar c2("c2","c2",1,-100,1e3);
       RooRealVar c3("c3","c3",1,-100,1e3);

      /* RooRealVar c0("c0","c0",10,0,1e5);
       RooRealVar c1("c1","c1",10,-1e3,1e5);
       RooRealVar c2("c2","c2",1,-100,1e3);
       RooRealVar c3("c3","c3",10,-100,1e3);*/

       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_M) + c2*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c3*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M))", RooArgSet(c0,c1,c2,c3,lj1_Xbb60_M)) ;
       gp.fitTo(*dh_mass);
        gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="4th_Poly")
    {  
       npars = 5; 
       /*RooRealVar c0("c0","c0",100,-1e5,1e7);
       RooRealVar c1("c1","c1",-10,-1e3,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,-10,10);*/
       
       

       RooRealVar c0("c0","c0",100,0,1e8);
       RooRealVar c1("c1","c1",-10,-1e5,1e7);
       RooRealVar c2("c2","c2",-10,-1e4,1e7);
       RooRealVar c3("c3","c3",-1,-100,1e7);
       RooRealVar c4("c4","c4",1,-10,100);

       /*RooRealVar c0("c0","c0",10,0,1e5);
       RooRealVar c1("c1","c1",10,-1e3,1e5);
       RooRealVar c2("c2","c2",1,-100,1e3);
       RooRealVar c3("c3","c3",10,-100,1e3);
       RooRealVar c4("c4","c4",1,-10,10);*/

       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_M) + c2*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c3*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c4*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M))", RooArgSet(c0,c1,c2,c3,c4,lj1_Xbb60_M)) ;
       gp.fitTo(*dh_mass);
        gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="5th_Poly")
    {  
       npars = 6; 
       RooRealVar c0("c0","c0",100,-1e5,1e7);
       RooRealVar c1("c1","c1",-10,-1e3,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,-10,10);
       RooRealVar c5("c5","c5",1,-10,10);

       /*RooRealVar c0("c0","c0",100,0,1e7);
       RooRealVar c1("c1","c1",-10,-1e3,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,0,10);
       RooRealVar c5("c5","c5",1,0,10);*/
       
       /*RooRealVar c0("c0","c0",100,0,1e8);
       RooRealVar c1("c1","c1",-10,-1e5,1e7);
       RooRealVar c2("c2","c2",-10,-1e4,1e7);
       RooRealVar c3("c3","c3",-1,-100,1e7);
       RooRealVar c4("c4","c4",1,-10,100);
       RooRealVar c5("c5","c5",1,-10,100);*/

       /*RooRealVar c0("c0","c0",100,0,1e7);
       RooRealVar c1("c1","c1",-10,-1e3,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,0,10);
       RooRealVar c5("c5","c5",1,0,10);*/
       
       

       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_M) + c2*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c3*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c4*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c5*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M))", RooArgSet(c0,c1,c2,c3,c4,c5,lj1_Xbb60_M)) ;
       gp.fitTo(*dh_mass);
       gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="6th_Poly")
    {  
       npars = 7; 
       RooRealVar c0("c0","c0",100,-1e5,1e7);
       RooRealVar c1("c1","c1",-10,-1e3,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,-10,10);
       RooRealVar c5("c5","c5",1,-10,10);
       RooRealVar c6("c6","c6",-10,-100,1e3);

       /*RooRealVar c0("c0","c0",100,0,1e7);
       RooRealVar c1("c1","c1",-10,-1e3,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,0,10);
       RooRealVar c5("c5","c5",1,0,10);*/
       
       /*RooRealVar c0("c0","c0",100,0,1e8);
       RooRealVar c1("c1","c1",-10,-1e5,1e7);
       RooRealVar c2("c2","c2",-10,-1e4,1e7);
       RooRealVar c3("c3","c3",-1,-100,1e7);
       RooRealVar c4("c4","c4",1,-10,100);
       RooRealVar c5("c5","c5",1,-10,100);*/

       /*RooRealVar c0("c0","c0",100,0,1e7);
       RooRealVar c1("c1","c1",-10,-1e3,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,0,10);
       RooRealVar c5("c5","c5",1,0,10);*/
       
       

       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_M) + c2*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c3*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c4*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c5*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M) + c5*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M)*(lj1_Xbb60_M))", RooArgSet(c0,c1,c2,c3,c4,c5,c6,lj1_Xbb60_M)) ;
       gp.fitTo(*dh_mass);
       gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    /*else if(FuncType[i]=="5th_Poly")
    {  
       npars = 6; 
       RooRealVar c0("c0","c0",100,0,1e7);
       RooRealVar c1("c1","c1",-10,-1e4,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",0.1,-1,10);
       RooRealVar c5("c5","c5",0.1,-1,10);
       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_mass) + c2*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c3*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c4*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c5*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass))", RooArgSet(c0,c1,c2,c3,c4,c5,lj1_Xbb60_mass)) ;
       gp.fitTo(dh_mass, Range("SL,SU"));
       gp.plotOn(frame3, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
       gp.plotOn(frame3, Range(50,150),LineStyle(kDashed));
       gp.paramOn(frame3);
    }*/

    else if(FuncType[i]=="1st_exp")
    {  
      npars = 1;   
      //RooRealVar c0("c0","c0",200,-1000,1000);
      //RooRealVar c0("c0","c0",-0.553,-1,1);
      //RooRealVar c0("c0","c0",20,-100,100);
      RooRealVar c0("c0","c0",-0.553,-10,10);
      //RooRealVar c0("c0","c0",20,-100,100);
      RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_M/100))", RooArgSet(c0,lj1_Xbb60_M)) ;
      gp.fitTo(*dh_mass);
      gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
      }

    else if(FuncType[i]=="2nd_exp")
    { 
      npars = 2;   
      /*RooRealVar c0("c0","c0",20,-100,100);
      RooRealVar c1("c1","c1",-20,-200,200);*/
      /*RooRealVar c0("c0","c0",10,-100,1e5);
      RooRealVar c1("c1","c1",-10,-100,1e5);*/
      /*RooRealVar c0("c0","c0",-0.54,-10,10);
      RooRealVar c1("c1","c1",0.55,-10,10);*/

      /*RooRealVar c0("c0","c0",-0.54,-1,1);
      RooRealVar c1("c1","c1",-0.4,-1,1);*/

     /* RooRealVar c0("c0","c0",-0.54,-10,10);
      RooRealVar c1("c1","c1",0.55,-10,10);*/

      RooRealVar c0("c0","c0",-0.553,-10,10);
      RooRealVar c1("c1","c1",0.55,-10,10);
      //RooRealVar c0("c0","c0",-20,-100,100);
      //RooRealVar c1("c1","c1",20,-100,100);
      /*RooRealVar c0("c0","c0",-0.54,-1,1);
      RooRealVar c1("c1","c1",-0.4,-1,1);*/
        /*RooRealVar c0("c0","c0",1000,-1e7,1e7);
      RooRealVar c1("c1","c1",-2000,-2e7,2e7);*/
        //RooRealVar c4("c4","c4",0.3,-1,1);
      RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_M/100)+c1*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100))", RooArgSet(c0,c1,lj1_Xbb60_M)) ;
      gp.fitTo(*dh_mass);
       gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
      }  

    else if(FuncType[i]=="3rd_exp")
    {   
        npars = 3;
        /*RooRealVar c0("c0","c0",20,-100,100);
        RooRealVar c1("c1","c1",-20,-200,200);
        RooRealVar c2("c2","c2",70,-200,200);*/

        
        /*RooRealVar c0("c0","c0",0.33,0,1);
        RooRealVar c1("c1","c1",-0.40,-1,1);
        RooRealVar c2("c2","c2",0.40,0,1);*/

        RooRealVar c0("c0","c0",-0.553,-10,10);
        RooRealVar c1("c1","c1",0.55,-10,10);
        RooRealVar c2("c2","c2",-0.54,-10,10);
        //RooRealVar c3("c3","c3",-0.4,-1,1);
        //RooRealVar c4("c4","c4",0.3,-1,1);
        RooGenericPdf gp("gp","gp","exp(c0*(lj1_Xbb60_M/100) + c1*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100) + c2*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100))", RooArgSet(c0,c1,c2,lj1_Xbb60_M)) ;
        gp.fitTo(*dh_mass);
        gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="4th_exp")
    {   
        npars = 4;
        
        /*RooRealVar c0("c0","c0",0.54,0,1);
        RooRealVar c1("c1","c1",0.33,0,1);
        RooRealVar c2("c2","c2",-0.54,-1,1);
        RooRealVar c3("c3","c3",-0.4,-1,1);*/

        /*RooRealVar c0("c0","c0",20,-100,100);
        RooRealVar c1("c1","c1",-20,-200,200);
        RooRealVar c2("c2","c2",70,-200,200);
        RooRealVar c3("c3","c3",70,-200,200);*/


        RooRealVar c0("c0","c0",-0.553,-10,10);
        RooRealVar c1("c1","c1",0.55,-10,10);
        RooRealVar c2("c2","c2",-0.54,-10,10);
        RooRealVar c3("c3","c3",-0.4,-10,10);
        //RooRealVar c4("c4","c4",0.3,-1,1);
        RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_M/100) + c1*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100) + c2*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)+c3*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100))", RooArgSet(c0,c1,c2,c3,lj1_Xbb60_M)) ;
        gp.fitTo(*dh_mass);
       gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);   
    }

    else if(FuncType[i]=="5th_exp")
    {   
        npars = 5;
    
        /*RooRealVar c0("c0","c0",-0.4,-10,10);
        RooRealVar c1("c1","c1",0.3,-10,10);
        RooRealVar c2("c2","c2",-0.54,-10,10);
        RooRealVar c3("c3","c3",-0.4,-10,10);
        RooRealVar c4("c4","c4",0.3,-10,10);*/

       /* RooRealVar c0("c0","c0",20,-100,100);
        RooRealVar c1("c1","c1",-20,-200,200);
        RooRealVar c2("c2","c2",70,-200,200);
        RooRealVar c3("c3","c3",70,-200,200);
        RooRealVar c4("c4","c4",-70,-200,200);*/

        RooRealVar c0("c0","c0",-0.553,-10,10);
        RooRealVar c1("c1","c1",0.55,-10,10);
        RooRealVar c2("c2","c2",-0.54,-10,10);
        RooRealVar c3("c3","c3",-0.4,-10,10);
        RooRealVar c4("c4","c4",0.3,-10,10);
       // RooRealVar c4("c5","c",-70,-200,200);
        RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_M/100) + c1*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100) + c2*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)+c3*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)+c4*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100))", RooArgSet(c0,c1,c2,c3,c4,lj1_Xbb60_M)) ;
        gp.fitTo(*dh_mass);
       gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="6th_exp")
    {   
        npars = 6;
        /*RooRealVar c0("c0","c0",20,-100,100);
        RooRealVar c1("c1","c1",-20,-200,200);
        RooRealVar c2("c2","c2",70,-200,200);
        RooRealVar c3("c3","c3",70,-200,200);
        RooRealVar c4("c4","c4",-70,-200,200);
        RooRealVar c5("c5","c5",70,-100,100);*/

       /* RooRealVar c0("c0","c0",-0.553,-10,10);
        RooRealVar c1("c1","c1",0.55,-10,10);
        RooRealVar c2("c2","c2",-0.54,-10,10);
        RooRealVar c3("c3","c3",-0.4,-10,10);
        RooRealVar c4("c4","c4",0.3,-10,10);
        RooRealVar c5("c5","c5",0.3,-10,10);*/

        RooRealVar c0("c0","c0",-0.553,-10,10);
        RooRealVar c1("c1","c1",0.55,-10,10);
        RooRealVar c2("c2","c2",-0.54,-10,10);
        RooRealVar c3("c3","c3",-0.4,-10,10);
        RooRealVar c4("c4","c4",0.3,-10,10);
        RooRealVar c5("c5","c5",0.3,-10,10);
        
        

        /*RooRealVar c0("c0","c0",-0.553,-10,10);
        RooRealVar c1("c1","c1",0.55,-10,10);
        RooRealVar c2("c2","c2",-0.54,-10,10);
        RooRealVar c3("c3","c3",-0.4,-10,10);
        RooRealVar c4("c4","c4",0.3,-10,10);*/
       // RooRealVar c4("c5","c",-70,-200,200);
        RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_M/100) + c1*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100) + c2*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)+c3*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)+c4*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)+c5*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100)*(lj1_Xbb60_M/100))", RooArgSet(c0,c1,c2,c3,c4,c5,lj1_Xbb60_M)) ;
        gp.fitTo(*dh_mass);
       gp.plotOn(frame, Range(50,150),LineStyle(kBlue));
        //gp.paramOn(frame);
    }

    


    /*
    if(FuncType[i]=="2nd_Poly")
    {   
        npars = 3;
        RooRealVar c0("c0","c0",100,0,1e5);
        RooRealVar c1("c1","c1",10,-10,1e3);
        RooRealVar c2("c2","c2",1,-100,1e3);
        RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_mass) + c2*(lj1_Xbb60_mass)*(lj1_Xbb60_mass))", RooArgSet(c0,c1,c2,lj1_Xbb60_mass)) ;
        gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="3rd_Poly")
    {  
       npars = 4; 
       RooRealVar c0("c0","c0",100,-10,1e5);
       RooRealVar c1("c1","c1",10,-10,1e5);
       RooRealVar c2("c2","c2",-1,-100,1e3);
       RooRealVar c3("c3","c3",10,-1,1e3);
       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_mass) + c2*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c3*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass))", RooArgSet(c0,c1,c2,c3,lj1_Xbb60_mass)) ;
       gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="4th_Poly")
    {  
       npars = 5; 
       RooRealVar c0("c0","c0",100,0,1e7);
       RooRealVar c1("c1","c1",-10,-1e4,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,-10,10);
       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_mass) + c2*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c3*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c4*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass))", RooArgSet(c0,c1,c2,c3,c4,lj1_Xbb60_mass)) ;
       gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="5th_Poly")
    {  
       npars = 6; 
       RooRealVar c0("c0","c0",100,0,1e7);
       RooRealVar c1("c1","c1",-10,-1e4,1e5);
       RooRealVar c2("c2","c2",-10,-1e4,1e5);
       RooRealVar c3("c3","c3",-1,-100,1e5);
       RooRealVar c4("c4","c4",1,0,10);
       RooRealVar c5("c5","c5",1,-10,10);
       RooGenericPdf gp("gp","GenericPDF","(c0 + c1*(lj1_Xbb60_mass) + c2*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c3*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c4*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass) + c5*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass)*(lj1_Xbb60_mass))", RooArgSet(c0,c1,c2,c3,c4,c5,lj1_Xbb60_mass)) ;
       gp.fitTo(dh_mass, Range("SL,SU"));
       gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
       gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
    }

    

    else if(FuncType[i]=="1st_exp")
    {  
      npars = 1;   
      RooRealVar c0("c0","c0",200,-1000,1000);
      RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_mass/100))", RooArgSet(c0,lj1_Xbb60_mass)) ;
      gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
      }

    else if(FuncType[i]=="2nd_exp")
    { 
      npars = 2;   
      RooRealVar c0("c0","c0",20,-100,100);
      RooRealVar c1("c1","c1",-20,-200,200);
      RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_mass/100)+c1*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100))", RooArgSet(c0,c1,lj1_Xbb60_mass)) ;
      gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
      }  

    else if(FuncType[i]=="3rd_exp")
    {   
        npars = 3;
        RooRealVar c0("c0","c0",20,-100,100);
        RooRealVar c1("c1","c1",-20,-200,200);
        RooRealVar c2("c2","c2",70,-200,200);
        RooGenericPdf gp("gp","gp","exp(c0*(lj1_Xbb60_mass/100) + c1*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100) + c2*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100))", RooArgSet(c0,c1,c2,lj1_Xbb60_mass)) ;
        gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
    }

    else if(FuncType[i]=="4th_exp")
    {   
        npars = 4;
        RooRealVar c0("c0","c0",20,0,100);
        RooRealVar c1("c1","c1",0,-200,200);
        RooRealVar c2("c2","c2",70,-200,200);
        RooRealVar c3("c3","c3",70,-200,200);
        RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_mass/100) + c1*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100) + c2*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)+c3*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100))", RooArgSet(c0,c1,c2,c3,lj1_Xbb60_mass)) ;
        gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);   
    }

    else if(FuncType[i]=="5th_exp")
    {   
        npars = 5;
        RooRealVar c0("c0","c0",20,-100,100);
        RooRealVar c1("c1","c1",-20,-200,200);
        RooRealVar c2("c2","c2",70,-200,200);
        RooRealVar c3("c3","c3",70,-200,200);
        RooRealVar c4("c4","c4",-70,-200,200);
        RooGenericPdf gp("gp","GenericPDF","exp(c0*(lj1_Xbb60_mass/100) + c1*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100) + c2*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)+c3*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)+c4*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100)*(lj1_Xbb60_mass/100))", RooArgSet(c0,c1,c2,c3,c4,lj1_Xbb60_mass)) ;
        gp.fitTo(dh_mass, Range("SL,SU"));
        gp.plotOn(frame, Range("SL,SU"), NormRange("SL,SU"), LineColor(kBlue)) ;
        gp.plotOn(frame, Range(50,150),LineStyle(kDashed));
        //gp.paramOn(frame);
    }
    */
     
     //lj1_mass.setRange("SL", 50, 80);
     //lj1_mass.setRange("SU", 110, 150);
     

     
  //RooPlot* frame3 = lj1_mass.frame(Title("Unbinned data shown in default frame binning")) ;
  
        
    Double_t chi2 = frame->chiSquare(npars);
    //Double_t chi2 = frame->chiSquare("model","data");
    //cout << " CHI2 = " << chi2 << endl;
    cout << " CHI2 = " << chi2<< endl;
    //Double_t red_chi2 = chi2/npars;
    double prob = ROOT::Math::chisquared_cdf_c(chi2*(npars), npars);
    cout << " CHI2 p-value = " << prob << endl;
    TString text_Chi2=Form("#chi2 = %0.3f",chi2*(100-npars));
    TString text_Chi2ndf=Form("#chi2/ndf = %0.3f",chi2);
    TString text_Chi2prob=Form("Percen. Prob = %0.3f",prob*100);
    TCanvas* c = new TCanvas("dataimport_check_TTree_Thist","dataimport_check_TTree_Thist",800,800) ;
    frame->SetTitle("");
   frame->SetXTitle("Mass/ GeV");
   frame->GetYaxis()->SetTitleOffset(1.5);
   frame->SetYTitle("Events");
   frame->Draw();
    double lx = 0.15;
   double ly = 0.25;
   TLatex *tex1 = new TLatex();
    if(FuncType[i]=="2nd_Poly")
    {
      tex1= new TLatex(lx,ly,"FF: #Sigma_{i=0 to 2} a_{i}*x^{i})");
    }
    else if(FuncType[i]=="3rd_Poly")
    {
      tex1= new TLatex(lx,ly,"FF: #Sigma_{i=0 to 3} a_{i}*x^{i})");
    }
    else if(FuncType[i]=="4th_Poly")
    {
      tex1= new TLatex(lx,ly,"FF: #Sigma_{i=0 to 4} a_{i}*x^{i})");
    }
    else if(FuncType[i]=="5th_Poly")
    {
      tex1= new TLatex(lx,ly,"FF: #Sigma_{i=0 to 5} a_{i}*x^{i})");
    }
    else if(FuncType[i]=="6th_Poly")
    {
      tex1= new TLatex(lx,ly,"FF: #Sigma_{i=0 to 6} a_{i}*x^{i})");
    }
    else if(FuncType[i]=="1st_exp")
    {
      tex1= new TLatex(lx,ly,"FF: exp(a_{1}*x/100)");
    }
    else if(FuncType[i]=="2nd_exp")
    {
      tex1= new TLatex(lx,ly,"dijet FF: exp(#Sigma_{i=1 to 2} a_{i}*(x/100)^{i})");
    }
    else if(FuncType[i]=="3rd_exp")
    {
      tex1= new TLatex(lx,ly,"dijet FF: exp(#Sigma_{i=1 to 3} a_{i}*(x/100)^{i})");
    }
    else if(FuncType[i]=="4th_exp")
    {
      tex1= new TLatex(lx,ly,"dijet FF: exp(#Sigma_{i=1 to 4} a_{i}*(x/100)^{i})");
    }
    else if(FuncType[i]=="5th_exp")
    {
      tex1= new TLatex(lx,ly,"dijet FF: exp(#Sigma_{i=1 to 5} a_{i}*(x/100)^{i})");
    }
    else if(FuncType[i]=="6th_exp")
    {
      tex1= new TLatex(lx,ly,"dijet FF: exp(#Sigma_{i=1 to 6} a_{i}*(x/100)^{i})");
    }
    tex1->SetNDC();
   tex1->SetTextSize(0.04);
   tex1->SetTextColor(1);
   tex1->SetTextFont(42);
   tex1->Draw();
    TLatex * Tag = new TLatex(0.15,0.45,"Xbb 60% Eff");
  Tag->SetNDC();
  Tag->SetTextFont(42);
  Tag->SetTextColor(1);
  Tag->SetTextSize(0.04);
  Tag->Draw();

  TLatex * pT = new TLatex(0.15,0.41,"pT=[600,1000] GeV");
  pT->SetNDC();
  pT->SetTextFont(42);
  pT->SetTextColor(1);
  pT->SetTextSize(0.04);
  pT->Draw();

  TLatex * chi2R = new TLatex(0.15,0.36,text_Chi2.Data());
  chi2R->SetNDC();
  chi2R->SetTextFont(42);
  chi2R->SetTextColor(1);
  chi2R->SetTextSize(0.04);
  chi2R->Draw();

  TLatex * chi2nd = new TLatex(0.15,0.32,text_Chi2ndf.Data());
  chi2nd->SetNDC();
  chi2nd->SetTextFont(42);
  chi2nd->SetTextColor(1);
  chi2nd->SetTextSize(0.04);
  chi2nd->Draw();

  TLatex * probi = new TLatex(0.15,0.28,text_Chi2prob.Data());
  probi->SetNDC();
  probi->SetTextFont(42);
  probi->SetTextColor(1);
  probi->SetTextSize(0.04);
  probi->Draw();
   gPad->SaveAs("./FTestFits/"+FuncType[i]+"UnBinned_lj1mass_Xbb60_pT610_SBRangeChange_QCD.pdf");
   //gPad->delete();
  }
  
}
