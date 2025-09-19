# ==========================================
# Systematic Smoothing for XML Model Making
#
# Original author:Migle Stankaityte
# 4 February 2020
# Refurbished for QT PhD 2019-2020
# ==========================================

#import
import sys, re

import sigfit_init as init

import ROOT
from ROOT import gROOT, gStyle, TH1F, TF1, TCanvas, TFile, TAxis
from ROOT import TSpline3
from array import array

import numpy as np

sys.path.insert(0, 'plotting/')
import plot_helpers

sys.path.insert(0, 'python/')
import helpers as other_helpers 

gROOT.SetBatch(True)


# ==========================================
# TODO add wider plots!!

def getSmoothedSysts(template,rebin,outpath,nbins,mmin,mmax,n_init,m_fitmin,m_fitmax,reg):

  #------------ coarse plots -------------
  # nominal
  hname = template['hist'].GetTitle()

  # variable binning 
  if template['smooth_th']:
    template['h_coarse_fixed'] = template['whist'].Clone(hname+'_coarse_fixed')
    template['h_coarse_fixed'].SetTitle(hname+'_coarse_fixed') 

    th = template['smooth_th'] # thershold
    template['h_coarse'], edges = setVarBins(template['h_coarse_fixed'],hname+'_coarse',n_init,th)
    nbins_ratio = len(edges)-1

  # fixed binning
  else :
    template['h_coarse'] = template['whist'].Clone(hname+'_coarse')
    template['h_coarse'].SetTitle(hname+'_coarse') 
    nbins_ratio = n_init

  # print
  print "Smoothing "+hname+" systematics using spline fits..."

  # systematics 
  template['hsyst_up_coarse'] = {}
  template['hsyst_down_coarse'] = {}
  if template['smooth_th']:
    template['hsyst_up_coarse_fixed'] = {}
    template['hsyst_down_coarse_fixed'] = {}

  template['hsyst_up_orig'] = {}
  template['hsyst_down_orig'] = {}

  # syst/nom ratios
  hup_ratio = {}
  hdown_ratio = {}

  # splines
  sup_ratio = {}
  sdown_ratio = {}
  
  # smoothed / fine bin histos
  hup_rsmooth = {}
  hdown_rsmooth = {}
  
  hup_smooth = {}
  hdown_smooth = {}

  print "\tRebinning, making ratios and splines..."
  for syst in template['whsyst_up'] :
    
    hup_name = template['hsyst_up'][syst].GetTitle()
    hdown_name = template['hsyst_down'][syst].GetTitle()

    # var-bin
    if template['smooth_th']:
      template['hsyst_up_coarse_fixed'][syst] = template['whsyst_up'][syst].Clone(hup_name+'_coarse_fixed')
      template['hsyst_down_coarse_fixed'][syst] = template['whsyst_down'][syst].Clone(hdown_name+'_coarse_fixed')
    
      template['hsyst_up_coarse_fixed'][syst].SetTitle(hup_name+'_coarse_fixed') 
      template['hsyst_down_coarse_fixed'][syst].SetTitle(hdown_name+'_coarse_fixed') 
    
      template['hsyst_up_coarse'][syst] = rebinToEdges(template['hsyst_up_coarse_fixed'][syst],hup_name+'_coarse',n_init,edges)
      template['hsyst_down_coarse'][syst] = rebinToEdges(template['hsyst_down_coarse_fixed'][syst],hdown_name+'_coarse',n_init,edges)

    # fixed-bin
    else : 
      template['hsyst_up_coarse'][syst] = template['whsyst_up'][syst].Clone(hup_name+'_coarse')
      template['hsyst_down_coarse'][syst] = template['whsyst_down'][syst].Clone(hdown_name+'_coarse')
    
      template['hsyst_up_coarse'][syst].SetTitle(hup_name+'_coarse') 
      template['hsyst_down_coarse'][syst].SetTitle(hdown_name+'_coarse') 
    
    # ratios
    hup_ratio[syst] = template['hsyst_up_coarse'][syst].Clone(hup_name+'_ratio')
    hdown_ratio[syst] = template['hsyst_down_coarse'][syst].Clone(hdown_name+'ratio')
    
    hup_ratio[syst].SetTitle(hup_name+'_ratio') 
    hdown_ratio[syst].SetTitle(hdown_name+'_ratio') 

    hup_ratio[syst].Divide(template['h_coarse'])
    hdown_ratio[syst].Divide(template['h_coarse'])

    # splines
    sup_ratio[syst], maxerr_up = makeSpline(hup_ratio[syst],nbins_ratio)
    sdown_ratio[syst], maxerr_down = makeSpline(hdown_ratio[syst],nbins_ratio)
    
    # initialising histos for smoothed ratios and temps
    hup_rsmooth[syst] = TH1F(hup_name+'_rsmooth',hup_name+'_rsmooth',nbins,mmin,mmax)
    hdown_rsmooth[syst] = TH1F(hdown_name+'_rsmooth',hdown_name+'_rsmooth',nbins,mmin,mmax)

    # smooth systematics (initially clones of nominal)
    hup_smooth[syst] = template['hist'].Clone(hup_name)
    hdown_smooth[syst] = template['hist'].Clone(hdown_name)
    
    hup_smooth[syst].SetTitle(hup_name) 
    hdown_smooth[syst].SetTitle(hdown_name)

    # save copies of originals
    template['hsyst_up_orig'][syst] = template['hsyst_up'][syst].Clone(hup_name+'_orig_fine')
    template['hsyst_down_orig'][syst] = template['hsyst_down'][syst].Clone(hdown_name+'_orig_fine')
    
    template['hsyst_up_orig'][syst].SetTitle(hup_name+'_orig_fine') 
    template['hsyst_down_orig'][syst].SetTitle(hdown_name+'_orig_fine') 

  #---------- draw slplines -----------
  
  print "\tDrawing spline plots.."
  for syst in template['whsyst_up'] : 
    drawSpline(hup_ratio[syst],sup_ratio[syst],maxerr_up,outpath,tname=hname,reg=reg)
    drawSpline(hdown_ratio[syst],sdown_ratio[syst],maxerr_down,outpath,tname=hname,reg=reg)

  #-------- smooth systematics --------
  
  print "\tMaking new smoothed systematics and comparison plots.."
  for syst in template['whsyst_up'] :

    # smooth ratio
    makeSmoothRatio(hup_rsmooth[syst],sup_ratio[syst],nbins)
    makeSmoothRatio(hdown_rsmooth[syst],sdown_ratio[syst],nbins)
    
    # multiply by nominal
    hup_smooth[syst].Multiply(hup_rsmooth[syst])   
    hdown_smooth[syst].Multiply(hdown_rsmooth[syst])   

    # set poisson uncertainties on new systematic plot
    setPoissonUnc(hup_smooth[syst],nbins)
    setPoissonUnc(hdown_smooth[syst],nbins)
    
    # make comparison plots
    drawComparison(template['hsyst_up'][syst],hup_smooth[syst],outpath,tname=hname,reg=reg)
    drawComparison(template['hsyst_down'][syst],hdown_smooth[syst],outpath,tname=hname,reg=reg)
    
    # make ratio comparison plots
    up_dummy = hup_smooth[syst].Clone()
    up_odummy = template['hsyst_up'][syst].Clone()
    down_dummy = hdown_smooth[syst].Clone()
    down_odummy = template['hsyst_down'][syst].Clone()
    up_dummy.Divide(template['hist'])
    up_odummy.Divide(template['hist'])
    down_dummy.Divide(template['hist'])
    down_odummy.Divide(template['hist'])
    drawComparison(up_odummy,up_dummy,outpath,ratio=True,tname=hname)
    drawComparison(down_odummy,down_dummy,outpath,ratio=True,tname=hname)

    # rewrite with new systematic
    template['hsyst_up'][syst] = hup_smooth[syst]
    template['hsyst_down'][syst] = hdown_smooth[syst]

  
  #------------- root output -------------
  
  print "\tMaking crosscheck root output: "+outpath+hname+"_systsmooth.root"
  
  # making output file
  fout = TFile(outpath+hname+"_systsmooth.root","RECREATE")
  
  # --- nominal ---
  template['hist'].Write()
  template['h_coarse'].Write()
  if template['smooth_th']:
    template['h_coarse_fixed'].Write()
  
  # --- systematics ---
  for syst in template['whsyst_up'] :
   
    # up
    template['hsyst_up_orig'][syst].Write()
    template['hsyst_up_coarse'][syst].Write()
    hup_ratio[syst].Write()
    sup_ratio[syst].Write()
    hup_rsmooth[syst].Write()
    hup_smooth[syst].Write()

    # down
    template['hsyst_down_orig'][syst].Write()
    template['hsyst_down_coarse'][syst].Write()
    hdown_ratio[syst].Write()
    sdown_ratio[syst].Write()
    hdown_rsmooth[syst].Write()
    hdown_smooth[syst].Write()

    # if var bins, store fixed crosscheck
    if template['smooth_th']:
      template['hsyst_up_coarse_fixed'][syst].Write()
      template['hsyst_up_coarse_fixed'][syst].Write()
    
  # closing
  fout.Close()
  print "\tFinished! All "+hname+" smoothing output in: "+outpath
  
  #sys.exit()
  #if "ttbar" in hname : sys.exit()

  # end
  return

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def makeSpline(h,nbins):

  xx = []
  yy = []
  
  # for spline plot
  maxerr = 0

  # get bin array
  ibin = 1
  
  while ibin<nbins+1 :

    # getting coordinates
    w = h.GetBinWidth(ibin)
    le = h.GetBinLowEdge(ibin)
    x = le + (w/2)
    y = h.GetBinContent(ibin)
    e = h.GetBinError(ibin)

    # first bin
    if ibin == 1 :
      xx.append(le)
      yy.append(y)

    # regular bins
    xx.append(x)
    yy.append(y)

    # first bin
    if ibin == nbins :
      xx.append(le+w)
      yy.append(y)

    # max error (last bin can be large
    if (maxerr < e) and not (ibin == nbins) : 
      maxerr = e

    ibin += 1

  spline = TSpline3(h.GetTitle()+"_spline",array('d',xx),array('d',yy),len(xx))

  return spline, maxerr

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def makeSmoothRatio(h,spline,nbins):
  # looping over bins
  ibin = 1  
  while ibin<nbins+1 :

    w = h.GetBinWidth(ibin)
    x = h.GetBinLowEdge(ibin)+(w/2)

    # filling histogram with mid values of spline
    y  = spline.Eval(x)
    #print('coords:\t',x,' ',y)

    h.SetBinContent(ibin,y)
    ibin += 1

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def setVarBins(h,hname,nbins,th):
  
  print "Merging bins with large stat unc. for: "+hname

  # new values
  values = []
  unc = [] # uncertainties
  edges = []

  # initial values for crosscheck
  ivalues = []
  iunc = [] 
  iedges = []

  # temporary dummy values for merging purposes 
  t_y = 0
  t_e2 = 0
  prev_fail = False

  # --- getting bins ---
  # looping over bins to get new values
  ibin = 1  
  while ibin<nbins+1 :

    # get info
    y  = h.GetBinContent(ibin) # value
    e = h.GetBinError(ibin) # error
    stat = abs(e/y) # size of uncert

    # zero check (will be absorbed so no need for warning)
    if y <= 0 : y = 0.000000000001

    w = h.GetBinWidth(ibin)
    le = h.GetBinLowEdge(ibin)
   
    # get first edge
    if ibin == 1 : 
     edges.append(le)
     iedges.append(le)
    
    # store crosschecks
    ivalues.append(y)
    iunc.append(e)
    iedges.append(le+w)

    # --- tests ---
    # did previous fail?
    if prev_fail :
      
      # add values
      t_y += y
      t_e2 += pow(e,2)
      
      # now pass? 
      if (abs(np.sqrt(t_e2)/t_y) <= th) or (ibin == nbins) :
        # store
        values.append(t_y)
        unc.append(np.sqrt(t_e2))
        edges.append(le+w)
        #clear dummys
        t_y = 0
        t_e2 = 0
        prev_fail = False
    
    # if previous bin stored, do check
    elif ((stat > th)  or (stat == 0)) and not (ibin == nbins):
      # store dummys
      t_y += y
      t_e2 += pow(e,2)
      prev_fail = True

    else : 
      values.append(y)
      unc.append(e)
      edges.append(le+w)

    ibin += 1
  
  # --- checking last bin ---
  
  l_y = values[-1]
  l_e = unc[-1]

  while abs(l_e/l_y) > th :
    print "Stat. unc. of last bin larger than threshold, merging with previous"
    # delete last bin and connecting edge
    del edges[-2]
    del values[-1]
    del unc[-1]
    
    # merge last bin
    values[-1] = values[-1]+l_y
    unc[-1] = np.sqrt(pow(unc[-1],2)+pow(l_e,2)) 

    # get new last bin values
    l_y = values[-1]
    l_e = unc[-1]
    
  
  # --- making new histo ---
  nbins_new = len(values)

  # init
  h_new = TH1F(hname,hname,nbins_new,array('d',edges))
  
  # looping over bins
  ibin = 1  
  ival = 0
  while ibin<nbins_new+1 :

    # get value
    y = values[ival]
    e = unc[ival] 

    # fill values
    h_new.SetBinContent(ibin,y)
    h_new.SetBinError(ibin,e)
    
    ibin += 1
    ival += 1

  # --- crosscheck ---
  #print "\nOriginal:"
  #for i,v in enumerate(ivalues) :
  #  if v <= 0 : v = 0.000000000001
  #  print "Low: ",iedges[i],"\tHigh: ",iedges[i+1],"\tValue: ",round(v,2),"\tError: ",round(iunc[i],2),"\tRel. Unc: ",round(abs(iunc[i]/v),3)
  #
  #print "\nNew:"
  #for i,v in enumerate(values) : 
  #  # zero check
  #  print "Low: ",edges[i],"\tHigh: ",edges[i+1],"\tValue: ",round(v,2),"\tError: ",round(unc[i],2),"\tRel. Unc: ",round(abs(unc[i]/v),3)
  
  # --- end ---
  return h_new, edges 

#TODO: check last bin
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def rebinToEdges(h,hname,nbins,edges):

  print "Merging bins according to given edges for: "+hname

  # new values
  values = []
  unc = [] # uncertainties

  # initial values for crosscheck
  ivalues = []
  iunc = [] 
  iedges = []

  # temporary dummy values for merging purposes 
  t_y = 0
  t_e2 = 0
  merge_prev = False

  # --- getting bins ---
  # looping over bins to get new values
  ibin = 1 
  ie = 0
  while ibin<nbins+1 :

    # get info
    y  = h.GetBinContent(ibin) # value
    e = h.GetBinError(ibin) # error

    w = h.GetBinWidth(ibin)
    le = h.GetBinLowEdge(ibin)
    te = le + w # top edge
    
    # store crosschecks
    ivalues.append(y)
    iunc.append(e)
    if ibin == 1 : 
     iedges.append(le)
    iedges.append(le+w)

    # --- tests ---
    # was it merged with previous? 
    if merge_prev :
      
      # add values
      t_y += y
      t_e2 += pow(e,2)
      
      # is this last merge bin?
      if te in edges :
        # store
        values.append(t_y)
        unc.append(np.sqrt(t_e2))
        #clear dummys
        t_y = 0
        t_e2 = 0
        merge_prev = False
    
    # if previous bin stored, check if this needs to merge
    elif te not in edges :
      # store dummys
      t_y += y
      t_e2 += pow(e,2)
      merge_prev = True

    else : 
      values.append(y)
      unc.append(e)

    ibin += 1
  
  # --- making new histo ---
  nbins_new = len(values)

  # init
  h_new = TH1F(hname,hname,nbins_new,array('d',edges))
  
  # looping over bins
  ibin = 1  
  ival = 0
  while ibin<nbins_new+1 :

    # get value
    y = values[ival]
    e = unc[ival] 

    # fill values
    h_new.SetBinContent(ibin,y)
    h_new.SetBinError(ibin,e)
    
    ibin += 1
    ival += 1

  # --- crosscheck ---
  #print "\nOriginal:"
  #for i,v in enumerate(ivalues) : 
  #  if v <= 0 : v = 0.000000000001
  #  print "Low: ",iedges[i],"\tHigh: ",iedges[i+1],"\tValue: ",round(v,2),"\tError: ",round(iunc[i],2),"\tRel. Unc: ",round(abs(iunc[i]/v),3)
  #
  #print "\nNew:"
  #for i,v in enumerate(values) : 
  #  print "Low: ",edges[i],"\tHigh: ",edges[i+1],"\tValue: ",round(v,2),"\tError: ",round(unc[i],2),"\tRel. Unc: ",round(abs(unc[i]/v),3)
  
  # --- end ---
  return h_new 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def setPoissonUnc(h,nbins):
  # looping over bins
  ibin = 1  
  while ibin<nbins+1 :

    # get value
    y  = h.GetBinContent(ibin)
    #print('value:\t',y)

    h.SetBinError(ibin,np.sqrt(y))
    ibin += 1

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def drawSpline(h,spline,maxerr,outpath,tname="",reg="") :

  # gStyle Options 
  gStyle.SetOptStat(0)
 
  # title
  hname = h.GetTitle()
  subtitle = other_helpers.getPrettySystName(h.GetTitle(),tname)
  subtitle = re.sub('ratio',' ',subtitle)
  tlabel=''
  if reg : tlabel = reg

  # legend
  #pos =[]
  pos=[0.68,0.75,0.88,0.87]

  l = plot_helpers.makeLegend([h,spline],
                              ['template','spline'],
                              pos=pos
                              #ncols=2
                              )

  # cosmetics
  if 'Higgs' in hname : h.SetMaximum(h.GetMaximum()+3*maxerr)
  else : h.SetMaximum(h.GetMaximum()+2*maxerr)
  h.SetMinimum(h.GetMinimum()-1.5*maxerr)
  h.GetXaxis().SetTitle('Large-R Jet Mass [GeV]')
  h.GetYaxis().SetTitle('Systematic / Nominal / 5 GeV')
  
  h.GetXaxis().SetLabelSize(0.045)
  h.GetYaxis().SetLabelSize(0.045)
  h.GetXaxis().SetTitleSize(0.045)
  h.GetYaxis().SetTitleSize(0.045)

  # making plot
  c = TCanvas("c"+hname,"c"+hname,1600,1000) 
  h.Draw()
  spline.Draw('samec')
  l.Draw()

  # atlas label
  plot_helpers.drawAtlasLabel(0.13,0.85,
                              p_type="Preliminary Internal",
                              s="13TeV",
                              slabel=subtitle,
                              tlabel=tlabel,
                              lumi="136 fb^{-1}")

  # printing plot
  c.Print(outpath+'spline_'+hname+".pdf")


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def drawComparison(h0,h,outpath,ratio=False,tag="",tname="",reg="") :

  # gStyle Options 
  gStyle.SetOptStat(0)
 
  # title
  hname = h.GetTitle()
  subtitle = other_helpers.getPrettySystName(h.GetTitle(),tname)
  tlabel=''
  if reg : tlabel = reg

  # colors
  plot_helpers.setHistoColors([h0,h],fill=False,colors=[3,5])


  # legend
  #pos =[]
  pos=[0.65,0.75,0.85,0.87]

  l = plot_helpers.makeLegend([h0,h],
                              ['Original','Smoothed'],
                              pos=pos
                              #ncols=2
                              )

  # cosmetics 
  h0.GetXaxis().SetTitle('Large-R Jet Mass [GeV]')
  if ratio : 
    h0.GetYaxis().SetTitle('Systematic / Nominal / 0.5 GeV')
    h0.SetMaximum(1.8)
    h0.SetMinimum(0.5)
  else : 
    h0.GetYaxis().SetTitle('Events / 0.5 GeV') #TODO: not hardcode this
    h0.SetMaximum(h.GetMaximum()*1.4)
  
  h0.GetXaxis().SetLabelSize(0.045)
  h0.GetYaxis().SetLabelSize(0.045)
  h0.GetXaxis().SetTitleSize(0.045)
  h0.GetYaxis().SetTitleSize(0.045)

  # making plot
  c = TCanvas("c"+hname,"c"+hname,1600,1000) 
  h0.Draw('hist')
  h.Draw('histsame')
  l.Draw()

  # atlas label
  plot_helpers.drawAtlasLabel(0.13,0.85,
                              p_type="Preliminary Internal",
                              s="13TeV",
                              slabel=subtitle,
                              tlabel=tlabel,
                              lumi="136 fb^{-1}")

  # printing plot
  if ratio : tag = tag+"ratio"

  c.Print(outpath+tag+'comparison_'+hname+".pdf")




