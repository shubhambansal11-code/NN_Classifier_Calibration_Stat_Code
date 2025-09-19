# ==========================================
# XML Model Creation Helpers:
#   Data Manipulation
#
# Original author: Migle Stankaityte
# 13 May 2019
# Refurbished for QT PhD 2019-2020
# ==========================================

#import
import sys #, os, re
#import json

import numpy as np
import random

import model_helpers as helpers

import ROOT
from ROOT import gROOT, TFile, TTree, TH1F, TKey
gROOT.SetBatch(True)

# ==========================================
# Getting data fistograms from Ntuple or Histo
# Depending on flags make histos as needed
# void

def getDataHistos(data,flags,tree,nbins,m_min,m_max,pt_min,pt_max,slice_calc,cand) :
  
  print '\n----------------------------------' 
  # adding candidate to region
  if data['region'] and cand :
    data['region'] += cand
    data['cand'] = cand
    print 'Filling data '+data['region']+' histogram:'
  else :
   print 'Filling data histogram:'
  
  # opening file
  f_data = TFile(str(data['path']),"READ")

  # ----------------------------------------
  # data from histogram input
  if data['hist_input'] : 

    print '\tHistogram input:'
    
    # getting histograms in input file
    histlist = f_data.GetListOfKeys()
       
    # single histogram input
    if data['hist_input'] == 1 :
      print '\t\t Single data histogram.'
      #assert(histlist.Contains("data")),'No full-data histogram in file.' 
      assert(histlist.Contains("data")),'No full-data histogram in file.' 
      data['hist'] = f_data.Get('data')
      data['hist'].GetXaxis().SetRangeUser(m_min,m_max)
      checkHistOpt([data['hist']],nbins,m_min,m_max)
    
    # slice input
    else :
      print '\t\tTaking data histogram slices:'
      data['shists'] = []
     
      dtitle = ""
      # loop over all histograms in file
      for key in histlist :
        dtitle = key.GetName()
        # getting slices
        if 'data_s' in dtitle : 
          print '\t\t\tGetting ',dtitle
          data['shists'].append(f_data.Get(dtitle))
      
      # getting total slice number
      assert ("_of" in dtitle),'Error! Total # of slices not in slice name!'
      data['n_slice'] = dtitle.split('_of',1)[1]
      print '\t\t\tNumber of slices in total: ',data['n_slice']

      # checking histogram input validity
      checkHistOpt(data['shists'],nbins,m_min,m_max)

  # ----------------------------------------
  # data from ntuple -- sliced
  elif data['slice'] :

    print '\tSlicing data from NTuple:'
    # setting number of slices in output
    data['n_out'] = setOutSlices(data['n_out'],data['n_slice'])
    
    # initialising histograms
    data['shists'] = []
    for i in range(data['n_out']) :
      dtitle = 'data_s'+str(i)+'_of'+str(data['n_slice'])
      data['shists'].append(TH1F(dtitle,dtitle,nbins,m_min,m_max))
    
    # create slices
    fillDataSlices(f_data,data['shists'],pt_min,pt_max,data['n_slice'],tree,data['region'],flags['conf'])
  
  # ----------------------------------------
  # data from ntuple -- not slices
  else :
    print '\tNtuple input:'
    data['hist'] = TH1F('data','data',nbins,m_min,m_max)
    helpers.fillHisto(f_data,data['hist'],-1,1e10,pt_min,pt_max,[],tree,data['region'],flags['conf'],'w',datamc=data['datamc'])

  # ----------------------------------------
  
  # closing file
  f_data.Close()
  
  # SR MC histogram container for slice calculation
  if slice_calc :
    data['srhist'] = TH1F('SR_qcd','SR_qcd',nbins,m_min,m_max)

  # if subtraction is performed, make a copy for cross-check
  if data['sub'] :
    if 'shists' in data :
      data['oshists'] = []
      for i, h, in enumerate(data['shists']) :
        otitle = h.GetTitle()+'_orig'
        oh = h.Clone(otitle)
        oh.SetTitle(otitle)
        data['oshists'].append(oh)
    else :
      data['ohist'] = data['hist'].Clone('data_orig')
      data['ohist'].SetTitle('data_orig')
  
  # ----------------------------------------

  # adding histogram name info to xml card info json
  data['dataname'] = []
  
  if flags['mk_xmljson'] : 
    if 'hist' in data : data['dataname'].append('data')
  
    # slices
    if 'shists' in data :
     
      # if number of slices to run = all or -1
      if (int(data['nrun']) < 0) or (int(data['nrun']) >= len(data['shists'])) :
        for h in data['shists']:
            data['dataname'].append(h.GetTitle())

      # if number of slices to run < all
      elif int(data['nrun']) < len(data['shists']) :
        random.seed(111)
        i=0
        while True:
          choice_name = random.choice(data['shists']).GetTitle()
          if choice_name not in data['dataname'] : 
            data['dataname'].append(random.choice(data['shists']).GetTitle())
            i += 1
          if i == int(data['nrun']) : break

      # if not specified - will do single random slice

  return

# end of: getDataHistos()
# ==========================================
# Template subtraction from data
# void

def subTempFromDataHist(data,templates) :

  print '\n----------------------------------' 
  print 'Template subtraction from data using gausian smearing:\n'
  
  # subtracted integrals
  sub_temp_int = []
  sub_gauss_int = []

  for i, st in enumerate(data['sub_t']):
    
    # is subtraction template found
    t_found = False

    # looping over templates
    for t in templates :
      if st == t['sample'] : 
        t_found = True
        print '\tSubtracting template ',t['sample']
        
        # subtract from single histogram
        if 'hist' in data : 
          data['hist'].Add(t['sub_hist'],-1)
          sub_temp_int.append(t['hist'].Integral())
        # subtract from slices 
        if 'shists' in data : 
          for h in data['shists'] : h.Add(t['sub_hist'],-1)
          sub_temp_int.append(t['hist'].Integral()/float(data['n_slice']))
        
        sub_gauss_int.append(t['sub_hist'].Integral())
    
    assert (t_found), 'Subtraction template '+st+' not found in template list.'
    
  # printing to screen
  print '\n * * * * * * * * * *'
  print ' Cross-check integrals:'
  
  if 'hist' in data : print '\tOrig. data:',data['ohist'].Integral() 
  if 'shists' in data : print '\tOrig. data (first slice): ',data['oshists'][0].Integral() 

  print '\tTemplates:'
  for i, gi in enumerate(sub_gauss_int) :
    tdiff = sub_temp_int[i]-gi
    tdiff_pc = round(tdiff/sub_temp_int[i],2)
    print '\t\t',data['sub_t'][i],'\torig: ',round(sub_temp_int[i],2),
    print '\tsub: ',round(gi,2),'\tdiff: ',round(tdiff,2),'\t=',tdiff_pc,'%'   
    sys.stdout.flush()

  if 'hist' in data : print '\tFinal data:',data['hist'].Integral() 
  if 'shists' in data : print '\tFinal data (first slice):',data['shists'][0].Integral() 
  
  if 'hist' in data :
    diff = data['ohist'].Integral()-sum(sub_gauss_int)-data['hist'].Integral()
    diff_pc = round(diff/data['ohist'].Integral(),2)
    print '\tDiff: ',round(diff,2),'\t = ',diff_pc,'%' 
  
  if 'shists' in data :
    diff = data['oshists'][0].Integral()-sum(sub_gauss_int)-data['shists'][0].Integral()
    diff_pc = round(diff/data['oshists'][0].Integral(),2)
    print '\tDiff (first slice): ',round(diff,2),'\t = ',diff_pc,'%' 

  sys.stdout.flush()
  return

# end of: subTempFromDataHist()
# ==========================================
# Data Slice Calculation using Histograms
# void

def dataSliceCalc(data,tree,conf) :

  print '\n----------------------------------' 
  print 'Data slice calculation using histograms:\n'
  
  # slice integrals
  slice_int = []

  #-----------------------------------------
  # calculating number of slices from MC
  
  print 'Getting number of events from QCD MC:'
  # getting SR MC tree
  sr_mc = TFile(data['sr_mc'],"READ")
  helpers.fillHisto(sr_mc,data['srhist'],-1,1e10,pt_min,pt_max,tree,'sr'+data['cand'],conf,'w')
   
  # getting number of events in SR and VR
  n_sr = data['srhist'].Integral() 
  n_cr = data['hist'].Integral()
  
  # calculating number of slices
  ratio = n_sr/n_cr
  data['n_slice'] = round(1/ratio)
 
  summary = {
      'SR' : n_sr,
      'VR' : n_cr,
      'ratio' : ratio, 
      'slices' : data['n_slice']
      }

  print 'Data slice calculation is done (summary at the end)!\n'
  
  return summary

# end of: dataSliceCalc()
# ==========================================

# Set Output Slices

def setOutSlices(nout, nslice) :

  # setting number or output slices
  
  # if more slices requested in output than produced, set out to max
  if nout > nslice :
    print '\nWARNING: Number of output slices requested (',nout,
    print ') is larger than slices made, all slices will be saved.\n'
    nout_new = nslice
  # else keep the same
  else : nout_new = nout
  print '\n\tSlices in output: ',nout_new
  
  return nout_new

# end of: setOutSlices()
# ==========================================
# Data Slicing
# Needs given number of slices as it is performed at ntuple level!
# void

def fillDataSlices(infile, hists, pt_min, pt_max, nslice, tree_name, dirname='',conf=False) :
  print '\tFilling slices:'
  
  # get tree
  if dirname : tree = infile.Get(str(dirname)+'/'+tree_name)
  else : tree = infile.Get(tree_name)

  # print
  nentries = tree.GetEntries()
  print '\t\tEntries in tree: ', nentries
  sys.stdout.flush()

  # filling array with random numbers
  r = []
  random.seed(111)
  for i in range(nentries):
    r.append(random.uniform(0.,nslice))

  # fill histograms
  for i, event in enumerate(tree) :
   
    if (nentries > 1000000) and (not i%1000000) : 
      print '\t\t\tFilled ',i,'events out of ',nentries
      sys.stdout.flush()
    for ns, h in enumerate(hists) :
      if (r[i] > ns) and (r[i] < ns+1) : 
        if conf :
          if (event.Hcand_p4.Pt()>pt_min) and (event.Hcand_p4.Pt()<pt_max):
            h.Fill(event.Hcand_p4.M())
        else :
          if (event.Hcand_pt>pt_min) and (event.Hcand_pt<pt_max):
            h.Fill(event.Hcand_m)



  # get total integral for consistency check
  totint = 0 
  for h in hists : 
    # reset directory
    h.SetDirectory(0)
    # debug -- print each histogram object and integral
    # print 'Integral of ',h,': ',h.Integral()
    totint += h.Integral()

  print'\t\tTotal Itegral of all slices: ',totint

  return

# end of: fillDataSlices()
# ==========================================
# Checking histogram input options to agree with global
# void

def checkHistOpt(datahists,nbins,m_min,m_max) :

  for h in datahists :
    
    # making sure input agrees with our settings
   # assert (h.GetNbinsX() == nbins), 'Number of bins in input data does not match run sequence'
   # assert (h.GetXaxis().GetXmin() == m_min), 'Input data x_min does not match run sequence'
   # assert (h.GetXaxis().GetXmax() == m_max), 'Input data x_max does not match run sequence'
  
    h.SetDirectory(0)

  return

# end of: checkHistOpt()
# ==========================================

# end of file
# ==========================================
# ==========================================
