# ==========================================
# XML Model Creation Helpers
#
# Original Author: Migle Stankaityte
# 18 March 2019
# Refurbished for QT PhD 2019-2020
# ==========================================

#import
import sys, os, re
import json

import numpy as np

import ROOT
from ROOT import gROOT, TFile, TTree, TH1F, TLorentzVector
gROOT.SetBatch(True)

# ==========================================
# Reading json file with model information
# returns: data path, samples paths and systematics

def readModelInfo(jpath): 
  
  data_path = ''
  templates = {}
  flags = {}
  
  print '\n----------------------------------'
  print 'Reading model config json file: '+jpath
  
  # --------------------------------
  with open(jpath, 'r') as f:
  
    # --- reading json file ---
    inf = json.load(f)

    # --- getting flags ---
    if 'flags' in inf :
      flags = inf['flags']
    else : flags = {}
   
    # initialising
    if 'nominal' not in flags : flags['nominal'] = False
    #if 'injection_test' not in flags : flags['injection_test'] = False
    if 'no_ws' not in flags : flags['no_ws'] = False
    if 'no_cc' not in flags : flags['no_cc'] = False
    if 'no_dataout' not in flags : flags['no_dataout'] = False
    if 'mkdatain' not in flags : flags['mkdatain'] = False
    if 'slice_calc' not in flags : flags['slice_calc'] = False
    if 'updown' not in flags : flags['updown'] = False
    if 'basicsymm' not in flags : flags['basicsymm'] = False
    if 'norm_systs' not in flags : flags['norm_systs'] = False
    if 'norm_all' not in flags : flags['norm_all'] = True
    if 'mk_xmljson' not in flags : flags['mk_xmljson'] = True # make json file xml card creation
    if 'slicetemp' not in flags : flags['slicetemp'] = False
    if 'conf' not in flags : flags['conf'] = False
    if 'sepuncerts' not in flags : flags['sepuncerts'] = True # build separate FTAG syst uncertainties

    # slice calculation - special run type (no general outputs)
    if flags['slice_calc'] : 
      flags['nominal'] = True 
      flags['no_ws'] = True 
      flags['no_dataout'] = True 
      flags['injection_test'] = False 
      flags['mk_xmljson'] = False 
    
    # --- reading data ---
    assert ('data' in inf), 'No data object provided, aborting.'
    data = inf['data']
    assert ('path' in data), 'No data path provided, aborting.'
   
    # hist input (0 - not, 1 - single, 2 - sliced)
    if 'hist_input' not in data : data['hist_input'] = 0
    assert ((type(data['hist_input']) == int) and (data['hist_input']<=2) and (data['hist_input']<=2) and (data['hist_input']>=0)), 'Invalid data[\'hist_input\'] value: '+str(data['n_out'])+' must be 0, 1 (full) or 2 (sliced)'

    # injection
    if 'injection' not in data : data['injection'] = None

    # is data actually MC (e.g. slice calc.)
    if 'datamc' not in data :
      data['datamc'] = False
    
    # region (not necessary, e.g. not specified in conf. note)
    if 'region' not in data :
      data['region'] = None

    # template subtraction from data
    if 'sub' not in data : 
      data['sub'] = False
    
    if data['sub']:
      # conditions
      assert (not flags['mkdatain']), 'Subtraction not possible when making a data input file for modelMaker.'
      assert ('sub_t' in data and data['sub_t']), 'Subtraction on, but templates not provided, aborting.'
      # if fitting of subtracted templates not defined - set all to false
      if 'sub_fit' not in data : data['sub_fit'] = [False for i in data['sub_t']]
    else : 
      data['sub_t'] = []
      data['sub_fit'] = []
    
    # data slicing and slice calculation
    # init
    if 'slice' not in data : 
      data['slice'] = False

    # slice calc.
    if flags['slice_calc'] :
      # conditions
      assert('sr_mc' in data),'SR MC file has to be given for slice calculation!'
      # slicing turned off
      if data['slice'] : 
        print 'WARNING! Slicing procedure is turned off for slice calculation.'
        data['slice'] = False
    # SR_mc only needed when slice calculation is done
    else : data['sr_mc'] = None
       
    # slicing
    if data['slice'] : 
      # conditions
      assert(not data['hist_input']),'Data slicing procedure is turned on with histogram input! Not possible!'
      assert ('n_slice' in data), 'Data slicing procedure turned on, but number of slices is not provided, aborting!'
      assert (not data['datamc'] ),'Cannot slice pseudo-data from MC'
      
      # if no. of slices in output not given, all will be saved by default
      if 'n_out' not in data or (data['n_out'] < 0): data['n_out'] = data['n_slice']
      elif data['n_out']:
        # check if integer
        assert (type(data['n_out']) == int), 'Number of output slices is not integer: '+str(data['n_out'])
    else :
      data['n_slice'] = None
      data['n_out'] = None

    if 'sr_mc' not in data : data['sr_mc'] = None
    if 'n_slice' not in data : data['n_slice'] = None
    
    # number of slices to run over in genxml( will only take affect if sliced)
    if 'nrun' not in data : data['nrun'] = 1

    sys.stdout.flush()
    
    # --- reading templates ---
    templates=[]

    if 'templates' in inf :
       
      for t in inf['templates'] : 
        
        assert ('sample' in t), 'No sample name provided, aborting.'
        assert ('path' in t), 'No sample path provided, aborting.'

        # dont record if template is a systematic and nominal flag is on
        if 'syst_for' not in t : t['syst_for'] = []
        if t['syst_for'] and flags['nominal'] : 
          print t['sample'],' is a systematic template, run set to nominal only, not recording.'
          continue

        # setting empty flags
        if 'syst' not in t or flags['nominal'] : t['syst'] = []
        if 'systpath' not in t or flags['nominal'] : t['systpath'] = ""
        if 'smoothsyst' not in t  or flags['nominal']: t['smoothsyst'] = False
        if 'Vjets_wsyst_up' not in t or flags['nominal'] : t['Vjets_wsyst_up'] = [] 
        if 'Vjets_wsyst_down' not in t or flags['nominal']: t['Vjets_wsyst_down'] = [] 
        if 'dsids' not in t: t['dsids'] = []
        if 'weightlist' not in t or flags['nominal'] : t['weightlist'] = []
        if 'fit' not in t : t['fit'] = False
        if 'smooth_th' not in t : t['smooth_th'] = 0
        if 'combine' not in t : t['combine'] = []
        if 'injection' not in t : t['injection'] = None
        if 'ignoreWS' not in t : 
          if t['syst_for'] : t['ignoreWS'] = True
          else : t['ignoreWS'] = False
        if 'region' not in t : t['region'] = None
        if 'spurious' not in t : t['spurious'] = False
        if 'fixed' not in t : t['fixed'] = False
        if 'truth_pt_min' not in t: t['truth_pt_min'] = -1
        if 'truth_pt_max' not in t: t['truth_pt_max'] = 1e10
        if 'mu' not in t: t['mu'] = None
    

        ## weight cuts (can't set to false/none as 0 can be a boundary)
        if 'w_min' not in t: t['w_min'] = 'unset'
        else : t['w_min'] = float(t['w_min']) # make sure it is a number

        if 'w_max' not in t: t['w_max'] = 'unset'
        else : t['w_max'] = float(t['w_max']) # make sure it is a number

        if 'truth_reweighting' not in t: t['truth_reweighting'] = 'unset'

        ## noting wether template will be subtracted
        if t['sample'] in data['sub_t'] : 
          t['sub'] = True
          t['sub_fit'] = data['sub_fit'][data['sub_t'].index(t['sample'])]
          if t['sub_fit'] and not t['fit']:
             t['fit'] = True
             print 'WARNING: Template to be fitted before subtraction from data, but template fitting is off!'
             print 'Teplate fitting is turned on!'
        else : 
          t['sub'] = False
          t['sub_fit'] = False

        templates.append(t)

  # --------------------------------
  # screen output

  # --- unusual run-modes ---
  print '\n**********************'
  if flags['slice_calc'] : 
    print 'Slice calculation mode!'
    print 'Some flags may be changed automatically (see flag output below).'
    print 'Only template fit output files will be made if turned on!'
  if flags['mkdatain'] : 
    print 'Making a data histogram input for faster running of modelMaker!' 
  
  # --- flags ---
  if flags :
    print '**********************\nFlags:'
    for flag, setto in flags.iteritems():
      print '\t',flag,': ',setto

  sys.stdout.flush()

  # --- data ---
  print '\n**********************'
  print '\nData: '+data['path']
  if data['hist_input'] : 'Input file is a histogram (not ntuple)!'
  if data['sub'] :
    print '\n\tTemplate subtraction will be done.\n\tSubtracted templates:'
    for i in range(len(data['sub_t'])): print '\t\t',data['sub_t'][i],'\t(Smoothed: ',data['sub_fit'][i],')'
  if data['slice'] : 
    print '\n\tData slicing will be performed'
    if data['sr_mc'] : print '\t\tSlice number will be defined by SR MC.\n\t\tPath: ',data['sr_mc']
    else : print '\t\tSlice number: ',data['n_slice'],' slices.'
    print '\t\tSlices to be saved: ',data['n_out']
    print '\t\tNote, if this number is larger than the total number of slices - all available slices will be saved.'
 
  sys.stdout.flush()
  print '\n**********************\n'

  # --- templates --- 
  i = 0
  print '\nTemplates:'
  for t in templates:
    i += 1
    
    print '%d. '%(i)+t['sample']+': '+t['path']
    
    print '\n\t'+t['sample']+' systematics:'
    if t['syst'] : 
      for syst in t['syst']:
        print '\t\t'+syst
    else : print '\t\tNominal Only!'
    if t['weightlist'] : 
      for w in t['weightlist']:
        print '\t\t'+w
    else : print '\t\tNo alternative btag uncertainty weights in sample'


    if t['fit']: print '\n\tParametric fit will be applied.'
    else : print '\n\tNo parametric fitting of template!' 
    
    if t['combine'] : 
      print '\n\tWill be combined into templates:'
      for c in t['combine'] : print '\t\t',c
      print '\n'
    else : print '\n\tTemplate is standalone - not to be combined.\n'
    
    if t['syst_for'] : 
      print '\tThis template is a shape systematic for:'
      for nt in t['syst_for'] : print '\t\t',nt
      #if 'truth_reweighting' in t: print '\twith truth reweighting given by: ',t['truth_reweighting']
      print '\n'
    else : print '\tTemplate is not a shape systematic.\n'

    if t['injection'] :
      print '\t'+str(t['sample'])+' injected mu:'+str(t['injection'])+'\n'
    else:
      print '\t'+str(t['sample'])+': no signal injection performed\n'
    
    if t['sub'] :
      print '\tTemplate will be subtracted from data',('smoothed\n' if t['sub_fit'] else 'without smoothing\n')

    if t['ignoreWS'] : print '\tTemplate will not be added to the workspace!\n'

  sys.stdout.flush()
  
  
  print '\n----------------------------------'
  return data, templates, flags

# end of: ReadModelInfo()
# ==========================================
# Printing paths to screen

def printOutPaths(flags,data_out,ws_pfx,model_name,cc,injtest_out,templates=[]) :

  print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
  print 'Output in:'
  if not flags['no_dataout'] : print 'data histogram: '+data_out+'.root'
  if not flags['no_ws'] :
    if templates :
      print 'template workspaces: '
      for t in templates :
        print '\t'+ws_pfx+'_'+t+'_combined_'+model_name+'_model.root'
    else :
      print 'template workspaces: '+ws_pfx+'_<template_name>_combined_'+model_name+'_model.root'
  if not flags['no_cc'] :
    print 'sanity check output: '+cc+'.root'
  if flags['injection_test'] : 
    print 'signal injection toy output: '+injtest_out+'.root\n'

  sys.stdout.flush()
  return

# ==========================================
# Adding output information to objects

def alocatePaths(flags,data,templates,asimov,pois,data_out,ws_pfx,model_name,injtest_out) :

  if not flags['no_dataout'] : 
    data_out_rel = re.sub('xmlAnaWSBuilder/','',data_out)
    data['outpath'] = data_out_rel+'.root'
  if not flags['no_ws'] :
    for t in templates :
      # skiping sub-templates and ignored templates
      if t['combine'] or t['ignoreWS'] : continue
      t['outpath'] = '../'+ws_pfx+'_'+t['sample']+'_combined_'+model_name+'_model.root'
  if flags['injection_test'] : 
    injtest_out_rel = re.sub('xmlAnaWSBuilder/','',injtest_out)
    asimov['outpath'] = injtest_out_rel+'.root'
    pois['outpath'] = injtest_out_rel+'.root'

  return

# end of: ReadModelInfo()
# ==========================================
# Writing Histograms from Ntuples
# void 

def fillHisto(infile, h, truth_pt_min, truth_pt_max, pt_min, pt_max, dsids, tree_name, dirname, conf, weight, w_min='unset', w_max='unset', truth_RW='unset', datamc=False):

  sample = h.GetTitle()
  # clone of h to store integral before capping of weights
  g = h.Clone()

  weight_name = weight
  if truth_RW != 'unset': weight_name += ' * '+truth_RW

  if ('data' in sample) and not datamc : print '\tFilling ',sample,' histogram with ',tree_name,' tree info.' 
  else : print '\tFilling ',sample,' histogram with ',tree_name,' tree info and weight ',weight_name,'.'

  # get tree
  if dirname : tree = infile.Get(str(dirname)+'/'+tree_name)
  else : tree = infile.Get(tree_name)

  # print
  nentries = tree.GetEntries()
  print '\t\tEntries in tree: ', nentries
  sys.stdout.flush()

  # fill histogram
  for i, event in enumerate(tree) :
   
    if (nentries > 1000000) and (not i%1000000) : 
      print '\t\t\tFilled ',i,'events out of ',nentries
      sys.stdout.flush()
    # data (no weight)
    if ('data' in sample) and not datamc :
      if conf :
        if (event.Hcand_p4.Pt()>pt_min) and (event.Hcand_p4.Pt()<pt_max):
          h.Fill(event.Hcand_p4.M())
      else :
        if (event.Hcand_pt>pt_min) and (event.Hcand_pt<pt_max):
          h.Fill(event.Hcand_m)


    # MC (templates or qcd MC)
    else :
      w = getattr(tree,str(weight))
      # is there any truth reweighting to be applied?
      if truth_RW!='unset' :
        rw = getattr(tree,str(truth_RW))
        w *= rw

      if conf :
        print '\tWARNING: conf inputs do not have truth pT, going to ignore truth_pt_min and truth_pt_max.'
        if (event.Hcand_p4.Pt()>pt_min) and (event.Hcand_p4.Pt()<pt_max):
          g.Fill(event.Hcand_p4.M(),w)
          if passWeightCut(w_min,w_max,w):
            h.Fill(event.Hcand_p4.M(),w)
      else :
        if (truth_pt_min!=-1) or (truth_pt_max!=1e10):
          if ( (dsids==[]) or (event.mcChannelNumber in dsids) ):
            if (event.Hcand_pt>pt_min) and (event.Hcand_pt<pt_max) and (event.truth_pt>truth_pt_min) and (event.truth_pt<truth_pt_max):
              g.Fill(event.Hcand_m,w)
              if passWeightCut(w_min,w_max,w):
                h.Fill(event.Hcand_m,w)
        else:
          if ( (dsids==[]) or (event.mcChannelNumber in dsids) ):
            if (event.Hcand_pt>pt_min) and (event.Hcand_pt<pt_max):
              g.Fill(event.Hcand_m,w)
              if passWeightCut(w_min,w_max,w):
                h.Fill(event.Hcand_m,w)
            
  if ('data' not in sample) :
    if (w_min != 'unset') or (w_max != 'unset') :
      print "\tEvent weight cuts were present"
      print "\t\tFull integral: ",g.Integral()
      print "\t\tCapped integral: ",h.Integral()
      print "\t\tSetting integral to original."
      h.Scale(g.Integral()/h.Integral())

  h.SetDirectory(0)

  print '\t'+sample+' histogram filled.\n'

  sys.stdout.flush()
  return 

# end of: fillHisto()
# ==========================================
# BELONGS TO : fillHisto
# returns true if passed

def passWeightCut(w_min, w_max, w):
  passCut = True

  if (w_min != 'unset') and (w < w_min) : 
    passCut = False
  elif (w_max != 'unset') and (w > w_max) : 
    passCut = False
  
  return  passCut

# ==========================================
# Check if item exists in file

def existsInFile(inpath,checklist,dirname) :
  
  exists = {}
  
  infile = TFile(inpath,"READ")
  # get tree
  if dirname : indir = infile.Get(str(dirname))
  else : indir = infile
  
  for item in checklist :
    exists[str(item)] = indir.GetListOfKeys().Contains(item)

  return exists

# end of: existsInFile()
# ==========================================
# Writing Histograms to file
# void 

def makeHistoRoot(outfile, hists):

  print 'Storing histograms in root file:'

  # making file
  outpath = outfile+'.root'
  f = TFile(outpath,"RECREATE")
  
  for h in hists :
    print '\tWriting '+h.GetTitle()
    h.Write()

  f.Close()

  print '\nDone. File at: '+outpath
  print '----------------------------------'

  return 

# end of: makeHistoRoot()
# ==========================================
# Making Directories
# void 

def mkdir(path):

  print 'Creating Directory: ',path
  if not os.path.exists(path):
    os.mkdir(path)
    print 'Directory ',path,' created.\n'
  else:    
    print 'Directory ',path,' already exists.\n'
  return 

# end of: mkdir()
# ==========================================
# Ensuring no negative bins
# void 

def zeroNegBins(h):

  i_int = h.Integral()
  nbins = h.GetNbinsX()

  # looping over bins
  ibin = 1  
  
  while ibin<nbins+1 :

    # getting bin content
    y = h.GetBinContent(ibin)
    
    # if negative set to infinitecimal ammount
    if y<=0 :
      
      x0 = h.GetBinLowEdge(ibin)
      x = x0 + h.GetBinWidth(ibin)
      
      y = 0.00000000001
      print 'WARNING!: in ',h.GetTitle(),' bin value <= 0 in range ',x0,' - ',x,', setting y = 0.00000000001!'

      h.SetBinContent(ibin,y)
    
    ibin += 1

  # rescaling to match the integral of original histogram
  h.Scale(i_int/h.Integral())
  
  return 

# end of: zeroNegBins()
# ==========================================
# Gaussian smearing of histogram
# Needed for template subtraction from data
# returns: smeared histogram

def smearGauss(h,tag):

  # initialising
  name = h.GetTitle()
  h_out = h.Clone(name+'_'+tag)
  h_out.SetTitle(name+'_'+tag)

  # looping over bins
  ibin = 1  
  nbins = h.GetNbinsX()

  #init. subtracted integral
  int_sub = 0.
  
  while ibin<nbins+1 :

    # getting bin content
    mu = h.GetBinContent(ibin)
    if mu < 0 : mu = 0
    
    # stat. uncertainty
    sigma = np.sqrt(mu)

    # gausian smearing 
    if mu : y = np.random.normal(mu,sigma)
    else : y = 0

    # only positive bins kept
    if y < 0 : y = 0

    h_out.SetBinContent(ibin,y)

    ibin += 1
 
  return h_out

# end of: smearGauss()
# ==========================================
# Printing all histogram objects for debug purposes
# void
def printDebugHists(data,templates):

  print 'Debug output:\n'
  # --- data ---
  if 'hist' in data : 
    print 'Data:'
    print data['hist']
  if 'ohist' in data : print data['ohist']
  
  if 'shists' in data:
    print 'Data slices:'
    for i, h in enumerate(data['shists']) : 
      print h
      if 'oshists' in data : print data['oshists'][i]
    
  # --- templates ---
  print '\nTemplates:'

  for template in templates:

    print '\n',template['sample']

    if 'hist' in template :
      
      if template['fit'] : 
        print template['ohist']
      print template['hist']
      
      if template['sub'] : print template['sub_hist']
     
      if template['syst'] :
        for sname in template['hsyst_up'] : 
          print '\t',sname
          print '\t',template['hsyst_up'][sname]
          print '\t',template['hsyst_down'][sname]
          if template['fit'] and (template['sample'] not in sname) : 
            print '\t',template['ohsyst_up'][sname]
            print '\t',template['ohsyst_down'][sname]
  
  return 

# end of: printDebugHists()
# ==========================================
# Making histogram container for various outputs
# returns : histogram list

def makeHistoContainer(data={},templates={},cc=False):

  histos = []

  # --- data ---
  if data :
    if 'hist' in data : histos.append(data['hist'])
    if 'shists' in data : histos += data['shists']
  if cc :
    if 'ohist' in data : histos.append(data['ohist'])
    if 'oshists' in data : histos += data['oshists']
    if 'srhist' in data : histos.append(data['srhist'])
  
  # --- templates ---
  if templates :
    for template in templates : 
      if 'hist' in template : 
        histos.append(template['hist'])
      if 'ohist' in template : 
        histos.append(template['ohist'])
      
      if template['syst']:
        for syst in template['hsyst_up'] : 
          histos.append(template['hsyst_up'][syst])
          histos.append(template['hsyst_down'][syst])
    
          if template['fit'] and (template['sample'] not in syst) : 
            histos.append(template['ohsyst_up'][syst])
            histos.append(template['ohsyst_down'][syst])
      
      if template['sub'] : 
        histos.append(template['sub_hist'])
 
  return histos


# end of: makeHistoContainer()
# ==========================================

# end of file
# ==========================================
# ==========================================
