# ==========================================
# Simple Automatic XML Model Creator
#
# Migle Stankaityte
# 18 March 2019
#
# changes for QT PhD project Z->bb calibration
# ==========================================

# import
import sys, os, re
import argparse
import copy
import json

import model_helpers as helpers
import data_helpers
import template_helpers
#import data_helpers
#import template_helpers
import sigfit
import syst_smooth

from  generate_asimov import *

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TFile, TH1F, RooStats

# ==========================================
# Start of procedure

# --------------------------------
# getting argumets

parser = argparse.ArgumentParser(description='Making HistFactory Model workspace for xmlAnaWSBuilder')
parser.add_argument('jpath', help='Path of JSON file with model info')
parser.add_argument('gevbins', help='GeV per bin') #NOTE:  now input is GeV and NOT # of bins!
parser.add_argument('model_name', help='Model title')
parser.add_argument('tag', help='Output tag')
parser.add_argument('-c','--cand', choices=['l','s',''], default = '', help='Candidate region: leading (l) / subleading (s)')
args = parser.parse_args()

# --------------------------------
# config

tree = 'outTree' 

# --- binning ---
m_fitmin = 50.
m_fitmax = 150.
m_min = 50.
m_max = 150.
pt_min = 450.
pt_max = 3000.
gevbins_fit = 2.5 #5 GeV always
gevbins = args.gevbins
#gevbins = 2.5

# --- threshold for gamma parameters ---
gammas_thr = 0.3

# --- thresholds for pruning --- 
prun_norm_thr = -1
prun_shape_thr = -1

# --------------------------------
# salutation
print ' ~~~~~~~~~~ XML Workspace Creator ~~~~~~~~~ \n'

# --------------------------------
# reading model info from json
data, templates, flags = helpers.readModelInfo(args.jpath)

# --------------------------------

# getting alternative range from json
if 'm_min' in data : m_min = float(data['m_min'])
if 'm_max' in data : m_max = float(data['m_max'])
if 'pt_min' in data : pt_min = float(data['pt_min'])
if 'pt_max' in data : pt_max = float(data['pt_max'])
if 'm_fitmin' in data : m_fitmin = float(data['m_fitmin'])
if 'm_fitmax' in data : m_fitmax = float(data['m_fitmax'])
if 'gammas_thr' in data : gammas_thr = float(data['gammas_thr'])
if 'prun_norm_thr' in data : prun_norm_thr = float(data['prun_norm_thr'])
if 'prun_shape_thr' in data : prun_shape_thr = float(data['prun_shape_thr'])

# setting number of bins
nbins_fit = int(round((m_fitmax-m_fitmin)/gevbins_fit))
nbins = int(round((m_max-m_min)/float(args.gevbins))) 
#nbins = int(round((m_max-m_min)/float(gevbins))) 
# --- output tags and paths ---
cand = args.cand
tag = args.tag+'_'+cand+'%d'%(nbins)
fittag = args.tag+'_'+cand+'%d'%(nbins_fit)
model_name = args.model_name
ws_pfx = 'workspaces/'+model_name+'/'+tag
ccpath = 'crosscheck/'+model_name+'_'+tag+'/'
ccname = model_name+'_'+tag
injtest_out = 'xmlAnaWSBuilder/Input/data/sig_inj/'+model_name+'_'+tag+'_injtest'
if flags['mkdatain'] : data_out = 'data_histos/'+model_name+'_'+tag+'_input_data'
else : data_out = 'xmlAnaWSBuilder/Input/data/histo/'+model_name+'_'+tag+'_data'

# --- json output for xml creation ---
if flags['mk_xmljson'] : 
  if cand : 
    jsonout = 'genxml/models/'+model_name+'_'+cand+'__'+args.tag+'.json'
  else : 
    jsonout = 'genxml/models/'+model_name+'__'+args.tag+'.json'

# --------------------------------
# region for plots 
reg = model_name
if cand == 'l' : reg = "Leading "+reg
elif cand == 's' : reg = "Subleading "+reg
if ('pt_min' in data) and ('pt_max' in data) : reg = reg+'p#_{T} ['+str(data['pt_min'])+'; '+str(data['pt_max'])+'] GeV'
print "Region: "+reg+'\n'

# --------------------------------
# print to screen
helpers.printOutPaths(flags,data_out,ws_pfx,model_name,ccpath+ccname,injtest_out)

print '\n----------------------------------' 

# creating folders
if not flags['slice_calc']:
  if not flags['no_ws'] : helpers.mkdir('./workspaces')
  if flags['mkdatain'] : helpers.mkdir('./data_histos')

if not flags['no_cc'] :
  helpers.mkdir('./crosscheck')
  helpers.mkdir(ccpath)

if flags['mk_xmljson'] : 
  helpers.mkdir('genxml/models/')


# --------------------------------
# smoothing binning
rebin = int(round(float(gevbins_fit)/float(gevbins)))
#print rebin

# ==========================================
# making and filling histograms

# ------ data -----
data_helpers.getDataHistos(data,flags,tree,nbins,m_min,m_max,pt_min,pt_max,flags['slice_calc'],cand)

# ----- templates -----

# init lists
# combination templates 
combtemps = []
# templates with systematics from other templates
nomtemps = []

# looping over templates
for template in templates:
  
  # set names
  temp_name = str(template['sample'])
  temp_path = str(template['path'])

  # in slice-calc mode ignore all templates which are not subtracted from data
  if flags['slice_calc'] and (temp_name not in data['sub_t']) :
    print 'Slice calculation mode: template ',temp_name,' not subtracted, therefore ignored.'
    continue

  print '\n----------------------------------' 
  print 'Filling '+temp_name+' histogram:'
  
  # adding candidate to region
  if template['region'] and cand :
    template['region']+= cand
    print 'Getting histograms for '+temp_name+' '+template['region']+':'
  else :
    print 'Getting histograms for '+temp_name+':'

  # --- getting histograms ---
  
  # procedure if the template is fit with parametric functions:
    # 1. Get histograms from input file with fit binning/range (add '_orig' to hist names)
    # 2. Run fitting procedure, output new histograms in output binning/range 
  if template['fit']: 
    
    # get histograms (fit range)
    template['ohist'], template['ohsyst_up'], template['ohsyst_down'] = template_helpers.getTemplateHistos(temp_name,
                                                                                                           template['syst'], template['weightlist'],
                                                                                                           template['Vjets_wsyst_up'], template['Vjets_wsyst_down'],
                                                                                                           temp_path,
                                                                                                           tree,
                                                                                                           nbins_fit,m_fitmin,m_fitmax,
                                                                                                           pt_min,pt_max,
                                                                                                           template['truth_pt_min'], template['truth_pt_max'],
                                                                                                           template['w_min'],template['w_max'],
                                                                                                           template['dsids'],
                                                                                                           template['region'],
                                                                                                           template['truth_reweighting'],
                                                                                                           flags['conf'], flags['sepuncerts'],
                                                                                                           fit=True,
                                                                                                           basicsymm=flags['basicsymm'],
                                                                                                           updown=flags['updown']) 
    # parametric fit template histograms
    fitpath = ccpath+'/fit_'+temp_name+'/'
    helpers.mkdir(fitpath)

    template['hist'],template['hsyst_up'],template['hsyst_down'] = sigfit.getFittedHistos(temp_name,
                                                                                          template['ohist'],
                                                                                          template['ohsyst_up'],
                                                                                          template['ohsyst_down'],
                                                                                          gevbins_fit,nbins_fit,m_fitmin,m_fitmax,
                                                                                          nbins,m_min,m_max,
                                                                                          fitpath,
                                                                                          flags['conf'],
                                                                                          cand)

  # if template not fit make histograms in output range/binning (no '_orig' extention to hist names)
  else :
    # get histograms (output range)
    template['hist'],template['hsyst_up'],template['hsyst_down'] = template_helpers.getTemplateHistos(temp_name,
                                                                                                      template['syst'], template['weightlist'],
                                                                                                      template['Vjets_wsyst_up'], template['Vjets_wsyst_down'],
                                                                                                      temp_path,
                                                                                                      tree,
                                                                                                      nbins,m_min,m_max,
                                                                                                      pt_min,pt_max,
                                                                                                      template['truth_pt_min'], template['truth_pt_max'],
                                                                                                      template['w_min'],template['w_max'],
                                                                                                      template['dsids'],
                                                                                                      template['region'],
                                                                                                      template['truth_reweighting'],
                                                                                                      flags['conf'], flags['sepuncerts'],
                                                                                                      basicsymm=flags['basicsymm'],
                                                                                                      updown=flags['updown'])

  if template['smoothsyst'] :
    # get histograms (fit range)
    template['whist'], template['whsyst_up'], template['whsyst_down'] = template_helpers.getTemplateHistos(temp_name,
                                                                                                           template['syst'], template['weightlist'],
                                                                                                           [], [],
                                                                                                           temp_path,
                                                                                                           tree,
                                                                                                           nbins_fit,m_fitmin,m_fitmax,
                                                                                                           pt_min,pt_max,
                                                                                                           template['truth_pt_min'], template['truth_pt_max'],
                                                                                                           template['w_min'],template['w_max'],
                                                                                                           template['dsids'],
                                                                                                           template['region'],
                                                                                                           template['truth_reweighting'],
                                                                                                           flags['conf'], flags['sepuncerts'],
                                                                                                           smooth=True,
                                                                                                           basicsymm=flags['basicsymm'],
                                                                                                           updown=flags['updown']) 
    if not template['combine'] :
      # smooth systematics (combination smoothing later, in combTemps())
      smoothpath = ccpath+'/smoothsysts_'+temp_name+'/'
      helpers.mkdir(smoothpath)
    
      syst_smooth.getSmoothedSysts(template,rebin,smoothpath,
                                  nbins,m_min,m_max,
                                  nbins_fit,m_fitmin,m_fitmax,reg)

    
  # --- subtraction templates ---     
  # storing histogram for subtraction from data 
  if template['sub'] :
    # if template is subtracted from data prefit - need nominal with correct binning
    if template['fit'] and not template['sub_fit'] : 
      hsub_dummy = template_helpers.getTemplateHistos(temp_name,
                                                      [], [],
                                                      [], [],
                                                      temp_path,
                                                      tree,
                                                      nbins,m_min,m_max,
                                                      pt_min,pt_max,
                                                      template['truth_pt_min'], template['truth_pt_max'],
                                                      template['w_min'],template['w_max'],
                                                      template['dsids'],
                                                      template['region'],
                                                      template['truth_reweighting'],
                                                      flags['conf'], flags['sepuncerts'],
                                                      tag='orig_rebin')[0]
      template['sub_hist'] = helpers.smearGauss(hsub_dummy,'smeared')
    else : template['sub_hist'] = helpers.smearGauss(template['hist'],'smeared')
    # if data is sliced - rescale subtraction histogram to slice size
    if 'oshists' in data : template['sub_hist'].Scale(1/float(data['n_slice']))  
 
  # --- combination and syst templates ---
  # ignore the rest if slice calculation mode
  if flags['slice_calc'] : continue
 
  # store all combination templates
  if template['combine'] : 
    for c in template['combine'] :
      if c not in combtemps : combtemps.append(c)
  
  # store all templates with shape syst. from other templates
  if template['syst_for'] : 
    for nt in template['syst_for'] :
      if nt not in nomtemps : nomtemps.append(nt)


# ==========================================
# DEBUG: check if all histograms are filled
#helpers.printDebugHists(data,templates)

# ==========================================
# Template subtraction from data histogram(s)

if data['sub'] :
  data_helpers.subTempFromDataHist(data,templates)

  # exiting 
  print '\nTemplate subtraction is done!'

# ==========================================
# Data slice calculation

#slice info
s_inf = {}

#runing procedure
if flags['slice_calc'] :

  # slice calculation
  s_inf = data_helpers.dataSliceCalc(data,tree,flags['conf']) 
  
  # exiting 
  print '\nSlice calculation is done!'

# ==========================================
# adding shape-systematic templates 

if nomtemps : 
  template_helpers.addTemplateSysts(templates,nomtemps,basicsymm=flags['basicsymm'],updown=flags['updown'])

# ==========================================
# making combination templates

if combtemps : 
  template_helpers.makeCombTemp(templates,combtemps,
                                flags['nominal'],
                                nbins,m_min,m_max,
                                nbins_fit,m_fitmin,m_fitmax,
                                ccpath,reg,rebin,flags['sepuncerts'])

print '--------------------------------------'

# ==========================================
# asimov data sets

# making containers
asimov = {}
pois = {}

# creating histograms
Asimov_data = TH1F("Asimov_data","Asimov_data",nbins,m_min,m_max)
# poisson sampled data sets
SAMPLED_Poiss_data = TH1F("SAMPLED_Poiss_data","SAMPLED_Poiss_data",nbins,m_min,m_max)
QCD = None
# filling histograms
if flags['injection_test'] == True:
  # data[injection] flag:
  # if flag is null -> no injection of QCD is performed
  # if is a list of parameters -> parameteric injection is performed,
  # if is a number -> the corresponding subtracted slice is injected
  if not data['injection']:
    print("NO QCD inj")
    pass
  elif type(data['injection']) == list:
    QCD = generate_asimov_QCD(nbins,m_min,m_max,data['injection'])
    Asimov_data.Add(QCD,1)
  else: 
    QCD = data['shists'][int(data['injection'])]
    Asimov_data.Add(QCD,1)

  print "\n Performing injection tests:"
  
  for template in templates:
    print("temp ",template['sample'],template['injection'])
    if not template['injection'] : continue
    print 'Adding ',template['sample'],' to Asimov dataset'
    print 'Integral',template['hist'].Integral()*template['injection']
    Asimov_data.Add(template['hist'],template['injection'])
    print(Asimov_data.Integral())
    

  trd = TRandom3(0)
  print("random seed ",trd)
  for _bin in range(1,nbins+1):
    eventsGenerated = Asimov_data.GetBinContent(_bin)
    Asimov_data.SetBinError(_bin,pow(eventsGenerated,0.5))
    # SAMPLED version of Asimov_data
    sampled_eventsGenerated = trd.Poisson(eventsGenerated)
    SAMPLED_Poiss_data.SetBinContent(_bin,sampled_eventsGenerated)
    SAMPLED_Poiss_data.SetBinError(_bin,pow(sampled_eventsGenerated,0.5))

  # filling containers with output info
  asimov['dataname'] = Asimov_data.GetTitle()
  pois['dataname'] = SAMPLED_Poiss_data.GetTitle()

# ==========================================
# making HistFactory model

# saving yields
yields = []
y_names = []

# template workspaces
ws_temps = []

print '\n===================================================='
if not flags['no_ws'] :

  print 'Making HistFactory Models\n' 
   
  # looping over templates
  # note: if you want to use the histogram later 
  #       you cannot pass it dirrectly to the HistFactory
  #       as it destroys it, hence make clones
  
  for template in templates: 
  
    # skiping combined sub-templates
    if template['combine'] : continue
    # skiping ingnored templates
    if template['ignoreWS'] : continue
 
    # names
    sname = str(template['sample'])
    
    # log for path output
    ws_temps.append(sname)

    print '\n----------------------------------'
    print 'Making HistFactory Model:\t',sname
  
    # create measurement
    meas = RooStats.HistFactory.Measurement(model_name,model_name)
    meas.SetOutputFilePrefix(ws_pfx+'_'+sname);
   # meas.AddPOI('mu');
    meas.AddPOI('mu_Zboson');
    meas.SetExportOnly(1);
  
    # scale histogram content, which already includes lumi, so set to 1
    meas.SetLumi(1.0);
    
    # making channel and sample
    # NEW: need different channel name per
    # each process, in order to have gamma
    # parameters with a different name
    tch = RooStats.HistFactory.Channel(sname)
    tsamp = RooStats.HistFactory.Sample(sname)
    # when avoiding parametric smoothing of 
    # templates -> activate gamma parameters
    if not template['fit'] and not flags['nominal'] :
      template['gammas'] = template_helpers.gammaBins(template['hist'],gammas_thr)
      tch.SetStatErrorConfig(gammas_thr, "Poisson")
      tsamp.ActivateStatError()
    else :
      template['gammas'] = {}
    
    # getting nominal histo
    thist = template['hist'].Clone(sname)
    
    # getting yield
    nom_yield = thist.Integral() if thist.Integral() else 1
    if flags['slicetemp'] and ('shists' in data) and not template['spurious']: 
      yields.append(nom_yield/float(data['n_slice']))
      template['yield'] = nom_yield/float(data['n_slice'])
    else:
      yields.append(nom_yield)
      template['yield'] = nom_yield
    y_names.append(sname)
    
    # normalise template
    if flags['norm_all'] :
      thist.Scale(1/nom_yield)
    elif flags['slicetemp'] and ('shists' in data): # sliced but not normalised
      thist.Scale(1/float(data['n_slice']))
      sys.exit()
    
    # setting histo
    tsamp.SetHisto(thist)
  
    #----------------
    # adding systematics
    template['syst_diffs'] = {}
    for syst in template['hsyst_up'] :
      print str(syst)
      syst_name = syst

      #decorrelate JET systematics? 
      if 'decorrelate_JET' in template: 
        if (template['decorrelate_JET']) and ('JET_CombMass_' in str(syst)) : 
          print '\tdecorrelate_JET=true, going to decorrelate ',str(syst),' for ',str(template['sample']),'\n'
          syst_name += '_'+template['sample']

      shape_syst = RooStats.HistFactory.HistoSys( str(syst_name) )
     
      hs_up = template['hsyst_up'][syst].Clone(sname+str(syst)+'__1up')
      hs_down = template['hsyst_down'][syst].Clone(sname+str(syst)+'__1down')
 
      # geting magnitude of difference
      diff_up = (hs_up.Integral()-nom_yield)/nom_yield
      diff_down = (hs_down.Integral()-nom_yield)/nom_yield
      y_names.append(sname+str(syst)+'__1up')
      y_names.append(sname+str(syst)+'__1down')
      if max(abs(diff_up),abs(diff_down)) < prun_norm_thr: 
        diff_up = 0
        diff_down = 0
        print '\tPruned norm variation from ',str(syst),' on ',str(template['sample']),'\n'
      yields.append(diff_down)
      yields.append(diff_up)

      # normalise syst templates
      if flags['norm_systs'] :
        hs_up.Scale(nom_yield/hs_up.Integral())
        hs_down.Scale(nom_yield/hs_down.Integral())

        # N.B. norm_all is the right thing to do
      elif flags['norm_all'] :
        if hs_up.Integral():
          hs_up.Scale(1/hs_up.Integral())
        if hs_down.Integral():
          hs_down.Scale(1/hs_down.Integral())
      elif flags['slicetemp'] and ('shists' in data): 
        hs_up.Scale(1/float(data['n_slice']))
        hs_down.Scale(1/float(data['n_slice']))
        sys.exit()
      
      # shape pruning 
      if template_helpers.maxShapeVariation(template['hist'],hs_up,hs_down) > prun_shape_thr : 
        shape_syst.SetHistoHigh(hs_up)
        shape_syst.SetHistoLow(hs_down)
    
        tsamp.AddHistoSys(shape_syst);
  
        # prevents attempts of double-destruction
        ROOT.SetOwnership(shape_syst,ROOT.kFALSE)
      else :
        print '\tPruned shape variation from ',str(syst),' on ',str(template['sample']),'\n'
        # append keyword to drop pdf variation in xml-card creation
        syst_name = syst_name+'_pruned_shape'

      template['syst_diffs'][str(syst_name)] = [diff_up,diff_down]


    #----------------
    # adding channel to sample
    
    tch.AddSample( tsamp )
    tsamp.Print()
    print '\t'+sname+' sample added to channel.\n'
    
    meas.AddChannel( tch )    
    print '\tChannel added to measurement.\n'
     
    #----------------
    # prevents attempts of double-destruction
    ROOT.SetOwnership(tsamp,ROOT.kFALSE)
    ROOT.SetOwnership(tch,ROOT.kFALSE)

    #----------------
    # creating the workspace

    print 'Creating workspace\n'
    ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast( meas )
    print 'Workspace created.\n'
    
    print 'Printing WS Model Info:'
    meas.PrintTree()
    
    # prevents attempts of double-destruction
    ROOT.SetOwnership(meas,ROOT.kFALSE)

else : 
  print '\nNo HistFactory model made! ',
  print '\n\'no_ws\' flag on in .json file.'

print '\n===================================================='

# ==========================================
# storing data histogram

if not flags['no_dataout'] :
  print 'Making data histogram file.'
  
  datahistos = helpers.makeHistoContainer(data=data)
  helpers.makeHistoRoot(data_out,datahistos)

# ==========================================
# making cross-check output
print 'Making sanity check output.'

if not flags['no_cc'] :
  
  histos = helpers.makeHistoContainer(data=data,
                                      templates=templates,
                                      cc=True)

  # adding asimov 
  if flags['injection_test']:
    histos.append(Asimov_data)
    if QCD == None:
      print("No QCD in injection data")
    else:
      histos.append(QCD)

  helpers.makeHistoRoot(ccpath+ccname,histos)

# ==========================================
# making signal injection histos: Asimov and Poiss
siginj_histos = []
if flags['injection_test'] == True:
  siginj_histos.append(Asimov_data)
  siginj_histos.append(SAMPLED_Poiss_data)
  helpers.makeHistoRoot(injtest_out,siginj_histos)

# ==========================================
# creating json for xml card creation

if flags['mk_xmljson'] : 
  
  # adding output info to containers
  helpers.alocatePaths(flags,data,templates,asimov,pois,
                       data_out,ws_pfx,model_name,injtest_out)

  # making new containers for printing
  p_glob = {}
  p_data = {}
  p_temp = []

  # data
  if not flags['no_dataout'] :
    p_data['outpath'] = data['outpath']
    p_data['dataname'] = data['dataname']
    p_glob['data'] = p_data
  
  # data
  if flags['injection_test'] :
    p_glob['asimov'] = asimov
    p_glob['pois'] = pois

  # templates
  for t in templates : 
    # skiping sub-templates and ignored templates
    if t['combine'] or t['ignoreWS'] : continue
    
    pt = {}
    pt['sample'] = t['sample']
    pt['yield'] = t['yield']
    pt['syst'] = t['syst_diffs']
    pt['gammas'] = t['gammas']
    pt['outpath'] = t['outpath']
    pt['spurious'] = t['spurious']
    pt['fixed'] = t['fixed']

    if t['mu'] : pt['mu'] = t['mu']

    p_temp.append(pt)
 
  p_glob['templates'] = p_temp

  # writing json
  with open(jsonout, 'w') as outfile:
    json.dump(p_glob, outfile,indent=2)

# ==========================================
# finishing

# out paths
helpers.printOutPaths(flags,data_out,ws_pfx,model_name,
                      ccpath+ccname,injtest_out,
                      templates=ws_temps)

# yields
if yields : 
  print '\nYields and Magnitudes: '
  if 'hist' in data : print '\tData:\t',data['hist'].Integral()
  if 'shists' in data : 
    for i, h in enumerate(data['shists']) : 
      print '\tData slice #',i,':\t',h.Integral()
  for i,t in enumerate(yields) : 
    print '\t',y_names[i],':\t',t

# injection yield
if flags['injection_test'] == True:
  print '\nAsimov dataset yield: ',Asimov_data.Integral()
  if QCD == None:
    print 'Asimov QCD yield: None'
  else:
    print 'Asimov QCD yield: ',QCD.Integral()

# slice number information
if s_inf : 
  print '\n Data slicing information:' 
  print '\tSR:\t',s_inf['SR']
  print '\tVR:\t',s_inf['VR']
  print '\tSR/VR:\t',s_inf['ratio']
  print '\tSlices:\t',s_inf['slices']

if flags['norm_systs'] :
  print "\nSystematic template yields were normalised to nominal."

if flags['mk_xmljson'] : 
  print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n'
  print 'Info file for xml card creation:'
  print '\t'+jsonout

print '\n~~~~~~~~~~ XML Workspace Creator ~~~~~~~~~ \n'


# end of file
# ==========================================
# ==========================================
