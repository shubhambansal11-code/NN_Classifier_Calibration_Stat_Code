# ==========================================
# xml card generation
#
# Original Author: Migle Stankaityte
# 24 September 2019
#
# Modifications made for Z->bb calibration
# QT PhD Project 2019-2020

# ==========================================

# import

import sys, os, re
import json
import copy
from lxml import etree

import argparse

import genxml_helpers as helpers

sys.path.insert(0, 'modelmaker')
import model_helpers

# ==========================================
# ----------------- config -----------------
# top level
o_path = 'config/Vbbj'
poi = 'mu_Z'
data_types = ['data']
qcd_fncs = ['None']
xml_binning=['32']
fit_range=['[50,240]']
clist = ['']
wsfolder = 'Vbbj'
tag = 'test'
qcd_float = ',0,1000000'
def_qcd_sy = '40000'
qcd_sy = []
# ------------- getting args ---------------

parser = argparse.ArgumentParser(description='Making xml cards for xmlaAnaWSBuilder')
parser.add_argument('inputs',   nargs='+',        help='Fit Region Inputs')
parser.add_argument('--title',  required = True,  help='Title of Fit')

parser.add_argument('--data',   nargs='+',    default=data_types,   choices =['data','asimov','pois'],  help='Data Types - data/asimov/pois')
parser.add_argument('--qcd',    nargs='+',    default=qcd_fncs,     help='QCD functions - None for no func.')
parser.add_argument('--fr',     nargs='+',    default=fit_range,    help='List of fit ranges eg: \'[70,230]\'')
parser.add_argument('--bins',   nargs='+',    default=xml_binning,  help='List of number of output bins, e.g.: 32')
parser.add_argument('--qcdsy',  nargs='+',    default=qcd_sy,       help='QCD starting yield list')
parser.add_argument('-c','--cand', nargs='+', default = clist,      help='List of candidate regions: \'l\' \'s\' \'ttbar\'')
                                                                              # only needed if multiple regions
                                                                              # can use any identifying string
                                                                              # will be appended to obs_x_channel!
parser.add_argument('--tag',        default=tag,        help='Tag')
parser.add_argument('--poi',        default=poi,        help='POIs e.g.: mu,mu_V')
parser.add_argument('--wsfolder',   default=wsfolder,   help='XMLAanWSBuilder output folder')
parser.add_argument('--output',     default=o_path,     help='xml card output folder')
parser.add_argument('--outtag',     default='',         help='Additional tag only on workspace output file')
parser.add_argument('--oldyield',   help='Add old yield systematics', action='store_true')
parser.add_argument('--fixt',       help='Fix Template Yields',       action='store_true')
parser.add_argument('--nominal',    help='Nominal',                   action='store_true')
parser.add_argument('-g','--gaus',  help='ttbar gaussian',            action='store_true')
parser.add_argument('--debug',      help='debug',                     action='store_true')
args = parser.parse_args()

debug = args.debug

title = args.title
tag = args.tag
outtag = args.outtag
fit_range = args.fr
poi = args.poi
fixt = args.fixt
xml_binning = args.bins
data_types = args.data
qcd_fncs = args.qcd
wsfolder = args.wsfolder
input_list = args.inputs
gaus = args.gaus
nominal = args.nominal
clist = args.cand

# no yield systs
oldyield = args.oldyield
if nominal : oldyield = False

# region list
# removing input tags for region list
region_list = []
print 'Regions given:'
for i in input_list:
  r,t = i.split('__')
  region_list.append(r)
  print '\t'+r

# checking lengths
nregions = len(region_list)
assert (nregions == len(data_types)),"Numer of regions and data types does not match!"
assert (nregions == len(qcd_fncs)),"Numer of regions and qcd func. indications does not match!"
assert (nregions == len(xml_binning)),"Numer of regions and binning inputs does not match!"
assert (nregions == len(fit_range)),"Numer of regions and fit ranges does not match!"
assert (nregions == len(clist)),"Numer of regions and given candidates does not match!"

# overwrite clist if only one region (no need as only used to distinguis obs_x_channels)
if nregions == 1 : clist = ['']

# data types (allowed/default values set in function)
#data_types = helpers.setRegAttributes(args.data,region_list,'data')

# qcd functions
#qcd_fncs = helpers.setRegAttributes(args.qcd,region_list,'qcd')

# --- qcd yields ---
qcd_sy = args.qcdsy
nqcds=len([i for i in qcd_fncs if (i != 'None')])

# if no custom staring yields, put default in
if not qcd_sy : 
  for i in range(nqcds): qcd_sy.append(def_qcd_sy)
# check # of custom yields
else :
  assert (nqcds == len(qcd_sy)),"Numer of qcd functions requested and init yields and does not match!"

# set qcd yields
qcd_yield = []
i = 0
for f in qcd_fncs :
  if f == 'None' : qcd_yield.append(None)
  else : 
    qcd_yield.append(str(qcd_sy[i])+qcd_float)
    i += 1

if debug : print "QCD yields: ",qcd_yield

# --------------- greeting -----------------
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'Creating xml cards for '+title+' fit\n' 
# --------------- attributes ---------------
# attributes
t_att_in = 'genxml/att/top_att.json'
cat_att_in = 'genxml/att/cat_att.json'
mod_att_in = 'genxml/att/mod_att.json'
qcd_att_in = 'genxml/att/qcd_att.json'

print '\n*** Reading card attributes ***'
t_att = helpers.readInfo(t_att_in)
cat_att = helpers.readInfo(cat_att_in)
mod_att = helpers.readInfo(mod_att_in)
qcd_att = helpers.readInfo(qcd_att_in)

# ------------- reading regions ------------
print '*** Reading region info ***'

regions = {}
runs = 1 # number of runs

# slicing
sliced = 0
sliced_reg = None 
sl_dtype = None

for i,r in enumerate(region_list) : 
  print 'Reading region:',r
  regions[r] = helpers.readInfo('genxml/models/'+input_list[i]+'.json')
  
  # adding selected data type
  assert (data_types[i] in regions[r]), 'No '+data_types[i]+' in genxml/models/'+r+'.json'
  regions[r]['data_type'] = data_types[i]
  print '\tData type: ',data_types[i]
  
  # adding selected qcd type
  if not (qcd_fncs[i] == 'None') : 
    regions[r]['qcd'] = qcd_fncs[i]
    print '\tQCD function: ',qcd_fncs[i]
    # checking qcd functions 
    assert( regions[r]['qcd'] in qcd_att ),'No '+regions[r]['qcd']+' function defined in '+qcd_att_in+'!'

  # checking if multiple slice slice data files are given
  # if so - create multiple xml cards
  # Note: only works with one-fold slicing (i.e. only one sliced region)

  # convert to list if not (for old files)
  if type(regions[r][data_types[i]]['dataname']) is not list :
    regions[r][data_types[i]]['dataname'] = [regions[r][data_types[i]]['dataname']]
  
  if len(regions[r][data_types[i]]['dataname']) > 1 :
    sliced += 1
    # over-slice warning
    assert(sliced in [0,1]),'Runing over multiple slices is only allowed for one region at a time, you have given '+sliced+' sliced regions!'  
    # sliced region becomes main
    sliced_reg = r
    # number of runs
    runs = len(regions[r][data_types[i]]['dataname'])
    # slice datatype
    sl_dtype = data_types[i]

# slice debug
if debug :
  print '********** '
  print 'Sliced region: ',sliced_reg
  print 'Number of slices: ',runs

# ------------------------------------------
# *********** RUNNNG CREATION **************

# making output folder
model_helpers.mkdir('./'+o_path)
# ---

print '\n*** RUNNING CREATTION ***'

for ir in range(runs):

  s_tag = '' # slice tag
  print '.......................................'
  # slice-specific tags
  if sliced :
    slicename = regions[sliced_reg][sl_dtype]['dataname'][ir]
    nslice = slicename.split('_of',1)[0].split('_s',1)[1]
    s_tag='_s'+nslice
    #greeting
    print 'Producing xml card set #',ir+1,' using ',slicename   
        
  # ----------------- output -----------------

  # output paths
  if outtag != '' : outtag = '_'+outtag
  xmloutput = 'workspace/'+wsfolder+'/'+title+'/'+title+s_tag+'_model_'+tag+outtag
  t_path = './'+o_path+'/'+title+s_tag+'_'+tag
  m_path = t_path+'/model'
  
  if debug :
    print 'xml output: ',xmloutput
    print 'xml top cards: ',t_path
    print 'xml model cards: ',m_path

  # making folders
  print '\n*** Making card output dirrectory ***'
  model_helpers.mkdir(t_path)
  model_helpers.mkdir(m_path)
  
  # creating symbolic links
  print '*** Creating links to AnaWSBuilder.dtd ***'
  helpers.mkLink(os.getcwd()+'/dtd',t_path,'AnaWSBuilder.dtd')
  helpers.mkLink(os.getcwd()+'/dtd',m_path,'AnaWSBuilder.dtd')
  
  # --- creating card paths ---
  # top card
  top_path = t_path+'/'+title+'.xml'
  
  # adding xml cards for each region / template
  for r in regions : 
    regions[r]['cat_path'] = t_path+'/'+r+'_category.xml'
    for t in regions[r]['templates'] : 
      t['xmlcard'] = m_path+'/'+r+'_'+t['sample']+'.xml'
    if 'qcd' in regions[r] : 
      regions[r]['qcd_path'] = m_path+'/'+r+'_QCD.xml'

  # ------------------------------------------
  # -------------- TOP LEVEL -----------------
  
  # creating top-level xml card
  print '-----------------------------------'
  print 'Creating TOP-LEVEL xml card..'
  
  # top element
  t_name = 'Combination'
  t_att[t_name]['OutputFile'] = xmloutput
  comb = etree.Element(t_name,t_att[t_name])
  
  # category input file
  for r in regions:
    t_region = etree.SubElement(comb, 'Input')
    t_region.text = regions[r]['cat_path']
  
  # poi's
  t_poi = etree.SubElement(comb, 'POI')
  t_poi.text = poi
  
  # fit description
  t_att['Asimov']['Name'] = title+'_fit'
  t_fit = etree.SubElement(comb, 'Asimov',t_att['Asimov'])
  
  # printing to file
  helpers.printToFile(comb,t_name,top_path)

  # --------------------------------
  
  for ireg,r in enumerate(regions) :
    
    # creating category xml card
    print '-----------------------------------'
    print 'Creating '+r+' CATEGORY xml card..'
    
    # getting attributes
    catt = copy.deepcopy(cat_att)
  
    # top element
    cat_name = 'Channel'
    catt[cat_name]['Name'] = r+' channel'
    chan = etree.Element(cat_name,catt[cat_name])
   
    # data
    cd = etree.Comment(' === Data === ')
    chan.insert(1, cd)
    data_type = regions[r]['data_type']
    catt['Data']['InputFile']=regions[r][data_type]['outpath']
    catt['Data']['HistName']=regions[r][data_type]['dataname'][ir]
    catt['Data']['Binning']=xml_binning[ireg]
    
    cand = ''
    if clist[ireg] != '' : cand = '_'+clist[ireg]
    catt['Data']['Observable']+=cand+fit_range[ireg]
    
    cat_data = etree.SubElement(chan,'Data',catt['Data'])
  
    # qcd
    if 'qcd' in regions[r] :
      fn = regions[r]['qcd']
      cq = etree.Comment(' === QCD ('+fn+') === ')  
      chan.insert(3, cq)
      cqcd_att = catt['QCD']
      cqcd_att['InputFile'] = regions[r]['qcd_path']
      cqcd = etree.SubElement(chan,"Sample",cqcd_att)
      # yield
      assert(qcd_yield[ireg]),'Something wrong! No qcd yield in region??'
      cqcd_yield = etree.SubElement(cqcd,"NormFactor",
          {"Name":'yield_QCD['+qcd_yield[ireg]+']'})
  
  
    # templates
    for i,t in enumerate(regions[r]['templates']) :
      name = t['sample']
      ct = etree.Comment(' === '+name+' === ')
      chan.insert((i*2)+5, ct)
      ctemp_att = catt['Template']
      ctemp_att["Name"] = name
      ctemp_att["InputFile"] = t['xmlcard']
      ctemp = etree.SubElement(chan,"Sample",ctemp_att)
      
      # systematics
      
      # old yield-systematics
      if oldyield : 
        helpers.addYieldSyst(name,r,ctemp,catt['YS'])
      
      # >> from templates
      if not nominal :
        for s in t['syst']:
          # comment out as new ttbar alternative generators will account for
	  # acceptance (thus, norm) variations and de-correlated syst variations
          # will contain the process name as a suffix
          #if name in s : continue # skip alt. gen.
          # get rid of shape-pruning suffix in category-level xml norm variation
          s_new = s
          catt['TS']['Name']='alpha_'+s_new.replace('_pruned_shape','')
          if t['syst'][s][0]==0 and t['syst'][s][1]==0:
            print 'Norm variation from '+s+' pruned or null for '+name
          else:
            catt['TS']['Mag']=helpers.getMagStr(s,name,t['syst'][s][0],t['syst'][s][1])
            csyst = etree.SubElement(ctemp,"Systematic",catt['TS'])
  
      # yield
      ctemp_yield = etree.SubElement(ctemp,"NormFactor",
          {"Name":'yield_'+name+'['+str(t['yield'])+']'})
     
      # for old files
      if 'spurious' not in t : t['spurious'] = False
      if 'fixed' not in t : t['fixed'] = False
      
      # other norm factors
      if ('mu' in t ) and not t['fixed'] and not (fixt and not spurious):
        if name == "Higgs" : mu_str = 'mu'+t['mu']
        if name == "Vboson" or name == "Zboson" or name == "Wboson" : mu_str = 'mu_'+name[0]+t['mu'] 
        else : mu_str = 'mu_'+name+t['mu']
        norm = etree.SubElement(ctemp,"NormFactor",{'Name':mu_str})
      else :
        if ('mu' in t ) and t['fixed'] : print 'WARNING: mu given but template is fixed'
        helpers.addDefaultNormFactors(name,r,ctemp,data_type,fixt,t['spurious'],t['fixed'],gaus)
  
    # printing to file
    helpers.printToFile(chan,cat_name,regions[r]['cat_path'])
  
    # --------------------------------
    # creating model xml cards for each template
    for i,t in enumerate(regions[r]['templates']) :
      name = t['sample']
      print 'Creating '+r+'_'+name+' model xml card..'
  
      # getting attributes
      matt = copy.deepcopy(mod_att)
      
      # top element
      mod_name = 'Model'
      matt[mod_name]['Input'] = t['outpath']
      matt[mod_name]['ModelName'] = name+matt[mod_name]['ModelName']
      matt[mod_name]['ObservableName'] = matt[mod_name]['ObservableName']+name
      mod = etree.Element(mod_name,matt[mod_name])
  
      # luminosity settings
      it = etree.SubElement(mod,"Item",matt['Item'])
      rn = etree.SubElement(mod,"Rename",matt['Rename'])
     
      # systematics
      if not nominal :
        for s in t['syst']:
          # drop pruned shape variations
          if '_pruned_shape' in s: continue
          ms_att = {}
          ms_att['ConstrName']='alpha_'+s+'Constraint'
          ms_att['NPName']='alpha_'+s
          ms_att['GOName']='nom_alpha_'+s
          msyst = etree.SubElement(mod,"ExtSyst",ms_att)
  
        # loop over gammas
        for s in t['gammas']:
          s = str(s)
          ms_att = {}
          ms_att['ConstrName']='gamma_stat_'+name+'_bin_'+s+'_constraint'
          ms_att['NPName']='gamma_stat_'+name+'_bin_'+s
          ms_att['GOName']='nom_gamma_stat_'+name+'_bin_'+s
          msyst = etree.SubElement(mod,"ExtSyst",ms_att)
  
      # printing to file
      helpers.printToFile(mod,mod_name,t['xmlcard'])
    # --------------------------------
    # creating model xml card for QCD
    
    if 'qcd' in regions[r] :
      fn = regions[r]['qcd']
      print 'Creating '+r+'_QCD model xml card using '+fn+' function..'
      qcd_name = 'Model'
      qcd = etree.Element(qcd_name,qcd_att[qcd_name])
      qcd_fn = etree.SubElement(qcd,'ModelItem',qcd_att[fn])
      
      # printing to file
      helpers.printToFile(qcd,qcd_name,regions[r]['qcd_path'])
  
  # --------------------------------
