'''
pyHXEXPRESS v0.1.0

Copyright 2025 Lisa M Tuttle

https://github.com/tuttlelm/pyHXExpress

26Mar2025
Various additions to handle ion_type, manually specify undeut_mz

22Aug2024
add rss normalization to calc_rss()
norm fit values before saving as fitparamsdf
add max_Int to datafits
'''

#%matplotlib widget 
import numpy as np, pandas as pd 
import scipy.stats as stats
import matplotlib.pyplot as plt 
import matplotlib.lines as mlines
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages
#import fitz #!pip install pymupdf 
#import matplotlib.colors as mcolors

import os
import sys
import importlib

from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
#from scipy.stats import rankdata, skew
#from scipy.special import comb
from math import gamma, lgamma, exp
from pyteomics import mass
#!pip install brain-isotopic-distribution
from brainpy import isotopic_variants
import random
from datetime import datetime
from collections import Counter
from Bio import SeqIO

import config
#from . import config

rng=np.random.default_rng(seed=config.Random_Seed)

now = datetime.now()
date = now.strftime("%d%b%Y")


mpl_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
mpl_colors_dark = ['#005794', '#df5f0e', '#0c800c', '#b60708', '#74479d', '#6c362b', '#b357a2', '#5f5f5f', '#9c9d02', '#079eaf']
mpl_colors_light = ['#2f97e4', '#ff9f2e', '#4cc04c', '#f64748', '#b487ed', '#ac766b', '#f397e2', '#9f9f9f', '#dcdd42', '#37deef']
mpl_colors_light2 = ['#4fb7f4', '#ffbf4e', '#6ce06c', '#f66768', '#d4a7fd', '#cc968b', '#f3b7f2', '#bfbfbf', '#fcfd62', '#57feff']
colors2d = [mpl_colors_dark, mpl_colors, mpl_colors_light,mpl_colors_light2]


import warnings
warnings.simplefilter (action='ignore', category=FutureWarning)

# class dotdict(dict):
#     """dot.notation access to dictionary attributes"""
#     __getattr__ = dict.get
#     __setattr__ = dict.__setitem__
#     __delattr__ = dict.__delitem__

# obsolete function after introduction of import config 
# def get_parameters(PARAMS_FILE):
#     params_dir = os.path.split(PARAMS_FILE)
#     curr_dir = os.getcwd()
#     os.chdir(params_dir[0])
#     para_file = params_dir[1].rsplit('.')[0]

#     #p,m = PARAMS_FILE.rsplit('.',1)
#     mod = importlib.import_module(para_file)
#     importlib.reload(mod)
#     values = {v: getattr(mod, v)
#                     for v in mod.__dict__
#                     if not v.startswith("_")}
#     globals().update({v: getattr(mod, v)
#                     for v in mod.__dict__
#                     if not v.startswith("_")})
#     os.chdir(curr_dir)
#     return values

def write_parameters(write_dir=os.getcwd(),overwrite=False):
    '''
    save current user settings to a python file
    this can be read in using: 
    > import pyhxexpress as hxex
    > import file_name.py as config
    > hxex.config = config #must run to update parameters to be used
        file needs to be in same file location as the current working directory set in the notebook (this is the default save location)
    '''

    all_params = ['Allow_Overwrite', 'BestFit_of_X', 'Binomial_dCorr', 'Data_DIR', 'Data_Type', 'DiffEvo_kwargs', 'DiffEvo_threshold','Dfrac', 'Env_limit', 
                  'Env_threshold', 'FullDeut_Time', 'Hide_Figure_Output', 'Keep_Raw', 'Limit_by_envelope', 'Max_Pops', 'Metadf_File', 'Min_Pops', 'Nboot', 
                  'Ncurve_p_accept', 'Nex_Max_Scale', 'Nterm_subtract', 'Output_DIR', 'Overlay_replicates', 'Peak_Resolution', 'Pop_Thresh', 'Preset_Pops', 
                  'Preset_Pops_File', 'process_ALL', 'Random_Seed', 'Read_Spectra_List', 'Residual_Cutoff','Save_Spectra', 'Scale_Y_Values', 'setNoise', 
                  'SVG', 'Test_Data', 'Use_DiffEvo', 'User_mutants', 'User_peptides', 'WRITE_PARAMS', 'Y_ERR', 'Zero_Filling']

    filename = "hdxms_params_"+config.date+".py"
    write_file = os.path.join(write_dir,filename)
    if not overwrite:
        add_string=""
        i = 0
        while os.path.exists(write_file):
            i += 1
            add_string=str(i)
            write_file = os.path.join(write_dir,"hdxms_params_"+config.date+"_"+add_string+".py")
    
    with open(write_file,"w") as p_file:
        print_name = os.path.join(os.path.basename(os.path.normpath(write_dir)),filename)
        print("Saving config parameters to "+print_name+"\n")
        #write header
        p_file.write("import os\n"+
                     "from datetime import datetime\n"+
                     "now = datetime.now()\n"+
                     "date = now.strftime('%d%b%Y')\n\n")
        p_file.write("#config parameters "+str(config.date)+"\n")
        for p in all_params:
            # #print(p,globals().get(p))
            # if type(vars().get(p)) == str:
            #     p_file.write(p+" = \""+str(globals().get(p))+"\"\n")
            # else: p_file.write(p+" = "+str(globals().get(p))+"\n")
            try:
                if type(getattr(config,p)) == str: 
                    #p_file.write(p+" = \""+str(getattr(config,p))+"\"\n")
                    p_file.write(p+" = r\""+os.path.normpath(getattr(config,p))+"\"\n")
                else:
                    p_file.write(p+" = "+str(getattr(config,p))+"\n")
            except: 
                p_file.write("# missing parameter "+p+"\n")

## Functions to read data in different formats

def get_hxexpress_meta(hx_file):
    '''
    this was written for a specific filename format: SAMPLE_START-END-PEPTIDE-zCHARGE-usertext.xlsx
    '''
    labels = hx_file.split('-')
    metadata = {}
    metadata['file']=hx_file
    metadata['sample'] = str(labels[0]).rsplit('_',1)[0] #split off the peptide start_seq
    start_seq = str(labels[0]).rsplit('_',1)[-1] 
    end_seq = labels[-4] #end sequence is fourth from end, allows for whatever sample name formatting
    metadata['start_seq'] = int(start_seq)
    metadata['end_seq'] = int(end_seq)
    metadata['peptide_range'] = '-'.join([start_seq,end_seq])
    metadata['charge']=float(labels[-2][1:])
    metadata['peptide'] = labels[-3]
    return metadata

def save_metadf(metadf,filename=None):
    '''
    save the metadf (or other dataframe) as csv file
    '''
    savedf = metadf.copy()
    int_cols = ['charge','start_seq','end_seq','rep']
    if not os.path.exists(config.Output_DIR): os.makedirs(config.Output_DIR)
    for ic in int_cols:
        if ic in savedf.columns:
            savedf[ic] = savedf[ic].astype(int)
    if filename:        
        saveto = os.path.join(config.Output_DIR,filename)
    else: 
        saveto = os.path.join(config.Output_DIR,'hdx_spectra_list_metadf_'+date+'.csv')
    savedf.to_csv(saveto,index_label='Index',mode='w')

def read_metadf(filename):
    '''
    read in the metadf dataframe from a csv file 
    '''
    try:
        readfrom = os.path.join(config.Data_DIR,filename)
        readmetadf = pd.read_csv(readfrom).drop('Index',axis=1)
        return readmetadf
    except: 
        print("No file found for spectra list")
        return

def get_metadf(quiet=False):    
    if config.Read_Spectra_List: 
        metadf = read_metadf(config.Metadf_File)        
    else: 
        metadf = pd.DataFrame()
        if config.Data_Type == 1:
            if config.Test_Data: 
                metadf = read_metadf(config.Metadf_File)
            else:
                hx_files = [ f for f in os.listdir(config.Data_DIR) if f[-5:]=='.xlsx'  ] 
                if not config.process_ALL:
                    if 'all' not in config.User_mutants[0].lower():
                        #available_mutants = set([um for hx in hx_files for um in User_mutants if um in hx])     
                        hx_files = [hx for um in config.User_mutants for hx in hx_files if um == hx.split('-')[0].rsplit('_',1)[0]]
                    if 'all' not in config.User_peptides[0].lower():
                        #available_charges = set([up for hx in hx_files for up in User_peptides if up in hx])
                        hx_files = [hx for up in config.User_peptides for hx in hx_files if up in hx]


                #HSPB1_B1B5_0001-0011-MTERRVPFSLL-z2-allspectra.xlsx
                #metadf = pd.DataFrame() #dataframe to hold filenames and sample/peptide/charge info
                for f in hx_files:
                    meta = get_hxexpress_meta(f)
                    metadf = metadf.append(meta,ignore_index=True)

        elif config.Data_Type == 2:
            fasta_files = [ f for f in os.listdir(config.Data_DIR) if f[-6:]=='.fasta'  ]
            if len(fasta_files)==0: print("no fasta files found")
            mutants = [ff.split('.')[0] for ff in fasta_files]
            mutant_dirs = [os.path.split(f)[-1] for f in os.scandir(os.path.join(config.Data_DIR)) if f.is_dir()]
            mutants = [o for o in mutant_dirs if o in mutants]  #ignore any extra fasta files

            if (not config.process_ALL) and (config.User_mutants[0].lower() != 'all'): 
                check = all(item in mutants for item in config.User_mutants)
                if not check:             
                    missing = list(set(config.User_mutants)-set(mutants))
                    remaining = list(set(config.User_mutants)-set(missing))
                    print("missing fasta files for: ", *list(set(config.User_mutants)-set(mutants)))      
                    mutants = [o for o in config.User_mutants if o in remaining]
                    if len(mutants): print("only processing",*mutants)
                else: mutants = config.User_mutants

            smeta = {}
            #metadf = pd.DataFrame()
            for mutant in mutants:
                peptide_dirs = [f.path for f in os.scandir(os.path.join(config.Data_DIR,mutant)) if f.is_dir()]
                peptide_ids = [os.path.split(ff)[-1] for ff in peptide_dirs]  
                peptide_ranges = [ff.rsplit('-',1)[0] for ff in peptide_ids]
                if (not config.process_ALL) and (config.User_peptides[0].lower() != 'all'): 
                    check = all(item in peptide_ranges for item in config.User_peptides)
                    if not check:
                        print("missing peptides for: ",mutant, *list(set(config.User_peptides)-set(peptide_ranges)))
                        missing = list(set(config.User_peptides)-set(peptide_ranges))
                        remaining = list(set(config.User_peptides)-set(missing))
                        peptide_ranges = [o for o in config.User_peptides if o in remaining]
                    else: 
                        peptide_ranges = config.User_peptides
                    peptide_dirs = list(filter(lambda x: any(userpep in x for userpep in config.User_peptides),peptide_dirs))
                
                fasta_sequence =  SeqIO.parse(open(os.path.join(config.Data_DIR,str(mutant)+'.fasta')),'fasta')
                for fasta in fasta_sequence:
                    sequence = str(fasta.seq)

                for spec_dir,peptide_range in zip(peptide_dirs,peptide_ranges):
                    csv_files = [ f for f in os.listdir(spec_dir) if f[-4:]=='.csv'  ]
                    charges = set([int(c[-5]) for c in csv_files if c[-5].isdigit() and c[-6]=='z'])
                    for charge in charges:
                        smeta['file']=os.path.split(spec_dir)[1]
                        smeta['sample']=os.path.split(os.path.split(spec_dir)[0])[1]
                        smeta['charge']=float(charge)
                        smeta['peptide_range']=peptide_range
                        startseq = int(peptide_range.split('-')[0])
                        endseq = int(peptide_range.split('-')[1])
                        smeta['start_seq'] = startseq
                        smeta['end_seq'] = endseq
                        smeta['peptide'] = sequence[startseq-1:endseq]
                        #metadf = metadf.append(smeta, ignore_index = True)   
                        metadf = pd.concat([metadf,pd.DataFrame([smeta])],ignore_index = True)
        else: print("Must specify Data_Type: 1 for HX-Express allspectra format or 2 for SpecExport")
    if metadf.empty: raise Exception("No data found. Check Data_Type specification and Data_DIR")
    else:
        metadf = metadf.sort_values(['peptide_range','sample','charge',],ignore_index=True) 
        if not quiet:
            print("Found",len(metadf['sample'].unique()),"sample types with",len(metadf),"total datasets to analyze.")
    return metadf


def makelist(thing):
    '''
    utility function to make any variable into a list
    '''
    thing = [thing] if not isinstance(thing, list) else thing
    return thing

def filter_df(metadf=pd.DataFrame(),samples=None,range=None,peptide_ranges=None,
                  charge=None,index=None,timept=None,timeidx=None,peptides=None,rep=None,data_ids=None,s_ids=None,quiet=True):
    ''' Utility function to Filter metadf (or any other dataframe) based on user specified values
        samples = ['sample1','sample2'] or 'sample1'
        range = [start,end]
        peptide_ranges = ['0001-0015','0030-0042'] or '0001-0015'
        charge = [1,2,3] or 1.0 
        index = [15,58,72] or [*range(0,50)] (Note: can just use filtered = metadf[0:50])
        timept = [0.0,5.0,1e6] or 0.0
        timeidx = [0,1,...,]
        peptides = ['PEPTIDESEQ','PEPITYPEPTIDE'] or 'PEPTIDESEQ'
        rep = [1,2,3] or 1
        data_ids = [0,1,...]
        s_ids = [123-4-1] #<data_id>-<time>-<rep>
    '''

    filtered = metadf.copy()
    if any(filtered.index.duplicated()): print("Warning, dataframe has duplicate indices. reset_index() to avoid selection errors")

    if samples:
        try: filtered = filtered[filtered['sample'].isin(makelist(samples))]
        except: 
            if not quiet: print("no column named sample")
    if not filtered.empty and range is not None:
        try: 
            if set(['start_seq','end_seq']).issubset(filtered.columns):
                filtered = filtered[(filtered['start_seq']>= range[0]) & (filtered['end_seq'] <= range[1])]
            elif 'peptide_range' in filtered.columns:
                filtered[['start_seq','end_seq']] = filtered['peptide_range'].str.split('-',expand=True).astype('int')
                filtered = filtered[(filtered['start_seq']>= range[0]) & (filtered['end_seq'] <= range[1])]
                filtered = filtered.drop(columns=['start_seq','end_seq'])
            else: print("Missing start/end seq or peptide_range column")
        except: 
            print("Filter error: specify range=[start,end]")
            return
    if not filtered.empty and charge is not None:
        try: filtered = filtered[filtered['charge'].isin(makelist(charge))]
        except: print("no column named charge")
    if not filtered.empty and (data_ids is not None or data_ids == 0):
        if 'data_id' in filtered.columns:
            filtered = filtered[filtered['data_id'].isin(makelist(data_ids))]
        else: index = data_ids
    if not filtered.empty and (index is not None or index == 0):
        filtered = filtered[filtered.index.isin(makelist(index))]
    
    if not filtered.empty and (timept is not None or timept == 0): #not in metadf but use to filter all_results_dataframe for Fixed_Pops
        try: filtered = filtered[filtered['time'].isin(makelist(timept))]
        except: print("no column named time")
    if not filtered.empty and (timeidx is not None or timeidx == 0): #not in metadf but use to filter all_results_dataframe for Fixed_Pops
        try: filtered = filtered[filtered['time_idx'].isin(makelist(timeidx))]
        except: print("no column named time_idx")
    if not filtered.empty and peptides is not None: #not in metadf but use to filter all_results_dataframe for Fixed_Pops
        try: filtered = filtered[filtered['peptide'].isin(makelist(peptides))]
        except: print("no column named peptide")
    if not filtered.empty and peptide_ranges is not None: 
        try: filtered = filtered[filtered['peptide_range'].isin(makelist(peptide_ranges))]
        except: print("no column named peptide_ranges")
    if not filtered.empty and (rep is not None or rep == 0): 
        try: filtered = filtered[filtered['rep'].isin(makelist(rep))]
        except: print("no column named rep")
    
    if not filtered.empty and (s_ids is not None):
        s_ids = makelist(s_ids)
        sidx = []
        for s_id in s_ids:
            d, t, r = [float(s) for s in s_id.split('-')]
            try: 
                sidx += filtered[(filtered['data_id'] == d) & (filtered['time'] == t) & (filtered['rep']==r)].index.tolist()
            except: print("no s_id",s_id)        
        filtered = filtered.loc[sidx]        
      
    if quiet == False: 
        print("Dataframe filtered to",len(filtered),"from",len(metadf),"total entries")
        if len(filtered) == 0: print("Warning: No datasets selected")
    return filtered


#safety function in case peptide sequence is bad
def goodseq(seq):
    try: 
        mass.most_probable_isotopic_composition(sequence=seq)
        return True
    except: 
        #print(f"Sequence {seq} is not defined")
        #exit()
        return False

def read_hexpress_data(f,dfrow,keep_raw = False,mod_dict={}):
    def pops(row): #just for test data
        pop = 0
        for col in ['p1','p2','p3']:
            if row[col] > 0: pop += 1
        return pop
    
    raw=[]
    all_raw = pd.DataFrame()
    dfs = pd.DataFrame()
    deutdata = pd.DataFrame()
    peaks=[]

    sample = dfrow['sample']
    start_seq = dfrow['start_seq']
    end_seq = dfrow['end_seq']
    peptide_range = dfrow['peptide_range']
    charge = dfrow['charge']
    peptide = dfrow['peptide']
    data_id = dfrow.name

    file=os.path.join(config.Data_DIR, f)
        
    try: 
        if (file[-4:] == 'xlsx') or (file[-3:] == 'xls'):
            ftype = 'excel'
            timepts = pd.read_excel(file,header=None,nrows=1) #get headers
        elif file[-3:] == 'csv':
            ftype = 'csv'
            timepts = pd.read_csv(file,header=None,nrows=1) #get headers
    except IOError as e:
        print (f"Could not read: {f} check path or if file is open")
        return deutdata
    #print(timepts)
    times =[x for x in timepts.values[0] if str(x) != 'nan']
    # print(times)

    delays = []
    for i,dtime in enumerate(times):
        rep = 1     
        if dtime[0:6] == 'undeut': delay = 0.0 
        elif dtime[0:2] == 'TD': delay = config.FullDeut_Time #1e6
        else: 
            tpstr = dtime.split(' ')[0].split('.')[0]
            if tpstr == '0':
                tp = float(dtime.split(' ')[0])  #leave as decimal, no unambiguous way to do replicates 
            else:
                tp = float(tpstr) 
            tunit = dtime.split(' ')[-1]
            #print(tp, tunit)
            delay = tp * np.power(60.0,'smh'.find(tunit[0]))
            # if Test_Data: 
            #     delay = tp * 60.0 #force to minutes, inconsistent dummy time units across sets
        delays += [delay]
        rep = Counter(delays)[delay]
        
        # #print(i,time)
        if ftype == 'excel':
            raw = pd.read_excel(file,skiprows=1,header=None,usecols=[i*2,i*2+1],names=['mz','Intensity']).dropna()
        elif ftype == 'csv':
            raw = pd.read_csv(file,skiprows=1,header=None,usecols=[i*2,i*2+1],names=['mz','Intensity']).dropna()
        raw = raw.sort_values('mz').reset_index()
        peaks = peak_picker( raw, peptide, charge, resolution=config.Peak_Resolution,count_sc=0,mod_dict=mod_dict)
        peaks['time']=delay
        peaks['data_id']=data_id
        peaks['sample']=sample
        peaks['peptide']=peptide
        peaks['charge']=charge
        peaks['rep']=rep
        peaks['peptide_range']=peptide_range
        peaks['start_seq']=start_seq
        peaks['end_seq']=end_seq
        peaks['file']=dfrow['file']
        if keep_raw:
            raw['time']=delay
            raw['data_id']=data_id
            raw['sample']=sample
            raw['peptide']=peptide
            raw['charge']=charge
            raw['rep']=rep
            raw['peptide_range']=peptide_range
            raw['start_seq']=start_seq
            raw['end_seq']=end_seq
            raw['file']=dfrow['file']
            all_raw = pd.concat([all_raw,raw],ignore_index=True)
        dfs = pd.concat([dfs,peaks],ignore_index=True)

    deutdata = dfs.copy()

    time_points = sorted(set(deutdata.time))
    #n_time_points = len(time_points)

    deutdata['time_idx'] = [ time_points.index(t) for t in deutdata.time ]

    ## read in solution file
    if config.Test_Data:
        file=os.path.join(config.Data_DIR, "bimodal_solutions2.txt")
        solution = pd.read_csv(file,delim_whitespace=True,)
        solution = solution.sort_values('time').reset_index()
        solution['npops'] = solution.apply(pops,axis=1)
        all_raw['time_idx'] = [ time_points.index(t) for t in all_raw.time ]
        return deutdata, all_raw, solution
    elif keep_raw: 
        all_raw['time_idx'] = [ time_points.index(t) for t in all_raw.time ]
        return deutdata, all_raw
    else: return deutdata

def read_specexport_data(csv_files,spec_path,row,keep_raw,mod_dict={}):
    '''
    read the HDExaminer exported data
    '''
    raw=[]
    dfs = []
    deutdata = pd.DataFrame()
    rawdata = pd.DataFrame()
    peaks=[]

    for f in csv_files:
        fileinfo = f.rsplit('.',1)[0].split('-')
        # print(fileinfo, len(fileinfo))
        rep = float(fileinfo[-2])
        if fileinfo[0] == 'Non': time = 0.0 
        elif fileinfo[0] == 'Full': time = config.FullDeut_Time #1e6
        else: time = float(fileinfo[0][:-1]) * np.power(60.0,'smh'.find(fileinfo[0][-1]))
        raw = pd.read_csv( os.path.join(spec_path,f),delimiter=",",header=None, names=["mz","Intensity"]).dropna()
        raw = raw.sort_values('mz').reset_index()
        peaks = peak_picker( raw, row['peptide'], row['charge'],resolution=config.Peak_Resolution,mod_dict=mod_dict )
        peaks['time']=time
        peaks['rep']=rep
        peaks['data_id']=row.name
        peaks['sample']=row['sample']
        peaks['charge']=row['charge']
        peaks['peptide']=row['peptide']
        peaks['peptide_range']=row['peptide_range']
        try:  peaks[['start_seq','end_seq']] = peaks['peptide_range'].str.split('-',expand=True).astype('int')
        except: pass
        peaks['file']=row['file']
        if peaks.Intensity.sum() > 0:
            dfs.append( peaks )         
        else: print (" File "+f+" contains no Intensity data at expected m/z values") #peaks['sample']+' '+peaks['peptide']+
        if keep_raw:
            raw['time']=time
            raw['data_id']=row.name
            raw['sample']=row['sample']
            raw['peptide']=row['peptide']
            raw['charge']=row['charge']
            raw['rep']=rep
            raw['peptide_range']=row['peptide_range']
            try: raw[['start_seq','end_seq']] = raw['peptide_range'].str.split('-',expand=True).astype('int')
            except: pass 
            raw['file']=row['file']
            rawdata = pd.concat([rawdata,raw],ignore_index=True)
    if len(dfs) > 0: 
        deutdata = pd.concat(dfs, ignore_index=True,)
        time_points = sorted(set(deutdata.time))
        #n_time_points = len(time_points)
        deutdata['time_idx'] = [ time_points.index(t) for t in deutdata.time ]
        deutdata['charge'] = row['charge']
        
    if keep_raw:
        time_points = sorted(set(rawdata.time)) 
        rawdata['time_idx'] = [ time_points.index(t) for t in rawdata.time ]
        return deutdata, rawdata
    else: return deutdata


def get_na_isotope(peptide,charge,npeaks=None,mod_dict={}):
    '''
    Get the Natural Abundance isotopic pattern for a given peptide,charge
    If npeaks=None the isotopic_variants routine will vary npeaks depending on composition
    e.g. 'YGGFL' will have 7 peaks but 'RDKVQKEYALFYKLD' has 15. 
    Letting value float depending on peptide length
    '''
    pepcomp = {}
    na_isotope=[]
    if goodseq(peptide): comp = mass.Composition(peptide) 
    else: comp = {}
    for key in list(comp):
        pepcomp[key] = comp[key]
    pepcomp['Hex']=count_amides(peptide,count_sc=0.0)

    #pepcomp = {'H': 53, 'C': 34, 'O': 15, 'N': 7}
    if mod_dict:
        for key in mod_dict:
            if key in pepcomp:
                pepcomp[key] = pepcomp[key] + mod_dict[key]
        pepcomp = {**pepcomp,**mod_dict}

    if set(['H','Hex']).issubset(pepcomp.keys()): 
        pepcomp['H'] = pepcomp['H'] - pepcomp['Hex']
    if 'Hex' in pepcomp.keys(): pepcomp.pop('Hex')
    if 'ion_type' in pepcomp.keys(): pepcomp.pop('ion_type')
    if 'undeut_mz' in pepcomp.keys(): pepcomp.pop('undeut_mz')

    theoretical_isotopic_cluster = isotopic_variants(pepcomp, npeaks=npeaks, charge=charge)

    for ipeak in theoretical_isotopic_cluster:
        na_isotope = np.append(na_isotope, ipeak.intensity) 

    return na_isotope

def choose_na(peptide,charge,mod_dict = {},deutdata=None,user_env=None):
    '''Get the user_specified NA_envelope for Current_Isotope 
        otherwise defaults to get_na_isotope()

        If the NA_envelop column is in the metadf file:
        1) '','normal','peptide' will use get_na_isotope()
        2) user specified array of values (may be string if dataframe is read from csv)
        3) 'mixed#' will fit # na_iso envelopes offset by 1 for each #
        4) 'infer','UnDeut' will use the peakpicked UnDeut spectra (avg over reps)
    '''

    def mixed_state(na_iso,*fracs):
        #assume offset by 1 bin for len(fracs)
        #e.g. c-2,c-1,c .. or z,z+1,z+2
        ysum = np.zeros(len(na_iso)+len(fracs))
        for i in range(len(fracs)):
            y = np.concatenate([np.zeros(i),na_iso,np.zeros(len(fracs)-i)])
            ysum += y*fracs[i]
        return ysum

    na_iso = []
       
    if user_env is not None:
        #user_env = row['NA_envelope']
        if str.lower(str(user_env)) in ['','normal','peptide']:
            na_iso = get_na_isotope(peptide,charge,mod_dict=mod_dict)
        elif isinstance(user_env,np.ndarray):
            na_iso = user_env
        elif ('[' in str(user_env)) and (']' in str(user_env)):
            if ',' in na_iso: #somehow have commas, most likely if input in session
                na_iso = [float(ue) for ue in str(user_env)[1:-1].split(',')]
            else: #more likely have spaces so can save as .csv 
                na_iso = [float(x) for x in str(user_env)[1:-1].split()]
            na_iso = np.array(na_iso)
        elif str.lower(str(user_env)).startswith('mixed'):
            n_mixed = int(user_env[5:])            
            fit_na_iso = []
            focal_data = deutdata.copy()[(deutdata.time==0)]  
            for rep in focal_data.rep.unique():
                focal_data = deutdata.copy()[(deutdata.time==0) & (deutdata.rep==rep)]           
                #mz=np.array(focal_data.mz.copy())
                y=np.array(focal_data.Intensity.copy())
                npeaks = len(y)-n_mixed  #y picked peaks should already be zero filled
                temp_na_iso = get_na_isotope(peptide,charge,mod_dict=mod_dict,npeaks=npeaks)
                temp_na_iso = np.concatenate([temp_na_iso,np.zeros(npeaks - len(temp_na_iso))])
                fit, covar = curve_fit( mixed_state, temp_na_iso, y, p0=(1,)*n_mixed, 
                                       maxfev=int(1e6), bounds = ([0]*n_mixed,[np.inf]*n_mixed)   )   
                yfit = mixed_state(temp_na_iso,*fit)   
                fit_na_iso += [list(yfit)]
            na_iso = list(map(lambda idx: sum(idx)/float(len(idx)),zip(*fit_na_iso)))
        elif str.lower(str(user_env)) in ['infer','UnDeut']:    
            focal_data = deutdata.copy()[(deutdata.time==0)]  
            na_iso = focal_data.groupby('n_deut')['Intensity'].mean().tolist()
        else: print("unknown NA_envelope specification")
        
    else: 
        na_iso = get_na_isotope(peptide,charge,mod_dict=mod_dict)
    if np.sum(na_iso) == 0:
        print("no data found for fit Isotopic envelope, using default")
        na_iso = get_na_isotope(peptide,charge,mod_dict=mod_dict)
 
    return na_iso/np.sum(na_iso)   #require normalization


def count_amides (peptide,count_sc=0.0):
    '''
    calculate the number of exchanging amides
    '''
    ex_sc = 0
    proline = peptide[config.Nterm_subtract:].count('P')
    for sidechain in 'STYCDEHW':
        ex_sc += peptide.count(sidechain)
    for sidechain in 'R':
        ex_sc += 2*peptide.count(sidechain)
    for sidechain in 'KQN':
        ex_sc += 2*peptide.count(sidechain)
    n_amides = max(0,len(peptide)-proline-config.Nterm_subtract+int(ex_sc*count_sc))
    return n_amides


def peak_picker(data, peptide,charge,resolution=50.0,count_sc=0.0,mod_dict={}):
    '''
    find the peaks at the expected m/z values (+/- resolution in ppm)
    '''
    padding = config.Zero_Filling 
    min_pts = config.Max_Pops*3+2 #minimum number of datapoints necessary for max pop fit
    na_buffer = len(get_na_isotope(peptide,charge,mod_dict=mod_dict))//2   
    n_amides = max(count_amides(peptide,count_sc=0.0),min_pts) + na_buffer + padding #include count from Isotopic Envelope

    if 'Hex' in mod_dict.keys():
        n_amides += mod_dict['Hex']  
    if 'ion_type' in mod_dict.keys(): 
        ion_type = mod_dict['ion_type']
        #mod_dict.pop('ion_type')
    else: ion_type = 'M'
    
    if 'undeut_mz' in mod_dict.keys(): #use user specified value and don't compute anything, ignore modification 
        undeut_mz = mod_dict['undeut_mz']
    else:
        undeut_mz = 0.0
        if len(peptide)>0:
            undeut_mz = mass.calculate_mass(sequence=peptide,show_unmodified_termini=True,charge=charge,ion_type=ion_type)
        if mod_dict:
            mod_comp = {}          
            comp_keys = set(list(mod_dict.keys())) - set(['Hex']) - set(['ion_type']) - set(['undeut_mz'])
            if len(comp_keys) > 0:
                for k in comp_keys:
                    mod_comp[k] = mod_dict[k]
                undeut_mz += mass.calculate_mass(composition=mod_comp,charge=charge,ion_type=ion_type)

    #print("undeut_mod",undeut_mz)
    #n_deut = np.arange(n_amides+1) #ExMS instead
    #pred_mzs = undeut_mz + (n_deut*1.006277)/charge #ExMS instead
    # Use ExMS method of expected mz values https://pubs.acs.org/doi/10.1007/s13361-011-0236-3
    delta_m = [1.003355, 1.003355, 1.004816, 1.004816, 1.006277] #dmC, dmC, dmHCavg, dmHCavg, dmH ...
    dmz = [0.0]
    for nd in range(0,n_amides):
        dm = delta_m[nd] if nd < 5 else delta_m[-1]
        dmz += [dmz[nd] + dm/charge]
    pred_mzs = np.array(dmz) + undeut_mz

    mz_mid = pred_mzs.mean() 

     #attempt to set threshold based on user noise value 
    if config.setNoise: 
        threshold = config.setNoise
        #print("config.SetNoise threshold",threshold)
    else: 
        threshold = config.Y_ERR/100.0 * data.Intensity.max()
        #print("config.Y_ERR threshold",threshold)

    peaks = []
    zeroes = 0

    # need to make sure the max_Int is a peak so we don't grab side values from overlapping peaks
    for i,pred_mz in enumerate(pred_mzs):
        #mass accuracy is in ppm, so default is 50ppm
        mz_range = pred_mz * (1.0+np.array([-1,1])*resolution/1e6)    
        #sort data in range, grab first entry which is the peak
        focal_data = data.copy()[data.mz.between(*mz_range)]
        focal_data['n_deut'] = i 
        focal_data = focal_data.sort_values('Intensity',ascending=False)#.reset_index(drop=True)
        
        if (len(focal_data) > 0):
            if len(focal_data) < 10: #assumes stick data or sparse data
                max_Int = focal_data['Intensity'].max()
            elif (focal_data.index[0] not in (focal_data.index.min(),focal_data.index.max())): #avoid shoulders
                max_Int = focal_data['Intensity'].max()
            else: max_Int = 0.0
            focal_data.reset_index(drop=True)
            
            #testing keeping value as threshold but counting as a zero
            # still getting 0 intenist peaks instead of low intensity values 
            intensity = max_Int
            if (intensity < threshold and pred_mz > mz_mid): zeroes += 1
            #previous 2 lines instead of following 3
            # if max_Int < threshold: intensity = 0.0
            # else: intensity = max_Int
            # if (intensity == 0.0 and pred_mz > mz_mid): zeroes += 1
        else:
            intensity = 0.0
            if (pred_mz > mz_mid): zeroes += 1
        #I don't remember what this next bit was meant to do but it zeroes out some real peaks >< 
        # if len(pred_mzs) - i < padding + 1:
        #     intensity = 0.0
        #     zeroes += 1
        peak = pd.DataFrame({'mz':[pred_mz],'Intensity':[intensity],'n_deut':[i]})
        peaks.append( peak )
        if (zeroes >= padding) and (i >= min_pts): 
            break #if we have enough zeroes and enough points, we're done peak_picking
    #print(n_amides, zeroes, mz_mid)
    peaks = pd.concat(peaks,ignore_index=True)

    # #clean up right side zero points
    # x = peaks['Intensity'].to_numpy()
    # cz = 0
    # for i in range(min_pts,len(x))[::-1]:
    #     if x[i] < threshold:
    #         cz += 1
    #     else: break
    # # print ("zero points:",cz)
    # # print ("points to cut:",(cz-padding))
    # # print ("last index:",len(new_x)-(cz-padding + 1))

    # l_idx = len(x) - (cz - padding + 1)
    # peaks = peaks[peaks['n_deut'] <= l_idx]

    ## adding this section to get X_features at the peakpick step
    try:
        envelope_height = peaks['Intensity'].max() * config.Env_threshold
        env, env_Int = get_mz_env(envelope_height,peaks,pts=True)
        y=np.array(peaks.Intensity.copy())
        env_symmetry_adj = 2.0 - (y.max() - env_Int)/y.max()
        peaks['env_width'] = charge*(env[1]-env[0])
        peaks['env_symm'] = env_symmetry_adj
        #peaks['skewness'] = skew(y_norm,bias=False)
    except:
        pass
    peaks['max_namides']=count_amides(peptide,count_sc=0.0)

    return peaks #pd.concat(peaks,ignore_index=True)

def get_mz_env(value, df, colname='Intensity',pts=False):
    '''
    get mz values at value = envelope_height (e.g. 0.1*maxIntensity)
    to define the envelope width, for assessment of expected polymodal fits
    '''
    if value == 0:
        if pts: return np.array([np.nan, np.nan]), np.nan
        else: return np.array([np.nan, np.nan])
    df = df.copy().reset_index()
    boolenv = df[colname].gt(value)
    #envpts = boolenv.value_counts()[True] # number of points in envelope
    loweridx = df[colname].where(boolenv).first_valid_index()
    upperidx = df[colname].where(boolenv).last_valid_index()
    if df[colname].iloc[0] > value:
        min_mz = df['mz'].iloc[0] #envelope starts with first datapoint
        left_Int = df[colname].iloc[0] 
    else:     
        x1a = df[colname].iloc[loweridx]
        x1b = df[colname].iloc[loweridx-1]
        y1a = df['mz'].iloc[loweridx]
        y1b = df['mz'].iloc[loweridx-1]
        min_mz = y1a + (y1b - y1a) * (value - x1a)/(x1b - x1a)
        left_Int = value
    
    if df[colname].iloc[-1] > value:
        max_mz = df['mz'].iloc[-1] #envelope goes to end of datapoints, poorly picked masspec?
    else:
        x2a = df[colname].iloc[upperidx]
        x2b = df[colname].iloc[upperidx+1]
        y2a = df['mz'].iloc[upperidx]
        y2b = df['mz'].iloc[upperidx+1]
        max_mz = y2a + (y2b - y2a) * (value - x2a)/(x2b - x2a)
    
    if pts: return np.array([min_mz, max_mz]), left_Int
    else: return np.array([min_mz, max_mz])


## Multi-binomial Functions
def nCk_real(n,k):
    '''
    real (vs integer) version of the NchooseK function
    '''
    #print("n,k:",n,k)
    if n == k: return 1.0
    elif n - k + 1 <= 0: return 0.0
    #if n - k <= 0: return 0.0
    elif n+k > 50: # https://stats.stackexchange.com/questions/72185/how-can-i-prevent-overflow-for-the-gamma-function
        log_nCk = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
        return exp(log_nCk)
    else: return gamma(n+1)/(gamma(k+1)*gamma(n-k+1))

def binom(bins, n, p):
    '''
    basic (but not simple) binomial function
    '''      

    k = np.arange(bins+1).astype('float64')
    nCk = [nCk_real(n,y) for y in k]

    with np.errstate(all='ignore'):
        # suppress errors: np.power may nan or overflow but we fix it after

        # this attempt was unstable near p = 1
        # fp1 = np.float_power(p,k)
        # fp2 = np.float_power(1-p,n-k) #fits, but np.power overflows     
        # binoval = nCk*fp1*fp2
        # binoval = np.nan_to_num(binoval,nan=0.0)
        
        
        t_fp1, t_fp2, t_binoval = [], [], []
        for ki in k.astype(int):
            if nCk[ki] == 0:
                t_binoval += [0.0]
                t_fp1 += [0.0]
                t_fp2  += [0.0]
            else: 
                t_fp1 += [np.power(p,ki)]
                fp2_val = np.clip(np.float_power(1-p,n-ki),sys.float_info.min,sys.float_info.max)
                t_fp2  += [fp2_val]
                t_binoval += [nCk[ki] * t_fp1[ki] * t_fp2[ki]]

    binoval = np.array(t_binoval/sum(t_binoval))
    binoval = np.nan_to_num(binoval,nan=0.0)
       
    return binoval

def binom_isotope(bins, n,p):
    '''
    binomial function using the Natural Abundance isotopic envelope
    '''
    bs = binom(bins,n,p)
    newbs=np.zeros(len(bs) + len(config.Current_Isotope)+1)
    for i in range(len(bs)):
        for j in range(len(config.Current_Isotope)):     
            newbs[i+j] += bs[i]*config.Current_Isotope[j]  
    return newbs[0:bins+1]

def n_binomials( bins, *params ): #allfracsversion
    '''
    basic binomial for n populations
    '''
    # params takes the form [scaler, n_1, ..., n_n, mu_1, ..., mu_n, frac_1, ..., frac_n] 
    n_curves = int( (len(params)+1) / 3.0 )
    log_scaler = params[0]
    n_array=np.array(params[1:n_curves+1])
    mu_array = np.array( params[n_curves+1:2*n_curves+1] )
    frac_array=[]
    frac_array = np.array( params[ -n_curves: ] )
    frac_array = frac_array/np.sum(frac_array)
    poissons = [ frac * binom( bins,n, mu ) for frac, n, mu in zip( frac_array, n_array, mu_array ) ]
    return np.power( 10.0, log_scaler ) * np.sum( poissons, axis=0, )


def n_binom_isotope( bins, *params ): #allfracsversion
    '''
    n population binomial using the Natural Abundance isotopic envelope
    '''
    # params takes the form [ scaler, mu_1, ..., mu_n, frac_1, ..., frac_n] 
    n_curves = int(( len(params) + 1) / 3.0 )
    log_scaler = params[0]
    n_array = np.array( params[1:n_curves+1] )
    mu_array = np.array( params[n_curves+1:2*n_curves+1] )
    frac_array = np.array( params[ -n_curves: ] )
    frac_array = frac_array/np.sum(frac_array)
    poissons = [ frac * binom_isotope( bins, n, mu ) for frac, n, mu in zip( frac_array, n_array, mu_array ) ]
    truncated = np.power( 10.0, log_scaler ) * np.sum( poissons, axis=0, )[0:bins+1]
    return truncated 

def calc_rss( true, pred,yerr_systematic=0.0,norm=True ):
    '''
    calculate the residual squared sum
    '''
    if norm: n = len(true)
    else: n = 1.0
    return np.sum( (pred-true)**2 + yerr_systematic**2 ) / n

def get_params(*fit, sort = False, norm = False, unpack = True, returnIndex = False):
    '''
    extract the fit parameters from the fits 
    ''' 
    # assuming all fracs
    # binom eqs of form: scaler, nex * n_curves, mu * n_curves, frac * n_curves
    # if sort then reorder nex,mu,frac in order of mu*nexs
    # if norm, then fracs will be scaled to sum to 1.0
        # need to be mindful of mixing up corresponding errors when sorting and norming

    num_curves = int((len(fit)-1)/3)

    scaler = fit[0]
    
    nexs = np.array(fit[1:num_curves+1])
    mus = np.array(fit[num_curves+1:num_curves*2+1])
    fracs = np.array(fit[-num_curves:])
    use_index = np.arange(len(mus))

    if sort:
        mn = nexs*mus
        use_index = mn.argsort()
        nexs = nexs[use_index]
        mus = mus[use_index]#[::-1]]
        fracs = fracs[use_index]#[::-1]]
    if norm: 
        sumfracs = np.sum(fracs)
        #print("fracs",fracs,sumfracs)
        if sumfracs > 0: fracs = fracs/sumfracs
    
    if unpack:
        if returnIndex:
            return scaler, nexs, mus, fracs, use_index
        else: return scaler, nexs, mus, fracs
    else:
        if returnIndex:
            return np.concatenate((np.array([scaler]),nexs,mus,fracs)), use_index
        else: return np.concatenate((np.array([scaler]),nexs,mus,fracs))

def init_params(n_curves,max_n_amides,seed=None):
    '''
    generate intial guess parameters for the fits 
    '''
    
    init_rng = np.random.default_rng(seed=seed)

    log_scaler_guess = 0.0  
    nex_guess = list(max_n_amides/config.Nex_Max_Scale*init_rng.random(n_curves)) #was same for all at max/2 in < ver3.0
    nex_low = 0.0 
    sampled = init_rng.random(n_curves*2) #array for both mus and fracs
    mu_guess =  list(sampled[0:n_curves]) #[0.9, 0.9, 0.9] #
    frac_guess = sampled[-n_curves:]#[0.70, 0.25, 0.05] #allfracs
    frac_guess = list(frac_guess/np.sum(frac_guess))
    frac_uppers = [1.0]*n_curves

    initial_estimate = [ log_scaler_guess ] + nex_guess + mu_guess[0:n_curves] + frac_guess[0:n_curves]
    lower_bounds = [ 0.0 ] + [nex_low]* n_curves + [0.0]* n_curves*2 
    upper_bounds = [ 1.0, ] + [max_n_amides]* n_curves + [ 1.0 ] *n_curves + frac_uppers[0:n_curves]
    bounds = ( lower_bounds, upper_bounds, )

    #print("intial seed, initial estimate",seed,initial_estimate)

    return initial_estimate, bounds

def minimize_func(params, *data):
    '''
    residual sum squared to be minimized using scipy.optimize.differential_evolution
    '''
    bins,y = data
    result = np.sum( (n_fitfunc(bins,*params)-y)**2 )

    return result #**0.5


def fit_bootstrap(p0_boot, bounds, datax, datay, sigma_res=None,yerr_systematic=0.0,nboot=100,
                  ax=None,yscale=1.0):
    '''
    perform a bootstrap type routine, nboot independent fits of the data + noise
    '''
    #p0_boot is a list of all p0 initial parameters with nboot entries
    #if ax != None: ax.plot( mz, datay, color = 'cyan', linestyle='solid',label='boot_datay' )

    p0 = p0_boot[0]
    num_curves = int((len(p0)-1)/3)
    
    if sigma_res==None:
        # Fit first time if no residuals
        pfit, perr = curve_fit( n_fitfunc, datax, datay, p0, maxfev=int(1e6), 
                                        bounds = bounds   )

        print("Ran initial bootstrap curve_fit to generate residuals")

        # Get the stdev of the residuals
        residuals = n_fitfunc(datax,*pfit) - datay
        sigma_res = np.std(residuals)
    sigma_err_total = np.sqrt(sigma_res**2 + yerr_systematic**2)

    # nboot random data sets are generated and fitted
    ps = []
    ps_cov = []
    centers = []
    boot_residuals = []
    for i in range(nboot):
        p0 = p0_boot[i]
        randomdataY = []
        randomDelta = []
        randomDelta = rng.normal(0., sigma_err_total, len(datay))
        randomDelta = [0 if a==0 else b for a,b in zip(datay,randomDelta)] #don't change y if y=0
        randomdataY = datay + randomDelta 
        #if datapoint is zero, leave as zero and don't let value go negative
        randomdataY = np.clip(randomdataY, 0.0, np.array(datay)*2)#np.inf) 

        #print(randomdataY)

        if (len(p0) > len(randomdataY)): 
                            print(f"attempting to fit more parameters than data points in bootstrap")
                            #should be able to exit at this point, haven't updated last fit parameters including p_err
                            break # exit the for n_curves loop
        try:
            #print(len(p0),len(bounds[0]),len(bounds[1]))
            randomfit, randomcov = curve_fit( n_fitfunc, datax, randomdataY, p0, maxfev=int(1e6), 
                                    bounds = bounds   )
        except RuntimeError:
            break
        
        
        randomfit[0] = np.power( 10.0, randomfit[0] )
        #rfit = randomfit
        rfit,use_index = get_params(*randomfit,sort=True,norm=True,unpack=False,returnIndex=True) #randomfit #
        ps.append(rfit)

        pcov = np.sqrt(np.diag(randomcov))
        rcscaler = pcov[0]  
        rcnexs = np.array(pcov[1:num_curves+1])[use_index]
        rcmus = np.array(pcov[num_curves+1:num_curves*2+1])[use_index]
        rcfracs = np.array(pcov[-num_curves:])[use_index]
        pcov_list = [rcscaler] + list(rcnexs) + list(rcmus) + list(rcfracs)

        ps_cov.append(pcov_list) ## this won't work if ever Sorting

        if ax != None: ax.plot( mz, randomdataY*yscale, color = 'orchid', linestyle='dashed', linewidth=2)
        # ax.plot( mz, randomdataY, color = 'green', linestyle='dashed', )

        tempr = rfit.copy()
        tempr[0] = np.log10(rfit[0])
        boot_y = n_fitfunc(datax, *tempr) 
        boot_residual = calc_rss(boot_y , datay)
        
        if ax != None: ax.plot( mz, boot_y*yscale, color = 'darkviolet', linestyle='solid', alpha=0.2 )
        s,n,m,f = get_params(*rfit,norm=True,unpack=True) #already did 10**s in rfit
        kcenter = [] #get centroids of each population
        for k in range( num_curves ):
            bfit_yk = s * f[k] * fitfunc( datax, n[k], m[k], ) * yscale
            kcenter += [sum(bfit_yk * mz)/sum(bfit_yk)]
            #plot_label = ('pop'+str(k+1)+' = '+format(frac,'.2f')+'\nNex'+str(k+1)+' = '+format(nex,'.1f'))
            if ax != None: ax.plot( mz, bfit_yk, color = 'green', linestyle='dashed',linewidth=3,alpha=0.2)#label=plot_label)
        centers.append(kcenter)
        boot_residuals.append(boot_residual/datax) #normalize by number of bins
    #print("boot_residuals:",len(boot_residuals),boot_residuals)
    ps = np.array(ps)
    ps_cov = np.array(ps_cov)
    # mean_pfit = np.mean(ps,axis=0)

    # You can choose the confidence interval that you want for your
    # parameter estimates: 
    Nsigma = 1. # 1sigma corresponds to 68.3% confidence interval
                # 2sigma corresponds to 95.44% confidence interval
    if ps.ndim == 2: ps[:,0] = np.log10(ps[:,0])
    # err_pfit = Nsigma * np.std(ps,axis=0)

    # mean_pfit[0] = np.log10(mean_pfit[0])

    # pfit_bootstrap = mean_pfit
    # perr_bootstrap = err_pfit
    return ps, boot_residuals, np.array(centers), ps_cov #return all the bootstrap fits 
    #return pfit_bootstrap, perr_bootstrap, np.array(centers)

def get_TDenv(datafits,mod_dict={}):
   '''
   haven't updated this to use choose_na()
   calculate the TD envelope to use as an X_feature in the ML prediction
   '''
   #global Current_Isotope
   df = datafits.copy()
   test_set = df[['peptide','max_namides']].drop_duplicates()
   test_peptides = dict(zip(test_set['peptide'],test_set['max_namides'])) 
   for peptide,namides in test_peptides.items():
      charge = 1
      config.Current_Isotope= get_na_isotope(peptide,charge,npeaks=None,mod_dict=mod_dict) 
      TD_max_spectrum = n_binom_isotope(namides+5,0.0, namides, config.Dfrac, 1.0) #use Dfrac for expected TDenv
      TD_spec = pd.DataFrame(zip(np.arange(len(TD_max_spectrum)),TD_max_spectrum),columns=['mz','Intensity'])
      [left,right] = get_mz_env(0.1*max(TD_max_spectrum),TD_spec,colname='Intensity')
      TD_env_width = (right - left)*charge
      peptide_idx = df[(df['peptide'] == peptide)].index
      df.loc[peptide_idx,'TD_env_width'] = TD_env_width
   return df

def predict_pops(trained_model,datafits):
   '''
   Predict 1 or more populations based on trained model and X_features
   '''
   import tensorflow as tf
   from keras.models import load_model

   df = datafits.copy()
   if 'TD_env_width' not in df.columns:
      df = get_TDenv(df)
   X_features = ['env_width','env_symm','max_namides','TD_env_width'] #"model4"
   Xtest = df[X_features].to_numpy()
   Xtest_pred = trained_model.predict(Xtest)
   ytest_pred = np.argmax(Xtest_pred,axis=1)
   ytest_pred_binary = np.array([min(y,1) for y in ytest_pred])+1
   df['pred_pops'] = ytest_pred_binary
   return df 


## Function to perform the fits on the metadf list of spectra
def run_hdx_fits(metadf,user_deutdata=pd.DataFrame(),user_rawdata=pd.DataFrame(),update_deutdata = False):
    '''
    Run the whole thing based on the metadf file: read in data, pick peaks, do fits, plot it all
    metadf: meta data dataframe with list of all sample/peptide/charge/filename information
    user_alldeutdata: optional peak picked data, to limit run (update_deutdata = False)
                    or to update a subset of deutdata (full set determined from metadf)
    user_allrawdata: user specified raw data instead of that found at filename
    '''
    
    GP = config.Generate_Plots
    global n_fitfunc, fitfunc, mz, now, date, deutdata, rawdata, reportdf #Current_Isotope
    global deutdata_all, rawdata_all, solution, data_fits, data_fit, config_df, fitparams_all

    def MinimizeStopper(xk=None,convergence=None):
        '''
        callback function to stop scipy.optimize.differential_evolution function
        when the minimize_func value is below a threshold
        use with config.DiffEvo_kwargs = {'callback':MinimizeStopper,'polish':True,'maxiter':100}
        '''
        return minimize_func(xk,*args) < config.DiffEvo_threshold

    now = datetime.now()
    date = now.strftime("%d%b%Y")

    os.makedirs(os.path.dirname(config.Output_DIR),exist_ok=True)
    #fout = open(os.path.join(config.Output_DIR,'output_v3allfracs_'+date+'.txt'),'w')
    data_output_file_inprogress = os.path.join(config.Output_DIR,"data_fits_asrun_"+date+".csv")

    
    if (not user_deutdata.empty) and (update_deutdata == False): #filter metadf to contain only user_deutdata sets 
        print("Setting metadf to data in user supplied deutdata")
        metadf = filter_df(metadf,data_ids=list(user_deutdata.data_id.unique()))
    
    reportdf = metadf.copy()
        
    save_metadf(reportdf,filename="metadf_asrun_"+date+".csv") #save the metadf for whatever is run to the Output_Dir

    deutdata_all = pd.DataFrame()
    rawdata_all = pd.DataFrame()
    data_fits = pd.DataFrame()
    fitparams_all = pd.DataFrame()

    sample_id_cols = ['data_id','sample', 'peptide', 'peptide_range','start_seq','end_seq','charge', 'time','time_idx', 'rep']
  
    # if config.USE_PARAMS_FILE:
    #     get_parameters()
    if config.WRITE_PARAMS:
        write_parameters(config.Output_DIR, config.Allow_Overwrite)


    Preset_Pops = False
    if config.Preset_Pops:
        try: 
            config_df = pd.read_csv(config.Preset_Pops_File).drop('Index',axis=1)
            print("Using user specified number of populations when available")
            Preset_Pops = True
            max_pops = config_df['max_pops'].max() #this is to determine columns in data_fit output            
        except: 
            print("Preset_Pops_File does not exist or improperly formatted, using default range")
            Preset_Pops = False
            max_pops = config.Max_Pops
    else: max_pops = config.Max_Pops

    n_fitfunc = n_binom_isotope # n_binomials #
    fitfunc = binom_isotope # binom #
    fitfunc_name=str(fitfunc).split()[1]  #fit function as string to add to filenames

    ## generalize to preset column headers corresponding to max_pops
    ## Columns may change still depending on X_features used for pops prediction
    data_fit_columns = [ 'data_id','sample', 'peptide', 'peptide_range','start_seq','end_seq','charge', 'time','time_idx', 'rep', 
                        'max_Int','centroid', 'env_width', 'env_symm', 'max_namides','UN_TD_corr','fit_pops','p-value']
    for imp in range(1,max_pops+1):
        data_fit_columns += ['centroid_'+str(imp), 'Dabs_'+str(imp),'Dabs_std_'+str(imp),'pop_'+str(imp), 'pop_std_'+str(imp), ]
    #                         'mu_'+str(imp), 'Nex_'+str(imp), 'Nex_std_'+str(imp)]
    if config.Test_Data: data_fit_columns += ['solution_npops']
    data_fit = pd.DataFrame(columns = data_fit_columns) #also reinitialize each rep
    data_fit.to_csv(data_output_file_inprogress,index=True,index_label='Index',header=True) #create new file 



    ### Process all the Data ###
    dataset_count = 0
    for index, row in metadf.iterrows():
        dataset_count += 1
        overlay_reps = config.Overlay_replicates
        rawdata = pd.DataFrame(None) 
        deutdata = pd.DataFrame(None)         

        mod_dic = {}
        if 'modification' in row.keys():
            mod = row['modification']
            mod = mod.split() #'H:1 Hex:1 ion_type:c undeut_mz'
            mod_dic = {x.split(':')[0]:(x.split(':')[1]) for x in mod}
            for k,v in mod_dic.items():
                if k == 'undeut_mz':
                    mod_dic[k] = float(v)                    
                elif k != 'ion_type':
                    mod_dic[k] = int(v)

        hdx_file = row['file']
        sample = row['sample']
        peptide = row['peptide']
        charge = row['charge']
        peptide_range = row['peptide_range']
        
        shortpeptide = peptide if len(peptide) < 22 else peptide[:9]+'...'+peptide[-9:] #for tidy plot labels
        
        print("\n\nDataset",index,"(",dataset_count,"of",len(metadf),")")
        print("Performing fits for "+sample+" "+peptide_range+": "+peptide+" z="+str(int(charge)))
        
        #update config.Max_Pops to pick the correct number of peaks 
        if config.Preset_Pops:
            try: 
                fixed_pop = filter_df(config_df,samples=sample,peptides=peptide,charge=charge,quiet=True)#['min/max_pops'][0]
                temp_max_pops = config.Max_Pops
                config.Max_Pops  = int(fixed_pop['max_pops'].values[0]) 
            except: pass 

        if config.Data_Type == 1:
            if config.Test_Data:
                config.Keep_Raw == True
                deutdata, rawdata, solution = read_hexpress_data(hdx_file,row,keep_raw=config.Keep_Raw)
            elif (config.Keep_Raw): 
                if (user_deutdata.empty) or (update_deutdata == True): #only peakpick if necessary
                    deutdata, rawdata = read_hexpress_data(hdx_file,row,keep_raw=config.Keep_Raw,mod_dict=mod_dic)     
                    # TODO need to add option if user_rawdata specified to just peak_pick user raw              
            else: 
                if (user_deutdata.empty) or (update_deutdata == True):
                    deutdata = read_hexpress_data(hdx_file,row,keep_raw=config.Keep_Raw,mod_dict=mod_dic)
                    # TODO need to add option if user_rawdata specified to just peak_pick user raw
        else: # Data_type == 2, already checked that it is 1 or 2            
            if config.Keep_Raw:
                if (user_deutdata.empty) or (update_deutdata == True):
                    spec_path = os.path.join(config.Data_DIR,row['sample'],row['file'])
                    csv_files = [ f for f in os.listdir(spec_path) if f[-5:]==str(int(charge))+'.csv'  ]
                    deutdata, rawdata = read_specexport_data(csv_files,spec_path,row,keep_raw=config.Keep_Raw,mod_dict=mod_dic)
                    # TODO need to add option if user_rawdata specified to just peak_pick user raw
            else: 
                if (user_deutdata.empty) or (update_deutdata == True):
                    spec_path = os.path.join(config.Data_DIR,row['sample'],row['file'])
                    csv_files = [ f for f in os.listdir(spec_path) if f[-5:]==str(int(charge))+'.csv'  ]
                    deutdata = read_specexport_data(csv_files,spec_path,row,keep_raw=False,mod=mod_dic)
                    # TODO need to add option if user_rawdata specified to just peak_pick user raw
        
        # deutdata['data_id'] = index  #should be in read_ functions
        # rawdata['data_id'] = index

        # update with user values ... 
        if not user_deutdata.empty:

            userdeut = filter_df(user_deutdata,data_ids=index,samples=sample,peptides=peptide,charge=charge,quiet=True)

            if not userdeut.empty:  # deutdata.loc[userdeut.index] = userdeut # 
                if update_deutdata == True:
                    update_deut_list = zip(*[userdeut[col] for col in ['data_id','sample','peptide_range','time','rep','charge']])
                    update_deut_list = list(dict.fromkeys(update_deut_list))
                    deut_filter_idx = []
                    for di,ss,pp,tt,rr,zz in update_deut_list:
                        spec = {'data_ids':di,'samples':ss,'peptide_ranges':pp,'timept':tt,'charge':zz,'rep':rr}
                        deut_filter_idx += list(filter_df(deutdata,**spec).index)
                    deutdata = pd.concat([deutdata.drop(index=deut_filter_idx),userdeut],ignore_index=True)
                else: deutdata = userdeut
        if not user_rawdata.empty:
            userraw = filter_df(user_rawdata,data_ids=index,samples=sample,peptides=peptide,charge=charge,quiet=True)
            if not userraw.empty:  # deutdata.loc[userdeut.index] = userdeut # 
                if update_deutdata == True:
                    update_raw_list = zip(*[userraw[col] for col in ['data_id','sample','peptide_range','time','rep','charge']])
                    update_raw_list = list(dict.fromkeys(update_raw_list))
                    raw_filter_idx = []
                    for di,ss,pp,tt,rr,zz in update_raw_list:
                        spec = {'data_ids':di,'samples':ss,'peptide_ranges':pp,'timept':tt,'charge':zz,'rep':rr}
                        raw_filter_idx += list(filter_df(rawdata,**spec).index)
                    rawdata = pd.concat([rawdata.drop(index=raw_filter_idx),userraw],ignore_index=True)
                else: rawdata = userraw

        if config.Preset_Pops:
            config.Max_Pops = temp_max_pops #put back to user set value   

        ## Now have deutdata, rawdata from any data format 
        
        if deutdata.empty:
            print("No intensity data for "+str(sample)+' peptide '+peptide_range+' z= '+str(charge))
            continue
        if deutdata.Intensity.sum() == 0: 
            print("No intensity data for "+str(sample)+' peptide '+peptide_range+' z= '+str(charge))
            continue
               
        time_points = sorted(set(deutdata.time))
        n_time_points = len(time_points)
        time_indexes = sorted(set(deutdata.time_idx)) 
        print("Found time points (s): "+', '.join(map(str, time_points)))
        if 'rep' not in deutdata.columns:
            deutdata['rep'] = 1
        if 'rep' not in rawdata.columns:
            rawdata['rep'] = 1
        max_time_reps = int(sorted(set(deutdata.rep))[-1])
        #config.Current_Isotope = get_na_isotope(peptide,charge,mod_dict=mod_dic)
        if 'NA_envelope' in row.keys():
            user_env = row['NA_envelope']
        else: user_env = None
        config.Current_Isotope = choose_na(peptide,charge,mod_dict=mod_dic,deutdata=deutdata,user_env=user_env)

        if GP:
            dax_legend_elements = []
            for irep in range(max_time_reps):
                dax_legend_elements += [ mlines.Line2D([0],[0], color='w',markerfacecolor = mpl_colors[irep],
                                                    marker='o',label="rep"+str(irep+1),markersize=10) ]
            if config.Test_Data:
                dax_legend_elements = [ mlines.Line2D([0],[0], color='w',markerfacecolor = mpl_colors[0],
                                                    marker='v',label="Fit Data",markersize=10),
                                        mlines.Line2D([0],[0], color='w',markerfacecolor = None,markeredgecolor='darkred',
                                                    marker='o',label="Solution",markersize=10) ]
        # print("\nDataset",dataset_count,"of",len(metadf))
        # print("Performing fits for "+sample+" "+peptide_range+": "+peptide+" z="+str(int(charge)))
        #print(n_time_points, max_time_reps) 

        nrows = n_time_points
        if max_time_reps == 1: overlay_reps = False
        if overlay_reps:
            ncols = max_time_reps + 1
        else: ncols = max_time_reps
        
        if GP:
            #figsize is width, height
            fig, ax = plt.subplots(figsize=(ncols*5+2, nrows*5), ncols=ncols, nrows = nrows, squeeze=False,num=1,clear=True)
            fig2, ax2 = plt.subplots(figsize=(ncols*5+2, nrows*5), ncols=ncols, nrows = nrows, squeeze=False,num=2,clear=True)
            if config.Test_Data: dfig,dax=plt.subplots(figsize=(nrows/3+9,6),num=3,clear=True) # numD vs time plot
            else: dfig,dax=plt.subplots(figsize=(9,6),num=3,clear=True)

        #time points are rows i, reps are columns j
        
        ## correction would be better based on n_curves = 1 fits of UN/TD
        ## get corrected deut values from centroids (not fit data) of Un and FullDeut
        ## Need these to compare to the 'solution' values
        n_amides = count_amides(peptide,count_sc=0.0)
        try: namide_scale = config.Nex_Max_Scale
        except: namide_scale = 1.2 #instead of using exchangable sidechain fraction to set max_n_amides, just use scale factor
        max_n_amides = count_amides(peptide,count_sc=0.0)*namide_scale #### Might be better way to set this

        if (max_n_amides < 1): 
            print(f"No exchangable amides found\n")
            continue # go to the next metadf row 
        
        undeut_mz = 0.0
        if 'ion_type' in mod_dic.keys(): 
            ion_type = mod_dic['ion_type']
            #mod_dic.pop('ion_type')
        else: ion_type = 'M'
        if len(peptide)>0:       
            undeut_mz = mass.calculate_mass(sequence=peptide,show_unmodified_termini=True,charge=charge,ion_type=ion_type) ####Adjust for mod_dic
        if mod_dic:
            mod_comp = {}
            if 'Hex' in mod_dic.keys():
                n_amides += mod_dic['Hex']
                max_n_amides += mod_dic['Hex']*namide_scale
            comp_keys = set(list(mod_dic.keys())) - set(['Hex']) - set(['ion_type']) - set(['undeut_mz'])
            if len(comp_keys) > 0:
                for k in comp_keys:
                    mod_comp[k] = mod_dic[k]
                undeut_mz += mass.calculate_mass(composition=mod_comp,charge=charge,ion_type=ion_type)

        Noise = 0.0
        d_corr = 1.0        

        if all(tp in time_points for tp in [0,config.FullDeut_Time]): #[0,1e6]
            time_reps = deutdata.rep[(deutdata.time==0.0)].unique().astype('int')
            n_time_reps = len(time_reps)
            centroidUD_all = []
            for r in time_reps:
                focal_data = deutdata.copy()[(deutdata.time==0.0) & (deutdata.rep==r)]
                mz=np.array(focal_data.mz.copy())
                y=np.array(focal_data.Intensity.copy())
                Noise += max(y)/n_time_reps
                if config.Binomial_dCorr: 
                    #use the binomial center for the d_corr calc, not the centroid which may capture impurity peaks
                    initial_estimate, bounds = init_params(1,max_n_amides,seed=None)#config.Random_Seed-1)
                    p0_UD = (0, 0.05, 0.05, 1.0 ) #scaler, nex, mu, frac
                    try:
                        fit, covar = curve_fit( n_fitfunc, len(y)-1, y/np.sum(y), p0=p0_UD, maxfev=int(1e6), 
                                                bounds = bounds   )
                        scaler,nexs,mus,fracs = get_params(*fit,sort=True,norm=True,unpack=True)
                        scaler = np.power( 10.0, scaler )  
                        fit_y = scaler * fitfunc( len(y)-1, nexs[0], mus[0], ) * np.sum(y)
                        #print("UN len(fit_y), sum(fit_y)",len(fit_y),sum(fit_y)) #TROUBLESHOOTING
                        cent_r =  [sum(mz*fit_y)/sum(fit_y)] 
                    except RuntimeError:
                        print (f"failed to fit: UnDeut data for d_corr calculation")
                        break
                else: cent_r = [sum(mz*y)/sum(y)]
                centroidUD_all += cent_r
            centroidUD = np.mean(centroidUD_all)         
            centroidUD_err = np.std(centroidUD_all)

            time_reps = deutdata.rep[(deutdata.time==config.FullDeut_Time)].unique().astype('int')
            n_time_reps = len(time_reps)
            centroidTD_all = []
            for r in time_reps:
                focal_data = deutdata.copy()[(deutdata.time==config.FullDeut_Time) & (deutdata.rep==r)]
                mz=np.array(focal_data.mz.copy())
                y=np.array(focal_data.Intensity.copy())
                Noise += max(y)/n_time_reps
                if config.Binomial_dCorr: 
                    #use the binomial center for the d_corr calc, not the centroid which may capture impurity peaks
                    initial_estimate, bounds = init_params(1,max_n_amides,seed=None)#config.Random_Seed-1)
                    p0_TD = (0, max(max_n_amides - 1,0), 0.8, 1.0 ) #
                    try:
                        fit, covar = curve_fit( n_fitfunc, len(y)-1, y/np.sum(y), p0=p0_TD, maxfev=int(1e6), 
                                                bounds = bounds   )
                        scaler,nexs,mus,fracs = get_params(*fit,sort=True,norm=True,unpack=True)
                        scaler = np.power( 10.0, scaler )  
                        fit_y = scaler * fitfunc( len(y)-1, nexs[0], mus[0], ) * np.sum(y)
                        #print("TD len(fit_y), sum(fit_y)",len(fit_y),sum(fit_y)) #TROUBLESHOOTING
                        cent_r = [sum(mz*fit_y)/sum(fit_y)]                      
                    except RuntimeError:
                        print (f"failed to fit: TD data for d_corr calculation")
                        break
                else: cent_r = [sum(mz*y)/sum(y)]
                centroidTD_all += cent_r
            centroidTD = np.mean(centroidTD_all)
            centroidTD_err = np.std(centroidTD_all)

            d_corr = (charge*(centroidTD - centroidUD)/n_amides)
            Noise = Noise/2.0 * config.Y_ERR/100.0 #noise is config.Y_ERR * avg maxInt of UnDeut and FullDeut
            dax_log = True #safety for plotting ndeut with log scale
        else: 
            print("Missing Un or Full Datasets\nUsing percentage of maxInt for Noise or setNoise if specified")
            dax_log = False
            Noise = config.Y_ERR/100.0 * deutdata.Intensity.max()
            #use pred undeut mz if no UN/TD data
            n_mz = np.arange(len(config.Current_Isotope)+1)
            mz = undeut_mz + (n_mz*1.006227)/charge

            centroidUD = sum(config.Current_Isotope*mz[0:len(config.Current_Isotope)])/sum(config.Current_Isotope)
            centroidTD = centroidUD + n_amides*config.Dfrac/charge
            d_corr = (charge*(centroidTD - centroidUD)/n_amides) # = config.Dfrac #units mass per amide

        
        #print ("centroids and delta mass:",centroidUD, centroidTD,(centroidTD-centroidUD)*charge,"charge, amides, dcorr",charge,n_amides,d_corr)

        if config.setNoise: Noise = config.setNoise  
        # print("namides, centroidUD, centroidTD",n_amides, centroidUD,centroidTD)
        # print("formula dcorr:",(charge*(centroidTD - centroidUD)/n_amides))
        # print("Noise:",Noise)
        # print("dcorr",d_corr)

        for i, ti in enumerate(time_indexes): #range(0,n_time_points): 
            n_time_rep = int(max(deutdata.rep[(deutdata.time_idx==ti)]))
            timept = max(deutdata.time[(deutdata.time_idx==ti)])
            if timept == config.FullDeut_Time: 
                timelabel = 'FullDeut'
            elif timept == 0: timelabel = 'UnDeut'
            else: timelabel = str(timept)+'s'
            if config.Test_Data: timelabel = 'Exp '+str(ti)
            ax_max_y_ov = 0 #for ylim in overlay plot
            ## TODO would like to test n_curves for all reps in time point, then do bootstrap with best_n_curves
            for j in range(1,n_time_rep+1):
                fitparamsdf = pd.DataFrame()

                lowermz = deutdata.mz[deutdata.rep==j].min()
                uppermz = deutdata.mz[deutdata.rep==j].max()
            
                focal_data = deutdata.copy()[(deutdata.time_idx == ti) & (deutdata.rep == j)]
                if (focal_data.empty) | (focal_data.Intensity.sum() == 0): continue
                if not rawdata.empty: 
                    if config.Keep_Raw: focal_raw = rawdata.copy()[(rawdata.time_idx == ti) & (rawdata.rep == j)]
                else: focal_raw = focal_data.copy()
                envelope_height = focal_data['Intensity'].max() * config.Env_threshold
                env, env_Int = get_mz_env(envelope_height,focal_data,pts=True)
                
                mz=np.array(focal_data.mz.copy())
                y=np.array(focal_data.Intensity.copy())
                
                #x=np.full(y.shape,len(y))

                env_symmetry_adj = 2.0 - (y.max() - env_Int)/y.max() # 0 -> assym, 1 -> symm  
                                                            # want 0 to be 2x and 1 to be 1x -> y = -1*x + 2

                max_y = np.max(y) #maximum Intensity, can filter after 
                ynorm_factor = np.sum(y)
                y_norm = y / ynorm_factor # 50secs with vs 2 mins without normalization
                

                if config.Scale_Y_Values: 
                    scale_y = ynorm_factor
                else:
                    if config.Keep_Raw: 
                        #print("duplicates",rawdata[rawdata.index.duplicated()])
                        focal_raw.Intensity = rawdata.Intensity / ynorm_factor  #norm raw instead
                    scale_y = 1.0

                n_bins = len(y)-1
                                                
                centroid_j = sum(mz*y)/sum(y) #centroid from unfit picked peaks

                #data_fit =pd.DataFrame({'time':[timept],'rep':[j],'centroid':[centroid_j]})
                data_fit = pd.DataFrame(columns = data_fit_columns)
                data_fit.loc[0,'data_id'] = index # index of metadf file used for fits 
                data_fit.loc[0,'time'] = timept
                data_fit.loc[0,'time_idx'] = ti
                data_fit.loc[0,'rep'] = j
                data_fit.loc[0,'max_Int'] = max_y
                data_fit.loc[0,'centroid'] = centroid_j
                data_fit.loc[0,'UN_TD_corr'] = d_corr 
                data_fit.loc[0,'sample'] = sample
                data_fit.loc[0,'peptide'] = peptide
                data_fit.loc[0,'peptide_range'] = peptide_range
                data_fit.loc[0,'start_seq'] = int(peptide_range.split('-')[0])
                data_fit.loc[0,'end_seq'] = int(peptide_range.split('-')[1])
                data_fit.loc[0,'charge'] = int(charge)
                
                        
                fstdev=[]     
                
                p_corr = 1.0 
                low_n = config.Min_Pops            
                if config.Limit_by_envelope: high_n = min(max(1,int(env_symmetry_adj * charge * (env[1]-env[0])/(3*config.Env_limit)-2/3)),config.Max_Pops)
                else: high_n = config.Max_Pops
                
                if Preset_Pops:
                    try:
                        fixed_pop = filter_df(config_df,samples=sample,peptides=peptide,charge=charge,rep=j,timept=timept,quiet=True)#['min/max_pops'][0]
                        low_n = int(fixed_pop['min_pops'].values[0])
                        high_n  = int(fixed_pop['max_pops'].values[0])                        
                        #print("Number of fit populations is fixed to",high_n)
                    except: 
                        #print("Failed to set the fixed population. Using,",low_n,"to",high_n,"populations")
                        pass
                high_n = max(low_n,high_n) #safety in case low_n > high_n

                best_n_curves = low_n
                for n_curves in range( low_n, high_n+1 ):  #[scaler] [n *n_curves] [mu *n_curves] [frac * (n_curves )] )]
                    print("Time point:",timelabel,"Rep:",j,"Npops:",n_curves,"          ",end='\r',flush=True) 
                    
                    initial_estimate, bounds = init_params(n_curves,max_n_amides,seed=config.Random_Seed)                 
                    if (len(initial_estimate) > n_bins): 
                        print(f"attempting to fit more parameters than data points: time {timelabel} rep {j} N={n_curves} curves")
                        #should be able to exit at this point, haven't updated last fit parameters including p_err
                        break # exit the for n_curves loop
                    #print("bounds",bounds)
                    try:                        
                        if config.Use_DiffEvo == True:
                            de_bounds = list(zip(bounds[0],bounds[1]))
                            args = (n_bins,y_norm)
                            config.DiffEvo_kwargs.update({'callback':MinimizeStopper,'init':'sobol'})
                            result = differential_evolution(minimize_func, de_bounds, args=args,**config.DiffEvo_kwargs)
                            fit =  get_params(*result.x,sort=True,norm=True,unpack=False)
                            #print(fit)
                            covar = np.zeros_like(fit)
                            fit_y = n_fitfunc( n_bins, *fit )
                            rss = calc_rss( y_norm, fit_y, )
                        else:
                            fit, covar = curve_fit( n_fitfunc, n_bins, y_norm, p0=initial_estimate, maxfev=int(1e6), 
                                                    bounds = bounds   )
                            fit_y = n_fitfunc( n_bins, *fit )
                            rss = calc_rss( y_norm, fit_y, )
                            #bestifit = 1
                            #print("trying fit 1 of ",config.BestFit_of_X," fits:",initial_estimate,rss)
                            for ifits in range(2,config.BestFit_of_X+1):
                                if config.Random_Seed: seed = config.Random_Seed+ifits*29
                                else: seed = None
                                initial_estimate, bounds = init_params(n_curves,max_n_amides,seed=seed)
                                
                                newfit, newcovar = curve_fit( n_fitfunc, n_bins, y_norm, p0=initial_estimate, maxfev=int(1e6), 
                                                    bounds = bounds   )
                                new_fit_y = n_fitfunc( n_bins, *newfit )
                                new_rss = calc_rss( y_norm, new_fit_y, )
                                #print("trying",ifits,"of ",config.BestFit_of_X," fits:",initial_estimate,new_rss)
                                if new_rss < rss:
                                    #bestifit = ifits 
                                    rss = new_rss
                                    fit = newfit
                                    fit_y = new_fit_y
                                    covar = newcovar
                    except RuntimeError:
                        print (f"failed to fit: time {timelabel} rep {j} N={n_curves} curves")
                        break
                    
                    #print("using bestifit number ",bestifit)
                    
                    # Perform statistical test, keep the best model
                    n_params = len( initial_estimate )
                  
                    ###IN PROGRESS - record all initial fit params
                    fitparamsdf[sample_id_cols] = data_fit[sample_id_cols] 
                    fitparamsdf['ncurves'] = n_curves
                    fitparamsdf['nboot'] = 0
                    fitparamsdf['rss'] = rss
                    fit = get_params(*fit,norm=True,unpack=False) #make sure populations are normalized before saving outputs 
                    fitparamsdf['Fit_Params'] = (' ').join(map(str,fit))
                    # scipy curve_fit does not produce same Covar matrix as Excel so this isn't useful info
                    # fstdev = np.sqrt(np.diag(covar))
                    # fitparamsdf['Fit_Var'] = (' ').join(map(str,fstdev))                    
                    # with np.errstate(divide='ignore'): ppar = stats.t.sf(fit/fstdev,n_bins+1-n_params)*2  #p-value from the covariance of linear regression fits
                    # fitparamsdf['Fit_Ppar'] =  (' ').join(map(str,ppar))
                    if config.Test_Data:
                                fitparamsdf['solution_npops'] = solution['npops'][solution.time==timept].to_numpy()[0] 

                    if n_curves == low_n:
                        best_fit = fit 
                        best_covar = covar 
                        # print( timelabel +' '+str(sample)+' '+peptide_range+' z'+str(int(charge))+' rep'+str(j)
                        #         +' N = ' + str(n_curves).ljust(5) + 'p = ' + format( p_corr, '.3e')
                        #         +' rss = '+format(rss,'.3e')+' '+np.array_str(np.array(fit),precision=5).replace('\n',''),file=fout)
                        fitparamsdf['p-value'] = p_corr
                        fitparams_all = pd.concat([fitparams_all,fitparamsdf],ignore_index=True)
                        save_metadf(fitparams_all,filename="fitparamsAll_asrun_"+date+".csv")
                                                
                    else: #if n_curves > low_n:
                        # F = (prev_rss - curr_rss)/(dof_prev - dof_curr) / (rss_curr / dof_curr)
                        #   = (prev_rss - rss)/(n_bins + 1 - (n_params - 3) - (nbins + 1 - n_params)) / (rss / (nbins + 1 - n-params))
                        #                      xnbinsx + x1x - xnparamsx + 3 - xnbinsx - x1x + xnparamsx  = 3 (1 x frac, mu, nex for each additional curve)
                        F = ( ( prev_rss - rss ) / (3)  ) / ( rss / ( n_bins + 1 - n_params ) )
                        p = 1.0 - stats.f.cdf( F, 3, n_bins + 1 - n_params )
                        p_corr = p * (n_curves-1) ### is this correct? this small scaling doesn't really affect the selection
                        
                        fitparamsdf['p-value'] = p_corr
                        fitparams_all = pd.concat([fitparams_all,fitparamsdf],ignore_index=True)
                        save_metadf(fitparams_all,filename="fitparamsAll_asrun_"+date+".csv")
                        
                        if p_corr >= config.Ncurve_p_accept:
                            p_corr = prev_pcorr 
                            break # exit the for n_curves loop; insufficient improvement in fit
                        ## fall back to n-1 curves if one of the populations <  Pop_Thresh
                        _,_,_,frac_check = get_params(*fit,norm=True,unpack=True)
                        if np.min(frac_check) < config.Pop_Thresh:
                            p_corr = prev_pcorr
                            print ("\nmin population below threshold: falling back to",(n_curves-1),"curve(s)")
                            break
                                       
                    prev_rss = rss   #only gets to these if N is better than N-1           
                    prev_pcorr = p_corr
                    best_fit = fit            
                    best_covar = covar 
                    best_n_curves = n_curves
                    if rss < config.Residual_Cutoff:
                        break # if the residual sum squared is sufficiently small, use that number of curves 
                      
                #end n_curves for loop

                data_fit.loc[0,'fit_pops'] = best_n_curves
                data_fit.loc[0,'p-value'] = p_corr
                fit_y = n_fitfunc( n_bins, *best_fit )
                fstdev = np.sqrt(np.diag(best_covar)) ### this error is not realistic and can be over and underestimate of error
                                                    ### Artifically small: fit can be very well converged but still miss the mark
                                                    ### since measurement contains error and fit is assuming perfect data
                                                    ### this is why adding bootstrap noise better samples around the 'solution'
                                                    ### Artifically large: essentially swapping between populations during curve_fit
                #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html see pcov
                # To compute one standard deviation errors on the parameters, use perr = np.sqrt(np.diag(pcov)) //when sigma=None
                # equivalent ppar for HX-Express:
                # ppar = stats.t.sf(best_fit/fstdev,n_bins+1-n_params)*2

                # do the bootstrap fits 
                if config.Nboot:            
                    p0_boot=[]
                    for boot in range(config.Nboot): #send different p0 for each Bootstrap iteration
                        if config.Random_Seed: bseed = boot+config.Random_Seed+1
                        else: bseed = None
                        p0, bbounds = init_params(best_n_curves,max_n_amides,seed=bseed) #best_n_curves
                        de_bounds = list(zip(bbounds[0],bbounds[1]))
                        
                        #randomize p0 from best_fit
                        #if config.Use_DiffEvo == True: 
                        p0 = best_fit 
                        randp0 = rng.normal(0., 0.05, len(p0)) # +/- fraction up to 0.1
                        p0_new = p0 + randp0 * [0.95*(u - l) for (l,u)  in (de_bounds)]
                        p0 = list(np.clip(p0_new, bbounds[0], bbounds[1]))
                        p0 = get_params(*p0,sort=True,norm=True,unpack=False)
                            
                        p0_boot.append(p0)
                    #p0_boot, bbounds = init_params(best_n_curves,n_amides,seed=config.Random_Seed)
                    #rss = 0.0 #calc_rss( y_norm, fit_y, ) #this doesn't seem appropriate for binom fits
                    if GP == False: boot_ax = None
                    else: boot_ax = ax[i,j-1]
                    
                    pfit, boot_rss, boot_centers, pvar = fit_bootstrap(p0_boot,bbounds,n_bins,y_norm,sigma_res=Noise/ynorm_factor, 
                                                nboot=config.Nboot,ax=boot_ax,yscale=scale_y)
                    #pfit = get_params(*pfit,sort=False,norm=True,unpack=False)
                    
                    # plot mus vs fracs for each dataset 
                    nm_alpha = [0.6]*int(config.Nboot) #default values if not scaling
                    nm_marker = [50.0]*int(config.Nboot) #default values if not scaling
                    p_array=np.array(pfit)
                    pvar_array = np.array(pvar) # pvar is already the diagonal of the covariance matrix of the fit
                    ns_array = p_array[:,1:best_n_curves+1]
                    mus_array = p_array[:,best_n_curves+1:best_n_curves*2+1]                
                    fracs_array = p_array[:,-best_n_curves:]
                    #nm_alpha = rankdata([-1*br for br in boot_rss],method='dense')/len(boot_rss)/2.0+0.2
                    
                    for pb in range(len(boot_rss)):
                        fitparamsdf = pd.DataFrame()
                        fitparamsdf[sample_id_cols] = data_fit[sample_id_cols] 
                        fitparamsdf['ncurves'] = best_n_curves
                        fitparamsdf['nboot'] = pb+1
                        fitparamsdf['rss'] = boot_rss[pb]
                        fitparamsdf['Fit_Params'] = (' ').join(map(str,p_array[pb]))
                        # fitparamsdf['Fit_Var'] = (' ').join(map(str,pvar_array[pb]))
                        # with np.errstate(divide='ignore'): ppar = stats.t.sf(p_array[pb]/pvar_array[pb],n_bins+1-len(p_array[pb]))*2  
                        # #p-value from the covariance of linear regression fits
                        # fitparamsdf['Fit_Ppar'] =  (' ').join(map(str,ppar))

                        fitparamsdf['p-value'] = p_corr # p-value from F-test of the n_curves versus n_curves - 1 fit. 
                        if config.Test_Data:
                            fitparamsdf['solution_npops'] = solution['npops'][solution.time==timept].to_numpy()[0]
                        fitparams_all = pd.concat([fitparams_all,fitparamsdf],ignore_index=True)
                    save_metadf(fitparams_all,filename="fitparamsAll_asrun_"+date+".csv")

                    bnm_avg, bnm_err = {}, {}
                    bfrac_avg,bfrac_err = {}, {}
                    for n in range(best_n_curves):
                        boot_centk = (boot_centers[:,n]-centroidUD)*charge*1/d_corr # like nm but discrete due to peak intervals 
                        #print(boot_centk)                            
                        nm = mus_array[:,n]*ns_array[:,n]*1/d_corr #apply correction based on TD-UN #may need to subtract nm_UNdeut?
                        nm_marker = [50.0]*len(nm)
                        nm_alpha = [0.6]*len(nm)                            
                        fracn = fracs_array[:,n]                            

                        #TODO this Pop_Thresh selection is reasonable for a 2 population fit, but is a bit off for 2+
                        # It would be better to analyze the fitparams dataframe afterwards 
                        # could add an ncurves_adj column to average fits of same values together.
                        # 
                        len_nm = len(nm[(fracn > config.Pop_Thresh) & (fracn < 1 - config.Pop_Thresh)])
                        if (best_n_curves == 2) & (len_nm > 0):
                            bnm_avg[n] = nm[(fracn > config.Pop_Thresh) & (fracn < 1 - config.Pop_Thresh)].mean()
                            bnm_err[n] = nm[(fracn > config.Pop_Thresh) & (fracn < 1 - config.Pop_Thresh)].std()
                            bfrac_avg[n] = fracn[(fracn > config.Pop_Thresh) & (fracn < 1 - config.Pop_Thresh)].mean()
                            bfrac_err[n] = fracn[(fracn > config.Pop_Thresh) & (fracn < 1 - config.Pop_Thresh)].std()
                        else: 
                            bnm_avg[n] = nm.mean()
                            bnm_err[n] = nm.std()
                            bfrac_avg[n] = fracn.mean()
                            bfrac_err[n] = fracn.std()
                        if GP:
                            #ax2[i,j-1].scatter(fracs_array[:,n],mus_array[:,n],label='pop'+str(n))
                            ax2[i,j-1].scatter(nm,fracn,label='pop'+str(n),alpha=nm_alpha,c=mpl_colors_light[n],s=nm_marker)
                            #ax2[i,j-1].scatter(boot_centk,fracn,label='pop'+str(n),alpha=0.8,c='purple',marker='x',zorder=0) 
                            ax2[i,j-1].errorbar(bnm_avg[n],bfrac_avg[n],xerr=bnm_err[n],yerr=bfrac_err[n],elinewidth=2,zorder=10,color=mpl_colors_dark[n])
                            if overlay_reps:
                                ax2[i,ncols-1].scatter(nm,fracs_array[:,n],alpha=0.6,label=str(j)+'_pop'+str(n),)
                                ax2[i,ncols-1].errorbar(bnm_avg[n],bfrac_avg[n],xerr=bnm_err[n],yerr=bfrac_err[n],elinewidth=2,zorder=0)
                            if config.Test_Data: #plot against time_set not the fake timept value
                                dax.errorbar(i,bnm_avg[n],yerr=bnm_err[n],alpha=1.0,color=mpl_colors[j-1])
                                dax.scatter(i,bnm_avg[n],marker='v',alpha=1.0,s=100.0*bfrac_avg[n]+20.0,
                                            color=mpl_colors[j-1],edgecolors=mpl_colors_dark[j-1])
                            else: 
                                dax.errorbar(timept,bnm_avg[n],yerr=bnm_err[n],alpha=1.0,color=mpl_colors[j-1])
                                dax.scatter(timept,bnm_avg[n],marker='v',alpha=1.0,s=50.0*bfrac_avg[n]+5.0,
                                            color=mpl_colors[j-1],edgecolors=mpl_colors_dark[j-1])
                if config.Test_Data:
                    for k in range(3):
                        s_nm = solution['fd'+str(k+1)][solution.time==timept].to_numpy()[0]/100*n_amides #no d_corr required
                        s_f = solution['p'+str(k+1)][solution.time==timept].to_numpy()[0]
                        if GP:
                            ax2[i,j-1].scatter(s_nm,s_f,edgecolors='darkred',alpha=1.0,marker='o',facecolor='none',s=60,linewidth=1.5)
                            dax.scatter(i,s_nm,edgecolors='darkred',alpha=1.0,marker='o',facecolor='none',s=100.0*s_f+20.0,)
                    data_fit.loc[0,'solution_npops'] = solution['npops'][solution.time==timept].to_numpy()[0]

                ### PLOTS ###
                y_plot = y_norm * scale_y
                fit_y = fit_y * scale_y
                scaled_env_height = max(y_plot)*config.Env_threshold
                env_resolution = env_symmetry_adj *charge * (env[1]-env[0]) / ( len(best_fit) + 1.0 ) #rough measure of whether there's enough information to do the n_curves fit
                
                #data_fit.loc[0,'env_res_1'] = env_symmetry_adj *charge * (env[1]-env[0]) / ( 2.0 ) #env_res for n_curves = 1 case 
                data_fit.loc[0,'env_width'] = charge*(env[1]-env[0])
                data_fit.loc[0,'env_symm'] = env_symmetry_adj
                #data_fit.loc[0,'skewness'] = skew(y_norm,bias=False)
                data_fit.loc[0,'max_namides']=n_amides

                if GP:
                    ax_max_y = max(max(y_plot),max(fit_y)) + Noise/ynorm_factor*scale_y #find max value on plot
                    
                    
                    env_label = "Env res: "+format(env_resolution,'0.2f')#+"/"+format(env_dof,'0.2f')
                    ax[i,j-1].plot(env,[scaled_env_height,scaled_env_height],label=env_label,color='darkorange')
                    if config.Keep_Raw:
                        ax_max_rawy = focal_raw[focal_raw['mz'].between(lowermz-3/charge,uppermz+9/charge)]['Intensity'].max() 
                        ax_max_y = max(ax_max_rawy,ax_max_y)
                        ax[i,j-1].plot( focal_raw.mz, focal_raw.Intensity, color='#999999' ) #ax[i,j-1]
                    ax[i,j-1].vlines( mz, 0.0, y_plot, color='#999999' ) #ax[i,j-1] # TROUBLESHOOTING #else: 
                    
                    ax[i,j-1].plot( mz, y_plot, 'ro', label='data '+ timelabel+", rep "+str(j), markersize='4')
                    ax[i,j-1].vlines( centroid_j, 0, max(y_plot), color='orange' ,label='m/z = '+format(centroid_j,'.2f'),linestyles='dashed',linewidth=2,zorder=10)

                    
                    ax[i,j-1].plot( mz, fit_y, '-', label='fit sum N='+str(best_n_curves))

                    if overlay_reps:
                        ax_max_y_ov = max(ax_max_y_ov,ax_max_y)
                        if config.Keep_Raw:
                            ax[i,ncols-1].plot( focal_raw.mz, focal_raw.Intensity, color=mpl_colors[j-1],alpha=0.5 ) 
                        else: ax[i,ncols-1].vlines( mz, 0.0, y_plot, color=mpl_colors[j-1], alpha=0.5 )
                        ax[i,ncols-1].plot( mz, y_plot, 'o',color=mpl_colors_dark[j-1], label="rep "+str(j), markersize='4')
                        ax[i,ncols-1].plot( mz, fit_y, '-',color=mpl_colors[j-1])# label='fit sum N='+str(best_n_curves))

                scaler,nexs,mus,fracs = get_params(*best_fit,sort=True,norm=False,unpack=True)
                scaler = np.power( 10.0, scaler )   
                fracsum = np.sum(fracs)

                #scaler at 0, n's at 1:best_n_curves+1, mu's at best_n_curves+1:2*best_ncurves+1, fracs at 2*best_n_curves+1:
                for k in range( best_n_curves ):
                    nex = nexs[k]
                    nex_err = fstdev[k+1]
                    mu = mus[k]
                    mu_err = fstdev[k+best_n_curves+1]
                    frac = fracs[k]/fracsum
                    frac_err = 1.0 #### issue with diffevo, min(1.0,fstdev[-1]/fracsum) #if best_n_curves > 1 else 0.
                    fit_yk = scaler * frac * fitfunc( n_bins, nex, mu, ) * scale_y
                    kindcent = sum(mz*fit_yk)/sum(fit_yk)
                    deut_corr_k = (kindcent - centroidUD)*charge * 1/d_corr
                    data_fit.loc[0,'centroid_'+str(k+1)]=kindcent
                    data_fit.loc[0,'Dabs_'+str(k+1)]=bnm_avg[k] if config.Nboot else deut_corr_k 
                    data_fit.loc[0,'Dabs_std_'+str(k+1)]=bnm_err[k] if config.Nboot else np.nan ##
                    data_fit.loc[0,'pop_'+str(k+1)]= bfrac_avg[k] if config.Nboot else frac
                    data_fit.loc[0,'pop_std_'+str(k+1)]=bfrac_err[k] if config.Nboot else frac_err
                    if GP: 
                        ax2[i,j-1].scatter(deut_corr_k,frac,marker='x',alpha=1.0,c='k',zorder=10)#mpl_colors[k])
                        #dax.errorbar(timept,bnm_avg[n],yerr=bnm_err[n],alpha=1.0,color=mpl_colors[j-1])
                        if config.Nboot == 0:
                            if config.Test_Data:
                                dax.scatter(i,deut_corr_k,marker='x',alpha=1.0,color=mpl_colors_dark[j-1])
                            else:
                                dax.scatter(timept,deut_corr_k,marker='x',alpha=1.0,color=mpl_colors_dark[j-1])
                        if overlay_reps:
                            ax2[i,ncols-1].scatter(deut_corr_k,frac,marker='x',alpha=1.0,c='k',zorder=10)#mpl_colors[k])

                        #data_fit.loc[0,'mu_'+str(k+1)]=mu  #don't include in fit report, keep with all params output
                        #data_fit.loc[0,'Nex_'+str(k+1)]=nex #don't include in fit report, keep with all params output
                        #data_fit.loc[0,'Nex_std_'+str(k+1)]=nex_err #don't include in fit report, keep with all params output

                        
                        nm_for_label = bnm_avg[k] if config.Nboot else deut_corr_k
                        frac_for_label = bfrac_avg[k] if config.Nboot else frac 
                        #print("nm:",bnm_avg[k],"deut_corr_k",deut_corr_k)
                        plot_label = ('pop'+str(k+1)+' = '+format(frac_for_label,'.2f')#+' ± '+format(frac_err,'.2f') 
                                            +'\nDabs_'+str(k+1)+' = '+format(nm_for_label,'.1f')#+' ± '+format(nex_err,'.2f')
                                            +'\nm/z = '+format(kindcent,'.2f'))
                        ax[i,j-1].plot( mz, fit_yk, color = 'black', linestyle='solid',linewidth=1.,label=plot_label)
                        ax[i,j-1].set(xlim=(lowermz-3/charge,uppermz+9/charge)) #gives rightside space for legend
                        try: ax[i,j-1].set(ylim=(-ax_max_y * 0.1, ax_max_y*1.1))
                        except: pass #issue with limits

                if GP:
                    ax[i,j-1].set_title(label=str(sample)+': '+peptide_range+" "+str(shortpeptide)+" z="+str(int(charge)),loc='center')
                    ax[i,j-1].legend(frameon=False,loc='upper right');
                    ax[i,j-1].set_xlabel("m/z")
                    ax[i,j-1].set_ylabel("Intensity")
                    
                    if overlay_reps: 
                        ax[i,ncols-1].set_title(label='Replicates Overlay')
                        ax[i,ncols-1].set(xlim=(lowermz-3/charge,uppermz+9/charge))
                        ax[i,ncols-1].set(ylim=(-ax_max_y_ov*0.1,ax_max_y_ov*1.1)) 
                        ax[i,ncols-1].legend(frameon=False,loc='upper right',title='data '+ timelabel);
                        ax[i,ncols-1].set_xlabel("m/z")
                        ax[i,ncols-1].set_ylabel("Intensity")
                #data_fit.loc[0,'scaler']=scaler #don't include in fit report, keep with all params output
                                
                data_fits = pd.concat([data_fits,data_fit],ignore_index = True).reset_index(drop = True)
                
                ddf = data_fits.copy()[data_fits['data_id']==index]
                polymodal = []
                polystr = ""
                polymodal = ddf[ddf['fit_pops']>1]['time'].unique().tolist()
                polystr = ' '.join(map(str,polymodal))
                reportdf.loc[index,'polymodal'] = polystr
                if any(tp in polymodal for tp in [0,config.FullDeut_Time]):
                    reportdf.loc[index,'comment'] = "Warning: UN/TD polymodal"
                reportdf.loc[index,'dataset_run'] = "Yes"
                if 'NA_envelope' in reportdf.columns:
                    reportdf.loc[index,'NA_envelope_calc'] = '['+(' ').join(map(str,config.Current_Isotope))+']'
                save_metadf(reportdf,filename="metadf_asrun_"+date+".csv")
                try:
                    data_fit.to_csv(data_output_file_inprogress,mode='a',index_label='Index',header=False) 
                except:
                    if (dataset_count == 1) and (i==0):
                        print("Could not save running results to",data_output_file_inprogress)
                
                centroid_j_corr = (centroid_j - centroidUD) * charge * 1/d_corr # = Dcorr = (m - m0)/(m100-m0)
                if GP:
                    ax2[i,j-1].scatter(centroid_j_corr,1.0,label='Centroid',alpha=1.0,c='darkorange',marker='*',zorder=0)
                    ax2[i,j-1].vlines(ymin=(0.0,0.0) ,ymax=(1.0,1.0) ,x=(0,n_amides),alpha=0.5,linestyles='dotted',color='grey')
                    ax2[i,j-1].set_title(label=str(sample)+': '+peptide_range+" "+str(shortpeptide)+" z="+str(int(charge)),loc='center')
                    ax2[i,j-1].legend(title=timelabel+", rep"+str(j),frameon=True,loc='upper right');
                    ax2[i,j-1].set(xlim=(-1.,max_n_amides+4),ylim=(-0.05,1.05))
                    #ax2[i,j-1].set(xlim=(-3,(uppermz-lowermz+9)))
                    ax2[i,j-1].set_xlabel("Deuterium Level (Da)") #N*mu")
                    ax2[i,j-1].set_ylabel("population") 
                    if overlay_reps:
                        ax2[i,ncols-1].scatter(centroid_j_corr,1.0,label='Centroid',alpha=0.8,c='k',marker='x',zorder=0)
                        ax2[i,ncols-1].vlines(ymin=(0.0,0.0) ,ymax=(1.0,1.0) ,x=(0,n_amides),alpha=0.5,linestyles='dotted',color='grey')
                        ax2[i,ncols-1].set(xlim=(-1.,max_n_amides+4),ylim=(-0.05,1.05))
                        #ax2[i,ncols-1].set(xlim=(-3,(uppermz-lowermz+9)))
                        ax2[i,ncols-1].legend(title=timelabel,frameon=True,loc='upper right');
                        ax2[i,ncols-1].set_title(label='Replicates Overlay')
                        ax2[i,ncols-1].set_xlabel("Deuterium Level (Da)")
                        ax2[i,ncols-1].set_ylabel("population") 

        if GP:
            dax.set_ylabel("Deuterium Level (Da)")
            dax.set_title(label=str(sample)+': '+peptide_range+" "+str(shortpeptide)+" z="+str(int(charge)),loc='center')
            dfig.tight_layout()
            if config.Test_Data: 
                yds = np.unique(solution[['fd1','fd2','fd3']].dropna().values)*n_amides/100.0
                dax.hlines(xmin=[-1.0] * len(yds),xmax=[n_time_points+1.0] * len(yds),y=yds,alpha=0.5,linestyles='dotted',color='grey')
                dax.set(xlim=(-0.5,n_time_points+0.5))
                dax.set_xlabel("Test Experiment")
                dax.legend(handles=dax_legend_elements,loc='upper left')
            else:
                dax.hlines(xmin=(-1.0,-1.0) ,xmax=(config.FullDeut_Time,config.FullDeut_Time) ,y=(0,n_amides),alpha=0.5,linestyles='dotted',color='grey')
                dax.set_xlabel("Time, s")
                if dax_log==True: dax.set_xscale('log')
                else: 
                    if len(time_points) > 1:
                        dax.set(xlim=(min(min(time_points),-max(time_points)*0.1),1.1*max(time_points)))
                dax.legend(handles=dax_legend_elements[0:max_time_reps],loc='lower right')
            fig.tight_layout()
            if config.Nboot: fig2.tight_layout()
            
            try:
                figfile = 'hxex_'+str(index)+'_'+sample+'_'+peptide_range+'_z'+str(int(charge))+'_IndFits_'+date
                print("saving figure as ",figfile)
                fig.savefig(os.path.join(config.Output_DIR,figfile+'.pdf'),format='pdf',dpi=600)
                if config.SVG: fig.savefig(os.path.join(config.Output_DIR,figfile+'.svg'),format='svg')
                if config.Hide_Figure_Output: plt.close(fig)
                else: plt.show()
            except IOError as e:
                print (f"Could not save: {figfile} file is open")   

            try:
                #fig2file = 'hdx_ms_hxex3_'+sample+peptide_range+fitfunc_name+'_p'+format(int(config.Ncurve_p_accept*100),'02d')+'_BootFits_'+yerr_name+date
                fig2file = 'hxex_'+str(index)+'_'+sample+'_'+peptide_range+'_z'+str(int(charge))+'_BootFits_'+date #details will be in separate file
                print("saving figure as ",fig2file)
                fig2.savefig(os.path.join(config.Output_DIR,fig2file+'.pdf'),format='pdf',dpi=600)
                if config.SVG: fig2.savefig(os.path.join(config.Output_DIR,fig2file+'.svg'),format='svg')
                if config.Hide_Figure_Output: plt.close(fig2)
                else: plt.show()
            except IOError as e:
                print (f"Could not save: {fig2file} file is open") 

            try:
                dfigfile = 'hxex_'+str(index)+'_'+sample+'_'+peptide_range+'_z'+str(int(charge))+'_ndeutBoot_'+date
                print("saving figure as ",dfigfile)
                dfig.savefig(os.path.join(config.Output_DIR,dfigfile+'.pdf'),format='pdf',dpi=600)
                if config.SVG: dfig.savefig(os.path.join(config.Output_DIR,dfigfile+'.svg'),format='svg')
                if config.Hide_Figure_Output: plt.close(dfig)
                else: plt.show()
            except IOError as e:
                print (f"Could not save: {dfigfile} file is open") 
            
        deutdata_all = pd.concat([deutdata_all,deutdata])
        rawdata_all = pd.concat([rawdata_all, rawdata])
    
    try: 
        data_output_file = os.path.join(config.Output_DIR,"data_fits"+date+".csv")
        print("Saving results table to",data_output_file)
        data_fits.to_csv(data_output_file,index_label='Index')
    except:
        print("Could not save results file")

    if config.Save_Spectra:
        try:
            print("Saving picked peaks as alldeutdata_ and Raw spectral data as allrawdata_")
            deutdata_all.to_csv(os.path.join(config.Output_DIR,'alldeutdata_'+date+'.csv'),index_label='Index')
            rawdata_all.to_csv(os.path.join(config.Output_DIR,'allrawdata_'+date+'.zip'),index_label='Index',compression={'method': 'zip', 'compresslevel': 9})
        except: 
            print("Could not save deutdata and rawdata to file")

    plt.close("all")
    #fout.close()

    return

def get_data(metadf):
    '''
    all_deutdata, all_rawdata = get_data(metadf)
    '''
    all_rawdata = pd.DataFrame(None) 
    all_deutdata = pd.DataFrame(None) 

    dataset_count = 0
    for index, row in metadf.iterrows():
        dataset_count += 1
        rawdata = pd.DataFrame(None) 
        deutdata = pd.DataFrame(None)

        mod_dic = {}
        if 'modification' in row.keys():
            mod = row['modification']
            mod = mod.split() #'H:1 Hex:1 ion_type:c'
            mod_dic = {x.split(':')[0]:(x.split(':')[1]) for x in mod}
            for k,v in mod_dic.items():
                if k == 'undeut_mz':
                    mod_dic[k] = float(v)                    
                elif k != 'ion_type':
                    mod_dic[k] = int(v)

        hdx_file = row['file']
        sample = row['sample']
        peptide = row['peptide']
        charge = row['charge']

        if config.Data_Type == 1:
            if config.Test_Data:
                deutdata, rawdata, solution = read_hexpress_data(hdx_file,row,keep_raw=True)
            else:
                deutdata, rawdata = read_hexpress_data(hdx_file,row,keep_raw=True,mod_dict=mod_dic)
        else: # Data_type == 2, already checked that it is 1 or 2
            spec_path = os.path.join(config.Data_DIR,sample,hdx_file)
            csv_files = [ f for f in os.listdir(spec_path) if f[-5:]==str(int(charge))+'.csv'  ]
            deutdata, rawdata = read_specexport_data(csv_files,spec_path,row,keep_raw=True,mod_dict=mod_dic)
        all_rawdata = pd.concat([all_rawdata,rawdata],ignore_index=True)
        all_deutdata = pd.concat([all_deutdata,deutdata],ignore_index=True)
    return all_deutdata, all_rawdata 

def export_to_hxexpress(rawdata,metadf,save_xls = False, removeUNTDreps = False):
    ''' 
    hxcols = export_to_hxexpress(all_rawdata,metadf, save_xls=True,removeUNTDreps=True)

    HX-Express format, use removeUNTDreps to remove undeut and TD replicates, not used in HX-Express
    [mz,Int] for each delay

    undeut          TD                  1 min               1.01 min
    456.347     98  456.347     3.5     456.347     28.25   456.347     29.10
    '''
    def time_col(row):
        if row['time'] == 0.0:
            if removeUNTDreps == True:
                if row['rep'] == 1: 
                    return "0.0_undeut"
                else: return "unused"
            x = "undeut"+"{}".format(int(row['rep']))
            order = '0.0'
        elif row['time'] == config.FullDeut_Time:
            if removeUNTDreps == True:
                if row['rep'] == 1: 
                    return "1.0_TD"
                else: return "unused"            
            x =  "TD"+"{}".format(int(row['rep']))
            order = '1.0'
        else:
            x = "{:0.2f}".format(row['time']+(row['rep'] - 1.0)/100.0) + ' sec'
            order = row['time_idx'] + 1 #offset to fit TD
        return str(order)+".{}".format(int(row['rep']))+'_'+x
    temp_raw = rawdata.copy()
    temp_raw['time_col'] = temp_raw.apply(time_col,axis=1)
    time_cols = list(temp_raw['time_col'].unique())
    time_cols.remove('unused')

    time_cols_dict = {}
    for tc in time_cols:
        i = float(tc.split('_')[0])
        time_cols_dict[i] = tc.split('_')[1]
    tc = list(time_cols_dict.keys())
    tc.sort()

    hxcols = pd.DataFrame()
    for i,k in enumerate(tc):
        hxcol = pd.DataFrame()
        v = time_cols_dict[k]
        x = temp_raw.mz[temp_raw['time_col']==str(k)+'_'+v].values
        y = temp_raw.Intensity[temp_raw['time_col']==str(k)+'_'+v].values
        hxcol[v] = x
        hxcol[' '*i] = y
        hxcols = pd.concat([hxcols,hxcol],axis=1,ignore_index=False)

    if save_xls: 
        sample = temp_raw['sample'].unique()[0]
        peprange = temp_raw['peptide_range'].unique()[0]
        pep = temp_raw['peptide'].unique()[0]
        z = str(int(temp_raw['charge'].unique()[0]))
        filename = sample+peprange+'-'+pep+'z'+z
        hxcols.to_excel(os.path.join(config.Data_DIR,filename+".xlsx"),index=None)
    return hxcols

def plot_spectrum(deutdata=pd.DataFrame(),rawdata=pd.DataFrame(),fit_params=pd.DataFrame(),norm=False,residual=False,
                  ax=None,rax=None,plt_kwargs={},rax_kwargs={},simfit=False,saveas=None,return_fit_data=False,user_env=None,mod_dict={}):
    '''
    intended to plot a single spectrum: raw data, picked peaks, and fits
    '''
    n_fitfunc = n_binom_isotope # n_binomials #
    fitfunc = binom_isotope # binom #
    #global Current_Isotope

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if (rax is None) and (residual == True):
            plt.close()
            fig, (ax,rax) =  plt.subplots(2, 1, figsize=(5,5),gridspec_kw={'height_ratios': [3, 1]})
    plt_raw = {'color':'#999999'}
    plt_raw.update(plt_kwargs)
    sample_id_cols = ['sample','peptide_range','time','rep','charge']
    deut_spectra = zip(*[deutdata[col] for col in sample_id_cols])
    deut_spectra = list(dict.fromkeys(deut_spectra))
    all_fit_data = pd.DataFrame()

    for s,p,t,r,z in deut_spectra:
        focal_data = filter_df(deutdata,quiet=True,samples=s,peptide_ranges=p,timept=t,charge=z,rep=r)
        focal_raw = filter_df(rawdata,quiet=True,samples=s,peptide_ranges=p,timept=t,charge=z,rep=r)
        focal_fit = filter_df(fit_params,quiet=True,samples=s,peptide_ranges=p,timept=t,charge=z,rep=r)

        fit_data = pd.DataFrame()

        peptide = focal_data['peptide'].values[0]
        peptide_range = focal_data['peptide_range'].values[0]
        shortpeptide = peptide if len(peptide) < 22 else peptide[:9]+'...'+peptide[-9:] #for tidy plot labels

        mz=np.array(focal_data.mz.copy())        
        if len(mz) == 0: continue
        y=np.array(focal_data.Intensity.copy())
        scale_y = np.sum(y)
        y_norm = y/scale_y
        rawint = focal_raw.Intensity
        scale_raw = rawint/scale_y
        if (len(deut_spectra) > 1): norm = True
        if (norm==True):
            y = y_norm
            rawint = scale_raw
            scale_y = 1.0
        n_bins = len(y)-1
                
        lowermz = min(mz)
        uppermz = max(mz)
        
        if len(deut_spectra)>1: pl = (',').join(map(str,[s,p,t,int(r),int(z)])) 
        elif config.Test_Data == True: pl='Exp ' + str(int(t/60)) #clunky conversion w/o grabbing time_idx
        else: pl='data '+str(t)+'s, rep '+str(r)
        ax.plot(mz, y, 'ro',  markersize='4',label=pl ,zorder=1)
        ax.plot(focal_raw.mz, rawint, **plt_raw ,zorder=0)

        fit_data = pd.DataFrame({'mz':mz,'Intensity':y})

        if not focal_fit.empty:
            fits = focal_fit.copy()
            #config.Current_Isotope = get_na_isotope(peptide,z,npeaks=None,mod_dict={})
            config.Current_Isotope = choose_na(peptide,z,npeaks=None,mod_dict=mod_dict,deutdata=deutdata,user_env=user_env)
            nboot_list = list(fits['nboot'].unique())
            if min(nboot_list) > 0:
                best_n_curves = fits[fits['nboot']==1]['ncurves'].values[0]
            else: best_n_curves = fits['ncurves'].values[-1]
            for nb in nboot_list:
                params_best_fit = fits.copy()[(fits['nboot']==nb) & (fits['ncurves']==best_n_curves)]['Fit_Params'].values[0]
                params_best_fit = [float(x) for x in params_best_fit.split()]
                scaler,nexs,mus,fracs = get_params(*params_best_fit,sort=True,norm=True,unpack=True)
                
                fit_y = n_fitfunc( n_bins, *params_best_fit ) * scale_y 
                if simfit: 
                    sim_scale_y = np.sum(y)/np.sum(fit_y)
                    fit_y = fit_y * sim_scale_y
                #### Issue with scale_y comparing actual simulated parameters, data ends up being y**2 ??
                # print("sum y: ", np.sum(y))
                # print("sumfit: ",np.sum(fit_y))
                # print("y/fit_y:", y/fit_y)

                ax.plot( mz, fit_y, '-', alpha=0.5)#label='fit sum N='+str(best_n_curves));
                fit_data['fit_sum_nb'+str(nb)] = fit_y

                if residual == True:
                    rlabel=''
                    rcolor = 'green'
                    y_residual = y - fit_y
                    if nb == 0:
                        rss = calc_rss(y,fit_y)
                        rlabel = 'rss = '+"{:.3e}".format(rss)
                        rcolor = 'blue'
                    rax.plot( mz, y_residual, '-', alpha=0.5,label=rlabel,color=rcolor)
                
                fracsum = np.sum(fracs)

                #scaler at 0, n's at 1:best_n_curves+1, mu's at best_n_curves+1:2*best_ncurves+1, fracs at 2*best_n_curves+1:
                for k in range( best_n_curves ):
                    nex = nexs[k]
                    mu = mus[k]                
                    frac = fracs[k]/fracsum                                    
                    fit_yk = np.power( 10.0, scaler )  * frac * fitfunc( n_bins, nex, mu, ) * scale_y
                    if simfit: fit_yk = fit_yk * sim_scale_y
                    kindcent = sum(mz*fit_yk)/sum(fit_yk)
                    plot_label = ('pop'+str(k+1)+' = '+format(frac,'.2f') 
                                        +'\nm/z = '+format(kindcent,'.2f'))
                    if nb==0: fit_kwds = {'color':'k','linestyle':'-','linewidth':1,'zorder':10}
                    else: fit_kwds = {'color':'green','linestyle':'dashed','linewidth':3,'alpha':0.2,'zorder':5}
                    if len(deut_spectra) == 1:
                        if (k==0) & (nb==1): fit_kwds.update({'label':'BootFits'})
                        if (nb==0): fit_kwds.update({'label':plot_label})
                    ax.plot( mz, fit_yk,**fit_kwds)
                    fit_data['fit_p'+str(k)+'_nb'+str(nb)] = fit_yk
        fit_data[sample_id_cols] =s,p,t,r,z
        all_fit_data = pd.concat([all_fit_data,fit_data],ignore_index=True)
    if len(deut_spectra) == 1:
        ax.set(xlim=(lowermz-3/z,uppermz+9/z))
        if residual: 
            rax.set_xlabel("m/z")
            rax.set_ylabel("residual")
            rax.set(xlim=(lowermz-3/z,uppermz+9/z))
            rax.set(ylim=(rax_kwargs.get('ymin',None),rax_kwargs.get('ymax',None)))
            rax.legend(frameon=False)
            fig.subplots_adjust(hspace=0.3)            
        ax.set_title(label=str(s)+': '+peptide_range+" "+str(shortpeptide)+" z="+str(int(z)),loc='center')
        ax.set_ylabel("Intensity")
    if norm==True: ax.set_ylabel("Normalized Intensity")
    ax.set_xlabel("m/z")
    
    if len(ax.get_legend_handles_labels()[0])>0:
        ax.legend(frameon=False,loc='upper right');
    
    fig.tight_layout()
    if saveas:
        try:
            figfile = saveas
            print("saving figure as ",os.path.join(config.Output_DIR,figfile))
            if figfile[:-4] == '.pdf':
                fig.savefig(os.path.join(config.Output_DIR,figfile),format='pdf',dpi=600)
            elif figfile[:-4] == '.svg':
                fig.savefig(os.path.join(config.Output_DIR,figfile),format='svg')
            else: fig.savefig(os.path.join(config.Output_DIR,figfile))
        except IOError as e:
            print (f"Could not save: {figfile} file is open")  

    if return_fit_data:
        return all_fit_data
    else: return


## Plot metadf as pdf table
# https://stackoverflow.com/questions/33155776/export-pandas-dataframe-into-a-pdf-file-using-python
# modified from answer by Lak
def dataframe_to_pdf(df, filename, numpages=(1, 1), pagesize=(11, 8.5),pagenos = False):
    '''
    #numpages=max(1,((len(metadf)+1)//35))
    dataframe_to_pdf(metadf, os.path.join(config.Output_DIR,'peptides_list_'+date+'.pdf'), numpages=(max(1,((len(metadf)+1)//35)), 1)) #,pagesize=(8.5, 11))
    '''
    def _draw_as_table(df, pagesize,col_widths,idx_width):
        alternating_colors = [['white'] * len(df.columns), ['lightgray'] * len(df.columns)] * len(df)
        alternating_colors = alternating_colors[:len(df)]
        fig, ax = plt.subplots(figsize=pagesize)
        ax.axis('tight')
        ax.axis('off')
        the_table = ax.table(cellText=df.values,
                            rowLabels=df.index.map(('{: >'+str(idx_width+3)+'d}').format),
                            colLabels=df.columns,
                            rowColours=['lightblue']*len(df),
                            colColours=['lightblue']*len(df.columns),
                            cellColours=alternating_colors,
                            colWidths=col_widths,
                            fontsize=18, 
                            loc='center')
        the_table.auto_set_font_size(False)
        the_table.scale(1,1.5)
        #the_table.auto_set_column_width(col=list(range(len(df.columns))))
        return fig  

    with PdfPages(filename) as pdf:
        nh, nv = numpages
        rows_per_page = len(df) // nh
        cols_per_page = len(df.columns) // nv
        col_widths = []
        for col in df.columns:
            col_widths += [df[col].astype('str').str.len().max()]
        header_widths = df.columns.str.len().to_numpy()
        col_widths = np.max(np.vstack((col_widths,header_widths)),axis=0) 
        col_widths = [max(x /sum(col_widths) * pagesize[0] *0.15,0.1) for x in col_widths] 
                    # frac_len * page_width *scaler
        idx_width = df.index.astype('str').str.len().max()
        for i in range(0, nh):
            for j in range(0, nv):
                page = df.iloc[(i*rows_per_page):min((i+1)*rows_per_page, len(df)),
                            (j*cols_per_page):min((j+1)*cols_per_page, len(df.columns))]
                fig = _draw_as_table(page, pagesize,col_widths,idx_width)
                if (nh > 1 or nv > 1) and pagenos:
                    # Add a part/page number at bottom-center of page
                    fig.text(0.5, 0.5/pagesize[0],
                            "Part-{}x{}: Page-{}".format(i+1, j+1, i*nv + j + 1),
                            ha='center', fontsize=8)
                pdf.savefig(fig, bbox_inches='tight')            
                plt.close()
    return



# ##############################################################################
# '''begin user input'''  these are in config.py, referenced as config.variable 
# ###############################################################################
# import os
# from datetime import datetime
# now = datetime.now()
# date = now.strftime("%d%b%Y")

# WRITE_PARAMS = True #save the params to hdxms_params_$.py file in Data_DIR, can then be read in as PARAMS_FILE 
# Allow_Overwrite = True #don't create a new filename if file already exists
# Save_Spectra = False
# Read_Spectra_List = False # Specify files to be run in a file, includes peptide/charge info. See example files.
#                 # To use this, Recommend setting to False to create and write 'metadf' to file with all availble datasets
#                 # then remove unwanted datasets from the file, and read it in with Read_Spectra_List = True
#                 # and Metadf_File set to the appropriate filename 
# Test_Data = False
# Data_Type = 2  #Test_Data = True has precedence over this, will make Data_Type = 1
#     #1: 'xlxs' , each file contains all timepoint reps at Data_DIR/onesample_onepeptide_alltimepts_allreps_onecharge.xlsx
#                 # current recognized format is e.g. HSPB1_B1B5_0001-0011-MTERRVPFSLL-z2-allspectra.xlsx
#                 # <sample name>_<peptide start>-<peptide end>-<peptide>-<zcharge>-<unused label>.xlsx
#                 # allows for replicats of UnDeut/TD even though HX-Express xlsm does not
#     #2: 'SpecExport', as exported from HDExaminer Data_DIR/samples/peptides/onetimept_onerep_onecharge.csv
#                 # this mode requires a sample.fasta file in the Data_DIR for each sample to be processed, with matching names
# if Data_Type == 1:
#     Data_DIR = 'c:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\ns_HSPB1_Bimodal_Peptide_Data'
#     #Data_DIR = 'C:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\pyHXExpress\\Bimodal_HDX_Data'
#     Metadf_File = "hdx_spectra_list_metadf_02Nov2023.csv" #only used if Read_Spectra_List = True; designates files to process
#     process_ALL = False # if True will assume all .xlsx files are HDX data, use with care
#     User_mutants = ['HSPB1only',]#'HSPB1_B1B6'] #['all'] #first element can be 'all' to include all mutants and/or peptides in directory
#     User_peptides = ['0001-0011',]#'0078-0094']
# if Data_Type == 2:
#     #Data_DIR = '/data/tuttle/HDX-MS/Pearl_SpecExport_30oct2023/SpecExport'
#     #Data_DIR = '/data/tuttle/HDX-MS/Pearl_FimHWTL34K_V6/SpecExport/'
#     Data_DIR = 'c:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\ns_HSPB1_Bimodal_Peptide_Data\\SpecExport'
#     Metadf_File = "hdx_spectra_list_metadf_02Nov2023.csv" #only used if Read_Spectra_List = True; designates files to process
#     process_ALL = True #process_all = True is limited to existing .fasta files, this setting overrides user_ settings
#     User_mutants = [''] #['WT','S19D','S45D','S59D','D3']#['All'] #
#     User_peptides =  ['']#'0049-0054']#['0034-0045'] #['0093-0116'] #['0090-0113']'0122-0166']#

# if Test_Data: 
#     Data_Type = 1
#     #Data_DIR = 'c:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\HX-Express3'
#     Data_DIR = 'C:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\pyHXExpress\\Bimodal_HDX_Data'
#     #Test_Sets = ['v3_Angiotensin_Bimodals.xlsx','v3_GluFib_Bimodals.xlsx']
#     Read_Spectra_List = True
#     Metadf_File = "hdxms_testsets_metadf.csv"
                
# Output_DIR = os.path.join(Data_DIR,'hdxms_analysis_'+str(date),'')
# #if not os.path.exists(Output_DIR): os.makedirs(Output_DIR)

# Preset_Pops = False #Use predetermined number of populations to fit curves, overrides Min/Max_Pops if given
# Preset_Pops_File = os.path.join(Output_DIR,"all_data_results.csv")


# Hide_Figure_Output = True #Recommended when processing lots of data. 
# SVG = False # also save figures as an svg file, slow, but better for making figures 

# Dfrac = 0.9 #fraction deuterated for "fully" deuterated exchange buffer 
# Nterm_subtract = 1 #number of N-term residues to remove from possible exchanging NH's (usually 1 or 2)
# Nboot = 20 # number of individual fits to perform, using n_best_curves from initial round of fits
# setNoise = 200 #if noise value is known, specify instead of estimating as Y_ERR % of avg Un+TD peaks
# Y_ERR = 1.0 #Percent random error applied during boot as y*+np.random.normal(0,yerr), 0.0 for NoNoise, ~0.5% for noise added
#             # the absolute Noise value is then Y_ERR * avg(maxInt of Un and TD)
#             # this is a very rough way to give a consistent Noise value throughout a dataset. 
# Binomial_dCorr = True #use the n=1 binomial fit to determine UN/TD centers for d_corr calculation
#                     # this is essential if dirty spectra; otherwise using sum(mz*y)/sum(y) which will be wrong if there are contaminating peaks
# Peak_Resolution = 50.0 #ppm, sensitivity of peak picker to expected m/z centers 
# Zero_Filling = 4 #number of zero points below noise threshold to include in peak picking
# Env_threshold = 0.1 #find envelope width at Env_threshold * Intensity_max
# Limit_by_envelope = False # only fit up to n = int(z*env/3*Env_limit - 2/3) 
# Env_limit = 1.0 #used if Limit_by_envelope = True, rough measure to constrain n_curves fit according to data width & num fit params
# Min_Pops = 1 # Min_Pops to Max_Pops sets the range of populations to fit, set to same value to force single population
# Max_Pops = 3 #maximum number of underlying populations to fit
# Pop_Thresh = 0.03 #fall back to n-1 curves if population is below this, does not apply to bootstrap fits, but does exclude from boot average values
# Ncurve_p_accept = 0.05 #stringency for accepting more fit populations      
# Random_Seed = 16 #used for parameter initialization
# Scale_Y_Values = True # if Scale_Y_Values = True, plots will be in original Intensity units
#                 # fit will always be on normalized Intensity as it is much faster               
# Keep_Raw = True # peak_picking will retain the Raw spectrum if True, if False will only keep peaks, auto True for Test_Data
# Overlay_replicates = True #add column to figures that is overlay of all available replicates

# #######################################
#  '''end user input''';
# #######################################
