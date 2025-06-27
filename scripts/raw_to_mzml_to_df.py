# copyright 2025
# Lisa M Tuttle, last updated 27 June 2025


import os
import subprocess

import pandas as pd, numpy as np
from pyteomics import mzml, auxiliary #auxiliary.print_tree(testdata)

import matplotlib.pyplot as plt
#mono ../ThermoRawFileParser.exe -i=TOPdownHDX_UN_test_01.raw -o=output -f=1 -m=1

def raw_to_mzml(raw_files,out_dir,out_suffix="",peakpick=False,
                trp_path='/home/tuttle/data/HDX-MS/ThermoRaw/ThermoRawFileParser.exe',use_mono=True,user_flags=None):
    
    if use_mono == True: trfp = ['mono',trp_path]
    else: trfp = [trp_path]

    flags = '-f=1 -m=1' #mzml, create meta txt file, don't peak pick
    if user_flags: flags = flags +' '+ user_flags
    if peakpick == False: flags = flags+' -p'

    for raw_file in raw_files:
        out_file = os.path.join(out_dir,os.path.basename(raw_file)[:-4]+out_suffix+".mzml")
        meta_file = os.path.join(out_dir,os.path.basename(raw_file)[:-4]+out_suffix+"_metafile.txt")
        args = "-i="+raw_file+" -b="+out_file+" -c="+meta_file+" "+flags
        args = args.split()

        subprocess.run(trfp+args)

    return



def read_mzml(mzml_files,scan_idxs=[],exp_dict = {},use_map=False,user_cols=[],user_dict={},quiet=True):
    '''
    reads in data from specified scans of a mzML file
    specified scan_idxs values are treated as indices so e.g. -1 would give the last index
    exp_dict is a dictionary that can map the MSLevel to some different value, e.g. exp_type = {1:'Global',2:'ETD'}
    user_cols should specify what subset of the scan_info_dict should be included in the output dataframe
    user_dict can specify additional columns with the path to the desired info, see scan_info_dict for guidance
        e.g. precursor_z = scan.get('precursorList').get('precursor')[0].get('selectedIonList').get('selectedIon')[0].get('charge state')
        so is specified in the dictionary as 'precursur_z':['precursorList',['precursor'],'selectedIonList',['selectedIon'],'charge state']
    '''
    
    scan_info_dict = {'filter string': ['scanList',['scan'],'filter string'],
            'collision energy':['precursorList',['precursor'],'activation','collision energy'],
            'precursor_mz' : ['precursorList',['precursor'],'selectedIonList',['selectedIon'],'selected ion m/z'],
            'precursor_z' : ['precursorList',['precursor'],'selectedIonList',['selectedIon'],'charge state']}
    scan_info_dict = {**scan_info_dict,**user_dict}

    if user_cols.lower() == 'all':
        user_cols = list(scan_info_dict.keys())
    else: 
        user_cols = set(user_cols).intersection(set(scan_info_dict.keys()))
            
    spectra = pd.DataFrame()
    mzml_files = [mzml_files] if not isinstance(mzml_files, list) else mzml_files
    scan_idxs = [scan_idxs] if not isinstance(scan_idxs, list) else scan_idxs
    user_cols = [user_cols] if not isinstance(user_cols, list) else user_cols
    for mzml_file in mzml_files:
        sample = os.path.basename(mzml_file)[:-5]
        f = mzml.MzML(mzml_file)
        spec_ids = []
        if (len(scan_idxs)==0) | (scan_idxs == 'ALL'):            
            if use_map == True:
                for spec in f.map():
                    spec_ids += [spec['id']]
            else: 
                for di in f.default_index:
                    spec_ids += [di]
            #scan_idxs = [int(r.rsplit('=',1)[1]) for r in spec_ids]
        else:
            for scani in scan_idxs:
                spec_ids += [f.get_by_index(scani)['id']]
        print(spec_ids)
        for id in spec_ids:
            scandf = pd.DataFrame()
            try:
                scan = f.get_by_id(id)
                rt = scan['scanList']['scan'][0]['scan start time']
                if rt.unit_info == 'millisecond': rt = rt/1e3
                else: rt = rt * np.power(60.0,'smh'.find(rt.unit_info[0])) #seconds
                if len(exp_dict) > 0: exp_type = exp_dict[scan['ms level']]         
                else: exp_type = 'MS'+str(scan['ms level'])
                scandf = pd.DataFrame.from_dict(
                    {'sample':sample,
                    'id':scan['id'],
                    'idx':scan['index'],
                    'scan start time, sec': rt, 
                    'mz':scan['m/z array'],
                    'Intensity':scan['intensity array'],
                    'Exp Type':exp_type})                              
            except:
                if quiet == False: print("scan_id not found: "+str(id))
            
            for k,col_path in scan_info_dict.items():
                try:        
                    part_get = scan.get(col_path[0])
                    for col in col_path[1:]:
                        if isinstance(col,list):
                            part_get = part_get.get(col[0])[0]
                        else: 
                            part_get = part_get.get(col)
                    #print(k,": ",part_get)
                    scandf[k] = part_get
                except: pass
            spectra = pd.concat([spectra,scandf])
   
    return spectra
