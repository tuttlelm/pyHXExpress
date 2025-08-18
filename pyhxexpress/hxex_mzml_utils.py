
import os
import sys
import subprocess
import pandas as pd, numpy as np
from pyteomics import mzml, auxiliary #auxiliary.print_tree(testdata)
import matplotlib.pyplot as plt

from datetime import datetime
now = datetime.now()
date = now.strftime("%d%b%Y")


hxex_path = os.path.join('/home/tuttle/data/HDX-MS/pyHXExpress')
sys.path.append(hxex_path)

#import hxex_kinetics_config as config
import hxex_updating as hxex
from hxex_utils import autoscale
#import hxex_fragments
#from hxex_fragments import get_fragments



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



def read_mzml(mzml_files,scan_idxs=[],use_map=False,user_cols=[],user_dict={},keep_mz_and_Int = False,quiet=True,keep_meta = pd.DataFrame()):
    ''' 21july2025
    reads in data from specified scans of a mzML file
    specified scan_idxs values are treated as indices so e.g. -1 would give the last index
    exp_dict is a dictionary that can map the MSLevel to some different value, e.g. exp_type = {1:'Global',2:'ETD'}
    user_cols should specify what subset of the scan_info_dict should be included in the output dataframe
    user_dict can specify additional columns with the path to the desired info, see scan_info_dict for guidance
        e.g. precursor_z = scan.get('precursorList').get('precursor')[0].get('selectedIonList').get('selectedIon')[0].get('charge state')
        so is specified in the dictionary as 'precursur_z':['precursorList',['precursor'],'selectedIonList',['selectedIon'],'charge state']
    m/z vs Intensity spectra will only be stored if keep_mz_and_Int=True (default is False)
    keep meta will keep the columns of a meta data type dataframe (which can contain additional info about each dataset), matched on 'file' 
    '''
    
    scan_info_dict = {'filter string': ['scanList',['scan'],'filter string'],
            'collision energy':['precursorList',['precursor'],'activation','collision energy'],
            'precursor_mz' : ['precursorList',['precursor'],'selectedIonList',['selectedIon'],'selected ion m/z'],
            'precursor_z' : ['precursorList',['precursor'],'selectedIonList',['selectedIon'],'charge state']}
    scan_info_dict = {**scan_info_dict,**user_dict}

    if 'all' in [str(x).lower() for x in user_cols]:
        user_cols = list(scan_info_dict.keys())
    else: 
        user_cols = set(user_cols).intersection(set(scan_info_dict.keys()))
            
    spectra = pd.DataFrame()
    mzml_files = [mzml_files] if not isinstance(mzml_files, list) else mzml_files
    scan_idxs = [scan_idxs] if not isinstance(scan_idxs, list) else scan_idxs
    user_cols = [user_cols] if not isinstance(user_cols, list) else user_cols
    for mzml_file in mzml_files:
        sample = os.path.basename(mzml_file)#[:-5] #keep extension
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
        #print(spec_ids)
        for id in spec_ids:
            scandf = pd.DataFrame()
            try:
                scan = f.get_by_id(id)
                rt = scan['scanList']['scan'][0]['scan start time']
                if rt.unit_info == 'millisecond': rt = rt/1e3
                else: rt = rt * np.power(60.0,'smh'.find(rt.unit_info[0])) #seconds
                ms_level = 'MS'+str(scan['ms level'])
                if keep_mz_and_Int == True:
                    scandf = pd.DataFrame.from_dict(
                        {'file':sample,
                        'id':scan['id'],
                        'idx':scan['index'],
                        'scan start time, sec': rt, 
                        'mz':scan['m/z array'],
                        'Intensity':scan['intensity array'],
                        'MS Level':ms_level})    
                else:
                    scandf = pd.DataFrame.from_dict(
                        {'file':[sample],
                        'id':[scan['id']],
                        'idx':[scan['index']],
                        'scan start time, sec': [rt], 
                        'MS Level':[ms_level]}) 
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
            spectra = pd.concat([spectra,scandf]).reset_index(drop=True)
    if len(keep_meta.index) > 0:
        spectra = pd.merge(spectra,keep_meta,on="file",how="left").set_index(spectra.index)
    return spectra



def categorize_exps(row):
    if row['MS Level'] == 'MS1':
        return 'MS1'
    fstr = row['filter string']
    exptype = ''
    if 'uvpd' in fstr:
        exptype = 'UVPD'
    elif 'hcd' in fstr:
        exptype = 'EThcD'
    elif 'etd' in fstr:
        exptype = 'ETD'
    if 'sid' in fstr:
        exptype = 'High'+exptype
    else: exptype = 'Low'+exptype
    if 'precursor_z' in row.keys():
        pre_z = '_z'+str(int(row['precursor_z']))
    else: pre_z = ''
    return exptype+pre_z


def spectra_peakpicker(metadf,spectradf,use_idx=[],index_column = 'idx',resolution=hxex.config.Peak_Resolution,
                       time_col='scan start time, sec',TrimRaw=True,ploteach=False):
    # also uses config.Zero_Filling within peak_picker, # should put it as variable in peak_picker() TODO

    peakpicked = pd.DataFrame()
    rawFull = pd.DataFrame()
    if len(use_idx) == 0: use_idx = spectradf[index_column].unique()
    for index,row in metadf.iterrows():
        seq = row['peptide'] 
        charge = int(row['charge']) 
        mod_dic = {x.split(':')[0]:(x.split(':')[1]) for x in row['modification'].split()}
        if 'ion_type' in mod_dic.keys():
            ion_type = mod_dic['ion_type']
        else: ion_type = 'M'
        undeut_mz = hxex.mass.calculate_mass(seq,charge=charge,ion_type=ion_type)
        n_amides = hxex.count_amides(seq)
        na_buffer = len(hxex.get_na_isotope(seq,charge,mod_dict=mod_dic))//2 
        delta_m = [1.003355, 1.003355, 1.004816, 1.004816, 1.006277] #dmC, dmC, dmHCavg, dmHCavg, dmH ...
        dmz = [0.0]
        for nd in range(0,n_amides+hxex.config.Zero_Filling+na_buffer):
            dm = delta_m[nd] if nd < 5 else delta_m[-1]
            dmz += [dmz[nd] + dm/charge]
        pred_mzs = np.array(dmz) + undeut_mz
        pred_width = max(pred_mzs) - min(pred_mzs)
        
        if ploteach:
            color = iter(plt.cm.rainbow(np.linspace(1,0,len(use_idx)+1)))
            c = next(color).reshape(1,-1)
            fig,ax = plt.subplots(figsize=(24,5),ncols=1,num=1,clear=True)
        for idx in use_idx:
            raw = spectradf.copy()[(spectradf[index_column]==idx) ]     
            #print("raw time A",raw[time_col].unique()) ## TROUBLESHOOTING
            rt = raw[time_col].values[0]  
            if TrimRaw == True: raw = raw[raw['mz'].between(min(pred_mzs)-pred_width/2,max(pred_mzs)+pred_width/2)]
            testpick = hxex.peak_picker(raw,seq,charge,resolution=resolution,mod_dict=mod_dic)
            testpick['time'] = rt
            raw['time'] = rt
            if ploteach:
                timelabel = "{:.1f}".format(rt/60.0)+" min"
                c = next(color).reshape(1,-1)
                raw.plot('mz','Intensity',label=timelabel,legend=True,ax=ax,color=c,alpha=0.5)
                ax.scatter(testpick['mz'].values,testpick['Intensity'].values,color=c)
            for k in row.keys():
                if k not in ['mz','Intensity','time',time_col,index_column]: #neeed to consider 'time'
                    testpick[k] = row[k]
                    raw[k] = row[k]
            testpick['data_id'] = index
            raw['data_id'] = index
            testpick['time_idx'] = idx #need to fix this for UnDeut cases 
            raw['time_idx'] = idx
            #print("raw time B",raw[time_col].unique()) ## TROUBLESHOOTING
            peakpicked = pd.concat([peakpicked,testpick],ignore_index = True).reset_index(drop=True)
            rawFull = pd.concat([rawFull,raw],ignore_index=True).reset_index(drop=True)
        if ploteach:
            max_int = rawFull[rawFull['mz'].between(min(pred_mzs)-pred_width/2,max(pred_mzs)+pred_width/2)].Intensity.max()
            ax.scatter(pred_mzs,np.zeros(len(pred_mzs))-0.1*max_int,marker=6,color='k')
            fig.suptitle(row['peptide_range']+" "+seq+" ion_type:"+ion_type+" z"+str(charge))
            ax.set_xlim(min(pred_mzs)-2,max(pred_mzs)+2)
            ax.legend(loc='upper left',bbox_to_anchor=[1.0,1.0])
            try: autoscale(ax,'y')
            except: pass
            plt.show()
    return peakpicked, rawFull