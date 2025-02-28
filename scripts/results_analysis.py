'''
(hxex_utils.py will be more updated)
Utility functions for analysis after running pyHXExpress 
Includes some functionality that requires pyHDX 
      (my github fork has modifications to allow for replicate data)

'''

import os
import importlib
import sys

# hxex_path = os.path.join('')
# sys.path.append(hxex_path)

import hxex_updating as hxex
import numpy as np, pandas as pd
import config as config
from datetime import datetime

pd.set_option('display.max_columns',None)

now = datetime.now()
date = now.strftime("%d%b%Y")

def hxex_reload():
    importlib.reload(hxex)
    importlib.reload(config)
    hxex.config = config

hxex_reload()

###     PyHDX documentation, Examples
###     https://pyhdx.readthedocs.io/en/stable/examples.html
###
#!pip install pyhdx

import proplot as pplt
from scipy.optimize import lsq_linear
from pathlib import Path
import yaml
from Bio import SeqIO
from collections import defaultdict

from pyhdx.models import Coverage
from pyhdx.plot import peptide_coverage
from pyhdx.batch_processing import StateParser
from pyhdx.fitting import fit_d_uptake
from pyhdx.config import cfg
from pyhdx import read_dynamx, HDXMeasurement
from pyhdx.fitting import fit_rates_half_time_interpolate, fit_rates_weighted_average, fit_gibbs_global
from pyhdx.process import filter_peptides, apply_control, correct_d_uptake
from dask.distributed import Client
from pyhdx.plot import dG_scatter_figure
import pyhdx.plot

from collections import defaultdict

# import warnings
# warnings.simplefilter (action='ignore', )#category=FutureWarning)
# warnings.filterwarnings (action='ignore',)# category=RuntimeWarning,)

# pd.set_option('display.max_columns',None)



def combine_batch_data(dirs):
    # Combine the batch outputs into single dataframe files for metadf, datafits, and fitparams
    metadf_comb = pd.DataFrame()
    datafits_comb = pd.DataFrame()
    fitparams_comb = pd.DataFrame()

    metadf_prefix = "metadf_asrun_"
    datafits_prefix = "data_fits_asrun"
    fitparams_prefix = "fitparamsAll_asrun_"

    for dir in dirs:
        csvfiles = [ f for f in os.listdir(dir) if f[-4:]=='.csv'  ]
        mf = [ m for m in csvfiles if m[:len(metadf_prefix)]==metadf_prefix]
        df = [ d for d in csvfiles if d[:len(datafits_prefix)]==datafits_prefix]
        ff = [ f for f in csvfiles if f[:len(fitparams_prefix)]==fitparams_prefix]

        for f in mf:
            md = pd.read_csv(os.path.join(dir,f))#.drop('Index',axis=1)
            md['original_file'] = f
            metadf_comb = pd.concat([metadf_comb,md])#,ignore_index = True)
        for f in df:
            dd = pd.read_csv(os.path.join(dir,f))
            dd['original_file'] = f
            datafits_comb = pd.concat([datafits_comb,dd])
        for f in ff:
            fd = pd.read_csv(os.path.join(dir,f)).drop('Index',axis=1)
            fd['original_file'] = f
            fitparams_comb = pd.concat([fitparams_comb,fd],ignore_index= True)
        

    metadf_comb = metadf_comb.set_index('Index').sort_values('Index')
    datafits_comb = datafits_comb.sort_values('data_id').reset_index(drop=True)
    fitparams_comb = fitparams_comb.sort_values(['data_id','time_idx','charge','rep','ncurves','nboot']).reset_index(drop=True)

    return metadf_comb, datafits_comb, fitparams_comb



def convert_to_uptakedf(datafit,proj=None):
    # convert datafits into a format that hdexa_to_pyhdx can recognize (for residue level avg uptake calcs in pyHDX)
    # there is an assumption in this that the more protected pop1 goes together and less protected pop2 goes together
    # probably mostly okay, but may need to try to account for actual % populations in deciding what goes together 
    # this needs more consideration if ever more than 2 pops

    TD_time = 1e6
    test = datafit.copy()

    test[['start','end']] = test['peptide_range'].str.split('-',expand=True).astype('int')

    try: 
        test['fracDeut_1'] = test['Dabs_1']/test['max_namides']
    except: pass
    try: 
            test['fracDeut_2'] = test['Dabs_2']/test['max_namides']
    except: test['fracDeut_2'] = np.nan
    #cleanup garbage output - accidental carryover of polymodal data when rolling back to previous number of populations, should be fixed now
    test.loc[test['pop_1'] == 1.0, ['centroid_2','Dabs_2','Dabs_std_2','pop_2','pop_std_2','fracDeut_2']] = np.nan

    pops = test.fit_pops.max()
    testpop = {}
    testpop[1] = test.copy()
    testpop[1]['Protein State']=testpop[1].apply(lambda row: row['sample']+'_pop1',axis=1)
    testpop[1]['Theor Uptake #D'] = testpop[1]['Dabs_1']
    testpop[1]['Conf Interval (#D)'] = testpop[1]['Dabs_std_1']
    testpop[1]['pop'] = testpop[1]['pop_1']
    testpop[1]['pop_std'] = testpop[1]['pop_std_1']
    for pop in range(2,pops+1):
        testpop[pop] = testpop[1].copy() 
        testpop[pop]['Protein State']=testpop[pop].apply(lambda row: row['sample']+'_pop'+str(pop),axis=1)
        testpop[pop].loc[testpop[pop]['fit_pops']==pop,['Theor Uptake #D']] = testpop[pop]['Dabs_'+str(pop)]
        testpop[pop].loc[testpop[pop]['fit_pops']==pop,['Conf Interval (#D)']] = testpop[pop]['Dabs_std_'+str(pop)]
        testpop[pop].loc[testpop[pop]['fit_pops']==pop,['pop']] = testpop[pop]['pop_'+str(pop)]
        testpop[pop].loc[testpop[pop]['fit_pops']==pop,['pop_std']] = testpop[pop]['pop_std_'+str(pop)]

    uptakedf = pd.DataFrame()
    for pop in range(1,pops+1):
        uptakedf = pd.concat([uptakedf,testpop[pop]])

    uptakedf['#D'] = uptakedf['Theor Uptake #D'] * uptakedf['UN_TD_corr'] 
    uptakedf['%D'] = uptakedf['Theor Uptake #D'] / uptakedf['max_namides'] * 100.0
    uptakedf['Deut Time (sec)'] = uptakedf['time']
    uptakedf.loc[uptakedf['time']==TD_time,['Deut Time (sec)']] = 'MAX'
    uptakedf['Sequence'] = uptakedf['peptide']
    uptakedf['Peptide Mass'] = uptakedf['centroid']
    uptakedf['maxD'] = uptakedf['max_namides']
    uptakedf = uptakedf.drop(columns=['start','end'])
    uptakedf['Start'] = uptakedf['start_seq']
    uptakedf['End'] = uptakedf['end_seq']
    if not proj:
         proj = test['sample'].unique()[0]
    uptakedf['Protein'] = proj

    #filter some garbage data based on TD-UN
    uptakedf = uptakedf[uptakedf['UN_TD_corr']<1.05]
    
    return uptakedf



def hdexa_to_pyhdx(data,d_percentage=0.85,protein='protein'):
    ### Function to convert data table exported from HDExaminer to the processed DynamX format pyHDX expects
    ### this will leave any extra columns, but will chop out the MAX time points after processing
    
    drop_first=2

    def _time_to_sec(tp,tpunit):
        return  tp * np.power(60.0,'smh'.find(tpunit[0])) 
    if '# Deut' in data.columns:
        data = data.rename(columns={"# Deut":"#D"})
        data['#D'] = data['#D'].fillna(0.0)
        data['#D'] = data['#D'].astype(float)
    if 'Deut %' in data.columns:
        data = data.rename(columns={"Deut %":"%D"})
        data['%D'] = data['%D'].fillna(0.0)
        data['%D'] = data['%D'].astype(float)
    if 'Deut Time' in data.columns:
        data.loc[data['Deut Time'] == 'FD','Deut Time'] = '1e6s'
        data['time unit'] = data['Deut Time'].str[-1]
        data['Deut Time (sec)'] = data['Deut Time'].str[:-1].astype(float)
        data['Deut Time (sec)'] = data.apply(lambda x: _time_to_sec(tp=x['Deut Time (sec)'],tpunit=x['time unit']),axis=1)
        data.loc[data['Deut Time (sec)'] == 1e6,'Deut Time (sec)'] = 'MAX'
    if 'Protein' not in data.columns:
        data['Protein'] = protein


    pyhdx_cols = ['start', 'end' ,'stop' ,'sequence', 'state', 'exposure' ,'uptake' ,'maxuptake',
                'fd_uptake' ,'fd_uptake_sd' ,'nd_uptake' ,'nd_uptake_sd' ,'rfu', 'protein',
                'modification', 'fragment', 'mhp' ,'center' ,'center_sd' ,'uptake_sd' ,'rt',
                'rt_sd' ,'rfu_sd' ,'_sequence' ,'_start' ,'_stop' ,'ex_residues',
                'uptake_corrected']
    data = data.rename(columns={
                "Protein State":"state",
                "Protein":"protein",
                "Start":"start",
                "End":"end",
                "Sequence":"_sequence",
                "Peptide Mass":"mhp",
                "RT (min)":"rt",
                "Deut Time (sec)":"exposure",
                "maxD":"maxuptake",
                "Theor Uptake #D":"uptake_corrected",
                "#D":"uptake",
                "%D":"rfu",
                "Conf Interval (#D)":"rfu_sd",
                "#Rep":"rep",
                "Confidence":"quality",
                "Stddev":"center_sd",
                #"p"
    })

    missing = list(set(pyhdx_cols)-set(data.columns))
    for mcol in missing:
        data[mcol] = np.nan
        if mcol == "rfu_sd": data[mcol] = 0.05 #set 5% error as dummy value

    data['rfu']=data['rfu']/100.
    data.loc[data['exposure']=="0",'rfu_sd']=0.0
    data['stop']=data['end']+1
    data['sequence']=data["_sequence"].copy()
    data['sequence']=[s.replace("P", "p") for s in data["sequence"]]
    # Find the total number of n terminal / c_terminal residues to remove from pyhdx/process.py
    n_term = np.array([len(seq) - len(seq[drop_first:].lstrip("p")) for seq in data["sequence"]])
    c_term = np.array([len(seq) - len(seq.rstrip("p")) for seq in data["sequence"]])
    data["sequence"] = ["x" * nt + s[nt:] for nt, s in zip(n_term, data["sequence"])]
    data["_start"] = data["start"] + n_term
    data["_stop"] = data["stop"] - c_term
    ex_residues = (np.array([len(s) - s.count("x") - s.count("p") for s in data["sequence"]])* d_percentage)
    data["ex_residues"] = ex_residues
    data["uptake_sd"]=data["center_sd"]
    data["nd_uptake"]=0.0
    data["nd_uptake_sd"]=0.0
    data["modification"]=float("nan")
    data["fragment"]=float("nan")
    # upeps = data[data["exposure"]=="0"]["_sequence"].unique()
    # fpeps = data[data["exposure"]=="MAX"]["_sequence"].unique()
    # good_peps = np.array(list(set(upeps) & set(fpeps)))
    #peps = data["_sequence"].unique()
    states = data["state"].unique()
    data["fd_uptake"]="novalue"
    data["fd_uptake_sd"]="novalue"

    for state in states:
        peps = data[data["state"]==state]["_sequence"].unique()
        for pep in peps:
            try: #may have gotten here without a TD measurement 
                fd_up = data[(data["_sequence"]==pep) & (data["exposure"]=="MAX")& (data["state"]==state)]['uptake'].iat[0]
                fd_up_sd = data[(data["_sequence"]==pep) & (data["exposure"]=="MAX")& (data["state"]==state)]['center_sd'].iat[0]
            except: 
                fd_up = data[(data["_sequence"]==pep) & (data["state"]==state)]['maxuptake'].iat[0]
                fd_up_sd = data[(data["_sequence"]==pep) & (data["state"]==state)]['maxuptake'].iat[0]
            data.loc[data["_sequence"]==pep, "fd_uptake"]=fd_up
            data.loc[data["_sequence"]==pep, "fd_uptake_sd"]=fd_up_sd
    data["center"]=data["mhp"]+data["uptake"]
    data["rt_sd"]=0.05 #dummy value

    data['uptake_corrected_orig'] = data['uptake_corrected']
    data['uptake_corrected'] = data["rfu"]*data['maxuptake']

    
    data = data[data["exposure"] != "MAX"]
    data = data[data["fd_uptake"] != 0]
    data = data[~data["uptake"].isna()]
    data["exposure"]=data["exposure"].astype(float)

    new_columns = [col for col in pyhdx_cols if col in data.columns] + [col for col in data.columns if col not in pyhdx_cols]
    return data[new_columns]


def prepare_kwargs(fit_result):
    """Prepare plot kwargs for fit result"""
    d = {
        "y": fit_result.d_uptake.mean(axis=0),
        "fadedata": np.percentile(fit_result.d_uptake, (5, 95), axis=0),
        "shadedata": np.percentile(fit_result.d_uptake, (25, 75), axis=0),
    }
    return d


## Tweaked from pyhdx.plot.single_linear_bar to leave levels in kwargs 
def plot_bar(ax, x, z, cmap, norm, height=1,**kwargs):
    """makes a linear bar plot on supplied axis with values z and corresponding x values x"""

    if isinstance(z, pd.Series):
        z = z.to_numpy()
    elif isinstance(z, pd.DataFrame):
        assert len(z.columns) == 1, "Can only plot dataframes with 1 column"
        z = z.to_numpy().squeeze()

    img = np.expand_dims(z, 0)
    
    collection = ax.pcolormesh(
        pplt.edges(x),
        np.array([0, height]),
        img,
        cmap=cmap,
        vmin=norm.vmin,
        vmax=norm.vmax,
        **kwargs,
    )
    ax.format(yticks=[])

    return collection

from pyhdx.plot import peptide_coverage

def plot_coverage(hdxm,states=None,times=None,savepath=None,peprange=None):
    ## Coverage plots
    use_hdxm = hdxm.copy()
    if states is None:
        states = list(use_hdxm.keys())
    if times is None:
        times = sorted(use_hdxm[states[0]].timepoints)

    fig, axes = pplt.subplots(nrows=len(states),ncols=len(times),  axwidth="200mm", sharey=False, refaspect=2,)

    for j,mutant in enumerate(states):
        for i, use_time in enumerate(times):  
            time_idx = np.where(hdxm[mutant].timepoints == use_time)[0][0]
            filtdf = use_hdxm[mutant][time_idx].data
            if peprange: filtdf = filter_range(filtdf,peprange)
            peptide_coverage(axes[j,i], filtdf,  cbar=True,linewidth=0.1)
            t = axes[j,i].set_title(f'{mutant} Peptides t = {int(use_time)}s')
            l = axes[j,i].set_xlabel('Residue number')
    
    if savepath:
        fig.savefig(savepath,format='pdf',dpi=600)
        
    return


def plot_rfu_residue(hdxm,states=None,times=None,colors=None,savepath=None,legendcols=5):
    # Prolines: 7,12,15,19,38,45,50,51,57,85,124,129,147,154,159,166,172,
    if colors is None:
        colors="#000000 #66C2A5 #56B4E9 #7570B3 #E7298A".split()
    if states is None:
        states = list(hdxm.keys())
    if times is None:
        times = sorted(hdxm[states[0]].timepoints)
    
    z_value ={'50':0.674,'68':1.0,'sd':1.0,'80':1.282,'90':1.645,'95':1.96,'98':2.326,'99':2.576} #confidence intervals multiplier 
  
    nset = 3 #rfu, redundancy, resolution 
    #ncols = 2
    nrows = len(times)*nset
    nfigs = nrows
    array=[]
    for i in range(1,nfigs+1):
        array += [[i,i]]
    hratios = [10,1,1]*len(times)

    fig, axes = pplt.subplots(array,  axheight="30mm",axwidth="180mm",  #refaspect=10,
                                wspace=(0), hspace=([0.7,0.3,3.0]*(len(times)-1)+[0.7,0.3]),
                                sharex=False, sharey=False, hratios=hratios)
    
    vmin = 0
    vmax = 20
    kwargs = dict(levels=pplt.arange(0, vmax*3, 1),) #extendsize=2.5, extendrect=True)
    res_kwargs = dict(levels=pplt.arange(0, vmax, 1),)


    for i, use_time in enumerate(times):

        for j,mutant in enumerate(states):
            time_idx = np.where(hdxm[mutant].timepoints == use_time)[0][0]
            norm = 1#hdxm_hxex[mutant][-1].rfu_residues.values[-1]
            x_res = hdxm[mutant][time_idx].r_number
            center = hdxm[mutant][time_idx].rfu_residues* 1/norm
            #high =  center + z_value['95']*hdxm[mutant][time_idx].rfu_residues_sd
            high50 = center + z_value['50']*hdxm[mutant][time_idx].rfu_residues_sd
            #low = center - z_value['95']*hdxm[mutant][time_idx].rfu_residues_sd
            low50 = center - z_value['50']*hdxm[mutant][time_idx].rfu_residues_sd

            #axes[i*nset].line(hdxm[mutant][time_idx].r_number,hdxm[mutant][time_idx].rfu_residues * 1/norm, label=str(mutant), color=colors[j])
            axes[i*nset].line(x_res,center, shadedata=(low50,high50),label=str(mutant), color=colors[j%len(colors)]) #,fadedata=(low,high)
            axes[i*nset].format(xlabel="",xtickrange=(-1,-1))
            axes[i*nset].format(title=str(use_time)+'s',titleloc= 'lower right',)
        
        redundancy = hdxm[states[0]].coverage.X.sum(axis=0).astype(float)
        resolution = np.repeat(hdxm[states[0]].coverage.block_length, hdxm[states[0]].coverage.block_length)
        resolution = resolution.astype(float)
        resolution[redundancy == 0] = np.nan
        redundancy[redundancy == 0] = np.nan
        red = plot_bar(axes[i*nset+1],hdxm[states[0]].coverage.r_number,redundancy,
                'blues',pplt.Norm("linear", vmin=vmin, vmax=vmax*3),height=1.0,**kwargs)
        axes[i*nset+1].format(xtickrange=(-2,-1))
        res = plot_bar(axes[i*nset+2],hdxm[states[0]].coverage.r_number,resolution, 
                'fire',pplt.Norm("linear", vmin=vmin, vmax=vmax),height=1.0,**res_kwargs)

        axes[i*nset].format(ylabel="D-uptake")
        axes[i*nset+1].set_ylabel("red.",rotation=0,)#loc='center',labelpad=10)
        axes[i*nset+1].yaxis.set_label_coords(-.02,-.05)
        axes[i*nset+1].format(xtickloc='none')
        axes[i*nset+2].set_ylabel("res.",rotation=0,)#loc='center',labelpad=10)
        axes[i*nset+2].yaxis.set_label_coords(-.02,-.05)

    #need to implement https://proplot.readthedocs.io/en/latest/subplots.html?highlight=spacing#Spacing-and-tight-layout
    
    fig.colorbar(red, label="Redundancy (peptides including replicates)",width=0.1,loc='b',length=0.5,col=1,ticks=np.arange(0, vmax*3+1, 5))# **kwargs)
    fig.colorbar(res, label="Resolution (residues)",width=0.1,loc='b',length=0.5,col=2)# **kwargs)
    
    axes[0].legend(loc="t", ncols=legendcols)
    axes.format(ylim=(0,1.1),)#xgrid=True,ygrid=True,ytickrange=(0,1.1))
    axes.format(xlim=(0,180))

    axes[i*nset+2].format(xlabel="Residue number")

    # save_path = os.path.join(project_dir,'B5_phospho_2popFiltShade_RFU_residueavg_plots_'+date+'.pdf')
    if savepath:
        fig.savefig(savepath,format='pdf',dpi=600)

    return #fig

def filter_range(hdxm,peprange):
    df = hdxm.copy()
    peprange = [peprange] if not isinstance(peprange,list) else peprange
    df = df[(df['start'] <= peprange[-1]) & (peprange[0] <= df['end'])]
    return df
