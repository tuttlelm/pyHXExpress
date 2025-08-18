import os
from collections import defaultdict
import numpy as np
import pandas as pd

# Combine the batch outputs into single dataframe files for metadf, datafits, and fitparams
def combine_batch_data(dirs):
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

    #filter some garbage data based on TD-UN ... leave this to baddata_id when reading in data
    #uptakedf = uptakedf[uptakedf['UN_TD_corr']<1.05]
    
    return uptakedf



def hdexa_to_pyhdx(data,d_percentage=0.85,protein='protein',dummyTD = '1e6s', keepTD = False):
    ### Function to convert data table exported from HDExaminer to the processed DynamX format pyHDX expects
    ### this will leave any extra columns
    ### If keepTD = false, this will chop out the MAX time points after processing
    ### if keepTD = True, TD timepoins will be kept with exposure = dummyTDs
    
    drop_first=2

    def _time_to_sec(tp,tpunit):
        return  tp * np.power(60.0,'smh'.find(tpunit[0])) 
    dummyTDs = _time_to_sec(float(dummyTD[:-1]),dummyTD[-1])
    if '# Deut' in data.columns:
        data = data.rename(columns={"# Deut":"#D"})
        data['#D'] = data['#D'].fillna(0.0)
        data['#D'] = data['#D'].astype(float)
    if 'Deut %' in data.columns:
        data = data.rename(columns={"Deut %":"%D"})
        data['%D'] = data['%D'].fillna(0.0)
        data['%D'] = data['%D'].astype(float)
    if 'Deut Time' in data.columns:
        data.loc[data['Deut Time'] == 'FD','Deut Time'] = dummyTD
        data['time unit'] = data['Deut Time'].str[-1]
        data['Deut Time (sec)'] = data['Deut Time'].str[:-1].astype(float)
        data['Deut Time (sec)'] = data.apply(lambda x: _time_to_sec(tp=x['Deut Time (sec)'],tpunit=x['time unit']),axis=1)
        data.loc[data['Deut Time (sec)'] == dummyTDs,'Deut Time (sec)'] = 'MAX'
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
                fd_up = data[(data["_sequence"]==pep) & (data["exposure"]=="MAX")& (data["state"]==state)]['uptake'].mean()#.iat[0]
                fd_up_sd = data[(data["_sequence"]==pep) & (data["exposure"]=="MAX")& (data["state"]==state)]['uptake'].std()#.iat[0]
                nd_up_sd = data[(data["_sequence"]==pep) & (data["exposure"]==0)& (data["state"]==state)]['uptake'].std()
            except: 
                fd_up = data[(data["_sequence"]==pep) & (data["state"]==state)]['maxuptake'].iat[0]
                fd_up_sd = 0.0
                nd_up_sd = 0.0
            data.loc[data["_sequence"]==pep, "fd_uptake"]=fd_up
            data.loc[data["_sequence"]==pep, "fd_uptake_sd"]=fd_up_sd
            data.loc[data["_sequence"]==pep, "nd_uptake_sd"]=nd_up_sd
    data["center"]=data["mhp"]+data["uptake"]
    data["rt_sd"]=0.05 #dummy value

    data['uptake_corrected_orig'] = data['uptake_corrected']
    data['uptake_corrected'] = data["rfu"]*data['maxuptake']

    
    if keepTD:
        data.loc[data["exposure"] == "MAX","exposure"] = dummyTDs
    else:
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
    import proplot as pplt

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



def plot_coverage(hdxm,states=None,times=None,savepath=None,svg=False,peprange=None,color_field='rfu'):
    ## Coverage plots
    from pyhdx.plot import peptide_coverage
    import proplot as pplt

    use_hdxm = hdxm.copy()
    if states is None:
        states = list(use_hdxm.keys())
    if times is None:
        times = sorted(use_hdxm[states[0]].timepoints)

    fig, axes = pplt.subplots(nrows=len(states),ncols=len(times),  axwidth="180mm", sharey=False, refaspect=2,)

    for j,mutant in enumerate(states):
        for i, use_time in enumerate(times):  
            time_idx = np.where(hdxm[mutant].timepoints == use_time)[0][0]
            filtdf = use_hdxm[mutant][time_idx].data
            if peprange: filtdf = filter_range(filtdf,peprange)
            peptide_coverage(axes[j,i], filtdf,  color_field=color_field, cbar=True,linewidth=0.1)
            t = axes[j,i].set_title(f'{mutant} Peptides t = {int(use_time)}s')
            l = axes[j,i].set_xlabel('Residue')
    
    if savepath:
        if svg:
            fig.savefig(savepath,format='svg')
        else:
            fig.savefig(savepath,format='pdf',dpi=600)
        
    return


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib as mpl

    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_rfu_residue(hdxm,states=None,times=None,seq=None,colors=None,savepath=None,svg=False,legendcols=5,plotZero=True,TD_time=1e6,UN_time=0):
    # Prolines: 7,12,15,19,38,45,50,51,57,85,124,129,147,154,159,166,172,

    import proplot as pplt
    import matplotlib as mpl
    import colorcet as cc

    if colors is None:
        colors="#000000 #66C2A5 #56B4E9 #7570B3 #E7298A".split()
    if states is None:
        states = list(hdxm.keys())
    if times is None:
        times = sorted(hdxm[states[0]].timepoints)
    if not plotZero:
        times.remove(0)
    
    z_value ={'50':0.674,'68':1.0,'sd':1.0,'80':1.282,'90':1.645,'95':1.96,'98':2.326,'99':2.576} #confidence intervals multiplier 
  
    nset = 3 #rfu, redundancy, resolution 
    ncols = 20
    nrows = len(times)*nset
    nfigs = nrows
    array=[]
    for i in range(1,nfigs+1):
        array += [np.repeat(i,ncols)]# [[i,i]]
    hratios = [10,1,1]*len(times)

    fig, axes = pplt.subplots(array,  axheight="30mm",axwidth="180mm",  #refaspect=10,
                                wspace=(0), hspace=([0.7,0.3,3.0]*(len(times)-1)+[0.7,0.3]),
                                sharex=False, sharey=False, hratios=hratios)
    
    vmin = 0
    vmax = 20
    red_scale = 2
    kwargs = dict(levels=pplt.arange(0, vmax, 1),) #extendsize=2.5, extendrect=True)
    res_kwargs = dict(levels=pplt.arange(0, vmax, 1),)


    for i, use_time in enumerate(times):
        if use_time == TD_time:
            time_label = 'FullDeut'
        elif use_time == UN_time:
            time_label = 'UnDeut'
        else:
            time_label = str(float(use_time))+'s'

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
            axes[i*nset].format(title=time_label,titleloc= 'lower right',)
        
        redundancy = hdxm[states[0]].coverage.X.sum(axis=0).astype(float)
        resolution = np.repeat(hdxm[states[0]].coverage.block_length, hdxm[states[0]].coverage.block_length)
        resolution = resolution.astype(float)
        resolution[redundancy == 0] = np.nan
        redundancy[redundancy == 0] = np.nan

        if seq is not None:
            if isinstance(seq,dict):
                seq = seq[states[0]]
            pro_index = [index+1 for index,s in enumerate(seq) if (s=="P")]
            pro_index = sorted(set(pro_index) & set(hdxm[states[0]].coverage.r_number))
            redundancy = pd.Series(redundancy,index = hdxm[states[0]].coverage.r_number)
            resolution = pd.Series(resolution,index = hdxm[states[0]].coverage.r_number)
            redundancy[pro_index] = np.nan
            resolution[pro_index] = np.nan

        #print("resolution",resolution)
        res_cmap = cc.cm['CET_L4_r'].copy() #.set_under() .set_over()
        res_cmap = truncate_colormap(res_cmap,0,0.7)
        # res_cmap = cc.cm['coolwarm'].copy()
        # res_cmap = truncate_colormap(res_cmap,0,1.0)
        res_cmap.set_bad(color='k')
        red_cmap = cc.cm['CET_D1A_r'].copy()
        red_cmap = truncate_colormap(red_cmap,0.1,0.9)
        red_cmap.set_bad(color='k')

        red = plot_bar(axes[i*nset+1],hdxm[states[0]].coverage.r_number,redundancy,
                red_cmap,pplt.Norm("linear", vmin=vmin, vmax=vmax*red_scale),height=1.0,**kwargs)
        axes[i*nset+1].format(xtickrange=(-2,-1))
        res = plot_bar(axes[i*nset+2],hdxm[states[0]].coverage.r_number,resolution, 
                res_cmap,pplt.Norm("linear", vmin=vmin, vmax=vmax),height=1.0,**res_kwargs)               

        axes[i*nset].format(ylabel="D-uptake")
        axes[i*nset+1].set_ylabel("red.",rotation=0,)#loc='center',labelpad=10)
        axes[i*nset+1].yaxis.set_label_coords(-.02,-.05)
        axes[i*nset+1].format(xtickloc='none')
        axes[i*nset+2].set_ylabel("res.",rotation=0,)#loc='center',labelpad=10)
        axes[i*nset+2].yaxis.set_label_coords(-.02,-.05)

    #need to implement https://proplot.readthedocs.io/en/latest/subplots.html?highlight=spacing#Spacing-and-tight-layout
    
    sm = mpl.cm.ScalarMappable(cmap= mpl.colors.ListedColormap([res_cmap.get_bad()]))

    cbar_width = 0.4
    gap = 1
    cols = [[1,int(ncols*cbar_width)], 
                [int(ncols*cbar_width)+1+gap,int(ncols*cbar_width*2)+1], 
                [int(ncols*cbar_width*2)+2+gap,ncols]]

    fig.colorbar(red, label="Redundancy (peptides including replicates)",width=0.1,pad=1.5,
                 loc='b',length=1.0,col=cols[0],ticks=np.arange(0, vmax+1, 5))# **kwargs)
    fig.colorbar(res, label="Resolution (residues)",width=0.1,
                 loc='b',length=1.0,col=cols[1],ticks=np.arange(0, vmax+1, 5))# **kwargs)
    pcb = fig.colorbar(sm,width=0.1,loc='b',col=cols[2],length=0.1,ticks=[],)
    pcb.set_label(label="No Coverage\nor Proline",labelpad=8)
    
    axes[0].legend(loc="t", ncols=legendcols)
    axes.format(ylim=(0,1.1),)#xgrid=True,ygrid=True,ytickrange=(0,1.1))
    
    #axes.format(xlim=(0,180))

    axes[i*nset+2].format(xlabel="Residue")
    fig.set_facecolor('white')
    # save_path = os.path.join(project_dir,'B5_phospho_2popFiltShade_RFU_residueavg_plots_'+date+'.pdf')
    if savepath:
        if svg:
            fig.savefig(savepath,format='svg',)
        else:
            fig.savefig(savepath,format='pdf',dpi=600)

    return #fig


def choose_fits(hdxm,state,times=None,avg_proj = 'fixed1pop',fixed2_proj = 'fixed2pops',z='50',
                    min_pop=0.05,max_sd=0.25,plot=True,plotpop=True,return_fits=False,TD_time=1e6,UN_time=0,savepath=None,svg=False):
    '''
    Based on the fixed1pop and fixed2pop fits and their errors, decide what the approprate number of fit populations is
    use bimodal if justified or fallback to the 'avg' value from teh fixed1pop (not a true average but what we can measure)
    Outputs are the plot and the fit values if return_fits = True

    hdxm:    the HDXMeasurement object from pyHDX (generated from other utility functions in hxed_utils)
    times:   the time value for the comparison
    z:       the z_value confidence % to use (e.g. value +/- z_value[z]*sd is solution with z% confidence)
             this is a little hand wavey since errors aren't likely guassian
    max_sd:  falls back to 'avg' fit if either pop1 or pop2 errors are larger than this cutoff 
    '''

    import proplot as pplt

    z_value ={'50':0.674,'68':1.0,'sd':1.0,'80':1.282,'90':1.645,'95':1.96,'98':2.326,'99':2.576} #confidence intervals multiplier 
    all_fits = defaultdict(dict)

    if times is None:
        times = sorted(hdxm[state].timepoints)
    times = [times] if not isinstance(times,list) else times

    if plotpop:
        nset = 2 #rfu, population
        hratios = [2,1]*(len(times))
        hspace = ([1.0,5.0]*(len(times)-1)+[1.0])
    else: 
        nset = 1
        hratios = [10]*len(times)
        hspace = ([0.7]*(len(times)-1))

    #print("hspace",hspace,len(hspace))
    #print("hratios",hratios,len(hratios))

    ncols = 20
    nrows = len(times)*nset
    nfigs = nrows
    array=[]
    for i in range(1,nfigs+1):
        array += [np.repeat(i,ncols)]# [[i,i]]
    
    if plot:
        fig, axes = pplt.subplots(array,  axheight="30mm",axwidth="180mm",  #refaspect=10,
                                wspace=(0), hspace=hspace, sharex=False, sharey=False, hratios=hratios)


    for i, use_time in enumerate(times):

        if use_time == TD_time:
            time_label = 'FullDeut'
        elif use_time == UN_time:
            time_label = 'UnDeut'
        else:
            time_label = str(float(use_time))+'s'        

        rfu = defaultdict(dict)

        rfu['center']['p1'] = hdxm[fixed2_proj][state+'_pop1'].rfu_residues[use_time]
        rfu['center']['p2'] = hdxm[fixed2_proj][state+'_pop2'].rfu_residues[use_time]
        rfu['center']['avg'] = hdxm[avg_proj][state+'_pop1'].rfu_residues[use_time]

        rfu['sd']['p1'] = hdxm[fixed2_proj][state+'_pop1'].rfu_residues_sd[use_time]
        rfu['sd']['p2'] = hdxm[fixed2_proj][state+'_pop2'].rfu_residues_sd[use_time]
        rfu['sd']['avg'] = hdxm[avg_proj][state+'_pop1'].rfu_residues_sd[use_time]

        resi = rfu['center']['p1'].index.union(rfu['center']['p2'].index).union(rfu['center']['avg'].index)
        for k in rfu['center'].keys():
            rfu['center'][k] = rfu['center'][k].reindex(resi)
            rfu['sd'][k] = rfu['sd'][k].reindex(resi)

        rfu['low']['p1'] = rfu['center']['p1'] - z_value[z]*rfu['sd']['p1']
        rfu['low']['p2'] = rfu['center']['p2'] - z_value[z]*rfu['sd']['p2']
        rfu['low']['avg'] = rfu['center']['avg'] - z_value[z]*rfu['sd']['avg']

        rfu['high']['p1'] = rfu['center']['p1'] + z_value[z]*rfu['sd']['p1']
        rfu['high']['p2'] = rfu['center']['p2'] + z_value[z]*rfu['sd']['p2']
        rfu['high']['avg'] = rfu['center']['avg'] + z_value[z]*rfu['sd']['avg']

        sep = (rfu['center']['p2'] - rfu['center']['p1']) - (rfu['sd']['p2'] + rfu['sd']['p1'])*z_value[z]

        pop1 = (rfu['center']['p2'] - rfu['center']['avg'])/(rfu['center']['p2'] - rfu['center']['p1'])
        #pop1 = np.where((sep < 0) | (rfu['sd']['p1'] > max_sd) |  (rfu['sd']['p2'] > max_sd), 1.0, pop1)
        pop1 = np.where(rfu['center']['p1'].isna(),np.nan,pop1)
        pop1 = np.where(rfu['center']['p2'].isna(),np.nan,pop1)
        pop1 = np.where((rfu['center']['p1'] < rfu['center']['avg']) & (rfu['center']['avg'] < rfu['center']['p2']),pop1,np.nan) #check avg between fits
        pop1 = pop1.clip(0,1)
        #pop1 = np.where(np.isnan(rfu_p1) & np.isnan(rfu_p2),np.nan,pop1) #if missing pop1/2 fits use avg

        rfu_p1 = np.where((rfu['center']['p1'] > 1.1) | (rfu['center']['p1'] < -0.1),np.nan,rfu['center']['p1']) #as it stands these still contribute to stdev
        rfu_p2 = np.where((rfu['center']['p2'] > 1.1) | (rfu['center']['p2'] < -0.1),np.nan,rfu['center']['p2'])
        rfu_p1 = np.where((sep < 0) | (rfu['sd']['p1'] > max_sd) |  (rfu['sd']['p2'] > max_sd) | (pop1 < min_pop), np.nan,rfu['center']['p1'])
        rfu_p2 = np.where((sep < 0) | (rfu['sd']['p1'] > max_sd) |  (rfu['sd']['p2'] > max_sd) | (pop1 > 1 - min_pop), np.nan,rfu['center']['p2'])
        rfu_p1 = np.where(np.isnan(rfu_p2),np.nan,rfu_p1) ## require both population fits to fit either
        rfu_p2 = np.where(np.isnan(rfu_p1),np.nan,rfu_p2) ## "" 
        rfu_p1 = np.where(np.isnan(pop1),np.nan,rfu_p1) ## require both population fits to fit either
        rfu_p2 = np.where(np.isnan(pop1),np.nan,rfu_p2) ## "" 
        rfu_avg = np.where((sep < 0) | (rfu['sd']['p1'] > max_sd) |  (rfu['sd']['p2'] > max_sd) | 
                        (pop1 < min_pop) | (pop1 > 1 - min_pop),rfu['center']['avg'],np.nan)
        rfu_avg = np.where(np.isnan(rfu_p1) & np.isnan(rfu_p2),rfu['center']['avg'],rfu_avg) #if missing pop1/2 fits use avg

        p1_sd = ((pop1 / (rfu_p2 - rfu_p1) )**2) * rfu['sd']['p1']**2 + \
                    ((rfu_avg - rfu_p1)/(rfu_p2 - rfu_p1)**2) * rfu['sd']['p2']**2 + \
                    ((1/(rfu_p2 - rfu_p1))**2) * rfu['sd']['avg']**2      #partial derivatives for error propogation
        

        if plot:
            axes[i*nset].line(resi,rfu['center']['p1'],fadedata=(rfu['low']['p1'],rfu['high']['p1']),labels="RFU pop1",color='cyan4') #more protected
            axes[i*nset].line(resi,rfu['center']['p2'],fadedata=(rfu['low']['p2'],rfu['high']['p2']),labels="RFU pop2",color='magenta') #more exchanged
            axes[i*nset].line(resi,rfu['center']['avg'],fadedata=(rfu['low']['avg'],rfu['high']['avg']),labels="RFU avg",color='dimgrey')

            marker_size = pop1.copy()
            marker_size[np.isnan(marker_size)] = 0
            
            if ~np.isnan(rfu_p1).all():
                axes[i*nset].scatter(resi,rfu_p1,labels="Fit pop1",color='green',s=(marker_size)**2,smin=5,smax=80,) #more protected
            else: axes[i*nset].scatter([-1],[-1],labels="Fit pop1",color='green',s=50) #need something to get label for legend
            if ~np.isnan(rfu_p2).all():
                axes[i*nset].scatter(resi,rfu_p2,labels="Fit pop2",color='darkred',s=2*((1-marker_size)**2),smin=5,smax=80,) #more exchanged
            else: axes[i*nset].scatter([-1],[-1],labels="Fit pop2",color='darkred',s=50) #need something to get label for legend
            if ~np.isnan(rfu_avg).all():
                axes[i*nset].scatter(resi,rfu_avg,labels="Fit AVG",edgecolors='black',facecolors='none',markeredgewidth=1.5,zorder=10)
            else: axes[i*nset].scatter([-1],[-1],labels="Fit AVG",color='black',facecolors='none',markeredgewidth=1.5,s=50) #need something to get label for legend
            
            axes[i*nset].format(ylim=(0,1.19),xlim=(0,max(resi)),ylocator=[0,0.25,0.5,0.75,1.0])
            axes[i*nset].format(ylabel="fractional D-uptake")
            axes[i*nset].format(xlabel='')
            #axes[i*nset].format(title=state+' '+time_label,titleloc= 'lower right',)
            axes[i*nset].text( max(resi)-3,1.25,state+' '+time_label,horizontalalignment='right',verticalalignment='bottom') #units of axes
            #axes[0].legend(loc='top',ncols=3);
            axes[i*nset].yaxis.set_label_coords(-0.049,0.5)

            if plotpop:
                #axes[i*nset].format(xtickloc='none')                
                axes[i*nset].format(xlabel="",xtickrange=(-1,-1))
                pop1_filt = np.where(np.isnan(rfu_avg),pop1,np.nan)
                axes[i*nset+1].scatter(resi,1-pop1,label="more exchanged",s=15,edgecolors='darkred',facecolors='white',linewidth=1)
                axes[i*nset+1].scatter(resi,1-pop1_filt,label="more exchanged",s=15,edgecolors='darkred',facecolors='magenta',linewidth=1)
                axes[i*nset+1].scatter(resi,pop1,label="more protected",s=15,edgecolors='green',facecolors='white',linewidth=1)
                axes[i*nset+1].scatter(resi,pop1_filt,label="more protected",s=15,edgecolors='green',facecolors='cyan4',linewidth=1)
                axes[i*nset+1].format(xlim=(0,max(resi)),ylim=(0,1.05),ylocator=[0,0.25,0.5,0.75,1.0])
                axes[i*nset+1].format(ylabel="Population")
                axes[i*nset+1].format(xlabel='Residue')
                axes[i*nset+1].yaxis.set_label_coords(-0.049,0.5)

        if return_fits:
            fits = pd.concat([pd.Series(rfu['center']['p1'],index=resi,name='RFU_pop1'),
                            pd.Series(rfu_p1,index=resi,name='Fit_pop1'),
                            pd.Series(rfu['sd']['p1'],index=resi,name='RFU1_sd'),
                            pd.Series(np.where(np.isnan(rfu_p1),np.nan,rfu['low']['p1']),index=resi,name='pop1_low'),
                            pd.Series(np.where(np.isnan(rfu_p1),np.nan,rfu['high']['p1']),index=resi,name='pop1_high'),
                            pd.Series(rfu['center']['p2'],index=resi,name='RFU_pop2'),
                            pd.Series(rfu_p2,index=resi,name='Fit_pop2'),
                            pd.Series(rfu['sd']['p2'],index=resi,name='RFU2_sd'),
                            pd.Series(np.where(np.isnan(rfu_p2),np.nan,rfu['low']['p2']),index=resi,name='pop2_low'),
                            pd.Series(np.where(np.isnan(rfu_p2),np.nan,rfu['high']['p2']),index=resi,name='pop2_high'),
                            pd.Series(rfu['center']['avg'],index=resi,name='RFU_avg'),
                            pd.Series(rfu_avg,index=resi,name='Fit_AVG'),
                            pd.Series(rfu['sd']['avg'],index=resi,name='RFUavg_sd'),
                            pd.Series(np.where(np.isnan(rfu_avg),np.nan,rfu['low']['avg']),index=resi,name='avg_low'),
                            pd.Series(np.where(np.isnan(rfu_avg),np.nan,rfu['high']['avg']),index=resi,name='avg_high'),
                            pd.Series(pop1,index=resi,name='Frac_pop1'),
                            pd.Series(p1_sd,index=resi,name='p1_sd')],axis=1)
            all_fits[use_time] = fits

    if plot:
        axes[-1].format(xlabel='Residue')
        axes[0].legend(ncols=3,edgecolor='white',bbox_to_anchor=(0.5,1.2),loc='center',); #loc='top',
        
        fig.set_facecolor('white')

        if savepath:
            #print("in savepath")
            if svg:
                fig.savefig(savepath,format='svg')
            else: 
                fig.savefig(savepath,format='pdf',dpi=600)
    if return_fits:
        return all_fits
    else: 
        return

def plot_comparison(fits,states=None,times=None,colors=None,savepath=None,svg=False,TD_time=1e6,UN_time=0):
    ''' 
    fits is a defaultdict as output from choose_fits(), fits[state][time]
    '''

    import proplot as pplt

    if states is None:
        states = [k for k in fits.keys()]
    if times is None:
        times = [k for k in fits[states[0]].keys()]
    if colors is None:
        colors="#000000 #66C2A5 #56B4E9 #7570B3 #E7298A".split()
    
    hspace = [1.0,6.0]*(len(times)-1)+[1.0]
    fig, axes = pplt.subplots(nrows=2*len(times),ncols=1,  axwidth="180mm", sharex=False, sharey=False, refaspect=5,hratios=[2,1]*len(times),hspace=hspace)

    for k,use_time in enumerate(times):

        if use_time == TD_time:
            time_label = 'FullDeut'
        elif use_time == UN_time:
            time_label = 'UnDeut'
        else:
            time_label = str(float(use_time))+'s'
        for i,state in enumerate(states):
            axes[2*k].line(fits[state][use_time].index,fits[state][use_time]['Fit_pop1'],fadedata=(fits[state][use_time]['pop1_low'],fits[state][use_time]['pop1_high']),label=state,color=colors[i])
            axes[2*k].line(fits[state][use_time].index,fits[state][use_time]['Fit_pop2'],fadedata=(fits[state][use_time]['pop2_low'],fits[state][use_time]['pop2_high']),label=None,color=colors[i])
            axes[2*k].scatter(fits[state][use_time].index,fits[state][use_time]['Fit_AVG'],fadedata=(fits[state][use_time]['avg_low'],fits[state][use_time]['avg_high']),s=10,label='Avg',edgecolors=colors[i],facecolors='none',linewidth=1)
            #paxes.scatter(fits[state][use_time].index,1-fits[state][use_time]['Frac_pop1'],label=state+" more exchanged",s=15,colors=colors[i],linewidth=1)
            pop1_filt = np.where(np.isnan(fits[state][use_time]['Fit_AVG']),fits[state][use_time]['Frac_pop1'],np.nan)
            axes[2*k+1].scatter(fits[state][use_time].index,fits[state][use_time]['Frac_pop1'],label=state+" less exchanged",s=15,edgecolors=colors[i],facecolors='#FFFFFF'+'20',linewidth=0.7,)#alpha=0.8)
            axes[2*k+1].scatter(fits[state][use_time].index,pop1_filt,label=state+" less exchanged",s=15,facecolors=colors[i],linewidth=0.7,)#alpha=0.8)
        #axes.text( 2,1.20,time_label,horizontalalignment='left',verticalalignment='bottom') #units of axes  ##above left
        axes[2*k].text( max(fits[state][use_time].index)-3,1.25,time_label,horizontalalignment='right',verticalalignment='bottom') #units of axes  ##above right
        #axes.text( max(fits[state][use_time].index)-3,1.12,time_label,horizontalalignment='right',verticalalignment='top') #units of axes ##inset top right
        axes[2*k].yaxis.set_label_coords(-0.048,0.5)
        axes[2*k+1].yaxis.set_label_coords(-0.048,0.5)
        axes[2*k].format(ylim=(0,1.19),xlim=(0,max(fits[state][use_time].index)),ylocator=[0,0.25,0.5,0.75,1.0])
        #paxes.format(title=time_label,titleloc= 'lower right',)
        axes[2*k].format(ylabel="fractional D-uptake",xlabel="",xtickrange=(-1,-1))
        #axes[2*k].format(xlabel='Residue')
        handles,labels=axes[2*k].get_legend_handles_labels()
        klabel = [kk*2 for kk in np.arange(len([k for k in states]))]
        # axes.legend(handles=[handles[i] for i in klabel],ncols=3, bbox_to_anchor=(1.0,1.2),loc='upper right', #units fraction of figure
        #             frame=True,edgecolor='white',facecolor='white',framealpha=0.8); #loc='top',    loc=(0,1)
        axes[2*k].legend(handles=[handles[j] for j in klabel],ncols=len(states), bbox_to_anchor=(0.5,1.1),loc='center', #units fraction of figure
                    frame=True,edgecolor='white',facecolor='white',framealpha=0.8); #loc='top',    loc=(0,1)
        
        axes[2*k+1].format(xlim=(0,max(fits[state][use_time].index)),ylim=(0.0,1.05))
        axes[2*k+1].format(ylabel="Population\n(protected)",ylocator=[0,0.25,0.5,0.75,1.0])
        axes[2*k+1].format(xlabel='Residue')
        #fig.suptitle("Big title woo")
        #paxes.format(title='Population of more \nexchanged state',titleloc='lr')
        #axes.legend(loc='top')
    fig.set_facecolor('white')

    if savepath:
        if svg:
            fig.savefig(savepath,format='svg')
        else:
            fig.savefig(savepath,format='pdf',dpi=600)


    return

def get_residues(rangetext):
    # assuming text is defined by ranges with '-' and comma separated values
    # not currently accepting negative residue values
    residues = []
    commagroups = rangetext.split(',')
    for i,c in enumerate(commagroups):      
        cc = [int(x) for x in c.split('-') if len(x) > 0]
        if len(cc) == 1: #single value specified
            residues += cc
        else: 
            residues += list(np.arange(cc[0],cc[1]+1))
    return residues

def group_fits(fitsdf,groups,max_sd = 0.5,min_pop = 0.05,z='50'):
    #cols = ['RFU_pop1','Fit_pop1','RFU1_sd','RFU_pop2','Fit_pop2','RFU2_sd','RFU_avg','Fit_AVG','RFUavg_sd']
    z_value ={'50':0.674,'68':1.0,'sd':1.0,'80':1.282,'90':1.645,'95':1.96,'98':2.326,'99':2.576}
    cols = ['RFU_pop1', 'RFU1_sd','RFU_pop2','RFU2_sd','RFU_avg','RFUavg_sd']
    groupfits_df = pd.DataFrame()
    for sample in fitsdf['sample'].unique():
        for group,v in groups.items():
            resgroup = get_residues(v)
            plot_data = fitsdf[(fitsdf['sample']==sample) & (fitsdf['r_number'].isin(resgroup))].copy()       
            #plot_data.loc[plot_data['exposure']==1e6,['RFU_pop1','RFU_pop2','Fit_pop1','Fit_pop2']] = np.nan #don't plot nonsense 2-state for TD
            #plot_data.loc[plot_data['exposure']==0,['RFU_avg','Fit_AVG']] = np.nan #don't draw line towards zero (looks misleading on log scale)
            use_data = plot_data.groupby(['exposure'])[cols].mean().reset_index() 

            use_data['pop1_low'] = use_data['RFU_pop1'] - z_value[z]*use_data['RFU1_sd']
            use_data['pop2_low'] = use_data['RFU_pop2'] - z_value[z]*use_data['RFU2_sd']
            use_data['avg_low'] = use_data['RFU_avg'] - z_value[z]*use_data['RFUavg_sd']

            use_data['pop1_high'] = use_data['RFU_pop1'] + z_value[z]*use_data['RFU1_sd']
            use_data['pop2_high'] = use_data['RFU_pop2'] + z_value[z]*use_data['RFU2_sd']
            use_data['avg_high'] = use_data['RFU_avg'] + z_value[z]*use_data['RFUavg_sd']

            use_data['sep'] = (use_data['RFU_pop2'] - use_data['RFU_pop1']) - (use_data['RFU2_sd'] + use_data['RFU1_sd'])*z_value[z]

            pop1 = (use_data['RFU_pop2'] - use_data['RFU_avg'])/(use_data['RFU_pop2'] - use_data['RFU_pop1'])
            pop1 = np.where(use_data['RFU_pop1'].isna(),np.nan,pop1)
            pop1 = np.where(use_data['RFU_pop2'].isna(),np.nan,pop1)
            pop1 = np.where((use_data['RFU_pop1'] < use_data['RFU_avg']) & (use_data['RFU_avg'] < use_data['RFU_pop2']),pop1,np.nan) #check avg between fits
            pop1 = pop1.clip(0,1)
            use_data['pop1'] = pop1
            use_data['pop2'] = 1 - pop1


            rfu_p1 = np.where((use_data['RFU_pop1'] > 1.1) | (use_data['RFU_pop1'] < -0.1),np.nan,use_data['RFU_pop1']) #excluded data still in stdev
            rfu_p2 = np.where((use_data['RFU_pop2'] > 1.1) | (use_data['RFU_pop2'] < -0.1),np.nan,use_data['RFU_pop2']) 

            rfu_p1 = np.where((use_data['sep'] < 0) | (use_data['RFU1_sd'] > max_sd) |  (use_data['RFU2_sd'] > max_sd) | (pop1 < min_pop), np.nan,use_data['RFU_pop1'])
            rfu_p2 = np.where((use_data['sep'] < 0) | (use_data['RFU1_sd'] > max_sd) |  (use_data['RFU2_sd'] > max_sd) | (pop1 > 1 - min_pop), np.nan,use_data['RFU_pop2'])
            rfu_p1 = np.where(np.isnan(rfu_p2),np.nan,rfu_p1) ## require both population fits to fit either
            rfu_p2 = np.where(np.isnan(rfu_p1),np.nan,rfu_p2) ## "" 
            rfu_p1 = np.where(np.isnan(pop1),np.nan,rfu_p1) ## require both population fits to fit either
            rfu_p2 = np.where(np.isnan(pop1),np.nan,rfu_p2) ## "" 
            rfu_avg = np.where((use_data['sep'] < 0) | (use_data['RFU1_sd'] > max_sd) |  (use_data['RFU2_sd'] > max_sd) | 
                            (pop1 < min_pop) | (pop1 > 1 - min_pop),use_data['RFU_avg'],np.nan)
            rfu_avg = np.where(np.isnan(rfu_p1) & np.isnan(rfu_p2),use_data['RFU_avg'],rfu_avg) #if missing pop1/2 fits use avg

            p1_sd = ((pop1 / (use_data['RFU_pop2'] - use_data['RFU_pop1']) )**2) * use_data['RFU1_sd']**2 + \
                    ((use_data['RFU_avg'] - use_data['RFU_pop1'])/(use_data['RFU_pop2'] - use_data['RFU_pop1'])**2) * use_data['RFU2_sd']**2 + \
                    ((1/(use_data['RFU_pop2'] - use_data['RFU_pop1']))**2) * use_data['RFUavg_sd']**2      #partial derivatives for error propogation

            use_data['p1_sd'] = np.sqrt(p1_sd)
            use_data['Fit_pop1'] = rfu_p1
            use_data['Fit_pop2'] = rfu_p2
            use_data['Fit_AVG'] = rfu_avg

            use_data['sample'] = sample
            use_data['Group'] = group 

            groupfits_df = pd.concat([groupfits_df,use_data]).reset_index(drop=True)
    return groupfits_df



def plot_resred(hdxm,states=None,times=None,seq=None,savepath=None,svg=False,plotZero=True,TD_time=1e6,UN_time=0):
    '''
    Just plot the redundancy and resolution color bars for the specified times and states
        (note that all times should have the same red/res if coming out of the hdxm without further peptide filtering)
    '''

    import proplot as pplt
    import matplotlib as mpl
    import colorcet as cc

    if states is None:
        states = list(hdxm.keys())
    if times is None:
        times = sorted(hdxm[states[0]].timepoints)
    if not plotZero:
        times.remove(0)
    
    z_value ={'50':0.674,'68':1.0,'sd':1.0,'80':1.282,'90':1.645,'95':1.96,'98':2.326,'99':2.576} #confidence intervals multiplier 
  
    nset = 2 # redundancy, resolution 
    ncols = 20
    nrows = len(times)*len(states)*nset
    nfigs = nrows
    array=[]
    for i in range(1,nfigs+1):
        array += [np.repeat(i,ncols)]# [[i,i]]
    hratios = [1,1]*len(times)*len(states)

    fig, axes = pplt.subplots(array,  axheight="4mm",axwidth="180mm",  #refaspect=10,
                                wspace=(0), hspace=([0.3,4.0]*(len(times)*len(states)-1)+[0.3]),
                                sharex=False, sharey=False, hratios=hratios)
    
    vmin = 0
    vmax = 20
    red_scale = 2
    kwargs = dict(levels=pplt.arange(0, vmax, 1),) #extendsize=2.5, extendrect=True)
    res_kwargs = dict(levels=pplt.arange(0, vmax, 1),)

    k=0
    for i, use_time in enumerate(times):
        if use_time == TD_time:
            time_label = 'FullDeut'
        elif use_time == UN_time:
            time_label = 'UnDeut'
        else:
            time_label = str(float(use_time))+'s'

        for j,mutant in enumerate(states):
            time_idx = np.where(hdxm[mutant].timepoints == use_time)[0][0]

            x_res = hdxm[mutant][time_idx].r_number

        
            redundancy = hdxm[mutant].coverage.X.sum(axis=0).astype(float)
            resolution = np.repeat(hdxm[mutant].coverage.block_length, hdxm[mutant].coverage.block_length)
            resolution = resolution.astype(float)
            resolution[redundancy == 0] = np.nan
            redundancy[redundancy == 0] = np.nan

            if seq is not None:
                if isinstance(seq,dict):
                    seq = seq[mutant]
                pro_index = [index+1 for index,s in enumerate(seq) if (s=="P")]
                pro_index = sorted(set(pro_index) & set(hdxm[mutant].coverage.r_number))
                redundancy = pd.Series(redundancy,index = hdxm[mutant].coverage.r_number)
                resolution = pd.Series(resolution,index = hdxm[mutant].coverage.r_number)
                redundancy[pro_index] = np.nan
                resolution[pro_index] = np.nan

            #print("resolution",resolution)
            res_cmap = cc.cm['CET_L4_r'].copy() #.set_under() .set_over()
            res_cmap = truncate_colormap(res_cmap,0,0.7)
            # res_cmap = cc.cm['coolwarm'].copy()
            # res_cmap = truncate_colormap(res_cmap,0,1.0)
            res_cmap.set_bad(color='k')
            red_cmap = cc.cm['CET_D1A_r'].copy()
            red_cmap = truncate_colormap(red_cmap,0.1,0.9)
            red_cmap.set_bad(color='k')

            red = plot_bar(axes[k*nset+0],hdxm[states[0]].coverage.r_number,redundancy,
                    red_cmap,pplt.Norm("linear", vmin=vmin, vmax=vmax*red_scale),height=1.0,**kwargs)
            axes[k*nset+0].format(xtickrange=(-2,-1),title=mutant+', '+time_label)
            res = plot_bar(axes[k*nset+1],hdxm[states[0]].coverage.r_number,resolution, 
                    res_cmap,pplt.Norm("linear", vmin=vmin, vmax=vmax),height=1.0,**res_kwargs)               

            axes[k*nset+0].set_ylabel("red.",rotation=0,)#loc='center',labelpad=10)
            axes[k*nset+0].yaxis.set_label_coords(-.02,-.05)
            axes[k*nset+0].format(xtickloc='none')
            axes[k*nset+1].set_ylabel("res.",rotation=0,)#loc='center',labelpad=10)
            axes[k*nset+1].yaxis.set_label_coords(-.02,-.05)
            #print("k time state",k,use_time,mutant)
            k += 1

    #need to implement https://proplot.readthedocs.io/en/latest/subplots.html?highlight=spacing#Spacing-and-tight-layout
    
    sm = mpl.cm.ScalarMappable(cmap= mpl.colors.ListedColormap([res_cmap.get_bad()]))

    cbar_width = 0.4
    gap = 1
    cols = [[1,int(ncols*cbar_width)], 
                [int(ncols*cbar_width)+1+gap,int(ncols*cbar_width*2)+1], 
                [int(ncols*cbar_width*2)+2+gap,ncols]]

    fig.colorbar(red, label="Redundancy (peptides including replicates)",width=0.1,pad=1.5,
                 loc='b',length=1.0,col=cols[0],ticks=np.arange(0, vmax+1, 5))# **kwargs)
    fig.colorbar(res, label="Resolution (residues)",width=0.1,
                 loc='b',length=1.0,col=cols[1],ticks=np.arange(0, vmax+1, 5))# **kwargs)
    pcb = fig.colorbar(sm,width=0.1,loc='b',col=cols[2],length=0.1,ticks=[],)
    pcb.set_label(label="No Coverage\nor Proline",labelpad=8)
    
    #axes[0].legend(loc="t", ncols=legendcols)
    axes.format(ylim=(0,1.1),)#xgrid=True,ygrid=True,ytickrange=(0,1.1))
    
    #axes.format(xlim=(0,180))
    fig.set_facecolor('white')
    axes[-1].format(xlabel="Residue")

    # save_path = os.path.join(project_dir,'B5_phospho_2popFiltShade_RFU_residueavg_plots_'+date+'.pdf')
    if savepath:
        if svg:
            fig.savefig(savepath,format='svg')
        else:
            fig.savefig(savepath,format='pdf',dpi=600)

    return #fig


def filter_range(hdxm,peprange,startcol='start',endcol='end',nterm_exch=2):
    df = hdxm.copy()
    peprange = [peprange] if not isinstance(peprange,list) else peprange
    df = df[(df[startcol] <= peprange[-1]-nterm_exch) & (peprange[0] <= df[endcol])]
    return df




def collate_pdfs(pdf_dir,project=None,samples=None,save_dir=None,addlabel=True):
    '''
    cobbles together a proper order to merge all the generated pdfs after the fact

    for proj in projects:
        for sample in samples[proj]:
            pdf_dir = os.path.join(batch_dir[proj],sample+'_run')
            #print(pdf_dir)
            collate_pdfs(pdf_dir,project='HSPB5_'+sample+'_'+proj,save_dir=batch_dir[proj])
    '''

    #!pip install PyMuPDF
    import pandas as pd
    import numpy as np
    import fitz
    import os
    from datetime import datetime

    now = datetime.now()
    date = now.strftime("%d%b%Y")

    def categorize(df,types,search_col,new_col):
        for i,type in enumerate(types):
            df[type] = np.where(df[search_col].str.contains(type),str(i)+type,'')
            df[new_col]=''
        # for type in types:
        #     df[new_col] += pdf_info[type]
        df[new_col] = df.apply(lambda x: max(x[types],key=len),axis=1)
        df = df.drop(types,axis=1)
        return df

    def add_header(mfile,header_text):
        for page in mfile:
            prect = page.rect
            header_rect = fitz.Rect(0, 10, prect.width-10,30)  # height 20 points
            page.insert_textbox(header_rect, header_text,
                                fontname="hebo", color= fitz.pdfcolor["blue"],    
                                align=fitz.TEXT_ALIGN_RIGHT)
        return mfile
    
    output_types = ['IndFits','BootFits','ndeut']

    pdfs = [ f for f in os.listdir(pdf_dir) if f[-4:]=='.pdf'  ]
    pdfs = [ p for p in pdfs if p[0:5]=='hxex_']
    pdf_info = pd.DataFrame({'file':pdfs})

    pdf_info = categorize(pdf_info,output_types,'file','type')
    if samples is not None: 
        samples = samples if isinstance(samples,list) else [samples]
        pdf_info = categorize(pdf_info,samples,'file','sample')
    else: pdf_info['sample'] = ''
    if project is not None: proj_str = project+'_'
    else: proj_str = ''
    if save_dir is None: save_dir = pdf_dir

    pdf_info[['start_res','end_res']] = ''
    #print("pdf_info['file']",pdf_info['file'])
    pdf_info['start_res'] = pdf_info['file'].str[:-4].str.split('-',n=1,).str[0].str[-4:].apply(pd.to_numeric,errors='coerce').astype('Int64')
    pdf_info['end_res'] = pdf_info['file'].str[:-4].str.split('-',n=1,).str[1].str[0:4].apply(pd.to_numeric,errors='coerce').astype('Int64')
    pdf_info['charge'] = pdf_info['file'].str[:-4].str.split('_').str[-3].str[1:].apply(pd.to_numeric,errors='coerce').astype('Int64')
    #pdf_info['data_id'].str[:-4].str.split('_').str[1].astype('Int64')
    pdf_info = pdf_info.sort_values(['start_res','end_res','sample','charge','type']).reset_index(drop=True)
    pdf_info['type'] = pdf_info['type'].str[1:]
    pdf_info['sample'] = pdf_info['sample'].str[1:]
    #display(pdf_info)

    sorted_pdfs = pdf_info['file'].values
    sorted_pdfs = [s for s in sorted_pdfs if '_all_' not in s]
    sorted_pdfs_Ind = [s for s in sorted_pdfs if 'Ind' in s]
    sorted_pdfs_BootFits = [s for s in sorted_pdfs if 'BootFits' in s]
    sorted_pdfs_ndeut = [s for s in sorted_pdfs if 'ndeut' in s]
    #sorted_pdfs

    

    result = fitz.open()
    for pdf in sorted_pdfs:        
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)
            result.insert_pdf(mfile)
    result.save(os.path.join(save_dir,'hxex_all_plots_'+proj_str+date+'.pdf'))

    result = fitz.open()
    for pdf in sorted_pdfs_Ind:
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)
            result.insert_pdf(mfile)            
    result.save(os.path.join(save_dir,'hxex_all_IndFits_'+proj_str+date+'.pdf'))

    result = fitz.open()
    for pdf in sorted_pdfs_BootFits:
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)            
            result.insert_pdf(mfile)
    result.save(os.path.join(save_dir,'hxex_all_BootFits_'+proj_str+date+'.pdf'))

    result = fitz.open()
    for pdf in sorted_pdfs_ndeut:
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)            
            result.insert_pdf(mfile)
    result.save(os.path.join(save_dir,'hxex_all_ndeutBoot_'+proj_str+date+'.pdf'))

    return


#Functions for zoom view of plots to autoscale y-axis based on xlims 
#https://stackoverflow.com/questions/29461608/fixing-x-axis-scale-and-autoscale-y-axis
def autoscale(ax=None, axis='y', margin=0.1):
    '''Autoscales the x or y axis of a given matplotlib ax object
    to fit the margins set by manually limits of the other axis,
    with margins in fraction of the width of the plot

    Defaults to current axes object if not specified.
    '''

    import matplotlib as mpl

    def calculate_new_limit(fixed, dependent, limit):
        '''Calculates the min/max of the dependent axis given 
        a fixed axis with limits
        '''
        if len(fixed) > 2:
            mask = (fixed>limit[0]) & (fixed < limit[1])
            window = dependent[mask]
            low, high = window.min(), window.max()
        else:
            low = dependent[0]
            high = dependent[-1]
            if low == 0.0 and high == 1.0:
                # This is a axhline in the autoscale direction
                low = np.inf
                high = -np.inf
        return low, high

    def get_xy(artist):
        '''Gets the xy coordinates of a given artist
        '''
        if "Collection" in str(artist):
            x, y = artist.get_offsets().T
        elif "Line" in str(artist):
            x, y = artist.get_xdata(), artist.get_ydata()
        else:
            raise ValueError("This type of object isn't implemented yet")
        return x, y


    if ax is None:
        ax = mpl.plt.gca()
    newlow, newhigh = np.inf, -np.inf

    for artist in ax.collections + ax.lines:
        x,y = get_xy(artist)
        if axis == 'y':
            setlim = ax.set_ylim
            lim = ax.get_xlim()
            fixed, dependent = x, y
        else:
            setlim = ax.set_xlim
            lim = ax.get_ylim()
            fixed, dependent = y, x

        low, high = calculate_new_limit(fixed, dependent, lim)
        newlow = low if low < newlow else newlow
        newhigh = high if high > newhigh else newhigh

    margin = margin*(newhigh - newlow)

    setlim(newlow-margin, newhigh+margin)