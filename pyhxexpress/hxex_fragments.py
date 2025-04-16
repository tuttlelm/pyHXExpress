''' 
supported ion_type in modification designation is from pyteomics

pyteomics std_ion_comp
https://github.com/levitsky/pyteomics/blob/master/pyteomics/mass/mass.py
'''

import pandas as pd
import numpy as np

ion_frag_dict = {
    'M':      'parent',  #Composition(formula=''),
    'M-H2O':  'parent',  #Composition(formula='H-2O-1'),
    'M-NH3':  'parent',  #Composition(formula='N-1H-3'),
    'a':      'N-term',  #Composition(formula='H-2O-1' + 'C-1O-1'),
    'a-H2O':  'N-term',  #Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':  'N-term',  #Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b':      'N-term',  #Composition(formula='H-2O-1'),
    'b-H2O':  'N-term',  #Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':  'N-term',  #Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':      'N-term',  #Composition(formula='H-2O-1' + 'NH3'),
    'c-1':    'N-term', #Composition(formula='H-2O-1' + 'NH3' + 'H-1'),
    'c-dot':  'N-term',  #Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+1':    'N-term',  #Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+2':    'N-term',  #Composition(formula='H-2O-1' + 'NH3' + 'H2'),
    'c-H2O':  'N-term',  #Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3':  'N-term',  #Composition(formula='H-2O-1'),
    'x':      'C-term',  #Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O':  'C-term',  #Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3':  'C-term',  #Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y':      'C-term',  #Composition(formula=''),
    'y-H2O':  'C-term',  #Composition(formula='H-2O-1'),
    'y-NH3':  'C-term',  #Composition(formula='N-1H-3'),
    'z':      'C-term',  #Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z-dot':  'C-term',  #Composition(formula='H-2O-1' + 'ON-1'),
    'z+1':    'C-term',  #Composition(formula='H-2O-1' + 'ON-1H1'),
    'z+2':    'C-term',  #Composition(formula='H-2O-1' + 'ON-1H2'),
    'z+3':    'C-term',  #Composition(formula='H-2O-1' + 'ON-1H3'),
    'z-H2O':  'C-term',  #Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3':  'C-term',  #Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
    }


def get_fragments(metadf,user_frags=None):
    '''
    if frags not specified, must be column named 'frags' in the (parent) metadf to generate fragments for each entry
    specify as 'ion_type:fragments:charges' where :charges is optional
    e.g. 'b-NH3:1,2,7 x:ALL y:1,2:1,2'
     '''
    def get_pepfrag(trunc_type):
        if (trunc_type == 'N-term'):
            start_idx = 0
            end_idx = trunc
        elif (trunc_type == 'C-term'):
            start_idx = len(peptide) - trunc 
            end_idx = len(peptide)
        else:
            start_idx = 0
            end_idx = len(peptide)
        return start_idx,end_idx

    frags_df = pd.DataFrame()

    for index, row in metadf.iterrows():
        if 'modification' in row.keys():
            mods = row['modification'].split()
            it = (list(filter(lambda x: x.startswith('ion_type'),mods)))
            new_mods = ' '.join(list(filter(lambda x: x not in it,mods)))+' '

        else: new_mods = ''
        peptide = row.peptide

        if user_frags is None:
            if 'frags' in metadf.columns:
                frags = row.frags.split()               
            else: 
                print("must specify fragments to generate")
        else: 
            frags = user_frags.split()

        for fx in frags:
            ionx = fx.split(':')[0]
            truncx = fx.split(':')[1].split(',') if len(fx.split(':'))>1 else ['ALL']
            chargex = (fx.split(':')[2].split(',')) if len(fx.split(':'))==3 else []
            #print(ionx,truncx,chargex)
            if ionx in ion_frag_dict.keys():
                trunc_type = ion_frag_dict[ionx]
            else:
                print("ion_type not recognized")
                continue
            max_charge = min(len(peptide)//2,row.charge)
            if len(chargex) > 0:
                if 'ALL' in chargex:
                    charges = 'ALL'
                else: 
                    try: charges = [int(c) for c in chargex if int(c) <= max_charge]
                    except: 
                        print("syntax error in charges, using ALL")
                        charges = 'ALL'
            else: charges = 'ALL'
            if charges == 'ALL':
                #max_charge = min(len(peptide)//2,row.charge)
                charges = np.arange(1,max_charge+1)

            for charge in charges:
                if len(truncx) > 0:
                    if 'ALL' in truncx:
                        truncs = 'ALL'
                    else: 
                        try: truncs = [int(t) for t in truncx]
                        except: 
                            print("syntax error in truncations, using ALL")
                            truncs = 'ALL'
                #else: truncs = 'ALL'
                if truncs == 'ALL':
                    min_trunc = 0 #
                    truncs = np.arange(min_trunc,len(row.peptide)-1)+1
                #print(truncs)
                for trunc in truncs:
                    frag_df = pd.DataFrame()
                    if trunc > len(peptide): 
                        print("specified fragment longer than peptide:",trunc)
                        continue
                    
                    sidx,eidx = get_pepfrag(trunc_type)
                    pep = peptide[sidx:eidx]
                    
                    frag_df = row.copy()
                    frag_df['parent'] = row['sample']
                    frag_df['sample'] = row['sample'] + '_' + ionx+'_'+str(trunc)+'_z'+str(charge)
                    frag_df['start_seq'] = row.start_seq + sidx
                    frag_df['end_seq'] = row.start_seq + eidx - 1
                    frag_df['parent_range'] = frag_df['peptide_range']
                    frag_df['peptide_range'] = '-'.join([str(row.start_seq+sidx),str(row.start_seq+eidx-1)])
                    frag_df['peptide'] = pep
                    frag_df['charge'] = charge
                    frag_df['modification'] = new_mods+'ion_type:'+ionx
                    frag_df['frags']=ionx+':'+str(trunc)
                    frag_df['parent_peptide']=peptide
                    frag_df['parent_charge']=row.charge
                    frag_df['parent_id']=index
                    
                    frags_df = pd.concat([frags_df,pd.DataFrame([frag_df])],ignore_index=True)

    return frags_df

    
