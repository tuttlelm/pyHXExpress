''' 
supported ion_type in modification designation is from pyteomics

pyteomics std_ion_comp
https://github.com/levitsky/pyteomics/blob/master/pyteomics/mass/mass.py
'''

from pyteomics import parser
from pyteomics import mass
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


def get_fragments(metadf,user_frags=None, limit_charge=True):
    '''
    if frags not specified, must be column named 'frags' in the (parent) metadf to generate fragments for each entry
    specify as 'ion_type:fragments:charges' where :charges is optional
    e.g. 'b-NH3:1,2,7 x:ALL y:1,2:1,2'
     '''
    def get_pepfrag(peptide,trunc_type,trunc):
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
    
    def expand_range(values,start,end):
        expanded_list = []
        if 'ALL' in values:
            expanded_list += list(np.arange(int(start),int(end)+1))
        else:
            for itrunc in values:
                if len(itrunc.split('-'))>1:
                    ltrunc = itrunc.split('-')[0]               
                    rtrunc = itrunc.split('-')[1]
                    ltrunc = int(ltrunc) if len(ltrunc) > 0 else int(start)
                    rtrunc = int(rtrunc) if len(rtrunc) > 0 else int(end)
                    expanded_list += list(np.arange(ltrunc,rtrunc+1))
                else:
                    expanded_list += [int(itrunc)]
        return set(expanded_list)

    frags_df = pd.DataFrame()

    for index, row in metadf.iterrows():
        if 'modification' in row.keys():
            mods = row['modification'].split()
            it = (list(filter(lambda x: x.startswith('ion_type'),mods)))
            new_mods = ' '.join(list(filter(lambda x: x not in it,mods)))+' '

        else: new_mods = ''
        sequence = row.peptide
        
        pep_parts = parser.parse(sequence)
        if '-' in pep_parts[0]:
            pep_nterm = pep_parts[0]
            pep_start = 1
        else: 
            pep_nterm = ''
            pep_start = 0
        if '-' in pep_parts[-1]:
            pep_cterm = pep_parts[-1]
            pep_end = -1
        else: 
            pep_cterm = ''
            pep_end = len(sequence)+1

        peptide = ''.join(pep_parts[pep_start:pep_end]) #ditch C/N modifications for truncations
        # print("sequence:",sequence)
        # print("peptide:",peptide)

        if user_frags is None:
            if 'frags' in metadf.columns:
                frags = row.frags.split()               
            else: 
                print("must specify fragments to generate")
        else: 
            frags = user_frags.split()

        for fx in frags:
            ionx = fx.split(':')[0]
            truncs = fx.split(':')[1].split(',') if len(fx.split(':'))>1 else ['ALL']
            charges = (fx.split(':')[2].split(',')) if len(fx.split(':'))==3 else ['ALL']
            truncx = expand_range(truncs,row.start_seq,row.end_seq)          

            if ionx in ion_frag_dict.keys():
                trunc_type = ion_frag_dict[ionx]
            else:
                print("ion_type not recognized")
                continue
            
            #lose any truncated N/C modifications
            if trunc_type == "N-term":  
                pep_c = ''
            else: pep_c = pep_cterm
            if trunc_type == "C-term":
                pep_n = ''
            else: pep_n = pep_nterm
            if trunc_type == "parent":
                truncx = [row.end_seq - row.start_seq + 1]
                max_charge = row.charge #     
            else:
                if limit_charge:
                    max_charge = row.charge-1 #parent charge - 1
                else:
                    max_charge = row.charge #                  
            chargex = expand_range(charges,1,max_charge)

            for charge in chargex:
                for trunc in truncx:
                    frag_df = pd.DataFrame()
                    if trunc > len(peptide): 
                        print("specified fragment longer than peptide:",trunc)
                        continue
                    
                    sidx,eidx = get_pepfrag(peptide,trunc_type,trunc)
                    pep = peptide[sidx:eidx]

                    frag_df = row.copy()
                    frag_df['parent'] = row['sample']
                    frag_df['sample'] = row['sample'] + '_' + ionx+'_'+str(trunc)+'_z'+str(charge)
                    frag_df['start_seq'] = row.start_seq + sidx
                    frag_df['end_seq'] = row.start_seq + eidx - 1
                    frag_df['parent_range'] = frag_df['peptide_range']
                    frag_df['peptide_range'] = '-'.join([str(row.start_seq+sidx),str(row.start_seq+eidx-1)])
                    frag_df['peptide'] = pep_n+pep+pep_c
                    frag_df['charge'] = charge
                    frag_df['modification'] = new_mods+'ion_type:'+ionx
                    frag_df['frags']=ionx+':'+str(trunc)
                    frag_df['parent_peptide']=sequence
                    frag_df['parent_charge']=row.charge
                    frag_df['parent_id']=index
                    
                    frags_df = pd.concat([frags_df,pd.DataFrame([frag_df])],ignore_index=True)

    return frags_df

    
