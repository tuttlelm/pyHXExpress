''' 
supported ion_type in modification designation is from pyteomics
get_fragments currently only supporting c and z-dot types 

pyteomics std_ion_comp
https://github.com/levitsky/pyteomics/blob/master/pyteomics/mass/mass.py

std_ion_comp.update({
    'M':        Composition(formula=''),
    'M-H2O':    Composition(formula='H-2O-1'),
    'M-NH3':    Composition(formula='N-1H-3'),
    'a':        Composition(formula='H-2O-1' + 'C-1O-1'),
    'a-H2O':    Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':    Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b':        Composition(formula='H-2O-1'),
    'b-H2O':    Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':    Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':        Composition(formula='H-2O-1' + 'NH3'),
    'c-1':      Composition(formula='H-2O-1' + 'NH3' + 'H-1'),
    'c-dot':    Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+1':      Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+2':      Composition(formula='H-2O-1' + 'NH3' + 'H2'),
    'c-H2O':    Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3':    Composition(formula='H-2O-1'),
    'x':        Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O':    Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3':    Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y':        Composition(formula=''),
    'y-H2O':    Composition(formula='H-2O-1'),
    'y-NH3':    Composition(formula='N-1H-3'),
    'z':        Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z-dot':    Composition(formula='H-2O-1' + 'ON-1'),
    'z+1':      Composition(formula='H-2O-1' + 'ON-1H1'),
    'z+2':      Composition(formula='H-2O-1' + 'ON-1H2'),
    'z+3':      Composition(formula='H-2O-1' + 'ON-1H3'),
    'z-H2O':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
    })

''';

import pandas as pd
import numpy as np

def get_fragments(metadf,frags=None):
    '''
    if frags not specified, must be column named 'frags' in the (parent) metadf to generate fragments for each entry
    formatting should be similar to the 'modification column' as a space deliminated dictionary
    e.g. 'c:ALL' or 'c:1,2,7 M:ALL' **note that ion types and chops need to be finalized

    ## need to work up including c1: or c2: or all possible charge subsets

    otherwise, if all entries should have the same fragment pattern can specify with frags
    frags is a dictionary of ion_type and a list of substates, e.g. frags = {'c':'ALL'} or {'c+1':'5,7'}
    'c':'ALL' will do all C-term truncations with ion_type = c
    alternatively specify {'c':'1,2,7']} for specific truncations
    '''
    frags_df = pd.DataFrame()

    for index, row in metadf.iterrows():
        if 'modification' in row.keys():
            mods = row['modification'].split()
            it = (list(filter(lambda x: x.startswith('ion_type'),mods)))
            new_mods = ' '.join(list(filter(lambda x: x not in it,mods)))+' '

        else: new_mods = ''
        peptide = row.peptide

        if frags is None:
            if 'frags' in metadf.columns:
                frags = {x.split(':')[0]:(x.split(':')[1]) for x in row.frags.split()}
            else: 
                print("must specify fragments to generate")
        #print("frags:",frags)
        c_keys = [k for k in frags.keys() if k.startswith('c')]
        z_keys = [k for k in frags.keys() if k.startswith('z')]
        #print(row.file,row.peptide)
        for c_key in c_keys:
            #print("c_key:",c_key)
            ion_type = 'c'
           

            c_charge = c_key.split('c')[1]
            #print("c_charge", c_charge)
            if len(c_charge) > 0:
                charges = [int(c_charge[1:]) ]
            else:
                max_charge = min(len(peptide)//2,row.charge)
                charges = np.arange(1,max_charge+1)

            for charge in charges:
                if frags[c_key] == 'ALL':
                    min_trunc = 0 # pepfrag uses charge*2-1-1 # but MS-Prospect starts at charge -1 (-1 is for index)# dunno
                    truncs = np.arange(min_trunc,len(row.peptide)-1)+1
                else:
                    truncs = list(map(int,frags[c_key].split(',')))
                for trunc in truncs:
                    frag_df = pd.DataFrame()
                    if trunc > len(peptide): print("specified fragment longer than peptide")  
                    pep = peptide[0:trunc]
                    
                    frag_df = row.copy()
                    frag_df['parent'] = row['sample']
                    frag_df['sample'] = row['sample']+'_z'+str(row['charge'])+ '_' + ion_type+str(trunc)
                    frag_df['end_seq'] = row.start_seq + trunc - 1
                    frag_df['parent_range'] = frag_df['peptide_range']
                    frag_df['peptide_range'] = '-'.join([str(row.start_seq),str(row.start_seq+trunc-1)])
                    frag_df['peptide'] = pep
                    frag_df['charge'] = charge
                    frag_df['modification'] = new_mods+'ion_type:'+ion_type
                    frag_df['frags']=ion_type+':'+str(trunc)
                    frag_df['parent_peptide']=peptide
                    frag_df['parent_charge']=row.charge
                    frag_df['parent_id']=index
                    
                    frags_df = pd.concat([frags_df,pd.DataFrame([frag_df])],ignore_index=True)
        for z_key in z_keys:
            #print("c_key:",c_key)
            ion_type = 'z-dot'
            

            z_charge = z_key.split('z')[1]
            #print("c_charge", c_charge)
            if len(z_charge) > 0:
                charges = [int(z_charge[1:]) ]
            else:
                max_charge = min(len(peptide)//2,row.charge)
                charges = np.arange(1,max_charge+1)

            for charge in charges:
                if frags[z_key] == 'ALL':
                    min_trunc = 0 # pepfrag uses charge*2-1-1 # but MS-Prospect starts at charge -1 (-1 is for index)# dunno
                    truncs = np.arange(min_trunc,len(row.peptide)-1)+1
                else:
                    truncs = list(map(int,frags[z_key].split(',')))
                for trunc in truncs:
                    frag_df = pd.DataFrame()
                    if trunc > len(peptide): print("specified fragment longer than peptide")
                    start_idx = len(peptide) - trunc  
                    pep = peptide[start_idx:]
                    
                    frag_df = row.copy()
                    frag_df['parent'] = row['sample']
                    frag_df['sample'] = row['sample']+'_z'+str(row['charge'])+ '_' + ion_type+str(trunc)
                    frag_df['start_seq'] = row.end_seq-trunc+1
                    frag_df['parent_range'] = frag_df['peptide_range']
                    frag_df['peptide_range'] = '-'.join([str(row.end_seq-trunc+1),str(row.end_seq)])
                    frag_df['peptide'] = pep
                    frag_df['charge'] = charge
                    frag_df['modification'] = new_mods+'ion_type:'+ion_type
                    frag_df['frags']=ion_type+':'+str(trunc)
                    frag_df['parent_peptide']=peptide
                    frag_df['parent_charge']=row.charge
                    frag_df['parent_id']=index
                    frags_df = pd.concat([frags_df,pd.DataFrame([frag_df])],ignore_index=True)
    return frags_df
