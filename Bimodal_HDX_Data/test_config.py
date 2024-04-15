import os
from datetime import datetime
now = datetime.now()
date = now.strftime("%d%b%Y")
##in notebook use: 
# import pyhxexpress
# import test_config as config
# pyhxexpress.config = config

##########################################
'''Settings for Test Data Sets'''
##########################################

WRITE_PARAMS = True #save the params to hdxms_params_$.py file in Data_DIR, can then be read in as PARAMS_FILE 
Allow_Overwrite = True #don't create a new filename if file already exists

Test_Data = True
Data_Type = 1 
Save_Spectra = False
process_ALL = True
Data_DIR = 'C:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\pyHXExpress\\Bimodal_HDX_Data'
#Data_DIR = '/home/tuttle/data/HDX-MS/pyHXExpress/Bimodal_HDX_Data'
Read_Spectra_List = True #get metadf from Metadf_File
Metadf_File = "hdxms_testsets_metadf.csv"              
Output_DIR = os.path.join(Data_DIR,'output_'+str(date),'')
#if not os.path.exists(Output_DIR): os.makedirs(Output_DIR)
Preset_Pops = False #Use predetermined number of populations to fit curves, overrides Min/Max_Pops if given
Preset_Pops_File = os.path.join(Data_DIR,"test_configdf_26feb2024.csv")
User_mutants = ['']
User_peptides = ['']

Generate_Plots = True
Hide_Figure_Output = False #True Recommended when processing lots of data. 
SVG = False # also save figures as an svg file, slow, but better for making figures 

BestFit_of_X = 3
Residual_Cutoff = 3.3e-5

Nboot = 20 # number of individual fits to perform, using n_best_curves from initial round of fits
setNoise = 5000 #if noise value is known, specify instead of estimating as Y_ERR % of avg Un+TD peaks, overrides Y_ERR
Y_ERR = 2.0 #Percent random error applied during boot as y*+np.random.normal(0,yerr), 0.0 for NoNoise, use ~0.5% for noise added
            # the absolute Noise value is then Y_ERR * avg(maxInt of Un and TD)
            # this is a very rough way to give a consistent Noise value throughout a dataset. 
Dfrac = 0.90
FullDeut_Time = 1e6 # dummy timept for fully deuterated (TD) sample
Nterm_subtract = 1 # number of N-term residues to remove from possible exchanging NH's (usually 1 or 2)
                   # This will mostly affect the corrected Dabs value as it is scaled to make TD = theoretical Nex
Zero_Filling = 3 # large value to try to fit 3 binomials on short peptides, not generally recommended.
Peak_Resolution = 70.0 #ppm, sensitivity of peak picker to expected m/z centers  #LeuEnk wants higher, pep122 wants lower ...
Binomial_dCorr = True # fit n=1 binomial for UN/TD to calculate d_corr and back exchange values
Env_threshold = 0.1 #find envelope width at Env_threshold * Intensity_max
Limit_by_envelope = False # only fit up to n = int(z*env/3*Env_limit - 2/3) 
Env_limit = 0.95 #used if Limit_by_envelope = True, rough measure to constrain n_curves fit according to data width & num fit params
Min_Pops = 1 # Min_Pops to Max_Pops sets the range of populations to fit, set to same value to force single population
Max_Pops = 3 # maximum number of underlying populations to fit
Nex_Max_Scale = 1.2 #multipler of how much to let Nex exceed the number predicted exchangable backbone NHs
Pop_Thresh = 0.03 #fall back to n-1 curves if population is below this, does not apply to bootstrap fits, but does exclude from boot average values
Ncurve_p_accept = 0.05 #stringency for accepting more fit populations, higher permits more populations, reasonable values are 0.01 to 0.05 
Random_Seed = 16 #used for parameter initialization

  
Scale_Y_Values = True # if Scale_Y_Values = True, plots will be in original Intensity units
                # fit will always be on normalized Intensity as it is much faster               
Keep_Raw = True # peak_picking will retain the Raw spectrum if True, if False will only keep peaks, auto True for Test_Data
Overlay_replicates = True #add column to figures that is overlay of all available replicates

Use_DiffEvo = False
DiffEvo_threshold = 0 #1e-15
DiffEvo_kwargs = {'polish':True,'maxiter':1000}

########################################
'''end user input''';
########################################
