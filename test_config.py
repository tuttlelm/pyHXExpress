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

USE_PARAMS_FILE = False  #### IF THIS IS TRUE ALL PARAMETERS ARE READ FROM PARAMS_FILE:
WRITE_PARAMS = True #save the params to hdxms_params_$.py file in Data_DIR, can then be read in as PARAMS_FILE 
Allow_Overwrite = True #don't create a new filename if file already exists

Read_Spectra_List = True 
Test_Data = True
Data_Type = 1 
process_ALL = True
Data_DIR = 'C:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\pyHXExpress\\Bimodal_HDX_Data'
Read_Spectra_List = True
Metadf_File = "hdxms_testsets_metadf.csv"              
Output_DIR = os.path.join(Data_DIR,'output_'+str(date),'')
if not os.path.exists(Output_DIR): os.makedirs(Output_DIR)
User_mutants = ['']
User_peptides = ['']

Hide_Figure_Output = False #Recommended when processing lots of data. 
SVG = False # also save figures as an svg file, slow, but better for making figures 

Bootstrap = True #False #
Full_boot=True #plot all the bootstrap fits, frac vs nex*mu

Nboot = 20 # number of individual fits to perform, using n_best_curves from initial round of fits
setNoise = 500 #if noise value is known, specify instead of estimating as Y_ERR % of avg Un+TD peaks
Y_ERR = 1.0 #Percent random error applied during boot as y*+np.random.normal(0,yerr), 0.0 for NoNoise, ~0.5% for noise added
            # the absolute Noise value is then Y_ERR * avg(maxInt of Un and TD)
            # this is a very rough way to give a consistent Noise value throughout a dataset. 

Zero_Filling = 6 # large value to try to fit 3 binomials on short peptides, not generally recommended.
Peak_Resolution = 100.0 #ppm, sensitivity of peak picker to expected m/z centers 
Env_threshold = 0.1 #find envelope width at Env_threshold * Intensity_max
Limit_by_envelope = False # only fit up to n = int(z*env/3*Env_limit - 2/3) 
Env_limit = 0.95 #used if Limit_by_envelope = True, rough measure to constrain n_curves fit according to data width & num fit params
Min_Pops = 1 # Min_Pops to Max_Pops sets the range of populations to fit, set to same value to force single population
Max_Pops = 3 # maximum number of underlying populations to fit
Pop_Thresh = 0.03 #fall back to n-1 curves if population is below this, does not apply to bootstrap fits, but does exclude from boot average values
Ncurve_p_accept = 1.0 #stringency for accepting more fit populations, higher permits more populations, reasonable values are 0.01 to 0.05 
Random_Seed = 16 #used for parameter initialization
Boot_Seed = True #if False, same seed as Random_Seed, 
                 #otherwise different seed for each boot iteration (0 to Nboot + Random_Seed + 1 to not repeat initial fit)   
Scale_Y_Values = True # if Scale_Y_Values = True, plots will be in original Intensity units
                # fit will always be on normalized Intensity as it is much faster               
Keep_Raw = True # peak_picking will retain the Raw spectrum if True, if False will only keep peaks, auto True for Test_Data
Overlay_replicates = True #add column to figures that is overlay of all available replicates

########################################
'''end user input''';
########################################
