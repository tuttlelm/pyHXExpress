import os
from datetime import datetime
now = datetime.now()
date = now.strftime('%d%b%Y')

#config parameters 02Mar2024
Allow_Overwrite = True
Binomial_dCorr = True
Boot_Seed = True
Bootstrap = True
Data_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data"
Data_Type = 1
Dfrac = 0.9
Env_limit = 0.95
Env_threshold = 0.1
Full_boot = True
Hide_Figure_Output = False
Keep_Raw = True
Limit_by_envelope = False
Max_Pops = 3
Metadf_File = r"hdxms_testsets_metadf.csv"
Min_Pops = 1
Nboot = 20
Ncurve_p_accept = 0.05
# missing parameter Nex_Max_ScaleNterm_subtract
Output_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data\output_02Mar2024"
Overlay_replicates = True
Peak_Resolution = 70.0
Pop_Thresh = 0.05
Preset_Pops = False
Preset_Pops_File = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data\test_configdf_26feb2024.csv"
process_ALL = True
Random_Seed = 16
Read_Spectra_List = True
Save_Spectra = False
Scale_Y_Values = True
setNoise = 100
SVG = False
Test_Data = True
User_mutants = ['']
User_peptides = ['']
WRITE_PARAMS = True
Y_ERR = 1.0
Zero_Filling = 6
