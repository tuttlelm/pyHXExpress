import os
from datetime import datetime
now = datetime.now()
date = now.strftime('%d%b%Y')

#config parameters 08Jan2024
Boot_Seed = True
Bootstrap = True
Data_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data"
Data_Type = 1
Env_limit = 0.95
Env_threshold = 0.1
Full_boot = True
Hide_Figure_Output = False
Keep_Raw = True
Limit_by_envelope = False
Max_Pops = 3
Metadf_File = r"hdxms_testsets_metadf.csv"
Nboot = 20
Ncurve_p_accept = 0.05
Overlay_replicates = True
Output_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data\output_fixedpops_08Jan2024"
Pop_Thresh = 0.03
process_ALL = True
setNoise = 100
Random_Seed = 16
Read_Spectra_List = True
SVG = False
Scale_Y_Values = True
Test_Data = True
User_mutants = ['']
User_peptides = ['']
Y_ERR = 1.0
WRITE_PARAMS = True
