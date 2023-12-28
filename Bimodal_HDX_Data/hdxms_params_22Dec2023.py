import os
from datetime import datetime
now = datetime.now()
date = now.strftime('%d%b%Y')
#config parameters 22Dec2023Boot_Seed = True
Bootstrap = True
Data_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data"
Data_Type = 1
Env_limit = 1.0
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
Output_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data\output_22Dec2023"
Pop_Thresh = 0.03
process_ALL = False
setNoise = 200
Random_Seed = 16
# missing parameter Read_meta
SVG = False
Scale_Y_Values = True
Test_Data = True
User_mutants = ['HSPB1only']
User_peptides = ['0001-0011']
Y_ERR = 1.0
WRITE_PARAMS = True
