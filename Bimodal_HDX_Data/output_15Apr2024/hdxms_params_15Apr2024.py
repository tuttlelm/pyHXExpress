import os
from datetime import datetime
now = datetime.now()
date = now.strftime('%d%b%Y')

#config parameters 15Apr2024
Allow_Overwrite = True
BestFit_of_X = 3
Binomial_dCorr = True
Data_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data"
Data_Type = 1
DiffEvo_kwargs = {'polish': True, 'maxiter': 1000}
DiffEvo_threshold = 0
Dfrac = 0.9
Env_limit = 0.95
Env_threshold = 0.1
FullDeut_Time = 1000000.0
Hide_Figure_Output = False
Keep_Raw = True
Limit_by_envelope = False
Max_Pops = 3
Metadf_File = r"hdxms_testsets_metadf.csv"
Min_Pops = 1
Nboot = 20
Ncurve_p_accept = 0.05
Nex_Max_Scale = 1.2
Nterm_subtract = 1
Output_DIR = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data\output_15Apr2024"
Overlay_replicates = True
Peak_Resolution = 70.0
Pop_Thresh = 0.05
Preset_Pops = False
Preset_Pops_File = r"C:\Users\tuttl\OneDrive\Documents\My Documents\KlevitHahn\hdx-ms\pyHXExpress\Bimodal_HDX_Data\test_configdf_26feb2024.csv"
process_ALL = True
Random_Seed = 16
Read_Spectra_List = True
Residual_Cutoff = 3.3e-05
Save_Spectra = False
Scale_Y_Values = True
setNoise = 5000
SVG = False
Test_Data = True
Use_DiffEvo = False
User_mutants = ['']
User_peptides = ['']
WRITE_PARAMS = True
Y_ERR = 2.0
Zero_Filling = 3
