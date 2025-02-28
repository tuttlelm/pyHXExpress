import os
from datetime import datetime
now = datetime.now()
date = now.strftime('%d%b%Y')

#config parameters 14Aug2024
Allow_Overwrite = True
BestFit_of_X = 5
Binomial_dCorr = True
Data_DIR = r"/data/tuttle/HDX-MS/Chris_B5_disease/SpecExport"
Data_Type = 2
DiffEvo_kwargs = {'polish': True, 'maxiter': 1000}
DiffEvo_threshold = 0
Dfrac = 0.9
Env_limit = 1.0
Env_threshold = 0.1
FullDeut_Time = 1000000.0
Hide_Figure_Output = True
Keep_Raw = True
Limit_by_envelope = False
Max_Pops = 2
Metadf_File = r"."
Min_Pops = 1
Nboot = 20
Ncurve_p_accept = 0.05
Nex_Max_Scale = 1.2
Nterm_subtract = 2
Output_DIR = r"/data/tuttle/HDX-MS/Chris_B5_disease/SpecExport/batchrun_14Aug2024/R120G_run"
Overlay_replicates = True
Peak_Resolution = 50.0
Pop_Thresh = 0.05
Preset_Pops = False
Preset_Pops_File = r"."
process_ALL = True
Random_Seed = 16
Read_Spectra_List = False
Residual_Cutoff = 1e-05
Save_Spectra = True
Scale_Y_Values = True
setNoise = 50
SVG = False
Test_Data = False
Use_DiffEvo = False
User_mutants = ['']
User_peptides = ['']
WRITE_PARAMS = True
Y_ERR = 1.0
Zero_Filling = 3
