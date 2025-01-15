import os
from datetime import datetime
now = datetime.now()
date = now.strftime("%d%b%Y")
USE_PARAMS_FILE = False


##########################################
'''begin user input'''
##########################################

WRITE_PARAMS = True #save the params to hdxms_params_$.py file in Data_DIR, can then be read in as PARAMS_FILE 
Allow_Overwrite = True #don't create a new filename if file already exists
Save_Spectra = False
Read_Spectra_List = False # Specify files to be run in a file, includes peptide/charge info. See example files.
                # To use this, Recommend setting to False to create and write 'metadf' to file with all availble datasets
                # then remove unwanted datasets from the file, and read it in with Read_Spectra_List = True
                # and Metadf_File set to the appropriate filename 
Test_Data = False
Data_Type = 2  #Test_Data = True has precedence over this, will make Data_Type = 1
    #1: 'xlxs' , each file contains all timepoint reps at Data_DIR/onesample_onepeptide_alltimepts_allreps_onecharge.xlsx
                # current recognized format is e.g. HSPB1_B1B5_0001-0011-MTERRVPFSLL-z2-allspectra.xlsx
                # <sample name>_<peptide start>-<peptide end>-<peptide>-<zcharge>-<unused label>.xlsx
                # allows for replicats of UnDeut/TD even though HX-Express xlsm does not
    #2: 'SpecExport', as exported from HDExaminer Data_DIR/samples/peptides/onetimept_onerep_onecharge.csv
                # this mode requires a sample.fasta file in the Data_DIR for each sample to be processed, with matching names
if Data_Type == 1:
    Data_DIR = os.getcwd()
    #Data_DIR = 'c:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\ns_HSPB1_Bimodal_Peptide_Data'
    #Data_DIR = 'C:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\pyHXExpress\\Bimodal_HDX_Data'
    Metadf_File = "hdx_spectra_list_metadf_02Nov2023.csv" #only used if Read_Spectra_List = True; designates files to process
    process_ALL = False # if True will assume all .xlsx files are HDX data, use with care
    User_mutants = ['HSPB1only',]#'HSPB1_B1B6'] #['all'] #first element can be 'all' to include all mutants and/or peptides in directory
    User_peptides = ['0001-0011',]#'0078-0094']
if Data_Type == 2:
    Data_DIR = os.getcwd()
    #Data_DIR = '/data/tuttle/HDX-MS/Pearl_SpecExport_30oct2023/SpecExport'
    #Data_DIR = '/data/tuttle/HDX-MS/Pearl_FimHWTL34K_V6/SpecExport/'
    #Data_DIR = 'c:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\ns_HSPB1_Bimodal_Peptide_Data\\SpecExport'
    Metadf_File = "hdx_spectra_list_metadf_02Nov2023.csv" #only used if Read_Spectra_List = True; designates files to process
    process_ALL = True #process_all = True is limited to existing .fasta files, this setting overrides user_ settings
    User_mutants = [''] #['WT','S19D','S45D','S59D','D3']#['All'] #
    User_peptides =  ['']#'0049-0054']#['0034-0045'] #['0093-0116'] #['0090-0113']'0122-0166']#

if Test_Data: 
    Data_Type = 1
    Data_DIR = os.getcwd()
    #Data_DIR = 'c:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\HX-Express3'
    #Data_DIR = 'C:\\Users\\tuttl\\OneDrive\\Documents\\My Documents\\KlevitHahn\\hdx-ms\\pyHXExpress\\Bimodal_HDX_Data'
    #Test_Sets = ['v3_Angiotensin_Bimodals.xlsx','v3_GluFib_Bimodals.xlsx']
    Read_Spectra_List = True
    Metadf_File = "hdxms_testsets_metadf.csv"
                
Output_DIR = os.path.join(Data_DIR,'hdxms_analysis_'+str(date),'')
#if not os.path.exists(Output_DIR): os.makedirs(Output_DIR)

Preset_Pops = False #Use predetermined number of populations to fit curves, overrides Min/Max_Pops if given
Preset_Pops_File = os.path.join(Output_DIR,"all_data_results.csv")

Generate_Plots = True # generate all plots and save as pdf files, can slow things down
Hide_Figure_Output = True #Recommended when processing lots of data. 
SVG = False # also save figures as an svg file, slow, but better for making figures 

BestFit_of_X = 3 # try this number of initial fits and use best fit to evaluate best_n_curves before nboot bootstrap fits
Residual_Cutoff = 1e-5 # if the residual sum squared is less than this value, use that number of curves, set to 0 to disable
Use_DiffEvo = False # use scipy.optimize.differential_evolution instead of curve_fit 
DiffEvo_threshold = 1e-11 # If using differential_evolution to fit data, stop when threshold satisfied, 0 to disable
Bootstrap = True #False #
Full_boot=True #plot all the bootstrap fits, frac vs nex*mu

FullDeut_Time = 1e6 # assumed time in seconds for the Total Deuteration sample for exchange plots
Dfrac = 0.9 #fraction deuterated for "fully" deuterated exchange buffer 
Nterm_subtract = 1 #number of N-term residues to remove from possible exchanging NH's (usually 1 or 2)
Nboot = 20 # number of individual fits to perform, using n_best_curves from initial round of fits
setNoise = 200 #if noise value is known, specify instead of estimating as Y_ERR % of avg Un+TD peaks
Y_ERR = 1.0 #Percent random error applied during boot as y*+np.random.normal(0,yerr), 0.0 for NoNoise, ~0.5% for noise added
            # the absolute Noise value is then Y_ERR * avg(maxInt of Un and TD)
            # this is a very rough way to give a consistent Noise value throughout a dataset. 
Binomial_dCorr = True #use the n=1 binomial fit to determine UN/TD centers for d_corr calculation
                    # this is essential if dirty spectra; otherwise using sum(mz*y)/sum(y) which will be wrong if there are contaminating peaks
Peak_Resolution = 50.0 #ppm, sensitivity of peak picker to expected m/z centers 
Zero_Filling = 4 #number of zero points below noise threshold to include in peak picking
Nex_Max_Scale = 1.2
Env_threshold = 0.1 #find envelope width at Env_threshold * Intensity_max
Limit_by_envelope = False # only fit up to n = int(z*env/3*Env_limit - 2/3) 
Env_limit = 1.0 #used if Limit_by_envelope = True, rough measure to constrain n_curves fit according to data width & num fit params
Min_Pops = 1 # Min_Pops to Max_Pops sets the range of populations to fit, set to same value to force single population
Max_Pops = 3 #maximum number of underlying populations to fit
Pop_Thresh = 0.03 #fall back to n-1 curves if population is below this, does not apply to bootstrap fits, but does exclude from boot average values
Ncurve_p_accept = 0.05 #stringency for accepting more fit populations      
Random_Seed = 16 #used for parameter initialization
Boot_Seed = True #if False, same seed as Random_Seed, 
                 #otherwise different seed for each boot iteration (0 to Nboot + Random_Seed + 1 to not repeat initial fit)   
Scale_Y_Values = True # if Scale_Y_Values = True, plots will be in original Intensity units
                # fit will always be on normalized Intensity as it is much faster               
Keep_Raw = True # peak_picking will retain the Raw spectrum if True, if False will only keep peaks, auto True for Test_Data
Overlay_replicates = True #add column to figures that is overlay of all available replicates

DiffEvo_kwargs = {'polish':True,'maxiter':100}

########################################
'''end user input''';
########################################
