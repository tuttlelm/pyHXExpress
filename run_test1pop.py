import os
import importlib
import pyhxexpress as hxex
import numpy as np, pandas as pd
#pd.set_option('display.max_columns',None)
#%pip install tensorflow==2.13.0
# import tensorflow as tf
# from keras.models import load_model
#import config  
#import hdxms_params_22Dec2023 as config
import test_config as config
#import config_shsp as config

def hxex_reload():
    importlib.reload(hxex)
    importlib.reload(config)
    hxex.config = config

hxex_reload()
#help(hxex)

hxex.config.Save_Spectra = True
hxex.config.Hide_Figure_Output = True
hxex.config.Min_Pops = 1
hxex.config.Max_Pops = 1

hxex.config.Output_DIR = os.path.join(config.Data_DIR,'output_fixed1pop_'+str(config.date),'')
if not os.path.exists(hxex.config.Output_DIR): os.makedirs(hxex.config.Output_DIR)
print("Output will be saved to:\n",hxex.config.Output_DIR,"\n")

metadf = hxex.get_metadf()
hxex.run_hdx_fits(metadf)

