
# run script for R120G

import os
import importlib
import sys

hxex_path = os.path.join('/home/tuttle/data/HDX-MS/pyHXExpress')
sys.path.append(hxex_path)

import hxex_updating as hxex
import numpy as np, pandas as pd
import config_b5disease as config

def hxex_reload():
    importlib.reload(hxex)
    importlib.reload(config)
    hxex.config = config

hxex_reload()

hxex.config.Output_DIR = os.path.join("/data/tuttle/HDX-MS/Chris_B5_disease/SpecExport/batchrun_14Aug2024/R120G_run")
if not os.path.exists(hxex.config.Output_DIR): os.makedirs(hxex.config.Output_DIR)
print("saving output to ",hxex.config.Output_DIR)

metadf = hxex.get_metadf()

filtered = hxex.filter_df(metadf,samples="R120G")

hxex.run_hdx_fits(filtered)
        
        