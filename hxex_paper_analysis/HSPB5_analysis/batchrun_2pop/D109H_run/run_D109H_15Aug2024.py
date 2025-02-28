
# run script for D109H

import os
import importlib
import sys

hxex_path = os.path.join('/home/tuttle/data/HDX-MS/pyHXExpress')
sys.path.append(hxex_path)

import hxex_updating as hxex
import numpy as np, pandas as pd
import config_b5disease_2pop as config

def hxex_reload():
    importlib.reload(hxex)
    importlib.reload(config)
    hxex.config = config

hxex_reload()

hxex.config.Output_DIR = os.path.join("/data/tuttle/HDX-MS/Chris_B5_disease/SpecExport/batchrun_2pop_15Aug2024/D109H_run")
if not os.path.exists(hxex.config.Output_DIR): os.makedirs(hxex.config.Output_DIR)
print("saving output to ",hxex.config.Output_DIR)

metadf = hxex.get_metadf()

filtered = hxex.filter_df(metadf,samples="D109H")

hxex.run_hdx_fits(filtered)
        
        