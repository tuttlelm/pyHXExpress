## Create python scripts to do a batch run for states in parallel 
## This is setup to do a separate run for each sample type found in the metadf file
## then run "python3 run_batch_<date>.py
## it will create a jumble of text output but will run in parellel speeding things up a bit

import os
import importlib
import sys

hxex_path = os.path.join('/home/tuttle/data/HDX-MS/pyHXExpress') ### <- modify as needed using pyhxexpress.py file
sys.path.append(hxex_path)

import pyhxexpress as hxex      ### <- modify as needed, this is set up for a pyhxexpress.py file
import numpy as np, pandas as pd
import config_fimh as config  ### <- modify as needed

config_file_name = "config_fimh" ### <- set this so it updates value in output script
python_name = "python3" # <- set to e.g. python or python3 as appropriate

def hxex_reload():
    importlib.reload(hxex)
    importlib.reload(config)
    hxex.config = config

hxex_reload()

hxex.config.Output_DIR = os.path.join(hxex.config.Data_DIR,'batchrun_'+str(config.date))
if not os.path.exists(hxex.config.Output_DIR): os.makedirs(hxex.config.Output_DIR)
print("saving output to ",hxex.config.Output_DIR)

metadf = hxex.get_metadf()
samples = metadf['sample'].unique()


### Create scripts for runs of each sample type
scripts = {}
for sample in samples:
    run_dir = os.path.join(hxex.config.Output_DIR,str(sample)+"_run")
    if not os.path.exists(run_dir): os.makedirs(run_dir)
    filename = "run_"+str(sample)+"_"+config.date+".py"
    write_file = os.path.join(run_dir,filename)
    scripts[sample] = write_file

    with open(write_file,"w") as p_file:

        script_text = f'''
# run script for {sample}

import os
import importlib
import sys

hxex_path = os.path.join('{hxex_path}')
sys.path.append(hxex_path)

import pyhxexpress as hxex
import numpy as np, pandas as pd
import {config_file_name} as config

def hxex_reload():
    importlib.reload(hxex)
    importlib.reload(config)
    hxex.config = config

hxex_reload()

hxex.config.Output_DIR = os.path.join("{run_dir}")
if not os.path.exists(hxex.config.Output_DIR): os.makedirs(hxex.config.Output_DIR)
print("saving output to ",hxex.config.Output_DIR)

metadf = hxex.get_metadf()

filtered = hxex.filter_df(metadf,samples="{sample}")

hxex.run_hdx_fits(filtered)
        
        '''

        #print(script_text)
        p_file.write(script_text)



### Create the multiprocessing script that will pool each of the sample scripts
## based on https://medium.com/@rohitobrai11/multithreading-in-python-running-2-scripts-in-parallel-8258c4635182

mp_script_text = f'''
import threading
import subprocess

def run_script(script_name):
    subprocess.run(["{python_name}", script_name])

scripts = {scripts}

if __name__ == "__main__":
    script_thread = {{}}
    for k,v in scripts.items():
        script_thread[k] = threading.Thread(target=run_script, args=[v])
        #print("value",v)

    for k,v in scripts.items():
        script_thread[k].start()

    for k,v in scripts.items():
        script_thread[k].join()

    print("Batch scripts have finished executing.")
'''

#print(mp_script_text)
write_file = os.path.join(hxex.config.Output_DIR,'run_batch_'+config.date+'.py')
with open(write_file,"w") as p_file:
    p_file.write(mp_script_text)
