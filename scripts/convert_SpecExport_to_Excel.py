''' 
Basic framework to convert read in SpecExport rawdata to Excel format
currently (pyHXExpress v0.0.100) the function export_to_hxexpress is expecting its arguments
to contain only the desired output data. Here a basic loop over entries in a filtered dataframe
takes care of this so a single excel file gets created for every data_id (metadf index).

To use this script, the user should update the Data_DIR to point to the desired SpecExport folder
and the filtered = hxex.filter_df() function should filter metadf to the desired subset of datasets
If you want to convert all of the data in SpecExport, you could just change this to filtered = metadf
'''

import os
import pyhxexpress as hxex
import numpy as np, pandas as pd

hxex.config.Data_Type = 2 #for SpecExport type data 

#Currently output is saved to Data_DIR, will update this in future versions
hxex.config.Data_DIR = '/home/tuttle/data/HDX-MS/sHSP_Heterooligomers/b5_hetero/SpecExport'
print("Output will be saved to:\n",hxex.config.Data_DIR,"\n")
metadf = hxex.get_metadf()

#filtered should contain the datasets to be converted
filtered = hxex.filter_df(metadf,index=[0,1],quiet=False)

deut,raw = hxex.get_data(filtered)
for id in filtered.index:

    print("converting data_id", id,filtered.loc[[id]]['sample'].values[0],
          filtered.loc[[id]]['file'].values[0],"z"+str(int(filtered.loc[[id]]['charge'].values[0])),
          "to HX-Express Excel format")
    hxcols = hxex.export_to_hxexpress(hxex.filter_df(raw,data_ids=id),filtered.loc[[id]],save_xls=True,removeUNTDreps=True)
