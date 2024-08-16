#!pip install PyMuPDF
import pandas as pd
import numpy as np
import fitz
import os
from datetime import datetime

def collate_pdfs(pdf_dir,project=None,samples=None,save_dir=None,addlabel=True):
    '''
    cobbles together a proper order to merge all the generated pdfs after the fact

    #example combining outputs from batch run (keeps samples separate)
    for proj in projects:
        for sample in samples[proj]:
            pdf_dir = os.path.join(batch_dir[proj],sample+'_run')
            #print(pdf_dir)
            collate_pdfs(pdf_dir,project='HSPB5_'+sample+'_'+proj,save_dir=batch_dir[proj])
    '''
    now = datetime.now()
    date = now.strftime("%d%b%Y")

    def categorize(df,types,search_col,new_col):
        for i,type in enumerate(types):
            df[type] = np.where(df[search_col].str.contains(type),str(i)+type,'')
            df[new_col]=''
        # for type in types:
        #     df[new_col] += pdf_info[type]
        df[new_col] = df.apply(lambda x: max(x[types],key=len),axis=1)
        df = df.drop(types,axis=1)
        return df

    def add_header(mfile,header_text):
        for page in mfile:
            prect = page.rect
            header_rect = fitz.Rect(0, 10, prect.width-10,30)  # height 20 points
            page.insert_textbox(header_rect, header_text,
                                fontname="hebo", color= fitz.pdfcolor["blue"],    
                                align=fitz.TEXT_ALIGN_RIGHT)
        return mfile
    
    output_types = ['IndFits','BootFits','ndeut']

    pdfs = [ f for f in os.listdir(pdf_dir) if f[-4:]=='.pdf'  ]
    pdfs = [ p for p in pdfs if p[0:5]=='hxex_']
    pdf_info = pd.DataFrame({'file':pdfs})

    pdf_info = categorize(pdf_info,output_types,'file','type')
    if samples is not None: 
        samples = samples if isinstance(samples,list) else [samples]
        pdf_info = categorize(pdf_info,samples,'file','sample')
    else: pdf_info['sample'] = ''
    if project is not None: proj_str = project+'_'
    else: proj_str = ''
    if save_dir is None: save_dir = pdf_dir

    pdf_info[['start_res','end_res']] = ''
    pdf_info['start_res'] = pdf_info['file'].str[:-4].str.split('-',n=1,).str[0].str[-4:].astype('Int64')
    pdf_info['end_res'] = pdf_info['file'].str[:-4].str.split('-',n=1,).str[1].str[0:4].astype('Int64')
    pdf_info['charge'] = pdf_info['file'].str[:-4].str.split('_').str[-3].str[1:].astype('Int64')
    #pdf_info['data_id'].str[:-4].str.split('_').str[1].astype('Int64')
    pdf_info = pdf_info.sort_values(['start_res','end_res','sample','charge','type']).reset_index(drop=True)
    pdf_info['type'] = pdf_info['type'].str[1:]
    pdf_info['sample'] = pdf_info['sample'].str[1:]
    #display(pdf_info)

    sorted_pdfs = pdf_info['file'].values
    sorted_pdfs = [s for s in sorted_pdfs if '_all_' not in s]
    sorted_pdfs_Ind = [s for s in sorted_pdfs if 'Ind' in s]
    sorted_pdfs_BootFits = [s for s in sorted_pdfs if 'BootFits' in s]
    sorted_pdfs_ndeut = [s for s in sorted_pdfs if 'ndeut' in s]
    #sorted_pdfs

    

    result = fitz.open()
    for pdf in sorted_pdfs:        
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)
            result.insert_pdf(mfile)
    result.save(os.path.join(save_dir,'hxex_all_plots_'+proj_str+date+'.pdf'))

    result = fitz.open()
    for pdf in sorted_pdfs_Ind:
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)
            result.insert_pdf(mfile)            
    result.save(os.path.join(save_dir,'hxex_all_IndFits_'+proj_str+date+'.pdf'))

    result = fitz.open()
    for pdf in sorted_pdfs_BootFits:
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)            
            result.insert_pdf(mfile)
    result.save(os.path.join(save_dir,'hxex_all_BootFits_'+proj_str+date+'.pdf'))

    result = fitz.open()
    for pdf in sorted_pdfs_ndeut:
        with fitz.open(os.path.join(pdf_dir,pdf)) as mfile:
            if addlabel: 
                data_id = pdf.split('_')[1]
                header_text = "data_id: "+data_id
                mfile = add_header(mfile,header_text)            
            result.insert_pdf(mfile)
    result.save(os.path.join(save_dir,'hxex_all_ndeutBoot_'+proj_str+date+'.pdf'))

    return
