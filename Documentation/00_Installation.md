
### Installation

You may consider creating a virtual environment before install. This can be helfpul if you have other packages with specific version requirements. See, for example, https://python.land/virtual-environments/virtualenv or https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

There are a few ways you can use pyHXExpress. This is not exhaustive but should generally work. If you encounter any issues, please create an issue or shoot me an email. 

#### 1. pip install

    
    pip install pyhxexpress
or

    pip3 install pyhxexpress

Using pip or pip3 depends on your python installation. This command should install pyhxexpress and all of the required modules. Depending on your particular python environment, some of these dependency installs may complain. You can try downloading the requirements.txt file and run 'pip install -r requirements.txt' and then repeat the 'pip install pyhxexpress' command. 

If you take this approach, you will load pyhxexpress in your python command line or Jupyter Notebook session by using 

    import pyhxexpress.hxex as hxex

#### 2. Just grab the files

It is sufficient to have they hxex.py and config.py files in the same folder as where you will run scripts or create a Jupyter Notebook. On github these are in the pyhxexpress subfolder. If you go to each file, there will be a download raw option. 

See the requirements.txt file for additional modules that may need to be installed (pip install -r requirements.txt). Generally, if you get going with this and see an error expressing some sense of 'module not found', you can just try pip install module to see if that sorts things out. 

If the files are in some different folder, you can just add the path to the files in your python script or Jupyter Notebook session. The particular syntax for your path will depend on your OS, example here is for Windows.


    import sys
    sys.path.append('c:\\Users\\user\\path\\to\\pyHXExpress_files\\')


If you take this approach, you will load pyhxexpress in your python command line or Jupyter Notebook session by using 

    import hxex

(note that if you have retained the \_\_init\_\_.py file in the directory you should use the same import command as for method (1): 'import pyhxexpress.hxex as hxex'). This assumes the directory name is pyhxexpress. 


###### Last edited by LMT on 14April2024