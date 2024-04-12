from setuptools import setup, find_packages

VERSION = '0.0.100'
DESCRIPTION = 'pyHXExpress'
LONG_DESCRIPTION = 'High-throughput polymodal analysis of protein HDX-MS spectra'

setup(
    name="pyhxexpress",
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    author="<Lisa M Tuttle>",
    author_email="<tuttlelm@uw.edu>",
    url='https://github.com/tuttlelm/pyHXExpress',
    license='GPLv3',
    packages=find_packages(),
    install_requires=['numpy','pandas','pymupdf','scipy','matplotlib','pyteomics',
              'brain-isotopic-distribution','biopython'],
    keywords='HDX-MS',

    py_modules = ['pyhxexpress'],
    
  
    classifiers= [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        "Programming Language :: Python :: 3",
    ]
)
