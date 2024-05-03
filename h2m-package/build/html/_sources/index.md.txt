% H2M documentation master file, created by
% sphinx-quickstart on Thu Feb  8 12:53:22 2024.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# H2M Documentation
```{image} figures/h2m-logo-final.png  
:width: 150px
:align: left
```
H2M is a Python package for high-throughput precision modeling of human variants in the mouse genome and vice cersa.   

H2M's main functions are:  

1. Reading and formatting mutation data from different pulic sources.  

2. Querying orthologous genes between mouse and human.  

3. Generating murine equivalents for human genetic variant input or vice versa. 

See more in the [the GitHub repository](https://github.com/kexindon/h2m-public.git).   

## Package Installation 
H2M is available through the python package index (PyPI). To install, use pip:  
 
```python
    pip install bioh2m
```
```{attention}
Python **3.9-3.12** are recommended since H2M has been tested compatible in them. 
```
```{hint}
H2M has `pysam` as a dependency. This is for a function that can read .vcf files. If you are experiencing installation problems due to pysam, you can download and install the wheel file in [the GitHub repository](https://github.com/kexindon/h2m-public/tree/main/install-wheels) without this function and the pysam dependency, which has been tested to solve most installation issues. The function rounded off in mini-h2m is also given in the repository.  
```

```{toctree}
:caption: 'Contents:'
:maxdepth: 6
quickstart
documentation.md
jupyter
downstream.md
about.md
```