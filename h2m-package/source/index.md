% H2M documentation master file, created by
% sphinx-quickstart on Thu Feb  8 12:53:22 2024.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# H2M Documentation
```{image} figures/h2m-logo-final.png  
:width: 150px
:align: left
```
H2M is a Python package for high-throughput precision modeling of human vairants in the mouse genome and vice cersa.   

H2M's main functions are:  

1. Reading and formatting mutation data from different pulic sources.  

2. Querying orthologous genes between mouse and human.  

3. Generating murine equivalents for human genetic variant input or vice versa. 

## Package Installation 
H2M is available through the python package index (PyPI). To install, use pip:  
 
```python
    pip install h2m
```
```{attention}
Python **3.9-3.12** are recommended since H2M has been tested compatible in them. 
```
```{hint}
Function vcf_reader() for reading .vcf files in H2M
```

```{toctree}
:caption: 'Contents:'
:maxdepth: 6
quickstart
jupyter
documentation.md
downstream.md
about.md
```