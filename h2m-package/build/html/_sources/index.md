% H2M documentation master file, created by
% sphinx-quickstart on Thu Feb  8 12:53:22 2024.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# H2M Documentation
```{image} figures/h2m-logo-final.png  
:width: 150px
:align: left
```
H2M is a python package for precise modeling of human vairants in the mouse genome.   

H2M's main functions are:  

1. Querying mouse orthologous genes with human genes input.  

2. Generating homologous mouse mutations with a list of input human mutations. The input format is extremely flexible, allowing for users to input a list of genome coordinates or sequences with desired edits.   

## Package Installation

H2M is available through the python package index. To install, use pip:  

```python
pip install h2m
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