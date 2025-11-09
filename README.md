## Preprint Announcement! (Feb 25, 2025)    
The manuscript of this work is now available on **bioRxiv**!  
Read the preprint: [Computational modeling of human genetic variants in mice](https://www.biorxiv.org/content/10.1101/2025.02.23.639784v1)

# H2M (Human-To-Mouse) 
Author: Kexin Dong  
Date: Jan 4, 2024  

<div align="left">
    <img src="readme_files/final.png"  style="width: 20%;">
</div>

-----

H2M is a python package for the precision modeling of human variants in mice and vice versa.    

H2M's main functions are:  

1. Reading and formatting mutation data from different pulic sources.  

2. Querying orthologous genes between mouse and human.  

3. Generating murine equivalents for human genetic variant input or vice versa. 

Documentation: https://h2m-public.readthedocs.io  

H2M Database Online Portal:  https://human2mouse.com/
## 🔧 Installation  

### 1. Via PyPI  


```python
pip install bioh2m
```

### 2. Download the **.whl** file from [the GitHub repository](https://github.com/kexindon/h2m-public.git) and:  


```python
pip install bioh2m-1.0.0-py3-non-any.whl  
```

ATTENTION: H2M has `pysam` as a dependency. This is for a function that can read .vcf files. If you are experiencing installation problems due to pysam, you can download and install the wheel file in [this GitHub repository](https://github.com/kexindon/h2m-public.git) without this function and the pysam dependency, which has been tested to solve most installation issues. The function rounded off in mini-h2m is also given in the repo. 

H2M has been tested in Python 3.9-3.12.  

## ⏬ Download dependency files  

Reference genome and genome annotation files for human and mouse are needed for this package. To download, visit [the Dropbox Folder](https://www.dropbox.com/scl/fo/1wtrnc9w6s9gemweuw2fv/h?rlkey=hli1z6tv096cjwit5oi6bwggg&dl=0). Some commonly used public data of human mutations are also available here.  

## 🚀 Quick Start     

### Load reference genome  


```python
import bioh2m as h2m
import pandas as pd

path_h_ref, path_m_ref = '.../GCF_000001405.25_GRCh37.p13_genomic.fna.gz', '.../GCF_000001635.27_GRCm39_genomic.fna.gz'
# remember to replace the paths with yours; for human, GRCh38 reference genome assembly is also provided  
records_h, index_list_h = h2m.genome_loader(path_h_ref)
records_m, index_list_m  = h2m.genome_loader(path_m_ref)

path_h_anno, path_m_anno = '.../gencode_v19_GRCh37.db', '.../gencode_vm33_GRCm39.db'
# remember to replace the paths with yours
db_h, db_m = h2m.anno_loader(path_h_anno), h2m.anno_loader(path_m_anno)
```

### 📊 Load human mutation data  


```python
filepath = '/Users/kexindong/Documents/GitHub/Database/PublicDatabase/ClinVar/GRCh37_clinvar_20240206.vcf.gz'
variation_ids = [32798013, 375926, 325626, 140953, 233866, 1796995, 17578, 573320]
df = h2m.clinvar_reader(filepath, variation_ids)
df = h2m.clinvar_to_maf(df)
df = df[['gene_name_h',	'start_h','end_h','ref_seq_h','alt_seq_h','type_h','format','ID']]
df = df.rename(columns={'ID':'index'})
```
<div align="left">
    <img src="readme_files/1.png"  style="width: 45%;">
</div>

### Get canonical transcript IDs for the human genes  


```python
df, df_fail = h2m.get_tx_batch(df, species='h', ver = 37)
```
<div align="left">
    <img src="readme_files/2.png"  style="width: 70%;">
</div>

### Query the gene orthologs in mouse  


```python
df_queried, df_fail = h2m.query_batch(df, direction='h2m')
```
<div align="left">
    <img src="readme_files/3.png"  style="width: 80%;">
</div>

### Get canonical transcript IDs for the murine genes  


```python
df_queried, df_fail = h2m.get_tx_batch(df_queried, species='m')
```
<div align="left">
    <img src="readme_files/4.png"  style="width: 93%;">
</div>
### Compute
the muerine variant equivalents  


```python
df_result, df_fail = h2m.model_batch(df_queried, records_h, index_list_h, records_m, index_list_m, db_h, db_m, 37)
```

![png](readme_files/5.png)

### Data visualization    


```python
df_vs = h2m.model(records_h, index_list_h, records_m, index_list_m, db_h, db_m, 
                  tx_id_h='ENST00000349496.5',tx_id_m='ENSMUST00000007130.15',
                  start=41266098, end=41266098, ref_seq='A', alt_seq='G', ty_h = 'SNP',
                  ver = 37, show_sequence=True)
h2m.visualization(df_vs, flank_size=2, print_size=6)
```
![png](readme_files/ctnnb1.png)

## 🪪 License

`h2m` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

## Cite H2M:  

## Release History  

### 1.0.0, May 3, 2024

First Release.

### 1.0.1, May 3, 2024

Updated the tutorial file.

### 1.0.2, May 7, 2024

Updated the class 5 input for reference and alternate sequences.  


### 1.0.3, May 13, 2024

Manually checked and updated the human gene symbol in the catalog of the human-mouse gene orthologs. 
