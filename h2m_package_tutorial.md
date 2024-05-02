# H2M Tutorial
Author: Kexin Dong  
Published: Jan 4, 2024    
Update: May 2, 2024  

H2M is a python package for the precision modeling of human vairants in mice and vice versa.  

H2M's main functions are:  

1. Reading and formatting mutation data from different pulic sources.  

2. Querying orthologous genes between mouse and human.  

3. Generating murine equivalents for human genetic variant input or vice versa.   

See more in the [the GitHub repository](https://github.com/kexindon/h2m-public.git).   

# Quick Start  

## Installation  


```python
pip install h2m  
```

### 2. Download the **.whl** file from [the GitHub repository](https://github.com/kexindon/h2m-public.git) and:  


```python
pip install h2m-1.0.0-py3-non-any.whl  
```

ATTENTION: H2M has `pysam` as a dependency. This is for a function that can read .vcf files. If you are experiencing installation problems due to pysam, you can download and install the wheel file named as mini-h2m in [the GitHub repository](https://github.com/kexindon/h2m-public.git) without this function and the pysam dependency, which has been tested to solve most installation issues. The function rounded off in mini-h2m is also given in the repo. 

### H2M has been tested in Python 3.9-3.12.  

## Importing packages


```python
# pip install h2m
import h2m
import pandas as pd
```

## Loading data  
We should upload reference genome and GENCODE annotation data for both human and mouse, which could be directly downloaded from a public [dropbox](https://www.dropbox.com/scl/fo/1wtrnc9w6s9gemweuw2fv/h?rlkey=hli1z6tv096cjwit5oi6bwggg&dl=0).  
Both GRCh37 and GRCh38 human reference genome assemblys are available. Upload the one that you are going to use.  


```python
path_h_ref, path_m_ref = '.../GCF_000001405.25_GRCh37.p13_genomic.fna.gz', '.../GCF_000001635.27_GRCm39_genomic.fna.gz'
# remember to replace the paths with yours; for human, GRCh38 reference genome assembly is also provided  
records_h, index_list_h = h2m.genome_loader(path_h_ref)
records_m, index_list_m  = h2m.genome_loader(path_m_ref)

path_h_anno, path_m_anno = '.../gencode_v19_GRCh37.db', '.../gencode_vm33_GRCm39.db'
# remember to replace the paths with yours
db_h, db_m = h2m.anno_loader(path_h_anno), h2m.anno_loader(path_m_anno)
```

It may take up to 3 minutes to load the ref genomes.

A notebook file about how to generate the db file: [1_prepare_gencode_annotation_file.ipynb]('1_prepare_gencode_annotation_file.ipynb')

## Batch Processing

### Input format  

Common mutation data formats include MAF (Mutation Annotation Format, used by cBioPortal), VCF (Variant Call Format, used by genomAD), and ClinVar (a modified VCF format, used by ClinVar). Mutation coordinates, reference and alternative sequences are recorded in slightly different ways between the three. 

![png](readme_files/format.png)

In batch processing, H2M accepts MAF (Mutation Annotation Format) input. More information about MAF format can be found at [GDC Documentation](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#:~:text=Mutation). For MAF files, you need to build a pandas dataframe with columns as the following example: 

![png](readme_files/1.png)  

For VCF and ClinVar files, you will need to convert the mutation coordinates and sequence information to MAF format after this. This can be achieved simply by using H2M built-in functions. Data used in this tutorial can be downloaded from [dropbox](https://www.dropbox.com/scl/fo/1wtrnc9w6s9gemweuw2fv/h?rlkey=hli1z6tv096cjwit5oi6bwggg&dl=0).  

#### Read from cBioPortal - MAF  
This format is compatible with all of the datasets in the cBioPortal, as well as TCGA and AACR-GENIE. Download the txt mutation data file from such public dataset and then load it as follows:  


```python
path_aacr = '/Users/kexindong/Documents/GitHub/Database/PublicDatabase/AACR-GENIE/v15.0/data_mutations_extended.txt'
df = h2m.cbio_reader(path_aacr)
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>index</th>
      <th>gene_name_h</th>
      <th>Entrez_Gene_Id</th>
      <th>Center</th>
      <th>NCBI_Build</th>
      <th>Chromosome</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>Strand</th>
      <th>Consequence</th>
      <th>...</th>
      <th>Polyphen_Prediction</th>
      <th>Polyphen_Score</th>
      <th>SIFT_Prediction</th>
      <th>SIFT_Score</th>
      <th>SWISSPROT</th>
      <th>n_depth</th>
      <th>t_depth</th>
      <th>Annotation_Status</th>
      <th>mutationInCis_Flag</th>
      <th>format</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>KRAS</td>
      <td>3845.0</td>
      <td>JHU</td>
      <td>GRCh37</td>
      <td>12</td>
      <td>25398285</td>
      <td>25398285</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>probably_damaging</td>
      <td>0.991</td>
      <td>deleterious</td>
      <td>0.04</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1623.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>BRAF</td>
      <td>673.0</td>
      <td>JHU</td>
      <td>GRCh37</td>
      <td>7</td>
      <td>140453136</td>
      <td>140453136</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>probably_damaging</td>
      <td>0.963</td>
      <td>deleterious</td>
      <td>0.00</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1031.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>EGFR</td>
      <td>1956.0</td>
      <td>JHU</td>
      <td>GRCh37</td>
      <td>7</td>
      <td>55249071</td>
      <td>55249071</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>probably_damaging</td>
      <td>1.000</td>
      <td>deleterious</td>
      <td>0.00</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>692.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>TP53</td>
      <td>7157.0</td>
      <td>JHU</td>
      <td>GRCh37</td>
      <td>17</td>
      <td>7577120</td>
      <td>7577120</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>possibly_damaging</td>
      <td>0.643</td>
      <td>tolerated</td>
      <td>0.13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>930.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>NRAS</td>
      <td>4893.0</td>
      <td>JHU</td>
      <td>GRCh37</td>
      <td>1</td>
      <td>115256529</td>
      <td>115256529</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>benign</td>
      <td>0.251</td>
      <td>tolerated</td>
      <td>0.06</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2277.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1806756</th>
      <td>1806756</td>
      <td>PHF6</td>
      <td>84295.0</td>
      <td>PROV</td>
      <td>GRCh37</td>
      <td>X</td>
      <td>133551330</td>
      <td>133551330</td>
      <td>+</td>
      <td>splice_region_variant,synonymous_variant</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>311.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>1806757</th>
      <td>1806757</td>
      <td>PHF6</td>
      <td>84295.0</td>
      <td>PROV</td>
      <td>GRCh37</td>
      <td>X</td>
      <td>133559236</td>
      <td>133559236</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>probably_damaging</td>
      <td>0.988</td>
      <td>deleterious</td>
      <td>0.00</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>381.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>1806758</th>
      <td>1806758</td>
      <td>PHF6</td>
      <td>84295.0</td>
      <td>PROV</td>
      <td>GRCh37</td>
      <td>X</td>
      <td>133559255</td>
      <td>133559256</td>
      <td>+</td>
      <td>frameshift_variant</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>509.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>1806759</th>
      <td>1806759</td>
      <td>PHF6</td>
      <td>84295.0</td>
      <td>PROV</td>
      <td>GRCh37</td>
      <td>X</td>
      <td>133559302</td>
      <td>133559302</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>benign</td>
      <td>0.066</td>
      <td>deleterious_low_confidence</td>
      <td>0.01</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>788.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>1806760</th>
      <td>1806760</td>
      <td>PHF6</td>
      <td>84295.0</td>
      <td>PROV</td>
      <td>GRCh37</td>
      <td>X</td>
      <td>133559304</td>
      <td>133559304</td>
      <td>+</td>
      <td>missense_variant</td>
      <td>...</td>
      <td>benign</td>
      <td>0.080</td>
      <td>tolerated_low_confidence</td>
      <td>0.24</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>624.0</td>
      <td>SUCCESS</td>
      <td>False</td>
      <td>MAF</td>
    </tr>
  </tbody>
</table>
<p>1806761 rows × 66 columns</p>
</div>



Set `keep = False` to keep H2M-required columns only.  


```python
df = h2m.cbio_reader(path_aacr,keep=False)
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>index</th>
      <th>gene_name_h</th>
      <th>tx_id_h</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>type_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>format</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>KRAS</td>
      <td>ENST00000256078.4</td>
      <td>25398285</td>
      <td>25398285</td>
      <td>SNP</td>
      <td>C</td>
      <td>A</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>BRAF</td>
      <td>ENST00000288602.6</td>
      <td>140453136</td>
      <td>140453136</td>
      <td>SNP</td>
      <td>A</td>
      <td>T</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>EGFR</td>
      <td>ENST00000275493.2</td>
      <td>55249071</td>
      <td>55249071</td>
      <td>SNP</td>
      <td>C</td>
      <td>T</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>TP53</td>
      <td>ENST00000269305.4</td>
      <td>7577120</td>
      <td>7577120</td>
      <td>SNP</td>
      <td>C</td>
      <td>T</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>NRAS</td>
      <td>ENST00000369535.4</td>
      <td>115256529</td>
      <td>115256529</td>
      <td>SNP</td>
      <td>T</td>
      <td>C</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>851083</th>
      <td>851083</td>
      <td>PHF6</td>
      <td>ENST00000332070.3</td>
      <td>133551330</td>
      <td>133551330</td>
      <td>SNP</td>
      <td>C</td>
      <td>T</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>851084</th>
      <td>851084</td>
      <td>PHF6</td>
      <td>ENST00000332070.3</td>
      <td>133559236</td>
      <td>133559236</td>
      <td>SNP</td>
      <td>A</td>
      <td>G</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>851085</th>
      <td>851085</td>
      <td>PHF6</td>
      <td>ENST00000332070.3</td>
      <td>133559255</td>
      <td>133559256</td>
      <td>INS</td>
      <td>-</td>
      <td>T</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>851086</th>
      <td>851086</td>
      <td>PHF6</td>
      <td>ENST00000332070.3</td>
      <td>133559302</td>
      <td>133559302</td>
      <td>SNP</td>
      <td>G</td>
      <td>T</td>
      <td>MAF</td>
    </tr>
    <tr>
      <th>851087</th>
      <td>851087</td>
      <td>PHF6</td>
      <td>ENST00000332070.3</td>
      <td>133559304</td>
      <td>133559304</td>
      <td>SNP</td>
      <td>G</td>
      <td>A</td>
      <td>MAF</td>
    </tr>
  </tbody>
</table>
<p>851088 rows × 9 columns</p>
</div>



#### Read from GenomAD  - VCF 
Search a specific gene in genomAD browser, and download the conlusion csv.  


```python
# downloaded TP53 variants from genomAD
df = h2m.vcf_reader('/Users/kexindong/Documents/GitHub/Database/PublicDatabase/genomAD/gnomAD_v4.1.0_ENSG00000141510_2024_04_29_12_41_42.csv',keep=False)
df['gene_name_h'] = 'TP53'
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>start_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>ID</th>
      <th>index</th>
      <th>format</th>
      <th>gene_name_h</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>7661882</td>
      <td>C</td>
      <td>T</td>
      <td>17-7661882-C-T</td>
      <td>0</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>1</th>
      <td>7661888</td>
      <td>T</td>
      <td>A</td>
      <td>17-7661888-T-A</td>
      <td>1</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2</th>
      <td>7661899</td>
      <td>A</td>
      <td>G</td>
      <td>17-7661899-A-G</td>
      <td>2</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>3</th>
      <td>7661904</td>
      <td>C</td>
      <td>T</td>
      <td>17-7661904-C-T</td>
      <td>3</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>4</th>
      <td>7661905</td>
      <td>G</td>
      <td>A</td>
      <td>17-7661905-G-A</td>
      <td>4</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2036</th>
      <td>7676665</td>
      <td>C</td>
      <td>T</td>
      <td>17-7676665-C-T</td>
      <td>2036</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2037</th>
      <td>7676667</td>
      <td>C</td>
      <td>A</td>
      <td>17-7676667-C-A</td>
      <td>2037</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2038</th>
      <td>7676667</td>
      <td>C</td>
      <td>T</td>
      <td>17-7676667-C-T</td>
      <td>2038</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2039</th>
      <td>7676668</td>
      <td>T</td>
      <td>G</td>
      <td>17-7676668-T-G</td>
      <td>2039</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2040</th>
      <td>7676669</td>
      <td>G</td>
      <td>T</td>
      <td>17-7676669-G-T</td>
      <td>2040</td>
      <td>VCF</td>
      <td>TP53</td>
    </tr>
  </tbody>
</table>
<p>2041 rows × 7 columns</p>
</div>




```python
df = h2m.vcf_to_maf(df)
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ID</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>type_h</th>
      <th>index</th>
      <th>format</th>
      <th>gene_name_h</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>17-7661882-C-T</td>
      <td>7661882</td>
      <td>7661882</td>
      <td>C</td>
      <td>T</td>
      <td>SNP</td>
      <td>0</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>1</th>
      <td>17-7661888-T-A</td>
      <td>7661888</td>
      <td>7661888</td>
      <td>T</td>
      <td>A</td>
      <td>SNP</td>
      <td>1</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2</th>
      <td>17-7661899-A-G</td>
      <td>7661899</td>
      <td>7661899</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>2</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>3</th>
      <td>17-7661904-C-T</td>
      <td>7661904</td>
      <td>7661904</td>
      <td>C</td>
      <td>T</td>
      <td>SNP</td>
      <td>3</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>4</th>
      <td>17-7661905-G-A</td>
      <td>7661905</td>
      <td>7661905</td>
      <td>G</td>
      <td>A</td>
      <td>SNP</td>
      <td>4</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2036</th>
      <td>17-7676665-C-T</td>
      <td>7676665</td>
      <td>7676665</td>
      <td>C</td>
      <td>T</td>
      <td>SNP</td>
      <td>2036</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2037</th>
      <td>17-7676667-C-A</td>
      <td>7676667</td>
      <td>7676667</td>
      <td>C</td>
      <td>A</td>
      <td>SNP</td>
      <td>2037</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2038</th>
      <td>17-7676667-C-T</td>
      <td>7676667</td>
      <td>7676667</td>
      <td>C</td>
      <td>T</td>
      <td>SNP</td>
      <td>2038</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2039</th>
      <td>17-7676668-T-G</td>
      <td>7676668</td>
      <td>7676668</td>
      <td>T</td>
      <td>G</td>
      <td>SNP</td>
      <td>2039</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
    <tr>
      <th>2040</th>
      <td>17-7676669-G-T</td>
      <td>7676669</td>
      <td>7676669</td>
      <td>G</td>
      <td>T</td>
      <td>SNP</td>
      <td>2040</td>
      <td>MAF</td>
      <td>TP53</td>
    </tr>
  </tbody>
</table>
<p>2041 rows × 9 columns</p>
</div>



#### Read from ClinVar
download a ClinVar vcf.gz file, and choose your desired Variation IDs that you wish to model. This can be accessed at: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/   


```python
filepath = '/Users/kexindong/Documents/GitHub/Database/PublicDatabase/ClinVar/GRCh37_clinvar_20240206.vcf.gz'
variation_ids = [32798013, 375926, 325626, 140953, 233866, 1796995, 17578, 573320]
df = h2m.clinvar_reader(filepath, variation_ids)
df = h2m.clinvar_to_maf(df)
df = df[['gene_name_h',	'start_h','end_h','ref_seq_h','alt_seq_h','type_h','format','ID']]
df = df.rename(columns={'ID':'index'})
```

    [E::idx_find_and_load] Could not retrieve index file for '/Users/kexindong/Documents/GitHub/Database/PublicDatabase/ClinVar/GRCh37_clinvar_20240206.vcf.gz'
    [W::vcf_parse] Contig '1' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '2' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '3' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '4' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '5' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '6' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '7' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '8' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '9' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '10' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '11' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '12' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '13' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '14' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '15' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '16' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '17' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '18' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '19' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '20' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '21' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig '22' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'X' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'Y' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'MT' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'NT_113889.1' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'NT_167222.1' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'NW_003315925.1' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'NW_003315947.1' is not defined in the header. (Quick workaround: index the file with tabix.)
    [W::vcf_parse] Contig 'NW_003315950.2' is not defined in the header. (Quick workaround: index the file with tabix.)



```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>type_h</th>
      <th>format</th>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CTNNB1</td>
      <td>41266098</td>
      <td>41266098</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>17578</td>
    </tr>
    <tr>
      <th>1</th>
      <td>KIT</td>
      <td>55599320</td>
      <td>55599321</td>
      <td>GA</td>
      <td>AT</td>
      <td>DNP</td>
      <td>MAF</td>
      <td>375926</td>
    </tr>
    <tr>
      <th>2</th>
      <td>APC</td>
      <td>112174150</td>
      <td>112174150</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>233866</td>
    </tr>
    <tr>
      <th>3</th>
      <td>APC</td>
      <td>112174154</td>
      <td>112174155</td>
      <td>-</td>
      <td>A</td>
      <td>INS</td>
      <td>MAF</td>
      <td>1796995</td>
    </tr>
    <tr>
      <th>4</th>
      <td>MRE11</td>
      <td>94200987</td>
      <td>94200987</td>
      <td>G</td>
      <td>A</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>140953</td>
    </tr>
    <tr>
      <th>5</th>
      <td>TP53</td>
      <td>7572148</td>
      <td>7572148</td>
      <td>G</td>
      <td>-</td>
      <td>DEL</td>
      <td>MAF</td>
      <td>325626</td>
    </tr>
    <tr>
      <th>6</th>
      <td>BRCA1</td>
      <td>41244325</td>
      <td>41244325</td>
      <td>T</td>
      <td>C</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>573320</td>
    </tr>
  </tbody>
</table>
</div>



### Get canonical transcript ID 

There will be returning two dataframes for success and failures.


```python
df, df_fail = h2m.get_tx_batch(df, species='h', ver = 37)
```

    No error occurs.



```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>type_h</th>
      <th>format</th>
      <th>index</th>
      <th>tx_id_h</th>
      <th>ref_genome_h</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CTNNB1</td>
      <td>41266098</td>
      <td>41266098</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>17578</td>
      <td>ENST00000349496.5</td>
      <td>GRCh37</td>
    </tr>
    <tr>
      <th>1</th>
      <td>KIT</td>
      <td>55599320</td>
      <td>55599321</td>
      <td>GA</td>
      <td>AT</td>
      <td>DNP</td>
      <td>MAF</td>
      <td>375926</td>
      <td>ENST00000288135.5</td>
      <td>GRCh37</td>
    </tr>
    <tr>
      <th>2</th>
      <td>APC</td>
      <td>112174150</td>
      <td>112174150</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>233866</td>
      <td>ENST00000457016.1</td>
      <td>GRCh37</td>
    </tr>
    <tr>
      <th>3</th>
      <td>APC</td>
      <td>112174154</td>
      <td>112174155</td>
      <td>-</td>
      <td>A</td>
      <td>INS</td>
      <td>MAF</td>
      <td>1796995</td>
      <td>ENST00000457016.1</td>
      <td>GRCh37</td>
    </tr>
    <tr>
      <th>4</th>
      <td>MRE11</td>
      <td>94200987</td>
      <td>94200987</td>
      <td>G</td>
      <td>A</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>140953</td>
      <td>None</td>
      <td>GRCh37</td>
    </tr>
    <tr>
      <th>5</th>
      <td>TP53</td>
      <td>7572148</td>
      <td>7572148</td>
      <td>G</td>
      <td>-</td>
      <td>DEL</td>
      <td>MAF</td>
      <td>325626</td>
      <td>ENST00000269305.4</td>
      <td>GRCh37</td>
    </tr>
    <tr>
      <th>6</th>
      <td>BRCA1</td>
      <td>41244325</td>
      <td>41244325</td>
      <td>T</td>
      <td>C</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>573320</td>
      <td>ENST00000471181.2</td>
      <td>GRCh37</td>
    </tr>
  </tbody>
</table>
</div>



### Query the gene orthologs in mouse  


```python
df_queried, df_fail = h2m.query_batch(df, direction='h2m')
```

    No error occurs.



```python
df_queried
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>type_h</th>
      <th>format</th>
      <th>index</th>
      <th>tx_id_h</th>
      <th>ref_genome_h</th>
      <th>gene_name_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>APC</td>
      <td>112174150</td>
      <td>112174150</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>233866</td>
      <td>ENST00000457016.1</td>
      <td>GRCh37</td>
      <td>Apc</td>
    </tr>
    <tr>
      <th>1</th>
      <td>APC</td>
      <td>112174154</td>
      <td>112174155</td>
      <td>-</td>
      <td>A</td>
      <td>INS</td>
      <td>MAF</td>
      <td>1796995</td>
      <td>ENST00000457016.1</td>
      <td>GRCh37</td>
      <td>Apc</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BRCA1</td>
      <td>41244325</td>
      <td>41244325</td>
      <td>T</td>
      <td>C</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>573320</td>
      <td>ENST00000471181.2</td>
      <td>GRCh37</td>
      <td>Brca1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CTNNB1</td>
      <td>41266098</td>
      <td>41266098</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>17578</td>
      <td>ENST00000349496.5</td>
      <td>GRCh37</td>
      <td>Ctnnb1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>KIT</td>
      <td>55599320</td>
      <td>55599321</td>
      <td>GA</td>
      <td>AT</td>
      <td>DNP</td>
      <td>MAF</td>
      <td>375926</td>
      <td>ENST00000288135.5</td>
      <td>GRCh37</td>
      <td>Kit</td>
    </tr>
    <tr>
      <th>5</th>
      <td>MRE11</td>
      <td>94200987</td>
      <td>94200987</td>
      <td>G</td>
      <td>A</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>140953</td>
      <td>None</td>
      <td>GRCh37</td>
      <td>Mre11a</td>
    </tr>
    <tr>
      <th>6</th>
      <td>TP53</td>
      <td>7572148</td>
      <td>7572148</td>
      <td>G</td>
      <td>-</td>
      <td>DEL</td>
      <td>MAF</td>
      <td>325626</td>
      <td>ENST00000269305.4</td>
      <td>GRCh37</td>
      <td>Trp53</td>
    </tr>
  </tbody>
</table>
</div>



### Get canonical transcript IDs for the murine genes  


```python
df_queried, df_fail = h2m.get_tx_batch(df_queried, species='m')
```

    No error occurs.



```python
df_queried
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>type_h</th>
      <th>format</th>
      <th>index</th>
      <th>tx_id_h</th>
      <th>ref_genome_h</th>
      <th>gene_name_m</th>
      <th>tx_id_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>APC</td>
      <td>112174150</td>
      <td>112174150</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>233866</td>
      <td>ENST00000457016.1</td>
      <td>GRCh37</td>
      <td>Apc</td>
      <td>ENSMUST00000079362.13</td>
    </tr>
    <tr>
      <th>1</th>
      <td>APC</td>
      <td>112174154</td>
      <td>112174155</td>
      <td>-</td>
      <td>A</td>
      <td>INS</td>
      <td>MAF</td>
      <td>1796995</td>
      <td>ENST00000457016.1</td>
      <td>GRCh37</td>
      <td>Apc</td>
      <td>ENSMUST00000079362.13</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BRCA1</td>
      <td>41244325</td>
      <td>41244325</td>
      <td>T</td>
      <td>C</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>573320</td>
      <td>ENST00000471181.2</td>
      <td>GRCh37</td>
      <td>Brca1</td>
      <td>ENSMUST00000017290.11</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CTNNB1</td>
      <td>41266098</td>
      <td>41266098</td>
      <td>A</td>
      <td>G</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>17578</td>
      <td>ENST00000349496.5</td>
      <td>GRCh37</td>
      <td>Ctnnb1</td>
      <td>ENSMUST00000007130.15</td>
    </tr>
    <tr>
      <th>4</th>
      <td>KIT</td>
      <td>55599320</td>
      <td>55599321</td>
      <td>GA</td>
      <td>AT</td>
      <td>DNP</td>
      <td>MAF</td>
      <td>375926</td>
      <td>ENST00000288135.5</td>
      <td>GRCh37</td>
      <td>Kit</td>
      <td>ENSMUST00000005815.7</td>
    </tr>
    <tr>
      <th>5</th>
      <td>MRE11</td>
      <td>94200987</td>
      <td>94200987</td>
      <td>G</td>
      <td>A</td>
      <td>SNP</td>
      <td>MAF</td>
      <td>140953</td>
      <td>None</td>
      <td>GRCh37</td>
      <td>Mre11a</td>
      <td>ENSMUST00000034405.11</td>
    </tr>
    <tr>
      <th>6</th>
      <td>TP53</td>
      <td>7572148</td>
      <td>7572148</td>
      <td>G</td>
      <td>-</td>
      <td>DEL</td>
      <td>MAF</td>
      <td>325626</td>
      <td>ENST00000269305.4</td>
      <td>GRCh37</td>
      <td>Trp53</td>
      <td>ENSMUST00000108658.10</td>
    </tr>
  </tbody>
</table>
</div>



### Compute the muerine variant equivalents  


```python
df_result, df_fail = h2m.model_batch(df_queried, records_h, index_list_h, records_m, index_list_m, db_h, db_m, 37)
```

    No error occurs.



```python
df_result
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>gene_id_h</th>
      <th>tx_id_h</th>
      <th>chr_h</th>
      <th>exon_num_h</th>
      <th>strand_h</th>
      <th>match</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>...</th>
      <th>alt_seq_m_ori</th>
      <th>HGVSc_m_ori</th>
      <th>HGVSp_m_ori</th>
      <th>start_m</th>
      <th>end_m</th>
      <th>ref_seq_m</th>
      <th>alt_seq_m</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
      <th>index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CTNNB1</td>
      <td>ENSG00000168036.12</td>
      <td>ENST00000349496.5</td>
      <td>chr3</td>
      <td>14</td>
      <td>+</td>
      <td>True</td>
      <td>41266098</td>
      <td>41266098</td>
      <td>A</td>
      <td>...</td>
      <td>G</td>
      <td>ENSMUST00000007130.15:c.95A&gt;G</td>
      <td>D32G</td>
      <td>120779669</td>
      <td>120779669</td>
      <td>A</td>
      <td>G</td>
      <td>ENSMUST00000007130.15:c.95A&gt;G</td>
      <td>D32G</td>
      <td>17578</td>
    </tr>
    <tr>
      <th>1</th>
      <td>APC</td>
      <td>ENSG00000134982.12</td>
      <td>ENST00000457016.1</td>
      <td>chr5</td>
      <td>15</td>
      <td>+</td>
      <td>True</td>
      <td>112174154</td>
      <td>112174155</td>
      <td>-</td>
      <td>...</td>
      <td>A</td>
      <td>ENSMUST00000079362.13:c.2857_2858&gt;A</td>
      <td>Y954Ifs*2</td>
      <td>34445962</td>
      <td>34445963</td>
      <td>-</td>
      <td>A</td>
      <td>ENSMUST00000079362.13:c.2857_2858&gt;A</td>
      <td>Y954Ifs*2</td>
      <td>1796995</td>
    </tr>
    <tr>
      <th>2</th>
      <td>APC</td>
      <td>ENSG00000134982.12</td>
      <td>ENST00000457016.1</td>
      <td>chr5</td>
      <td>15</td>
      <td>+</td>
      <td>True</td>
      <td>112174150</td>
      <td>112174150</td>
      <td>A</td>
      <td>...</td>
      <td>G</td>
      <td>ENSMUST00000079362.13:c.2853A&gt;G</td>
      <td>K951delinsK</td>
      <td>34445958</td>
      <td>34445958</td>
      <td>A</td>
      <td>G</td>
      <td>ENSMUST00000079362.13:c.2853A&gt;G</td>
      <td>K951delinsK</td>
      <td>233866</td>
    </tr>
    <tr>
      <th>3</th>
      <td>TP53</td>
      <td>ENSG00000141510.11</td>
      <td>ENST00000269305.4</td>
      <td>chr17</td>
      <td>10</td>
      <td>-</td>
      <td>True</td>
      <td>7572148</td>
      <td>7572148</td>
      <td>G</td>
      <td>...</td>
      <td>-</td>
      <td>ENSMUST00000108658.10:c.279C&gt;</td>
      <td>None</td>
      <td>69482536</td>
      <td>69482536</td>
      <td>C</td>
      <td>-</td>
      <td>ENSMUST00000108658.10:c.279C&gt;</td>
      <td>None</td>
      <td>325626</td>
    </tr>
    <tr>
      <th>4</th>
      <td>KIT</td>
      <td>ENSG00000157404.11</td>
      <td>ENST00000288135.5</td>
      <td>chr4</td>
      <td>21</td>
      <td>+</td>
      <td>True</td>
      <td>55599320</td>
      <td>55599321</td>
      <td>GA</td>
      <td>...</td>
      <td>AT</td>
      <td>ENSMUST00000005815.7:c.2452_2453GA&gt;AT</td>
      <td>D818I</td>
      <td>75810291</td>
      <td>75810292</td>
      <td>GA</td>
      <td>AT</td>
      <td>ENSMUST00000005815.7:c.2452_2453GA&gt;AT</td>
      <td>D818I</td>
      <td>375926</td>
    </tr>
    <tr>
      <th>5</th>
      <td>BRCA1</td>
      <td>ENSG00000012048.15</td>
      <td>ENST00000471181.2</td>
      <td>chr17</td>
      <td>23</td>
      <td>-</td>
      <td>True</td>
      <td>41244325</td>
      <td>41244325</td>
      <td>T</td>
      <td>...</td>
      <td>C</td>
      <td>ENSMUST00000017290.11:c.3130A&gt;G</td>
      <td>N1044D</td>
      <td>101415003</td>
      <td>101415003</td>
      <td>T</td>
      <td>C</td>
      <td>ENSMUST00000017290.11:c.3130A&gt;G</td>
      <td>N1044D</td>
      <td>573320</td>
    </tr>
  </tbody>
</table>
<p>6 rows × 43 columns</p>
</div>



## Single variant input  

### Query orthologous genes
First of all, you can use H2M to query a human gene for the presence of mouse homologs and vice versa.  


```python
query_result = h2m.query('TP53')
```

    Query human gene: TP53;
    Mouse ortholog(s): Trp53;
    Homology type: one2one;
    Sequence Simalarity(%):77.3537.



```python
query_result = h2m.query('Trp53', direction='m2h')
```

    Query human gene: Trp53;
    Mouse ortholog(s): TP53;
    Homology type: one2one;
    Sequence Simalarity(%):77.3537.


The output is a list of information for all the mouse ortholog(s) (if have; sometimes more than one).  
Each element is a dictionary of **mouse gene name**, **mouse gene id**, **homology type** (one to one/one to multiple/many to many), and **similarity of human and mouse gene in percentage**.


```python
h2m.query('U1')
```

    Query human gene: U1;
    Mouse ortholog(s): Gm22866,Gm25938;
    Homology type: one2many;
    Sequence Simalarity(%):68.75, 62.3457.





    [{'gene_name_m': 'Gm22866',
      'gene_id_m': 'ENSMUSG00000065881',
      'homology_type': 'ortholog_one2many',
      'similarity': 68.75},
     {'gene_name_m': 'Gm25938',
      'gene_id_m': 'ENSMUSG00000077327',
      'homology_type': 'ortholog_one2many',
      'similarity': 62.3457}]




```python
h2m.query('TPT1P6')
```

    The query human gene: TPT1P6 has no mouse ortholog or this gene id is not included in the database. Please check the input format.





    [None]



A print output is helpful in interactive tasks like a Jupyter Notebook. You can also turn off this.

Except for gene names, both ENSEMBL gene id and transcript id are accepted to identify a human gene. You can use the **ty** parameter ('tx_id','gene_id' or 'name') to specify your input type, but this is totally optional.

Using gene id:


```python
query_result = h2m.query('ENSG00000141510')
```

    Query human gene: TP53;
    Mouse ortholog(s): Trp53;
    Homology type: one2one;
    Sequence Simalarity(%):77.3537.


Using transcript id. Should include a db annotation file with the same ref genome version.


```python
query_result = h2m.query('ENST00000269305.4', db=db_h, ty='tx_id')
```

    Query human gene: TP53;
    Mouse ortholog(s): Trp53;
    Homology type: one2one;
    Sequence Simalarity(%):77.3537.


The query result of all human genes, as well as corresponding transcript IDs, is also available as [a csv file]('https://www.dropbox.com/scl/fi/o6735wok25t5dstvpz9kd/Supp_Table_1_Homo_Genes.csv?rlkey=skvxlcfv4r8ksspiq5itxxjc9&dl=0').

### Get transcript ID (Internet connection needed)

One gene may have different transcripts. For mutation modeling, it is important to specify one transcript. If you do not have this information in hand, you can use H2M to get it.

Again, both gene IDs and gene names are accepted as identificaitons for human and mouse genes.


```python
list_tx_id_h = h2m.get_tx_id('TP53', 'h', ver=37)
```

    Genome assembly: GRCh37;
    The canonical transcript is: ENST00000269305.4;
    You can choose from the 17 transcripts below for further analysis:
    (1)ENST00000269305.4 (2)ENST00000413465.2 (3)ENST00000359597.4 (4)ENST00000504290.1 (5)ENST00000510385.1 (6)ENST00000504937.1 (7)ENST00000455263.2 (8)ENST00000420246.2 (9)ENST00000445888.2 (10)ENST00000576024.1 (11)ENST00000509690.1 (12)ENST00000514944.1 (13)ENST00000574684.1 (14)ENST00000505014.1 (15)ENST00000508793.1 (16)ENST00000604348.1 (17)ENST00000503591.1
    



```python
list_tx_id_m = h2m.get_tx_id('ENSMUSG00000059552', 'm')
```

    Genome assembly: GRCm39;
    The canonical transcript is: ENSMUST00000108658.10;
    You can choose from the 6 transcripts below for further analysis:
    (1)ENSMUST00000108658.10 (2)ENSMUST00000171247.8 (3)ENSMUST00000005371.12 (4)ENSMUST00000147512.2 (5)ENSMUST00000108657.4 (6)ENSMUST00000130540.2
    


More information is offered except for a complete list of all transcripts.

- the chromosome


```python
list_tx_id_h[0]
```




    17



- the start and end location of the gene on the chromosome


```python
list_tx_id_h[1:3]
```




    [7565097, 7590856]



- the canonical transcript annotated by ENSEMBL database and used by major clinical datasets (e.g. AACR-GENIE)


```python
list_tx_id_h[3]
```




    'ENST00000269305.4'



- the complete list of transcripts starting with the canonical one


```python
list_tx_id_h[4]
```




    ['ENST00000269305.4',
     'ENST00000413465.2',
     'ENST00000359597.4',
     'ENST00000504290.1',
     'ENST00000510385.1',
     'ENST00000504937.1',
     'ENST00000455263.2',
     'ENST00000420246.2',
     'ENST00000445888.2',
     'ENST00000576024.1',
     'ENST00000509690.1',
     'ENST00000514944.1',
     'ENST00000574684.1',
     'ENST00000505014.1',
     'ENST00000508793.1',
     'ENST00000604348.1',
     'ENST00000503591.1']



The output is a list of information for all the mouse ortholog(s) (if have; sometimes more than one).  
Each element is a dictionary of **mouse gene name**, **mouse gene id**, **homology type** (one to one/one to multiple), and **similarity of human and mouse gene in percentage**.


### Modeling human variants in the mouse genome

#### Input Parameters Guidance  

Now you can use H2M to model your human mutations of interest.  
You should have at least such information in hand:  
1. transcript id of the human gene
2. transcript id of the mouse gene    

With MAF file:   

3. start location of human variants *on the chromosome*  
4. end location of human variants *on the chromosome*  
5. the reference and alternative sequence *on the positive strand of the chromosome*  
6. mutation type in `'SNP','DNP','TNP','ONP','INS','DEL'` 
7. the version number of human ref genome

#### Usage

Taking *TP53* R273H (ENST00000269305.4:c.818G>A) as an example.


```python
tx_id_h, tx_id_m = list_tx_id_h[3], list_tx_id_m[3]
# use the canonical transcript
```

Another non-coding example.


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7578291,7578291, 'C','T', ty_h = 'SNP', ver = 37)
pd.DataFrame(model_result)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>gene_id_h</th>
      <th>tx_id_h</th>
      <th>chr_h</th>
      <th>exon_num_h</th>
      <th>strand_h</th>
      <th>match</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>...</th>
      <th>ref_seq_m_ori</th>
      <th>alt_seq_m_ori</th>
      <th>HGVSc_m_ori</th>
      <th>HGVSp_m_ori</th>
      <th>start_m</th>
      <th>end_m</th>
      <th>ref_seq_m</th>
      <th>alt_seq_m</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>TP53</td>
      <td>ENSG00000141510.11</td>
      <td>ENST00000269305.4</td>
      <td>chr17</td>
      <td>10</td>
      <td>-</td>
      <td>False</td>
      <td>7578291</td>
      <td>7578291</td>
      <td>T</td>
      <td>...</td>
      <td>A</td>
      <td>A</td>
      <td>ENSMUST00000108658.10:c.551-2A&gt;A</td>
      <td>X183_splice</td>
      <td>69479450</td>
      <td>69479450</td>
      <td>A</td>
      <td>A</td>
      <td>ENSMUST00000108658.10:c.551-2A&gt;A</td>
      <td>X183_splice</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 42 columns</p>
</div>



We can see that this human mutaton can be originally modeled by introducing the same neucleotide alteration.

By setting `show_sequence = True`, we can output the sequences of the wild-type and mutated human gene, wild-type, originally-modeled, and alternatively-modeled (if exsist) mouse gene. (If it can be originally modeled, the new_seq_m_alt would be the same as new_seq_m_ori)


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577120, 7577120, 'C','T', ty_h = 'SNP', ver = 37, show_sequence=True)
pd.DataFrame(model_result)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>gene_id_h</th>
      <th>tx_id_h</th>
      <th>chr_h</th>
      <th>exon_num_h</th>
      <th>strand_h</th>
      <th>match</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>...</th>
      <th>seq_m</th>
      <th>new_seq_m_ori</th>
      <th>mouse_tx_idx_ori</th>
      <th>mouse_p_idx_ori</th>
      <th>mouse_new_p_idx_ori</th>
      <th>dist_m</th>
      <th>new_seq_m</th>
      <th>mouse_tx_idx</th>
      <th>mouse_p_idx</th>
      <th>mouse_new_p_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>TP53</td>
      <td>ENSG00000141510.11</td>
      <td>ENST00000269305.4</td>
      <td>chr17</td>
      <td>10</td>
      <td>-</td>
      <td>True</td>
      <td>7577120</td>
      <td>7577120</td>
      <td>C</td>
      <td>...</td>
      <td>ATGACTGCCATGGAGGAGTCACAGTCGGATATCAGCCTCGAGCTCC...</td>
      <td>ATGACTGCCATGGAGGAGTCACAGTCGGATATCAGCCTCGAGCTCC...</td>
      <td>[808]</td>
      <td>[269]</td>
      <td>[269]</td>
      <td>None</td>
      <td>ATGACTGCCATGGAGGAGTCACAGTCGGATATCAGCCTCGAGCTCC...</td>
      <td>[808]</td>
      <td>[269]</td>
      <td>[269]</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 58 columns</p>
</div>



Modeling results with `show_sequence = True` can be directly visulaized by `h2m.visulization`.  


```python
h2m.visualization(model_result, flank_size=4, print_size=2)
```


    
![png](h2m_package_tutorial_files/h2m_package_tutorial_81_0.png)
    


The length of the identical sequences between human and mouse on teh left/right side of the mutation is provided in order to give you a sense of the local homology and how confident you should be in the fidelity of this modeling.  


```python
pd.DataFrame(model_result)[['flank_size_left','flank_size_right']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>flank_size_left</th>
      <th>flank_size_right</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>4aa</td>
      <td>15aa</td>
    </tr>
  </tbody>
</table>
</div>



### Alternative modeling

Sometimes the human mutation cannot be originally modeled in the mouse genome by using the same neucleotide alteration. Under this circumsatance, some alternative modeling strategies may be found by searching the codon list of the target amino acids. 

Taking TP53 R306Q as an example. 


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577021, 7577021, 'C','T', ty_h = 'SNP', ver = 37)
pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                            'HGVSc_m_ori','HGVSp_m_ori',
                            'HGVSc_m','HGVSp_m']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>HGVSc_h</th>
      <th>HGVSp_h</th>
      <th>HGVSc_m_ori</th>
      <th>HGVSp_m_ori</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENST00000269305.4:c.917G&gt;A</td>
      <td>R306Q</td>
      <td>ENSMUST00000108658.10:c.908G&gt;A</td>
      <td>R303K</td>
      <td>ENSMUST00000108658.10:c.907_908AG&gt;CA</td>
      <td>R303Q</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENST00000269305.4:c.917G&gt;A</td>
      <td>R306Q</td>
      <td>ENSMUST00000108658.10:c.908G&gt;A</td>
      <td>R303K</td>
      <td>ENSMUST00000108658.10:c.907_909AGA&gt;CAG</td>
      <td>R303Q</td>
    </tr>
  </tbody>
</table>
</div>




Taking TP53 R249_T253delinsS as an example.


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577523, 7577534, 'GTGAGGATGGGC', '-', ty_h = 'DEL', ver = 37)
pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                            'HGVSc_m_ori','HGVSp_m_ori',
                            'HGVSc_m','HGVSp_m']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>HGVSc_h</th>
      <th>HGVSp_h</th>
      <th>HGVSc_m_ori</th>
      <th>HGVSp_m_ori</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENST00000269305.4:c.747_758GCCCATCCTCAC&gt;</td>
      <td>R249_T253delinsS</td>
      <td>ENSMUST00000108658.10:c.738_749ACCTATCCTTAC&gt;</td>
      <td>P247_T250del</td>
      <td>ENSMUST00000108658.10:c.736_748CGACCTATCCTTA&gt;T</td>
      <td>R246_T250delinsS</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENST00000269305.4:c.747_758GCCCATCCTCAC&gt;</td>
      <td>R249_T253delinsS</td>
      <td>ENSMUST00000108658.10:c.738_749ACCTATCCTTAC&gt;</td>
      <td>P247_T250del</td>
      <td>ENSMUST00000108658.10:c.736_749CGACCTATCCTTAC&gt;AG</td>
      <td>R246_T250delinsS</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENST00000269305.4:c.747_758GCCCATCCTCAC&gt;</td>
      <td>R249_T253delinsS</td>
      <td>ENSMUST00000108658.10:c.738_749ACCTATCCTTAC&gt;</td>
      <td>P247_T250del</td>
      <td>ENSMUST00000108658.10:c.736_750CGACCTATCCTTACC...</td>
      <td>R246_T250delinsS</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENST00000269305.4:c.747_758GCCCATCCTCAC&gt;</td>
      <td>R249_T253delinsS</td>
      <td>ENSMUST00000108658.10:c.738_749ACCTATCCTTAC&gt;</td>
      <td>P247_T250del</td>
      <td>ENSMUST00000108658.10:c.736_750CGACCTATCCTTACC...</td>
      <td>R246_T250delinsS</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENST00000269305.4:c.747_758GCCCATCCTCAC&gt;</td>
      <td>R249_T253delinsS</td>
      <td>ENSMUST00000108658.10:c.738_749ACCTATCCTTAC&gt;</td>
      <td>P247_T250del</td>
      <td>ENSMUST00000108658.10:c.736_750CGACCTATCCTTACC...</td>
      <td>R246_T250delinsS</td>
    </tr>
  </tbody>
</table>
</div>



The default maximum number of output alternatives is 5. You can definitly change that by the parameter **max_alternative**.


```python
model_result_long = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577523, 7577534, 'GTGAGGATGGGC', '-', ty_h = 'DEL', ver = 37, max_alternative=10)
len(model_result), len(model_result_long)
```




    (5, 6)



If you do not want to alternatively model variants, you can set **search_alternatve** to False.


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577523, 7577534, 'GTGAGGATGGGC', '-', ty_h = 'DEL', ver = 37, search_alternative= False)
model_result[0]['statement']
```




    'Class 6: This mutation cannot be originally modeled.'



### Original modeling with uncertain effects

For frame-shifting mutations and mutations in the non-coding region, we cannot find such alternative modeling strategies with the same protein change effects. H2M will only offer the original modeling and its effect.

- Example 1: *TP53* C275Lfs*31


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577115, 7577116, '','A', ty_h = 'INS', ver = 37)
pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                            'HGVSc_m','HGVSp_m',
                            'statement']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>HGVSc_h</th>
      <th>HGVSp_h</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
      <th>statement</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENST00000269305.4:c.822_823&gt;T</td>
      <td>C275Lfs*31</td>
      <td>ENSMUST00000108658.10:c.813_814&gt;T</td>
      <td>C272Lfs*24</td>
      <td>Class 2: This mutation can be modeled, but the...</td>
    </tr>
  </tbody>
</table>
</div>



- Example 2: TP53 splice site mutation


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7578555, 7578555, 'C', 'T', ty_h = 'SNP', ver = 37)
pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                            'HGVSc_m','HGVSp_m',
                            'statement']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>HGVSc_h</th>
      <th>HGVSp_h</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
      <th>statement</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENST00000269305.4:c.376-1G&gt;A</td>
      <td>X125_splice</td>
      <td>ENSMUST00000108658.10:c.367-1G&gt;A</td>
      <td>X122_splice</td>
      <td>Class 2: This mutation can be modeled, but the...</td>
    </tr>
  </tbody>
</table>
</div>



## Additional Usage Hint   

### Additional function 1: modeling M2H

Replace human variant coordinates and sequences with murine ones, and set `direction = 'm2h'`.  Use TP53 R273H as an example.  

#### H2M:  


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577120, 7577120, 'C','T', ty_h = 'SNP', ver = 37)
pd.DataFrame(model_result)[['start_h','end_h','ref_seq_h','alt_seq_h','HGVSp_h','start_m','end_m','ref_seq_m','alt_seq_m','HGVSp_m']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>HGVSp_h</th>
      <th>start_m</th>
      <th>end_m</th>
      <th>ref_seq_m</th>
      <th>alt_seq_m</th>
      <th>HGVSp_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>7577120</td>
      <td>7577120</td>
      <td>C</td>
      <td>T</td>
      <td>R273H</td>
      <td>69480434</td>
      <td>69480434</td>
      <td>G</td>
      <td>A</td>
      <td>R270H</td>
    </tr>
  </tbody>
</table>
</div>



#### M2H:  


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 
                         69480434, 69480434, 'G', 'A', ty_h = 'SNP', ver = 37, 
                         direction='m2h')
pd.DataFrame(model_result)[['start_h','end_h','ref_seq_h','alt_seq_h','HGVSp_h','start_m','end_m','ref_seq_m','alt_seq_m','HGVSp_m']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>HGVSp_h</th>
      <th>start_m</th>
      <th>end_m</th>
      <th>ref_seq_m</th>
      <th>alt_seq_m</th>
      <th>HGVSp_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>7577120</td>
      <td>7577120</td>
      <td>C</td>
      <td>T</td>
      <td>R273H</td>
      <td>69480434</td>
      <td>69480434</td>
      <td>G</td>
      <td>A</td>
      <td>R270H</td>
    </tr>
  </tbody>
</table>
</div>



### Additional function 2: modeling H2H/M2M paralogs  

Replace the reference genome and gencode annotation database input parameter to do so.  Take human IDH1 R172G as an example.  


```python
df = df[df['class']==1].reset_index(drop=True)
```


```python
tx_id_1_h, tx_id_2_h = h2m.get_tx_id('SMARCA2','h',ver=37)[3],h2m.get_tx_id('SMARCA4','h',ver=37)[3]
```

    Genome assembly: GRCh37;
    The canonical transcript is: ENST00000382203.1;
    You can choose from the 17 transcripts below for further analysis:
    (1)ENST00000382203.1 (2)ENST00000450198.1 (3)ENST00000457226.1 (4)ENST00000439732.1 (5)ENST00000382194.1 (6)ENST00000491574.1 (7)ENST00000452193.1 (8)ENST00000302401.3 (9)ENST00000423555.1 (10)ENST00000382186.1 (11)ENST00000417599.1 (12)ENST00000382185.1 (13)ENST00000382183.1 (14)ENST00000416751.1 (15)ENST00000349721.2 (16)ENST00000357248.2 (17)ENST00000324954.5
    
    Genome assembly: GRCh37;
    The canonical transcript is: ENST00000429416.3;
    You can choose from the 20 transcripts below for further analysis:
    (1)ENST00000429416.3 (2)ENST00000344626.4 (3)ENST00000541122.2 (4)ENST00000589677.1 (5)ENST00000444061.3 (6)ENST00000590574.1 (7)ENST00000591545.1 (8)ENST00000592604.1 (9)ENST00000586122.1 (10)ENST00000587988.1 (11)ENST00000591595.1 (12)ENST00000585799.1 (13)ENST00000592158.1 (14)ENST00000586892.1 (15)ENST00000538456.3 (16)ENST00000586985.1 (17)ENST00000586921.1 (18)ENST00000358026.2 (19)ENST00000413806.3 (20)ENST00000450717.3
    



```python
pd.DataFrame(model_result)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_m</th>
      <th>gene_id_m</th>
      <th>tx_id_m</th>
      <th>chr_m</th>
      <th>exon_num_m</th>
      <th>strand_m</th>
      <th>matcm</th>
      <th>start_m</th>
      <th>end_m</th>
      <th>ref_seq_m</th>
      <th>...</th>
      <th>ref_seq_h_ori</th>
      <th>alt_seq_h_ori</th>
      <th>HGVSc_h_ori</th>
      <th>HGVSp_h_ori</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>alt_seq_h</th>
      <th>HGVSc_h</th>
      <th>HGVSp_h</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>IDH1</td>
      <td>ENSG00000138413.9</td>
      <td>ENST00000415913.1</td>
      <td>chr2</td>
      <td>8</td>
      <td>-</td>
      <td>True</td>
      <td>90631839</td>
      <td>90631839</td>
      <td>A</td>
      <td>...</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 42 columns</p>
</div>




```python
model_result = h2m.model(records_h,index_list_h, records_h, index_list_h, db_h, db_h, tx_id_1_h, tx_id_2_h, 
                        2115855, 2115855, 'G', 'A', ty_h = 'SNP', ver = 37,
                        direction='h2h')
pd.DataFrame(model_result)[['gene_name_h_1','start_h_1','end_h_1','ref_seq_h_1','alt_seq_h_1','HGVSp_h_1','gene_name_h_2','start_h_2','end_h_2','ref_seq_h_2','alt_seq_h_2','HGVSp_h_2']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h_1</th>
      <th>start_h_1</th>
      <th>end_h_1</th>
      <th>ref_seq_h_1</th>
      <th>alt_seq_h_1</th>
      <th>HGVSp_h_1</th>
      <th>gene_name_h_2</th>
      <th>start_h_2</th>
      <th>end_h_2</th>
      <th>ref_seq_h_2</th>
      <th>alt_seq_h_2</th>
      <th>HGVSp_h_2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>SMARCA2</td>
      <td>2115855</td>
      <td>2115855</td>
      <td>G</td>
      <td>A</td>
      <td>G1164R</td>
      <td>SMARCA4</td>
      <td>11143999</td>
      <td>11143999</td>
      <td>G</td>
      <td>A</td>
      <td>G1194R</td>
    </tr>
  </tbody>
</table>
</div>



### Additional function 3: modeling for base editing

When you set **param = 'BE'**, you will get modeling results that can be modeled by base editing (A->G, G->A, C->T, T->C, AA->GG, ...etc.). If one mutation can be originally modeled in the mouse genome but not in a BE style, alternative BE modeling strategies will be returned too.

Taking *KEAP1* F221L as an example.


```python
h2m.query('KEAP1')
```

    Query human gene: KEAP1;
    Mouse ortholog(s): Keap1;
    Homology type: one2one;
    Sequence Simalarity(%):94.0705.





    [{'gene_name_m': 'Keap1',
      'gene_id_m': 'ENSMUSG00000003308',
      'homology_type': 'ortholog_one2one',
      'similarity': 94.0705}]




```python
tx_id_h_2, tx_id_m_2 = h2m.get_tx_id('KEAP1','h',ver=37, show=False)[3], h2m.get_tx_id('Keap1','m', show=False)[3]
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h_2, tx_id_m_2, 10602915, 10602915, 'G','T', ty_h = 'SNP', ver = 37, param='BE')
```


```python
pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h','HGVSc_m_ori','HGVSp_m_ori','statement','HGVSc_m','HGVSp_m']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>HGVSc_h</th>
      <th>HGVSp_h</th>
      <th>HGVSc_m_ori</th>
      <th>HGVSp_m_ori</th>
      <th>statement</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENST00000171111.5:c.663C&gt;A</td>
      <td>F221L</td>
      <td>ENSMUST00000164812.8:c.663C&gt;A</td>
      <td>F221L</td>
      <td>Class 1: This mutation can be alternatively mo...</td>
      <td>ENSMUST00000164812.8:c.661T&gt;C</td>
      <td>F221L</td>
    </tr>
  </tbody>
</table>
</div>



### Additional function 4: modeling by amino acid change input  

Set **coor = 'aa'** and modeling variants by amino acid change input. Use TP53 R175H as an example.


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 175, 175, 'R', 'H', coor = 'aa', ty_h = 'SNP', ver = 37)
pd.DataFrame(model_result)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_name_h</th>
      <th>gene_id_h</th>
      <th>tx_id_h</th>
      <th>chr_h</th>
      <th>exon_num_h</th>
      <th>strand_h</th>
      <th>match</th>
      <th>start_h</th>
      <th>end_h</th>
      <th>ref_seq_h</th>
      <th>...</th>
      <th>ref_seq_m_ori</th>
      <th>alt_seq_m_ori</th>
      <th>HGVSc_m_ori</th>
      <th>HGVSp_m_ori</th>
      <th>start_m</th>
      <th>end_m</th>
      <th>ref_seq_m</th>
      <th>alt_seq_m</th>
      <th>HGVSc_m</th>
      <th>HGVSp_m</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>TP53</td>
      <td>ENSG00000141510.11</td>
      <td>ENST00000269305.4</td>
      <td>chr17</td>
      <td>10</td>
      <td>-</td>
      <td>True</td>
      <td>7578405</td>
      <td>7578407</td>
      <td>GCG</td>
      <td>...</td>
      <td>CGC</td>
      <td>CAC</td>
      <td>ENSMUST00000108658.10:c.514_516CGC&gt;CAC</td>
      <td>R172H</td>
      <td>69479338</td>
      <td>69479338</td>
      <td>G</td>
      <td>A</td>
      <td>ENSMUST00000108658.10:c.515G&gt;A</td>
      <td>R172H</td>
    </tr>
    <tr>
      <th>1</th>
      <td>TP53</td>
      <td>ENSG00000141510.11</td>
      <td>ENST00000269305.4</td>
      <td>chr17</td>
      <td>10</td>
      <td>-</td>
      <td>True</td>
      <td>7578405</td>
      <td>7578407</td>
      <td>GCG</td>
      <td>...</td>
      <td>CGC</td>
      <td>CAC</td>
      <td>ENSMUST00000108658.10:c.514_516CGC&gt;CAC</td>
      <td>R172H</td>
      <td>69479338</td>
      <td>69479339</td>
      <td>GC</td>
      <td>AT</td>
      <td>ENSMUST00000108658.10:c.515_516GC&gt;AT</td>
      <td>R172H</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 42 columns</p>
</div>



All of these can also be done in a batch-processing style by using `h2m.model_batch`.   

## Appendix  


```python
list(model_result[0].keys())
```




    ['gene_name_h',
     'gene_id_h',
     'tx_id_h',
     'chr_h',
     'exon_num_h',
     'strand_h',
     'match',
     'start_h',
     'end_h',
     'ref_seq_h',
     'alt_seq_h',
     'HGVSc_h',
     'HGVSp_h',
     'classification_h',
     'exon_h',
     'type_h',
     'status',
     'class',
     'statement',
     'flank_size_left',
     'flank_size_right',
     'gene_name_m',
     'gene_id_m',
     'tx_id_m',
     'chr_m',
     'exon_num_m',
     'strand_m',
     'type_m',
     'classification_m',
     'exon_m',
     'start_m_ori',
     'end_m_ori',
     'ref_seq_m_ori',
     'alt_seq_m_ori',
     'HGVSc_m_ori',
     'HGVSp_m_ori',
     'start_m',
     'end_m',
     'ref_seq_m',
     'alt_seq_m',
     'HGVSc_m',
     'HGVSp_m']




```python
list_of_name = ['gene_name_h',
 'gene_id_h',
 'tx_id_h',
 'chr_h',
 'exon_num_h',
 'strand_h',
 'match',
 'start_h | end_h',
 'ref_seq_h | alt_seq_h',
 'HGVSc_h | HGVSp_h',
 'classification_h',
 'exon_h',
 'type_h',
 'status',
 'class',
 'statement',
 'flank_size_left | flank_size_right',
 'gene_name_m',
 'gene_id_m',
 'tx_id_m',
 'chr_m',
 'exon_num_m',
 'strand_m',
 'type_m',
 'classification_m',
 'exon_m',
 'start_m_ori | end_m_ori',
 'ref_seq_m_ori | alt_seq_m_ori',
 'HGVSc_m_ori | HGVSp_m_ori',
 'start_m | end_m',
 'ref_seq_m | alt_seq_m',
 'HGVSc_m | HGVSp_m']
```


```python
annotation_of_column = [
    'Human gene name',
    'Human gene ID',
    'Human transcript ID',
    'Human chromosome number',
    'Total number of exons of the human transcript',
    '+ or - strand of the human transcript on the chromosome',
    'The computed reference sequence by given coordinate is matched with the input reference sequence or not',
    'Start and end position of the human variant on the chromosome in MAF format',
    'Reference and alternate sequence of the human variant on the chromosome in MAF format',
    'HGVSc and HGVSp expression of the human variant',
    'Human variant effect classification, including missense/nonsense/in-frame indel/fram-shift indel/intron, etc.',
    'Exon/Intron location of the given human mutation, for example, E_7/I_5',
    'Human variant type in MAF format, including SNP/DNP/TNO/ONP/INS/DEL',
    'This mutation can be modeled in the given target transcript or not, True or False',
    'H2M modeling result class, 0-5',
    'Statement of the H2M result class',
    'Length of the identical sequences between human and mouse on teh left/right side of the mutation',
    'Mouse gene name',
    'Mouse gene ID',
    'Mouse transcript ID',
    'Mouse chromosome number',
    'Total number of exons of the mouse transcript',
    '+ or - strand of the mouse transcript on the chromosome',
    'Mouse variant type in MAF format, including SNP/DNP/TNO/ONP/INS/DEL',
    'Mouse variant effect classification',
    'Exon/Intron location of the murine mutation',
    'Start and end position of the mouse variant (with exactly the same DNA change) on the chromosome in MAF format',
    'Reference and alternate sequence of the mouse variant (with exactly the same DNA change) on the chromosome in MAF format',
    'HGVSc and HGVSp expression of the mouse variant (with exactly the same DNA change)',
    'Start and end position of the mouse variant (with the same amino acid change) on the chromosome in MAF format',
    'Reference and alternate sequence of the mouse variant (with the same amino acid change) on the chromosome in MAF format',
    'HGVSc and HGVSp expression of the mouse variant (with the same amino acid change)'
]
```


```python
df = pd.DataFrame(zip(list_of_name, annotation_of_column))
df.columns = ['Column','Annotation']
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Column</th>
      <th>Annotation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>gene_name_h</td>
      <td>Human gene name</td>
    </tr>
    <tr>
      <th>1</th>
      <td>gene_id_h</td>
      <td>Human gene ID</td>
    </tr>
    <tr>
      <th>2</th>
      <td>tx_id_h</td>
      <td>Human transcript ID</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr_h</td>
      <td>Human chromosome number</td>
    </tr>
    <tr>
      <th>4</th>
      <td>exon_num_h</td>
      <td>Total number of exons of the human transcript</td>
    </tr>
    <tr>
      <th>5</th>
      <td>strand_h</td>
      <td>+ or - strand of the human transcript on the c...</td>
    </tr>
    <tr>
      <th>6</th>
      <td>match</td>
      <td>The computed reference sequence by given coord...</td>
    </tr>
    <tr>
      <th>7</th>
      <td>start_h | end_h</td>
      <td>Start and end position of the human variant on...</td>
    </tr>
    <tr>
      <th>8</th>
      <td>ref_seq_h | alt_seq_h</td>
      <td>Reference and alternate sequence of the human ...</td>
    </tr>
    <tr>
      <th>9</th>
      <td>HGVSc_h | HGVSp_h</td>
      <td>HGVSc and HGVSp expression of the human variant</td>
    </tr>
    <tr>
      <th>10</th>
      <td>classification_h</td>
      <td>Human variant effect classification, including...</td>
    </tr>
    <tr>
      <th>11</th>
      <td>exon_h</td>
      <td>Exon/Intron location of the given human mutati...</td>
    </tr>
    <tr>
      <th>12</th>
      <td>type_h</td>
      <td>Human variant type in MAF format, including SN...</td>
    </tr>
    <tr>
      <th>13</th>
      <td>status</td>
      <td>This mutation can be modeled in the given targ...</td>
    </tr>
    <tr>
      <th>14</th>
      <td>class</td>
      <td>H2M modeling result class, 0-5</td>
    </tr>
    <tr>
      <th>15</th>
      <td>statement</td>
      <td>Statement of the H2M result class</td>
    </tr>
    <tr>
      <th>16</th>
      <td>flank_size_left | flank_size_right</td>
      <td>Length of the identical sequences between huma...</td>
    </tr>
    <tr>
      <th>17</th>
      <td>gene_name_m</td>
      <td>Mouse gene name</td>
    </tr>
    <tr>
      <th>18</th>
      <td>gene_id_m</td>
      <td>Mouse gene ID</td>
    </tr>
    <tr>
      <th>19</th>
      <td>tx_id_m</td>
      <td>Mouse transcript ID</td>
    </tr>
    <tr>
      <th>20</th>
      <td>chr_m</td>
      <td>Mouse chromosome number</td>
    </tr>
    <tr>
      <th>21</th>
      <td>exon_num_m</td>
      <td>Total number of exons of the mouse transcript</td>
    </tr>
    <tr>
      <th>22</th>
      <td>strand_m</td>
      <td>+ or - strand of the mouse transcript on the c...</td>
    </tr>
    <tr>
      <th>23</th>
      <td>type_m</td>
      <td>Mouse variant type in MAF format, including SN...</td>
    </tr>
    <tr>
      <th>24</th>
      <td>classification_m</td>
      <td>Mouse variant effect classification</td>
    </tr>
    <tr>
      <th>25</th>
      <td>exon_m</td>
      <td>Exon/Intron location of the murine mutation</td>
    </tr>
    <tr>
      <th>26</th>
      <td>start_m_ori | end_m_ori</td>
      <td>Start and end position of the mouse variant (w...</td>
    </tr>
    <tr>
      <th>27</th>
      <td>ref_seq_m_ori | alt_seq_m_ori</td>
      <td>Reference and alternate sequence of the mouse ...</td>
    </tr>
    <tr>
      <th>28</th>
      <td>HGVSc_m_ori | HGVSp_m_ori</td>
      <td>HGVSc and HGVSp expression of the mouse varian...</td>
    </tr>
    <tr>
      <th>29</th>
      <td>start_m | end_m</td>
      <td>Start and end position of the mouse variant (w...</td>
    </tr>
    <tr>
      <th>30</th>
      <td>ref_seq_m | alt_seq_m</td>
      <td>Reference and alternate sequence of the mouse ...</td>
    </tr>
    <tr>
      <th>31</th>
      <td>HGVSc_m | HGVSp_m</td>
      <td>HGVSc and HGVSp expression of the mouse varian...</td>
    </tr>
  </tbody>
</table>
</div>




```python
list_of_class = [
    'Class 0',
    'Class 1',
    'Class 2',
    'Class 3',
    'Class 4',
    'Class 5',
    'Class 6'
]

list_of_staterment = [
    'This mutation can be originally modeled.',
    'This mutation can be alternatively modeled.',
    'This mutation can be modeled, but the effect may not be consistent.',
    'This mutation cannot be originally modeled and no alternative is found.',
    'Mutated sequences are not identical.',
    'Coordinate error. This mutation is not in the query gene.',
    'This mutation cannot be originally modeled.'
]

df = pd.DataFrame(zip(list_of_class, list_of_staterment))
df.columns = ['Class','Statement']
```


```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Class</th>
      <th>Statement</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Class 0</td>
      <td>This mutation can be originally modeled.</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Class 1</td>
      <td>This mutation can be alternatively modeled.</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Class 2</td>
      <td>This mutation can be modeled, but the effect m...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Class 3</td>
      <td>This mutation cannot be originally modeled and...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Class 4</td>
      <td>Mutated sequences are not identical.</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Class 5</td>
      <td>Coordinate error. This mutation is not in the ...</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Class 6</td>
      <td>This mutation cannot be originally modeled.</td>
    </tr>
  </tbody>
</table>
</div>


