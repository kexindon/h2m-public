# Quick Start  

```{image} figures/h2m-logo-final.png  
:width: 150px
:align: left
```

```{note}
Data used in this tutorial can be downloaded from <a href='https://www.dropbox.com/scl/fo/1wtrnc9w6s9gemweuw2fv/h?rlkey=hli1z6tv096cjwit5oi6bwggg&dl=0' target="_blank" rel="noopener noreferrer">this Dropbox folder</a>.  
```
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

## Importing packages

```python
    import bioh2m as h2m
    import pandas as pd
```

## Loading data  
We should upload reference genome and GENCODE annotation data for both human and mouse, which could be directly downloaded from <a href='https://www.dropbox.com/scl/fo/1wtrnc9w6s9gemweuw2fv/h?rlkey=hli1z6tv096cjwit5oi6bwggg&dl=0' target="_blank" rel="noopener noreferrer">this Dropbox folder</a>.  
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

## Batch Processing

### Input format  

Common mutation data formats include **MAF** (Mutation Annotation Format, used by cBioPortal), **VCF** (Variant Call Format, used by genomAD), and **ClinVar** (a modified VCF format, used by ClinVar). Mutation coordinates, reference and alternative sequences are recorded in slightly different ways between the three. 

```{image} figures/format.png  
:width: 600px
:align: center
```

In batch processing, **H2M accepts MAF input**. More information about MAF format can be found at <a href='https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#:~:text=Mutation' target="_blank" rel="noopener noreferrer">GDC Documentation</a>. 

For MAF files, you need to build a pandas dataframe with columns as the following example: 

```{image} figures/1.png  
:width: 600px
:align: center
```

For VCF and ClinVar files, you will need to convert the mutation coordinates and sequence information to MAF format after this. 

This can be achieved simply by using H2M built-in functions. 

#### Read from cBioPortal - MAF  
This format is compatible with all of the datasets in the <a href='https://www.cbioportal.org' target="_blank" rel="noopener noreferrer">cBioPortal</a>, as well as <a href='https://www.cancer.gov/ccg/research/genome-sequencing/tcga' target="_blank" rel="noopener noreferrer">TCGA</a> and <a href='https://www.synapse.org/#!Synapse:syn7222066/wiki/405659' target="_blank" rel="noopener noreferrer">AACR-GENIE</a>. Download the txt mutation data file from such public dataset and then load it as follows:  


```python
    path_aacr = '.../data_mutations_extended.txt'
    df = h2m.cbio_reader(path_aacr)
    df
```

```{image} figures/cbio_reader.png
:width: 600px
:align: center
```


#### Read from GenomAD  - VCF 
Search a specific gene in <a href='https://gnomad.broadinstitute.org' target="_blank" rel="noopener noreferrer">GenomAD browser</a>, and download the conluson csv.  

```{image} figures/genomad.png
:width: 600px
:align: center
```

```python
    # downloaded TP53 variants from genomAD
    df = h2m.vcf_reader('.../gnomAD_v4.1.0_ENSG00000141510.csv',keep=False)
    df['gene_name_h'] = 'TP53'
    df
```

```{image} figures/vcf_reader.png
:width: 600px
:align: center
```

Then convert it to **MAF** format.  

```python
    df = h2m.vcf_to_maf(df)
    df
```
```{image} figures/vcf_to_maf.png
:width: 600px
:align: center
```


#### Read from ClinVar
download a <a href='https://www.ncbi.nlm.nih.gov/clinvar/' target="_blank" rel="noopener noreferrer">ClinVar</a> vcf.gz file, and choose your desired Variation IDs that you wish to model. These vcf.gz files are available in <a href='https://www.dropbox.com/scl/fo/1wtrnc9w6s9gemweuw2fv/h?rlkey=hli1z6tv096cjwit5oi6bwggg&dl=0' target="_blank" rel="noopener noreferrer">this Dropbox folder</a>.

```{image} figures/clinvar.png
:width: 600px
:align: center
```

```python
    filepath = '/Users/kexindong/Documents/GitHub/Database/PublicDatabase/ClinVar/GRCh37_clinvar_20240206.vcf.gz'
    variation_ids = [32798013, 375926, 325626, 140953, 233866, 1796995, 17578, 573320]
    df = h2m.clinvar_reader(filepath, variation_ids)
    df = h2m.clinvar_to_maf(df)
    df = df[['gene_name_h',	'start_h','end_h','ref_seq_h','alt_seq_h','type_h','format','ID']]
    df = df.rename(columns={'ID':'index'})
    df
```
```{image} figures/clinvar_result.png
:width: 600px
:align: center
```

### Get canonical transcript IDs for human 

There will be returning two dataframes for success and failures.


```python
    df, df_fail = h2m.get_tx_batch(df, species='h', ver = 37)
    df
```
```{image} figures/2.png
:width: 600px
:align: center
```

### Query orthologous genes


```python
    df_queried, df_fail = h2m.query_batch(df, direction='h2m')
    df_queried
```

```{image} figures/3.png  
:width: 700px
:align: center
```
### Get canonical transcript IDs for mouse


```python
    df_queried, df_fail = h2m.get_tx_batch(df_queried, species='m')
    df_queried
```

```{image} figures/4.png  
:width: 700px
:align: center
```


### Compute the muerine variant equivalents  

```python
    df_result, df_fail = h2m.model_batch(df_queried, records_h, index_list_h, records_m, index_list_m, db_h, db_m, 37)
```


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


```python
    h2m.query('TPT1P6')
```

        The query human gene: TPT1P6 has no mouse ortholog or this gene id is not included in the database. Please check the input format.


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


The query result of all human genes, as well as corresponding transcript IDs, is also available as a csv file in the <a href='https://www.dropbox.com/scl/fo/1wtrnc9w6s9gemweuw2fv/h?rlkey=hli1z6tv096cjwit5oi6bwggg&dl=0' target="_blank" rel="noopener noreferrer">this Dropbox folder</a>.

### Get transcript ID

```{note}
Internet connection needed for this function
```

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
    
Now you can use H2M to model your human mutations of interest.  
You should have at least such information in hand:  
1. Transcript ID of the human gene
2. Transcript ID of the mouse gene    

Also, multiple infomation for the huaman variant in **MAF** format:   

3. Start postion
4. End position
5. Reference sequence
6. Alternate sequence
6. Type in `'SNP','DNP','TNP','ONP','INS','DEL'` 
7. The version number of human ref genome `'37','38'` 


### Modeling human variants in the mouse genome

#### Basic usage  

Taking *TP53* R273H (ENST00000269305.4:c.818G>A) as an example.


```python
    tx_id_h, tx_id_m = list_tx_id_h[3], list_tx_id_m[3]
    # use the canonical transcript
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577120, 7577120, 'C','T', ty_h = 'SNP', ver = 37)
    model_result
```

| Key                   | Value                          |
|-----------------------|--------------------------------|
| gene_name_h           | TP53                           |
| gene_id_h             | ENSG00000141510.11             |
| tx_id_h               | ENST00000269305.4              |
| chr_h                 | chr17                          |
| exon_num_h            | 10                             |
| strand_h              | -                              |
| match                 | True                           |
| start_h               | 7577120                        |
| end_h                 | 7577120                        |
| ref_seq_h             | C                              |
| alt_seq_h             | T                              |
| HGVSc_h               | ENST00000269305.4:c.818G>A     |
| HGVSp_h               | R273H                          |
| classification_h      | Missense                       |
| exon_h                | E_7                            |
| type_h                | SNP                            |
| status                | True                           |
| class                 | 0                              |
| statement             | Class 0: This mutation can be originally modeled. |
| flank_size_left       | 4aa                            |
| flank_size_right      | 15aa                           |
| gene_name_m           | Trp53                          |
| gene_id_m             | ENSMUSG00000059552.14          |
| tx_id_m               | ENSMUST00000108658.10          |
| chr_m                 | chr11                          |
| exon_num_m            | 10                             |
| strand_m              | +                              |
| type_m                | SNP                            |
| classification_m      | Missense                       |
| exon_m                | E_7                            |
| start_m_ori           | 69480434                       |
| end_m_ori             | 69480434                       |
| ref_seq_m_ori         | G                              |
| alt_seq_m_ori         | A                              |
| HGVSc_m_ori           | ENSMUST00000108658.10:c.809G>A |
| HGVSp_m_ori           | R270H                          |
| start_m               | 69480434                       |
| end_m                 | 69480434                       |
| ref_seq_m             | G                              |
| alt_seq_m             | A                              |
| HGVSc_m               | ENSMUST00000108658.10:c.809G>A |
| HGVSp_m               | R270H                          |


We can see that this human mutaton can be originally modeled by introducing the same neucleotide alteration.

##### Flank Size    
The length of the identical sequences between human and mouse on teh left/right side of the mutation is provided in order to give you a sense of the local homology and how confident you should be in the fidelity of this modeling.  

```python
    pd.DataFrame(model_result)[['flank_size_left','flank_size_right']]
```

| flank_size_left       | flank_size_left                |
|-----------------------|--------------------------------|
| 4aa           | 15aa                           |


##### Result visualization  

By setting `show_sequence = True`, we can output the sequences of the wild-type and mutated human gene, wild-type, originally-modeled, and alternatively-modeled (if exsist) mouse gene. Modeling results with `show_sequence = True` can be directly visulaized by `h2m.visulization`.  

```python
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577120, 7577120, 'C','T', ty_h = 'SNP', ver = 37, show_sequence=True)
    pd.DataFrame(model_result)
    h2m.visualization(model_result, flank_size=4, print_size=2)
```

```{image} figures/h2m_visual.png
:width: 800px
:align: center
```

#### Alternative modeling

Sometimes the human mutation cannot be originally modeled in the mouse genome by using the same neucleotide alteration. Under this circumsatance, some alternative modeling strategies may be found by searching the codon list of the target amino acids. 

- Example 1: TP53 R306Q. 


```python
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577021, 7577021, 'C','T', ty_h = 'SNP', ver = 37)
    pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                                'HGVSc_m_ori','HGVSp_m_ori',
                                'HGVSc_m','HGVSp_m']]
```
|     |               |             |
|---------------|--------------------------------------------|---------------------------------------------|
| HGVSc_h       | ENST00000269305.4:c.917G>A                 | ENST00000269305.4:c.917G>A                  |
| HGVSp_h       | R306Q                                      | R306Q                                       |
| HGVSc_m_ori   | ENSMUST00000108658.10:c.908G>A             | ENSMUST00000108658.10:c.908G>A              |
| HGVSp_m_ori   | R303K                                      | R303K                                       |
| HGVSc_m       | ENSMUST00000108658.10:c.907_908AG>CA       | ENSMUST00000108658.10:c.907_909AGA>CAG      |
| HGVSp_m       | R303Q                                      | R303Q                                       |

- Example 2: TP53 R249_T253delinsS.

```python
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577523, 7577534, 'GTGAGGATGGGC', '-', ty_h = 'DEL', ver = 37)
    pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                                'HGVSc_m_ori','HGVSp_m_ori',
                                'HGVSc_m','HGVSp_m']]
```

```{image} figures/delins.png
:width: 800px
:align: center
```

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



#### Original modeling with uncertain effects

For frame-shifting mutations and mutations in the non-coding region, we cannot find such alternative modeling strategies with the same protein change effects. H2M will only offer the original modeling and its effect.

- Example 1: *TP53* C275Lfs*31


```python
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, 
                                db_h, db_m, tx_id_h, tx_id_m, 
                                7577115, 7577116, '','A', ty_h = 'INS', ver = 37)
    pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                                'HGVSc_m','HGVSp_m']]
```

|    |                                             |
|----------|---------------------------------------------|
| HGVSc_h  | ENST00000269305.4:c.822_823>T               |
| HGVSp_h  | C275Lfs*31                                  |
| HGVSc_m  | ENSMUST00000108658.10:c.813_814>T           |
| HGVSp_m  | C272Lfs*24                                  |




- Example 2: TP53 splice site mutation

```python
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, 
                            db_h, db_m, tx_id_h, tx_id_m, 7578555, 7578555, 
                            'C', 'T', ty_h = 'SNP', ver = 37)
    pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h',
                            'HGVSc_m','HGVSp_m']]
```

|         |                                        |
|-----------|--------------------------------------------|
| HGVSc_h   | ENST00000269305.4:c.376-1G>A               |
| HGVSp_h   | X125_splice                                |
| HGVSc_m   | ENSMUST00000108658.10:c.367-1G>A           |
| HGVSp_m   | X122_splice                                |



## Additional Usage Hint   

### Additional function 1: modeling M2H

Replace human variant coordinates and sequences with murine ones, and set `direction = 'm2h'`.  Use TP53 R273H as an example.  

H2M:  
```python
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577120, 7577120, 'C','T', ty_h = 'SNP', ver = 37)
    pd.DataFrame(model_result)[['start_h','end_h','ref_seq_h','alt_seq_h','HGVSp_h','start_m','end_m','ref_seq_m','alt_seq_m','HGVSp_m']]
```

|          |       |
|------------|-----------|
| start_h    | 7577120   |
| end_h      | 7577120   |
| ref_seq_h  | C         |
| alt_seq_h  | T         |
| HGVSp_h    | R273H     |
| start_m    | 69480434  |
| end_m      | 69480434  |
| ref_seq_m  | G         |
| alt_seq_m  | A         |
| HGVSp_m    | R270H     |


M2H:  
```python
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 
                            69480434, 69480434, 'G', 'A', ty_h = 'SNP', ver = 37, 
                            direction='m2h')
    pd.DataFrame(model_result)[['start_h','end_h','ref_seq_h','alt_seq_h','HGVSp_h','start_m','end_m','ref_seq_m','alt_seq_m','HGVSp_m']]
```

|           |       |
|-------------|-----------|
| start_h     | 7577120   |
| end_h       | 7577120   |
| ref_seq_h   | C         |
| alt_seq_h   | T         |
| HGVSp_h     | R273H     |
| start_m     | 69480434  |
| end_m       | 69480434  |
| ref_seq_m   | G         |
| alt_seq_m   | A         |
| HGVSp_m     | R270H     |


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
model_result = h2m.model(records_h,index_list_h, records_h, index_list_h, db_h, db_h, tx_id_1_h, tx_id_2_h, 
                        2115855, 2115855, 'G', 'A', ty_h = 'SNP', ver = 37,
                        direction='h2h')
pd.DataFrame(model_result)[['gene_name_h_1','start_h_1','end_h_1','ref_seq_h_1','alt_seq_h_1','HGVSp_h_1','gene_name_h_2','start_h_2','end_h_2','ref_seq_h_2','alt_seq_h_2','HGVSp_h_2']]
```

|             |         |
|---------------|-------------|
| gene_name_h_1 | SMARCA2     |
| start_h_1     | 2115855     |
| end_h_1       | 2115855     |
| ref_seq_h_1   | G           |
| alt_seq_h_1   | A           |
| HGVSp_h_1     | G1164R      |
| gene_name_h_2 | SMARCA4     |
| start_h_2     | 11143999    |
| end_h_2       | 11143999    |
| ref_seq_h_2   | G           |
| alt_seq_h_2   | A           |
| HGVSp_h_2     | G1194R      |


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





```python
    tx_id_h_2, tx_id_m_2 = h2m.get_tx_id('KEAP1','h',ver=37, show=False)[3], h2m.get_tx_id('Keap1','m', show=False)[3]
    model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h_2, tx_id_m_2, 10602915, 10602915, 'G','T', ty_h = 'SNP', ver = 37, param='BE')
```


```python
    pd.DataFrame(model_result)[['HGVSc_h','HGVSp_h','HGVSc_m_ori','HGVSp_m_ori','statement','HGVSc_m','HGVSp_m']]
```

|            |                                                   |
|--------------|-------------------------------------------------------|
| HGVSc_h      | ENST00000171111.5:c.663C>A                            |
| HGVSp_h      | F221L                                                 |
| HGVSc_m_ori  | ENSMUST00000164812.8:c.663C>A                         |
| HGVSp_m_ori  | F221L                                                 |
| statement    | Class 1: This mutation can be alternatively modeled.  |
| HGVSc_m      | ENSMUST00000164812.8:c.661T>C                         |
| HGVSp_m      | F221L                                                 |


### Additional function 4: modeling by amino acid change input  

Set **coor = 'aa'** and modeling variants by amino acid change input. Use TP53 R175H as an example.


```python
model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 175, 175, 'R', 'H', coor = 'aa', ty_h = 'SNP', ver = 37)
pd.DataFrame(model_result)
```
| Key              | 0                                             | 1                                             |
|------------------|-----------------------------------------------|-----------------------------------------------|
| gene_name_h      | TP53                                          | TP53                                          |
| gene_id_h        | ENSG00000141510.11                            | ENSG00000141510.11                            |
| tx_id_h          | ENST00000269305.4                             | ENST00000269305.4                             |
| chr_h            | chr17                                         | chr17                                         |
| exon_num_h       | 10                                            | 10                                            |
| strand_h         | -                                             | -                                             |
| match            | True                                          | True                                          |
| start_h          | 7578405                                       | 7578405                                       |
| end_h            | 7578407                                       | 7578407                                       |
| ref_seq_m_ori    | CGC                                           | CGC                                           |
| alt_seq_m_ori    | CAC                                           | CAC                                           |
| HGVSc_m_ori      | ENSMUST00000108658.10:c.514_516CGC>CAC        | ENSMUST00000108658.10:c.514_516CGC>CAC        |
| HGVSp_m_ori      | R172H                                         | R172H                                         |
| start_m          | 69479338                                      | 69479338                                      |
| end_m            | 69479338                                      | 69479339                                      |
| ref_seq_m        | G                                             | GC                                            |
| alt_seq_m        | A                                             | AT                                            |
| HGVSc_m          | ENSMUST00000108658.10:c.515G>A                | ENSMUST00000108658.10:c.515_516GC>AT          |
| HGVSp_m          | R172H                                         | R172H                                         |


All of these can also be done in a batch-processing style by using `h2m.model_batch`.   