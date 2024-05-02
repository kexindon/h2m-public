# Complete Documentation 

```{image} figures/h2m-logo-final.png  
:width: 150px
:align: left
```
## H2M Output Data Description

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: middle;
    }

    .dataframe thead th {
        text-align: middle;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th>Column</th>
      <th>Description</th>
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
      <td>Positive or Negative strand of the human transcript on the chromosome</td>
    </tr>
    <tr>
      <th>6</th>
      <td>match</td>
      <td>The computed reference sequence by given coordinate is matched with the input reference sequence or not</td>
    </tr>
    <tr>
      <th>7</th>
      <td>start_h | end_h</td>
      <td>Start and end position of the human variant on the chromosome in MAF format</td>
    </tr>
    <tr>
      <th>8</th>
      <td>ref_seq_h | alt_seq_h</td>
      <td>Reference and alternate sequence of the human variant on the chromosome in MAF format</td>
    </tr>
    <tr>
      <th>9</th>
      <td>HGVSc_h | HGVSp_h</td>
      <td>HGVSc and HGVSp expression of the human variant</td>
    </tr>
    <tr>
      <th>10</th>
      <td>classification_h</td>
      <td>Human variant effect classification, including missense/nonsense/in-frame indel/fram-shift indel/intron, etc.</td>
    </tr>
    <tr>
      <th>11</th>
      <td>exon_h</td>
      <td>Exon/Intron location of the given human mutation, for example, E_7/I_5</td>
    </tr>
    <tr>
      <th>12</th>
      <td>type_h</td>
      <td>Human variant type in MAF format, including SNP/DNP/TNO/ONP/INS/DEL</td>
    </tr>
    <tr>
      <th>13</th>
      <td>status</td>
      <td>This mutation can be modeled in the given target transcript or not, True or False</td>
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
      <td>Length of the identical sequences between human and mouse on the left/right side of the mutation</td>
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
      <td>Positive or negative strand of the mouse transcript on the chromosome</td>
    </tr>
    <tr>
      <th>23</th>
      <td>type_m</td>
      <td>Mouse variant type in MAF format, including SNP/DNP/TNO/ONP/INS/DEL</td>
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
      <td>Start and end position of the mouse variant (with exactly the same DNA change) on the chromosome in MAF format</td>
    </tr>
    <tr>
      <th>27</th>
      <td>ref_seq_m_ori | alt_seq_m_ori</td>
      <td>Reference and alternate sequence of the mouse variant (with exactly the same DNA change) on the chromosome in MAF format</td>
    </tr>
    <tr>
      <th>28</th>
      <td>HGVSc_m_ori | HGVSp_m_ori</td>
      <td>HGVSc and HGVSp expression of the mouse variant (with exactly the same DNA change)</td>
    </tr>
    <tr>
      <th>29</th>
      <td>start_m | end_m</td>
      <td>Start and end position of the mouse variant (with the same amino acid change) on the chromosome in MAF format</td>
    </tr>
    <tr>
      <th>30</th>
      <td>ref_seq_m | alt_seq_m</td>
      <td>Reference and alternate sequence of the mouse variant (with the same amino acid change) on the chromosome in MAF format</td>
    </tr>
    <tr>
      <th>31</th>
      <td>HGVSc_m | HGVSp_m</td>
      <td>HGVSc and HGVSp expression of the mouse variant (with the same amino acid change)</td>
    </tr>
  </tbody>
</table>
</div>


## H2M Modeling Class Description
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: middle;
    }

    .dataframe thead th {
        text-align: left;
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
      <td>This mutation can be modeled, but the effect may not be consistent.</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Class 3</td>
      <td>This mutation cannot be originally modeled and no alternative is found.</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Class 4</td>
      <td>Mutated sequences are not identical.</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Class 5</td>
      <td>Coordinate error. This mutation is not in the query gene.</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Class 6</td>
      <td>This mutation cannot be originally modeled.</td>
    </tr>
  </tbody>
</table>
</div>


## Functions
```{eval-rst}
.. autofunction:: h2m.genome_loader
.. autofunction:: h2m.anno_loader
.. autofunction:: h2m.cbio_reader
.. autofunction:: h2m.clinvar_reader
.. autofunction:: h2m.clinvar_to_maf
.. autofunction:: h2m.vcf_reader
.. autofunction:: h2m.vcf_to_maf
.. autofunction:: h2m.get_variant_type
.. autofunction:: h2m.get_tx_id
.. autofunction:: h2m.get_tx_batch
.. autofunction:: h2m.query
.. autofunction:: h2m.query_batch
.. autofunction:: h2m.model
.. autofunction:: h2m.model_batch
.. autofunction:: h2m.visualization
```