import pysam

def clinvar_reader(path, list_of_ids = None, keep = True):
    """
    Generate h2m input from ClinVar data.  

    Parameter:
        - path (str): the path of clinvar renference vcf.gz data.
        - list_of_ids (list): the list of variation ids. If no value, the function would output all entries in the ClinVar data file.  
        - keep (bool): True: keep all the original columns in the dataframe/ False: keep the necesssary columns for h2m only. Default to False.  

    Output: 
        An input dataframe for h2m modeling.

    Example:   
        >>> filepath = '.../GrCh37_clinvar_20230923.vcf.gz'
        >>> variation_ids = [925574, 925434, 926695, 925707, 325626, 1191613, 308061, 361149, 1205375, 208043]
        >>> df = h2m.clinvar_reader(filepath, variation_ids)
    """

    def clinvar_VCF_translator(filepath, variation_ids=None):
        vcf_file = pysam.VariantFile(filepath)

        records = []
        for record in vcf_file.fetch():
            if variation_ids is None or int(record.id) in variation_ids:
                if record.alts is not None:
                    for alt in record.alts:
                        records.append({
                            'Hugo_Symbol': record.info['GENEINFO'].split(':')[0] if 'GENEINFO' in record.info else None,
                            'Chromosome': record.chrom,
                            'Start_Position': record.pos,
                            'End_Position': record.stop,
                            'Reference_Allele': record.ref,
                            'Tumor_Seq_Allele2': alt,
                            'Variant_Type': record.info['CLNVC'] if 'CLNVC' in record.info else None,
                            'Variation_ID': int(record.id),
                            'Allele_ID': record.info['ALLELEID'] if 'ALLELEID' in record.info else None,
                            'CLNSIG': record.info['CLNSIG'] if 'CLNSIG' in record.info else None,
                            'CLNHGVS': record.info['CLNHGVS'] if 'CLNHGVS' in record.info else None,
                            'CLNDN': record.info['CLNDN'] if 'CLNDN' in record.info else None,
                            'ID':record.id
                        })

        vcf_file.close()
        return pd.DataFrame(records)

    # Translate VCF to DataFrame
    df = clinvar_VCF_translator(path, list_of_ids)

    # Rename columns
    col_name_old = ['Hugo_Symbol', 'Start_Position', 'End_Position', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2']
    col_name_new = ['gene_name_h', 'start_h', 'end_h', 'type_h', 'ref_seq_h', 'alt_seq_h']
    df = df.rename(columns=dict(zip(col_name_old, col_name_new)))

    # Select columns if not keeping all
    if not keep:
        df = df.loc[:, col_name_new]

    # Drop duplicates and reset index
    df = df.drop_duplicates().reset_index(drop=True)
    
    # Add index column
    df['index'] = range(len(df))
    
    # Rearrange columns to put 'index' at the front
    columns = ['index'] + [col for col in df.columns if col != 'index']
    df = df[columns]
    df['format'] = 'ClinVar'
    df['ref_seq_h'] = df['ref_seq_h'].fillna('')
    df['alt_seq_h'] = df['alt_seq_h'].fillna('')
    return df