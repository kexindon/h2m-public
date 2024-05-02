import h2m
import pegg
import pandas as pd
import os
import sys
from pathlib import Path

#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 3:
    print("Usage: python3 task_aacr.py <input_df> <output_file_name>")
    sys.exit(1)

#usage
path_h_ref, path_m_ref = '/Users/kexindong/Documents/GitHub/Database/RefGenome/ncbi-2023-09-12/GCF_000001405.25_GRCh37.p13_genomic.fna.gz', '/Users/kexindong/Documents/GitHub/Database/RefGenome/mouse-2023-09-13/GCF_000001635.27_GRCm39_genomic.fna.gz'
chrom_dict_h, i = pegg.prime.genome_loader(path_h_ref)
chrom_dict_m, i = pegg.prime.genome_loader(path_m_ref)

df_all = df_all[df_all['Chromosome']!='M'].reset_index(drop=True)
df_all['Reference_Allele'] = [str(x) for x in df_all['Reference_Allele']]
df_all['Tumor_Seq_Allele2'] = [str(x) for x in df_all['Tumor_Seq_Allele2']]
df_all['Start_Position'] = [int(x) for x in df_all['Start_Position']]
df_all['End_Position'] = [int(x) for x in df_all['End_Position']]

df_all = pd.read_csv('/Users/kexindong/Documents/GitHub/Output/h2m_database/PEGG/pe_human.csv')
# save result
df_be_h = pegg.base.run_base(df_result_h, 'cBioPortal', chrom_dict_h, PAM='NGN',filtration='ABE+CBE', ideal_edit_window=[4, 8], auto_SNP_filter=True, proto_size=19,context_size=120, RE_sites=None, polyT_threshold=4, before_proto_context=5,sensor_length=40, sensor_orientation='reverse-complement', sensor=True)
df_be_h

result[0].to_csv(f'h2m_db_result_{sample_name}.csv', index=False)
# save failed samples
result[1].to_csv(f'h2m_db_left_over_{sample_name}.csv',index=False)