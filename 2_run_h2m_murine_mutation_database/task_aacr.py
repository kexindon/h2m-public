import h2m
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
df_run = pd.read_csv(Path(sys.argv[1]))
sample_name = str(sys.argv[2])

path_h_ref, path_m_ref = 'GCF_000001405.25_GRCh37.p13_genomic.fna.gz', 'GCF_000001635.27_GRCm39_genomic.fna.gz'
records_h, index_list_h = h2m.genome_loader(path_h_ref)
records_m, index_list_m  = h2m.genome_loader(path_m_ref)
path_h_anno, path_m_anno = 'gencode_v19_GRCh37.db', 'gencode_vm33_GRCm39.db'
db_h, db_m = h2m.anno_loader(path_h_anno), h2m.anno_loader(path_m_anno)

result_list = []
for f in [0,2,5,10]:
    result = h2m.model_batch(df_run, records_h, index_list_h, records_m, index_list_m, db_h, db_m, 37, flank_size = f, memory_size=2000)
    # save result
    result[0].to_csv(f'{sample_name}_result_f{f}.csv', index=False)
    # save failed samples
    result[1].to_csv(f'{sample_name}_f{f}_left_over.csv',index=False)