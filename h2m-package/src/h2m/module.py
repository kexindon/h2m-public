# import packages
import os
import json
import requests
import gffutils
import Bio.Align
from Bio import pairwise2
from Bio.Data import CodonTable
from itertools import product
import pandas as pd
import warnings
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pysam
from Bio import SeqIO
import gzip

warnings.filterwarnings("ignore")

h2m_dir = os.path.dirname(__file__)
json_path = os.path.join(h2m_dir, 'Data', 'human_name_dict.json')
with open(json_path, 'r') as json_file:
    name_dict = json.load(json_file)

json_path = os.path.join(h2m_dir, 'Data', 'homo_dict.json')
with open(json_path, 'r') as json_file:
    homo_dict = json.load(json_file)


json_path = os.path.join(h2m_dir, 'Data', 'homo_dict_mouse.json')
with open(json_path, 'r') as json_file:
    homo_dict_mouse = json.load(json_file)

json_path = os.path.join(h2m_dir, 'Data', 'dict_tx_h.json')
with open(json_path, 'r') as file:
    dict_tx_h = json.load(file)

json_path = os.path.join(h2m_dir, 'Data', 'dict_tx_m.json')
with open(json_path, 'r') as file:
    dict_tx_m = json.load(file)

# external data dictionary
## 1. refgenome 
## 2. genome annotation
## 3. nonstop_size
## 4. flank_size

# how to load?
def genome_loader(path):
    """
    Load the refernce genome file.

    Parameters:
        - path (str): path of the genome file.  

    Return:
        - reference genome records and the index list of chromosomes.

    Example:
        >>> records_h, index_list_h = h2m.genome_loader(path_h_ref)
    """
    #------loading in reference genome and organizing it into a 2-d list by chromosome---------------------

    with gzip.open(path, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta")) #about 4 Gb in  memory
        #records = list that contains sequences split up by chromosome (and intrachromosome splits up to some size)

    #filtering out alternative sequences to only select consensus chromosome sequences
    wrong = ["alternate", "unplaced", "unlocalized", "patch"]
    badlist = []
    for key in wrong:
        for i in records:
            ii = i.description
            if key in ii:
                badlist.append(ii)

    #creating an 
    filtered = []
    index_list = []
    for idx, i in enumerate(records):
        ii = i.description
        if ii not in badlist:
            filtered.append(ii)
            index_list.append(idx)
            
    return records, index_list

def anno_loader(path):
    """
    Load the GENCODE annotation file.

    Parameters:
        - path (str): path of the annotation file.  

    Return:
        - a FeatureDB

    Example:
        >>> db_h = h2m.anno_loader(path_h_anno)
    """
    return gffutils.FeatureDB(path)

# dict_parameter = {'non_stopsize': 300, 'flank_size' :2}

# server for different ref-genomes
server_dic = {37:'http://grch37.rest.ensembl.org/',
              38:'https://rest.ensembl.org/'}
species_dic = {'h':'homo_sapiens','m':'mouse'}

standard_table = CodonTable.unambiguous_dna_by_id[1]
be_tuple = [('A','G'),('T','C'),('C','T'),('G','A')]

list_order_of_condon = [
    'TTC', 'TTT', 'TTG', 'TTA', 'TAC', 'TAT', 'TGA', 'TAA', 'TAG',
    'CTG', 'CTC', 'CTT', 'CTA', 'CAC', 'CAT', 'CAG', 'CAA',
    'ATG', 'ATC', 'ATT', 'ATA', 'AAG', 'AAA', 'AAC', 'AAT',
    'GTG', 'GTC', 'GTT', 'GTA', 'GAG', 'GAC', 'GAT',  'GAA',
    'TCC', 'TCT', 'TCA', 'TCG', 'TGT', 'TGC', 'TGG',
    'CCT', 'CCC', 'CCA', 'CCG',
    'AGG', 'AGA', 'CGG', 'CGC', 'CGT', 'CGA',
    'ACC', 'ACA', 'ACT', 'ACG', 'AGC', 'AGT',
    'GCC', 'GCT', 'GCA', 'GCG', 'GGC', 'GGA', 'GGG', 'GGT'
]

list_order_of_coding_condon = [
    'TTC', 'TTT', 'TTG', 'TTA', 'TAC', 'TAT', 
    'CTG', 'CTC', 'CTT', 'CTA', 'CAC', 'CAT', 'CAG', 'CAA',
    'ATG', 'ATC', 'ATT', 'ATA', 'AAG', 'AAA', 'AAC', 'AAT',
    'GTG', 'GTC', 'GTT', 'GTA', 'GAG', 'GAC', 'GAT', 'GAA',
    'TCC', 'TCT', 'TCA', 'TCG', 'TGT', 'TGC', 'TGG',
    'CCT', 'CCC', 'CCA', 'CCG',
    'AGG', 'AGA', 'CGG', 'CGC', 'CGT', 'CGA',
    'ACC', 'ACA', 'ACT', 'ACG', 'AGC', 'AGT',
    'GCC', 'GCT', 'GCA', 'GCG', 'GGC', 'GGA', 'GGG', 'GGT'
]

#######query, Main function 1, local data based ########
######################################
def query(id, db = None , direction = 'h2m', ty = 'default',show = True): 
    """
    Query homologous mouse genes of human genes.

    Parameters:
        - id (str): name/gene_id/tx_id of human.
        - direction (str): OPTIONAL. query from human gene to the mouse gene ('h2m') or vise versa ('m2h').
        - db (FeatureDB): OPITONAL. The transcript annotation database of specific version.  
        - ty (str): OPITONAL. Specify the id type. one of 'gene_id'/'tx_id'/'name'.  
    
    Return:
        - a list of human gene name, mouse gene name, mapping type and sequence similarity.  

    Example:
        >>> h2m.query('TP53')
    """
    def get_what(lists, name):
        x = None
        for i in lists:
            if i[0] == name:
                x = i[1][0]
                break
        return x

    # - genome assembly version optional except for using tx_id for query
    if ty == 'default':
        if len(id)>4:
            if id[:4] == 'ENSG': # gene_id, 37 and 38 both applicable
                ty = 'gene_id'
            elif id[:4] == 'ENST': # transcript id
                ty = 'tx_id'
            else:
                ty = 'name'
        else:
            ty = 'name'

    if direction == 'h2m':
    ################################ Part1, unify the input format ################################ 
        # check the input type: name/gene_id/tx_id/...
        if ty == 'gene_id': # if user has specified the type as gene_id/ changed to gene_id above
            if '.' in id:
                # only keep the stable id
                id = id.split('.')[0]
            if id not in name_dict.keys(): # if the gene id is not included in the dictionary
                web_query = 'Error 1: This gene id is not included in the database. Please check the input format.'
                raise ValueError('Error 1: This gene id is not included in the database. Please check the input format.')
            else:
                name = name_dict[id]
        
        if ty == 'tx_id':
            if db is None:
                web_query = 'Error 2: Please include the transcript annotation database of specific version since you are querying with transcript id.'
                raise ValueError('Error 2: Please include the transcript annotation database of specific version since you are querying with transcript id.')
            lists = db[id].attributes.items()
            name = get_what(lists, 'gene_name')
            # name = db[id].attributes.items()[4][1][0]
        if ty == 'name':
            name = id
            if name not in name_dict.values():
                web_query = f'The query human gene: {name} has no mouse ortholog or this gene id is not included in the database. Please check the input format.'
                raise ValueError(f'The query human gene: {name} has no mouse ortholog or this gene id is not included in the database. Please check the input format.')

        ##############################################################################################
        ################################ Part2, check homology status ################################ 
        if name not in homo_dict.keys():
            output = [None]
            web_query = f'The query human gene: {name} has no mouse ortholog or this gene id is not included in the database. Please check the input format.'
        else:
            output = homo_dict[name]
            homo_type = list(output[0].values())[2]
            if homo_type is not None:
                homo_type = homo_type.split('_')[1]
            if len(output)>1:
                name_m = ','.join([list(output[x].values())[0] for x in range(len(output))])
                similarity = ', '.join([str(list(output[x].values())[3]) for x in range(len(output))])
            else:
                name_m, similarity = list(output[0].values())[0], str(list(output[0].values())[3])
            web_query = f'Query human gene: {name};\nMouse ortholog(s): {name_m};\nHomology type: {homo_type};\nSequence Simalarity(%):{similarity}.'
        if show:
            print(web_query)
        return output
    
    else:
        if ty != 'name':
            raise ValueError('Error 2: Please use gene name for m2h query.')
        else:
            output = homo_dict_mouse[id]
            homo_type = list(output[0].values())[1]
            if homo_type is not None:
                homo_type = homo_type.split('_')[1]
            if len(output)>1:
                name_h = ','.join([list(output[x].values())[0] for x in range(len(output))])
                similarity = ', '.join([str(list(output[x].values())[2]) for x in range(len(output))])
            else:
                name_h, similarity = list(output[0].values())[0], str(list(output[0].values())[2])
            web_query = f'Query human gene: {id};\nMouse ortholog(s): {name_h};\nHomology type: {homo_type};\nSequence Simalarity(%):{similarity}.'
        if show:
            print(web_query)
        return output


def query_batch(df, direction = 'h2m'):
    """
    Batch query of orthologous mouse gene of given human genes.
    
    Parameters:  
        - df (Pandas DataFrame): Must include a column of gene names named 'gene_name_h'. An index column is recommended.
        - direction: OPTIONAL. query from human gene to the mouse gene ('h2m') or vise versa ('m2h').
    Return:  
        - Two dataframes. The first dataframe is the processed original dataframe with canonical transcirpt id attached in the column named 'gene_name_m'. The second dataframe contains all rows that are not successfully processed.

    Example:
        >>> h2m.query_batch(df)
    """
    if direction == 'h2m':
        df_cor = df[[x in homo_dict.keys() for x in df['gene_name_h']]].reset_index(drop=True)
        df_wrong =  df[[x not in homo_dict.keys() for x in df['gene_name_h']]].reset_index(drop=True)
        if len(df_wrong) == 0:
            print('No error occurs.')
        else:
            print('There were rows that could not be processed.')

        list_of_gene = list(set(df_cor['gene_name_h']))
        list_homo = []
        for name_h in list_of_gene:
            list_m = [x['gene_name_m'] for x in homo_dict[name_h]]
            for name_m in list_m:
                list_homo.append({'gene_name_h':name_h,
                                    'gene_name_m':name_m})
        df_add = pd.DataFrame(list_homo)
        df_cor = pd.merge(df_cor,df_add,on='gene_name_h', how='outer')
        return([df_cor, df_wrong])
    else:
        df_cor = df[[x in homo_dict_mouse.keys() for x in df['gene_name_m']]].reset_index(drop=True)
        df_wrong =  df[[x not in homo_dict_mouse.keys() for x in df['gene_name_m']]].reset_index(drop=True)
        if len(df_wrong) == 0:
            print('No error occurs.')
        else:
            print('There were rows that could not be processed.')

        list_of_gene = list(set(df_cor['gene_name_m']))
        list_homo = []
        for name_m in list_of_gene:
            list_h = [x['gene_name_h'] for x in homo_dict_mouse[name_m]]
            for name_h in list_h:
                list_homo.append({'gene_name_m':name_m,
                                    'gene_name_h':name_h})
        df_add = pd.DataFrame(list_homo)
        df_cor = pd.merge(df_cor,df_add,on='gene_name_m', how='outer')
        return([df_cor, df_wrong])

#######get_tx_id, Main function 2, API based ########
######################################
def get_tx_id(id,  species, ver=None,ty = 'default', show = True):
    """ 
    Query a human or mouse gene for coordinate and information of all its transcripts. Internet needed.

    Parameters:  
        - id (str):, identification of a human gene. Multiple input forms are accepted, including gene name, stable ensembl gene id with or without version number.  
        - species (str): 'h' for human or 'm' for mouse.  
        - ver (int): specify the version of human, one of 37/38. It is a necessary parameter.  
        - ty (str): OPTIONAL. type of your input id. string, one of 'name'/'gene_id'.
        - show (bool): OPTIONAL. print summary of output or not.      

    Return:  
        - A list [chromosome, start location(of gene), end location(of gene), canonical transcript id, list of all transcript id (the canonical one included and always at the first place), a list of additional information of each transcript]  

    Example:
        >>> h2m.get_tx_id('TP53','h',ver=37)
    """
    ls = species_dic[species]
    if species not in ['h','m']:
        raise ValueError('Error 0: Please include species (h or m).')
    # choose server
    elif species == 'h':
        if (ver is None):
            raise ValueError('Error 1: Please include genome assembly version according to your transcript ID by adding ver=37 or ver=38.')
        elif ver not in [37,38]:
            raise ValueError('Error 1: Please include genome assembly version accfording to your transcript ID by adding ver=37 or ver=38.')
        else:
        # URL
            server = server_dic[ver]
    else:
        server = server_dic[38]
    
    if ty == 'default':
        if len(id)>4:
            if (id[:4] == 'ENSG') or (id[:4] == 'ENSM'): # gene_id, 37 and 38 both applicable
                ty = 'gene_id'
            else:
                ty = 'name'
        else:
            ty = 'name'
    
    if ty == 'name':
        ext = f'lookup/symbol/{ls}/{id}?expand=1'
    else:
        if '.' in id:
            # only keep the stable id
            id = id.split('.')[0]
        ext = f'lookup/id/{id}?expand=1'

    #output format
    #json (JavaScript Object Notation), the most a readable one 
    headers = {"Content-Type" : "application/json"} 
    #request
    response = requests.get(server+ext, headers=headers)
    # check the response status
    if not response.ok:
            response.raise_for_status()
            return None
    gene_annotation = response.json()
    # chr, main_tx
    chr, main_tx, start, end = gene_annotation['seq_region_name'], gene_annotation['canonical_transcript'], gene_annotation['start'], gene_annotation['end']

    if species == 'h': 
        if chr =='X':
            chr = 22
        elif chr =='Y':
            chr = 23
    elif species == 'm':
        if chr =='X':
            chr = 19
        elif chr =='Y':
            chr = 20
    chr = int(chr)
    # tx
    tx_dic_list = gene_annotation['Transcript']
    # put main_tx at the first location
    list_of_tx = []
    for tx_dic in tx_dic_list:
        head, tail = tx_dic['id'], tx_dic['version']
        list_of_tx.append(f'{head}.{tail}')

    loc = list_of_tx.index(main_tx)
    if (loc!=0):
        if loc==len(list_of_tx)-1:
            list_of_tx = ([list_of_tx[loc]]+list_of_tx[0:loc]).copy()
        else:
            list_of_tx = ([list_of_tx[loc]]+list_of_tx[0:loc] + list_of_tx[loc+1:]).copy()
    print_tx = ' '.join([''.join([f'({i+1})',x]) for i,x in enumerate(list_of_tx)])
    if show:
        if species == 'h':
            print(f'Genome assembly: GRCh{ver};\nThe canonical transcript is: {main_tx};\nYou can choose from the {len(list_of_tx)} transcripts below for further analysis:\n{print_tx}\n')
        else:
            print(f'Genome assembly: GRCm39;\nThe canonical transcript is: {main_tx};\nYou can choose from the {len(list_of_tx)} transcripts below for further analysis:\n{print_tx}\n')
    return([chr, start, end, main_tx, list_of_tx, tx_dic_list])

def get_tx_batch(df, species, ver=None):
    """
    Batch query of canonical transcript IDs of human or mouse genes.

    Parameters:  
        - df (Pandas DataFrame): Must include a column of gene names named 'gene_name_h'/'gene_name_m', depending on the species. An index column is recommended.
        - species (str): 'h' for human or 'm' for mouse.  
        - ver (int): specify the version of human, one of 37/38.

    Return:  
        - Two dataframes. The first dataframe is the processed original dataframe with canonical transcirpt id attached in the column named 'tx_id_h'/'tx_id_m'. The second dataframe contains all rows that are not successfully processed.

    Example:
        >>> h2m.get_tx_batch(df,'h',ver=37)
    """
    if species not in ['h','m']:
        raise ValueError('Error 0: Please include species (h or m).')
    # choose server
    elif species == 'h':
        if (ver is None):
            raise ValueError('Error 1: Please include genome assembly version according to your transcript ID by adding ver=37 or ver=38.')
        elif ver not in [37,38]:
            raise ValueError('Error 1: Please include genome assembly version accfording to your transcript ID by adding ver=37 or ver=38.')
        df_cor = df[[x in dict_tx_h.keys() for x in df['gene_name_h']]].reset_index(drop=True)
        df_wrong =  df[[x not in dict_tx_h.keys() for x in df['gene_name_h']]].reset_index(drop=True)
        if ver == 37:
            df_cor['tx_id_h'] = [dict_tx_h[x][0] for x in df_cor['gene_name_h']]
            df_cor['ref_genome_h'] = 'GRCh37'
        else:
            df_cor['tx_id_h'] = [dict_tx_h[x][1] for x in df_cor['gene_name_h']]
            df_cor['ref_genome_h'] = 'GRCh38'

        # should add an API here
    elif species == 'm':
        df_cor = df[[x in dict_tx_m.keys() for x in df['gene_name_m']]].reset_index(drop=True)
        df_wrong =  df[[x not in dict_tx_m.keys() for x in df['gene_name_m']]].reset_index(drop=True)

        df_cor['tx_id_m'] = [dict_tx_m[x] for x in df_cor['gene_name_m']]

    if len(df_wrong) == 0:
        print('No error occurs.')
    else:
        print('There were rows that could not be processed.')
    return([df_cor, df_wrong])
########## inner functions
# for model
def check_type(ty_h):
    if (ty_h is None):
        raise ValueError('Alarm 2: Please include mutation type (SNP/DNP/TNP/ONP/INS/DEL).')
    elif (ty_h not in ['SNP','DNP','TNP','ONP','INS','DEL']):
        raise ValueError('Alarm 2: Please include mutation type as SNP/DNP/ONP/INS/DEL.')

def get_non_coding(set_of_coding):
    coding = []
    if 'E' in set_of_coding:
        coding.append('Exon')
    if '5' in set_of_coding:
        coding.append("5'")
    if '3' in set_of_coding:
        coding.append("3'")
    if 'I' in set_of_coding:
        coding.append('Intron')
    return(','.join(coding))

def GateTx(range_idx, start_h, end_h):
    global status,error,statement
    area = range(start_h, end_h+1)
    # target human mutation is not in human CDS
    inter = set(range_idx).intersection(set(area))
    inclusion = len(inter)==len(area)
    if not inclusion:
        status = False
        statement = f'Class 5: Coordinate error. This mutation is not in the query gene.'
        error = 5

def GateFlankH(start_h, end_h, align_idx_h, align_idx_m, mut_tx_idx_h):
    global status, error, statement
    identity = set(list(range(start_h, end_h+1))).issubset(set(align_idx_h))
    def is_continuous(list_of_idx):
        if len(list_of_idx) == len(list(range(list_of_idx[0],list_of_idx[-1]+1))):
            return True
        else:
            return False
    if identity == False:
        status = False
        statement = 'Class 4: Mutated sequences are not identical.'
        error = 4
    else:
        mut_tx_idx_m_test = transfer_aligned_idx(mut_tx_idx_h, align_idx_h, align_idx_m)
        if is_continuous(mut_tx_idx_m_test) == False:
            status = False
            statement = 'Class 4: Mutated sequences are not identical.'
            error = 4

def get_flank(mut_tx_idx_h, align_idx_h, mut_tx_idx_m, align_idx_m,ty):
    def longest_continuous_both_sides(coords, target):
        target_left = target[0]
        target_right = target[-1]
        coord_set = set(coords)  # 将列表转换为集合
        longest_left, longest_right = 0,0

        longest_left = 0
        current = target_left - 1
        while current in coord_set:
            longest_left += 1
            current -= 1
        longest_right = 0

        current = target_right + 1
        while current in coord_set:
            longest_right += 1
            current += 1
        return longest_left, longest_right
    longest_left_h,longest_right_h = longest_continuous_both_sides(align_idx_h, mut_tx_idx_h)
    longest_left_m,longest_right_m = longest_continuous_both_sides(align_idx_m, mut_tx_idx_m)
    left, right = int(min(longest_left_h,longest_left_m)),int(min(longest_right_h,longest_right_m))
    return f'{left}{ty}', f'{right}{ty}'

def GateNonstop(new_p_seq_m, stop_loc_p_m):
    global status, error, statement
    if '*' in new_p_seq_m:
        if new_p_seq_m.index('*') == stop_loc_p_m:
            status, error, statement = False, 6, 'Class 6: This mutation cannot be originally modeled.'
        elif new_p_seq_m.index('*') < stop_loc_p_m:
            status, error, statement = False, 6, 'Class 6: This mutation cannot be originally modeled.'
        else:
            status, error, statement = True, 0, 'Class 0: This mutation can be originally modeled.'
    else:
        status, error, statement = True, 0, 'Class 0: This mutation can be originally modeled.'

def GateSilentstop(new_p_seq_m, stop_loc_p_m):
    global status, error, statement
    if '*' in new_p_seq_m:
        if new_p_seq_m.index('*') == stop_loc_p_m:
            status, error, statement = True, 0, 'Class 0: This mutation can be originally modeled.'
        else:
            status, error, statement = False, 6, 'Class 6: This mutation cannot be originally modeled.'
    else:
        status, error, statement = False, 6, 'Class 6: This mutation cannot be originally modeled.'

def GateNonsense(new_p_seq_m, loc):
    global status, error, statement
    if '*' in new_p_seq_m:
        if new_p_seq_m.index('*') == loc:
            status, error, statement = True, 0, 'Class 0: This mutation can be originally modeled.'
        else:
            status, error, statement = False, 6, 'Class 6: This mutation cannot be originally modeled.'
    else:
        status, error, statement =False, 6, 'Class 6: This mutation cannot be originally modeled.'

# for model
def transfer_aligned_idx(list_of_idx, tx_align_idx_from, tx_align_idx_to):
    aloc = [tx_align_idx_from.index(x) for x in list_of_idx]
    return([tx_align_idx_to[x] for x in aloc])

def transfer_stop_index(list_of_idx, stop_loc_tx_h, stop_loc_tx_m):
    mut_tx_idx_m = [x + stop_loc_tx_m - stop_loc_tx_h for x in list_of_idx]
    return mut_tx_idx_m

def GetAltM(strand_m, mut_dna_idx_m, alt_seq_m, ref_seq_m):
    if strand_m == "+":
        return([mut_dna_idx_m[0], mut_dna_idx_m[-1], alt_seq_m, ref_seq_m])
    if strand_m == "-":
        alt_r = str(Bio.Seq.Seq(alt_seq_m).reverse_complement())
        ref_r = str(Bio.Seq.Seq(ref_seq_m).reverse_complement())
        return([mut_dna_idx_m[-1], mut_dna_idx_m[0], alt_r, ref_r])

def GetAltH(strand_h, alt_seq_m, ref_seq_m):
    if strand_h == "+":
        return([alt_seq_m, ref_seq_m])
    if strand_h == "-":
        alt_r = str(Bio.Seq.Seq(alt_seq_m).reverse_complement())
        ref_r = str(Bio.Seq.Seq(ref_seq_m).reverse_complement())
        return([alt_r, ref_r])
#######sequencing and splicing, Main function set 3, mainly based on GENCODE local db and genome sequencing db########
############################################################
# 3.1, get sequence from chr
def GetChrSeq(records,index_list,chrom, species):
    if chrom[:3] == 'chr':
        chrom = chrom[3:]
    if species == 'h': 
        if chrom =='X':
            chrom = 22
        elif chrom =='Y':
            chrom = 23
        else:
            chrom = int(chrom)-1
        seq = records[index_list[chrom]].seq.upper()

    elif species == 'm':
        if chrom =='X':
            chrom = 19
        elif chrom =='Y':
            chrom = 20
        else:
            chrom = int(chrom)-1
        seq = records[index_list[chrom]].seq.upper()     
    return seq

def GetTx(records,index_list, db_input, tx_id, species, ver, nonstop_size):
    db = db_input[tx_id]
    cds = list(db_input.children(tx_id, order_by='+end', featuretype=['CDS']))
    def get_what(lists, name):
        x = None
        for i in lists:
            if i[0] == name:
                x = i[1][0]
                break
        return x
    gene_id,  gene_name = get_what(db.attributes.items(),'gene_id'), get_what(db.attributes.items(),'gene_name')
    gene_start, gene_end = db.start, db.end

    strand = db.strand
    chr = db.chrom
    chr_seq = GetChrSeq(records,index_list,chr,species)

    exonlist= [[i.start, i.end] for i in cds]
    intronlist = [[i.end +1, j.start-1] for i,j in zip(cds[0:-1],cds[1:])]
    genelist, utr_3_list, utr_5_list = [],[],[]

    if exonlist == []:
        if strand == '-':
            utr_5_list = []
            utr_3_list = []
            genelist = [gene_end, gene_start]
        else:
            utr_5_list = []
            utr_3_list = []
            genelist = [gene_start, gene_end]
    else:
        if strand == '-':
            utr_3_list = [min(gene_start,max(1,exonlist[0][0]-nonstop_size)),exonlist[0][0]-1]
            utr_5_list = [exonlist[-1][-1]+1, gene_end]
            genelist = [utr_3_list[0], utr_5_list[-1]]

        else:
            utr_5_list = [gene_start,exonlist[0][0]-1]
            utr_3_list = [exonlist[-1][-1]+1, max(min(exonlist[-1][-1]+nonstop_size, len(chr_seq)),gene_end)]
            genelist = [utr_5_list[0], utr_3_list[-1]]
        
    chr_index = list(range(genelist[0],genelist[1]+1)) # in the order of chromosome
    gene_seq = ''
    gene_seq += chr_seq[(genelist[0]-1):genelist[1]]

    if strand == '-':
        chr_index, exonlist, intronlist = chr_index[::-1], exonlist[::-1], intronlist[::-1] # in the order of stranded_chromosome
        gene_seq = Bio.Seq.Seq(gene_seq).reverse_complement()

    # annotation of intro, extro and nonstop region
    # get annotation of extron and intron number
    gene_region_anno, gene_region_idx = [None] * len(chr_index), [None] * len(chr_index)
    dict_map_dna_to_tx = dict(zip(chr_index, range(len(chr_index))))

    n_of_exon = len(exonlist)

    if exonlist == []:
        for l in range(genelist[0], genelist[1]):
            gene_region_anno[dict_map_dna_to_tx[l]] = 'I'
            gene_region_idx[dict_map_dna_to_tx[l]] = 'I'
            gene_region_idx[-1]
            stop_loc = gene_region_anno[::-1][:3]
            for l in stop_loc:
                gene_region_anno[dict_map_dna_to_tx[l]] = 'S'
                gene_region_idx[dict_map_dna_to_tx[l]] = -1
    else:
        for i,k in enumerate(exonlist):
            for l in range(k[0],k[1]+1):
                gene_region_anno[dict_map_dna_to_tx[l]] = 'E'
                gene_region_idx[dict_map_dna_to_tx[l]] = f'E_{i+1}'

        for i,k in enumerate(intronlist):
            for l in range(k[0],k[1]+1):
                gene_region_anno[dict_map_dna_to_tx[l]] = 'I'
                gene_region_idx[dict_map_dna_to_tx[l]] = f'I_{i+1}'

        for l in range(utr_5_list[0], utr_5_list[1]+1):
            gene_region_anno[dict_map_dna_to_tx[l]] = '5'
            gene_region_idx[dict_map_dna_to_tx[l]] = -5
        
        for l in range(utr_3_list[0], utr_3_list[1]+1):
            gene_region_anno[dict_map_dna_to_tx[l]] = '3'
            gene_region_idx[dict_map_dna_to_tx[l]] = -3

        if strand == '-':
            stop_loc = [exonlist[-1][0]-1, exonlist[-1][0]-2, exonlist[-1][0]-3]
        else:
            stop_loc = [exonlist[-1][-1]+1, exonlist[-1][-1]+2, exonlist[-1][-1]+3]
        
        for l in stop_loc:
            gene_region_anno[dict_map_dna_to_tx[l]] = 'S'
            gene_region_idx[dict_map_dna_to_tx[l]] = -1

    return([strand, gene_id, gene_name, chr, n_of_exon, gene_seq, chr_index, gene_region_anno, gene_region_idx, exonlist])

def Splicing(seq, Set):
    new_seq = Bio.Seq.Seq(''.join([seq[x] for x in range(len(seq)) if x in Set]))
    return(new_seq)

def Splicing_idx(idx, Set):
    new_index = [idx[x] for x in range(len(idx)) if x in Set]
    return(new_index)

def PrintAlign(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    print(pairwise2.format_alignment(*alignments[0]))

# for web, to generate visualized alignment, for both transcript and protein

def match_index(align1, align2, score, begin, end): 
    s = [] 
    s.append("%s\n" % align1) 
    s.append(" " * begin) 
    for a, b in zip(align1[begin:end], align2[begin:end]): 
        if a == b: 
            s.append("|")  # match 
        elif a == "-" or b == "-": 
            s.append(" ")  # gap 
        else: 
            s.append(".")  # mismatch 
    s.append("\n") 
    s.append("%s\n" % align2) 
    s.append("  Score=%g\n" % score)
    c = []
    for pos, char in enumerate(s):
        pipe = "|"
        if char == pipe:
            c.append(pos-2)
    return(c)

#######alignment, Main function set 4, mainly based on GENCODE local db and genome sequencing db########
############################################################
def align_idx(human_prot, mouse_prot, glo = True):
    "Returns indexes of HOMOLOGOUS regions in mouse and human"
    if glo:
        alignments = pairwise2.align.globalxx(human_prot, mouse_prot)
    else:
        alignments = pairwise2.align.localxx(human_prot, mouse_prot)

    match_idx = match_index(*alignments[0])
    
    human_aln_seq = alignments[0][0] #human alignment sequence
    mouse_aln_seq = alignments[0][1] #human alignment sequence

    human_idx = []
    mouse_idx = []
    
    for index in match_idx:
        #calculate human original idx
        gaps = human_aln_seq[:index].count('-')
        original_index = len(human_aln_seq[:index]) - gaps
        human_idx.append(original_index)
        
        #calculate mouse original idx
        gaps = mouse_aln_seq[:index].count('-')
        original_index = len(mouse_aln_seq[:index]) - gaps
        mouse_idx.append(original_index)

    return human_idx, mouse_idx
    
def GetAlign(tx_seq_1, tx_seq_2):
    alignments = pairwise2.align.globalxx(Translate(tx_seq_1), Translate(tx_seq_2))
    a_p_seq_1, a_p_seq_2 = alignments[0], alignments[1]

    # produce tx sequence
    match_idx = match_index(*alignments[0])
    human_aln_seq = alignments[0][0] #human alignment sequence
    mouse_aln_seq = alignments[0][1] #human alignment sequence
    human_idx = []
    mouse_idx = []
    for index in match_idx:
        #calculate human original idx
        gaps = human_aln_seq[:index].count('-')
        original_index = len(human_aln_seq[:index]) - gaps
        human_idx.append(original_index)
        
        #calculate mouse original idx
        gaps = mouse_aln_seq[:index].count('-')
        original_index = len(mouse_aln_seq[:index]) - gaps
        mouse_idx.append(original_index)
    return(human_idx, mouse_idx)

def get_type(ref,alt):
    """
    Get MAF style variant type.

    Parameters:
        - ref (str): reference seuquence in MAF/VCF style. 
        - alt (str): reference seuquence in MAF/VCF style. 
    """
    if ref == '-':
        r = 'INS'
    elif alt == '-':
        r = 'DEL'
    elif len(ref) == len(alt):
        if len(ref) ==1:
            r = 'SNP'
        elif len(ref) == 2:
            r = 'DNP'
        elif len(ref) ==3:
            r = 'TNP'
        elif len(ref) >=4:
            r = 'ONP'
        else:
            r = None
    elif len(ref)<len(alt):
        r = 'INS'
    else:
        r = 'DEL'
    return r

def ModifySeq(s,e,old_seq, alt_seq):
    # do not include s!
    # do not include e!
    new_seq = old_seq[:s] + alt_seq + old_seq[(e+1):]
    return(new_seq)

# for all of the 
def get_modify_seq(ty, mut_tx_idx, tx_seq, alt_seq):
    if ty in ['SNP','DNP','TNP','ONP']:
        new_seq = ModifySeq(mut_tx_idx[0], mut_tx_idx[-1], tx_seq, alt_seq)
    elif 'DEL' in ty:
        new_seq = ModifySeq(mut_tx_idx[0], mut_tx_idx[-1], tx_seq, alt_seq)
    elif 'INS' in ty:
        new_seq = ModifySeq(mut_tx_idx[0]+1, mut_tx_idx[-1]-1, tx_seq, alt_seq)
    return new_seq
    # for web
def get_tx_seq(p_seq, tx_seq):
    new_tx_seq = []
    for i, s in p_seq:
        if s == '-':
            new_tx_seq.append('---')
        else:
            new_tx_seq = new_tx_seq + [tx_seq[x] for x in p_to_tx_3([i])]
    return(new_tx_seq)
            
#a_tx_seq_1, a_tx_seq_2 = get_tx_seq(a_p_seq_1, tx_seq_1), get_tx_seq(a_p_seq_2, tx_seq_2)
#return a_tx_seq_1, a_p_seq_1, a_tx_seq_2, a_p_seq_2

#######Transcribe related, Main function set 5 ########
#######################################################
def Translate(seq):
    return Bio.Seq.Seq(seq).transcribe().translate()

def tx_to_p_idx(tx):
    return([x//3 for x in tx])

def tx_to_p_con(tx):
    return([x%3 for x in tx])

def p_to_tx_3(p):
    list_of_condn = []
    for i in p:
        list_of_condn.extend([3*i, 3*i+1, 3*i+2])
    return(list_of_condn)

def p_to_tx_con(idx, con):
    return([x*3+y for x,y in zip(idx,con)])

def GetCondon(seq, condon_list):
    return ''.join([seq[x] for x in p_to_tx_3(condon_list)])

def unique_ordered(lst):
    seen = set()
    result = []
    for item in lst:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result

def IfAlternative(ty_h):
    global status, error, statement
    if ty_h in ['INS','DEL']:
        status, error, statement = False, 2, 'Class 2: This mutation can be modeled, but the effect may not be consistent.'

def getCombiList(lists):
    all_combinations = list(product(*lists))
    combined_strings = [''.join(combination) for combination in all_combinations]
    return combined_strings

def GetConList(p):
    ll = []
    for x in list(p):
        if x == '*':
            sub_list_condon = standard_table.stop_codons
        else:
            sub_list_condon = StandardCondon(str(x))
            sub_list_condon = [x for x in list_order_of_condon if x in sub_list_condon]
        ll.append(sub_list_condon)
    return(getCombiList(ll))

def compare_lists(list1, list2):
    if len(list1) != len(list2):
        raise ValueError("Different Length")
    
    unequal_indices = []

    for i in range(len(list1)):
        if list1[i] != list2[i]:
            unequal_indices.append(i)

    return unequal_indices

def StandardCondon(a):
    list_of_condon = [codon for codon, aa in standard_table.forward_table.items() if aa == a]
    return list_of_condon

def all_characters_same(s):
    if not s:
        return True
    first_char = s[0]
    for char in s:
        if char != first_char:
            return False
    return True

def BE_check(a,b):
    result = False
    if ((all_characters_same(a)) and all_characters_same(b)):
        if len(a) == len(b):
            if (a[0],b[0]) in be_tuple:
                result = True
    return result

def GateBE(ref_seq_m, alt_seq_m):
    global status, error, statement
    if not BE_check(ref_seq_m,alt_seq_m):
        status, error, statement = False, 6, 'Class 6: This mutation cannot be originally modeled.'
def get_tx_change(tx_id,start,end,ref,alt,c=0,len=0, direct = '-', dist = 0):
    if c == '5':
        start_number, end_number = len-start, len-end
        if start==end:
            string = f'{tx_id}:c.-{start_number}{ref}>{alt}'
        else:
            string = f'{tx_id}:c.-{start_number}_{end_number}{ref}>{alt}'
    if c == '3':
        if start==end:
            string = f'{tx_id}:c.*{start+1}{ref}>{alt}'
        else:
            string = f'{tx_id}:c.*{start+1}_{end+1}{ref}>{alt}'
    if c == 'I':
        if start==end:
            string = f'{tx_id}:c.{len+1}{direct}{dist}{ref}>{alt}'
        else:
            if direct == '-':
                x = -1
            else:
                x = 1
            pair = [abs(x*dist+end-start), dist]
            string = f'{tx_id}:c.{len+1}{direct}{min(pair)}_{max(pair)}{ref}>{alt}'
    else:
        if start==end:
            string = f'{tx_id}:c.{start+1}{ref}>{alt}'
        else:
            string = f'{tx_id}:c.{start+1}_{end+1}{ref}>{alt}'
    return string

def longest_non_matching_substring(str1, str2):
    min_len = min(len(str1), len(str2))
    front_match = 0
    while (front_match <= min_len - 1) and (str1[front_match] == str2[front_match]):
        front_match += 1
    start1 = start2 = front_match

    end1 = len(str1) 
    end2 = len(str2) 
    back_match = 0
    while (back_match < (min_len - front_match)) and str1[-1 - back_match] == str2[-1 - back_match]:
        back_match += 1
        end1 -= 1
        end2 -= 1
    #len1, len2 = end1 - start1 + 1, end2 - start2 + 1
    if '*' in str2:
        if (end2-1>str2.index('*') and end1-1>str2.index('*')):
            end1 = end2 = str2.index('*')+1
    elif '*' in str1:
        if (end2-1>str1.index('*') and end1-1>str1.index('*')):
            end1 = end2 = str1.index('*')+1
    return start1, end1, start2, end2

def get_new(mut_tx_idx, ref_seq, alt_seq):

    s = unique_ordered(tx_to_p_idx(mut_tx_idx))
    if ((len(alt_seq)-len(ref_seq)) % 3) == 0:
        len_h = (len(alt_seq)-len(ref_seq)) // 3
    else:
        len_h = ((len(alt_seq)-len(ref_seq)) // 3) +1
    len_final = len(s) + len_h
    if len_final <= 0:
        s = []
    else:
        s = list(range(s[0], s[0]+len_final))
    return s
    #if (len(alt_seq) - len(ref_seq)) % 3 ==0:
    #    if mut_tx_idx[0] //3 != 0:
    #        lists = list(range(mut_tx_idx[0],(mut_tx_idx[-1]+len(alt_seq)-len(ref_seq)+1)))
    #       x = unique_ordered(tx_to_p_idx(lists))
    #else:
    #    lists = list(range(mut_tx_idx[0],(mut_tx_idx[-1]+len(alt_seq)-len(ref_seq)+1)))
    #    if len(lists) == 0:
    #        lists = [mut_tx_idx[0]]
    #    return(unique_ordered(tx_to_p_idx(lists)))

def get_cor(mut_p_idx, new_mut_p_idx, mut_p, new_mut_p):
    a,b,c,d = longest_non_matching_substring(mut_p, new_mut_p)
    return mut_p_idx[a:b], new_mut_p_idx[c:d], mut_p[a:b], new_mut_p[c:d]

def get_pro_ref(seq, idx):
    if len(seq)==1:
        str1 = f'{seq}{idx[0]+1}'
    elif len(seq)>1:
        str1 = f'{seq[0]}{idx[0]+1}_{seq[-1]}{idx[-1]+1}'
    return str1

def if_in_frame(ref_seq_h, alt_seq_h):
    if ((len(alt_seq_h) - len(ref_seq_h)) % 3 == 0):
        classification_h = 'In_Frame'
    else:
        classification_h = 'Frame_Shifting'
    return classification_h

# check classification 
    # non-sense, could find alternative
        # non-sense-stop: the same
    # missense, could find alternative
    # silent, could find alternative
        # non-sense-stop: the same
    
        # if_in_frame() (may be different for ins and del)

    # will not effect splicing!
    # frame-shifting: keep the same transcript as the non-stop one!
    # in-frame-complete or in-frame-disturb: keep the same complete transcript!
def first_mismatch_position(list1, list2):
    min_len = min(len(list1), len(list2))
    
    for i in range(min_len):
        if list1[i] != list2[i]:
            return i  

    return None

def type_snop(p_seq, new_p_seq, stop_loc_p, mut_p, mut_p_idx, cor_mut_p_idx, cor_new_mut_p_idx, cor_mut_p, cor_new_mut_p, category, ty, ref_seq, alt_seq, nonstop_size):

    def hgvsp_non_stop(p_seq, new_p_seq, stop_loc_p_original):
        p_ori, p_mut = p_seq[stop_loc_p_original], new_p_seq[stop_loc_p_original]
        if '*' not in new_p_seq:
            ext = ''.join(['>',str(nonstop_size//3)])
        else:
            ext = str(new_p_seq.index('*') - stop_loc_p_original)
            pro_change = ''.join([str(p_ori), str(stop_loc_p_original+1), str(p_mut), 'ext*', ext])
        return pro_change

    if category == 3:
        if ty in ['SNP','DNP','ONP','TNP']:
            if (len(cor_mut_p) == len(cor_new_mut_p) == 0):
                classification = 'Silent'
                pro_change = f'{get_pro_ref(mut_p, mut_p_idx)}delins{mut_p}'
            elif '*' in cor_new_mut_p:
                classification = 'Nonsense'
                if len(cor_mut_p)>0:
                    if cor_new_mut_p == '*':
                        pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}{cor_new_mut_p}'
                    else:
                        pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}delins{cor_new_mut_p}'
                else:
                    s = cor_mut_p_idx[0]-1
                    pro_change = f'{p_seq[s]}{s+1}_{p_seq[s+1]}{s+2}ins{cor_new_mut_p}'
            else:
                classification = 'Missense'
                if(len(cor_mut_p)==1):
                    pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}{cor_new_mut_p}'
                else:
                    pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}delins{cor_new_mut_p}'

        elif ty in ['INS','DEL']:
            classification = if_in_frame(ref_seq, alt_seq)
            if classification == 'Frame_Shifting':
                loc = first_mismatch_position(p_seq, new_p_seq)
                if '*' not in new_p_seq:
                    pro_change = f'{p_seq[loc]}{loc+1}{new_p_seq[loc]}fs*>{nonstop_size}'
                else:
                    stop_loc = new_p_seq.index('*')
                    if cor_new_mut_p_idx == []:
                        x = stop_loc - loc + 1
                    else:
                        x = stop_loc - cor_new_mut_p_idx[0]+1
                    if x == 1:
                        pro_change = f'{p_seq[loc]}{loc+1}{new_p_seq[loc]}'
                    else:
                        pro_change = f'{p_seq[loc]}{loc+1}{new_p_seq[loc]}fs*{x}'
            elif classification == 'In_Frame':
                if cor_new_mut_p == '': # del
                    s, e = cor_mut_p_idx[0], cor_mut_p_idx[-1]
                    if s == e:
                        pro_change = f'{p_seq[s]}{s+1}del'
                    else:
                        pro_change = f'{p_seq[s]}{s+1}_{p_seq[e]}{e+1}del'
                elif cor_mut_p == '': # ins
                    s = cor_new_mut_p_idx[0]-1
                    pro_change = f'{p_seq[s]}{s+1}_{p_seq[s+1]}{s+2}ins{cor_new_mut_p}'
                else:
                    pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}delins{cor_new_mut_p}'
                
    elif category == 2:
            if '*' not in new_p_seq:
                classification = 'Non_Stop'
                pro_change = hgvsp_non_stop(p_seq, new_p_seq, stop_loc_p)
            else:
                new_stop_loc_p = new_p_seq.index('*')
                if stop_loc_p < new_stop_loc_p:   
                    classification = 'Non_Stop'
                    pro_change = hgvsp_non_stop(p_seq, new_p_seq, stop_loc_p)
                elif stop_loc_p == new_stop_loc_p:
                    if cor_mut_p == cor_new_mut_p:
                        classification = 'Silent'
                        pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}delins{cor_new_mut_p}'
                    else:
                        classification = 'Missense'
                        pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}delins{cor_new_mut_p}'
                elif stop_loc_p > new_stop_loc_p:
                    classification = 'Nonsense' 
                    if len(cor_mut_p)>0:
                        pro_change = f'{get_pro_ref(cor_mut_p, cor_mut_p_idx)}delins{cor_new_mut_p}'
                    else:
                        s = cor_mut_p_idx[0]-1
                        pro_change = f'{p_seq[s]}{s}_{p_seq[s+1]}{s+1}ins{cor_new_mut_p}'
    return classification, pro_change

def get_seq(seq, idx):
    if len(idx) == 0:
        return('')
    else:
        return(seq[idx[0]:(idx[-1]+1)])

########## Main function 3, get align

def get_tx_and_align(records_h, index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, ver, nonstop_size = 300, align = True, memory_protect = True, memory_size = 10000):
    def get_splicing(category):
        def condition(anno_list, condition_list):
            return [x for x in range(len(anno_list)) if anno_list[x] in condition_list]
        if category == 2: # for mutations that cover the stop codon
            set_of_mut_idx_idx_h = set(condition(gene_region_anno_h, ['E','S','3']))
            tx_seq_h, tx_idx_h, region_idx_h = Splicing(gene_seq_h,set_of_mut_idx_idx_h), Splicing_idx(gene_index_h,set_of_mut_idx_idx_h), Splicing_idx(gene_region_idx_h,set_of_mut_idx_idx_h)
            set_of_mut_idx_idx_m = set(condition(gene_region_anno_m, ['E','S','3']))
            tx_seq_m, tx_idx_m, region_idx_m = Splicing(gene_seq_m,set_of_mut_idx_idx_m), Splicing_idx(gene_index_m,set_of_mut_idx_idx_m), Splicing_idx(gene_region_idx_m,set_of_mut_idx_idx_m)
            return [tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m]
    
    strand_h, gene_id_h, gene_name_h, chr_h, n_exon_h, gene_seq_h, gene_index_h, gene_region_anno_h, gene_region_idx_h, extron_list_h = GetTx(records_h, index_list_h, db_h, tx_id_h, 'h', ver, nonstop_size)
    strand_m, gene_id_m, gene_name_m, chr_m, n_exon_m, gene_seq_m, gene_index_m, gene_region_anno_m, gene_region_idx_m, extron_list_m  = GetTx(records_m, index_list_m, db_m, tx_id_m, 'm', ver, nonstop_size)

    tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m = get_splicing(category = 2)
    p_seq_h, p_seq_m = Translate(tx_seq_h), Translate(tx_seq_m)
    stop_loc_p_h, stop_loc_p_m = p_seq_h.index('*'), p_seq_m.index('*')
    if align:
        gate_memory(memory_protect, memory_size, stop_loc_p_h, stop_loc_p_m)
        align_p_idx_h, align_p_idx_m = align_idx(p_seq_h[:stop_loc_p_h+1], p_seq_m[:stop_loc_p_m+1])
        #align_tx_idx_h, align_tx_idx_m = align_idx(gene_seq_h, gene_seq_m)
        return [[strand_h, gene_id_h, gene_name_h, chr_h, n_exon_h, gene_seq_h, gene_index_h, gene_region_anno_h, gene_region_idx_h, extron_list_h],[strand_m, gene_id_m, gene_name_m, chr_m, n_exon_m, gene_seq_m, gene_index_m, gene_region_anno_m, gene_region_idx_m, extron_list_m],[tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m], [align_p_idx_h, align_p_idx_m]]
    else:
        return [[strand_h, gene_id_h, gene_name_h, chr_h, n_exon_h, gene_seq_h, gene_index_h, gene_region_anno_h, gene_region_idx_h, extron_list_h],[strand_m, gene_id_m, gene_name_m, chr_m, n_exon_m, gene_seq_m, gene_index_m, gene_region_anno_m, gene_region_idx_m, extron_list_m],[tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m], None]

def gate_memory(memory_protect, memory_size, len_1, len_2):
    if memory_protect:
        if (len_1>memory_size or len_2>memory_size):
            raise ValueError(f'Alarm 3: To perform this protection, you need to do the sequence alignment that involves at least one sequence longer than {memory_size}. Please make sure you have enough memory and then set memory_protect = False. Otherwise, your kernel may crush.')

# Main function 4, model
def model(records_h, index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, start, end, ref_seq, alt_seq, ty_h = None, ver = None, direction = 'h2m', param = 'default', coor = 'nc',search_alternative = True, max_alternative = 5, nonstop_size = 300, splicing_size = 30, batch = False, show_sequence = False, align_input = None, memory_protect = True, memory_size = 10000):
    """
    Model human variants in the mouse genome.  

    Parameters:
        - records_h, index_list_h, records_m, index_list_m: human and mouse reference genome.  
        - db_h, db_m: human and mouse GENCODE annotation.  
        - tx_id_h, tx_id_m: human and mouse transcript id (could get by **h2m.get_tx_id()**). Transcript ids of input and output variants if use direction = 'h2h' or  direction = 'm2m'.
        - start_h, end_h: int, start and end location of the mutation on the chromosome. 
        - ref_seq_h: str, human mutation reference sequence. Reference sequence of the input variant if use direction = 'm2h'/'h2h'/'m2m'.
        - alt_seq_h: str, human mutation alternate sequence. Alternate sequence of the input variant if use direction = 'm2h'/'h2h'/'m2m'.
        - ty_h: str, human variantion type in MAF format. One of ['SNP', 'DNP', 'TNP','ONP', 'INS', 'DEL'].  
        - ver: int, human ref genome number. 37 or 38.
        - direction (optional): str, set the modeling direction by 'h2m' (default) or 'm2h', 'h2h', 'm2m'.
        - param (optional): set param = 'BE' and will only output base editing modelable results.  
        - coor (optional): default = 'nc'. set input = 'aa' and will be compatable with input of amino acid variants.
        - search_alternative (optional): set search_alternative = False and will only output original modeling results.  
        - max_alternative (optional): the maximum number of output alternatives of one human variants.  
        - nonstop_size (optional): the length of neucleotides that are included after the stop codon for alignment and translation in case of the nonstop mutations or frame shifting mutations.  
        - splicing_size (optional): the number of amino acids or neucleotides (for non-coding mutations) that are included after the top codon for the consideration of frame-shifting effect.  
        - batch (optional): set batch = True and will use input align_dict to save time in batch processing.  
        - show_sequence (optional): set batch = True and will output the whole sequences.
        - align_dict (optional): input a prepared dictionary of alignment indexes to save time in batch processing.  
        - memory_protect (optional): default True. Break long alignments that may lead to death of the kernel.    
        - memory_size (optional): maxlength of aligned sequence when memory_protect == True.

    Other rules:
        1. If the mutation falls in the coding and non-coding regions at the same time, it would be considered and processed as a ORIGIAL-MODELING ONLY mutation.
        3. The alt_seq input should be in the positive strand and the start_h coordinate should be smaller than or equal the end_h coordinate.
        4. If the ref-seq or alt-see has no length, it could be input as '' or '-'.

    Example:
        >>> h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577120, 7577120, 'C', 'T', ty_h = 'SNP', ver = 37)
    """
    start_h, end_h, ref_seq_h, alt_seq_h = int(start), int(end), ref_seq, alt_seq
    global status, error, statement
    if (ver is None):
        raise ValueError('Alarm 1: Please include genome assembly version according to your transcript ID by adding ver=37 or ver=38.')
    elif ver not in [37,38]:
        raise ValueError('Alarm 1: Please include genome assembly version according to your transcript ID by adding ver=37 or ver=38.')
    status,error,statement = True, 0, "'Class 0: This mutation can be originally modeled." # initialization
    tx_change_h=tx_change_m_ori=tx_change_m=extron_list_h = extron_list_m = classification_m_ori = start_m_ori = end_m_ori = ref_seq_m_ori = alt_seq_m_ori = pro_change_m_ori = new_tx_seq_h = new_tx_seq_m_ori = new_tx_seq_m = exon_m = exon_h = ref_seq_m = category = start_m = end_m = pro_change_m = alt_seq_m = mut_tx_idx_h = strand_m = strand_h = gene_name_h = gene_id_h = chr_h = n_exon_h = n_exon_m = pro_change_h  = classification_h = classification_m = gene_name_m = gene_id_m = chr_m = n_exon_m  = ty_m = None
    mut_p_idx_h=new_mut_p_h=mut_p_idx_m= new_mut_p_m=new_mut_p_idx_h=new_mut_p_idx_m=None
    tx_seq_h = tx_seq_m = new_tx_seq_h = new_tx_seq_m_ori = new_tx_seq_m = mut_tx_idx_h = mut_tx_idx_m = dist_h_return = dist_m_return = dist_h = dist_m = None
    match_h=True
    flank_size_left = flank_size_right = None

    # check the intergrity of input
    def check_mut_info():
        nonlocal alt_seq_h, ref_seq_h, start_h, end_h, match_h

        if ref_seq_h == '-':
            ref_seq_h = ''
        ref_seq_h_input = ref_seq_h

        if alt_seq_h == '-':
            alt_seq_h = ''
        if ty_h == 'INS' and len(ref_seq_h) >= 1:
            start_h = start_h - 1
            end_h = end_h + 1

        if coor == 'aa':

            if alt_seq_h != '':
                alt_seq_h = GetConList(alt_seq_h)[0]
            
            if ty_h == 'INS':
                ref_seq_h = ''.join([gene_seq_h[gene_index_h.index(x)] for x in range(start_h+1, end_h)])
            else:
                ref_seq_h = ''.join([gene_seq_h[gene_index_h.index(x)] for x in range(start_h, end_h+1)])
            
            if strand_h == '-':
                ref_seq_h = ref_seq_h[::-1]
            
            if ref_seq_h_input is not None:
                if ref_seq_h == '':
                    if ref_seq_h_input != '':
                        match_h = False 
                if str(Translate(ref_seq_h)) != ref_seq_h_input:
                    match_h = False

        else:
            if ty_h == 'INS':
                ref_seq_h = ''.join([gene_seq_h[gene_index_h.index(x)] for x in range(start_h+1, end_h)])
            else:
                ref_seq_h = ''.join([gene_seq_h[gene_index_h.index(x)] for x in range(start_h, end_h+1)])
        
            # keep consistent with the transcript
            if strand_h == '-':
                ref_seq_h = ref_seq_h[::-1]
                alt_seq_h = str(Bio.Seq.Seq(alt_seq_h).reverse_complement())
                ref_seq_h_input = str(Bio.Seq.Seq(ref_seq_h_input).reverse_complement())

            if ref_seq_h_input is not None:
                if ref_seq_h != ref_seq_h_input:
                    match_h = False    
          
    check_type(ty_h)
    # with type, get the peptide effect of the human mutation. should be the same in the mouse genome, aka what we need to model.
    # get strand, cds, and sequencing
    if batch:
        if align_input is None:
            raise ValueError('Alarm 2: Please include align input first since you have selected batch processing.')
        else:
            result_1 = align_input[0]
            result_2 = align_input[1]
    else:
        if direction == 'h2h':
            result_1 = GetTx(records_h, index_list_h, db_h, tx_id_h, 'h', ver, nonstop_size)
            result_2 = GetTx(records_m, index_list_m, db_m, tx_id_m, 'h', ver, nonstop_size)
        elif direction == 'm2m':
            result_1 = GetTx(records_h, index_list_h, db_h, tx_id_h, 'm', ver, nonstop_size)
            result_2 = GetTx(records_m, index_list_m, db_m, tx_id_m, 'm', ver, nonstop_size)
        else:
            result_1 = GetTx(records_h, index_list_h, db_h, tx_id_h, 'h', ver, nonstop_size)
            result_2 = GetTx(records_m, index_list_m, db_m, tx_id_m, 'm', ver, nonstop_size)

    if direction in ['h2m','h2h','m2m']:
        strand_h, gene_id_h, gene_name_h, chr_h, n_exon_h, gene_seq_h, gene_index_h, gene_region_anno_h, gene_region_idx_h, extron_list_h = result_1
        strand_m, gene_id_m, gene_name_m, chr_m, n_exon_m, gene_seq_m, gene_index_m, gene_region_anno_m, gene_region_idx_m, extron_list_m  = result_2
    else:
        strand_h, gene_id_h, gene_name_h, chr_h, n_exon_h, gene_seq_h, gene_index_h, gene_region_anno_h, gene_region_idx_h, extron_list_h = result_2
        strand_m, gene_id_m, gene_name_m, chr_m, n_exon_m, gene_seq_m, gene_index_m, gene_region_anno_m, gene_region_idx_m, extron_list_m  = result_1
        temp = tx_id_h
        tx_id_h = tx_id_m
        tx_id_m = temp
    
    ######################
    # inner functions for splicing (transcript selection)
    def replace_keys(d, direct):
        if direct == 'h2m':
            return d
        elif direct == 'm2h':
            new_dict = {}
            for key, value in d.items():
                if key.startswith('human'):
                    new_key = 'mouse' + key[5:]  
                elif key.startswith('mouse'):
                    new_key = 'human' + key[5:]  
                else:
                    new_key = key

                if new_key.endswith('m'):
                    new_key = new_key[:-1] + 'h'
                elif new_key.endswith('h'):
                    new_key = new_key[:-1] + 'm' 

                if new_key.endswith('m_ori'):
                    new_key = new_key[:-5] + 'h_ori'
                elif new_key.endswith('h_ori'):
                    new_key = new_key[:-5] + 'm_ori' 

                new_dict[new_key] = value
            return new_dict
        elif direct == 'h2h':
            new_dict = {}
            for key, value in d.items():
                if key.startswith('human'):
                    new_key = 'human_1' + key[5:]  
                elif key.startswith('mouse'):
                    new_key = 'human_2' + key[5:]  
                else:
                    new_key = key

                if new_key.endswith('m'):
                    new_key = new_key[:-1] + 'h_2'
                elif new_key.endswith('h'):
                    new_key = new_key[:-1] + 'h_1' 

                if new_key.endswith('m_ori'):
                    new_key = new_key[:-5] + 'h_2_ori'
                elif new_key.endswith('h_ori'):
                    new_key = new_key[:-5] + 'h_1_ori' 
                new_dict[new_key] = value
            return new_dict
        elif direct == 'm2m':
            new_dict = {}
            for key, value in d.items():
                if key.startswith('human'):
                    new_key = 'mouse_1' + key[5:]  
                elif key.startswith('mouse'):
                    new_key = 'mouse_2' + key[5:]  
                else:
                    new_key = key

                if new_key.endswith('m'):
                    new_key = new_key[:-1] + 'm_2'
                elif new_key.endswith('h'):
                    new_key = new_key[:-1] + 'm_1' 

                if new_key.endswith('m_ori'):
                    new_key = new_key[:-5] + 'm_2_ori'
                elif new_key.endswith('h_ori'):
                    new_key = new_key[:-5] + 'm_1_ori' 
                new_dict[new_key] = value
            return new_dict

    def final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m, flank_size_left, flank_size_right, dist_h_return = None, dist_m_return = None):
        # nonlocal start_m, end_m, alt_seq_m, ref_seq_m, alt_seq_h, ref_seq_h, alt_seq_m_ori, ref_seq_m_ori, mut_dna_idx_m_ori, alt_seq_m_ori, ref_seq_m_ori, start_m_ori, end_m_ori
        alt_seq_h, ref_seq_h = GetAltH(strand_h, alt_seq_h, ref_seq_h)

        if alt_seq_h == '':
            alt_seq_h = '-'
            
        if ref_seq_h == '':
            ref_seq_h = '-'

        if (ty_h == 'INS') and (ref_seq_h != '-'):
                start_h = start_h+1
                end_h = end_h-1

        if status:
            start_m, end_m, alt_seq_m, ref_seq_m = GetAltM(strand_m, mut_dna_idx_m, alt_seq_m, ref_seq_m)
            start_m_ori, end_m_ori, alt_seq_m_ori, ref_seq_m_ori = GetAltM(strand_m, mut_dna_idx_m_ori, alt_seq_m_ori, ref_seq_m_ori)

            if alt_seq_m == '':
                alt_seq_m = '-'
            
            if ref_seq_m == '':
                ref_seq_m = '-'

            if alt_seq_m_ori == '':
                alt_seq_m_ori = '-'
            
            if ref_seq_m_ori == '':
                ref_seq_m_ori = '-'

            if (ty_m == 'INS') and (ref_seq_m != '-'):
                start_m = start_m+1
                end_m = end_m-1
                start_m_ori = start_m_ori+1
                end_m_ori = end_m_ori-1

        else:
            ty_m = classification_m = exon_m = start_m_ori = end_m_ori = ref_seq_m_ori = alt_seq_m_ori = pro_change_m = start_m = end_m = ref_seq_m = alt_seq_m = pro_change_m_ori = tx_change_m_ori = tx_change_m =  None

        if classification_h is not None:
            classification_h = str(classification_h)
        if classification_m is not None:
            classification_m = str(classification_m)
            
        return_dict = {'gene_name_h':gene_name_h,'gene_id_h': gene_id_h, 'tx_id_h': tx_id_h,'chr_h': chr_h, 'exon_num_h':n_exon_h, 'strand_h':strand_h,# human gene
                        'match':match_h,
                        'start_h': start_h, 'end_h': end_h, 'ref_seq_h': ref_seq_h, 'alt_seq_h':alt_seq_h, 'HGVSc_h':tx_change_h,'HGVSp_h': pro_change_h, # human mutation
                        'classification_h': classification_h, 'exon_h': exon_h, 'type_h': ty_h,    # human mutation
                        'status':status, 'class':error,'statement':statement,
                        'flank_size_left':flank_size_left, 'flank_size_right':flank_size_right, # status
                        'gene_name_m':gene_name_m,'gene_id_m': gene_id_m,'tx_id_m': tx_id_m,'chr_m': chr_m, 'exon_num_m':n_exon_m, 'strand_m':strand_m, # mouse gene
                        'type_m':ty_m, 'classification_m': classification_m, 'exon_m':exon_m,                                                           # mouse mutation
                        'start_m_ori': start_m_ori, 'end_m_ori': end_m_ori, 'ref_seq_m_ori':ref_seq_m_ori, 'alt_seq_m_ori': alt_seq_m_ori,'HGVSc_m_ori':tx_change_m_ori, 'HGVSp_m_ori': pro_change_m_ori,  # mouse mutation ori                                                                                                    
                        'start_m': start_m, 'end_m': end_m, 'ref_seq_m':ref_seq_m, 'alt_seq_m': alt_seq_m, 'HGVSc_m':tx_change_m,'HGVSp_m': pro_change_m}                             # mouse mutation alt
        if show_sequence:
            #  
            additional_dict = {
                # human
                'seq_h':str(tx_seq_h),
                'new_seq_h': str(new_tx_seq_h),
                'human_tx_idx':mut_tx_idx_h,
                'human_p_idx':mut_p_idx_h,
                'human_new_p_idx':new_mut_p_idx_h,
                'dist_h':dist_h_return,
                # mouse
                'seq_m':str(tx_seq_m),
                'new_seq_m_ori':str(new_tx_seq_m_ori),
                'mouse_tx_idx_ori':mut_tx_idx_m,
                'mouse_p_idx_ori':mut_p_idx_m,
                'mouse_new_p_idx_ori':new_mut_p_idx_m,
                'dist_m':dist_m_return,
                # mouse-alt
                'new_seq_m': str(new_tx_seq_m),
                'mouse_tx_idx':mut_tx_idx_m,
                'mouse_p_idx':mut_p_idx_m,
                'mouse_new_p_idx':new_mut_p_idx_m
            }
            return_dict.update(additional_dict)

        # change the appendix of h and m depending on the direction
        return [replace_keys(return_dict, direction).copy()]
    
    def get_splicing(category):
        def condition(anno_list, condition_list):
            return [x for x in range(len(anno_list)) if anno_list[x] in condition_list]
        if category == 4:
            # for 5' and 3', not influenced by # of exons
            # get segments that are covered by the mutated region
            set_of_mut_idx = set([gene_region_idx_h[gene_index_h.index(x)] for x in range(start_h, end_h+1)])# get the index (distinct for each segment)
            set_of_mut_idx_idx_h = set(condition(gene_region_idx_h, set_of_mut_idx))
            tx_seq_h, tx_idx_h = Splicing(gene_seq_h,set_of_mut_idx_idx_h), Splicing_idx(gene_index_h,set_of_mut_idx_idx_h)
            set_of_mut_idx_idx_m = set(condition(gene_region_idx_m, set_of_mut_idx))
            tx_seq_m, tx_idx_m = Splicing(gene_seq_m,set_of_mut_idx_idx_m), Splicing_idx(gene_index_m,set_of_mut_idx_idx_m)
            tx_region_anno_m = Splicing_idx(gene_region_anno_m,set_of_mut_idx_idx_m)
            return [tx_seq_h, tx_idx_h, tx_seq_m, tx_idx_m, tx_region_anno_m]
        
        elif category == 1:
            set_of_mut_idx = set(range(max(min(gene_index_h),start_h-splicing_size), min(max(gene_index_h),end_h+splicing_size+1))) # 'I_13'
            # get the start and end location of the expanded region
            temp_mut_h_s, temp_mut_h_e = min(set_of_mut_idx), max(set_of_mut_idx)
            exons_h_distance = [abs(start-start_h) for start,_ in extron_list_h]
            # This part means to make a shortest alignment segment in order to save time. Expand the human non-coding region to a 300-nt-long segments, 
            # and the attach it to the nearest exon for alignment.
            # for the mouse, if they have the same numbers of exons, then do the same thing. Otherwise, keep the whole seuqence.
            ##################
            # dist_min_idx_h
            ##################
            dist_min_idx_h = exons_h_distance.index(min(exons_h_distance)) # get the closest exon number
            temp_h_s, temp_h_e = gene_index_h.index(extron_list_h[dist_min_idx_h][0]), gene_index_h.index(extron_list_h[dist_min_idx_h][1])

            temp_mut_h_s, temp_mut_h_e = gene_index_h.index(temp_mut_h_s), gene_index_h.index(temp_mut_h_e)
            set_of_mut_idx_idx_h = set(range(min(temp_mut_h_s, temp_mut_h_e, temp_h_s, temp_h_e), max(temp_mut_h_s, temp_mut_h_e, temp_h_s, temp_h_e)+1))
            # set_of_mut_idx_idx_h = set(range(min(temp_mut_h_s, temp_mut_h_e, temp_h_s, temp_h_e), max(temp_mut_h_s, temp_mut_h_e, temp_h_s, temp_h_e)+1))
            tx_seq_h, tx_idx_h = Splicing(gene_seq_h,set_of_mut_idx_idx_h), Splicing_idx(gene_index_h,set_of_mut_idx_idx_h)
            tx_region_anno_h_inner = Splicing_idx(gene_region_anno_h,set_of_mut_idx_idx_h)
            region_idx_h = Splicing_idx(gene_region_idx_h,set_of_mut_idx_idx_h)
            if n_exon_h == n_exon_m:
                # get segments that are covered by the mutated region
                if strand_m == '+':
                    dist_min_idx_m = dist_min_idx_h
                else: 
                    dist_min_idx_m = n_exon_m - dist_min_idx_h - 1

                temp_mut_m_s, temp_mut_m_e = gene_index_m.index(extron_list_m[dist_min_idx_m][0]), gene_index_m.index(extron_list_m[dist_min_idx_m][1])
                min_m, max_m = min(temp_mut_m_s, temp_mut_m_e), max(temp_mut_m_s, temp_mut_m_e)

                len_m = extron_list_m[dist_min_idx_m][1] - extron_list_m[dist_min_idx_m][0]
                len_h = extron_list_h[dist_min_idx_h][1] - extron_list_h[dist_min_idx_h][0]

                if max(set_of_mut_idx_idx_h) in [temp_h_s, temp_h_e]: # end as exon
                    set_of_mut_idx_idx_m = set(range(max_m-(max(set_of_mut_idx_idx_h)-min(set_of_mut_idx_idx_h))-len_m+len_h,max_m+1))
                else:
                    set_of_mut_idx_idx_m = set(range(min_m, min_m+1+(max(set_of_mut_idx_idx_h)-min(set_of_mut_idx_idx_h)+len_m-len_h)))
                tx_seq_m, tx_idx_m = Splicing(gene_seq_m,set_of_mut_idx_idx_m), Splicing_idx(gene_index_m,set_of_mut_idx_idx_m)
                tx_region_anno_m_inner = Splicing_idx(gene_region_anno_m,set_of_mut_idx_idx_m)
                region_idx_h_m = Splicing_idx(gene_region_idx_m,set_of_mut_idx_idx_m)
            else:
                tx_seq_m, tx_idx_m = gene_seq_m, gene_index_m
                tx_region_anno_m_inner = gene_region_anno_m
                region_idx_h = gene_region_idx_h
            return [tx_seq_h, tx_idx_h, tx_region_anno_h_inner, region_idx_h, tx_seq_m, tx_idx_m, tx_region_anno_m_inner, region_idx_h_m]
        
        elif category == 2: # for mutations that cover the stop codon, and coding-only
            if batch:
                if direction in ['h2m','h2h','m2m']:
                    return align_input[2]
                else:
                    list_of_output = align_input[2][3:]+align_input[2][:3]
                    return list_of_output
            else:
                set_of_mut_idx_idx_h = set(condition(gene_region_anno_h, ['E','S','3']))
                tx_seq_h, tx_idx_h, region_idx_h = Splicing(gene_seq_h,set_of_mut_idx_idx_h), Splicing_idx(gene_index_h,set_of_mut_idx_idx_h), Splicing_idx(gene_region_idx_h,set_of_mut_idx_idx_h)
                set_of_mut_idx_idx_m = set(condition(gene_region_anno_m, ['E','S','3']))
                tx_seq_m, tx_idx_m, region_idx_m = Splicing(gene_seq_m,set_of_mut_idx_idx_m), Splicing_idx(gene_index_m,set_of_mut_idx_idx_m), Splicing_idx(gene_region_idx_m,set_of_mut_idx_idx_m)
                return [tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m]

    def FindAlternative(max_alternative, mut_tx_idx_m_ori, gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h,mut_p_idx_m,new_mut_p_idx_m, mut_p_idx_h, new_mut_p_idx_h, align_p_idx_h, align_p_idx_m):
        global status, error, statement
        new_mut_p_idx_m_ori = new_mut_p_idx_m
        mut_p_idx_m_ori = mut_p_idx_m
        def snop(ref_seq, alt_seq):
            delta_len = len(ref_seq) - len(alt_seq)
            if delta_len>0:
                ty_m = 'DEL'
            elif delta_len<0:
                ty_m = 'INS'
            else:
                if len(ref_seq) == 1:
                    ty_m = 'SNP'
                elif len(ref_seq) == 2:
                    ty_m = 'DNP'
                elif len(ref_seq) == 3:
                    ty_m = 'TNP'
                elif len(ref_seq) >=4:
                    ty_m = 'ONP'
            return ty_m
        #, condon_m_mut, p_m_mut, mut_tx_idx_m, ref_seq_m, alt_seq_m, mut_p_idx_m, mut_p_idx_m_short, condon_m_ori
        if classification_h == 'Non_Stop':
            list_condon = list_order_of_coding_condon
            condon_m_ori = GetCondon(tx_seq_m, [stop_loc_p_m])
            mut_p_idx_m = [stop_loc_p_m]
            p_m_ori = '*'
        else:
            mut_p_idx_m = transfer_aligned_idx(cor_mut_p_idx_h, align_p_idx_h, align_p_idx_m)
            mut_p_m = ''.join(p_seq_m[x] for x in mut_p_idx_m)
            list_condon = GetConList(cor_new_mut_p_h)
            condon_m_ori = GetCondon(tx_seq_m, mut_p_idx_m)
            p_m_ori = Translate(condon_m_ori)

        original_condon_index = p_to_tx_3(mut_p_idx_m)
        # numbers of different base pairs for all the condons
        loc_list = [longest_non_matching_substring(std,str(condon_m_ori)) for std in list_condon]
        loc_list_num = [loc_list[x][1] - loc_list[x][0] for x in range(len(loc_list))]
        list_return = []
        mut_tx_idx_m = mut_tx_idx_m_ori
        for num in range(0, max(loc_list_num)):
            num = num + 1
            if num in loc_list_num:
                indexes = [index for index, element in enumerate(loc_list_num) if element == num]
                for idx in indexes:
                    if len(list_return) == max_alternative:
                        break
                    loc, std = loc_list[idx], list_condon[idx]
                    p_m_mut = Bio.Seq.Seq(std).transcribe().translate()
                    loc_alt = list(range(loc[0],loc[1]))
                    loc_ref = list(range(loc[2],loc[3]))
                    #mut_p_idx_m_short = list(set(mut_p_idx_m))
                    ref_seq_m = ''.join([condon_m_ori[x] for x in loc_ref])
                    alt_seq_m = ''.join([std[x] for x in loc_alt])
                    if loc_ref != []: # for ins
                        mut_tx_idx_m = [original_condon_index[x] for x in loc_ref]
                    else:
                        mut_tx_idx_m = mut_tx_idx_m_ori
                    mut_p_idx_m = tx_to_p_idx(mut_tx_idx_m)
                    flank_size_left, flank_size_right = get_flank(mut_p_idx_h, align_p_idx_h, mut_p_idx_m, align_p_idx_m,'aa')
                    # BE
                    check = False
                    if param is not None:
                        if param == 'BE':
                            if BE_check(ref_seq_m, alt_seq_m):
                                check = True
                        else:
                            check = True
                    else:
                        check = True
                    
                    if check:
                        status, error, statement = True, 1, 'Class 1: This mutation can be alternatively modeled.'
                        mut_dna_idx_m = [tx_idx_m[x] for x in mut_tx_idx_m]
                        classification_m = classification_h
                        ty_m = snop(ref_seq_m, alt_seq_m)
                        new_tx_seq_m = get_modify_seq(ty_h, mut_tx_idx_m, tx_seq_m, alt_seq_m)
                        new_p_seq_m = Translate(new_tx_seq_m)
                        mut_p_m, new_mut_p_m = get_seq(p_seq_m,mut_p_idx_m), get_seq(new_p_seq_m,new_mut_p_idx_m)
                        new_mut_p_idx_m  = get_new(mut_tx_idx_m, ref_seq_m, alt_seq_m)
                        tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m,alt_seq_m)
                        classification_m, pro_change_m = type_snop(p_seq_m, new_p_seq_m, stop_loc_p_m, mut_p_m, mut_p_idx_m, mut_p_idx_m, cor_new_mut_p_idx_m, p_m_ori, p_m_mut, category, ty_m, ref_seq_m, alt_seq_m, nonstop_size)
                        start_m, end_m, alt_seq_m, ref_seq_m = GetAltM(strand_m, mut_dna_idx_m, alt_seq_m, ref_seq_m)
                        alt_seq_h_cor, ref_seq_h_cor = GetAltH(strand_h, alt_seq_h, ref_seq_h)
                        start_m_ori, end_m_ori, alt_seq_m_ori_cor, ref_seq_m_ori_cor = GetAltM(strand_m, mut_dna_idx_m_ori, alt_seq_m_ori, ref_seq_m_ori)
                        return_dict = {'gene_name_h':gene_name_h,'gene_id_h': gene_id_h, 'tx_id_h': tx_id_h,'chr_h': chr_h, 'exon_num_h':n_exon_h, 'strand_h':strand_h,# human gene
                                        'match':match_h,
                                        'start_h': start_h, 'end_h': end_h, 'ref_seq_h': ref_seq_h_cor, 'alt_seq_h':alt_seq_h_cor, 'HGVSc_h':tx_change_h,'HGVSp_h': pro_change_h, # human mutation
                                        'classification_h': classification_h, 'exon_h': exon_h, 'type_h': ty_h,    # human mutation
                                        'status':status, 'class':error,'statement':statement,
                                        'flank_size_left':flank_size_left, 'flank_size_right':flank_size_right, # status                                                                                                                                                       # status
                                        'gene_name_m':gene_name_m,'gene_id_m': gene_id_m,'tx_id_m': tx_id_m,'chr_m': chr_m, 'exon_num_m':n_exon_m, 'strand_m':strand_m, # mouse gene
                                        'type_m':ty_m, 'classification_m': classification_m, 'exon_m':exon_m,                                                           # mouse mutation
                                        'start_m_ori': start_m_ori, 'end_m_ori': end_m_ori, 'ref_seq_m_ori':ref_seq_m_ori_cor, 'alt_seq_m_ori': alt_seq_m_ori_cor,'HGVSc_m_ori':tx_change_m_ori, 'HGVSp_m_ori': pro_change_m_ori,  # mouse mutation ori                                                                                                    
                                        'start_m': start_m, 'end_m': end_m, 'ref_seq_m':ref_seq_m, 'alt_seq_m': alt_seq_m, 'HGVSc_m':tx_change_m,'HGVSp_m': pro_change_m}                    # mouse mutation alt
                        if show_sequence:
                            additional_dict = {
                                    # human
                                    'seq_h':str(tx_seq_h),
                                    'new_seq_h': str(new_tx_seq_h),
                                    'human_tx_idx':mut_tx_idx_h,
                                    'human_p_idx':mut_p_idx_h,
                                    'human_new_p_idx':new_mut_p_idx_h,

                                    # mouse
                                    'seq_m':str(tx_seq_m),
                                    'new_seq_m_ori':str(new_tx_seq_m_ori),
                                    'mouse_tx_idx_ori':mut_tx_idx_m_ori,
                                    'mouse_p_idx_ori':mut_p_idx_m_ori,
                                    'mouse_new_p_idx_ori':new_mut_p_idx_m_ori,

                                    # mouse-alt
                                    'new_seq_m': str(new_tx_seq_m),
                                    'mouse_tx_idx':mut_tx_idx_m,
                                    'mouse_p_idx':mut_p_idx_m,
                                    'mouse_new_p_idx':new_mut_p_idx_m
                                }
                            return_dict.update(additional_dict)

                        return_dict = replace_keys(return_dict,direction).copy()

                        list_return.append(return_dict)

                if len(list_return) == max_alternative:
                        break
        if len(list_return)>=1:
            return list_return
        else:
            status, error, statement = False, 3, 'Class 3: This mutation cannot be originally modeled and no alternative is found.'
            x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m,  flank_size_left, flank_size_right)
            return x
    ######################
    # transcript selection
    # check coordinate
    if coor == 'nc':
        GateTx(gene_index_h, start_h, end_h)
        # strand reversed
        if not status:
            x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m, flank_size_left, flank_size_right)
            return x
        else:
            check_mut_info()
            set_of_mut = set([gene_region_anno_h[gene_index_h.index(x)] for x in range(start_h, end_h+1)])
            # if not cds only
            if 'S' in set_of_mut:
                category = 2
            elif (('5' in set_of_mut) or ('3' in set_of_mut)):        
                category = 4
            elif set_of_mut != {'E'}: #intron and splicing
                category = 1
            elif set_of_mut == {'E'}:
                category = 3

            ################################################category 4: 5'3' non-coding region ################################################
            if  category ==4: # 5' 3'
                t = ','.join(set_of_mut)
                classification_h = f'Non_Coding,{t}'
                tx_seq_h, tx_idx_h, tx_seq_m, tx_idx_m, tx_region_anno_m = get_splicing(category = category)
                gate_memory(memory_protect=memory_protect, memory_size=memory_size, len_1 = len(tx_seq_h), len_2 = len(tx_seq_m))
                if n_exon_h == n_exon_m:
                    align_idx_h, align_idx_m = align_idx(tx_seq_h, tx_seq_m) # make alignment
                else:
                    align_idx_h, align_idx_m = align_idx(tx_seq_h, tx_seq_m, glo=False)
                # check if mutation with flank size within the index
                mut_tx_idx_h = [tx_idx_h.index(x) for x in tx_idx_h if x in range(start_h, end_h+1)]
                new_tx_seq_h = get_modify_seq(ty_h, mut_tx_idx_h, tx_seq_h, alt_seq_h)
                tx_change_h = get_tx_change(tx_id_h,mut_tx_idx_h[0],mut_tx_idx_h[-1],ref_seq_h,alt_seq_h)
                GateFlankH(mut_tx_idx_h[0], mut_tx_idx_h[-1], align_idx_h, align_idx_m, mut_tx_idx_h)
                if status: # if so, get it to the nt change
                    mut_tx_idx_m = transfer_aligned_idx(mut_tx_idx_h, align_idx_h, align_idx_m)
                    flank_size_left, flank_size_right = get_flank(mut_tx_idx_h, align_idx_h, mut_tx_idx_m, align_idx_m,'nt')
                    mut_dna_idx_m_ori = mut_dna_idx_m = [tx_idx_m[x] for x in mut_tx_idx_m]
                    if ty_h == 'INS':
                        ref_seq_m_ori = ref_seq_m =''.join(tx_seq_m[x] for x in range(mut_tx_idx_m[0]+1, mut_tx_idx_m[-1]))
                    else:
                        ref_seq_m_ori = ref_seq_m = ''.join(tx_seq_m[x] for x in mut_tx_idx_m)
                    alt_seq_m_ori = alt_seq_m = alt_seq_h
                    tx_change_m_ori = tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m,alt_seq_m)
                    classification_m = get_non_coding(set([tx_region_anno_m[x] for x in mut_tx_idx_m]))
                    ty_m = ty_h
                    error, statement = 2, 'Class 2: This mutation can be modeled, but the effect may not be consistent.'
                # no alternative
                x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m, flank_size_left, flank_size_right)
                return x
            
            ################################################category 1: for intron################################################
            elif  category == 1:  # intron
                t = ','.join(set_of_mut)
                classification_h = f'Non_Coding,{t}'
                tx_seq_h, tx_idx_h, tx_region_anno_h, region_idx_h, tx_seq_m, tx_idx_m, tx_region_anno_m, region_idx_m = get_splicing(category = category)
                gate_memory(memory_protect=memory_protect, memory_size=memory_size, len_1 = len(tx_seq_h), len_2 = len(tx_seq_m))
                if n_exon_h == n_exon_m:
                    align_idx_h, align_idx_m = align_idx(tx_seq_h, tx_seq_m) # make alignment
                else:
                    align_idx_h, align_idx_m = align_idx(tx_seq_h, tx_seq_m, glo=False)
                # check if mutation with flank size within the index
                mut_tx_idx_h = [tx_idx_h.index(x) for x in tx_idx_h if x in range(start_h, end_h+1)]
                new_tx_seq_h = get_modify_seq(ty_h, mut_tx_idx_h, tx_seq_h, alt_seq_h)

                tx_idx_h_for_hgvsc, tx_idx_m_for_hgvsc = get_splicing(category =2)[1], get_splicing(category =2)[4]
                def find_nearest_E(sequence, tx_idx_for_hgvsc, tx_idx, coord):
                    left = right = coord
                    len_seq = len(sequence)
                    while left >= 0 or right < len_seq:
                        # left search
                        if left >= 0:
                            if sequence[left] == 'E':
                                loc = tx_idx_for_hgvsc.index(tx_idx[left])
                                distance = coord - left
                                return loc, '+', distance
                            left -= 1
                        # right search
                        if right < len_seq:
                            if sequence[right] == 'E':
                                loc = tx_idx_for_hgvsc.index(tx_idx[right])
                                distance = right - coord
                                return loc, '-', distance
                            right += 1
                    return None, None
                
                loc_h, direct_h, dist_h = find_nearest_E(tx_region_anno_h, tx_idx_h_for_hgvsc, tx_idx_h, mut_tx_idx_h[0])
                exon_h = region_idx_h[mut_tx_idx_h[0]]

                if direct_h == '-':
                    x = -1
                    dist_h_return = dist_h*(-1)
                else:
                    x = 1
                    dist_h_return = dist_h
                    
                if (dist_h <= 2) or (abs(x*dist_h+len(mut_tx_idx_h)-1) <=2):
                    classification_h = 'Splice_Site'
                    pro_change_h= f'X{loc_h//3}_splice'
                # customize tx_change
                if '5' in t:
                    tx_change_h = get_tx_change(tx_id_h,mut_tx_idx_h[0],mut_tx_idx_h[-1],ref_seq_h,alt_seq_h, '5')
                elif '3' in t:
                    tx_change_h = get_tx_change(tx_id_h,mut_tx_idx_h[0],mut_tx_idx_h[-1],ref_seq_h,alt_seq_h, '3', len(tx_seq_h))
                else:
                    tx_change_h = get_tx_change(tx_id_h,mut_tx_idx_h[0],mut_tx_idx_h[-1],ref_seq_h,alt_seq_h, 'I', loc_h, direct_h, dist_h)

                GateFlankH(mut_tx_idx_h[0], mut_tx_idx_h[-1], align_idx_h, align_idx_m, mut_tx_idx_h)
                
                if status: # if so, get it to the nt change
                    mut_tx_idx_m = transfer_aligned_idx(mut_tx_idx_h, align_idx_h, align_idx_m)
                    flank_size_left, flank_size_right = get_flank(mut_tx_idx_h, align_idx_h, mut_tx_idx_m, align_idx_m,'nt')
                    mut_dna_idx_m_ori = mut_dna_idx_m = [tx_idx_m[x] for x in mut_tx_idx_m]
                    if ty_h == 'INS':
                        ref_seq_m_ori = ref_seq_m =''.join(tx_seq_m[x] for x in range(mut_tx_idx_m[0]+1, mut_tx_idx_m[-1]))
                    else:
                        ref_seq_m_ori = ref_seq_m = ''.join(tx_seq_m[x] for x in mut_tx_idx_m)
                    new_tx_seq_m_ori = new_tx_seq_m = get_modify_seq(ty_h, mut_tx_idx_m, tx_seq_m, alt_seq_h)
                    alt_seq_m_ori = alt_seq_m = alt_seq_h
                    loc_m, direct_m, dist_m = find_nearest_E(tx_region_anno_m, tx_idx_m_for_hgvsc, tx_idx_m, mut_tx_idx_m[0])
                    if '5' in t:
                        tx_change_m_ori = tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m,alt_seq_m, '5')
                    elif '3' in t:
                        tx_change_m_ori = tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m,alt_seq_m, '3', len(tx_seq_m))
                    else:
                        tx_change_m_ori = tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m,alt_seq_m, 'I', loc_m, direct_m, dist_m)
                    if direct_m == '-':
                        dist_m_return = dist_m*(-1)
                        x = -1
                    else:
                        dist_m_return = dist_m
                        x = 1

                    if (dist_m <= 2) or (abs(x*loc_m+len(mut_tx_idx_m)-1) <=2):
                        classification_m = 'Splice_Site'
                        pro_change_m_ori = pro_change_m = f'X{loc_m//3}_splice'
                        # customize tx_change
                    else:
                        classification_m = get_non_coding(set([tx_region_anno_m[x] for x in mut_tx_idx_m]))

                    exon_m = region_idx_m[mut_tx_idx_m[0]]

                    ty_m = ty_h
                    error, statement = 2, 'Class 2: This mutation can be modeled, but the effect may not be consistent.'
                # no alternative
                x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m, flank_size_left, flank_size_right, dist_h_return, dist_m_return)
                return x
            ################################################category 2: for stop codon################################################
            elif category == 2:
                tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m = get_splicing(category = 2)
                p_seq_h, p_seq_m = Translate(tx_seq_h), Translate(tx_seq_m)
                stop_loc_tx_h, stop_loc_tx_m = region_idx_h.index(-1), region_idx_m.index(-1)
                stop_loc_p_h, stop_loc_p_m = p_seq_h.index('*'), p_seq_m.index('*')

                mut_tx_idx_h = [tx_idx_h.index(x) for x in tx_idx_h if x in range(start_h, end_h+1)]
                mut_tx_idx_m = transfer_stop_index(mut_tx_idx_h, stop_loc_tx_h, stop_loc_tx_m)
                tx_change_h = get_tx_change(tx_id_h,mut_tx_idx_h[0],mut_tx_idx_h[-1],ref_seq_h,alt_seq_h)
                if ty_h == 'INS':
                    ref_seq_m_ori = ref_seq_m =''.join(tx_seq_m[x] for x in range(mut_tx_idx_m[0]+1, mut_tx_idx_m[-1]))
                else:
                    ref_seq_m_ori = ref_seq_m = ''.join(tx_seq_m[x] for x in mut_tx_idx_m)

                alt_seq_m_ori = alt_seq_m = alt_seq_h
                mut_dna_idx_m_ori = mut_dna_idx_m = [tx_idx_m[x] for x in mut_tx_idx_m]

                new_tx_seq_h = get_modify_seq(ty_h, mut_tx_idx_h, tx_seq_h, alt_seq_h)
                new_tx_seq_m_ori = new_tx_seq_m = get_modify_seq(ty_h, mut_tx_idx_m, tx_seq_m, alt_seq_h)
                tx_change_m_ori = tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m_ori,alt_seq_m_ori)
                new_p_seq_h , new_p_seq_m = Translate(new_tx_seq_h), Translate(new_tx_seq_m)

                mut_p_idx_h = unique_ordered(tx_to_p_idx(mut_tx_idx_h))
                mut_p_idx_m = unique_ordered(tx_to_p_idx(mut_tx_idx_m))
                new_mut_p_idx_h, new_mut_p_idx_m = get_new(mut_tx_idx_h, ref_seq_h, alt_seq_h), get_new(mut_tx_idx_m, ref_seq_m, alt_seq_m)
                mut_p_h, new_mut_p_h = get_seq(p_seq_h,mut_p_idx_h), get_seq(new_p_seq_h,new_mut_p_idx_h)
                mut_p_m, new_mut_p_m = get_seq(p_seq_m,mut_p_idx_m), get_seq(new_p_seq_m,new_mut_p_idx_m)
                cor_mut_p_idx_h, cor_new_mut_p_idx_h, cor_mut_p_h, cor_new_mut_p_h = get_cor(mut_p_idx_h, new_mut_p_idx_h, mut_p_h, new_mut_p_h)
                cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m = get_cor(mut_p_idx_h, new_mut_p_idx_h, mut_p_h, new_mut_p_h)
                # to see if it is a non_stop or a silent mutation
                # non_stop, silent, or nonsense
                # return p_seq_h, new_p_seq_h, stop_loc_p_h, mut_p_h, mut_p_idx_h, cor_mut_p_idx_h, cor_new_mut_p_idx_h, cor_mut_p_h, cor_new_mut_p_h, category, ty_h, ref_seq_h, alt_seq_h
                classification_h, pro_change_h = type_snop(p_seq_h, new_p_seq_h, stop_loc_p_h, mut_p_h, mut_p_idx_h, cor_mut_p_idx_h, cor_new_mut_p_idx_h, cor_mut_p_h, cor_new_mut_p_h, category, ty_h, ref_seq_h, alt_seq_h, nonstop_size)
                classification_m_ori, pro_change_m_ori = type_snop(p_seq_m, new_p_seq_m, stop_loc_p_m, mut_p_m, mut_p_idx_m, cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m, category, ty_h, ref_seq_m, alt_seq_m, nonstop_size)

            ################################################category 3: for cds mutations################################################
            elif category == 3:
                # splicing for the coding sequence
                tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m = get_splicing(category = 2)
                p_seq_h, p_seq_m = Translate(tx_seq_h), Translate(tx_seq_m)
                stop_loc_tx_h, stop_loc_tx_m = region_idx_h.index(-1), region_idx_m.index(-1)
                stop_loc_p_h, stop_loc_p_m = p_seq_h.index('*'), p_seq_m.index('*')
                # Translate
                # All info for human
                mut_tx_idx_h = [tx_idx_h.index(x) for x in tx_idx_h if x in range(start_h, end_h+1)]
                tx_change_h = get_tx_change(tx_id_h,mut_tx_idx_h[0],mut_tx_idx_h[-1],ref_seq_h,alt_seq_h)
                new_tx_seq_h = get_modify_seq(ty_h, mut_tx_idx_h, tx_seq_h, alt_seq_h)
                exon_h = ','.join(list(set([region_idx_h[x] for x in mut_tx_idx_h])))
                new_p_seq_h = Translate(new_tx_seq_h)
                mut_p_idx_h = unique_ordered(tx_to_p_idx(mut_tx_idx_h))
                new_mut_p_idx_h = get_new(mut_tx_idx_h, ref_seq_h, alt_seq_h)
                mut_p_h, new_mut_p_h = get_seq(p_seq_h,mut_p_idx_h), get_seq(new_p_seq_h,new_mut_p_idx_h)
                cor_mut_p_idx_h, cor_new_mut_p_idx_h, cor_mut_p_h, cor_new_mut_p_h = get_cor(mut_p_idx_h, new_mut_p_idx_h, mut_p_h, new_mut_p_h)
                
                # check classification 
                    # non-sense, could find alternative
                        # non-sense-stop: the same
                    # missense, could find alternative
                    # silent, could find alternative
                        # non-sense-stop: the same
                    
                        # if_in_frame() (may be different for ins and del)

                    # will not effect splicing!
                    # frame-shifting: keep the same transcript as the non-stop one!
                    # in-frame-complete or in-frame-disturb: keep the same complete transcript!
                
                classification_h, pro_change_h = type_snop(p_seq_h, new_p_seq_h, stop_loc_p_h, mut_p_h, mut_p_idx_h, cor_mut_p_idx_h, cor_new_mut_p_idx_h, cor_mut_p_h, cor_new_mut_p_h, category, ty_h, ref_seq_h, alt_seq_h, nonstop_size)

                # alignment
                if batch:
                    if direction in ['h2m','h2h','m2m']:
                        align_p_idx_h, align_p_idx_m = align_input[3]
                    else:
                        align_p_idx_h, align_p_idx_m = [align_input[3][-1],align_input[3][0]]
                else:
                    gate_memory(memory_protect=memory_protect, memory_size=memory_size, len_1 = stop_loc_p_h, len_2 = stop_loc_p_m)
                    align_p_idx_h, align_p_idx_m = align_idx(p_seq_h[:stop_loc_p_h+1], p_seq_m[:stop_loc_p_m+1]) # make alignment
                align_tx_idx_h, align_tx_idx_m = p_to_tx_3(align_p_idx_h), p_to_tx_3(align_p_idx_m)

                # flank size gating
                GateFlankH(mut_p_idx_h[0], mut_p_idx_h[-1], align_p_idx_h, align_p_idx_m, mut_p_idx_h)
                
                if not status: # if so, get it to the nt change
                    x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m,  flank_size_left, flank_size_right) # should have human information before this
                    return x
                else:
                    mut_tx_idx_m = transfer_aligned_idx(mut_tx_idx_h, align_tx_idx_h, align_tx_idx_m)
                    mut_p_idx_m = transfer_aligned_idx(mut_p_idx_h, align_p_idx_h, align_p_idx_m)
                    flank_size_left, flank_size_right = get_flank(mut_p_idx_h, align_p_idx_h, mut_p_idx_m, align_p_idx_m,'aa')
                    if ty_h == 'INS':
                        ref_seq_m_ori = ref_seq_m =''.join(tx_seq_m[x] for x in range(mut_tx_idx_m[0]+1, mut_tx_idx_m[-1]))
                    else:
                        ref_seq_m_ori = ref_seq_m = ''.join(tx_seq_m[x] for x in mut_tx_idx_m)
                    alt_seq_m_ori = alt_seq_m = alt_seq_h
                    
                    tx_change_m_ori = tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m_ori,alt_seq_m_ori)

                    mut_dna_idx_m_ori = mut_dna_idx_m = [tx_idx_m[x] for x in mut_tx_idx_m]
                    exon_m = ','.join(list(set([region_idx_m[x] for x in mut_tx_idx_m])))
                    new_tx_seq_m_ori = new_tx_seq_m = get_modify_seq(ty_h, mut_tx_idx_m, tx_seq_m, alt_seq_h)
                    new_p_seq_m = Translate(new_tx_seq_m)
                    mut_p_idx_m = unique_ordered(tx_to_p_idx(mut_tx_idx_m))
                    new_mut_p_idx_m = get_new(mut_tx_idx_m, ref_seq_m, alt_seq_m)
                    mut_p_m, new_mut_p_m = get_seq(p_seq_m,mut_p_idx_m), get_seq(new_p_seq_m,new_mut_p_idx_m)
                    cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m = get_cor(mut_p_idx_m, new_mut_p_idx_m, mut_p_m, new_mut_p_m)
                    
                    # return p_seq_m, new_p_seq_m, stop_loc_p_m, mut_p_m, mut_p_idx_m, cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m, category, ty_h, ref_seq_m, alt_seq_m
                    classification_m_ori, pro_change_m_ori = type_snop(p_seq_m, new_p_seq_m, stop_loc_p_m, mut_p_m, mut_p_idx_m, cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m, category, ty_h, ref_seq_m, alt_seq_m, nonstop_size)
                    
                    
            # gate same mutation
            if classification_h == 'Frame_Shifting':
                start_m, end_m, alt_seq_m, ref_seq_m = start_m_ori, end_m_ori, alt_seq_m_ori, ref_seq_m_ori
                classification_m, pro_change_m  = classification_m_ori, pro_change_m_ori
                ty_m = ty_h
                error, statement = 2, 'Class 2: This mutation can be modeled, but the effect may not be consistent.'
                if param is not None:
                    if param == 'BE':
                        GateBE(ref_seq_m, alt_seq_m)
                x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m,  mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m,  flank_size_left, flank_size_right)
                return x
                          
            if classification_h == 'Non_Stop':
                # get ori m info and gating
                GateNonstop(new_p_seq_m, stop_loc_p_m)
                if param is not None:
                    if param == 'BE':
                        GateBE(ref_seq_m, alt_seq_m)
                if status:
                    start_m, end_m, alt_seq_m, ref_seq_m = start_m_ori, end_m_ori, alt_seq_m_ori, ref_seq_m_ori
                    classification_m, pro_change_m  = classification_m_ori, pro_change_m_ori
                    ty_m = ty_h
                    x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m,  flank_size_left, flank_size_right)
                    return x
                else:
                    IfAlternative() # del and ins: do not find alternative.
                    if error != 2:
                        if search_alternative:
                            x =  FindAlternative(max_alternative, mut_tx_idx_m, gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_p_idx_m,new_mut_p_idx_m, mut_p_idx_h, new_mut_p_idx_h, align_p_idx_h, align_p_idx_m)
                            
                            return x
                        else:
                            x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m,mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m,  flank_size_left, flank_size_right)
                            return x
            
            else:
                if (cor_mut_p_m == cor_mut_p_h) and (cor_new_mut_p_m == cor_new_mut_p_h):
                    status, error, statement = status,error,statement = True, 0, "Class 0: This mutation can be originally modeled."
                    if param is not None:
                        if param == 'BE':
                            GateBE(ref_seq_m, alt_seq_m)
                else:
                    status, error, statement = False, 6, 'Class 6: This mutation cannot be originally modeled.'

                if status:
                    start_m, end_m, alt_seq_m, ref_seq_m = start_m_ori, end_m_ori, alt_seq_m_ori, ref_seq_m_ori
                    classification_m, pro_change_m  = classification_m_ori, pro_change_m_ori
                    ty_m = ty_h
                    x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m,  flank_size_left, flank_size_right)
                    return x
                else:
                    if search_alternative:
                        x =  FindAlternative(max_alternative, mut_tx_idx_m, gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_p_idx_m,new_mut_p_idx_m, mut_p_idx_h, new_mut_p_idx_h, align_p_idx_h, align_p_idx_m)
                        return x
                    else:
                        x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m, new_mut_p_idx_h, new_mut_p_idx_m, flank_size_left, flank_size_right)
                        return x
                          
    elif (coor == 'aa'):
        category = 3
        tx_seq_h, tx_idx_h, region_idx_h, tx_seq_m, tx_idx_m, region_idx_m = get_splicing(category = 2)
        mut_list = [tx_idx_h[x] for x in p_to_tx_3(list(range(start_h-1, end_h-1+1)))]
        if strand_h == '+':
            start_h, end_h = mut_list[0],mut_list[-1]
        else:
            start_h, end_h = mut_list[-1],mut_list[0]
        check_mut_info()
        p_seq_h, p_seq_m = Translate(tx_seq_h), Translate(tx_seq_m)
        stop_loc_tx_h, stop_loc_tx_m = region_idx_h.index(-1), region_idx_m.index(-1)
        stop_loc_p_h, stop_loc_p_m = p_seq_h.index('*'), p_seq_m.index('*')
        # Translate
        # All info for human
        mut_tx_idx_h = [tx_idx_h.index(x) for x in tx_idx_h if x in range(start_h, end_h+1)]
        tx_change_h = get_tx_change(tx_id_h,mut_tx_idx_h[0],mut_tx_idx_h[-1],ref_seq_h,alt_seq_h)
        new_tx_seq_h = get_modify_seq(ty_h, mut_tx_idx_h, tx_seq_h, alt_seq_h)
        exon_h = ','.join(list(set([region_idx_h[x] for x in mut_tx_idx_h])))
        new_p_seq_h = Translate(new_tx_seq_h)
        mut_p_idx_h = unique_ordered(tx_to_p_idx(mut_tx_idx_h))
        new_mut_p_idx_h = get_new(mut_tx_idx_h, ref_seq_h, alt_seq_h)
        mut_p_h, new_mut_p_h = get_seq(p_seq_h,mut_p_idx_h), get_seq(new_p_seq_h,new_mut_p_idx_h)
        cor_mut_p_idx_h, cor_new_mut_p_idx_h, cor_mut_p_h, cor_new_mut_p_h = get_cor(mut_p_idx_h, new_mut_p_idx_h, mut_p_h, new_mut_p_h)
        
        # check classification 
            # non-sense, could find alternative
                # non-sense-stop: the same
            # missense, could find alternative
            # silent, could find alternative
                # non-sense-stop: the same
            
                # if_in_frame() (may be different for ins and del)

            # will not effect splicing!
            # frame-shifting: keep the same transcript as the non-stop one!
            # in-frame-complete or in-frame-disturb: keep the same complete transcript!
        
        classification_h, pro_change_h = type_snop(p_seq_h, new_p_seq_h, stop_loc_p_h, mut_p_h, mut_p_idx_h, cor_mut_p_idx_h, cor_new_mut_p_idx_h, cor_mut_p_h, cor_new_mut_p_h, category, ty_h, ref_seq_h, alt_seq_h, nonstop_size)

        # alignment
        if batch:
            if direction in ['h2m','h2h','m2m']:
                align_p_idx_h, align_p_idx_m = align_input[3]
            else:
                align_p_idx_h, align_p_idx_m = [align_input[3][-1],align_input[3][0]]
        else:
            gate_memory(memory_protect=memory_protect, memory_size=memory_size, len_1 = stop_loc_p_h, len_2 = stop_loc_p_m)
            align_p_idx_h, align_p_idx_m = align_idx(p_seq_h[:stop_loc_p_h+1], p_seq_m[:stop_loc_p_m+1]) # make alignment
        align_tx_idx_h, align_tx_idx_m = p_to_tx_3(align_p_idx_h), p_to_tx_3(align_p_idx_m)

        # flank size gating
        GateFlankH(mut_p_idx_h[0], mut_p_idx_h[-1], align_p_idx_h, align_p_idx_m, mut_p_idx_h)
        if not status: # if so, get it to the nt change
            x = final_return(gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_tx_idx_m, mut_p_idx_h, mut_p_idx_m,new_mut_p_idx_h, new_mut_p_idx_m, flank_size_left, flank_size_right) # should have human information before this
            if status == False:
                ty_m = classification_m = exon_m = start_m_ori = end_m_ori = ref_seq_m_ori = alt_seq_m_ori = pro_change_m_ori = start_m = end_m = ref_seq_m = alt_seq_m = pro_change_m = None
            return x
        else:
            mut_tx_idx_m = transfer_aligned_idx(mut_tx_idx_h, align_tx_idx_h, align_tx_idx_m)
            mut_p_idx_m = transfer_aligned_idx(mut_p_idx_h, align_p_idx_h, align_p_idx_m)
            flank_size_left, flank_size_right = get_flank(mut_p_idx_h, align_p_idx_h, mut_p_idx_m, align_p_idx_m,'aa')
            if ty_h == 'INS':
                ref_seq_m_ori = ref_seq_m =''.join(tx_seq_m[x] for x in range(mut_tx_idx_m[0]+1, mut_tx_idx_m[-1]))
            else:
                ref_seq_m_ori = ref_seq_m = ''.join(tx_seq_m[x] for x in mut_tx_idx_m)
            alt_seq_m_ori = alt_seq_m = alt_seq_h
            
            tx_change_m_ori = tx_change_m = get_tx_change(tx_id_m,mut_tx_idx_m[0],mut_tx_idx_m[-1],ref_seq_m_ori,alt_seq_m_ori)

            mut_dna_idx_m_ori = mut_dna_idx_m = [tx_idx_m[x] for x in mut_tx_idx_m]
            exon_m = ','.join(list(set([region_idx_m[x] for x in mut_tx_idx_m])))
            new_tx_seq_m_ori = new_tx_seq_m = get_modify_seq(ty_h, mut_tx_idx_m, tx_seq_m, alt_seq_h)
            new_p_seq_m = Translate(new_tx_seq_m)
            mut_p_idx_m = unique_ordered(tx_to_p_idx(mut_tx_idx_m))
            new_mut_p_idx_m = get_new(mut_tx_idx_m, ref_seq_m, alt_seq_m)
            mut_p_m, new_mut_p_m = get_seq(p_seq_m,mut_p_idx_m), get_seq(new_p_seq_m,new_mut_p_idx_m)
            cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m = get_cor(mut_p_idx_m, new_mut_p_idx_m, mut_p_m, new_mut_p_m)
            
            # return p_seq_m, new_p_seq_m, stop_loc_p_m, mut_p_m, mut_p_idx_m, cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m, category, ty_h, ref_seq_m, alt_seq_m
            classification_m_ori, pro_change_m_ori = type_snop(p_seq_m, new_p_seq_m, stop_loc_p_m, mut_p_m, mut_p_idx_m, cor_mut_p_idx_m, cor_new_mut_p_idx_m, cor_mut_p_m, cor_new_mut_p_m, category, ty_h, ref_seq_m, alt_seq_m, nonstop_size)

            x = FindAlternative(max_alternative, mut_tx_idx_m, gene_name_h, gene_id_h, tx_id_h, chr_h, n_exon_h, strand_h, match_h, start_h, end_h, ref_seq_h, alt_seq_h, tx_change_h, pro_change_h, classification_h, exon_h, ty_h, gene_name_m, gene_id_m, tx_id_m, chr_m, n_exon_m, strand_m, ty_m, classification_m, exon_m, start_m_ori, end_m_ori, ref_seq_m_ori, alt_seq_m_ori, tx_change_m_ori, pro_change_m_ori, start_m, end_m, ref_seq_m, alt_seq_m, tx_change_m, pro_change_m, tx_seq_h, tx_seq_m, new_tx_seq_h, new_tx_seq_m_ori, new_tx_seq_m, mut_tx_idx_h, mut_p_idx_m,new_mut_p_idx_m, mut_p_idx_h, new_mut_p_idx_h,align_p_idx_h, align_p_idx_m)
            return x

def model_batch(df, records_h, index_list_h, records_m, index_list_m, db_h, db_m, ver, param = 'default', direction = 'h2m',coor = 'nc', search_alternative = True, max_alternative = 5, nonstop_size = 300, splicing_size = 30, show_sequence = False, align_input = None, memory_protect = True, memory_size = 10000, bind = False):
    """
    Batch modeling of human variants in the mouse genome.  

    Parameters:  
        - df (Pandas DataFrame): Must include columns {'start_h','end_h','type_h','ref_seq_h','alt_seq_h','tx_id_h','tx_id_m','index'}.
        - ecords_h, index_list_h, records_m, index_list_m: reference genome files
        - db_h, db_m: genomic annotation files
        - ver (int): specify the version of human, one of 37/38.
        - param (optional): set param = 'BE' and will only output base editing modelable results.  
        - direction (optional): str, set the modeling direction by 'h2m' (default) or 'm2h'.
        - coor (optional): default = 'nc'. set input = 'aa' and will be compatable with input of amino acid variants.
        - search_alternative (optional): set search_alternative = False and will only output original modeling results.  
        - max_alternative (optional): the maximum number of output alternatives of one human variants.  
        - nonstop_size (optional): the length of neucleotides that are included after the stop codon for alignment and translation in case of the nonstop mutations or frame shifting mutations.  
        - splicing_size (optional): 
        - batch (optional): set batch = True and will use input align_dict to save time in batch processing.  
        - show_sequence (optional): set batch = True and will output the whole sequences.  
        - align_dict (optional): input a prepared dictionary of alignment indexes to save time in batch processing.  
        - memory_protect (optional): default True. Break long alignments that may lead to death of the kernel.    
        - memory_size (optional): maxlength of aligned sequence when memory_protect == True.  
        - bind (optional): to bind the output dataframe with the original input or not.  

    Return:  
        - Two dataframes. The first dataframe is the processed original dataframe. The second dataframe contains all rows that are not successfully processed.

    Example:
        >>> h2m.model_batch(df, records_h, index_list_h, records_m, index_list_m, db_h, db_m, ver = 37, param = 'BE')
    """

    result_list = []
    error_list = []
    df_list = []
    r = 0
    if (ver is None):
        raise ValueError('Error 1: Please include genome assembly version according to your transcript ID by adding ver=37 or ver=38.')
    elif ver not in [37,38]:
        raise ValueError('Error 1: Please include genome assembly version accfording to your transcript ID by adding ver=37 or ver=38.')

    col_list = list(df.columns)
    if 'ref_seq_h' in col_list:
        input_name_list = ['start_h','end_h','type_h','ref_seq_h','alt_seq_h','tx_id_h','tx_id_m','index']
    elif 'ref_seq_m' in col_list:
        input_name_list = ['start_m','end_m','type_m','ref_seq_m','alt_seq_m','tx_id_h','tx_id_m','index']
    else:
        input_name_list = ['start','end','type','ref_seq','alt_seq','tx_id_h','tx_id_m','index']
        
    if not all(x in col_list for x in input_name_list):
        raise ValueError("Error 2: Please include columns of 'start','end','type','ref_seq','alt_seq','tx_id_h','tx_id_m','index'.")
    else:
        col_idx = [col_list.index(x) for x in input_name_list]

    df_mutation = df.sort_values(by=['tx_id_h','tx_id_m']).reset_index(drop=True)
    list_of_counts = df_mutation[['tx_id_h','tx_id_m']].value_counts(sort=False).to_list()

    for c_index, c_count in enumerate(list_of_counts):
        print(f'{c_index+1}/{len(list_of_counts)}',end='\r')
        tx_id_h, tx_id_m = df_mutation.iloc[r,col_idx[5:7]].to_list()
        #################################### align ###################################
        try:
            align_input = get_tx_and_align(records_h, index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, ver=ver, align = True)
        except:
            # if align goes wrong
            error_list = error_list + list(range(r,r+c_count))
            r = r+c_count
        else:
            # if align doesn't go wrong
            #################################### model ###################################
            for _ in range(c_count):
                start_h, end_h, ty_h, ref_seq_h, alt_seq_h = df_mutation.iloc[r,col_idx[:5]].to_list()
                unique_index = df_mutation.iloc[r,col_idx[7]]
                
                try:
                    result = model(records_h, index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, start_h, end_h, ref_seq_h, alt_seq_h, ty_h, ver = ver,direction = direction, batch=True, align_input=align_input, param = param, coor = coor, search_alternative = search_alternative, max_alternative = max_alternative, nonstop_size = nonstop_size, splicing_size = splicing_size, show_sequence = show_sequence, memory_protect = memory_protect, memory_size = memory_size)
                except:
                    # if goes wrong
                    error_list.append(r)
                    r = r+1
                else:
                    # if doesn't go wrong
                    r = r+1
                    if type(result) == list: # multiple alternative modeling
                        for k in result:
                            k['index'] = unique_index
                            result_list.append(k)
                    else:
                        result['index'] = unique_index
                        result_list.append(result)
            
        df_list.append(pd.DataFrame(result_list))
        result_list = []


    df_result = pd.concat(df_list, ignore_index=True)
    
    col = [col for col in df.columns if col not in ['start','end','type','ref_seq','alt_seq','tx_id_h','tx_id_m']]
    if bind == True:
        df_result = pd.merge(df_result, df[col], how = 'left', on = 'index')

    if len(df_result)>1:
        df_result = df_result.sort_values(by='index',ascending=True).reset_index(drop=True)
    df_wro = df_mutation.iloc[error_list,:].reset_index(drop=True)
    if len(error_list)==0:
        print('No error occurs.')
    else:
        print('There were rows that could not be processed.')
    return(df_result, df_wro)

def vcf_reader(path, keep = True):
    """
    Generate MAF-style h2m input from VCF-format data, including genomAD data.  

    Parameter:
        - path (str): the path of input csv data.
        - keep (bool): True: keep all the original columns in the dataframe/ False: keep the necesssary columns for h2m only. Default to False.  

    Output: 
        An input dataframe for h2m modeling.

    Example:
        >>> filepath = '.../gnomAD_v4.0.0_ENSG00000141510_2024_02_07_11_36_03.csv'
        >>> df = h2m.vcf_reader('','TP53')
    """
    
    df = pd.read_csv(path)
    col_name_old = ['Position','Reference','Alternate','gnomAD ID']
    col_name_new = ['start_h', 'ref_seq_h', 'alt_seq_h','ID']
    df = df.rename(columns = dict(zip(col_name_old,col_name_new)))
    col_index = [list(df.columns).index(x) for x in ['start_h','ref_seq_h','alt_seq_h','ID']]
    if keep == False:
        df = df.iloc[:,col_index]
    df['index'] = list(range(len(df)))
    df['format'] = 'VCF'
    return df

def vcf_to_maf(df):
    """
    Convert vcf-style input into maf-style input, while keeping other columns intact.
    """
    def iter_row(row):
        # Extract the relevant columns
        start, ref, alt = row['start_h'], row['ref_seq_h'], row['alt_seq_h']
        # Perform transformations
        if len(ref) == len(alt):
            if len(ref) == 1:
                type_ = 'SNP'
            elif len(ref) == 2:
                type_ = 'DNP'
            elif len(ref) == 3:
                type_ = 'TNP'
            elif len(ref) >= 4:
                type_ = 'ONP'
            end = start + len(ref) - 1
        elif (len(ref) == 1) and (len(alt) >= 2):
            end = start + 1
            alt = alt[1:]
            ref = '-'
            type_ = 'INS'
        elif (len(ref) >= 2) and (len(alt) == 1):
            start += 1
            ref = ref[1:]
            alt = '-'
            end = start + len(ref) - 1
            type_ = 'DEL'
        elif len(alt) > len(ref):
            start += 1
            ref = ref[1:]
            alt = alt[1:]
            end = start + len(ref) - 1
            type_ = 'INS'
        elif len(alt) < len(ref):
            start += 1
            ref = ref[1:]
            alt = alt[1:]
            end = start + len(ref) - 1
            type_ = 'DEL'
        row['start_h'], row['end_h'], row['type_h'], row['ref_seq_h'], row['alt_seq_h'] = start, end, type_, ref, alt
        return row

    df['end_h'] = df['start_h']
    df['type_h'] = 'SNP'
    df = df.apply(iter_row, axis=1)
    priority_columns = ['ID','start_h','end_h','ref_seq_h','alt_seq_h','type_h']
    remaining_columns = [col for col in df.columns if col not in priority_columns]
    new_column_order = priority_columns + remaining_columns
    df = df[new_column_order]
    df['format'] = 'MAF'    
    return df

def cbio_reader(path = None, df = None, keep = True):
    """
    Generate MAF-sty h2m input from cbioportal data.  

    Parameter:
        - path (str): the path of mutation data in txt format. 
        - keep (bool): True: keep all the original columns in the dataframe/ False: keep the necesssary columns for h2m only. Default to False.  

    Output: 
        An input dataframe for h2m modeling.

    Example:   
        >>> h2m.cbio_reader('.../data_mutations.txt', keep=False)
    """
    df = pd.read_csv(path, header=0, sep='\t', comment="#", na_values = 'Not Applicable')
    df = df.dropna(subset = 'HGVSc')
    df['HGVSc'] = [x.split(':')[0] for x in df['HGVSc']]
    col_name_old = ['Hugo_Symbol','HGVSc','Start_Position','End_Position','Variant_Type','Reference_Allele','Tumor_Seq_Allele2']
    col_name_new = ['gene_name_h', 'tx_id_h', 'start_h', 'end_h', 'type_h', 'ref_seq_h', 'alt_seq_h']
    df = df.rename(columns = dict(zip(col_name_old,col_name_new)))
    col_index = [list(df.columns).index(x) for x in col_name_new]
    if keep == False:
        df = df.iloc[:,col_index]
    df = df.drop_duplicates().reset_index(drop=True)
    df['index'] = list(range(len(df)))
    col = len(list(df.columns))
    df = df.iloc[:,[col-1]+list(range(col-1))]
    df['format'] = 'MAF'
    return df

def get_variant_type(df, ref_col, alt_col, col_name = 'type_h'):
    """
    Generate h2m/cbio-style neucleotide variant type annotations.

    Parameter:
        - df (str): a dataframe of mutations including columns for both reference and alternate sequences
        - ref_col (str): column name for reference sequences.
        - alt_col (str): column name for alternate sequences.

    Output: 
        The input dataframe but with a column of vairant type added.

    Example:   
        >>> h2m.get_variant_type(df, 'ref_seq_h','alt_seq_h')
    """
    df_change = df
    df_change[ref_col] = df_change[ref_col].fillna('')
    df_change[alt_col] = df_change[alt_col].fillna('')

    df_change[ref_col] = [str(x) for x in df_change[ref_col]]
    df_change[alt_col] = [str(x) for x in df_change[alt_col]]

    ref = df_change[ref_col]
    alt = df_change[alt_col]

    def change(ref, alt):
        if len(ref) == len(alt):
            if len(ref) ==1:
                r = 'SNP'
            elif len(ref) == 2:
                r = 'DNP'
            elif len(ref) == 3:
                r = 'TNP'
            elif len(ref) >3:
                r = 'ONP'
            else:
                r = None
        elif len(ref)<len(alt):
            r = 'INS'
        else:
            r = 'DEL'
        return r

    df_change[col_name] = [change(x,y) for x,y in zip(ref, alt)]
    return df_change

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

def clinvar_to_maf(df):
    """
    Convert clinvar-style input into maf-style input, while keeping other columns intact.
    """
    def iter_row(row):
        # Extract the relevant columns
        start, end, type_, ref, alt = row['start_h'], row['end_h'], row['type_h'], row['ref_seq_h'], row['alt_seq_h']
        # Perform transformations
        if len(ref) == len(alt):
            if len(ref) == 1:
                type_ = 'SNP'
            elif len(ref) == 2:
                type_ = 'DNP'
            elif len(ref) == 3:
                type_ = 'TNP'
            elif len(ref) >= 4:
                type_ = 'ONP'
        elif (len(ref) == 1) and (len(alt) >= 2):
            end += 1
            alt = alt[1:]
            ref = '-'
            type_ = 'INS'
        elif (len(ref) >= 2) and (len(alt) == 1):
            start += 1
            ref = ref[1:]
            alt = '-'
            type_ = 'DEL'
        elif len(alt) > len(ref):
            type_ = 'INS'
        elif len(alt) < len(ref):
            type_ = 'DEL'

        row['start_h'], row['end_h'], row['type_h'], row['ref_seq_h'], row['alt_seq_h'] = start, end, type_, ref, alt
        return row

    # Apply transformation to each row
    df = df.apply(iter_row, axis=1)
    df['format'] = 'MAF'
    return df

def visualization(model_result, flank_size = 0, print_size = 6):
    """
    Visualize h2m modeling results.  

    Parameter:
        - model_result (list): the output of `h2m.model(show_sequence = True)` function.
        - flank_size (int) (de). 
        - print_size (int): lenth of neucleotide/peptide included on both sides of the flank region.  

    Output: 
        A visualization plot.

    Example:   
        >>> model_result = h2m.model(records_h,index_list_h, records_m, index_list_m, db_h, db_m, tx_id_h, tx_id_m, 7577120, 7577120, 'C','T', ty_h = 'SNP', ver = 37, show_sequence=True)
        >>> h2m.visualization(model_result, flank_size = 2, print_size = 4)
    """
    def Translate(seq):
        return Bio.Seq.Seq(seq).transcribe().translate()
    
    def p_to_tx_3(p):
        list_of_condn = []
        for i in p:
            list_of_condn.extend([3*i, 3*i+1, 3*i+2])
        return(list_of_condn)

    dict_color_row = {'WT':'#F4E285',
                    'MUT':'#B2D3A8',
                    'ORI':'#F5D3C8',
                    'ALT':'#8E4162'}

    def GenerateAlignIdx(print_seq):
        list_idx = []
        sum = 0
        for i in print_seq:
            if i == '-':
                list_idx.append('-')
            else:
                list_idx.append(sum)
                sum = sum+1
        return list_idx

    def get_flank_idx(start_h, end_h, len_seq_h, flank_size):
        flank_left = min([start_h, flank_size])
        flank_right = min([len_seq_h-end_h-1, flank_size])
        s_flank, e_flank = start_h-flank_left, end_h+flank_right
        return s_flank, e_flank

    def print_tx_for_pep(list_idx, tx_seq, idx_ori, color):
        def change_to_tx(x):
                if x == '-':
                    return ['---','---']
                else:
                    idx = p_to_tx_3([x])
                    return [idx, tx_seq[idx[0]:idx[-1]+1]]
        seq_output = [change_to_tx(x)[1] for x in list_idx]
        list_output = [change_to_tx(x)[0] for x in list_idx]
        seq_output = ''.join(seq_output)
        list_output = [item for inner_list in list_output for item in inner_list]
        dict_color = dict()
        if len(idx_ori)>0:
            dict_color = dict(zip(idx_ori, [color]*len(idx_ori)))
        return list_output, seq_output, dict_color

    def get_color_dict(list_idx, mut_index, flank_size, color_mut):
        output_dict = dict()

        for x in list_idx:
            if x in mut_index:
                color = color_mut
            elif x in range(mut_index[0]-flank_size, mut_index[-1]+flank_size+1):
                color = 'lightblue'
            else:
                color = 'lightgrey'
            output_dict.update({x:color})

        return output_dict

    def get_color_dict_del(list_idx, mut_index, flank_size):
        output_dict = dict()
        for x in list_idx:
            if x in range(mut_index-flank_size, mut_index+flank_size):
                color = 'lightblue'
            else:
                color = 'lightgrey'
            output_dict.update({x:color})
        return output_dict

    def add_delta(list_idx, point, delta):
        if type(list_idx) == str:
            return list_idx[:point+1] + '-'*delta + list_idx[point+1:]
        elif type(list_idx) == list:
            return list_idx[:point+1] + ['-']*delta + list_idx[point+1:]

    def get_print_info_ori(tx_idx_ori, len_seq, len_mut_seq, print_ori_seq, print_index, flank_size, print_size):
        sh, eh = get_flank_idx(tx_idx_ori[0], tx_idx_ori[-1],len_seq, flank_size)
        shp, ehp = print_index.index(sh), print_index.index(eh)
        print_s, print_e = get_flank_idx(shp, ehp,len(print_ori_seq), print_size)
        list_idx = print_index[print_s:print_e+1] # whcih could have '-' inside
        list_seq = print_ori_seq[print_s:print_e+1]
        # list_idx_new = []
        # for idx in list_idx:
        #     if idx != '-':
        #         idx = str(int(idx)+1)
        #     list_idx_new.append(idx)
        color_tx = get_color_dict(list_idx, tx_idx_ori, flank_size, 'pink')
        delta = len_mut_seq - len_seq
        # if the mutated seq is longer than the original seq, add extra gaps
        if delta>0:
            list_idx = add_delta(list_idx, flank_size+print_size, delta)
            list_seq = add_delta(list_seq, flank_size+print_size, delta)
        return list_idx, list_seq, color_tx

    def get_align_tx_print(seq1, seq2):
        alignments = pairwise2.align.globalxx(seq1, seq2)
        x = pairwise2.format_alignment(*alignments[0]).split('\n')
        return x[0],x[2]

    def get_print_elements_ori(seq_h, seq_m, seq_mut_h, seq_mut_m, idx_h, idx_m, flank_size, print_size):
        # align
        # print alignment of h&m WT seqs
        print_ori_seq_h, print_ori_seq_m = get_align_tx_print(seq_h, seq_m)
        print_index_h, print_index_m = GenerateAlignIdx(print_ori_seq_h), GenerateAlignIdx(print_ori_seq_m)
        list_idx_h, print_ori_seq_h_short, color_tx_h_ori = get_print_info_ori(idx_h, len(seq_h), len(seq_mut_h), print_ori_seq_h, print_index_h, flank_size, print_size)
        list_idx_m, print_ori_seq_m_short, color_tx_m_ori = get_print_info_ori(idx_m, len(seq_m), len(seq_mut_m), print_ori_seq_m, print_index_m, flank_size, print_size)
        return list_idx_h, print_ori_seq_h_short, color_tx_h_ori, list_idx_m, print_ori_seq_m_short, color_tx_m_ori

    def get_print_elements_ori_non_stop(seq_h, seq_m, seq_mut_h, seq_mut_m, idx_h, idx_m, flank_size, print_size):
        def get_align_tx_non_stop(seq_h, seq_m):
            print1, print2 = '',''
            stop_1, stop_2 = seq_h.index('*'), seq_m.index('*')
            delta = stop_1 - stop_2
            if delta>0:
                print2 = print2 + '-'*delta
            elif delta<0:
                print1 = print1 + '-'*delta
            print1 = print1+seq_h
            print2 = print2+seq_m
            delta = len(print1) - len(print2)
            if delta>0:
                print2 = print2 + '-'*delta
            elif delta<0:
                print1 = print1 + '-'*delta
            return print1, print2
        print_ori_seq_h, print_ori_seq_m = get_align_tx_non_stop(seq_h, seq_m)
        print_index_h, print_index_m = GenerateAlignIdx(print_ori_seq_h), GenerateAlignIdx(print_ori_seq_m)
        list_idx_h, print_ori_seq_h_short, color_tx_h_ori = get_print_info_ori(idx_h, len(seq_h), len(seq_h), print_ori_seq_h, print_index_h, flank_size, print_size)
        list_idx_m, print_ori_seq_m_short, color_tx_m_ori = get_print_info_ori(idx_m, len(seq_m), len(seq_m), print_ori_seq_m, print_index_m, flank_size, print_size)
        return list_idx_h, print_ori_seq_h_short, color_tx_h_ori, list_idx_m, print_ori_seq_m_short, color_tx_m_ori
    
    def get_print_elements_mut(list_idx_ori, print_ori_seq, color_ori, tx_idx, ref_seq, alt_seq, ty, flank_size, print_size):
        def plus_k(l,k):
            l_o = []
            for x in l:
                if x == '-':
                    l_o.append('-')
                else:
                    l_o.append(int(x)+k)
            return l_o
        def minus_k(l,k):
            l_o = []
            for x in l:
                if x == '-':
                    l_o.append('-')
                else:
                    l_o.append(int(x)-k)
            return l_o
        
        if (ty != 'INS') and (ty != 'DEL'):
            print_ori_seq = print_ori_seq[:print_size+flank_size]+alt_seq+print_ori_seq[print_size+flank_size+len(ref_seq):]
            list_idx_ori = list_idx_ori
            color_ori_output = color_ori.copy()
            for nucleotide in tx_idx:
                color_ori_output.update({
                    nucleotide:'#A26EC4' #purple
                })
        elif ty == 'INS':
            delta = -len(ref_seq)+len(alt_seq)
            new_tx_idx = list(range(tx_idx[0], tx_idx[0]+2+delta))
            print_ori_seq = print_ori_seq[:print_size+flank_size+1]+alt_seq+print_ori_seq[print_size+flank_size+delta+1:]
            list_idx_ori = list_idx_ori[:print_size+flank_size] + new_tx_idx + plus_k(list_idx_ori[print_size+flank_size+delta+2:],delta)
            color_ori_output = get_color_dict(list_idx_ori, new_tx_idx, flank_size, '#A26EC4')
            color_ori_output.update({new_tx_idx[0]:'pink',new_tx_idx[-1]:'pink'})
        elif ty == 'DEL':
            delta = -len(alt_seq)+len(ref_seq)
            new_tx_idx = list(range(tx_idx[0], tx_idx[0]+len(alt_seq_h)))
            # add additional gaps
            print_ori_seq = print_ori_seq[:print_size+flank_size]+alt_seq+'-'*delta+print_ori_seq[print_size+flank_size+len(ref_seq):]
            
            list_idx_ori = list_idx_ori[:print_size+flank_size]+new_tx_idx+['-']*delta+minus_k(list_idx_ori[print_size+flank_size+len(ref_seq):],delta)
            if new_tx_idx == []: 
                color_ori_output = get_color_dict_del(list_idx_ori, tx_idx[0], flank_size)
            else:
                color_ori_output = get_color_dict(list_idx_ori, new_tx_idx, flank_size, '#A26EC4')
        return list_idx_ori, print_ori_seq, color_ori_output

    def get_print_elements_mut_pep(list_idx_h, print_ori_seq_h_short, color_p_h_ori, mut_p_seq_h ,p_idx_h_new, ty_h, class_h, flank_size, print_size):
        if ty_h in ['SNP','DNP','ONP']:
            list_idx_mut = list_idx_h
            print_seq_h_short_mut = get_seq_from_idx(list_idx_h, mut_p_seq_h)
            color_p_h_output = color_p_h_ori.copy()
            for key in color_p_h_output.keys():
                if color_p_h_output[key] == 'pink':
                    color_p_h_output[key] = '#A26EC4'
            return list_idx_h, print_seq_h_short_mut, color_p_h_output
        if ty_h == 'INS' or ty_h == 'DEL':
            # for before p_idx_h_new[0], keep everything as the ori
            # for after that, splice the same length
            stop_point = list_idx_h.index(p_idx_h_new[0])
            list_idx_mut = list_idx_h[:stop_point]
            q = p_idx_h_new[0]
            while len(list_idx_mut)<len(list_idx_h):
                list_idx_mut.append(q)
                q = q + 1
            print_seq_h_short_mut = get_seq_from_idx(list_idx_mut, mut_p_seq_h)

            if class_h == 'Frame_Shifting':
                list_key = [key for key in list_idx_mut if (key !='-') and (key >= p_idx_h_new[0])]
                color_p_h_output = dict(zip(list_key, len(list_key)*['#A26EC4']))
            else:
                color_p_h_output = dict(zip(p_idx_h_new, len(p_idx_h_new)*['#A26EC4']))
            return list_idx_mut, print_seq_h_short_mut, color_p_h_output

    def get_seq_from_idx(idx, seq):
        def change(x):
            if x =='-':
                return '-'
            else:
                return seq[x]
        return ''.join([change(x) for x in idx])

    def change_to_mutated(tx_idx_h, ty, alt_seq, ref_seq):
        if ty == 'INS':
            return list(range(tx_idx_h[0],tx_idx_h[0]+len(alt_seq)+2))
        elif ty == 'DEL':
            return list(range(tx_idx_h[0],tx_idx_h[0]+len(alt_seq)))
        else:
            return tx_idx_h

    def GetAlignIdxTx(tx_seq_h,tx_seq_m, align_idx_tx_h, align_idx_tx_m):
        align_tx_seq_h, align_tx_seq_m = '',''
        i,j = 0,0
        while ((i<len(tx_seq_h)) or (j<len(tx_seq_m))):
            if ((i<len(tx_seq_h)) and (j<len(tx_seq_m))):
                if i in align_idx_tx_h: #h配
                    if j in align_idx_tx_m:
                        align_tx_seq_h = align_tx_seq_h + tx_seq_h[i]
                        align_tx_seq_m = align_tx_seq_m + tx_seq_m[j]
                        i,j = i+1, j+1
                    else: #m不配
                        align_tx_seq_h = align_tx_seq_h + '-'
                        align_tx_seq_m = align_tx_seq_m + tx_seq_m[j]
                        j = j+1
                else:
                    align_tx_seq_h = align_tx_seq_h + tx_seq_h[i]
                    align_tx_seq_m = align_tx_seq_m + '-'
                    i = i + 1
            elif i>=len(tx_seq_h):
                align_tx_seq_m = align_tx_seq_m + tx_seq_m[j:]
                align_tx_seq_h = align_tx_seq_h + '-'*(len(tx_seq_m)-j)
                break
            else:
                align_tx_seq_h = align_tx_seq_h + tx_seq_h[i:]
                align_tx_seq_m = align_tx_seq_m + '-'*(len(tx_seq_h)-i)
                break
        return(str(align_tx_seq_h),str(align_tx_seq_m))    
    
    list_of_keys = ['gene_name_h',
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
 'HGVSp_m',
 'seq_h',
 'new_seq_h',
 'human_tx_idx',
 'human_p_idx',
 'human_new_p_idx',
 'dist_h',
 'seq_m',
 'new_seq_m_ori',
 'mouse_tx_idx_ori',
 'mouse_p_idx_ori',
 'mouse_new_p_idx_ori',
 'dist_m',
 'new_seq_m',
 'mouse_tx_idx',
 'mouse_p_idx',
 'mouse_new_p_idx']

    model_result[0] = dict(zip(list_of_keys,list(model_result[0].values())))
    
    if ('Non_Coding' in model_result[0]['classification_h']) or ('Splice_Site' in model_result[0]['classification_h']):
        # if Non_Coding
        # there is no alternative modeling
        # for h, h_m, m, m_m, output flank-sized mut seq with print-size

        flank_size = flank_size
        print_size = print_size

        #######read in
        # get all seq
        strand_h, strand_m = model_result[0]['strand_h'], model_result[0]['strand_m']
        tx_seq_h, tx_seq_m = model_result[0]['seq_h'], model_result[0]['seq_m']
        mut_seq_h, mut_seq_m = model_result[0]['new_seq_h'], model_result[0]['new_seq_m']
        # get ori idex
        tx_idx_h_ori, tx_idx_m_ori = model_result[0]['human_tx_idx'], model_result[0]['mouse_tx_idx']
        ref_seq_h, alt_seq_h = model_result[0]['ref_seq_h'], model_result[0]['alt_seq_h']
        ref_seq_m, alt_seq_m = model_result[0]['ref_seq_m'], model_result[0]['alt_seq_m']
        dist_h, dist_m = model_result[0]['dist_h'], model_result[0]['dist_m']

        def relative_loc(distance, ref_loc, x):
            return x+distance-ref_loc
        
        def GetAlt(strand, alt_seq, ref_seq):
            if strand == "+":
                return([alt_seq, ref_seq])
            if strand == "-":
                alt_r = str(Bio.Seq.Seq(alt_seq).reverse_complement())
                ref_r = str(Bio.Seq.Seq(ref_seq).reverse_complement())
                return([alt_r, ref_r])
        alt_seq_h, ref_seq_h = GetAlt(strand_h, alt_seq_h, ref_seq_h)
        alt_seq_m, ref_seq_m = GetAlt(strand_m, alt_seq_m, ref_seq_m)
        ty_h, ty_m = model_result[0]['type_h'], model_result[0]['type_m']
        name_h,name_m, tx_change_h,  tx_change_m = model_result[0]['gene_name_h'], model_result[0]['gene_name_m'], model_result[0]['HGVSc_h'], model_result[0]['HGVSc_m']
        class_h, class_m = model_result[0]['classification_h'], model_result[0]['classification_m']
        ######## get the original seqs
        # print alignment of h&m WT seqs
        list_idx_h, print_ori_seq_h_short, color_tx_h_ori, list_idx_m, print_ori_seq_m_short, color_tx_m_ori = get_print_elements_ori(tx_seq_h, tx_seq_m, mut_seq_h, mut_seq_m,
                                                                                                                                    tx_idx_h_ori, tx_idx_m_ori,
                                                                                                                                    flank_size, print_size)

        list_idx_h_mut, print_mut_seq_h_short, color_tx_h_mut = get_print_elements_mut(list_idx_h, print_ori_seq_h_short, color_tx_h_ori, tx_idx_h_ori, ref_seq_h, alt_seq_h, ty_h, flank_size, print_size)
        list_idx_m_mut, print_mut_seq_m_short, color_tx_m_mut = get_print_elements_mut(list_idx_m, print_ori_seq_m_short, color_tx_m_ori, tx_idx_m_ori, ref_seq_m, alt_seq_m, ty_m, flank_size, print_size)

        # plot for non-coding
        fig, ax = plt.subplots(figsize=(0.7*len(print_ori_seq_h_short)+12.6, 8), dpi = 400)

        # ratio of height
        cell_ratio = 4/3

        # 绘制蛋白质序列的颜色块和文本
        # for i, residue in enumerate(protein_seq):
        #     rect = patches.Rectangle((i*cell_width, 0), cell_width, cell_height, linewidth=15,
        #                              edgecolor='white', facecolor=color_tx_h_ori.get(residue, 'lightgrey'))
        #     ax.add_patch(rect)
        #     ax.text(i*cell_width + cell_width/2, cell_height/2, r'{}$^{{35{}}}$'.format(residue, i+3),
        #             horizontalalignment='center', verticalalignment='center', fontsize=90, fontname='Avenir')

        # for i, residue in enumerate(print_ori_seq_h_short):
        #     ax.text(i*cell_width + cell_width/2, cell_height/2, r'{}$^{{35{}}}$'.format(residue, i+3),
        #                 horizontalalignment='center', verticalalignment='center', fontsize=90, fontname='Avenir')

                
        # human original seq
        for i, nucleotide in enumerate(print_ori_seq_h_short):
            idx = list_idx_h[i]
            rect = patches.Rectangle((i, .5*cell_ratio), 1, cell_ratio, linewidth=5,
                                    edgecolor='white', facecolor=color_tx_h_ori.get(idx, 'lightgrey'))
            ax.add_patch(rect)
            ax.text((i + 0.5), cell_ratio, nucleotide,
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
            if idx!='-':
                idx =str(int(relative_loc(dist_h, tx_idx_h_ori[0], int(idx))))
            ax.text((i + 0.5), cell_ratio*1.75, idx,
                    horizontalalignment='center', verticalalignment='center', fontsize=20, fontname='arial')
            
        # human mutated seq
        for i, nucleotide in enumerate(print_mut_seq_h_short):
            idx = list_idx_h_mut[i]
            rect = patches.Rectangle((i, 2*cell_ratio), 1, cell_ratio, linewidth=5,
                                    edgecolor='white', facecolor=color_tx_h_mut.get(idx, 'lightgrey'))
            ax.add_patch(rect)
            ax.text((i + 0.5), cell_ratio*2.5, nucleotide,
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
            if idx!='-':
               idx =str(int(relative_loc(dist_h, tx_idx_h_ori[0], int(idx))))
            ax.text((i + 0.5), cell_ratio*3.25, idx,
                    horizontalalignment='center', verticalalignment='center', fontsize=20, fontname='arial')
            
        # mouse
        for i, nucleotide in enumerate(print_ori_seq_m_short):
            idx = list_idx_m[i]

            rect = patches.Rectangle((i, -1.5*cell_ratio), 1, cell_ratio, linewidth=5,
                                    edgecolor='white', 
                                    facecolor=color_tx_m_ori.get(idx, 'lightgrey'))
            ax.add_patch(rect)
            ax.text((i + 0.5), -cell_ratio, nucleotide,
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
            if idx!='-':
                idx =str(int(relative_loc(dist_m, tx_idx_m_ori[0], int(idx))))
            ax.text((i + 0.5), -cell_ratio*1.75, idx,
                    horizontalalignment='center', verticalalignment='center', fontsize=20, fontname='arial')
            
        # mouse mutated seq
        for i, nucleotide in enumerate(print_mut_seq_m_short):
            idx = list_idx_m_mut[i]

            rect = patches.Rectangle((i, -3*cell_ratio), 1, cell_ratio, linewidth=5,
                                    edgecolor='white', facecolor=color_tx_m_mut.get(idx, 'lightgrey'))
            ax.add_patch(rect)
            ax.text((i + 0.5), -cell_ratio*2.5, nucleotide,
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
            if idx!='-':
                idx =str(int(relative_loc(dist_m, tx_idx_m_ori[0], int(idx))))
            ax.text((i + 0.5), -cell_ratio*3.25, idx,
                    horizontalalignment='center', verticalalignment='center', fontsize=20, fontname='arial')

        i = i+1
        ### row annotation
        #human WT
        rect = patches.Rectangle((i+1, .5*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                    edgecolor='white', facecolor=dict_color_row['WT'])
        ax.add_patch(rect)
        ax.text((i + 2.5), cell_ratio, 'WT',
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')


        #human MUT
        rect = patches.Rectangle((i+1, 2*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                    edgecolor='white', facecolor=dict_color_row['MUT'])
        ax.add_patch(rect)
        ax.text((i + 2.5), cell_ratio*2.5, 'MUT',
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')

        #mouse WT
        rect = patches.Rectangle((i+1, -1.5*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                    edgecolor='white', facecolor=dict_color_row['WT'])
        ax.add_patch(rect)
        ax.text((i + 2.5), -cell_ratio, 'WT',
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')

        # mouse MUT
        rect = patches.Rectangle((i+1, -3*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                    edgecolor='white', facecolor=dict_color_row['MUT'])
        ax.add_patch(rect)
        ax.text((i + 2.5), -cell_ratio*2.5, 'MUT',
                    horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')

        # human annotation
        ##name
        ax.text(-6, cell_ratio*1.75, f'Human {name_h}\n\n{tx_change_h}',fontweight='bold', color = '#33608B',
                    horizontalalignment='center', verticalalignment='center', fontsize=25, fontname='arial')

        # mouse annotation
        ##name & HGVsc
        ax.text(-6, -cell_ratio*1.75, f'Mouse {name_m}\n\n{tx_change_m}',fontweight='bold', color = '#E1854C',
                    horizontalalignment='center', verticalalignment='center', fontsize=25, fontname='arial')

        # legend
        # flank size
        for loc, color, txt in zip([0, 6, 14], ['lightblue', 'pink','#A26EC4'],['Flank Size','Reference Sequence','Alternate Sequence']):
                rect = patches.Rectangle((loc, -4.5*cell_ratio), 1/2*cell_ratio, 1/2*cell_ratio, linewidth=15,
                                        edgecolor='white', facecolor=color)
                ax.add_patch(rect)
                ax.text(loc+1, -4.25*cell_ratio, txt,
                        horizontalalignment='left', verticalalignment='center', fontsize=20, fontname='arial')

        ax.axis('off')

        ax.set_xlim(-12, max(i+1+5, 18))
        ax.set_ylim(-4.8*cell_ratio, 3.8*cell_ratio)

        plt.tight_layout()
        plt.show()

    else:
        # flank_size, print_size = flank_size, print_size
        if (model_result[0]['class'] == 1) or (len(model_result)>1): # alternative modeling
            tx_seq_h, tx_seq_m = model_result[0]['seq_h'], model_result[0]['seq_m']
            mut_seq_h, mut_seq_m = model_result[0]['new_seq_h'], model_result[0]['new_seq_m_ori']
            class_h, class_m = model_result[0]['classification_h'], model_result[0]['classification_m']
            # translate
            p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m = str(Translate(tx_seq_h)), str(Translate(tx_seq_m)), str(Translate(mut_seq_h)), str(Translate(mut_seq_m))
            stop_loc_p_h, stop_loc_p_m, stop_mut_loc_p_h, stop_mut_loc_p_m = p_seq_h.index('*'), p_seq_m.index('*'), mut_p_seq_h.index('*'), mut_p_seq_m.index('*')
            if class_h != 'Non_stop':
                p_seq_h, p_seq_m = p_seq_h[:stop_loc_p_h+1], p_seq_m[:stop_loc_p_m+1]
                if class_h == 'Nonsense':
                    mut_p_seq_h, mut_p_seq_m = mut_p_seq_h[:stop_loc_p_h+1], mut_p_seq_m[:stop_loc_p_m+1]
                else:
                    mut_p_seq_h, mut_p_seq_m = mut_p_seq_h[:stop_mut_loc_p_h+1], mut_p_seq_m[:stop_mut_loc_p_m+1]
            # get ori idex
            tx_idx_h_ori, tx_idx_m_ori = model_result[0]['human_tx_idx'], model_result[0]['mouse_tx_idx_ori']
            p_idx_h_ori, p_idx_h_new, p_idx_m_ori, p_idx_m_new = model_result[0]['human_p_idx'], model_result[0]['human_new_p_idx'], model_result[0]['mouse_p_idx_ori'], model_result[0]['mouse_new_p_idx_ori']

            ref_seq_h, alt_seq_h = model_result[0]['ref_seq_h'], model_result[0]['alt_seq_h']
            ref_seq_m, alt_seq_m = model_result[0]['ref_seq_m_ori'], model_result[0]['alt_seq_m_ori']
            ty_h, ty_m = model_result[0]['type_h'], model_result[0]['type_m']
            name_h,name_m, tx_change_h,  tx_change_m = model_result[0]['gene_name_h'], model_result[0]['gene_name_m'], model_result[0]['HGVSc_h'], model_result[0]['HGVSc_m_ori']
            p_change_h, p_change_m = model_result[0]['HGVSp_h'], model_result[0]['HGVSp_m_ori']

            if class_h != 'Non_Stop':
                list_idx_h, print_ori_seq_h_short, color_p_h_ori, list_idx_m, print_ori_seq_m_short, color_p_m_ori = get_print_elements_ori(p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m,
                                                                                                                                        p_idx_h_ori, p_idx_m_ori,
                                                                                                                                        flank_size, print_size)
            else:
                list_idx_h, print_ori_seq_h_short, color_p_h_ori, list_idx_m, print_ori_seq_m_short, color_p_m_ori = get_print_elements_ori_non_stop(p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m,
                                                                                                                                        p_idx_h_ori, p_idx_m_ori,
                                                                                                                                        flank_size, print_size)
            # print transcript seqs for peptides
            list_idx_h_tx, print_ori_tx_h, color_tx_h_ori = print_tx_for_pep(list_idx_h, tx_seq_h, tx_idx_h_ori, '#A3BEBC')
            list_idx_m_tx, print_ori_tx_m, color_tx_m_ori = print_tx_for_pep(list_idx_m, tx_seq_m, tx_idx_m_ori, '#A3BEBC')

            ##### get mutated ones
            # class 1: SNV/Nonsense
            list_idx_h_mut, print_mut_seq_h_short, color_p_h_mut = get_print_elements_mut_pep(list_idx_h, print_ori_seq_h_short, color_p_h_ori, mut_p_seq_h ,p_idx_h_new, ty_h, class_h, flank_size, print_size)
            list_idx_m_mut, print_mut_seq_m_short, color_p_m_mut = get_print_elements_mut_pep(list_idx_m, print_ori_seq_m_short, color_p_m_ori, mut_p_seq_m, p_idx_m_new, ty_m, class_m, flank_size, print_size)

            
            # print transcript seqs for peptides
            list_idx_h_tx_mut, print_mut_tx_h, color_tx_h_mut = print_tx_for_pep(list_idx_h_mut, mut_seq_h, change_to_mutated(tx_idx_h_ori, ty_h, alt_seq_h, ref_seq_h), '#C05252')
            list_idx_m_tx_mut, print_mut_tx_m, color_tx_m_mut = print_tx_for_pep(list_idx_m_mut, mut_seq_m, change_to_mutated(tx_idx_m_ori, ty_m, alt_seq_m, ref_seq_m), '#C05252')

            fig, ax = plt.subplots(figsize=(1.3*len(print_ori_seq_h_short)+14, 9.6+4.7*len(model_result)), dpi = 400)
            # ratio of height
            cell_ratio = 4/3

            ################ peps
            # in the order of human-mut, huaman-wt, mouse-wt, mouse-mut
            for list_idx, seq, color_dict, height in zip(
            # list_idx
                [list_idx_h_mut, list_idx_h, list_idx_m, list_idx_m_mut],
                [print_mut_seq_h_short, print_ori_seq_h_short, print_ori_seq_m_short, print_mut_seq_m_short],
                [color_p_h_mut, color_p_h_ori, color_p_m_ori, color_p_m_mut],
                [7/2, 1/2, -3/2, -9/2]
            ):
                # human-wt-pep

                for i, pep in enumerate(seq):
                    idx = list_idx[i]
                    rect = patches.Rectangle((3*i, height*cell_ratio), 3, cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=color_dict.get(idx, 'lightgrey'))
                    ax.add_patch(rect)
                    if idx!='-':
                        idx =str(int(idx)+1)
                        txt = r'{}$^{{{}}}$'.format(pep, idx)
                    else:
                        txt = '-'
                    ax.text(3*(i + 0.5), (height+.4)*cell_ratio, txt,
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                

            ################ transcripts
            #human-mut-tx
            for list_idx, seq, color_dict, height in zip(
            # list_idx
                [list_idx_h_tx_mut, list_idx_h_tx, list_idx_m_tx, list_idx_m_tx_mut],
                [print_mut_tx_h, print_ori_tx_h, print_ori_tx_m, print_mut_tx_m],
                [color_tx_h_mut, color_tx_h_ori, color_tx_m_ori, color_tx_m_mut],
                [5/2, 3/2, -5/2, -7/2]
            ):

                for i, neucleotide in enumerate(seq):
                    idx = list_idx[i]
                    rect = patches.Rectangle((i, height*cell_ratio), 1, cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=color_dict.get(idx, 'white'))
                    ax.add_patch(rect)
                    ax.text((i + 0.5), (height+.5)*cell_ratio, neucleotide,
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')

        # human annotation
            ##name
            ax.text(-9, cell_ratio*2.5, f'Human {name_h}\n\n{tx_change_h}\n\n{p_change_h}',fontweight='bold', color = '#33608B',
                        horizontalalignment='center', verticalalignment='center', fontsize=25, fontname='arial')

            # mouse annotation
            ##name & HGVsc
            ax.text(-9, -cell_ratio*2.5, f'Mouse {name_m}\n\n{tx_change_m}\n\n{p_change_m}',fontweight='bold', color = '#E1854C',
                        horizontalalignment='center', verticalalignment='center', fontsize=25, fontname='arial')
                    
            # WT
            for height in [1/2, -3/2]:
                rect = patches.Rectangle((i+1, height*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=dict_color_row['WT'])
                ax.text((i + 2.5), (height+.5)*cell_ratio, 'WT',
                        horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                ax.add_patch(rect)

            # MUT
            for height in [7/2,-9/2]:
                rect = patches.Rectangle((i+1, height*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=dict_color_row['ORI'])
                ax.add_patch(rect)
                ax.text((i + 2.5), (height+.5)*cell_ratio, 'ORI',
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')

            h = -9/2 -1

            for r in range(len(model_result)):
                tx_seq_h, tx_seq_m = model_result[r]['seq_h'], model_result[r]['seq_m']
                mut_seq_h, mut_seq_m = model_result[r]['new_seq_h'], model_result[r]['new_seq_m']
                class_h, class_m = model_result[r]['classification_h'], model_result[r]['classification_m']
                # translate
                p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m = str(Translate(tx_seq_h)), str(Translate(tx_seq_m)), str(Translate(mut_seq_h)), str(Translate(mut_seq_m))
                stop_loc_p_h, stop_loc_p_m, stop_mut_loc_p_h, stop_mut_loc_p_m = p_seq_h.index('*'), p_seq_m.index('*'), mut_p_seq_h.index('*'), mut_p_seq_m.index('*')
                if class_h != 'Non_stop':
                    p_seq_h, p_seq_m = p_seq_h[:stop_loc_p_h+1], p_seq_m[:stop_loc_p_m+1]
                    if class_h == 'Nonsense':
                        mut_p_seq_h, mut_p_seq_m = mut_p_seq_h[:stop_loc_p_h+1], mut_p_seq_m[:stop_loc_p_m+1]
                    else:
                        mut_p_seq_h, mut_p_seq_m = mut_p_seq_h[:stop_mut_loc_p_h+1], mut_p_seq_m[:stop_mut_loc_p_m+1]
                # get ori idex
                tx_idx_h_ori, tx_idx_m_ori = model_result[r]['human_tx_idx'], model_result[r]['mouse_tx_idx']
                p_idx_h_ori, p_idx_h_new, p_idx_m_ori, p_idx_m_new = model_result[r]['human_p_idx'], model_result[r]['human_new_p_idx'], model_result[r]['mouse_p_idx'], model_result[r]['mouse_new_p_idx']

                ref_seq_h, alt_seq_h = model_result[r]['ref_seq_h'], model_result[r]['alt_seq_h']
                ref_seq_m, alt_seq_m = model_result[r]['ref_seq_m_ori'], model_result[r]['alt_seq_m']
                ty_h, ty_m = model_result[r]['type_h'], model_result[r]['type_m']
                name_h,name_m, tx_change_h,  tx_change_m = model_result[r]['gene_name_h'], model_result[r]['gene_name_m'], model_result[r]['HGVSc_h'], model_result[r]['HGVSc_m']
                p_change_h, p_change_m = model_result[r]['HGVSp_h'], model_result[r]['HGVSp_m']

                if class_h != 'Non_Stop':
                    list_idx_h, print_ori_seq_h_short, color_p_h_ori, list_idx_m, print_ori_seq_m_short, color_p_m_ori = get_print_elements_ori(p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m,
                                                                                                                                            p_idx_h_ori, p_idx_m_ori,
                                                                                                                                            flank_size,print_size)
                else:
                    list_idx_h, print_ori_seq_h_short, color_p_h_ori, list_idx_m, print_ori_seq_m_short, color_p_m_ori = get_print_elements_ori_non_stop(p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m,
                                                                                                                                            p_idx_h_ori, p_idx_m_ori,
                                                                                                                                            flank_size,print_size)
                # print transcript seqs for peptides
                list_idx_h_tx, print_ori_tx_h, color_tx_h_ori = print_tx_for_pep(list_idx_h, tx_seq_h, tx_idx_h_ori, '#A3BEBC')
                list_idx_m_tx, print_ori_tx_m, color_tx_m_ori = print_tx_for_pep(list_idx_m, tx_seq_m, tx_idx_m_ori, '#A3BEBC')

                ##### get mutated ones
                # class 1: SNV/Nonsense
                list_idx_h_mut, print_mut_seq_h_short, color_p_h_mut = get_print_elements_mut_pep(list_idx_h, print_ori_seq_h_short, color_p_h_ori, mut_p_seq_h ,p_idx_h_new, ty_h, class_h, flank_size, print_size)
                list_idx_m_mut, print_mut_seq_m_short, color_p_m_mut = get_print_elements_mut_pep(list_idx_m, print_ori_seq_m_short, color_p_m_ori, mut_p_seq_m, p_idx_m_new, ty_m, class_m, flank_size, print_size)

                # print transcript seqs for peptides
                list_idx_h_tx_mut, print_mut_tx_h, color_tx_h_mut = print_tx_for_pep(list_idx_h_mut, mut_seq_h, change_to_mutated(tx_idx_h_ori, ty_h, alt_seq_h, ref_seq_h), '#C05252')
                list_idx_m_tx_mut, print_mut_tx_m, color_tx_m_mut = print_tx_for_pep(list_idx_m_mut, mut_seq_m, change_to_mutated(tx_idx_m_ori, ty_m, alt_seq_m, ref_seq_m), '#C05252')

                ################ peps
                # in the order of human-mut, huaman-wt, mouse-wt, mouse-mut
                for list_idx, seq, color_dict, height in zip(
                # list_idx
                    [list_idx_m, list_idx_m_mut],
                    [print_ori_seq_m_short, print_mut_seq_m_short],
                    [color_p_m_ori, color_p_m_mut],
                    [h, h-3]
                ):
                    # human-wt-pep

                    for i, pep in enumerate(seq):
                        idx = list_idx[i]
                        rect = patches.Rectangle((3*i, height*cell_ratio), 3, cell_ratio, linewidth=5,
                                                edgecolor='white', facecolor=color_dict.get(idx, 'lightgrey'))
                        ax.add_patch(rect)
                        if idx!='-':
                            idx =str(int(idx)+1)
                            txt = r'{}$^{{{}}}$'.format(pep, idx)
                        else:
                            txt = '-'
                        ax.text(3*(i + 0.5), (height+.4)*cell_ratio, txt,
                                horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                

                ################ transcripts
                #human-mut-tx
                for list_idx, seq, color_dict, height in zip(
                # list_idx
                    [list_idx_m_tx, list_idx_m_tx_mut],
                    [print_ori_tx_m, print_mut_tx_m],
                    [color_tx_m_ori, color_tx_m_mut],
                    [h-1,h-2]
                ):

                    for i, neucleotide in enumerate(seq):
                        idx = list_idx[i]
                        rect = patches.Rectangle((i, height*cell_ratio), 1, cell_ratio, linewidth=5,
                                                edgecolor='white', facecolor=color_dict.get(idx, 'white'))
                        ax.add_patch(rect)
                        ax.text((i + 0.5), (height+.5)*cell_ratio, neucleotide,
                                horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                

                rect = patches.Rectangle((i+1, h*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                                edgecolor='white', facecolor=dict_color_row['WT'])
                ax.add_patch(rect)
                ax.text((i + 2.5), (h+.5)*cell_ratio, 'WT',
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                

                rect = patches.Rectangle((i+1, (h-3)*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                                edgecolor='white', facecolor=dict_color_row['ALT'])
                ax.add_patch(rect)
                ax.text((i + 2.5), (h-3+.5)*cell_ratio, 'ALT',
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                
                ##name & HGVsc
                ax.text(-9, cell_ratio*(h-1.5), f'Mouse {name_m}\n\n{tx_change_m}\n\n{p_change_m}',fontweight='bold', color = '#E1854C',
                        horizontalalignment='center', verticalalignment='center', fontsize=25, fontname='arial')
                h = h - 4

            h = h - 1
            # legend
            for loc, color, txt in zip([0,9,18,27,38], ['lightblue', 'pink','#A26EC4', '#A3BEBC', '#C05252'],['Flank Size','Reference AA','Alternate AA', 'Reference Sequence', 'Alternate Sequence']):
                    rect = patches.Rectangle((loc, (h+.5)*cell_ratio), 1/2*cell_ratio, 1/2*cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=color)
                    ax.add_patch(rect)
                    ax.text(loc+cell_ratio, (h+.75)*cell_ratio, txt,
                            horizontalalignment='left', verticalalignment='center', fontsize=20, fontname='arial')


            ax.axis('off')
            ax.set_xlim(-18, max(3*len(print_ori_seq_h_short)+1+5, 40))
            ax.set_ylim((h-.5)*cell_ratio, 6*cell_ratio)

            plt.tight_layout()
            plt.show()

        else:

            tx_seq_h, tx_seq_m = model_result[0]['seq_h'], model_result[0]['seq_m']
            mut_seq_h, mut_seq_m = model_result[0]['new_seq_h'], model_result[0]['new_seq_m']
            class_h, class_m = model_result[0]['classification_h'], model_result[0]['classification_m']
            # translate
            p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m = str(Translate(tx_seq_h)), str(Translate(tx_seq_m)), str(Translate(mut_seq_h)), str(Translate(mut_seq_m))
            stop_loc_p_h, stop_loc_p_m, stop_mut_loc_p_h, stop_mut_loc_p_m = p_seq_h.index('*'), p_seq_m.index('*'), mut_p_seq_h.index('*'), mut_p_seq_m.index('*')
            if class_h != 'Non_stop':
                p_seq_h, p_seq_m = p_seq_h[:stop_loc_p_h+1], p_seq_m[:stop_loc_p_m+1]
                if class_h == 'Nonsense':
                    mut_p_seq_h, mut_p_seq_m = mut_p_seq_h[:stop_loc_p_h+1], mut_p_seq_m[:stop_loc_p_m+1]
                else:
                    mut_p_seq_h, mut_p_seq_m = mut_p_seq_h[:stop_mut_loc_p_h+1], mut_p_seq_m[:stop_mut_loc_p_m+1]
            # get ori idex
            tx_idx_h_ori, tx_idx_m_ori = model_result[0]['human_tx_idx'], model_result[0]['mouse_tx_idx']
            p_idx_h_ori, p_idx_h_new, p_idx_m_ori, p_idx_m_new = model_result[0]['human_p_idx'], model_result[0]['human_new_p_idx'], model_result[0]['mouse_p_idx'], model_result[0]['mouse_new_p_idx']

            ref_seq_h, alt_seq_h = model_result[0]['ref_seq_h'], model_result[0]['alt_seq_h']
            ref_seq_m, alt_seq_m = model_result[0]['ref_seq_m'], model_result[0]['alt_seq_m']
            ty_h, ty_m = model_result[0]['type_h'], model_result[0]['type_m']
            name_h,name_m, tx_change_h,  tx_change_m = model_result[0]['gene_name_h'], model_result[0]['gene_name_m'], model_result[0]['HGVSc_h'], model_result[0]['HGVSc_m']
            p_change_h, p_change_m = model_result[0]['HGVSp_h'], model_result[0]['HGVSp_m']

            if class_h != 'Non_Stop':
                list_idx_h, print_ori_seq_h_short, color_p_h_ori, list_idx_m, print_ori_seq_m_short, color_p_m_ori = get_print_elements_ori(p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m,
                                                                                                                                        p_idx_h_ori, p_idx_m_ori,
                                                                                                                                        flank_size,print_size)
            else:
                list_idx_h, print_ori_seq_h_short, color_p_h_ori, list_idx_m, print_ori_seq_m_short, color_p_m_ori = get_print_elements_ori_non_stop(p_seq_h, p_seq_m, mut_p_seq_h, mut_p_seq_m,
                                                                                                                                        p_idx_h_ori, p_idx_m_ori,
                                                                                                                                        flank_size,print_size)
            # print transcript seqs for peptides
            list_idx_h_tx, print_ori_tx_h, color_tx_h_ori = print_tx_for_pep(list_idx_h, tx_seq_h, tx_idx_h_ori, '#A3BEBC')
            list_idx_m_tx, print_ori_tx_m, color_tx_m_ori = print_tx_for_pep(list_idx_m, tx_seq_m, tx_idx_m_ori, '#A3BEBC')

            ##### get mutated ones
            # class 1: SNV/Nonsense
            list_idx_h_mut, print_mut_seq_h_short, color_p_h_mut = get_print_elements_mut_pep(list_idx_h, print_ori_seq_h_short, color_p_h_ori, mut_p_seq_h ,p_idx_h_new, ty_h, class_h, flank_size, print_size)
            list_idx_m_mut, print_mut_seq_m_short, color_p_m_mut = get_print_elements_mut_pep(list_idx_m, print_ori_seq_m_short, color_p_m_ori, mut_p_seq_m, p_idx_m_new, ty_m, class_m, flank_size, print_size)

            # print transcript seqs for peptides
            list_idx_h_tx_mut, print_mut_tx_h, color_tx_h_mut = print_tx_for_pep(list_idx_h_mut, mut_seq_h, change_to_mutated(tx_idx_h_ori, ty_h, alt_seq_h, ref_seq_h), '#C05252')
            list_idx_m_tx_mut, print_mut_tx_m, color_tx_m_mut = print_tx_for_pep(list_idx_m_mut, mut_seq_m, change_to_mutated(tx_idx_m_ori, ty_m, alt_seq_m, ref_seq_m), '#C05252')

            fig, ax = plt.subplots(figsize=(1.3*len(print_ori_seq_h_short)+10, 9.6), dpi = 400)
            # ratio of height
            cell_ratio = 4/3

            ################ peps
            # in the order of human-mut, huaman-wt, mouse-wt, mouse-mut
            for list_idx, seq, color_dict, height in zip(
            # list_idx
                [list_idx_h_mut, list_idx_h, list_idx_m, list_idx_m_mut],
                [print_mut_seq_h_short, print_ori_seq_h_short, print_ori_seq_m_short, print_mut_seq_m_short],
                [color_p_h_mut, color_p_h_ori, color_p_m_ori, color_p_m_mut],
                [7/2, 1/2, -3/2, -9/2]
            ):
                # human-wt-pep

                for i, pep in enumerate(seq):
                    idx = list_idx[i]
                    rect = patches.Rectangle((3*i, height*cell_ratio), 3, cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=color_dict.get(idx, 'lightgrey'))
                    ax.add_patch(rect)
                    if idx!='-':
                        idx =str(int(idx)+1)
                        txt = r'{}$^{{{}}}$'.format(pep, idx)
                    else:
                        txt = '-'
                    ax.text(3*(i + 0.5), (height+.4)*cell_ratio, txt,
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                

            ################ transcripts
            #human-mut-tx
            for list_idx, seq, color_dict, height in zip(
            # list_idx
                [list_idx_h_tx_mut, list_idx_h_tx, list_idx_m_tx, list_idx_m_tx_mut],
                [print_mut_tx_h, print_ori_tx_h, print_ori_tx_m, print_mut_tx_m],
                [color_tx_h_mut, color_tx_h_ori, color_tx_m_ori, color_tx_m_mut],
                [5/2, 3/2, -5/2, -7/2]
            ):

                for i, neucleotide in enumerate(seq):
                    idx = list_idx[i]
                    rect = patches.Rectangle((i, height*cell_ratio), 1, cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=color_dict.get(idx, 'white'))
                    ax.add_patch(rect)
                    ax.text((i + 0.5), (height+.5)*cell_ratio, neucleotide,
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')

        # human annotation
        ##name
            ax.text(-9, cell_ratio*2.5, f'Human {name_h}\n\n{tx_change_h}\n\n{p_change_h}',fontweight='bold', color = '#33608B',
                        horizontalalignment='center', verticalalignment='center', fontsize=25, fontname='arial')

            # mouse annotation
            ##name & HGVsc
            ax.text(-9, -cell_ratio*2.5, f'Mouse {name_m}\n\n{tx_change_m}\n\n{p_change_m}',fontweight='bold', color = '#E1854C',
                        horizontalalignment='center', verticalalignment='center', fontsize=25, fontname='arial')


            # legend
            for loc, color, txt in zip([0,9,18,27,38], ['lightblue', 'pink','#A26EC4', '#A3BEBC', '#C05252'],['Flank Size','Reference AA','Alternate AA', 'Reference Sequence', 'Alternate Sequence']):
                    rect = patches.Rectangle((loc, -6*cell_ratio), 1/2*cell_ratio, 1/2*cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=color)
                    ax.add_patch(rect)
                    ax.text(loc+cell_ratio, -5.75*cell_ratio, txt,
                            horizontalalignment='left', verticalalignment='center', fontsize=20, fontname='arial')

            i = len(print_ori_seq_h_short)*3
                    
            # WT
            for height in [1/2, -3/2]:
                rect = patches.Rectangle((i+1, height*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=dict_color_row['WT'])
                ax.text((i + 2.5), (height+.5)*cell_ratio, 'WT',
                        horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')
                ax.add_patch(rect)

            # MUT
            for height in [7/2,-9/2]:
                rect = patches.Rectangle((i+1, height*cell_ratio), 3, 1*cell_ratio, linewidth=5,
                                            edgecolor='white', facecolor=dict_color_row['MUT'])
                ax.add_patch(rect)
                ax.text((i + 2.5), (height+.5)*cell_ratio, 'MUT',
                            horizontalalignment='center', verticalalignment='center', fontsize=30, fontname='arial')

            ax.axis('off')
            ax.set_xlim(-18, max(3*len(print_ori_seq_h_short)+1+5, 40))
            ax.set_ylim(-6.5*cell_ratio, 6*cell_ratio)

            plt.tight_layout()
            plt.show()
