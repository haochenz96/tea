from pathlib import Path
import time
# from missionbio.h5.data import Assay
from h5.data import Assay
import numpy as np
import pandas as pd
from tea.parse import *
import sys
from tea.format import isNaN, check_matrix_format, CONDENSED_SNV_FORMAT
import time


NONFUNC_SO = ['2kb_upstream_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'intron_variant', 'synonymous_variant', ]

# ----- (1) functions to query from opencravat.org -----
def write_cravat_input(variants, cravat_input_path, sample_name, *additional_var_info):
    '''
    Writes a cravat input file for a given sample name
    
        inputs:
        - varaints: list of variants to be analyzed (in Tapestri output format)
        - cravat_input_path: path to write cravat input file
        - sample_name: sample name to be used for CRAVAT query
        - *additional_var_info: one-column dataframes with variant_name as index, and additional information in the first column

        outputs:

    '''

    # unpack additional_var_info
    additional_var_info = list(additional_var_info)

    if not isinstance(cravat_input_path, Path):
        cravat_input_path = Path(cravat_input_path)

    print(f"[INFO] Writing CRAVAT input file to -- {cravat_input_path}...")
    cravat_input_path.parents[0].mkdir(
        parents=True, 
        exist_ok=True
    )
    if not cravat_input_path.is_file():
        with open(cravat_input_path, 'w') as f:

            var_count=0
            for var in variants:

                var_count += 1
                chrom = var.split(':')[0]
                pos = var.split(':')[1]
                ref = var.split(':')[2].split('/')[0]
                alt = var.split(':')[2].split('/')[1]
                if len(additional_var_info) > 0:
                    f.write(f'{chrom} {pos} + {ref} {alt} var_{var_count} {" ".join([str(additional_var_info[i][var]) for i in range(len(additional_var_info))])}\n')
                else:
                    f.write(f'{chrom} {pos} + {ref} {alt} var_{var_count}\n')
                    
        print(f"[INFO] Finished writing CRAVAT input file!")
    else:
        print(f"[WARNING] CRAVAT input file already exists!")

# sys.setrecursionlimit(100) # <---- will cause issue for loading plotly
# print('[INFO] Recursion limit set to 100')

def get_cravat_output(session, job_id, output_path, max_query_time = 120, query_times = 0):
    '''
    Query the current session to check CRAVAT job status. If job is finished, download the output file to output_path.
    
    Params:
    ----------
    - session: requests.session where the cravat query is in
    - job_id: can be feteched by post.json()['id']
    - output_path: path to write the CRAVAT output file
    - max_query_time: maximum time to query for the output file (unit: seconds)
    - query_times: number of times the query has been attempted; for keeping track of recursion.

    Returns:
    ----------
    None

    '''
    if query_times >= int(max_query_time / 10):
        raise RuntimeError(f'[ERROR] Maximum query time of {max_query_time} seconds exceeded!')

    if not isinstance(output_path, Path):
        output_path = Path(output_path)
    output_path.parents[0].mkdir(parents=True, exist_ok=True)
    response = session.get(
        'https://run.opencravat.org/submit/jobs/' + job_id + '/status')
    if response.json()['status'] == 'error':
        print(f'[WARNING] CRAVAT run failed! Check input file')
    elif response.json()['status'] == 'Finished':
        output = session.get('https://run.opencravat.org/submit/jobs/' + job_id + '/reports/text')

        with open(output_path, 'w') as f:
            f.write(output.text)    
    elif response.json()['status'] != 'Finished':
        time.sleep(10) # wait 10 seconds before checking again
        get_cravat_output(session, job_id, output_path, max_query_time, query_times + 1)
    else:
        print(f'[WARNING] CRAVAT run failed! Unknown error')

def clean_and_format_cravat_df(cravat_df, fill_na = False):
    '''
    Clean and format the cravat output dataframe
    '''

    # --- format multi-layer index
    cravat_df = adjust_multiindex(cravat_df)

    # --- create Tapestri format
    for i, row in cravat_df.iterrows():
        tap = row['Original Input']['Chrom'] + ':' + str(row['Original Input']['Pos']) + ':' + row['Original Input']['Reference allele'] + '/' + row['Original Input']['Alternate allele']
        cravat_df.loc[i, ('Tapestri_format', 'index')] = tap
        
    cravat_df.set_index(('Tapestri_format', 'index'), inplace=True)
    cravat_df.index.rename('condensed_format', inplace=True)

    # --- rearrange annotation order
    original_input_cols = cravat_df.pop('Original Input')
    for i in range(original_input_cols.shape[1]):
        cravat_df.insert(i, ('Original Input', original_input_cols.columns[i]), original_input_cols.iloc[:,i])

    # 2. drop useless columns
    # (1) drop the hg38 annotations
    for i in ['UID','Chrom','Position','Ref Base','Alt Base','Note', 'All Mappings', 'Sample Count', 'Samples']:
        if i in cravat_df['Variant Annotation'].columns:
            try:
                cravat_df.drop(columns=('Variant Annotation', i), inplace=True)
            except KeyError:
                print(f'[WARNING] Variant Annotation - {i} column does not exist!')
    # (2) OncoKB is not functional yet, so drop
    try:
        cravat_df.drop(columns='OncoKB', inplace=True)
    except KeyError:
        print(f'[WARNING] OncoKB column does not exist!')

    # (3) CHASM's Score, Transcript and trailing info 
    for i in ['Score', 'Transcript', 'All Annotations']:
        CHASM_cols = [j for j in cravat_df.columns.levels[0] if 'CHASM' in j]
        for k in CHASM_cols:
            if i in cravat_df[k].columns:
                cravat_df.drop(columns=(k, i), inplace=True)

    # 3. rename the single-cell mutational prevalence column
    Tapestri_sc_mut_prev = cravat_df.pop(('Variant Annotation', 'Tags'))
    cravat_df.insert(4, ('Tapestri_result', 'sc_mut_prev'), Tapestri_sc_mut_prev)

    # 4. fill N/A's
    if fill_na:
        cravat_df = cravat_df.fillna('')

    return cravat_df

def create_ann_map_from_cravat_df(cravat_df, voi=None):
    '''
    Create a dictionary of SNV from the cravat output dataframe based on the following rule:
        if protein change is available, use protein change
        if protein change is not available, use sequence ontology

    inputs:
    - cravat_df: cleaned, formatted cravat output dataframe

    outputs:
    - ann_map: dictionary of annotations
    - ann_map_df: dataframe of annotations

    '''
    ann_map = {}
    if voi is None:
        voi = cravat_df.index
    for var_i, row in cravat_df.iterrows():
        if not var_i in voi:
            continue
        if not var_i in ann_map:
            if not isNaN(row['Variant Annotation']['Protein Change']): # coding variant
                ann_map[var_i] = row['Variant Annotation']['Gene'] + ' ' + row['Variant Annotation']['Protein Change']
            else:
                ann_map[var_i] = row['Variant Annotation']['Gene'] + ' ' + row['Variant Annotation']['Sequence Ontology']
    ann_map_df = pd.DataFrame.from_dict(ann_map, orient='index')
    ann_map_df.columns = ['variant annotation']
    ann_map_df.index.rename('condensed_format', inplace=True)
    
    return ann_map, ann_map_df

def get_technical_artifact_mask(
    voi,
    cravat_df, 
    num_cells, 
    mut_prev_series,
    bq_prev_threshold = 0.005, 
    normals_pon_occurence=4, 
    rescue_1000genome_af = 0.01, 
    filter_broad_wes_pon = False,
    ado_threshold = None,
    ngt_df = None,
    filter_TtoC_artifact = None,
    blacklist = None
    ):
    '''
    Create a mask for technical artifacts

    inputs:
    - voi: list of variants of interest
    - cravat_df: cleaned, formatted cravat output dataframe
    - mut_prev_series: series of mutational prevalence for each SNV
    - num_cells: number of cells in the sample
    - bq_prev_threshold: base quality threshold for filtering. If none, will not filter. If yes, make sure that it matches up with mut_prev_series
    - normals_pon_occurence: number of normals that the SNV must be present in
    - rescue_1000genome: whether to rescue recurrent PoN SNVs that are in 1000 genome
    - filter_broad_wes_pon: whether to filter based on broad_wes_pon
    - ado_threshold: filter SNVs that have 
    - ngt_df: for getting high ADO SNVs. Should be single-cells x SNVs
    - filter_TtoC_artifact: filter T>C artifacts
    - blacklist: list of variants to be blacklisted

    outputs:
    - mask: boolean mask for technical artifacts

    '''
    try:
        cravat_df = cravat_df.loc[voi]
    except KeyError:
        raise KeyError('Some SNVs are not found in CRAVAT df!')
    # ----- base quality mask -----
    if bq_prev_threshold is None:
        print('[WARNING] No base quality threshold is given. Skipping filtering...')
        bq_mask = pd.Series(True, index=voi)
    else:
        try:
            bq_mask = ~(cravat_df[('blacklist_comparison', 'blacklist-base_qual-sc_prev')] >= mut_prev_series) 
            if num_cells is not None:
                print(f'[INFO] Filtering out variants with base quality flag in >= {bq_prev_threshold*num_cells}/{num_cells} single cells')
                bq_mask = bq_mask & ~(cravat_df[('blacklist_comparison', 'blacklist-base_qual-sc_prev')] >= bq_prev_threshold*num_cells) 
            else:
                raise ValueError('Number of cells is required to filter based on base quality threshold')

        except KeyError:
            print('[WARNING] Base quality threshold given but no blacklist-base_qual-sc_prev column is found in the CRAVAT file. Skipping filtering...')
            bq_mask = pd.Series(True, index=voi)
            
    # ----- PoN mask, potentially rescuing 1000 genome variants -----
    # find the normals_occurence column
    normals_pon_occurence_col = None
    for pon_col in cravat_df.columns.levels[1]:
        if 'normals-occurence' in pon_col:
            normals_pon_occurence_col = pon_col
            break

    if normals_pon_occurence_col is None:
        print('[WARNING] No normals-occurence column is found. Skipping PoN filtering...')
        pon_mask = pd.Series(True, index=voi)
    else:
        pon_mask = ~(cravat_df[('PoN_comparison',normals_pon_occurence_col)] >= normals_pon_occurence)
        if rescue_1000genome_af is not None:
            try:
                pon_mask = pon_mask | (cravat_df[('1000 Genomes', 'AF')] > rescue_1000genome_af)
            except KeyError:
                print('[WARNING] 1000 Genomes column not found! Skipping rescuing...')
            
        if filter_broad_wes_pon:
            try:
                pon_mask = pon_mask & isNaN(cravat_df[('bulk_comparison', 'bulk-broad_wes_pon-AF')])
            except KeyError:
                print('[WARNING] broad_wes_pon column not found! Skipping filtering...')

    # ----- ADO mask -----
    ado_mask = pd.Series(True, index=voi)
    if ado_threshold is not None:
        try:
            ngt_df = ngt_df[voi]
        except KeyError:
            raise KeyError('Some SNVs are not found in ngt_df!')
        ado_high_vars = voi[
            (ngt_df == 3).sum(axis=0) > (ado_threshold*num_cells)
        ]

        print(f"[INFO] Filtering {len([v for v in voi if v in ado_high_vars])} SNVs due to ADO > {ado_threshold}")
        ado_mask[ado_high_vars] = False

    # ----- T>C variants mask -----
    tc_mask = pd.Series(True, index=voi)
    if filter_TtoC_artifact is not None:
        try:
            filter_TtoC_artifact_binary = filter_TtoC_artifact['filter'] 
            filter_TtoC_artifact_lower_thres = filter_TtoC_artifact['lower_thres']
            filter_TtoC_artifact_upper_thres = filter_TtoC_artifact['upper_thres']
        except KeyError:
            raise ValueError(f"[ERROR] filter_TtoC_artifact must have keys ['filter', 'lower_thres', 'upper_thres']")

        if filter_TtoC_artifact_binary:
            tc_mask = (cravat_df.index.str.endswith('T/C')) & \
                    (mut_prev_series >= filter_TtoC_artifact_lower_thres * num_cells) & \
                    (mut_prev_series <= filter_TtoC_artifact_upper_thres * num_cells)
            print(f"[INFO] There are {tc_mask.sum()} T>C artifacts that are rare. They will be filtered out.")
            tc_mask = ~tc_mask

    # ----- manual blacklist mask -----
    blacklist_mask = pd.Series(True, index=voi)
    if blacklist is not None:
        if type(blacklist) is list:
            blacklist_snvs = set(blacklist)
        else:
            try:
                blacklist_snv_df = pd.read_csv(blacklist, index_col=0)
                if not check_matrix_format(blacklist_snv_df, CONDENSED_SNV_FORMAT):
                    raise ValueError(f"[ERROR] blacklist_snvs file not in the correct format (index column needs to be in condensed SNV format).")
                blacklist_snvs = set(blacklist_snv_df.index)
            except:
                raise ValueError(f"[ERROR] blacklist_snvs should either be a list of SNVs or a path to a CSV file whose index is the list of SNVs.")
        blacklist_mask.loc[[v for v in blacklist_snvs if v in blacklist_mask.index]] = False

    return bq_mask & pon_mask & ado_mask & tc_mask & blacklist_mask