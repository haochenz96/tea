from pathlib import Path
import time
# from missionbio.h5.data import Assay
from h5.data import Assay
import numpy as np
import pandas as pd
from tea.parse import *
import sys
# ----- (1) functions to query from opencravat.org -----
def write_cravat_input(variants, cravat_input_path, sample_name, *additional_var_info):
    '''
    Writes a cravat input file for a given sample name
    
        inputs:
        - varaints: list of variants to be analyzed (in Tapestri output format)
        - cravat_input_path: path to write cravat input file
        - sample_name: sample name to be used for CRAVAT query

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

                f.write(f'{chrom} {pos} + {ref} {alt} var_{var_count} {" ".join([str(additional_var_info[i][var]) for i in range(len(additional_var_info))])}\n')
        print(f"[INFO] Finished writing CRAVAT input file!")
    else:
        print(f"[WARNING] CRAVAT input file already exists!")

sys.setrecursionlimit(100)
print('[INFO] Recursion limit set to 100')
def get_cravat_output(session, job_id, output_path):
    '''
    Writes a cravat input file for a given sample name
    
        inputs:
        - session: requests.session where the cravat query is in
        - job_id: can be feteched by post.json()['id']
        - output_path: path to write the CRAVAT output file

        outputs:

    '''
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
        get_cravat_output(session, job_id, output_path)
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

# def create_cravat_assay(cravat_df) -> Assay:
#     """
#     Create CRAVAT assay from CRAVAT dataframe

#     The first column should be SNV's in condensed format (e.g. "chr12:25398284:C/A")
#     The first 2 rows should be CRAVAT output annotation column names
#         1st row name: 'annotator'
#         2nd row name: 'info'

#     Inputs:
#     - cravat_df: a pandas dataframe with the CRAVAT output

#     Outputs:
#     - assay: a missionbio.h5.data.Assay object
#     """

#     # sanity check cravat_df:
#     condensed_format_regex = r'chr(\d+|X|Y):\d+:(A|T|G|C)+/(A|T|G|C)+'
#     if not cravat_df.index.str.contains(condensed_format_regex, regex=True).all():
#         print('[WARNING] The first column should be SNV\'s in condensed format (e.g. "chr12:25398284:C/A")')
    
#     if not cravat.df.columns.names == ['annotator', 'info']:
#         print('[WARNING] The first 2 rows should be CRAVAT output annotation column names')
#         print('    1st row name: \'annotator\'')
#         print('    2nd row name: \'info\'')