# parse CRAVAT output
import io
from multiprocessing.sharedctypes import Value
import pandas as pd
import numpy as np
from pathlib import Path
from pysam import VariantFile
import logging
logger = logging.getLogger(__name__)

def read_until_match(in_file, match_pattern, num_occurences=1, include_last_match=False):

    '''
    Read a file until a match pattern is found.

    Parameters
    ----------
    in_file : str
        path to input file object
    match_pattern : function that maps a line to a boolean

    num_occurences : int

    include_last_match : bool
        if True, include the last match in the output
    
    Returns
    -------
    out_file : file object
        file object ending at the matched pattern

    '''
    out_file = io.StringIO()
    with open(in_file) as in_f:
        count = 0
        while count < num_occurences:
            line = in_f.readline()
            #print(line)
            if match_pattern(line):
                #print(line)
                count += 1
                if count == num_occurences:
                    if include_last_match:
                        out_file.write(line)
                        break
                    else:
                        break
                else:
                    out_file.write(line)
            else:
                out_file.write(line)
    out_file.seek(0)
    return out_file

def read_between_matches(in_file, match_pattern, num_start, num_end, include_last_match=False):
    out_file = io.StringIO()
    with open(in_file) as in_f:
        match_count = 0
        for line in in_f:
            if match_pattern(line):
                match_count += 1

            if match_count < num_start:
                continue
            elif match_count >= num_start and match_count < num_end:
                out_file.write(line)
            elif match_count == num_end and include_last_match == True:
                out_file.write(line)
                break
            else:
                break
    out_file.seek(0)
    return out_file

def adjust_multiindex(df) -> pd.DataFrame:
    """
    Rename unamed columns name for a multi-index Pandas DataFrame

    See https://stackoverflow.com/questions/41221079/rename-multiindex-columns-in-pandas

    Parameters
    ----------
    df : pd.DataFrame object
        Input dataframe

    Returns
    -------
    pd.DataFrame
        Output dataframe

    """
    new_index = []

    for index_pair_i in df.columns:
        if 'Unnamed' in index_pair_i[0]:
            try: 
                last_seen
            except NameError: 
                last_seen = None
            if last_seen is None:
                raise ValueError('Unnamed column found before named column')
            new_index_pair = (last_seen, index_pair_i[1])
        else:
            last_seen = index_pair_i[0]
            new_index_pair = index_pair_i
        new_index.append(new_index_pair)

    df.columns = pd.MultiIndex.from_tuples(new_index, names=["annotator", "info"])
    return df

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

# get pan-cohort variant protein annotation
def write_cohort_wide_variant_protein_annotation(cravat_dfs, cohort_name, output_dir):
    """
    Given multiple cravat dataframes, get all unique variants' protein annotation and write to a txt file

    Parameters
    ----------
    cravat_dfs : dict
        A dictionary mapping sample name to its cravat df
        key: str
            The sample name
        value: pandas.DataFrame
            cravat_df, formatted by clean_and_format_cravat_df()

    """
    ann_map = {}

    for sample_i in cravat_dfs:
        cravat_df = cravat_dfs[sample_i]

        for var_i, row in cravat_df.iterrows():
            if not var_i in ann_map:
                ann = row['Variant Annotation']['Gene'] + ' ' + row['Variant Annotation']['Protein Change']
                ann_map[var_i] = ann

    ann_map_df = pd.DataFrame.from_dict(ann_map, orient='index')
    ann_map_df.columns = ['variant_annotation']
    ann_map_df.index.rename('condensed_format', inplace=True)

    if not isinstance(output_dir, Path):
        output_dir = Path(working_dir)
    ann_map_df.to_csv( (output_dir / f'{cohort_name}_all_var_annotation.txt'), sep = '\t' )

def highlight_for_col_condition(s, queries):
    '''
    Given a dataframe, highlight the rows that meet the condition in the specified columns.

    Parameters
    ----------
    cravat_df : pandas.DataFrame
        THe dataframe to be highlighted.
    
    queries : dict
        a dictionary mapping each column of interest to condition and style.
        key: str
            The column name.
        value: tuple
            (condition, style)
    '''
    for col_of_interest, action in queries.items():
        idx = pd.IndexSlice
        _slice = idx[:, col_of_interest] # this allows for multiindexxing
        condition = action[0]
        pos_style = action[1]
        s = s.applymap(lambda cell: __style_cell(cell, condition, pos_style), subset=_slice)

    return s

def __style_cell(cell, condition, pos_style = 'background-color:lightgreen;', neg_style = ''):
    """
    Takes a scalar and returns a string with
    the css property `'color: red'` for negative
    strings, black otherwise.
    """
    if condition(cell):
        return pos_style
    else:
        return neg_style

def read_vcf_to_df(in_bcf, sample_name=None):
    '''
    Reads a vcf file and returns a pandas dataframe with the fields specified in fields_to_keep

    inputs:
    - in_bcf: path to vcf file
    - sample_name: sample name to be used as index (only one at a time)
        if None:
            the dataframe will only have the variant info

    outputs:
    - out_df: parsed pandas dataframe 
    '''
    if not isinstance(in_bcf, VariantFile):
        in_bcf = VariantFile(in_bcf)
    out_df = pd.DataFrame(columns = ['chr', 'start', 'end', 'ref_base', 'alt_base'])

    idx = 0
    for rec in in_bcf.fetch():

        chr = rec.chrom
        if not chr.startswith('chr'):
            chr = 'chr' + chr # add chr prefix for b37 chromosome notation
        start = rec.pos
        end = rec.stop

        alleles = rec.alleles
        ref_base=alleles[0]
        alt_bases=alleles[1:]

        if sample_name is not None:
            # get ref and alt counts
            AD_field_idx = np.where([i == 'AD' for i in rec.format.keys()])[0][0]
            AD_field = rec.samples[sample_name].items()[AD_field_idx]
            ref_read_count = AD_field[1][0]
            alt_read_counts = AD_field[1][1:]

            if len(alt_bases) > 1:
                logger.debug(f'[read_vcf_to_df][warning] position - {chr}:{start} has more than 2 alleles')
                # consider multiallelic variant site
                for sub_idx in range(len(alt_bases)):
                    out_df.loc[idx + sub_idx, ['chr', 'start', 'end', 'ref_base', 'alt_base','ref_read_count', 'alt_read_count']] = [chr, start, end, ref_base, alt_bases[sub_idx], ref_read_count, alt_read_counts[sub_idx]]
                    
                    out_df.loc[idx + sub_idx, 'AF'] = alt_read_counts[sub_idx] / (ref_read_count + sum(alt_read_counts))
            
                idx += len(alt_bases)

            else:
                if len(alt_read_counts) != 0:
                    try:
                        AF = alt_read_counts[0] / (ref_read_count + alt_read_counts[0])
                    except ZeroDivisionError:
                        logger.debug(f'[read_vcf_to_df][warning] position - {chr}:{start} has 0 depth')
                        idx += 1
                        continue
                    if AF != 0:
                        out_df.loc[idx, ['chr', 'start', 'end', 'ref_base', 'alt_base','ref_read_count', 'alt_read_count']] = [chr, start, end, ref_base, alt_bases[0], ref_read_count, alt_read_counts[0]]
                        out_df.loc[idx, 'AF'] = AF

                # skip if no AD info is available
                else:
                    logger.debug(f'[read_vcf_to_df][warning] position - {chr}:{start} no alternate allele information')
                idx += 1
          
        else: # when no sample name is provided, only keep variant info
            if len(alt_bases) > 1:
                logger.debug(f'[read_vcf_to_df][warning] position - {chr}:{start} has more than 2 alleles')
                # consider multiallelic variant site
                for sub_idx in range(len(alt_bases)):
                    out_df.loc[idx + sub_idx, ['chr', 'start', 'end', 'ref_base', 'alt_base']] = [chr, start, end, ref_base, alt_bases[sub_idx]]            
                idx += len(alt_bases)

            else:    
                out_df.loc[idx, ['chr', 'start', 'end', 'ref_base', 'alt_base']] = [chr, start, end, ref_base, alt_bases[0]]
                idx += 1

    out_df['condensed_format'] = out_df['chr'].astype(str) + ':' + out_df['start'].astype(str) + ':' + out_df['ref_base'].astype(str) + '/' + out_df['alt_base'].astype(str)
    out_df.set_index('condensed_format', inplace=True)
    return out_df