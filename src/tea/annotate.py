from pybedtools import BedTool
import pandas as pd
from pathlib import Path
import mosaic.io as mio
import numpy as np

def annotate_snv_amplicon_coverage(sample, working_dir: str or pathlib.Path, insert_bed: str or pathlib.Path = None, ) -> pd.DataFrame:
    '''
    Annotate SNV dataframe with amplicon coverage information

    Inputs:
        sample:
            MissionBio sample object.
        insert_bed: 
            amplicon insert bed file path
        working_dir:
            working directory path

    Outputs:
        snv_panel_bed_isec_df: 
            annotated SNV dataframe
    '''
    working_dir = Path(working_dir)
    
    # infer amplicon bed from sample if needed
    if insert_bed is None:
        panel_version = sample.cnv.metadata['panel_version']
        insert_bed_path = working_dir / f'{panel_version}-amplicon.bed'
        if not insert_bed_path.exists():
            raise ValueError(f'amplicon bed file not found: {insert_bed_path}')
            exit(1)
        insert_bed_df = pd.DataFrame()
        insert_bed_df['CHROM'] = sample.cnv.col_attrs['CHROM']
        if not insert_bed_df['CHROM'].str.startswith('chr').all():
            insert_bed_df['CHROM'] = 'chr' + insert_bed_df['CHROM']
        insert_bed_df['start_pos'] = sample.cnv.col_attrs['start_pos']
        insert_bed_df['end_pos'] = sample.cnv.col_attrs['end_pos']
        insert_bed_df['amplicon_id'] = sample.cnv.col_attrs['id']
        insert_bed_df.to_csv(insert_bed_path, sep='\t', index=False)
        insert_bed = BedTool(insert_bed_path)
    else:
        insert_bed = BedTool(insert_bed)

    # parse condensed foramt
    snv_bed_df = pd.DataFrame(index = sample.dna.ids())
    snv_bed_df['chrom'] = pd.Series(snv_bed_df.index).str.split(':', expand=True)[0].values
    snv_bed_df['start'] = pd.Series(snv_bed_df.index).str.split(':', expand=True)[1].values
    snv_bed_df['end'] = snv_bed_df['start'] 
    snv_bed_df['condensed_format'] = snv_bed_df.index.values
    snv_bed_df[['chrom', 'start', 'end', 'condensed_format']].to_csv(working_dir / 'snv_bed.bed', sep='\t', index=False, header=False) # save it to a bed file; order of columns matters

    snv_bed = BedTool(working_dir / 'snv_bed.bed')
    snv_panel_bed_isec = snv_bed.intersect(insert_bed, wao = True)
    snv_panel_bed_isec.saveas(working_dir / 'snv_bed_isec_with_panel_insert.bed')

    snv_panel_bed_isec_df = pd.read_csv((working_dir / 'snv_bed_isec_with_panel_insert.bed'), sep='\t', usecols = [3,4,5,6,7], names = ['condensed_format','chrom', 'start', 'end', 'amplicon_id']).set_index('condensed_format')

    # check to make sure the isec_df is the same shape as the snv_df
    assert snv_panel_bed_isec_df.shape[0] == sample.dna.shape[1], '[ERROR] The number of SNVs in the intersected DF does not match the number of SNVs in the original Tapestri data. Check if file formats are incorrect!'

    # fetch amplicon mean rc
    amp_median_rc = sample.cnv.get_attribute('read_counts', constraint='row').median(axis=0)
    snv_panel_bed_isec_df['median_rc'] = snv_panel_bed_isec_df['amplicon_id'].apply(lambda x: amp_median_rc[x] if x != '.' else np.nan)

    return snv_panel_bed_isec_df
