#@HZ
import numpy as np
import pandas as pd
from pysam import VariantFile

def get_ann(samples): # Obsolete
    '''
    create mapping between genomic variant and HGVSp short
    input: 
    - samples: can be a single missionbio.mosaic.sample.Sample instance, 
               or a list of that instance
    return:
    - ann_map: maps genomic info to protein info
               e.g. 'chr12:25398284:C/T': '(PATH) KRAS:p.G12D '
    '''
    if isinstance(samples, dict):
        samples = list(samples.values())
        
    if isinstance(samples, list):
        map_list = []
        for i in range(len(samples)):
            annotation_map = {}
            for j, k in zip(samples[i].dna.ids(), samples[i].dna.get_annotations()):
                annotation_map[j] = k.split('-')[0]
                map_list.append(annotation_map)
        ann_map = {}
        for d in map_list:
            ann_map.update(d)
        return ann_map
    
    elif isinstance(samples, mio.Sample):
        annotation_map = {}
        for j, k in zip(samples.dna.ids(), samples.dna.get_annotations()):
            annotation_map[j] = k.split('-')[0]
            return annotation_map
    else:
        print("input error")


def sort_for_var(dna, vars, attribute, method='hier'):

    '''
    function to sort a particular set of barcodes for DNA values
    
    input:
    - dna: a single missionbio.mosaic.sample.Sample.dna instance
    - vars: a list of variants of the form 'chr12:25398284:C/T'
    - attribute: "AF_MISSING", "NGT" etc.
    - method: 'hier' (multiple vars) or 'single_var' (single var, ascending value)
    
    return:
      ordered cell barcodes (hierarhchically)
    '''

    if not isinstance(vars, list):
        vars = [vars]

    # check if every variant is present in the provided matrix
    good_vars = []
    for i in vars:
        if not i in dna.ids():
            print(f'[WARNING] --- {i} not in dna.ids')
        else:
            good_vars.append(i)

    if method == 'hier':
        assay = dna[dna.barcodes(), good_vars]
        return assay.clustered_barcodes(orderby = attribute)
        
    elif method == 'single_var':
        data = dna.get_attribute(attribute, constraint = 'row')
        return data.sort_values(by=good_vars).index

def label_for_var(sample, vars_of_interest, AF_threshold=20, min_mut_var=0):
    '''
    function to create label based on the VAF of a particular set of variants, with the option to filter based on genomic variant quality (GQ). 
    '''

    if not isinstance(vars_of_interest, list):
        vars_of_interest = [vars_of_interest]
    
    # get mutant cells: 
    AF_df = sample.dna.get_attribute('AF_MISSING', constraint='row+col')
    mut_cells = AF_df[vars_of_interest].index[((AF_df[vars_of_interest] > AF_threshold).sum(axis=1) > min_mut_var) == True].tolist()
    
    sample.dna.set_labels(clusters = {'mutant_cells': mut_cells},
                          others_label='others')

def rand_3split_normal_cells(sample):

    '''
    randomly splits the normal cells into 3 parts

    ---
    input:

    - sample: missionbio.mosaic.sample object, which should already have a proportion of cells labeled 'normal'

    '''

    normal_idx = np.argwhere(sample.dna.get_labels() == 'normal')
    #print(normal_idx)
    normal_bars = sample.dna.barcodes()[sample.dna.get_labels() == 'normal']

    ran_idx = np.random.choice(a= normal_idx.shape[0],\
                     size = (int(normal_idx.shape[0]/3), 2), replace = False)

    norm_a = sample.dna.barcodes()[normal_idx[ran_idx[:,0]]]
    norm_b = sample.dna.barcodes()[normal_idx[ran_idx[:,1]]]
    norm_c = normal_bars[~np.in1d(
        normal_bars,
        np.concatenate([norm_a,norm_b])
        )
    ]

    return norm_a.T[0], norm_b.T[0], norm_c

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
    out_df = pd.DataFrame(columns = ['chr', 'start', 'end', 'ref_base', 'alt_base','ref_read_count', 'alt_read_count', 'AF'])

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
                print(f'[WARNING] position - {chr}:{start} has more than 2 alleles')
                # consider multiallelic variant site
                for sub_idx in range(len(alt_bases)):
                    out_df.loc[idx + sub_idx, ['chr', 'start', 'end', 'ref_base', 'alt_base','ref_read_count', 'alt_read_count']] = [chr, start, end, ref_base, alt_bases[sub_idx], ref_read_count, alt_read_counts[sub_idx]]
                    
                    out_df.loc[idx + sub_idx, 'AF'] = alt_read_counts[sub_idx] / (ref_read_count + sum(alt_read_counts))
            
                idx += len(alt_bases)

            else:
                if len(alt_read_counts) != 0:
                    AF = alt_read_counts[0] / (ref_read_count + alt_read_counts[0])
                    if AF == 0:
                        print(f'[WARNING] position - {chr}:{start} AF is 0')
                        # skip if AF is 0
                    else:
                        out_df.loc[idx, ['chr', 'start', 'end', 'ref_base', 'alt_base','ref_read_count', 'alt_read_count']] = [chr, start, end, ref_base, alt_bases[0], ref_read_count, alt_read_counts[0]]
                        out_df.loc[idx, 'AF'] = AF

                # skip if no AD info is available
                else:
                    print(f'[WARNING] position - {chr}:{start} no alternate allele information')
                idx += 1

                
        else: # when no sample name is provided, only keep variant info
            if len(alt_bases) > 1:
                print(f'[WARNING] position - {chr}:{start} has more than 2 alleles')
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

def isNaN(num):
    # check for float('nan') values
    return num != num
