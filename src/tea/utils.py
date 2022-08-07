#@HZ
import numpy as np
import pandas as pd
from pysam import VariantFile
import logging
logger = logging.getLogger(__name__)
from datetime import datetime
import pytz as tz

# ----- get current date and time
def get_simple_timestamp(timezone = 'US/Eastern') -> str:
    '''
    function to get current date and time in a particular timezone
    '''
    tz_obj = tz.timezone(timezone) # <------ uses US eastern time by default
    now = datetime.now(tz.utc).astimezone(tz_obj) 
    timestamp = now.strftime("%a %m %d %H:%M:%S %Z %Y")
    # datetime_simple = now.strftime("%Y-%m-%d--%H_%M")
    return f'[{timezone}]] {timestamp}'


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
