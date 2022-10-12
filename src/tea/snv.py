import mosaic.io as mio
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px

# snv_pair_correlation
import scipy
import scipy.cluster.hierarchy as sch
from scipy.stats import betabinom
import matplotlib.pyplot as plt
from scipy.stats import binom
from scipy.stats import beta
from scipy.special import rel_entr
from scipy.stats import entropy

# __reduce_concat
from functools import reduce

from distutils.command.sdist import sdist

import logging
logger = logging.getLogger(__name__)

def mut_likelihood(f, a, b):
    '''
    f: float
        VAF value to evaluate
    a: np.array or pd.DataFrame
        number of alt reads
    b: np.array or pd.DataFrame
        number of ref reads (DP - alt)
    '''
    out_mat = 1 - beta.cdf(f, a, b) # inverse beta CDF

    # if ALT == 0, assign 0
    out_mat = np.where(a == 0, 0, out_mat)
    
    # if ALT == 0 and DP == 0, assign nan
    out_mat = np.where((a == 0) & (b == 0), np.nan, out_mat) 

    # if ALT == DP != 0, assign 1
    out_mat = np.where((a > 0) & (b == 0), 1, out_mat)
    
    return out_mat

def aggregate_likelihoods(p, q, verbose=False, likelihood_threshold = 0.5):
    # For cells where the variant has DP = 0, its likelihood is 0
    p = np.nan_to_num(p)
    q = np.nan_to_num(q)
    
    diff = abs(p - q)
    
    p_indices = p > likelihood_threshold
    num_p_indices = np.sum(p_indices)
    q_indices = q > likelihood_threshold
    num_q_indices = np.sum(q_indices)
    common_indices = p_indices & q_indices

    if num_p_indices == 0 or num_q_indices == 0:
        # if either SNV has any cell that has a likelihood > 0.5, return 0
        return 0

    forward_diff = diff[p_indices].sum() / num_p_indices

    reverse_diff = diff[q_indices].sum() / num_q_indices

    if verbose == True:
        logging.info(f'number of p indices: {num_p_indices}')
        logging.info(f'number of q indices: {num_q_indices}')
        logging.info(f'number of common indices: {common_indices.sum()}')
        logging.info(f'forward diff: {forward_diff}')
        logging.info(f'reverse diff: {reverse_diff}')

    # normed_p = p[common_indices] / p[common_indices].sum()
    # normed_q = q[common_indices] / q[common_indices].sum()

    solution = min(
        forward_diff,
        reverse_diff,
    )
    # if solution == np.inf:
    #     # solution is infinity, which means that there is no common index where both SNVs have positive DP. Usually the case for ultrarare SNVs. The mutual correlation score is assigned as 0.
    #     return 0
    # else:
    return 1 - solution


    # num = ( ~np.isnan(entropy_vector) ).sum()
    # return np.nansum( entropy_vector ) / num

    # return -np.nansum( entropy_vector ) / -np.nansum( p * np.log(p) )

def snv_pair_mutual_correlation(alt_mat, dp_mat, voi=None, f = 0.2):

    '''
    f controls how "sensitive" the test is to the difference between the two SNVs. The higher the more distinction spotted.
    '''

    if voi is None:
        # sanity check that both matricies have the same column names:
        assert (alt_mat.columns == dp_mat.columns).all()
    else:
        try:
            alt_mat = alt_mat.loc[:, voi]
            dp_mat = dp_mat.loc[:, voi]
        except KeyError:
            raise KeyError('some of voi not found in alt_mat or dp_mat')
    
    # for each snv in each single cell, calculate inverse beta CDF
    raw_inv_beta_cdf_array = pd.DataFrame(
        data = mut_likelihood(f, alt_mat, dp_mat-alt_mat),
        columns = alt_mat.columns,
    )

    out_df = pd.DataFrame(data = raw_inv_beta_cdf_array, index = alt_mat.columns, columns=alt_mat.columns)
    for i in out_df.index:
        for j in out_df.columns:
            out_df.loc[i,j] = aggregate_likelihoods(
                raw_inv_beta_cdf_array.loc[:, i], 
                raw_inv_beta_cdf_array.loc[:, j]
                )
    return out_df

def plot_snv_correlation_heatmap(
    me_mat, 
    color_scale = ["blue", "yellow", "red"],
    cluster = True,
    ):
    if cluster:
        me_mat = cluster_corr(me_mat)
    fig = px.imshow(me_mat, color_continuous_scale=color_scale, zmin=0, zmax=1)
    fig.update_layout(height = 800, width = 800)
    
    return fig

def translate(trinuc_seq):
       
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    aa_seq =""
    if len(trinuc_seq)%3 == 0:
        for i in range(0, len(trinuc_seq), 3):
            codon = trinuc_seq[i:i + 3]
            aa_seq+= table[codon]
    return aa_seq

def __reduce_concat(x, sep=""):
    return reduce(lambda x, y: str(x) + sep + str(y), x)

def paste(*lists, sep=" ", collapse=None):
    result = map(lambda x: __reduce_concat(x, sep=sep), zip(*lists))
    if collapse is not None:
        return __reduce_concat(result, sep=collapse)
    return list(result)

def __mutate(codon, pos, nt):
    '''
    Mutate a codon at a given position to a given nucleotide

    Parameters
    ----------
    codon : str
        The codon to mutate
    pos : int
        The position to mutate (1, 2, or 3)
    nt : str
        The nucleotide to mutate to

    Returns
    -------
    str
    '''
    if not len(codon) == 3:
        raise ValueError
    if not pos in [0, 1, 2]:
        raise ValueError
    if not nt in ['A', 'C', 'G', 'T']:
        raise ValueError
    mutated_codon = codon[:pos] + nt + codon[pos+1:]
    return mutated_codon

nts = np.array(['A', 'C', 'G', 'T'])
trinuc_list = paste(
    np.repeat(nts, 16),
    np.repeat(np.tile(nts, 4), 4),
    np.tile(nts, 16),
    sep = '',
)

# single-nt mutation in all possible trinucleotide contexts
IMPACT_MATRIX = pd.DataFrame(np.zeros((64,64)))
IMPACT_MATRIX.columns = IMPACT_MATRIX.index = trinuc_list
for trinuc_i in IMPACT_MATRIX.index:
    for trinuc_j in IMPACT_MATRIX.columns:
        from_aa = translate(trinuc_i)
        to_aa = translate(trinuc_j)
        if from_aa == to_aa:
            IMPACT_MATRIX.loc[trinuc_i, trinuc_j] = 1
        elif to_aa == '*':
            IMPACT_MATRIX.loc[trinuc_i, trinuc_j] = 3
        elif to_aa != '*' and from_aa != '*' and to_aa != from_aa:
            IMPACT_MATRIX.loc[trinuc_i, trinuc_j] = 2
        else: # from _aa == '*'
            IMPACT_MATRIX.loc[trinuc_i, trinuc_j] = np.nan

trinuc_ind = pd.Series(range(64), index = trinuc_list) # all the trinucleotide mutation patterns. Should be 4 * 3 * (4^2) = 192 
TRINUC_SUBS = np.array(
    paste(
        np.repeat(nts, 48),
        np.tile(np.repeat(nts, 3), 16),
        np.tile(np.repeat(nts, 12), 4),
        np.repeat('>', 192),
        np.repeat(nts, 48),
        np.array([nts[nts != j].tolist() for j in np.tile(nts, 16)]).flatten(),
        np.tile(np.repeat(nts, 12), 4),
        sep = '',
    )
)
trinuc_subsind = pd.Series(range(192), index=TRINUC_SUBS)
# sanity check
# np.unique(np.array(trinuc_subs)).shape

def get_all_possible_mutation_mat(seq, orf_start, overlap_start, orf_end, overlap_end, strand = '+', log_level = logging.INFO):
    """
    seq: a string of nucleotides
    (all the positions need to be 0 based)
    orf_start: the position of the coding exon sequence
    overlap_start: the position of the overlap sequence
    orf_end: the position of the coding exon sequence
    overlap_end: the position of the overlap sequence
    strand: the strand of the sequence

    """
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    # by default, set orf_end to the second to last position
    if orf_start < 1: 
        raise ValueError('ORF start position must be >= 1.')
    if overlap_start < orf_start:
        raise ValueError('Overlap start must be >= ORF start.')
    if orf_end <= overlap_start:
        raise ValueError('ORF end must be > overlap start.')
    if overlap_end > orf_end:
        raise ValueError('Overlap end must be <= ORF end.')
    if orf_end > len(seq):
        raise ValueError('End position must be less than the length of the sequence.')
    # if (orf_end-orf_start+1) % 3 != 0 and 3 - ( (orf_end-orf_start+1) % 3 ) > (len(seq) - orf_end): # if the last codon cannot be completed
    #     raise ValueError('last codon cannot be completed with current setting. Please adjust the sequence and orf_start/orf_end positions.')
    
    seq = np.array([i for i in seq.upper()])
    if strand == '-':
        # reverse strand
        seq = seq[::-1]
    seq_1up = np.array( seq[overlap_start-1:overlap_end] ) # 1-bp upstream e.g. ATG
    seq_mid = np.array( seq[overlap_start:overlap_end+1] ) # middle e.g. TGC
    seq_1down = np.array( seq[overlap_start+1:overlap_end+2] ) # 1-bp downstream e.g. GCA
    seq_length = len(seq_mid)

    ##### ----- trinucleotide context ----- #####
    indices = np.repeat([i for i in range(seq_length)], 3) # times 3 because for each single nucleotide, there are 3 ways to mutate it; e.g. array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5])
    old_trinuc = np.array(
        paste(seq_1up[indices], seq_mid[indices], seq_1down[indices], sep = "") 
        ) # array(['CAT', 'CAT', 'CAT', 'ATG', 'ATG', 'ATG', 'TGT', 'TGT', 'TGT','GTC', 'GTC', 'GTC', 'TCC', 'TCC', 'TCC', 'CCA', 'CCA', 'CCA'])
    new_bases = np.array(
        [nts[nts != x].tolist() for x in seq_mid] 
        ).flatten() # generate all possible new bases at each position e.g. array(['C', 'G', 'T', 'A', 'C', 'G', 'A', 'C', 'T', 'A', 'C', 'G', 'A', 'G', 'T', 'A', 'G', 'T'])
    new_trinuc = np.array(
        paste(seq_1up[indices], new_bases, seq_1down[indices], sep="") 
        ) # array(['CCT', 'CGT', 'CTT', 'AAG', 'ACG', 'AGG', 'TAT', 'TCT', 'TTT', 'GAC', 'GCC', 'GGC', 'TAC', 'TGC', 'TTC', 'CAA', 'CGA', 'CTA'])

    ##### ----- codon ----- #####
    start_skip = (3-overlap_start+orf_start) % 3 # number of bps to skip at the beginning
    end_skip = (3-orf_end+overlap_end) % 3 # number of bps to skip at the end
    num_codon = (seq_length - start_skip - end_skip) / 3
    assert num_codon.is_integer(), f'the input ORF is not a multiple of 3; seq_length {seq_length}, start_skip {start_skip}, end_skip {end_skip}, num_codon {num_codon}'
    num_codon = int(num_codon)

    codon_start_ind = np.repeat([i for i in range(start_skip, seq_length-end_skip, 3)], 9) # relative indices for codon start in the sequence; multiplies by 9 here because for each codon, there are 9 possible ways to mutate it
    # e.g. array([0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3])
    pos_in_codon_to_mutate = np.tile(np.repeat([i for i in range(3)], 3), num_codon) # e.g. array([0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2])

    if start_skip > 0 or end_skip > 0:
        logger.debug('the first/last codon would be partially considered for mutation.')

        codon_start_ind = np.concatenate([
            [start_skip - 3] * (start_skip) * 3, 
            codon_start_ind, 
            [start_skip + num_codon*3] * (end_skip) * 3
            ]).astype(int)
        pos_in_codon_to_mutate = np.concatenate([
            np.repeat([i for i in range(start_skip)], 3),
            pos_in_codon_to_mutate, 
            np.repeat([i for i in range(end_skip)], 3) 
            ]).astype(int) # the last codon will only have 3 or 6 possible mutations

    old_codon = np.array(
        paste(seq[overlap_start + codon_start_ind], seq[overlap_start + codon_start_ind + 1], seq[overlap_start + codon_start_ind + 2], sep="")
        ) # array(['ATG', 'ATG', 'ATG', 'ATG', 'ATG', 'ATG', 'ATG', 'ATG', 'ATG', 'TCC', 'TCC', 'TCC', 'TCC', 'TCC', 'TCC', 'TCC', 'TCC', 'TCC'])
    
    new_codon = np.array(
        [__mutate(
            old_codon[x], 
            pos_in_codon_to_mutate[x], 
            new_bases[x]
            ) for x in range(len(pos_in_codon_to_mutate))]) # e.g. array(['CTG', 'GTG', 'TTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT', 'ACC', 'CCC', 'GCC', 'TAC', 'TGC', 'TTC', 'TCA', 'TCG', 'TCT']

    imp = np.array([ 
        IMPACT_MATRIX.loc[old_codon[x], new_codon[x]] 
        for x in range(len(old_codon)) 
        ]) # array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1])

    ##### ----- create L matrix ----- #####
    trinuc_subs_oi = pd.DataFrame(
        data = imp,
        index = paste(old_trinuc, new_trinuc, sep=">"),
        columns = ['IMPACT']
        )

    L = pd.DataFrame(index = TRINUC_SUBS, columns = ['Synonymous', 'Missense', 'Nonsense'])
    # Synonymous
    val_counts = trinuc_subs_oi[imp==1].index.value_counts()
    L.loc[ val_counts.index, 'Synonymous' ] = val_counts.values

    # Missense
    val_counts = trinuc_subs_oi[imp==2].index.value_counts()
    L.loc[ val_counts.index, 'Missense' ] = val_counts.values

    # Nonsense
    val_counts = trinuc_subs_oi[imp==3].index.value_counts()
    L.loc[ val_counts.index, 'Nonsense' ] = val_counts.values
    return L

def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold, 
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    
    if not inplace:
        corr_array = corr_array.copy()
    
    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]
    return corr_array[idx, :][:, idx]