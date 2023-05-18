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

from distutils.command.sdist import sdist

import logging
logger = logging.getLogger(__name__)

##############################################
#         1. snv_pair_correlation            #
##############################################

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
