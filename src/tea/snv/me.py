# import necessary libraries
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import betabinom, binom
from .poibin import PoiBin
import time
import logging

logging.basicConfig(level=logging.INFO)
logging.getLogger().setLevel(logging.INFO)

# %% [NEW]ME functions
def __get_presence_and_null_prob(
    df_total, df_alt, 
    ado_precision = 15, fp = 0.001,
    rm_irrelevant_cells = True,
    ):
    """
    Given a matrix of read depths and alt read counts, compute (1) presence probabilities and (2) null probabilities for each mutation.

    Inputs
    ----------
    depth_mat: named DataFrame ncells x nmutations. Should already be subset to mutations of interest
    alt_mat: named DataFrame ncells x nmutations. Should already be subset to mutations of interest

	Paramters
	----------
	ado_precision: precision parameter
	fp: false positive rate
	rm_irrelevant_cells: remove cells with zero alt count for all mutations

    Returns
    -------
	normalized_presence_prob: ncells x nmutations
	null_prob_vector: nmutations

    """
    # assert (df_total.columns == df_alt.columns).all()
    # assert (df_total.index == df_alt.index).all()
    # if rm_irrelevant_cells:
    #     # remove cells with zero total depth/alt count for all mutations
    #     coi = df_total.index[(df_alt > 0).any(axis=1)]
    #     logging.info(f"Subsetting to {len(coi)} / {len(df_total)} cells with positive read count for any one mutation")
    #     df_total = df_total.loc[coi]
    #     df_alt = df_alt.loc[coi]
    depth_mat = df_total.values
    alt_mat = df_alt.values
    ncells = depth_mat.shape[0]

    bb_alpha = fp * ado_precision
    bb_beta = (1 - fp) * ado_precision

    presence_coeff_mat = betabinom.logpmf(alt_mat, depth_mat, 1, 1)
    absence_coeff_mat = betabinom.logpmf(alt_mat, depth_mat, bb_alpha, bb_beta)

    presence_absence_tensor = np.stack([presence_coeff_mat, absence_coeff_mat], axis=2)
    total_coeff_mat = logsumexp(presence_absence_tensor, axis=2)

    normalized_presence_prob = pd.DataFrame(np.exp(presence_coeff_mat - total_coeff_mat), index = df_total.index, columns = df_total.columns)

    if rm_irrelevant_cells:
        # if alt_read_count is zero, set the probability to 0, so that it is not used in the ME analysis
        normalized_presence_prob[alt_mat == 0] = 0
    else:
        normalized_presence_prob[depth_mat == 0] = 0
    # @HZ move this to later, because cells of interest depend on which mutation pair is being analyzed
    # null_prob_vector = pd.Series(np.exp(logsumexp(np.log(normalized_presence_prob), axis=0) - np.log(ncells)), index=df_total.columns)

    return normalized_presence_prob # , null_prob_vector

def __compute_mutual_correlation_pval(
    npp_mut_1, npp_mut_2, 
    # null_prob_mut_1, null_prob_mut_2, 
    gametes = [(0,1), (1,0)],
    rm_irrelevant_cells = True
    ):
    """
    Given normalized presence probabilities (npp) for one mutation and N other mutations, compute the mutual exclusivity p-value for the pair(s).
    Note that npp_mut_1 must be 1D while npp_mut_2 can be 1D (1 mutation) or 2D (multiple mutations). 
    
    Inputs
    ----------
    npp_mut_1: 1D array of normalized presence probabilities for mutation 1
    npp_mut_2: 1D or 2D array of normalized presence probabilities for other mutation(s) to compare with mutation 1. If 2D, each column corresponds to a mutation.
    
    Paramters
    ----------
    gametes: list of tuples, each tuple corresponds to a gamete (g00, g01, g10, g11). For example, g01 means mutation 1 is absent in the cell and mutation 2 is present.
    rm_irrelevant_cells: remove cells with zero npp for both mutations. If enabled, npp_mut_2 must be 1D.
    
    Returns
    -------
    me_inv_pval: inverted colocalization p-value for each mutation pair (greater the value, more likely that the pair's colocalization pattern DID NOT happen by chance) If npp_mut_2 is 1D, then a float number; if npp_mut_2 is 2D, then 1D array. 
    
    """
    if npp_mut_1.ndim != 1:
        raise ValueError("npp_mut_2 must be 1D")
    if rm_irrelevant_cells:
        if npp_mut_2.ndim != 1:
            raise ValueError("npp_mut_2 must be 1D if rm_irrelevant_cells is enabled")
        else:
            # select cells with positive read count for either mutation
            coi = npp_mut_1.index[(npp_mut_1 > 0) | (npp_mut_2 > 0)]
            npp_mut_1 = npp_mut_1.loc[coi].values
            npp_mut_2 = npp_mut_2.loc[coi].values
            ncells = len(coi)
    else:
        npp_mut_1 = npp_mut_1.values
        npp_mut_2 = npp_mut_2.values
        ncells = npp_mut_1.shape[0]
    null_prob_mut_1 = np.exp(logsumexp(np.log(npp_mut_1), axis=0) - np.log(ncells))
    null_prob_mut_2 = np.exp(logsumexp(np.log(npp_mut_2), axis=0) - np.log(ncells))

    # get individual gametes:
    g00 = gametes[0][0]
    g01 = gametes[0][1]
    g10 = gametes[1][0]
    g11 = gametes[1][1]

    # null_prob = ( (-1)**(g00+1) * null_prob_mut_1 + 1 - g00)  * ( (-1)**(g01+1) * null_prob_mut_2 + 1 - g01) + ( (-1)**(g10+1) * null_prob_mut_1 + 1 - g10)  * ( (-1)**(g11+1) * null_prob_mut_2 + 1 - g11)

    null_prob = g00 * (1-null_prob_mut_1) * (1-null_prob_mut_2) + g01 * (1-null_prob_mut_1) * (null_prob_mut_2) + g10 * (null_prob_mut_1) * (1 - null_prob_mut_2) + g11 * (null_prob_mut_1) * (null_prob_mut_2)
    if npp_mut_2.ndim == 2:
        p_vector = g00 * (1-npp_mut_1[:, np.newaxis]) * (1-npp_mut_2) + g01 * (1-npp_mut_1[:, np.newaxis]) * (npp_mut_2) + g10 * (npp_mut_1[:, np.newaxis]) * (1 - npp_mut_2) + g11 * (npp_mut_1[:, np.newaxis]) * (npp_mut_2)
    
        # get PoiBin per column
        pb = [PoiBin(p_vector[:, i]) for i in range(p_vector.shape[1])]
    elif npp_mut_2.ndim == 1:
        p_vector = g00 * (1-npp_mut_1) * (1-npp_mut_2) + g01 * (1-npp_mut_1) * (npp_mut_2) + g10 * (npp_mut_1) * (1 - npp_mut_2) + g11 * (npp_mut_1) * (npp_mut_2)
        pb = [PoiBin(p_vector)]
    else:
        raise ValueError("npp_mut_2 must be 1D or 2D")

    ncell_array = np.arange(ncells)
    
    me_inv_pval = (
        np.array( [pb_i.pval(ncell_array) for pb_i in pb] ) * np.array( [binom.pmf(i, ncells, null_prob) for i in ncell_array] ).T
        ).sum(axis=1)

    # me_inv_pval = 0
    # for k in range(ncells):
    #     me_inv_pval += np.array([pb_i.pval(k) for pb_i in pb]) * binom.pmf(k, ncells, null_prob)
        
    return me_inv_pval

def get_exclusivity(
    depth_df: pd.DataFrame,
    alt_df: pd.DataFrame, 
    ado_precision = 15, # precision parameter
    fp = 0.001, # false positive rate
    gametes = [(0,1), (1,0)],
    rm_irrelevant_cells = True,
    ):
    """
    Given a matrix of read depths and alt read counts, compute the exclusivity p-value for each pair of mutations.

    Parameters
    ----------

    depth_mat: ncells x nmutations
    alt_mat: ncells x nmutations

    Returns
    -------
    exclusivity_mat: nmutations x nmutations

    """

    # time each step
    start = time.time()
    normalized_presence_prob = __get_presence_and_null_prob(depth_df, alt_df, ado_precision, fp, rm_irrelevant_cells)
    
    # # @HZ DEBUG
    # normalized_presence_prob[normalized_presence_prob >= 0.9] = 0.9
    # normalized_presence_prob[normalized_presence_prob < 0.9] = 0.1
    # ####################
    
    print("Time to compute presence probabilities: {}".format(time.time() - start))

    nmutations = depth_df.shape[1]
    exclusivity_mat = pd.DataFrame(np.zeros((nmutations, nmutations)), columns = depth_df.columns, index = depth_df.columns)

    if rm_irrelevant_cells:
        # element-wise iteration
        for i in range(nmutations):
            for j in range(i+1, nmutations):
                exclusivity_mat.iloc[i,j] = __compute_mutual_correlation_pval(
                    normalized_presence_prob.iloc[:,i], 
                    normalized_presence_prob.iloc[:,j], 
                    gametes,
                    rm_irrelevant_cells,
                )
        # mirrow the upper triangle
        exclusivity_mat = exclusivity_mat + exclusivity_mat.T

    else:# row-wise iteration
        exclusivity_mat = exclusivity_mat.apply(
            lambda x: __compute_mutual_correlation_pval(
                normalized_presence_prob[x.name],
                normalized_presence_prob,
                gametes,
                rm_irrelevant_cells,
            ), axis=0
        )
    print("Time to compute mutual exclusivity: {}".format(time.time() - start))

    return exclusivity_mat
