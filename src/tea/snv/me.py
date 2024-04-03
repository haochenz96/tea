# import necessary libraries
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import betabinom, binom
from .poibin import PoiBin
import time

def __get_presence_and_null_prob(df_total, df_alt, ado_precision = 15, fp = 0.001):
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
    # this part needs work -- how do you go from read counts to probabilities of presence of mutations?
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

    # if depth is zero, then the probability of presence is zero
    normalized_presence_prob[depth_mat == 0] = 0

    null_prob_vector = pd.Series(np.exp(logsumexp(np.log(normalized_presence_prob), axis=0) - np.log(ncells)), index=df_total.columns)

    return normalized_presence_prob, null_prob_vector

def __compute_mutual_correlation_pval(npp_mut_1, npp_mut_2, null_prob_mut_1, null_prob_mut_2, gametes = [(0,1), (1,0)]):
    """
    Given normalized presence probabilities and null probabilities for one mutation and N other mutations, compute the mutual exclusivity p-value for the pair(s).

    Note that npp_mut_1 must be 1D while npp_mut_2 can be 1D (1 mutation) or 2D (multiple mutations). 

    """
    if npp_mut_1.ndim != 1:
        raise ValueError("npp_mut_2 must be 1D")

    ncells = npp_mut_1.shape[0]

    # get individual gametes:
    g00 = gametes[0][0]
    g01 = gametes[0][1]
    g10 = gametes[1][0]
    g11 = gametes[1][1]

    # null_prob = ( (-1)**(g00+1) * null_prob_mut_1 + 1 - g00)  * ( (-1)**(g01+1) * null_prob_mut_2 + 1 - g01) + ( (-1)**(g10+1) * null_prob_mut_1 + 1 - g10)  * ( (-1)**(g11+1) * null_prob_mut_2 + 1 - g11)

    null_prob = g00 * (1-null_prob_mut_1) * (1-null_prob_mut_2) + g01 * (1-null_prob_mut_1) * (null_prob_mut_2) + g10 * (null_prob_mut_1) * (1 - null_prob_mut_2) + g11 * (null_prob_mut_1) * (null_prob_mut_2)
    if npp_mut_2.ndim == 2:
        # p_vector = (
        #     npp_mut_1[:, np.newaxis] * (1 - npp_mut_2) + 
        #     npp_mut_2 * (1 - npp_mut_1[:, np.newaxis])
        #     )
        
        # p_vector = ( (-1)**(g00+1) * npp_mut_1[:, np.newaxis] + 1 - g00)  * ( (-1)**(g01+1) * npp_mut_2 + 1 - g01) + ( (-1)**(g10+1) * npp_mut_1[:, np.newaxis] + 1 - g10)  * ( (-1)**(g11+1) * npp_mut_2 + 1 - g11)

        p_vector = g00 * (1-npp_mut_1[:, np.newaxis]) * (1-npp_mut_2) + g01 * (1-npp_mut_1[:, np.newaxis]) * (npp_mut_2) + g10 * (npp_mut_1[:, np.newaxis]) * (1 - npp_mut_2) + g11 * (npp_mut_1[:, np.newaxis]) * (npp_mut_2)
    
        # get PoiBin per column
        pb = [PoiBin(p_vector[:, i]) for i in range(p_vector.shape[1])]
    elif npp_mut_2.ndim == 1:
        p_vector = g00 * (1-npp_mut_1) * (1-npp_mut_2) + g01 * (1-npp_mut_1) * (npp_mut_2) + g10 * (npp_mut_1) * (1 - npp_mut_2) + g11 * (npp_mut_1) * (npp_mut_2)
        pb = [PoiBin(p_vector)]
    else:
        raise ValueError("npp_mut_2 must be 1D or 2D")

    ncell_array = np.arange(ncells)

    # me_inv_pval = 0
    # for k in range(ncells):
    #     me_inv_pval += np.array([pb_i.pval(k) for pb_i in pb]) * binom.pmf(k, ncells, null_prob)
    # parallelize above    
    me_inv_pval = (np.array([pb_i.pval(ncell_array) for pb_i in pb]) * np.array([binom.pmf(i, ncells, null_prob) for i in ncell_array]).T).sum(axis=1)

    # @HZ TODO -- sometimes when me_inv_pval is too small, it becomes negative. In that case, set it to 0.
    me_inv_pval[me_inv_pval < 0] = 0
        
    return me_inv_pval

def get_exclusivity(
    depth_df: pd.DataFrame,
    alt_df: pd.DataFrame, 
    ado_precision = 15, # precision parameter
    fp = 0.001, # false positive rate
    gametes = [(0,1), (1,0)],
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
    normalized_presence_prob, null_prob_vector = __get_presence_and_null_prob(depth_df, alt_df, ado_precision, fp)
    
    # # @HZ DEBUG
    # normalized_presence_prob[normalized_presence_prob >= 0.9] = 0.9
    # normalized_presence_prob[normalized_presence_prob < 0.9] = 0.1
    # ####################
    
    print("Time to compute presence probabilities: {}".format(time.time() - start))

    nmutations = depth_df.shape[1]
    exclusivity_mat = pd.DataFrame(np.zeros((nmutations, nmutations)), columns = depth_df.columns, index = depth_df.columns)

    # row-wise iteration
    exclusivity_mat = exclusivity_mat.apply(
        lambda x: __compute_mutual_correlation_pval(
            normalized_presence_prob[x.name].values,
            normalized_presence_prob.values,
            null_prob_vector[x.name],
            null_prob_vector,
            gametes,
        ), axis=0
    )
    print("Time to compute mutual exclusivity: {}".format(time.time() - start))

    # # element-wise iteration
    # exclusivity_mat = exclusivity_mat.apply(
    #     lambda x: exclusivity_mat.apply(
    #         lambda y: __compute_mutual_correlation_pval(
    #             normalized_presence_prob[x.name].values,
    #             normalized_presence_prob[y.name].values,
    #             null_prob_vector[x.name],
    #             null_prob_vector[y.name]
    #         ), axis=0
    #     ), axis=1
    # )

    return exclusivity_mat