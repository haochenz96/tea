# __reduce_concat
from functools import reduce

##############################################
#                  3. dnds                   #
##############################################

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