import pysam

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt

def config_params(font_size=7):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'

import logging
logger = logging.getLogger(__name__)

##############################################
#            mutational signature            #
##############################################

def _get_snv_trinuc_context(snv, fasta: pysam.FastaFile) -> str:
    """
    Fetch trinucleotide context for a given snv, specified in the form 'chr:pos:ref/alt'

    Parameters
    ----------
    snv : str
        snv in the form 'chr:pos:ref/alt'
    
    fasta : pysam.FastaFile
        pysam FastaFile object for the reference genome
        it is to be loaded like this: fasta = pysam.FastaFile('/Users/haochen/Desktop/Tapestri_analysis/resources/genome/ucsc_hg19.fa')
    
    Returns
    -------
    str
        Trinucleotide context for the snv. E.g. 'AAT>ATT'
    """
    # contig = snv.split(':')[0]
    # pos = int(snv.split(':')[1])
    ref = snv.split(':')[2].split('/')[0]
    alt = snv.split(':')[2].split('/')[1]
    if len(ref) != 1 or len(alt) != 1:
        # print('[WARNING] Only snvs are supported! Returning None.')
        return None

    # the region needs to be specified in samtools region string format (https://pysam.readthedocs.io/en/latest/glossary.html#term-region). E.g. 'chr1:100-200' (both ends would be 1-based and included)
    region = f'{snv.split(":")[0]}:{int(snv.split(":")[1])-1}-{int(snv.split(":")[1])+1}'

    ref_trinuc = fasta.fetch(region = region).upper()
    alt_trinuc = ref_trinuc[0] + alt + ref_trinuc[2]
    return ref_trinuc + '>' + alt_trinuc

# 05/18/2023 adapted from TRACERx Signatures-Code.ipynb
def create_snv_class(snv_trinuc_context: str) -> str:
    """
    Create the channels, in pyrimidine reference
    Args:
        snv in the form ref_trinuc + '>' + alt_trinuc (e.g. AAA>ATA)
    Returns:
        str: Fixed channels
    """
    
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N':'N'}
    if snv_trinuc_context[1] in pyr:
        out = '{}[{}>{}]{}'.format(snv_trinuc_context[0], snv_trinuc_context[1], snv_trinuc_context[5], snv_trinuc_context[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[snv_trinuc_context[2]], rev[snv_trinuc_context[1]], rev[snv_trinuc_context[5]], rev[snv_trinuc_context[0]])
    
    return out

def snvs_order():
    """
    Get correct order for plotting
    Returns:
        list: 
    """
    
    order = []
    first = ['A', 'C', 'G', 'T']
    pyr = ['C', 'T']
    for p in pyr:
        for mut in first:
            if mut != p:
                for f in first:
                    for f2 in first:
                        comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                        order.append(comb)
    return order

def chunks(l, n):
    """
    Split into even chunks
    Args:
        l (list): 
        n (int): number of chunks
    Yields:
        [list]: 
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]

def plot_snvs(sig, title, outpath='', fontsize = 12):
    """
    Plotting function
    Args:
        sig (list): vector, size 96 
        title (str): [description]
        outpath (str, optional): Defaults to ''.
        fontsize (int, optional): Defaults to 12.
    """

    config_params(fontsize)

    fig, axs = plt.subplots(
        nrows=2, ncols=1, figsize=(10.2, 3), gridspec_kw={'height_ratios': [1, 9]}
    )
    
    plt.title(title)
    order_plot = snvs_order()

    vals = []
    colors = []
    colors_mut = [
        '#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'
    ]
    bot = -0.5
    for ix, c in enumerate(chunks(sig, 16)):
        colors.extend([colors_mut[ix] for s in c])
        axs[0].barh(1, 16, left=bot, color=colors_mut[ix])
        bot += 16
        vals.extend(c)

    axs[0].set_xlim(-1, 96)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(
        ['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
        verticalalignment="center", ha='center', rotation=90, fontsize=6,
        color='grey'
    )

    plt.tight_layout()
    plt.xlim(-1, 96)

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Mutational Count')
    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    if len(outpath)>0:
        plt.savefig(outpath)
    plt.show()
