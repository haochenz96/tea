
from utils import *

label_for_var(samples[sample_i], ['chr12:25398284:C/A'])
samples[sample_i].cnv.normalize_reads()
samples[sample_i].cnv.compute_ploidy(diploid_cells = samples[sample_i].dna.barcodes('others'))

# set up labels and colors
samples[sample_i].cnv.set_labels(samples[sample_i].dna.get_labels())
samples[sample_i].cnv.set_palette(samples[sample_i].dna.get_palette())