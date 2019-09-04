import numpy as np
import collections


GeneInterval = collections.namedtuple(
    'GeneInterval', ['start', 'end', 'sequence', 'strand']
)

def get_seq(x):
    if 'translation' in x.qualifiers:
        return x.qualifiers['translation'][0]
    else:
        return ''

def get_context(genes, idx, window_size):
    """ Retrieves context genes

    Parameters
    ----------
    genes : list of GeneInterval
        List of genes in genome
    idx : int
        Index for gene of interest
    """
    lidx, ridx = idx - 1, idx + 1
    context = []

    # only grab gene if it is in the same strand
    while lidx >= 0 and ((genes[idx].start - genes[lidx].end) < window_size):
        if genes[lidx].strand == genes[idx].strand:
            context.append(genes[lidx])
        lidx = lidx - 1
    while ridx < len(genes) and ((genes[ridx].start - genes[idx].end) < window_size):
        if genes[ridx].strand == genes[idx].strand:
            context.append(genes[ridx])
        ridx = ridx + 1
    return context

def draw_exclusive(n, idx):
    if n <= 1:
        raise ValueError('Cannot draw exclusively, n<=1')
    j = np.random.randint(n)
    while j == idx:
        j = np.random.randint(n)
    return j
