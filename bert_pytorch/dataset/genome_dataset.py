from torch.utils.data import Dataset
import tqdm
import torch
import random
from Bio import SeqIO
import pandas as pd
import glob
import collections


GeneInterval = collections.namedtuple(
    'GeneInterval', ['start', 'end', 'sequence']
)

def get_seq(x):
    if 'translation' in x.qualifiers:
        return x.qualifiers['translation'][0]
    else:
        return ''

def get_operon(genes, idx, window_size):
    """ Retrieves genes within the window size surrounding gene

    Parameters
    ----------
    genes : list of GeneInterval
        List of genes
    idx : int
        Index of gene of interest
    window_size : int
        Size of window surrounding gene.
    """
    s, e = genes[idx].start, genes[idx].end
    coord = (s + e) / 2
    lidx = max(0, idx - 1)
    ridx = min(idx + 1, len(genes))

    while (coord - genes[lidx].end) < window_size and lidx > 0:
        lidx = lidx - 1

    while genes[ridx].start - coord) < window_size and ridx < len(genes):
        ridx = ridx + 1

    return genes[lidx : idx] + genes[idx + 1 : ridx]


class GenomeDataset(Dataset):
    def __init__(self, genbank_directory, vocab, genbank_ext='.gb',
                 within_operon_prob=0.5, mask_prob=0.8, skipgram_size=10000):

        self.genbank_files = glob.glob(genbank_directory, '*' + genbank_ext)
        self.vocab = vocab.mask
        self.within_operon_prob = within_operon_prob
        self.mask_prob = mask_prob
        self.skipgram_size = skipgram_size

    def __len__(self):
        return len(self.genbank_files)

    def __getitem__(self, item):
        # return a list of items for each genome
        gb_file = self.genbank_files[item]
        record = GenomeDataset.read_genbank(gb_file)

        # get random gene and another context gene
        g1, g2, is_in_context = self.random.gene(genes)

        # get random peptide within gene


    def __iter__(self):
        start = 0
        end = len(self)
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None:  # single-process data loading
            for i in range(end):
                res = self[i]
                for k, v in res:
                    yield k, v

    def random_peptide(self, gene):
        seq = list(gene.sequence)
        output_label = []
        for i, pep in enumerate(seq):
            prob = random.random()
            if prob < self.mask_prob:
                seq[i] = self.vocab.mask
                output_label.append(1)
            else:
                output_label.append(0)
        return seq, output_label

    def random_gene(self, genes):
        """ Retrieve random gene and a pair

        Parameters
        ----------
        genes : list of GeneInterval
            List of genes
        """
        idx = random.randint(len(genes))
        operon = get_operon(genes, idx, window_size=self.window_size)
        g1 = genes[idx]
        if random.random() < within_operon_prob:
            jdx = random.randint(len(operon))
            return genes[idx], genes[jdx], 1
        else:
            jdx = random.randint(len(genes))
            return genes[idx], genes[jdx], 0

    @staticmethod
    def read_genbank(gb_file):
        """ """
        gb_record = SeqIO.read(open(gb_file, "r"), "genbank")
        cds = list(filter(lambda x: x.type == 'CDS', gb_record.features))
        starts = list(map(lambda x: int(x.location.start), cds))
        ends = list(map(lambda x: int(x.location.end), cds))
        strand = list(map(lambda x: x.strand, cds))
        seqs = list(map(get_seq, cds))
        res = zip(starts, ends, seqs)

        # TODO: will also need to factor in reverse complement

        # sequences with start, end and position
        res = list(filter(lambda x: len(x) > 0))
        res = list(map(lambda x: GeneInterval(
            start=x[0], end=x[1], sequence=x[2]
        )))
        return res


class PeptideVocab():

    def __init__(self, peptide_size=1):
        # enumerate all peptides
        peptides = ['A', 'C', 'D', 'E', 'F', 'G', 'H',
                    'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                    'R', 'S', 'T', 'V', 'W', 'Y']
        # this is the missing peptide by IUPAC standards
        # see https://www.bioinformatics.org/sms2/iupac.html
        self.mask = 'X'
        self.peptide_lookup = pd.Series(
            np.arange(len(peptides))
            index=peptides,
        )
        self.peptides = set(peptides)

    def __getitem__(self, x):
        if x in self.peptides:
            return self.peptide_lookup.loc[x]
        else:
            return self.mask

    def __contains__(self, x):
        return x in self.peptides
