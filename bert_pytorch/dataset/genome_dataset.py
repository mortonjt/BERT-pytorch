from torch.utils.data import Dataset
import tqdm
import torch
import numpy as np
from Bio import SeqIO
import pandas as pd
import glob
import collections
from util import draw_exclusive, get_context, get_seq, GeneInterval


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

    while (genes[ridx].start - coord) < window_size and ridx < len(genes):
        ridx = ridx + 1

    return genes[lidx : idx] + genes[idx + 1 : ridx]


class ExtractIntervals(object):
    def __init__(self, window_size=10000):
        self.window_size = window_size

    def __call__(self, gb_file):
        gb_record = SeqIO.read(open(gb_file, "r"), "genbank")
        cds = list(filter(lambda x: x.type == 'CDS', gb_record.features))
        starts = list(map(lambda x: int(x.location.start), cds))
        ends = list(map(lambda x: int(x.location.end), cds))
        strand = list(map(lambda x: x.strand, cds))
        seqs = list(map(get_seq, cds))
        res = zip(starts, ends, seqs, strand)

        # sequences with start, end and position
        res = list(filter(lambda x: len(x) > 0, res))
        res = list(map(lambda x: GeneInterval(
            start=x[0], end=x[1], sequence=x[2], strand=x[3]
        ), res))
        return {'gene_intervals': res}


class SampleGenes(object):
    def __init__(self, num_sampled, within_prob=0.5, window_size=10000):
        """ Randomly samples genes within genome

        Parameters
        ----------
        num_sampled : int
            Number of genes sampled per genome
        within_prob : float
            The probability of drawing a gene within the same operon
        window_size : int
            The radius of the operon (typically around 10kb)
        """
        self.num_sampled = num_sampled
        self.within_prob = within_prob
        self.window_size = window_size

    def __call__(self, record):
        """ Randomly samples genes are their paired genes

        Parameters
        ----------
        record: dict
           key : 'gene_intervals'
           values : list of GeneInterval objects
        """
        gis = record['gene_intervals']
        # draw genes
        idx = np.random.randint(0, len(gis), size=self.num_sampled)
        draws = np.random.random(size=self.num_sampled)

        def draw_operon(x):
            operon = get_context(gis, x, self.window_size)
            # TODO: check this assumption carefully
            mid = len(operon) // 2
            print(operon)
            i = draw_exclusive(len(operon), mid)
            return gis[i]

        # draw context genes
        rand_operons = list(map(draw_operon, idx[draws<self.within_prob]))
        rand_pairs = list(map(lambda x: draw_exclusive(len(gis), x),
                              idx[draws>self.within_prob]))
        rand_pairs = list(map(lambda i : gis[i], rand_pairs))
        genes = list(map(lambda i: gis[i], idx))
        rand_pairs = rand_operons + rand_pairs


        return {'genes' : genes, 'next_genes' : rand_pairs}


class MaskPeptides(object):
    def __init__(self, mask_prob=0.8, swap_prob=0.25):
        self.mask_prob = mask_prob
        self.swap_prob = swap_prob

    def __call__(self):
        pass

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
        res = list(filter(lambda x: len(x) > 0, res))
        res = list(map(lambda x: GeneInterval(
            start=x[0], end=x[1], sequence=x[2]
        ), res))
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
            np.arange(len(peptides)),
            index=peptides
        )
        self.peptides = set(peptides)

    def __getitem__(self, x):
        if x in self.peptides:
            return self.peptide_lookup.loc[x]
        else:
            return self.mask

    def __contains__(self, x):
        return x in self.peptides
