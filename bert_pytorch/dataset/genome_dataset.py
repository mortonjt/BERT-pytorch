from torch.utils.data import Dataset
import tqdm
import torch
import random
from Bio import SeqIO
import pandas as pd
import glob

def get_seq(x):
    if 'translation' in x.qualifiers:
        return x.qualifiers['translation'][0]
    else:
        return ''




class GenomeDataset(Dataset):
    def __init__(self, genbank_directory, vocab, genbank_ext='.gb',
                 within_operon_prob=0.5, mask_prob=0.8):
        self.genbank_files = glob.glob(genbank_directory, '*' + genbank_ext)

        self.within_operon_prob = within_operon_prob
        self.mask_prob = mask_prob

    def __len__(self):
        return len(self.genbank_files)

    def __getitem__(self, item):
        # return a list of items for each genome
        gb_file = self.genbank_files[item]


        # get random gene
        i = random.randint(len(res))
        s, e, seq = res[i]

        # get random next gene

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
        pass

    def random_gene(self, genome):
        pass

    @staticmethod
    def read_genbank(f):
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
