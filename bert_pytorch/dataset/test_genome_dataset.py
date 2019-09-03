import os
import random
import numpy as np
from genome_dataset import GeneInterval, ExtractIntervals, SampleGenes

import unittest


class TestExtractIntervals(unittest.TestCase):

    def setUp(self):
        self.record_name = os.path.join('./data/ichnovirus.gb')

    def test_extract(self):
        self.maxDiff = None
        int_tfm = ExtractIntervals()
        res = int_tfm(self.record_name)['gene_intervals']

        exp = [
            GeneInterval(250, 796,
                         'MEIFPMDRLFKKNTMGNNIFHEIAIEGSLLMLRRVRDNVNEQMDTYLSD'
                         'TNDQGETCIVIVADRHRGHLAIELIEIFVGLGADINGTDNGGNTALHYT'
                         'VFNGDHALAEWLCQQPGINLNAANHDELTPLGLAIQLNIQSMKALLEAT'
                         'GAIFHDIESNDSDNDDDDDDDDDDDDDVSTRRHG', -1),
            GeneInterval(2240, 2810,
                         'MVNSCVLELFEGNTSAGNNIFHEIAMKGSLALLLEIRDKFDRPTDHALR'
                         'EWNGHGETCLHLVALMNRGQNAIRMIDILVELGADLNAKNHLGHTLLHY'
                         'ALENDDCELINWLLLHPEMNLSVRDYYDMQTDDDCFVEESEEEQEEEET'
                         'EETEEEEKTRVSFSAFSDDLMDFESDEFDDIPRWIDELVSIL', -1)
        ]
        self.assertListEqual(exp, res)


class TestSampleGenes(unittest.TestCase):

    def setUp(self):
        self.record_name = os.path.join('./data/polyoma.gb')

    def test_sample_genes(self):
        np.random.seed(0)
        samp_tfm = SampleGenes(num_sampled=2, within_prob = 0.9, window_size=100)

        records = [
            GeneInterval(
                268, 469,
                ('MVLRQLSRQASVRVSKTWTGTKRRAQRIFIFILELLLEFCRGEDSVDGKNKSTTALPA'
                 'VKDSVKDS'), 1
            ),
            GeneInterval(
                504, 1560,
                ('MGAALALLGDLVASVSEAAAATGFSVAEIAAGEAAAAIEVQIAS'
                 'LATVEGITSTSEAIAAIGLTPQTYAVIAGAPGAIAGFAALIQTVTGISSLAQVGYRFF'
                 'SDWDHKVSTVGLYQQSGMALELFNPDEYYDILFPGVNTFVNNIQYLDPRHWGPSLFAT'
                 'ISQALWHVIRDDIPAITSQELQRRTERFFRDSLARFLEETTWTIVNAPVNFYNYIQDY'
                 'YSNLSPIRPSMVRQVAEREGTQVNFGHTYRIDDADSIQEVTQRMELRNKENVHSGEFI'
                 'EKTIAPGGANQRTAPQWMLPLLLGLYGTVTPALEAYEDGPNQKKRRVSRGSSQKAKGT'
                 'RASAKTTNKRRSRSSRS'), 1
            ),
            GeneInterval(
                861, 1560,
                ('MALELFNPDEYYDILFPGVNTFVNNIQYLDPRHWGPSLFATISQ'
                 'ALWHVIRDDIPAITSQELQRRTERFFRDSLARFLEETTWTIVNAPVNFYNYIQDYYSN'
                 'LSPIRPSMVRQVAEREGTQVNFGHTYRIDDADSIQEVTQRMELRNKENVHSGEFIEKT'
                 'IAPGGANQRTAPQWMLPLLLGLYGTVTPALEAYEDGPNQKKRRVSRGSSQKAKGTRAS'
                 'AKTTNKRRSRSSRS'), 1
            )
        ]

        neighbor = GeneInterval(
            504, 1560,
            ('MGAALALLGDLVASVSEAAAATGFSVAEIAAGEAAAAIEVQIAS'
             'LATVEGITSTSEAIAAIGLTPQTYAVIAGAPGAIAGFAALIQTVTGISSLAQVGYRFF'
             'SDWDHKVSTVGLYQQSGMALELFNPDEYYDILFPGVNTFVNNIQYLDPRHWGPSLFAT'
             'ISQALWHVIRDDIPAITSQELQRRTERFFRDSLARFLEETTWTIVNAPVNFYNYIQDY'
             'YSNLSPIRPSMVRQVAEREGTQVNFGHTYRIDDADSIQEVTQRMELRNKENVHSGEFI'
             'EKTIAPGGANQRTAPQWMLPLLLGLYGTVTPALEAYEDGPNQKKRRVSRGSSQKAKGT'
             'RASAKTTNKRRSRSSRS'), 1
        )

        res = samp_tfm({'gene_intervals': records})
        self.assertIn(neighbor in res['genes'])
        self.assertIn(neighbor in res['next_genes'])

    def test_avoid_clash(self):
        pass


class TestGenomeDataset(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_operon(self):
        pass

    def test_get_item(self):
        pass

    def test_iter(self):
        pass

    def test_random_peptide(self):
        pass

    def test_random_gene(self):
        pass

    def test_read_genbank(self):
        pass


if __name__ == '__main__':
    unittest.main()
