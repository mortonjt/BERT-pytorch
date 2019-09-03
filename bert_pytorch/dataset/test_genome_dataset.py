import os
import unittest
from genome_dataset import GeneInterval, ExtractIntervals


class TestExtractIntervals(unittest.TestCase):
    def setUp(self):
        self.record_name = os.path.join('./data/ichnovirus.gb')
    def test_extract(self):
        self.maxDiff = None
        int_tfm = ExtractIntervals()
        res = int_tfm(self.record_name)

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
