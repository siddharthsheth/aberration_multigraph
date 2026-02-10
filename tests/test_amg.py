import unittest, warnings
from aberration_multigraph.amg import AberrationMultigraph
import numpy as np

class TestAMG(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)

    def test_simple(self):
        amg = AberrationMultigraph([(3,4), (2,1), (5,6)],
                                   [(4,5), (3,2)],
                                   [(3,4), (5,2)])

        self.assertEqual(amg.chromatins, ((1,2), (3,4), (5,6)))
        self.assertEqual(amg.dsbs, ((2,3), (4,5)))
        self.assertEqual(amg.rejoins, ((2,5), (3,4)))

    def test_labeled_vertex(self):
        amg = AberrationMultigraph([('A','B'), ('C','D'), ('E','F')],
                                   [('B','C'), ('D','E')],
                                   [('B','E'), ('C','D')])

        self.assertEqual(amg.chromatins, (('A','B'), ('C','D'), ('E','F')))
        self.assertEqual(amg.dsbs, (('B','C'), ('D','E')))
        self.assertEqual(amg.rejoins, (('B','E'), ('C','D')))
    
    def test_num_chromosome(self):
        chromatin = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16)]
        dsb = [(2,3),(6,7),(10,11),(14,15)]
        rejoin = [(2,6),(3,11),(7,14),(10,15)]
        amg = AberrationMultigraph(chromatin, dsb, rejoin)

        self.assertEqual(amg.num_chromosome, 4)

    def test_diameter_connected(self):
        amg = AberrationMultigraph((('A','B'), ('C','D'), ('E','F'), ('G','H')),
                            (('B','C'), ('D','E'), ('F','G')),
                            (('B','E'), ('C','G'), ('D','F')))
        
        self.assertEqual(amg.diameter(), 4)

    def test_diameter_inf(self):
        chrom = (('A','B'), ('C','D'), ('E','F'), ('G','H'), ('I','J'), ('K','L'))
        dsb = (('B','C'), ('D','E'), ('H','I'), ('J','K'))
        rejoin = (('B', 'E'), ('C', 'D'), ('H', 'J'), ('I', 'K'))

        amg = AberrationMultigraph(chrom, dsb, rejoin)

        self.assertEqual(amg.diameter(), np.inf)

    def test_girth(self):
        amg = AberrationMultigraph((('A','B'), ('C','D'), ('E','F'), ('G','H')),
                            (('B','C'), ('D','E'), ('F','G')),
                            (('B','E'), ('C','G'), ('D','F')))
        
        self.assertEqual(amg.girth(), 3)

    def test_total_swap(self):
        chromatin = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16)]
        dsb = [(2,3),(6,7),(10,11),(14,15)]
        rejoin = [(2,6),(3,11),(7,14),(10,15)]

        amg = AberrationMultigraph(chromatin, dsb, rejoin)
        amg_2 = amg.total_swap(0,1)

        self.assertEqual(amg_2, AberrationMultigraph(chromatin,
                                                     dsb,
                                                     ((2,6), (3,14), (7,11), (10,15))))

    def test_total_twist(self):
        chromatin = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16)]
        dsb = [(2,3),(6,7),(10,11),(14,15)]
        rejoin = [(2,6),(3,11),(7,14),(10,15)]

        amg = AberrationMultigraph(chromatin, dsb, rejoin)
        amg_2 = amg.total_twist(1)

        self.assertEqual(amg_2, AberrationMultigraph(chromatin,
                                                     dsb,
                                                     ((2,7), (3,11), (6,14), (10,15))))

