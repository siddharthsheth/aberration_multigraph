import unittest

from aberration_multigraph.incomplete_amg import IncompleteAMG
from aberration_multigraph.generator import AMGGenerator


class TestIncompleteAMGConstruction(unittest.TestCase):
    """Basic construction and sanity checks."""

    def test_empty_rejoin_is_allowed(self):
        chromatin = [(0, 1), (2, 3)]
        dsbs = [(1, 2)]
        rejoin = []

        inc = IncompleteAMG(chromatin, dsbs, rejoin)
        self.assertIsNotNone(inc)

    def test_even_number_of_free_vertices_required(self):
        chromatin = [(0, 1), (2, 3), (4, 5)]
        dsbs = [(1, 2), (3, 4)]
        rejoin = [(1, 2)]  # leaves 3,4 free

        inc = IncompleteAMG(chromatin, dsbs, rejoin)
        self.assertIsNotNone(inc)


class TestIncompleteAMGCounting(unittest.TestCase):
    """Tests for counting completions of incomplete AMGs."""

    def test_single_dsb_completion_count(self):
        chromatin = [(0, 1), (2, 3)]
        dsbs = [(1, 2)]
        rejoin = []

        inc = IncompleteAMG(chromatin, dsbs, rejoin)
        self.assertEqual(inc.count_amgs(), 1)

    def test_partial_rejoin_reduces_count(self):
        chromatin = [(0, 1), (2, 3), (4, 5)]
        dsbs = [(1, 2), (3, 4)]
        rejoin = [(1, 4)]

        inc = IncompleteAMG(chromatin, dsbs, rejoin)
        self.assertEqual(inc.count_amgs(), 1)


class TestIncompleteAMGEnumeration(unittest.TestCase):
    """Tests for enumerating complete AMGs from an incomplete AMG."""

    def test_complete_amgs_yields_connected_amgs(self):
        chromatin = [(0, 1), (2, 3), (4, 5)]
        dsbs = [(1, 2), (3, 4)]
        rejoin = []

        inc = IncompleteAMG(chromatin, dsbs, rejoin)

        for amg in inc.complete_amgs():
            self.assertTrue(amg.is_connected())

    def test_all_completions_are_distinct(self):
        chromatin = [(0, 1), (2, 3), (4, 5)]
        dsbs = [(1, 2), (3, 4)]
        rejoin = []

        inc = IncompleteAMG(chromatin, dsbs, rejoin)

        seen = set()
        for amg in inc.complete_amgs():
            key = frozenset(amg.rejoins)
            self.assertNotIn(key, seen)
            seen.add(key)


class TestIncompleteAMGConsistencyWithGenerator(unittest.TestCase):
    """Cross-check incomplete AMG results against AMGGenerator."""

    def test_count_matches_generator(self):
        num_chromosomes = 2
        num_dsbs = [1, 1]

        gen = AMGGenerator(num_chromosomes, num_dsbs)
        amgs = list(gen.generate_amgs())

        # Build an equivalent incomplete AMG
        chromatin = gen.chromatins
        dsbs = gen.dsbs
        rejoin = []

        inc = IncompleteAMG(chromatin, dsbs, rejoin)

        self.assertEqual(inc.count_amgs(), len(amgs))

    def test_enumeration_matches_generator(self):
        num_chromosomes = 2
        num_dsbs = [1, 1]

        gen = AMGGenerator(num_chromosomes, num_dsbs)
        gen_amgs = {
            frozenset(amg.rejoins) for amg in gen.generate_amgs()
        }

        inc = IncompleteAMG(gen.chromatins, gen.dsbs, [])

        inc_amgs = {
            frozenset(amg.rejoins) for amg in inc.complete_amgs()
        }

        self.assertEqual(gen_amgs, inc_amgs)


class TestIncompleteAMGInvariants(unittest.TestCase):
    """Invariant checks for all generated completions."""

    def test_rejoins_form_perfect_matching(self):
        chromatin = [(0, 1), (2, 3), (4, 5)]
        dsbs = [(1, 2), (3, 4)]
        rejoin = []

        inc = IncompleteAMG(chromatin, dsbs, rejoin)

        dsb_vertices = set()
        for u, v in dsbs:
            dsb_vertices.add(u)
            dsb_vertices.add(v)

        for amg in inc.complete_amgs():
            used = set()
            for u, v in amg.rejoins:
                self.assertNotIn(u, used)
                self.assertNotIn(v, used)
                used.add(u)
                used.add(v)

            self.assertEqual(used, dsb_vertices)


if __name__ == "__main__":
    unittest.main()