import unittest

from aberration_multigraph.generator import AMGGenerator


class TestAMGGeneratorStructure(unittest.TestCase):
    """Structural and deterministic tests for AMGGenerator."""

    def test_label_count(self):
        gen = AMGGenerator(2, [1, 2])
        expected = 2 * (1 + 2) + 2 * 2
        self.assertEqual(len(gen.vertices), expected)

    def test_dsb_pairs_are_symmetric(self):
        gen = AMGGenerator(2, [2, 1])
        for u, v in gen.dsbs:
            self.assertEqual(gen.dsb_pair[u], v)
            self.assertEqual(gen.dsb_pair[v], u)

    def test_chromatin_edge_count(self):
        num_dsbs = [1, 0, 2]
        gen = AMGGenerator(3, num_dsbs)
        expected = sum(1 + d for d in num_dsbs)
        self.assertEqual(len(gen.chromatins), expected)


class TestAMGGeneratorGeneration(unittest.TestCase):
    """Tests for small, fully enumerable AMG instances."""

    def test_single_chromosome_single_dsb(self):
        gen = AMGGenerator(1, [1])
        amgs = list(gen.generate_amgs())

        self.assertEqual(len(amgs), 1)

        amg = amgs[0]
        self.assertEqual(len(amg.rejoins), 1)
        self.assertTrue(amg.is_connected())

    def test_no_rejoin_equals_dsb(self):
        gen = AMGGenerator(2, [1, 1])
        dsbs = {tuple(sorted(e)) for e in gen.dsbs}

        for amg in gen.generate_amgs():
            for u, v in amg.rejoins:
                self.assertNotIn(tuple(sorted((u, v))), dsbs)

    def test_all_generated_amgs_are_connected(self):
        gen = AMGGenerator(2, [2, 1])
        for amg in gen.generate_amgs():
            self.assertTrue(amg.is_connected())


class TestAMGGeneratorInvariants(unittest.TestCase):
    """Tests that fundamental AMG invariants always hold."""

    def test_rejoins_form_perfect_matching(self):
        gen = AMGGenerator(2, [2, 1])

        dsb_vertices = set()
        for u, v in gen.dsbs:
            dsb_vertices.add(u)
            dsb_vertices.add(v)

        for amg in gen.generate_amgs():
            used = set()
            for u, v in amg.rejoins:
                self.assertNotIn(u, used)
                self.assertNotIn(v, used)
                used.add(u)
                used.add(v)

            self.assertEqual(used, dsb_vertices)


class TestAMGGeneratorKnownCounts(unittest.TestCase):
    """Regression tests for very small known cases."""

    def test_known_small_counts(self):
        cases = [
            ((1, [0]), 0),
            ((1, [1]), 1),
            ((2, [1, 1]), 2),
            ((1, [3]), 8),
            ((3, [1, 1, 1]), 8),
            ((1, [4]), 60),
            ((2, [2, 2]), 56),
            ((4, [1, 1, 1, 1]), 48),
        ]

        for (n, dsbs), expected in cases:
            with self.subTest(num_chromosomes=n, num_dsbs=dsbs):
                gen = AMGGenerator(n, dsbs)
                count = sum(1 for _ in gen.generate_amgs())
                self.assertEqual(count, expected)


if __name__ == "__main__":
    unittest.main()