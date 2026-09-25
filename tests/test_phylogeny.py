import unittest
from treeio import (Tree, build_tree, distance_matrix, neighbor_joining)


class _SeqRecord:
    """Duck-typed stand-in for a Biopython SeqRecord."""
    def __init__(self, sid, seq):
        self.id = sid
        self.seq = type("S", (), {"__str__": lambda self: seq})()


class TestPhylogeny(unittest.TestCase):
    def test_distance_matrix_p(self):
        labels, D = distance_matrix({"A": "AAAA", "B": "AAAT", "C": "AATT"}, model="p")
        self.assertEqual(labels, ["A", "B", "C"])
        self.assertAlmostEqual(D[0][1], 0.25)
        self.assertAlmostEqual(D[0][2], 0.5)
        self.assertAlmostEqual(D[1][2], 0.25)
        self.assertAlmostEqual(D[0][0], 0.0)
        self.assertTrue(all(abs(D[i][j] - D[j][i]) < 1e-9 for i in range(3) for j in range(3)))

    def test_distance_matrix_jc(self):
        _, D = distance_matrix({"A": "AAAA", "B": "AAAT"}, model="jc")
        # p = 0.25 -> jukes-cantor = -0.75 ln(1 - 4*0.25/3) = -0.75 ln(2/3)
        import math
        self.assertAlmostEqual(D[0][1], -0.75 * math.log(2.0 / 3.0), places=6)

    def test_build_tree_topology(self):
        seqs = {"A": "AAAA", "B": "AAAT", "C": "AATT"}
        t = build_tree(seqs, model="p")
        self.assertEqual(set(t.tip_names), {"A", "B", "C"})
        self.assertTrue(t.is_binary())
        self.assertEqual(t.nnodes, 5)
        # A-B are the closest pair -> they form a clade excl. C
        self.assertTrue(t.is_monophyletic("A", "B"))
        self.assertFalse(t.is_monophyletic("A", "C"))
        # branch lengths present
        self.assertGreater(t.get_tree_length(), 0)
        self.assertIn("(A", t.as_newick())

    def test_tree_from_alignment(self):
        t = Tree.from_alignment({"A": "AA", "B": "AT", "C": "TT"}, model="p")
        self.assertEqual(set(t.tip_names), {"A", "B", "C"})

    def test_seq_record_input(self):
        seqs = [
            _SeqRecord("human", "ATG"),
            _SeqRecord("chimp", "ATG"),
            _SeqRecord("mouse", "ATC"),
        ]
        t = build_tree(seqs)
        self.assertEqual(set(t.tip_names), {"human", "chimp", "mouse"})

    def test_pair_input(self):
        t = build_tree([("a", "AA"), ("b", "AA"), ("c", "CC")], model="p")
        self.assertEqual(set(t.tip_names), {"a", "b", "c"})

    def test_bad_model(self):
        with self.assertRaises(ValueError):
            distance_matrix({"A": "AA", "B": "AA"}, model="nope")

    def test_bad_matrix_size(self):
        with self.assertRaises(ValueError):
            neighbor_joining(["A", "B"], [[0.0, 1.0]])

    def test_single_taxon(self):
        t = build_tree({"A": "AA"})
        self.assertEqual(t.tip_names, ["A"])

    def test_gaps_ignored(self):
        _, D = distance_matrix({"A": "A-", "B": "A-", "C": "AT"}, model="p")
        # A vs C: only the first site is comparable (A==A), the gap site is ignored
        self.assertAlmostEqual(D[0][2], 0.0)
        # a real difference: A=AA vs C=AT -> 1 of 2 sites differ -> 0.5
        _, D2 = distance_matrix({"A": "AA", "B": "AA", "C": "AT"}, model="p")
        self.assertAlmostEqual(D2[0][2], 0.5)


if __name__ == "__main__":
    unittest.main()
