import math
import unittest

import matplotlib

matplotlib.use("Agg")

from treeio import (
    Tree,
    neighbor_joining,
    read_phyloxml,
    tree_coords,
    write_phyloxml,
)


class TestFoundBugs(unittest.TestCase):
    # Bug 1: PhyloXML reader injects a spurious anonymous ``unknown`` root node,
    # mangling topology (node count +1, real root demoted to a child).
    # Reproduces at treeio/phyloxml.py:50 (_parse_clade called on <phylogeny>)
    def test_phyloxml_no_extra_root(self):
        t = Tree("root")
        inner = Tree("inner")
        inner.append_child(Tree("A"))
        inner.append_child(Tree("B"))
        t.append_child(inner)
        t.append_child(Tree("C"))
        self.assertEqual(t.nnodes, 5)

        xml = write_phyloxml(t)
        # correct writer must not itself produce a bogus root; sanity check
        self.assertEqual(xml.count("<clade"), 5)

        t2 = read_phyloxml(xml)
        # EXPECTED: exactly the same 5-node topology
        self.assertEqual(t2.nnodes, 5)
        self.assertNotEqual(t2.name, "unknown")
        self.assertEqual(t2.name, "root")
        self.assertNotIn("unknown", {n.name for n in t2.traverse()})

    # Bug 2: NJ clamps negative limb lengths to 0 (phylogeny.py:181-182), so the
    # reconstructed tree no longer reproduces the input distance matrix.
    # For this matrix a valid NJ tree must reproduce all three patristic distances.
    def test_nj_reproduces_distance_matrix(self):
        labels = ["A", "B", "C"]
        D = [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 3.0],
            [1.0, 3.0, 0.0],
        ]
        t = neighbor_joining(labels, D)
        got = t.get_cophenetic_distance()
        # EXPECTED: each pair reproduces the input distance exactly
        for i in range(3):
            for j in range(i + 1, 3):
                self.assertAlmostEqual(
                    got[i][j],
                    D[i][j],
                    places=9,
                    msg=f"dist({labels[i]},{labels[j]}) not reproduced",
                )

    # Bug 3: ``daylight`` layout collapses balanced-tree tips onto the same
    # coordinate (plot.py:358-381, _daylight_coords). Tips must be distinct.
    def test_daylight_tips_distinct(self):
        t = Tree("r")
        t.append_child(Tree("A"))
        t.append_child(Tree("B"))
        t.append_child(Tree("C"))
        t.append_child(Tree("D"))
        c = tree_coords(t, layout="daylight")
        tips = [tip for tip in t.get_tips()]
        for i in range(len(tips)):
            for j in range(i + 1, len(tips)):
                x1, y1 = c[tips[i]]
                x2, y2 = c[tips[j]]
                self.assertGreater(
                    math.hypot(x1 - x2, y1 - y2),
                    1e-9,
                    msg=f"tips {tips[i].name} and {tips[j].name} overlap",
                )


if __name__ == "__main__":
    unittest.main()
