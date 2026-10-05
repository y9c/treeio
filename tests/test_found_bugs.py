import math
import unittest

import matplotlib

matplotlib.use("Agg")

from treeio import (
    Tree,
    detect_format,
    neighbor_joining,
    read,
    read_many,
    read_newick,
    read_phyloxml,
    tree_coords,
    write,
    write_nexus,
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

    # -- fourth bug-hunt pass (fresh multi-agent sweep on fixed code) --
    # root() must not change each tip-to-tip patristic distance (reroot is a
    # root-placement change on a FIXED unrooted topology).  Reproduces at
    # tree.py:790-821 (path-collapse dropsinternal nodes, flipping the split).
    def test_root_preserves_tip_distances(self):
        root = Tree("root")
        X = Tree("X", 3.0)
        A = Tree("A", 4.0)
        B = Tree("B", 1.0)
        C = Tree("C", 2.0)
        root.children = [X, A]
        X.children = [B, C]

        tips = root.tip_names
        pre = root.get_cophenetic_distance()
        idx = {name: i for i, name in enumerate(tips)}
        pre_dist = {(a, b): pre[idx[a]][idx[b]] for a in tips for b in tips}

        root.root("B")

        post_tips = root.tip_names
        post = root.get_cophenetic_distance()
        pidx = {name: i for i, name in enumerate(post_tips)}
        post_dist = {
            (a, b): post[pidx[a]][pidx[b]] for a in post_tips for b in post_tips
        }

        for a, b in sorted({tuple(sorted(p)) for p in pre_dist}):
            self.assertAlmostEqual(
                post_dist[(a, b)],
                pre_dist[(a, b)],
                places=9,
                msg=f"dist({a},{b}) changed after reroot",
            )

    # NeXML writer must preserve per-node set_data annotations; it currently
    # emits only the 'name' and 'support' metas (nexml.py:145-193), so trait
    # data is silently lost on a write->read round-trip.
    def test_nexml_preserves_annotations(self):
        t = read_newick("(Homo:0.1,Pan:0.2):0.3;")
        leaf = t.get_tips()[0]
        leaf.set_data("rate", 0.5)
        leaf.set_data("height", 2.0)
        leaf.set_data("is_model", True)
        internal = next(n for n in t.traverse("preorder") if not n.is_leaf())
        internal.set_data("sci_name", "Hominidae")

        t2 = read(write(t, format="nexml"), format="nexml")
        leaf2 = t2.get_tips()[0]
        internal2 = next(n for n in t2.traverse("preorder") if not n.is_leaf())

        self.assertEqual(leaf2.get_data("rate"), 0.5)
        self.assertEqual(leaf2.get_data("height"), 2.0)
        self.assertIs(leaf2.get_data("is_model"), True)
        self.assertEqual(internal2.get_data("sci_name"), "Hominidae")

    # _tip_angles hard-codes the circular base angle for every non-fan layout
    # (plot.py:1038-1046), so for 'unrooted' the ring/heatmap wedges are rotated
    # ~90 deg from where the tips actually sit.  Each tip's wedge must align
    # with the tip's own angle (as 'circular' does).
    def test_unrooted_ring_aligns_with_tips(self):
        from treeio import tree_coords
        from treeio.plot import _tip_angles

        t = read_newick("(" + ",".join(f"t{i}:1" for i in range(6)) + ");")
        n = t.nleaves
        tips = t.get_tips()
        half = 180.0 / n  # circular convention: tip at wedge start
        coords = tree_coords(t, layout="unrooted")
        for i, tip in enumerate(tips):
            x, y = coords[tip]
            ta = math.degrees(math.atan2(y, x)) % 360.0
            a0, a1 = _tip_angles(t, i, n, "unrooted")
            center = (a0 + ((a1 - a0) % 360.0) / 2.0) % 360.0
            rot = (center - ta - half) % 360.0
            # For 'unrooted', the ring must align like 'circular' (rot ~ 0).
            self.assertLess(abs(min(rot, rot - 360.0)), 1.0)

    # -- fifth pass: secondary defects fixed by the focused sweep --
    # detect_format must not treat a one-node JSON tree named "tree" as jplace,
    # nor a multi-tree JSON array as newick (io.py:118-136 now parses structure).
    def test_json_not_misdetected(self):
        one = '{"name":"tree","branch_length":0,"children":[]}'
        self.assertEqual(detect_format(one), "json")
        arr = "[" + one + "," + one + "]"
        self.assertEqual(detect_format(arr), "json")
        self.assertEqual(len(read_many(arr)), 2)
        # jplace still detected by its "tree" key
        self.assertEqual(
            detect_format('{"tree":"(A,B);","fields":[],"placements":[]}'), "jplace"
        )

    # write_nexus must not crash on a None annotation value (nexus.py:169 now
    # skips None, matching newick.py's guard).
    def test_nexus_none_annotation_ok(self):
        t = read_newick("(A:1,B:1);")
        t.get_tip_by_label("A").set_data("foo", None)
        out = write_nexus(t)
        self.assertIsInstance(out, str)
        self.assertIn("BEGIN TREES", out)

    # ladderize / is_binary / prune must not overflow the stack on a deep tree
    # (tree.py rewrote them iteratively).
    def test_deep_tree_methods_iterative(self):
        root = Tree("r")
        cur = root
        for i in range(2200):
            tip = Tree(f"t{i}")
            nxt = Tree()
            cur.children = [tip, nxt]
            cur = nxt
        tip_deep = Tree("deep")
        cur.append_child(tip_deep)
        root.ladderize()
        self.assertIsInstance(root.is_binary(), bool)
        pruned = root.copy()
        pruned.prune(["deep"])
        self.assertEqual(pruned.nnodes, 1)

    # reverse_x must actually affect equal_angle/daylight coords (plot.py wired it up)
    def test_eq_angle_reverse_x_used(self):
        from treeio import tree_coords

        t = read_newick("(A:1,B:1,C:1,D:1);")
        c1 = tree_coords(t, layout="equal_angle")
        c2 = tree_coords(t, layout="equal_angle", reverse_x=True)
        tips = t.get_tips()
        self.assertTrue(
            any(c1[tip] != c2[tip] for tip in tips), msg="reverse_x had no effect"
        )

    # ellipse layout must be genuinely distinct from circular
    def test_ellipse_distinct_from_circular(self):
        from treeio import tree_coords

        t = read_newick("(A:1,B:1,C:1,D:1);")
        circ = tree_coords(t, layout="circular")
        ell = tree_coords(t, layout="ellipse")
        tips = t.get_tips()
        self.assertTrue(
            any(circ[tip] != ell[tip] for tip in tips), msg="ellipse==circular"
        )


if __name__ == "__main__":
    unittest.main()
