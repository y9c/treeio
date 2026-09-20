#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Comprehensive tests for the core treeio engine, I/O and plotting."""

import unittest
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

from treeio import (
    Tree, read, read_many, write, detect_format,
    read_newick, write_newick, read_nexus, write_nexus,
    read_phyloxml, write_phyloxml, read_json, write_json,
    draw, draw_tree, render, treeplot, save,
    color_by_value, color_map, named_palette,
    attach, get_tipdata, get_nodedata,
)


N = "((raccoon:19.19959,bear:6.80041):0.84600," \
    "((sea_lion:11.99700,seal:12.00300):7.52973," \
    "((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382,dog:25.46154);"


def _small():
    # root with children: (A,(B,C))
    a = Tree("A", branch_length=0.2)
    b = Tree("B", branch_length=0.3)
    c = Tree("C", branch_length=0.4)
    r = Tree("root")
    r.children = [a, b, c]  # actually 3 children here; build proper below
    return r


def _animals():
    return read_newick(N)


class TestNewick(unittest.TestCase):
    def test_roundtrip(self):
        t = _animals()
        s = write_newick(t)
        self.assertEqual(read_newick(s), t)

    def test_tips_and_counts(self):
        t = _animals()
        self.assertEqual(t.tip_names, ["raccoon", "bear", "sea_lion", "seal", "monkey", "cat", "weasel", "dog"])
        self.assertEqual(t.nleaves, 8)
        self.assertEqual(t.ntips, 8)
        self.assertGreater(t.nnodes, 8)

    def test_tree_length(self):
        t = _animals()
        self.assertAlmostEqual(t.get_tree_length(), 277.27722, places=3)

    def test_quotes(self):
        t = read_newick("('A x':1,'B:2':3,C);")
        self.assertEqual(t.tip_names, ["A x", "B:2", "C"])

    def test_support(self):
        t = read_newick("((A:0.1,B:0.2)95:0.3,C:0.4);")
        supps = [n.support for n in t.get_internal_nodes()]
        self.assertIn(95.0, supps)

    def test_multiple_trees(self):
        trees = read_newick("(A,B);") and read("(A,B);(C,D);")
        # read returns first tree
        self.assertEqual(trees.tip_names, ["A", "B"])

    def test_multifurcation_not_binary(self):
        t = read_newick("(A,B,C);")
        self.assertFalse(t.is_binary())


class TestTreeOps(unittest.TestCase):
    def test_parent_setter(self):
        a, b = Tree("A"), Tree("B")
        b.parent = a
        self.assertIs(b.parent, a)
        self.assertIn(b, a.children)

    def test_append_remove(self):
        a, b = Tree("A"), Tree("B")
        a.append_child(b)
        self.assertIn(b, a.children)
        a.remove_child(b)
        self.assertNotIn(b, a.children)
        self.assertIsNone(b.parent)

    def test_traverse(self):
        t = read_newick("((A,B),(C,D));")
        pre = [n.name for n in t.traverse("preorder")]
        self.assertEqual(pre, ["unknown", "unknown", "A", "B", "unknown", "C", "D"])

    def test_mrca(self):
        t = read_newick("((A,B),(C,D));")
        n = t.get_mrca("A", "B")
        below = {x.name for x in n}
        self.assertEqual(below, {"unknown", "A", "B"})

    def test_cophenetic_distance(self):
        t = read_newick("((A:1,B:1):0.5,C:1);")
        mat = t.get_cophenetic_distance()
        tips = t.tip_names
        i, j = tips.index("A"), tips.index("B")
        self.assertAlmostEqual(mat[i][j], 2.0)
        ia, ic = tips.index("A"), tips.index("C")
        self.assertAlmostEqual(mat[ia][ic], 2.5)

    def test_ladderize(self):
        t = read_newick("((A,B),(C,D,E));")
        t.ladderize()
        # larger clade first
        self.assertEqual(len(t.children[0]), 3)

    def test_drop_tips(self):
        t = read_newick("((A,B),(C,D));")
        t.drop_tips(["D"])
        self.assertNotIn("D", t.tip_names)
        self.assertIn("A", t.tip_names)

    def test_prune(self):
        t = read_newick("((A,B),(C,D));")
        t.prune(["A", "B"])
        self.assertEqual(set(t.tip_names), {"A", "B"})

    def test_root(self):
        t = read_newick("((A,B),(C,D));")
        t.root("A")
        self.assertIn("A", t.children[0].tip_names if not t.children[0].is_leaf() else t.children[0].name)

    def test_resolve_polytomies(self):
        t = read_newick("(A,B,C,D);")
        t.resolve_polytomies()
        self.assertTrue(t.is_binary())

    def test_copy(self):
        t = read_newick("(A:1,(B:2,C:3));")
        c = t.copy()
        self.assertEqual(c, t)
        # independent
        c.get_tip_by_label("A").branch_length = 99
        self.assertNotEqual(t.get_tip_by_label("A").branch_length, 99)

    def test_is_monophyletic(self):
        t = read_newick("((A,B),(C,D));")
        self.assertTrue(t.is_monophyletic("A", "B"))
        self.assertFalse(t.is_monophyletic("A", "C"))

    def test_edge_vector(self):
        t = read_newick("(A,B);")
        edges = t.edge_vector()
        self.assertEqual(len(edges), 2)

    def test_deep_tree_no_recursion_error(self):
        # a very deep (ladder) tree must not hit the recursion limit
        root = Tree("r")
        cur = root
        for i in range(3000):
            tip = Tree(f"t{i}", branch_length=1.0)
            nxt = Tree()
            cur.children = [tip, nxt]
            cur = nxt
        cur.children = [Tree("last", branch_length=1.0)]
        self.assertEqual(root.nnodes, 3000 * 2 + 2)
        self.assertEqual(len(root.get_tips()), 3001)
        # deep copy + coordinate calc must also survive
        self.assertEqual(root.copy().nnodes, root.nnodes)
        root.get_y_and_x()
        self.assertIsNotNone(root)

    def test_label_index(self):
        t = read_newick("((A,B),(C,D));")
        idx = t.label_index()
        self.assertIs(idx["A"], t.get_tip_by_label("A"))
        self.assertIn("A", idx)

    def test_slots_and_annotations(self):
        t = read_newick("(A,B);")
        # nodes use __slots__ (no per-node __dict__)
        self.assertFalse(hasattr(t, "__dict__"))
        # annotation key that collides with a method name must still read back
        node = t.get_tip_by_label("A")
        node.set_data("root", "R")
        self.assertEqual(node.get_data("root"), "R")
        self.assertIsInstance(node.get_data("root"), str)
        # deep copy keeps the annotation
        self.assertEqual(t.copy().get_tip_by_label("A").get_data("root"), "R")
        self.assertIsNone(node.get_data("nonexistent"))


class TestIO(unittest.TestCase):
    def test_detect(self):
        self.assertEqual(detect_format(N), "newick")
        self.assertEqual(detect_format("#NEXUS\n"), "nexus")
        self.assertEqual(detect_format("{\"a\":1}"), "json")
        self.assertEqual(detect_format("<?xml..."), "phyloxml")

    def test_read_write_newick_file(self):
        import tempfile, os
        d = tempfile.mkdtemp()
        src = os.path.join(d, "in.nwk")
        with open(src, "w") as fh:
            fh.write(N)
        t = read(src)
        self.assertEqual(t.tip_names, _animals().tip_names)

    def test_convert_json(self):
        t = _animals()
        s = write(t, format="json")
        t2 = read_json(s)
        self.assertIsInstance(t2, list)
        self.assertEqual(t2[0].tip_names, t.tip_names)
        # and json->newick
        t3 = read(s, format="json")
        self.assertEqual(t3.tip_names, t.tip_names)

    def test_convert_nexus(self):
        t = _animals()
        s = write_nexus(t)
        t2 = read_nexus(s)
        self.assertEqual(t2.tip_names, t.tip_names)
        # topology + branch lengths preserved through the TRANSLATE table
        self.assertEqual(t2, t)
        self.assertAlmostEqual(t2.get_tree_length(), t.get_tree_length(), places=3)

    def test_convert_phyloxml(self):
        t = _animals()
        s = write_phyloxml(t)
        t2 = read_phyloxml(s)
        self.assertEqual(t2.tip_names, t.tip_names)

    def test_annotations_roundtrip(self):
        ann = read_newick("((A:1[&rate=0.5],B:1[&rate=0.7])[&height=2.3],C:1);", annotations=True)
        inner = ann.get_internal_nodes()
        self.assertAlmostEqual(inner[0].height, 2.3)
        out = write_newick(ann, annotations=["rate", "height"])
        ann2 = read_newick(out, annotations=True)
        self.assertAlmostEqual(ann2.get_internal_nodes()[0].height, 2.3)


class TestAttach(unittest.TestCase):
    def test_attach(self):
        t = _animals()
        attach(t, {"raccoon": "mammal", "dog": "mammal", "sea_lion": "mammal"}, key="class")
        self.assertEqual(get_tipdata(t, "class")["raccoon"], "mammal")
        self.assertNotIn("A", get_tipdata(t, "class"))

    def test_attach_rows(self):
        t = read_newick("(A,B);")
        attach(t, [{"name": "A", "size": 10}, {"name": "B", "size": 20}])
        self.assertEqual(get_tipdata(t, "size")["A"], 10)
        self.assertEqual(get_tipdata(t, "size")["B"], 20)


class TestPlot(unittest.TestCase):
    def setUp(self):
        import matplotlib
        matplotlib.use("Agg")

    def _ax_ok(self, ax):
        self.assertIsNotNone(ax)
        self.assertTrue(len(ax.collections) >= 1)
        return True

    def test_rectangular(self):
        t = _animals()
        ax = draw(t, layout="rectangular", tip_labels=True, node_support=True)
        self.assertTrue(self._ax_ok(ax))

    def test_circular(self):
        t = _animals()
        self.assertTrue(self._ax_ok(draw(t, layout="circular")))

    def test_slanted(self):
        t = _animals()
        self.assertTrue(self._ax_ok(draw(t, layout="slanted")))

    def test_no_branch(self):
        t = read_newick("(A,B,(C,D));")
        self.assertTrue(self._ax_ok(draw(t)))

    def test_color_by_value(self):
        cmap = color_by_value({"A": 1, "B": 2, "C": 3})
        self.assertEqual(set(cmap), {"A", "B", "C"})
        self.assertTrue(all(v.startswith("#") for v in cmap.values()))

    def test_layouts_valid(self):
        t = _animals()
        for layout in ("rectangular", "roundrect", "circular", "fan", "radial", "slanted", "unrooted", "time"):
            self.assertTrue(self._ax_ok(draw_tree(t, layout=layout)))

    def test_render_save_default(self):
        import tempfile, os
        t = _animals()
        ax = render(t, backend="mpl", layout="circular")
        self.assertTrue(self._ax_ok(ax))
        p = os.path.join(tempfile.mkdtemp(), "t.png")
        self.assertEqual(render(t, backend="mpl", path=p, layout="rectangular"), p)
        self.assertTrue(os.path.exists(p))

    def test_treeplot_fluent(self):
        t = _animals()
        tp = treeplot(t, layout="rectangular")
        tp.plot_edges()
        tp.scatter_tips()
        tp.plot_tip_labels()
        tp.plot_node_labels()
        self.assertTrue(self._ax_ok(tp.ax))

    def test_color_map_and_palette(self):
        self.assertEqual(color_map({"A": 1}), color_by_value({"A": 1}))
        pal = named_palette(["A", "B", "C"])
        self.assertEqual(set(pal), {"A", "B", "C"})
        self.assertEqual(pal["A"], pal["A"])

    def test_tree_save(self):
        import tempfile, os
        t = _animals()
        p = os.path.join(tempfile.mkdtemp(), "t.png")
        self.assertEqual(t.save(p, layout="circular"), p)
        self.assertTrue(os.path.exists(p))

    def test_branch_continuous_color(self):
        t = _animals()
        for n in t.get_nodes():
            n.set_data("rate", n.branch_length if n.branch_length is not None else 0.0)
        ax = draw_tree(t, layout="rectangular", tip_labels=True,
                       branch_color_field="rate", show_colorbar=True, colorbar_label="rate")
        self.assertIsNotNone(ax)
        # a colorbar is added as an extra axes on the figure
        self.assertGreaterEqual(len(ax.figure.axes), 2)

    def test_treeplot_colorbar_and_theme(self):
        t = _animals()
        for n in t.get_nodes():
            n.set_data("rate", n.branch_length if n.branch_length is not None else 0.0)
        tp = treeplot(t, layout="rectangular")
        tp.theme("clean")
        tp.plot_edges(aes="rate", cmap="viridis")
        tp.scatter_tips()
        tp.plot_tip_labels()
        cb = tp.add_colorbar("rate")
        self.assertIsNotNone(cb)
        self.assertGreaterEqual(len(tp.ax.figure.axes), 2)

    def test_scalebar(self):
        t = _animals()
        tp = treeplot(t, layout="rectangular")
        tp.plot_edges()
        tp.scale_bar(label="substitutions")
        texts = [x.get_text() for x in tp.ax.texts]
        self.assertTrue(any("substitutions" in s for s in texts))

    def test_theme_styles(self):
        import matplotlib.pyplot as plt
        from treeio import tree_theme
        ax = plt.subplots()[1]
        for style in ("clean", "void", "minimal", "plain"):
            tree_theme(ax, style=style)
        plt.close(ax.figure)

    def test_facet_grid(self):
        import matplotlib.pyplot as plt
        from treeio import facet_grid
        t = _animals()
        attach(t, {n.name: i for i, n in enumerate(t.get_tips())}, key="size")
        attach(t, {"raccoon": "m", "bear": "m", "monkey": "p", "dog": "p"}, key="grp")
        fig = facet_grid(t, "size", "grp", tip_labels=True)
        self.assertGreaterEqual(len(fig.axes), 3)  # 3 panels (+ optional colorbar)
        plt.close(fig)

    def test_grid_of_trees_shared_y(self):
        import matplotlib.pyplot as plt
        from treeio import grid_of_trees
        t = _animals()
        fig = grid_of_trees([t, t.copy()], labels=["a", "b"], layout="rectangular", tip_labels=True)
        self.assertEqual(len(fig.axes), 2)
        self.assertIsNotNone(fig.axes[0].get_xlim())
        plt.close(fig)

    def test_tip_color_field_categorical(self):
        import matplotlib.pyplot as plt
        t = _animals()
        attach(t, {"raccoon": "m", "monkey": "p", "dog": "p"}, key="class")
        ax = draw_tree(t, layout="rectangular", tip_labels=True, tip_color_field="class")
        # legend added as a separate axes (bbox_to_anchor) or legend object
        self.assertTrue(ax.get_legend() is not None)
        plt.close(ax.figure)

    def test_tip_color_field_continuous(self):
        import matplotlib.pyplot as plt
        t = _animals()
        attach(t, {n.name: (n.branch_length or 0) for n in t.get_tips()}, key="size")
        ax = draw_tree(t, layout="rectangular", tip_points=True, tip_labels=False, tip_color_field="size")
        self.assertGreaterEqual(len(ax.figure.axes), 2)  # colorbar axes
        plt.close(ax.figure)

    def test_polar_tip_label_rotation(self):
        t = _animals()
        ax = draw_tree(t, layout="circular", tip_labels=True)
        # labels are rotated on a radial tree
        self.assertTrue(any(text.get_rotation() for text in ax.texts))
        plt.close(ax.figure)

    def test_grid_of_trees_multirow_layouts(self):
        import matplotlib.pyplot as plt
        from treeio import grid_of_trees
        t = _animals()
        fig = grid_of_trees([t, t.copy(), t.copy()], labels=["a", "b", "c"],
                            layouts=["rectangular", "circular", "slanted"], ncols=2)
        # 2 cols x 2 rows (3 trees + 1 hidden)
        self.assertEqual(len(fig.axes), 4)
        plt.close(fig)

    def test_save_publication_dpi(self):
        import tempfile, os
        t = _animals()
        p = os.path.join(tempfile.mkdtemp(), "t.png")
        ret = render(t, layout="rectangular", path=p, dpi=150)
        self.assertEqual(ret, p)
        self.assertTrue(os.path.exists(p))
        # Tree.save passes dpi through
        p2 = os.path.join(tempfile.mkdtemp(), "t.svg")
        t.save(p2, layout="circular", format="svg", dpi=150)
        self.assertTrue(os.path.exists(p2))


if __name__ == "__main__":
    unittest.main()
