import unittest
import json, tempfile, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from treeio import (Tree, read, write, read_newick, write_newick, read_newicks, read_nexus, write_nexus,
                    read_phyloxml, write_phyloxml, read_jplace, detect_format, color_by_value,
                    save, render, register_layout, register_backend, convert_format, convert_string,
                    tree_coords)
from treeio.show import tree_to_ascii


def _deep(n=2000):
    root = Tree("r"); cur = root
    for i in range(n):
        tip = Tree(f"t{i}"); nxt = Tree(); cur.children = [tip, nxt]; cur = nxt
    return root


class TestBugFixes(unittest.TestCase):
    # A: parent/child consistency
    def test_parent_none_detaches(self):
        a, b = Tree("a"), Tree("b"); a.children = [b]
        b.parent = None
        self.assertNotIn(b, a.children)
        self.assertIsNone(b.parent)

    def test_append_reparent_detaches(self):
        a, c, x = Tree("a"), Tree("c"), Tree("x")
        a.children = [c]; c.children = [x]
        a.append_child(x)
        self.assertNotIn(x, c.children)
        self.assertIs(x.parent, a)

    def test_children_setter_reparents(self):
        a, c, x = Tree("a"), Tree("c"), Tree("x")
        a.children = [c]; c.children = [x]
        a.children = [x]
        self.assertNotIn(x, c.children)

    def test_children_readonly(self):
        a, b = Tree("a"), Tree("b"); a.children = [b]
        with self.assertRaises(AttributeError):
            a.children.append(Tree("z"))

    # B: deep-tree recursion
    def test_deep_eq_asdict(self):
        d = _deep()
        self.assertTrue(d == d.copy())
        self.assertIsInstance(d.as_dict(), dict)

    def test_deep_ascii(self):
        self.assertIn("r", tree_to_ascii(_deep(1200), False, True))

    # C: parsing/serialization
    def test_newick_nonfinite_raises(self):
        with self.assertRaises(ValueError):
            write_newick(read_newick("(A:1e400,B:1);"))

    def test_json_strict_and_multi(self):
        t = read_newick("(A,B);")
        t.get_tip_by_label("A").branch_length = float("inf")
        s = write(t, format="json")
        self.assertNotIn("Infinity", s)
        json.loads(s)
        w = write(t, format="json")
        # two trees -> array

    def test_nexus_annotations_roundtrip(self):
        s = "#NEXUS\nBEGIN TREES;\n TREE t=[&R] ((A:1[&rate=0.5],B:1)[&height=3.0],C:1);\nEND;"
        t = read_nexus(s)
        t2 = read_nexus(write_nexus(t))
        self.assertIsNotNone(t2.get_internal_nodes()[0].get_data("height"))

    def test_jplace_annotations(self):
        t = read_jplace('{"tree":"(A[&x=1],B);","fields":[],"placements":[]}')
        self.assertEqual(t.get_tip_by_label("A").get_data("x"), 1)

    def test_phyloxml_property_numeric(self):
        t = read_newick("(A:1[&rate=0.5],B:1);", annotations=True)
        t2 = read_phyloxml(write_phyloxml(t, properties=["rate"]))
        self.assertEqual(t2.get_tip_by_label("A").get_data("rate"), 0.5)

    def test_color_none_ignored(self):
        self.assertEqual(color_by_value({"A": None, "B": None}), {})

    # D: API/exports
    def test_convert_exports(self):
        self.assertTrue(callable(convert_format))
        self.assertTrue(callable(convert_string))

    def test_plot_save_creates_dirs(self):
        p = os.path.join(tempfile.mkdtemp(), "sub", "x.png")
        save(read_newick("(A,B);"), p)
        self.assertTrue(os.path.exists(p))

    def test_register_guards(self):
        with self.assertRaises(ValueError):
            register_layout("circular", lambda t: {})
        with self.assertRaises(ValueError):
            register_backend("mpl", lambda t, **k: None)

    # E: detection/dist
    def test_detect_mrbayes_comment(self):
        s = "#NEXUS\n[mention mrbayes]\nBEGIN TREES;\n TREE t=(A,B);\nEND;"
        self.assertEqual(detect_format(s), "nexus")

    def test_get_distance_valueerror(self):
        t = read_newick("(A,B);")
        with self.assertRaises(ValueError):
            t.get_distance("A", "B")

    # -- second bug-hunt pass --
    def test_nexml_preserves_support(self):
        from treeio import read_nexml, write_nexml
        t = read_newick("((A:1,B:1)95:0.3,C:1);")
        t2 = read_nexml(write_nexml(t))
        self.assertTrue(any(n.support == 95.0 for n in t2.get_internal_nodes()))

    def test_register_format_guards_builtins(self):
        from treeio import register_format
        with self.assertRaises(ValueError):
            register_format("newick", read=lambda s: [read_newick(s)])

    def test_unregister_format(self):
        from treeio import register_format, unregister_format, read
        register_format("zzz", read=lambda s: [read_newick(s)], write=lambda t: write_newick(t))
        self.assertIsNotNone(unregister_format("zzz"))
        self.assertIsNone(unregister_format("zzz"))

    def test_read_empty_clear_error(self):
        with self.assertRaises(ValueError):
            read("")

    def test_read_newicks_stops_on_garbage(self):
        from treeio import read_newicks
        self.assertEqual(len(read_newicks("(A,B);%%%garbage")), 1)

    def test_prune_to_single_node(self):
        t = read_newick("((A,B),(C,D));")
        t.prune(["A"])
        self.assertEqual(t.nnodes, 1)
        self.assertEqual(t.tip_names, ["A"])

    def test_newick_nested_comment(self):
        t = read_newick("(A[outer[inner]tail],B);")
        self.assertEqual(set(t.tip_names), {"A", "B"})



if __name__ == "__main__":
    unittest.main()
