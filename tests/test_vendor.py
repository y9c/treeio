import unittest
import tempfile, os
from treeio import (read, write,
                    read_nexml, write_nexml, read_mrbayes, read_iqtree, read_raxml, Tree)


class TestVendorFormats(unittest.TestCase):
    def test_nexml_roundtrip(self):
        t = read("(A:1,B:2):3;")
        xml = write_nexml(t)
        t2 = read_nexml(xml)
        self.assertEqual(set(t2.tip_names), {"A", "B"})

    def test_nexml_generic(self):
        xml = ("<nexml xmlns=\"http://www.nexml.org/2009\">"
               "<otus id=\"o\"><otu id=\"t1\" label=\"A\"/><otu id=\"t2\" label=\"B\"/></otus>"
               "<trees id=\"ts\"><tree id=\"tr\"><node id=\"n1\" otu=\"t1\"/>"
               "<node id=\"n2\" otu=\"t2\"/><node id=\"n0\"/>"
               "<root id=\"n0\"/><edge id=\"e1\" source=\"n0\" target=\"n1\"/>"
               "<edge id=\"e2\" source=\"n0\" target=\"n2\"/></tree></trees></nexml>")
        t = read(xml)
        self.assertEqual(set(t.tip_names), {"A", "B"})

    def test_mrbayes(self):
        text = ("#NEXUS\nBEGIN TREES;\n TREE t = ((A,B),C);\nEND;\n"
                "BEGIN MRBAYES;\n lset nst=6 rates=gamma;\n mcmc ngen=1000;\nEND;")
        t = read_mrbayes(text)
        self.assertEqual(set(t.tip_names), {"A", "B", "C"})
        self.assertIn("lset", t.get_data("mrbayes"))

    def test_iqtree_raxml(self):
        nwk = "((A:1,B:1)95:0.3,C:1);"
        self.assertEqual(read_iqtree(nwk).tip_names[:2], ["A", "B"])
        self.assertEqual(read_raxml("(A,B);").tip_names, ["A", "B"])


class TestStreaming(unittest.TestCase):
    def setUp(self):
        import matplotlib
        matplotlib.use("Agg")

    def _tree(self):
        return read("(A:1,B:2,C:3);")

    def test_render_save_file(self):
        import matplotlib
        from treeio import render
        d = tempfile.mkdtemp()
        p = os.path.join(d, "t.png")
        ret = render(self._tree(), backend="mpl", layout="rectangular", path=p)
        self.assertEqual(ret, p)
        self.assertTrue(os.path.exists(p))

    def test_tree_save(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "t.png")
        self._tree().save(p, layout="circular")
        self.assertTrue(os.path.exists(p))


if __name__ == "__main__":
    unittest.main()
