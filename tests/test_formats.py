import unittest
from treeio import (read, write, read_phylip, write_phylip, read_jplace,
                    read_newick, write_newick, detect_format)


class TestFormats(unittest.TestCase):
    def test_phylip_roundtrip(self):
        t = read_phylip("3 1\n((A:0.1,B:0.2):0.3,C:0.4);")
        self.assertEqual(t.tip_names, ["A", "B", "C"])
        s = write_phylip(t)
        self.assertTrue(s.startswith("3 1\n"))
        self.assertEqual(read_phylip(s), t)

    def test_phylip_no_header(self):
        t = read_phylip("(A,B);")
        self.assertEqual(t.tip_names, ["A", "B"])

    def test_jplace(self):
        doc = '{"tree": "(A,B);", "fields": ["like_weight_ratio"], "placements": [{"p": [[1,0.9]], "n": ["q1"]}]}'
        t = read_jplace(doc)
        self.assertEqual(t.tip_names, ["A", "B"])
        self.assertEqual(t.placements["fields"], ["like_weight_ratio"])
        self.assertEqual(len(t.placements["placements"]), 1)

    def test_nhx_annotations(self):
        t = read_newick("(A:1[&&NHX:S=human,D=Y],B:1[&&NHX:S=chimp,D=Y])[&&NHX:S=Hominidae];",
                        annotations=True)
        tips = {x.name: x for x in t.get_tips()}
        self.assertEqual(tips["A"].get_data("S"), "human")
        self.assertEqual(tips["B"].get_data("S"), "chimp")
        # write back with annotations preserved
        s = write_newick(t, annotations=["S", "D"])
        self.assertIn("[&", s)

    def test_detect_phylip_and_jplace(self):
        self.assertEqual(detect_format("5 1\n((A,B));"), "phylip")
        self.assertEqual(detect_format('{"tree":"(A,B);","placements":[]}'), "jplace")

    def test_generic_read_phylip(self):
        t = read("4 1\n((A:1,B:0.2):0.1,(C:0.4,D:0.3));")
        self.assertEqual(len(t.tip_names), 4)


if __name__ == "__main__":
    unittest.main()
