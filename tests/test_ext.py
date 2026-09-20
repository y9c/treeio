import unittest
from treeio import (read, write, register_format, register_backend, register_layout,
                    tree_coords, render, layouts)
from treeio import read_newick, write_newick


class TestExtensibility(unittest.TestCase):
    def test_register_format(self):
        # a trivial "myfmt" = newick with a marker
        def myread(text):
            return [read_newick(text.split("\n", 1)[1])]
        def mywrite(tree):
            return "MYFMT\n" + write_newick(tree)
        register_format("myfmt", read=myread, write=mywrite, extensions=[".myfmt"])
        t = read("MYFMT\n(A:1,B:2);", format="myfmt")
        self.assertEqual(t.tip_names, ["A", "B"])
        s = write(t, format="myfmt")
        self.assertTrue(s.startswith("MYFMT\n"))
        # extension routing
        import tempfile, os
        p = os.path.join(tempfile.mkdtemp(), "x.myfmt")
        with open(p, "w") as fh:
            fh.write("MYFMT\n(A,B);")
        self.assertEqual(read(p).tip_names, ["A", "B"])

    def test_register_layout(self):
        def spiral(tree):
            base = tree_coords(tree, layout="rectangular")
            import math
            return {n: (x * math.cos(y * 0.5), x * math.sin(y * 0.5)) for n, (x, y) in base.items()}
        register_layout("spiral", spiral)
        self.assertIn("spiral", layouts())
        t = read_newick("(A,B,C);")
        coords = tree_coords(t, layout="spiral")
        self.assertTrue(all(n in coords for n in t.traverse()))

    def test_register_backend(self):
        def mybackend(tree, layout="rectangular", **kw):
            return f"<custom {layout} {len(tree.tip_names)}/>"
        register_backend("customxml", mybackend)
        t = read_newick("(A,B,C);")
        out = render(t, backend="customxml", layout="circular")
        self.assertEqual(out, "<custom circular 3/>")


if __name__ == "__main__":
    unittest.main()
