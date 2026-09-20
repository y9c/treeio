import unittest, tempfile, os
from treeio import read, render
from treeio.show import tree_to_ascii


class TestUnifiedRender(unittest.TestCase):
    def setUp(self):
        import matplotlib
        matplotlib.use("Agg")
        self.t = read("(A:1,B:2,C:3);")

    def test_mpl_default_backend(self):
        import matplotlib
        ax = render(self.t, layout="rectangular")
        self.assertIsNotNone(ax)
        self.assertGreaterEqual(len(ax.collections), 1)
        # tree.render sugar
        self.assertIsNotNone(self.t.render(backend="mpl", layout="circular"))

    def test_ascii_backend(self):
        out = render(self.t, backend="ascii")
        self.assertIn("A", out)
        self.assertEqual(out, tree_to_ascii(self.t, False, True))
        self.assertEqual(self.t.render(backend="ascii"), out)

    def test_mpl_save_path(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "t.png")
        ret = render(self.t, backend="mpl", layout="rectangular", path=p)
        self.assertEqual(ret, p)
        self.assertTrue(os.path.exists(p))

    def test_bad_backend(self):
        with self.assertRaises(ValueError):
            render(self.t, backend="svg")


if __name__ == "__main__":
    unittest.main()
