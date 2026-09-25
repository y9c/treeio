treeio
======

Why is there no tree-visualization package as powerful as ggtree in R?

Build a tree straight from an alignment (pure Python, no Biopython needed):

.. code-block:: python

   from treeio import build_tree
   tree = build_tree({"human": "ATG", "chimp": "ATG", "mouse": "ATC"}, model="jc")

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
