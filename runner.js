/*
 * treeio live-runner (Pyodide / WASM, no backend)
 *
 * Loads CPython-in-WASM (Pyodide), installs numpy + matplotlib + a locally
 * hosted treeio wheel, and exposes `run(code)` which executes user Python and
 * returns the resulting figure as a PNG data URL.  Everything happens in the
 * browser; GitHub Pages static hosting is sufficient.
 *
 * Public API
 * ----------
 *   window.treeioLive = {
 *     run: async (pythonCode) => ({ image: dataURL|null, stdout: '' }),
 *     status: () => 'idle' | 'loading' | 'ready' | 'error',
 *     ready: Promise,
 *     onReady: (cb) => void,
 *     warmup: () => void,
 *   }
 */
(function () {
  'use strict';

  var PYODIDE_BASE =
    window.TREEIO_PYODIDE_URL ||
    'https://cdn.jsdelivr.net/pyodide/v0.26.2/full/';

  function wheelUrl() {
    return window.TREEIO_WHEEL_URL ||
      './packages/treeio-0.0.0.dev10-py3-none-any.whl';
  }

  var status = 'idle';
  var kernel = null;
  var bootPromise = null;
  var callbacks = [];

  function setStatus(s) {
    status = s;
    callbacks.forEach(function (cb) { try { cb(s); } catch (e) {} });
  }

  async function boot() {
    if (kernel) return kernel;
    if (bootPromise) return bootPromise;
    bootPromise = (async () => {
      setStatus('loading');
      const { loadPyodide } = await import(PYODIDE_BASE + 'pyodide.mjs');
      const py = await loadPyodide({ indexURL: PYODIDE_BASE });

      // Prebuilt binary packages -> loadPackage (lower overhead than micropip).
      await py.loadPackage(['numpy', 'matplotlib', 'micropip']);

      const micropip = py.pyimport('micropip');
      await micropip.install(wheelUrl());

      // One-time prelude: Agg renderer + pre-import treeio so each run() is
      // only user code.  Kept in a Python module so imports are cached.
      await py.runPythonAsync(`
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import treeio
from treeio import *
plt.ioff()

def __treeio_render(user_code):
    _ns = {'treeio': treeio, 'plt': plt, 'Tree': Tree, 'draw_tree': draw_tree,
           'highlight_clade': highlight_clade, 'add_ring': add_ring,
           'add_rings': add_rings, 'gheatmap': gheatmap, 'attach': attach,
           'read_newick': read_newick, 'named_palette': named_palette,
           'tree_coords': tree_coords, 'add_boxed_labels': add_boxed_labels,
           'box_label': box_label}
    exec(compile(user_code, '<live>', 'exec'), _ns)
    figs = [plt.figure(n) for n in plt.get_fignums()]
    if not figs:
        result = {'image': None, 'stdout': ''}
    else:
        f = figs[-1]
        buf = BytesIO()
        f.savefig(buf, format='png', dpi=110, bbox_inches='tight')
        data = base64.b64encode(buf.getvalue()).decode('ascii')
        plt.close('all')
        result = {'image': 'data:image/png;base64,' + data, 'stdout': ''}
    import json
    return json.dumps(result)
`);
      kernel = py;
      setStatus('ready');
      return py;
    })();
    return bootPromise;
  }

  async function run(code, opts) {
    const py = await boot();
    // __treeio_render returns a JSON *string*; parse it to a JS object.
    const raw = await py.runPythonAsync(
      '__treeio_render(' + JSON.stringify(code) + ')'
    );
    if (typeof raw === 'string') {
      return JSON.parse(raw);
    }
    // conservatively coerce a dict that pyodide already converted to a JS object
    return { image: (raw && raw.image) || null, stdout: (raw && raw.stdout) || '' };
  }

  window.treeioLive = {
    run: run,
    status: function () { return status; },
    ready: boot(),
    onReady: function (cb) {
      callbacks.push(cb);
      boot().then(function () { cb('ready'); }).catch(function (e) { cb('error'); });
    },
    warmup: function () { boot().catch(function () {}); },
  };

  boot().catch(function (e) {
    setStatus('error');
    if (window.console) console.error('treeio live boot failed:', e);
  });
})();
