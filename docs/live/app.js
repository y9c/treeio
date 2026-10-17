/*
 * treeio live gallery app
 *
 * - First paint is instant: each example shadows a static PNG preview.
 * - The Pyodide kernel is warmed up in the background (toast shows progress).
 * - Once ready, the user can edit the snippet and re-run it live (no backend).
 */
(function () {
  'use strict';

  var examples = window.TREEIO_EXAMPLES || [];
  var current = null;
  var renderedLive = false;

  // runner.js is an ES module (loadPyodide is ESM-only); import it lazily.
  // It assigns window.treeioLive once loaded, so poll briefly for it.
  var _liveResolve = null;
  var _liveReady = new Promise(function (res) { _liveResolve = res; });
  function ensureRunner() {
    if (window.treeioLive) return Promise.resolve(window.treeioLive);
    return _liveReady.then(function () { return window.treeioLive; });
  }
  (function loadRunner() {
    if (!window.treeioLive) {
      import('./runner.js').then(function () {
        _liveResolve(window.treeioLive);
      }).catch(function (e) {
        if (window.console) console.error('failed to import runner.js', e);
      });
    }
  })();

  var listEl = document.getElementById('list');
  var titleEl = document.getElementById('title');
  var descEl = document.getElementById('desc');
  var codeEl = document.getElementById('code');
  var runBtn = document.getElementById('runBtn');
  var canvas = document.getElementById('canvas');
  var placeholder = document.getElementById('placeholder');
  var liveEl = document.getElementById('live');
  var toastEl = document.getElementById('toast');

  // ---------- toast (bottom-right progress) -------------------------------
  function toast(msg, state, ttl) {
    var el = document.createElement('div');
    el.className = 'toast' + (state ? ' ' + state : '');
    el.innerHTML = '<span class="spinner"></span><span class="msg"></span>';
    el.querySelector('.msg').textContent = msg;
    toastEl.appendChild(el);
    if (ttl) {
      setTimeout(function () {
        el.classList.add('hide');
        setTimeout(function () { el.remove(); }, 320);
      }, ttl);
    }
    return el;
  }

  function setLive(state, label) {
    liveEl.className = 'live ' + state;
    liveEl.textContent = label;
  }

  // ---------- build the sidebar -------------------------------------------
  examples.forEach(function (ex, i) {
    var btn = document.createElement('button');
    btn.className = 'ex'; btn.dataset.id = ex.id;
    btn.innerHTML = '<span>' + ex.title + '</span><span class="d">' + (ex.desc || '') + '</span>';
    btn.addEventListener('click', function () { select(i); });
    listEl.appendChild(btn);
  });

  // ---------- render static preview (instant) -----------------------------
  function showStatic(ex) {
    canvas.innerHTML = '';
    var img = document.createElement('img');
    img.src = ex.preview;
    img.alt = ex.title;
    canvas.appendChild(img);
    renderedLive = false;
  }

  // ---------- select an example -------------------------------------------
  function select(i) {
    current = examples[i];
    Array.prototype.forEach.call(listEl.querySelectorAll('.ex'), function (b) {
      b.classList.toggle('active', b.dataset.id === current.id);
    });
    titleEl.textContent = current.title;
    descEl.textContent = current.desc || '';
    codeEl.value = (current.code || '').replace(/^\n/, '');
    showStatic(current);
    // status pill reflects the async engine state; default is 'loading'.
    if (window.treeioLive && window.treeioLive.status) {
      var st = window.treeioLive.status();
      if (st === 'ready') setLive('ready', 'engine ready');
      else if (st === 'error') setLive('error', 'engine error');
      else setLive('loading', 'engine loading…');
    } else {
      setLive('loading', 'engine loading…');
    }
  }

  // ---------- run user code via WASM ---------------------------------------
  async function run() {
    if (!current) return;
    // indicate a live render is starting
    placeholder.textContent = 'Rendering…';
    runBtn.disabled = true;
    try {
      var r = await window.treeioLive.run(codeEl.value);
      if (r && r.image) {
        canvas.innerHTML = '';
        var img = document.createElement('img');
        img.src = r.image;
        canvas.appendChild(img);
        renderedLive = true;
        toast('Figure rendered live', 'done', 2600);
      } else {
        placeholder.textContent = 'No figure was produced (did you create a plot?)';
        canvas.innerHTML = '';
        canvas.appendChild(placeholder);
      }
    } catch (e) {
      placeholder.textContent = 'Error: ' + e.message;
      canvas.innerHTML = '';
      canvas.appendChild(placeholder);
      toast('Run error', 'error', 4000);
      if (window.console) console.error(e);
    } finally {
      runBtn.disabled = window.treeioLive.status() !== 'ready';
    }
  }

  runBtn.addEventListener('click', run);
  // Ctrl/Cmd+Enter runs too
  codeEl.addEventListener('keydown', function (e) {
    if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) {
      e.preventDefault(); run();
    }
  });

  // ---------- select the first example immediately ------------------------
  select(0);

  // ---------- warm the kernel in the background ---------------------------
  var loadingToast = toast('Loading the in-browser engine (Pyodide + matplotlib)…', '');
  setLive('loading', 'engine loading…');

  ensureRunner().then(function (live) {
    live.onReady(function () {
      loadingToast.classList.add('done');
      loadingToast.querySelector('.spinner').style.borderTopColor = '';
      loadingToast.querySelector('.msg').textContent = 'Engine ready — you can edit & run now.';
      setTimeout(function () { loadingToast.classList.add('hide'); }, 3500);
      setLive('ready', 'engine ready');
      runBtn.disabled = false;
    });
    live.ready.catch(function (e) {
      loadingToast.classList.add('error');
      loadingToast.querySelector('.msg').textContent = 'Engine failed to load: ' + (e && e.message);
      setLive('error', 'engine error');
    });
  }).catch(function (e) {
    setLive('error', 'engine error');
  });
})();
