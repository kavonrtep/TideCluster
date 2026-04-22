/* TideCluster report v2 — shared client-side behaviour.
   Runs after jQuery + DataTables have loaded. Requires nothing else. */

(function () {
  'use strict';

  // ---------- DataTables init -------------------------------------------------
  // Every <table class="tc-table tc-datatable"> becomes a DataTable with the
  // TideCluster defaults. Opt out by omitting the tc-datatable class.
  function initDataTables() {
    if (typeof jQuery === 'undefined' || !jQuery.fn.DataTable) return;
    jQuery('table.tc-datatable').each(function () {
      var $t = jQuery(this);
      if ($t.hasClass('dataTable')) return; // already wired
      var opts = {
        pageLength: parseInt($t.data('page-length') || '25', 10),
        lengthMenu: [10, 25, 50, 100, 250, 1000],
        order: [],
        deferRender: true,
        autoWidth: false,
        stripeClasses: ['odd', 'even']
      };
      // Respect an explicit column ordering via data-order="[[col,dir]]".
      var dataOrder = $t.data('order');
      if (dataOrder) { try { opts.order = JSON.parse(dataOrder); } catch (e) {} }
      $t.DataTable(opts);
    });
  }

  // ---------- Consensus copy + download button --------------------------------
  // Buttons are rendered by the Python side as:
  //   <button class="tc-btn" data-copy-target="#some-pre">Copy</button>
  //   <a      class="tc-btn" data-download="TRC_1.fasta"
  //                          data-download-source="#some-pre">Download FASTA</a>
  function initConsensusTools() {
    document.addEventListener('click', function (ev) {
      var copyBtn = ev.target.closest('[data-copy-target]');
      if (copyBtn) {
        var target = document.querySelector(copyBtn.dataset.copyTarget);
        if (!target) return;
        var text = target.textContent;
        navigator.clipboard.writeText(text).then(function () {
          var original = copyBtn.textContent;
          copyBtn.textContent = 'copied';
          setTimeout(function () { copyBtn.textContent = original; }, 1200);
        });
      }
      var dl = ev.target.closest('[data-download-source]');
      if (dl && dl.tagName === 'A' && !dl.href) {
        ev.preventDefault();
        var src = document.querySelector(dl.dataset.downloadSource);
        if (!src) return;
        var blob = new Blob([src.textContent + '\n'], { type: 'text/plain' });
        var url = URL.createObjectURL(blob);
        dl.href = url;
        dl.setAttribute('download', dl.dataset.download || 'consensus.fasta');
        dl.click();
        setTimeout(function () { URL.revokeObjectURL(url); dl.removeAttribute('href'); }, 500);
      }
    });
  }

  // ---------- TRC jump dropdown ----------------------------------------------
  // A <select class="tc-trc-jumper" data-base="trc/"> navigates on change.
  function initTrcJumper() {
    document.addEventListener('change', function (ev) {
      var sel = ev.target.closest('select.tc-trc-jumper');
      if (!sel) return;
      var target = sel.value;
      if (!target) return;
      var base = sel.dataset.base || '';
      window.location.href = base + target + '.html';
    });
  }

  // ---------- TRA rectangle tooltip ------------------------------------------
  // Browsers render <title> tooltips on SVG elements with a ~1 s delay
  // that made users think hovering "didn't work". We show an instant
  // floating tooltip instead, populated from the element's
  // data-title attribute.
  function ensureTooltip() {
    var t = document.getElementById('tc-tooltip');
    if (t) return t;
    t = document.createElement('div');
    t.id = 'tc-tooltip';
    t.className = 'tc-tooltip';
    t.setAttribute('role', 'tooltip');
    document.body.appendChild(t);
    return t;
  }
  function initTraTooltip() {
    var tip = null;
    document.addEventListener('mouseover', function (ev) {
      var r = ev.target.closest('rect.tc-tra');
      if (!r) return;
      var text = r.getAttribute('data-title');
      if (!text) return;
      tip = tip || ensureTooltip();
      tip.textContent = text;
      tip.classList.add('tc-tooltip-visible');
    });
    document.addEventListener('mousemove', function (ev) {
      if (!tip || !tip.classList.contains('tc-tooltip-visible')) return;
      var x = ev.pageX + 12, y = ev.pageY - 8;
      // Nudge into viewport if the tooltip would overflow on the right.
      var vw = document.documentElement.clientWidth;
      if (x + tip.offsetWidth > vw - 4) x = ev.pageX - tip.offsetWidth - 12;
      tip.style.left = x + 'px';
      tip.style.top  = y + 'px';
    });
    document.addEventListener('mouseout', function (ev) {
      var r = ev.target.closest('rect.tc-tra');
      if (!r || !tip) return;
      // Only hide when we truly leave the rect; related target check
      // avoids flicker when moving between adjacent rectangles.
      if (ev.relatedTarget && ev.relatedTarget.closest &&
          ev.relatedTarget.closest('rect.tc-tra')) return;
      tip.classList.remove('tc-tooltip-visible');
    });
  }

  // ---------- init on DOMContentLoaded ---------------------------------------
  function ready(fn) {
    if (document.readyState !== 'loading') { fn(); }
    else { document.addEventListener('DOMContentLoaded', fn); }
  }

  ready(function () {
    initDataTables();
    initConsensusTools();
    initTrcJumper();
    initTraTooltip();
  });
})();
