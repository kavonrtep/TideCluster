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

  // ---------- init on DOMContentLoaded ---------------------------------------
  function ready(fn) {
    if (document.readyState !== 'loading') { fn(); }
    else { document.addEventListener('DOMContentLoaded', fn); }
  }

  ready(function () {
    initDataTables();
    initConsensusTools();
    initTrcJumper();
  });
})();
