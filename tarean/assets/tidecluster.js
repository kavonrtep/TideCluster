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
  // Shared by the TRA ideogram rectangles and the cluster-overview
  // scatter circles. Elements opt in via class + a data-title attribute;
  // data-html="1" means data-title is HTML (used by the scatter points
  // for a multi-line tooltip), otherwise it is plain text.
  var TOOLTIP_SEL = 'rect.tc-tra, rect.tc-idx-tra, circle.tc-point';
  function initTraTooltip() {
    var tip = null;
    document.addEventListener('mouseover', function (ev) {
      var r = ev.target.closest(TOOLTIP_SEL);
      if (!r) return;
      var text = r.getAttribute('data-title');
      if (!text) return;
      tip = tip || ensureTooltip();
      if (r.getAttribute('data-html')) tip.innerHTML = text;
      else tip.textContent = text;
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
      var r = ev.target.closest(TOOLTIP_SEL);
      if (!r || !tip) return;
      // Only hide when we truly leave the element; related target check
      // avoids flicker when moving between adjacent shapes.
      if (ev.relatedTarget && ev.relatedTarget.closest &&
          ev.relatedTarget.closest(TOOLTIP_SEL)) return;
      tip.classList.remove('tc-tooltip-visible');
    });
  }

  // ---------- init on DOMContentLoaded ---------------------------------------
  function ready(fn) {
    if (document.readyState !== 'loading') { fn(); }
    else { document.addEventListener('DOMContentLoaded', fn); }
  }

  // ---------- per-array Details child row (DataTables row.child) ------------
  // The arrays table carries only a tiny `data-array-key` per <tr>; the heavy
  // per-peak diagnostics are NOT inlined. They live in one embedded per-TRC
  // JSON block (`#tc-peaks`, with column dictionary + tier map in
  // `#tc-coldict`), and the child-row HTML is BUILT LAZILY in JS on first
  // unfold — so satellite-rich dashboards stay light and nothing is rendered
  // until a row is opened. Self-contained: no fetch(), works from file://.
  var _peaksData = null, _colDict = null, _legendHtml = null;
  function _data() {
    if (_peaksData === null) {
      var el = document.getElementById('tc-peaks');
      try { _peaksData = el ? JSON.parse(el.textContent) : {}; }
      catch (e) { _peaksData = {}; }
      var cd = document.getElementById('tc-coldict');
      try { _colDict = cd ? JSON.parse(cd.textContent) : {cols: [], tiers: {}}; }
      catch (e) { _colDict = {cols: [], tiers: {}}; }
    }
    return _peaksData;
  }
  function _esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }
  function _fmt(v, fmt) {
    if (v === null || v === undefined || v === '') return '';
    if (fmt === 'int') return String(v);
    if (fmt === 'f3') return Number(v).toFixed(3);
    if (fmt === 'f4') return Number(v).toFixed(4);
    if (fmt === 'yn') return v ? '✓' : '';
    return _esc(v);
  }
  function _idstr(v) { return (v === null || v === undefined) ? '' : 'id_med ' + Number(v).toFixed(3); }
  // Shared "column meanings" legend block — built once, revealed from any
  // child row (NOT duplicated per row).
  function _legend() {
    if (_legendHtml !== null) return _legendHtml;
    var cols = (_colDict && _colDict.cols) || [];
    var rows = cols.map(function (c) {
      // c = [key, header, fmt, description]
      return '<dt><code>' + _esc(c[1]) + '</code></dt><dd>' + _esc(c[3]) + '</dd>';
    }).join('');
    _legendHtml = '<details class="tc-col-legend"><summary>ⓘ column meanings</summary>'
      + '<dl class="tc-legend-list tc-legend-cols">' + rows + '</dl></details>';
    return _legendHtml;
  }
  function _buildDetail(key) {
    var d = _data()[key];
    if (!d) return '<div class="tc-dim" style="margin-top:6px">no rescored peaks for this array</div>';
    var cols = (_colDict && _colDict.cols) || [];
    var tiers = (_colDict && _colDict.tiers) || {};
    // ----- header summary -----
    var bits = [];
    if (d.f !== null && d.f !== undefined) {
      var fb = d.fb ? '<sup title="rescore fallback to top-scored kite peak">*</sup>' : '';
      var cr = d.cr ? '<sup title="founder = median of summed-score peak cluster (Pass 5 rescue; rescore returned NA)">‡</sup>' : '';
      bits.push('<strong>Founder:</strong> ' + _esc(d.f) + fb + cr
        + ' <span class="tc-dim">' + _idstr(d.fid) + '</span>');
    } else { bits.push('<strong>Founder:</strong> NA'); }
    if (d.s !== null && d.s !== undefined) {
      bits.push('<strong>Strongest:</strong> ' + _esc(d.s)
        + ' <span class="tc-dim">' + _idstr(d.sid) + '</span>');
    }
    if (d.m > 1) {
      var kraw = (d.mr !== null && d.mr !== undefined)
        ? ' <span class="tc-dim">(k = ' + Number(d.mr).toFixed(3) + ')</span>' : '';
      bits.push('<strong>&times;' + d.m + '</strong>'
        + ' <span class="hor-badge hor-tier-' + _esc(d.ht) + '" title="HOR-order confidence: '
        + _esc(d.ht) + '">' + _esc(d.htl || d.ht) + '</span>' + kraw);
    } else { bits.push('<strong>&times;' + d.m + '</strong>'); }
    if (d.d !== null && d.d !== undefined) {
      var dv = Number(d.d); bits.push('<span class="tc-dim">&Delta;id ' + (dv >= 0 ? '+' : '') + dv.toFixed(2) + '&nbsp;pp</span>');
    }
    if (d.ws) {
      var altxt = (d.al !== null && d.al !== undefined)
        ? ' &nbsp;<span class="tc-dim">(longer alternative &asymp; ' + _esc(d.al) + '&nbsp;bp)</span>'
        : ' <span class="tc-dim">(no credible longer alternative)</span>';
      bits.push('<span class="tc-weakshort" title="short founder kept on weak kite support (founder peak kite score &lt; 0.20)">⚠ weak short founder</span>' + altxt);
    }
    var header = '<div class="tc-details-summary">' + bits.join(' &nbsp;·&nbsp; ') + '</div>';
    // ----- alternative periodicities -----
    var alt = '';
    if (d.alt && d.alt.length) {
      var ab = d.alt.map(function (c) {  // c = [period, score_sum, n_peaks]
        var meta = [];
        if (c[1] !== null && c[1] !== undefined) meta.push('&Sigma;&nbsp;' + Number(c[1]).toFixed(4));
        if (c[2] !== null && c[2] !== undefined) meta.push('n=' + c[2]);
        var m = meta.length ? ' <span class="tc-meta">(' + meta.join(' &middot; ') + ')</span>' : '';
        return _esc(c[0]) + '&nbsp;bp' + m;
      }).join(' &nbsp;&middot;&nbsp; ');
      alt = '<div class="tc-details-alt tc-dim" style="margin-top:2px;font-size:0.9em">Other periodicities: ' + ab + '</div>';
    }
    // ----- SSR scan -----
    var ssr = '';
    if (d.ssr) {  // [motif, dom_pct, tot_pct, top]
      ssr = '<div class="tc-details-ssr">SSR scan: <strong>' + _esc(d.ssr[0]) + '</strong>';
      if (d.ssr[1] !== null && d.ssr[1] !== undefined) ssr += ' ' + Number(d.ssr[1]).toFixed(1) + '% dominant';
      if (d.ssr[2] !== null && d.ssr[2] !== undefined && d.ssr[2] !== d.ssr[1]) ssr += ' &nbsp;·&nbsp; ' + Number(d.ssr[2]).toFixed(1) + '% total';
      if (d.ssr[3] && d.ssr[3] !== d.ssr[0]) ssr += ' &nbsp;·&nbsp; top: ' + _esc(d.ssr[3]);
      ssr += '</div>';
    }
    // ----- peak table -----
    var table = '';
    if (d.pk && d.pk.length) {
      var thead = cols.map(function (c) {
        return '<th title="' + _esc(c[3]) + '">' + _esc(c[1]) + '</th>';
      }).join('');
      var body = d.pk.map(function (row) {
        var tds = cols.map(function (c, i) {
          if (c[0] === 'tier') {
            var t = row[i]; if (!t) return '<td></td>';
            var ts = tiers[t] || ['?', 'rej'];
            return '<td><span class="tc-tier tc-tier-' + _esc(ts[1]) + '">' + ts[0] + '</span></td>';
          }
          return '<td>' + _fmt(row[i], c[2]) + '</td>';
        }).join('');
        return '<tr>' + tds + '</tr>';
      }).join('');
      table = '<div class="tc-details-table-wrap"><table class="tc-details-table"><thead><tr>'
        + thead + '</tr></thead><tbody>' + body + '</tbody></table></div>' + _legend();
    } else {
      table = '<div class="tc-dim" style="margin-top:6px">no rescored peaks for this array</div>';
    }
    return header + alt + ssr + table;
  }

  function initDetailsExpand() {
    document.addEventListener('click', function (ev) {
      var ctrl = ev.target.closest('.tc-details-control');
      if (!ctrl || ctrl.classList.contains('tc-details-control-static')) return;
      ev.preventDefault();
      var tr = ctrl.closest('tr');
      if (!tr) return;
      var key = tr.getAttribute('data-array-key');
      if (!key) return;
      var html = _buildDetail(key);
      // Try the DataTables path first.
      if (window.jQuery) {
        var $tr = jQuery(tr);
        var $table = $tr.closest('table.tc-datatable');
        if ($table.length && jQuery.fn.DataTable.isDataTable($table[0])) {
          var dt = $table.DataTable();
          var row = dt.row($tr);
          if (row.child.isShown()) {
            row.child.hide();
            tr.classList.remove('tc-details-open');
            ctrl.classList.remove('open');
          } else {
            row.child(html, 'tc-details-row').show();
            tr.classList.add('tc-details-open');
            ctrl.classList.add('open');
          }
          return;
        }
      }
      // Fallback: toggle a sibling <tr>.
      var next = tr.nextElementSibling;
      if (next && next.classList.contains('tc-details-row')) {
        next.parentNode.removeChild(next);
        tr.classList.remove('tc-details-open');
        ctrl.classList.remove('open');
      } else {
        var nrow = document.createElement('tr');
        nrow.className = 'tc-details-row';
        var nc = document.createElement('td');
        nc.colSpan = tr.children.length;
        nc.innerHTML = html;
        nrow.appendChild(nc);
        tr.parentNode.insertBefore(nrow, tr.nextSibling);
        tr.classList.add('tc-details-open');
        ctrl.classList.add('open');
      }
    });
  }

  // ---------- index page: TRC distribution selector ---------------------
  // The cross-TRC ideogram on the index page lets the reader pick one
  // TRC from the left-side list to highlight just its arrays on the
  // chromosome bars. Default state shows every TRA in muted grey.
  //
  // - Click a TRC item        → make it active (others fade via CSS).
  // - Click the active item   → clear (back to "show all").
  // - "Show all" button       → same as clicking the active item.
  // - Search input            → live-filter the TRC list by substring.
  // - The TRC name itself is a <a href> link to the TRC dashboard; we
  //   keep that link working — clicks on the link area navigate, clicks
  //   on the surrounding row toggle the highlight.
  function initIndexDistribution() {
    var dist = document.querySelector('.tc-idx-dist');
    if (!dist) return;
    var svg = dist.querySelector('.tc-idx-svg');
    var list = dist.querySelector('.tc-idx-trc-list');
    var search = dist.querySelector('.tc-idx-trc-search');
    var reset = dist.querySelector('.tc-idx-trc-reset');
    if (!svg || !list) return;

    function applyActive(trcId) {
      svg.setAttribute('data-active-trc', trcId || '');
      // Toggle `.tc-idx-active` on every TRA rect so the per-rect --sf
      // CSS variable can take effect on the matching ones only.
      var rects = svg.querySelectorAll('rect.tc-idx-tra');
      rects.forEach(function (r) {
        if (trcId && r.getAttribute('data-trc') === trcId) {
          r.classList.add('tc-idx-active');
        } else {
          r.classList.remove('tc-idx-active');
        }
      });
      // Highlight the corresponding sidebar row.
      list.querySelectorAll('.tc-idx-trc-item').forEach(function (li) {
        li.classList.toggle('tc-idx-trc-active',
          !!trcId && li.getAttribute('data-trc') === trcId);
      });
    }

    list.addEventListener('click', function (ev) {
      // Let real link clicks navigate.
      if (ev.target.closest('a.tc-idx-trc-link')) return;
      var item = ev.target.closest('.tc-idx-trc-item');
      if (!item) return;
      var trcId = item.getAttribute('data-trc');
      var current = svg.getAttribute('data-active-trc') || '';
      applyActive(current === trcId ? '' : trcId);
    });

    if (reset) {
      reset.addEventListener('click', function () { applyActive(''); });
    }
    if (search) {
      search.addEventListener('input', function () {
        var q = search.value.trim().toLowerCase();
        list.querySelectorAll('.tc-idx-trc-item').forEach(function (li) {
          var id = (li.getAttribute('data-trc') || '').toLowerCase();
          li.classList.toggle('tc-idx-trc-hidden', !!q && id.indexOf(q) === -1);
        });
      });
    }
  }

  ready(function () {
    initDataTables();
    initConsensusTools();
    initTrcJumper();
    initTraTooltip();
    initDetailsExpand();
    initIndexDistribution();
  });
})();
