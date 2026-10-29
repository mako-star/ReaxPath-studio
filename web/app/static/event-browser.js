"use strict";

// Event records come from the analysis index; every atom label is a simulation ID.
window.EventBrowser = (() => {
  const byId = (id) => document.getElementById(id);
  const escape = (value) => String(value ?? "").replace(/[&<>"']/g, (c) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  }[c]));
  const number = (value) => value === null || value === undefined ? "未提供" : Number(value).toLocaleString("zh-CN", { maximumFractionDigits: 8 });
  const bondKey = (pair) => [...pair].sort((a, b) => a - b).join("-");
  const model = { jobId: null, offset: 0, limit: 50, total: 0, selected: null, detail: null, atom: null, request: 0, detailRequest: 0 };

  async function get(url) {
    const response = await fetch(url);
    if (!response.ok) {
      let message = `读取失败（${response.status}）`;
      try {
        const data = await response.json();
        message = typeof data.detail === "string" ? data.detail : message;
      } catch (_) { /* The HTTP status remains useful for a non-JSON error. */ }
      throw new Error(message);
    }
    return response.json();
  }

  function error(message) {
    byId("event-error").textContent = message || "";
    byId("event-error").classList.toggle("hidden", !message);
  }

  function clearDetail(message = "选择事件，查看反应前后的连接图。") {
    model.detailRequest += 1;
    model.selected = null;
    model.detail = null;
    model.atom = null;
    byId("event-detail").setAttribute("aria-busy", "false");
    byId("event-detail").innerHTML = `<p class="detail-placeholder">${escape(message)}</p>`;
  }

  function reset() {
    model.request += 1;
    model.jobId = null;
    model.offset = 0;
    model.total = 0;
    byId("event-filters").reset();
    byId("events-table").replaceChildren();
    byId("events-count-label").textContent = "";
    byId("event-page-label").textContent = "";
    byId("event-prev").disabled = true;
    byId("event-next").disabled = true;
    clearDetail();
    error("");
  }

  function params() {
    const query = new URLSearchParams({ offset: model.offset, limit: model.limit, status: "all" });
    const fields = { q: "event-query", formula: "event-formula", start_step: "event-start", end_step: "event-end" };
    for (const [key, id] of Object.entries(fields)) {
      const value = byId(id).value.trim();
      if (value) query.set(key, value);
    }
    if (query.has("start_step") && query.has("end_step") && Number(query.get("start_step")) > Number(query.get("end_step"))) {
      throw new Error("起始 step 不能大于结束 step");
    }
    return query;
  }

  async function load(jobId = model.jobId) {
    if (!jobId) return;
    if (model.jobId !== jobId) {
      const url = new URL(window.location.href);
      const eventId = url.searchParams.get("event");
      if (url.searchParams.get("job") === jobId && /^E\d+$/.test(eventId || "")) byId("event-query").value = eventId;
    }
    model.jobId = jobId;
    const request = ++model.request;
    error("");
    byId("events-count-label").textContent = "读取事件…";
    byId("event-prev").disabled = true;
    byId("event-next").disabled = true;
    try {
      const query = params();
      const payload = await get(`/api/jobs/${encodeURIComponent(jobId)}/events?${query}`);
      if (request !== model.request) return;
      model.total = payload.total;
      const rows = payload.rows || [];
      const first = rows.length ? model.offset + 1 : 0;
      byId("events-count-label").textContent = `共 ${number(payload.total)} 条 · 当前 ${number(first)}–${number(model.offset + rows.length)}`;
      byId("event-page-label").textContent = payload.total ? `${Math.floor(model.offset / model.limit) + 1} / ${Math.ceil(payload.total / model.limit)} 页` : "无匹配事件";
      byId("event-prev").disabled = model.offset === 0;
      byId("event-next").disabled = !payload.has_more;
      renderRows(rows);
      if (!rows.length) clearDetail("没有符合条件的候选事件。可清空筛选查看全部原始事件。");
      else if (!rows.some((row) => row.event_id === model.selected)) await select(rows[0].event_id);
    } catch (cause) {
      if (request !== model.request) return;
      error(cause.message);
      byId("events-count-label").textContent = "事件读取失败";
      byId("event-page-label").textContent = "";
      byId("events-table").replaceChildren();
      clearDetail("请检查筛选条件或重新分析任务。");
    }
  }

  function renderRows(rows) {
    byId("events-table").innerHTML = `<thead><tr><th scope="col">事件 / 帧间隔</th><th scope="col">候选反应式</th></tr></thead><tbody>${rows.map((row) => `
      <tr data-event-row="${escape(row.event_id)}" class="${row.event_id === model.selected ? "selected-event" : ""}">
        <td><button type="button" class="event-select" data-event-id="${escape(row.event_id)}" aria-current="${row.event_id === model.selected}" aria-label="查看事件 ${escape(row.event_id)}">${escape(row.event_id)}</button><span class="table-subline">${number(row.previous_timestep)} → ${number(row.timestep)} step</span></td>
        <td class="reaction-cell">${escape(row.reaction)}<span class="table-subline">${number(row.time_start_ps)} → ${number(row.time_end_ps)} ps</span></td>
      </tr>`).join("") || '<tr><td colspan="2" class="empty-events">没有符合条件的事件</td></tr>'}</tbody>`;
  }

  async function select(eventId) {
    model.selected = eventId;
    model.detail = null;
    model.atom = null;
    const request = ++model.detailRequest;
    const jobId = model.jobId;
    byId("events-table").querySelectorAll("[data-event-row]").forEach((row) => {
      const selected = row.dataset.eventRow === eventId;
      row.classList.toggle("selected-event", selected);
      row.querySelector(".event-select").setAttribute("aria-current", String(selected));
    });
    byId("event-detail").setAttribute("aria-busy", "true");
    byId("event-detail").innerHTML = `<p class="detail-placeholder">读取事件 ${escape(eventId)} 的证据…</p>`;
    try {
      const event = await get(`/api/jobs/${encodeURIComponent(jobId)}/events/${encodeURIComponent(eventId)}`);
      if (request !== model.detailRequest || jobId !== model.jobId) return;
      model.detail = event;
      model.atom = null;
      renderDetail(event);
    } catch (cause) {
      if (request !== model.detailRequest) return;
      byId("event-detail").innerHTML = `<p class="inline-error" role="alert">${escape(cause.message)}</p>`;
    } finally {
      if (request === model.detailRequest && jobId === model.jobId) byId("event-detail").setAttribute("aria-busy", "false");
    }
  }

  function moleculeCard(molecule, side, index) {
    const names = { ok: "标识已生成", partial: "部分标识可用", failed: "标识未生成", unavailable: "无标识后端", disabled: "标识已关闭" };
    const status = molecule.identity_status || "unavailable";
    return `<article class="molecule-card">
      <div class="molecule-heading"><strong>${escape(molecule.structure_id || molecule.formula)}</strong><span>${escape(molecule.formula)}</span></div>
      <div class="molecule-graph" data-molecule-side="${side}" data-molecule-index="${index}"></div>
      <span class="molecule-caption">${molecule.atoms.length} 个原子 · ${molecule.connections.length} 条连接</span>
      <details class="molecule-identifiers"><summary>${escape(names[status] || status)}</summary>
        <dl><dt>实例</dt><dd>${escape(molecule.instance_id)}</dd><dt>SMILES</dt><dd>${escape(molecule.smiles || "未生成")}</dd><dt>InChI</dt><dd>${escape(molecule.inchi || "未生成")}</dd><dt>说明</dt><dd>${escape(molecule.message || "无附加说明")}</dd></dl>
      </details>
    </article>`;
  }

  function renderDetail(event) {
    const conserved = event.conservation?.atom_ids === true && event.conservation?.elements === true;
    const sideCards = (side) => {
      const molecules = event[side];
      const visible = molecules.slice(0, 24);
      return visible.map((m, i) => moleculeCard(m, side, i)).join("") +
        (molecules.length > visible.length ? `<p class="topology-note">当前展示 ${visible.length} / ${number(molecules.length)} 个分子实例；完整参与分子见事件 JSON。</p>` : "");
    };
    byId("event-detail").innerHTML = `
      <div class="detail-header"><div><span class="detail-kicker">候选事件 ${escape(event.event_id)}</span><h3>${escape(event.reaction.length > 2000 ? event.reaction.slice(0, 2000) + "…（完整反应见 JSON）" : event.reaction)}</h3></div><a class="small-button" href="/api/jobs/${encodeURIComponent(model.jobId)}/events/${encodeURIComponent(event.event_id)}/export" download>导出此事件 JSON</a></div>
      <div class="detail-meta"><span>${number(event.previous_timestep)} → ${number(event.timestep)} step</span><span>${number(event.time_start_ps)} → ${number(event.time_end_ps)} ps</span><span class="${conserved ? "conserved" : "inline-error"}">${conserved ? "原子与元素守恒已核对" : "守恒尚未确认"}</span></div>
      <p class="topology-note">二维连接拓扑示意；连线不表示已确认的化学键型，位置不代表真实坐标。点击原子编号核对去向。</p>
      <div class="bond-legend"><span class="bond-sample broken">消失的连接 · 虚线</span><span class="bond-sample formed">新增的连接 · 实线</span></div>
      <div class="reaction-comparison">
        <section aria-label="反应前分子"><h4>反应前 · ${number(event.previous_timestep)} step</h4>${sideCards("reactants")}</section>
        <section aria-label="反应后分子"><h4>反应后 · ${number(event.timestep)} step</h4>${sideCards("products")}</section>
      </div>
      <div class="connection-changes"><div><strong>消失的连接</strong><span>${escape(bondSummary(event.broken_bonds))}</span></div><div><strong>新增的连接</strong><span>${escape(bondSummary(event.formed_bonds))}</span></div></div>
      <p id="atom-destination" class="atom-destination" role="status">选择图中的原子，或下方原子去向表中的编号。</p>
      <details class="atom-mapping"><summary>原子去向表（${number(event.involved_atom_ids.length)} 个原子）</summary><div id="atom-map-table" class="table-scroll"></div></details>`;
    paintGraphs(event);
    const rows = event.involved_atom_ids.slice(0, 250).map((atomId) => {
      const before = event.reactants.find((m) => m.atoms.some((a) => a.atom_id === atomId));
      const after = event.products.find((m) => m.atoms.some((a) => a.atom_id === atomId));
      const atom = before?.atoms.find((a) => a.atom_id === atomId);
      return `<tr><td><button type="button" class="atom-table-button" data-select-atom="${atomId}">${escape(atom?.element)} ${atomId}</button></td><td>${escape(before?.structure_id || before?.formula)}</td><td>${escape(after?.structure_id || after?.formula)}</td></tr>`;
    }).join("");
    byId("atom-map-table").innerHTML = `<table><thead><tr><th>原子 ID</th><th>反应前结构</th><th>反应后结构</th></tr></thead><tbody>${rows}</tbody></table>${event.involved_atom_ids.length > 250 ? '<p class="topology-note">页面列出前 250 个原子；完整原子映射可导出此事件 JSON。</p>' : ""}`;
  }

  function bondSummary(edges) {
    return (edges || []).slice(0, 80).map(bondKey).join("；") +
      ((edges || []).length > 80 ? `；…共 ${number(edges.length)} 条，完整列表见 JSON` : "") || "无";
  }

  // Deterministic graph layout, for inspection only. No geometry is inferred.
  function layout(atoms, edges) {
    const ordered = [...atoms].sort((a, b) => a.atom_id - b.atom_id);
    const points = new Map(ordered.map((a, i) => [a.atom_id, { x: Math.cos(i * 2 * Math.PI / ordered.length) * 100, y: Math.sin(i * 2 * Math.PI / ordered.length) * 100 }]));
    if (ordered.length === 1) return new Map([[ordered[0].atom_id, { x: 0, y: 0 }]]);
    for (let iteration = 0; iteration < 95; iteration += 1) {
      const delta = new Map(ordered.map((a) => [a.atom_id, { x: 0, y: 0 }]));
      for (let i = 0; i < ordered.length; i += 1) {
        for (let j = i + 1; j < ordered.length; j += 1) {
          const a = ordered[i].atom_id, b = ordered[j].atom_id;
          const dx = points.get(a).x - points.get(b).x, dy = points.get(a).y - points.get(b).y;
          const distance = Math.max(1, Math.hypot(dx, dy));
          const force = 2000 / (distance * distance);
          delta.get(a).x += dx / distance * force;
          delta.get(a).y += dy / distance * force;
          delta.get(b).x -= dx / distance * force;
          delta.get(b).y -= dy / distance * force;
        }
      }
      for (const [a, b] of edges) {
        if (!points.has(a) || !points.has(b)) continue;
        const dx = points.get(b).x - points.get(a).x, dy = points.get(b).y - points.get(a).y;
        const distance = Math.max(1, Math.hypot(dx, dy));
        const force = (distance - 65) * 0.035;
        delta.get(a).x += dx / distance * force;
        delta.get(a).y += dy / distance * force;
        delta.get(b).x -= dx / distance * force;
        delta.get(b).y -= dy / distance * force;
      }
      for (const { atom_id: id } of ordered) {
        const p = points.get(id), d = delta.get(id);
        p.x += Math.max(-5, Math.min(5, d.x)) - p.x * 0.002;
        p.y += Math.max(-5, Math.min(5, d.y)) - p.y * 0.002;
      }
    }
    return points;
  }

  function paintGraphs(event) {
    const comparison = byId("event-detail").querySelector(".reaction-comparison");
    const hasComplexGraph = [...event.reactants, ...event.products].some((m) => m.atoms.length > 10);
    comparison.classList.toggle("stacked-topologies", hasComplexGraph && byId("event-detail").clientWidth < 760);
    const atoms = new Map();
    const edges = new Map();
    [...event.reactants, ...event.products].forEach((m) => {
      m.atoms.forEach((a) => atoms.set(a.atom_id, a));
      m.connections.forEach((e) => edges.set(bondKey(e), e));
    });
    const sharedLayout = atoms.size <= 80 ? layout([...atoms.values()], [...edges.values()]) : null;
    const broken = new Set(event.broken_bonds.map(bondKey));
    const formed = new Set(event.formed_bonds.map(bondKey));
    const changed = new Set(event.changed_atom_ids);
    byId("event-detail").querySelectorAll(".molecule-graph").forEach((container) => {
      const side = container.dataset.moleculeSide;
      const full = event[side][Number(container.dataset.moleculeIndex)];
      let molecule = full;
      let scope = "";
      if (full.atoms.length > 40) {
        const focus = full.atoms.filter((a) => changed.has(a.atom_id));
        const selected = new Set(focus.slice(0, 40).map((a) => a.atom_id));
        for (const [a, b] of full.connections) {
          if (selected.size >= 40) break;
          if (changed.has(a) && selected.has(a)) selected.add(b);
          else if (changed.has(b) && selected.has(b)) selected.add(a);
        }
        if (!selected.size) full.atoms.slice(0, 40).forEach((a) => selected.add(a.atom_id));
        molecule = { ...full, atoms: full.atoms.filter((a) => selected.has(a.atom_id)), connections: full.connections.filter(([a, b]) => selected.has(a) && selected.has(b)) };
        scope = `<p class="topology-note graph-scope">局部连接图：${number(molecule.atoms.length)} / ${number(full.atoms.length)} 个原子，优先展示发生连接变化的原子及邻居${focus.length > 40 ? "；变化原子超过显示上限" : ""}。完整结构见 JSON。</p>`;
      }
      const width = Math.max(200, Math.floor(container.clientWidth));
      const height = molecule.atoms.length > 24 ? 500 : molecule.atoms.length > 12 ? 380 : width < 340 && molecule.atoms.length > 8 ? 320 : 210;
      const positions = sharedLayout || layout(molecule.atoms, molecule.connections);
      const ps = molecule.atoms.map((a) => positions.get(a.atom_id));
      const minX = Math.min(...ps.map((p) => p.x)), maxX = Math.max(...ps.map((p) => p.x));
      const minY = Math.min(...ps.map((p) => p.y)), maxY = Math.max(...ps.map((p) => p.y));
      const scale = Math.min((width - 68) / Math.max(65, maxX - minX), (height - 70) / Math.max(65, maxY - minY));
      const screen = new Map(molecule.atoms.map(({ atom_id: id }) => [id, { x: width / 2 + (positions.get(id).x - (minX + maxX) / 2) * scale, y: height / 2 + (positions.get(id).y - (minY + maxY) / 2) * scale }]));
      // Resolve glyph collisions in screen pixels so labels stay readable at any width.
      for (let iteration = 0; iteration < 90; iteration += 1) {
        let overlap = 0;
        for (let i = 0; i < molecule.atoms.length; i += 1) {
          for (let j = i + 1; j < molecule.atoms.length; j += 1) {
            const a = screen.get(molecule.atoms[i].atom_id), b = screen.get(molecule.atoms[j].atom_id);
            let dx = b.x - a.x, dy = b.y - a.y;
            let distance = Math.hypot(dx, dy);
            if (distance < 0.001) { dx = 1; dy = 0; distance = 1; }
            if (distance >= 49) continue;
            const push = (49 - distance) / 2;
            overlap = Math.max(overlap, push);
            a.x -= dx / distance * push; a.y -= dy / distance * push;
            b.x += dx / distance * push; b.y += dy / distance * push;
          }
        }
        for (const p of screen.values()) {
          p.x = Math.max(24, Math.min(width - 24, p.x));
          p.y = Math.max(24, Math.min(height - 24, p.y));
        }
        if (overlap < 0.1) break;
      }
      const coord = (id) => screen.get(id);
      const lines = molecule.connections.map((pair) => {
        const a = coord(pair[0]), b = coord(pair[1]);
        const kind = side === "reactants" && broken.has(bondKey(pair)) ? "broken-edge" : side === "products" && formed.has(bondKey(pair)) ? "formed-edge" : "";
        return `<line class="topology-edge ${kind}" x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}"></line>`;
      }).join("");
      const nodes = molecule.atoms.map((atom) => {
        const p = coord(atom.atom_id);
        return `<g class="topology-atom ${changed.has(atom.atom_id) ? "changed-atom" : ""}" data-select-atom="${atom.atom_id}" data-atom-id="${atom.atom_id}" role="button" tabindex="0" aria-label="查看 ${escape(atom.element)} 原子 ${atom.atom_id} 的去向" transform="translate(${p.x} ${p.y})"><circle r="21"></circle><text class="element-label" y="-3">${escape(atom.element)}</text><text class="atom-id-label" y="13">${atom.atom_id}</text></g>`;
      }).join("");
      container.innerHTML = `${scope}<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="group" aria-label="${escape(molecule.structure_id || molecule.formula)} 的${scope ? "局部" : "完整"}连接拓扑">${lines}${nodes}</svg>`;
    });
  }

  function selectAtom(id) {
    const event = model.detail;
    if (!event) return;
    model.atom = Number(id);
    const before = event.reactants.find((m) => m.atoms.some((a) => a.atom_id === model.atom));
    const after = event.products.find((m) => m.atoms.some((a) => a.atom_id === model.atom));
    const atom = before?.atoms.find((a) => a.atom_id === model.atom);
    byId("atom-destination").textContent = `${atom?.element || ""} 原子 ${model.atom}：${before?.structure_id || before?.formula || "未知"} → ${after?.structure_id || after?.formula || "未知"}`;
    byId("event-detail").querySelectorAll("[data-atom-id]").forEach((node) => node.classList.toggle("selected-atom", Number(node.dataset.atomId) === model.atom));
    byId("event-detail").querySelectorAll(".atom-table-button").forEach((node) => node.setAttribute("aria-pressed", String(Number(node.dataset.selectAtom) === model.atom)));
  }

  function init() {
    byId("event-filters").addEventListener("submit", (event) => { event.preventDefault(); model.offset = 0; clearDetail(); load(); });
    byId("event-reset").addEventListener("click", () => { byId("event-filters").reset(); model.offset = 0; clearDetail(); load(); });
    byId("event-page-size").addEventListener("change", () => { model.limit = Number(byId("event-page-size").value); model.offset = 0; clearDetail(); load(); });
    byId("event-prev").addEventListener("click", () => { model.offset = Math.max(0, model.offset - model.limit); clearDetail(); load(); });
    byId("event-next").addEventListener("click", () => { model.offset += model.limit; clearDetail(); load(); });
    byId("events-table").addEventListener("click", (event) => { const row = event.target.closest("[data-event-row]"); if (row) select(row.dataset.eventRow); });
    byId("event-detail").addEventListener("click", (event) => { const node = event.target.closest("[data-select-atom]"); if (node) selectAtom(node.dataset.selectAtom); });
    byId("event-detail").addEventListener("keydown", (event) => {
      const node = event.target.closest("g[data-select-atom]");
      if (node && (event.key === "Enter" || event.key === " ")) { event.preventDefault(); selectAtom(node.dataset.selectAtom); }
    });
    let resizeTimer;
    window.addEventListener("resize", () => {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(() => { if (model.detail) { paintGraphs(model.detail); if (model.atom !== null) selectAtom(model.atom); } }, 150);
    });
  }

  return { init, load, reset };
})();
