"use strict";

// Counts are read from the per-frame CSV outputs. Rendering never changes them.
window.CountCharts = (() => {
  const palette = ["#E7E5CC", "#9CC184", "#3C7D3F", "#1A2717"];
  const dash = ["", "9 4", "3 3", "12 4 2 4"];
  const names = { hydroxyl: "Hydroxyl candidate", ether: "Ether candidate", aldehyde: "Aldehyde candidate", ketone: "Ketone-like carbonyl", carboxyl: "Carboxyl candidate", ester: "Ester candidate" };
  const models = {};
  const escape = (text) => String(text ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const numeric = (value) => value !== null && value !== undefined && value !== "" && Number.isFinite(Number(value));
  const number = (value) => numeric(value) ? Number(value).toLocaleString("zh-CN", { maximumFractionDigits: 5 }) : "—";
  const el = (kind, name) => document.getElementById(`counts-${kind}-${name}`);
  const figureLabel = (category) => names[category.key] || category.label || category.key.replaceAll("-", "–");

  function shell(kind) {
    const functional = kind === "functional";
    return `
      <div class="table-heading"><div><h3>${functional ? "官能团候选数量" : "元素对连接数量"}</h3><p>${functional ? "依据局部连接拓扑逐帧匹配，不代表已确认的化学键型或化学身份。" : "直接统计 bonds 文件中的无向连接，每对原子只计一次；C–O 等标签不区分单、双键。"}</p></div></div>
      ${functional ? '<p class="counts-basis-note">羧酸中的 OH 归入羧基，酯中的氧归入酯基，不再重复计入普通羟基和醚；“酮型羰基”不表示全部 C=O。水和 OH 自由基不计羟基。</p>' : ''}
      <p id="counts-${kind}-coverage-note" class="counts-basis-note hidden"></p>
      <div id="counts-${kind}-status" class="counts-status" role="status">打开此页时读取统计。</div>
      <button id="counts-${kind}-backfill" type="button" class="small-button primary-small hidden">补算当前范围统计</button>
      <div id="counts-${kind}-content" class="hidden">
        <form id="counts-${kind}-controls" class="counts-controls">
          <div><label for="counts-${kind}-metric">统计指标</label><select id="counts-${kind}-metric">${functional ? '<option value="count">官能团候选个数</option><option value="molecules">含候选的分子数</option>' : '<option value="count">现存连接数</option><option value="formed">相邻帧新增数</option><option value="broken">相邻帧消失数</option>'}</select></div>
          <div><label for="counts-${kind}-from">起始时间 / ps</label><input id="counts-${kind}-from" type="number" min="0" step="any" placeholder="全部"></div>
          <div><label for="counts-${kind}-to">结束时间 / ps</label><input id="counts-${kind}-to" type="number" min="0" step="any" placeholder="全部"></div>
          <div><label for="counts-${kind}-smooth">移动平均窗口 / ps</label><input id="counts-${kind}-smooth" type="number" min="0" step="any" value="${functional ? 5 : 0}"></div>
          <div class="counts-control-actions"><button type="submit" class="small-button">应用范围</button><button id="counts-${kind}-range-reset" type="button" class="small-button">全部时间</button></div>
        </form>
        <p class="counts-help">可同时比较最多 4 类。窗口为 0 时显示原始逐帧线；平滑时保留浅色原始线，粗线表示居中时间窗均值。下载 CSV 始终保留完整原始统计。</p>
        <fieldset class="counts-series"><legend>选择曲线</legend><div id="counts-${kind}-series"></div></fieldset>
        <p id="counts-${kind}-error" class="inline-error hidden" role="alert"></p>
        <div class="counts-figure-heading"><span id="counts-${kind}-coverage"></span><div class="counts-downloads"><a id="counts-${kind}-csv" class="small-button" download>完整 CSV</a><button id="counts-${kind}-svg" type="button" class="small-button">导出 SVG</button><button id="counts-${kind}-png" type="button" class="small-button">导出 PNG</button></div></div>
        <div id="counts-${kind}-plot" class="counts-plot"></div>
        <label class="counts-frame-label" for="counts-${kind}-frame">查看具体保存帧（也可在图上移动鼠标）</label>
        <input id="counts-${kind}-frame" class="counts-frame-slider" type="range" min="0" max="0" value="0" step="1" aria-label="查看统计保存帧">
        <p id="counts-${kind}-position" class="counts-position"></p>
        <div id="counts-${kind}-values" class="table-scroll counts-values"></div>
        <details class="counts-provenance"><summary>统计口径与来源</summary><pre id="counts-${kind}-metadata"></pre></details>
      </div>`;
  }

  function init() {
    for (const kind of ["bonds", "functional"]) {
      document.getElementById(`count-${kind}`).innerHTML = shell(kind);
      models[kind] = { jobId: null, request: 0, timer: null, payload: null, selected: new Set(), visibleRows: [], chart: null };
      el(kind, "controls").addEventListener("submit", (event) => { event.preventDefault(); render(kind); });
      el(kind, "metric").addEventListener("change", () => {
        if (kind === "bonds") el(kind, "smooth").value = el(kind, "metric").value === "count" ? "0" : "5";
        load(models[kind].jobId, kind, true);
      });
      el(kind, "smooth").addEventListener("change", () => render(kind));
      el(kind, "range-reset").addEventListener("click", () => { el(kind, "from").value = ""; el(kind, "to").value = ""; render(kind); });
      el(kind, "series").addEventListener("change", (event) => {
        const input = event.target.closest("[data-count-series]");
        if (!input) return;
        const selected = models[kind].selected;
        if (input.checked && selected.size >= 4) {
          input.checked = false;
          showError(kind, "最多同时显示 4 类；先取消一类，再选择其他曲线。");
          return;
        }
        if (input.checked) selected.add(input.dataset.countSeries); else selected.delete(input.dataset.countSeries);
        render(kind);
      });
      el(kind, "frame").addEventListener("input", () => selectFrame(kind, Number(el(kind, "frame").value)));
      el(kind, "backfill").addEventListener("click", () => backfill(kind));
      el(kind, "svg").addEventListener("click", () => exportFigure(kind, false));
      el(kind, "png").addEventListener("click", () => exportFigure(kind, true));
    }
  }

  function reset() {
    for (const [kind, model] of Object.entries(models)) {
      model.request += 1;
      clearTimeout(model.timer);
      model.jobId = null;
      model.payload = null;
      model.selected.clear();
      model.visibleRows = [];
      model.chart = null;
      el(kind, "content").classList.add("hidden");
      el(kind, "backfill").classList.add("hidden");
      el(kind, "status").textContent = "打开此页时读取统计。";
      el(kind, "coverage-note").classList.add("hidden");
      el(kind, "controls").reset();
      el(kind, "plot").replaceChildren();
    }
  }

  async function request(url, options) {
    const response = await fetch(url, options);
    let payload;
    try { payload = await response.json(); } catch (_) { payload = {}; }
    if (!response.ok) throw new Error(typeof payload.detail === "string" ? payload.detail : `读取统计失败（${response.status}）`);
    return payload;
  }

  function showError(kind, message = "") {
    el(kind, "error").textContent = message;
    el(kind, "error").classList.toggle("hidden", !message);
  }

  async function load(jobId, kind, force = false) {
    if (!jobId || !models[kind]) return;
    const model = models[kind];
    const metric = el(kind, "metric").value;
    if (!force && model.jobId === jobId && model.payload?.status === "completed" && model.payload.metric === metric) { render(kind); return; }
    clearTimeout(model.timer);
    model.jobId = jobId;
    const ticket = ++model.request;
    el(kind, "status").textContent = "读取逐帧统计…";
    el(kind, "content").classList.add("hidden");
    el(kind, "backfill").classList.add("hidden");
    showError(kind);
    try {
      const payload = await request(`/api/jobs/${encodeURIComponent(jobId)}/counts?${new URLSearchParams({ kind, metric, max_points: 5000 })}`);
      if (ticket !== model.request || jobId !== model.jobId) return;
      if (payload.status === "running") {
        const progress = payload.progress;
        el(kind, "status").textContent = `${payload.message || "正在补算逐帧统计"}${progress ? ` · ${typeof progress === "number" ? number(progress) + "%" : number(progress.frames_processed ?? progress.frame_count) + " 帧"}` : ""}。反应结果仍可查看。`;
        model.timer = setTimeout(() => load(jobId, kind, true), 1500);
        return;
      }
      if (payload.status !== "completed") {
        model.payload = null;
        el(kind, "status").textContent = payload.error || payload.reason || payload.message || (payload.status === "failed" ? "数量统计失败，可重新补算。" : "此任务尚无数量统计。补算只读取本任务对应范围，不重做反应分析。");
        el(kind, "backfill").classList.remove("hidden");
        el(kind, "backfill").disabled = false;
        return;
      }
      const wasComplete = model.payload?.status === "completed";
      model.payload = payload;
      const keys = new Set((payload.categories || []).map((c) => c.key));
      model.selected = new Set([...model.selected].filter((key) => keys.has(key)));
      if (!model.selected.size) (payload.categories || []).slice(0, 4).forEach((c) => model.selected.add(c.key));
      el(kind, "series").innerHTML = (payload.categories || []).map((category) => `<label><input type="checkbox" data-count-series="${escape(category.key)}" ${model.selected.has(category.key) ? "checked" : ""}><span>${escape(category.label || category.key)}</span></label>`).join("");
      el(kind, "content").classList.remove("hidden");
      el(kind, "status").textContent = `${number(payload.frame_count)} 个保存帧的统计已读取。${kind === "functional" ? "官能团结果为指定规则下的拓扑候选。" : "连接数包含文件中所有已列出的连接。"}`;
      const csv = new URL(payload.csv_url || "", window.location.href);
      const csvAvailable = csv.origin === window.location.origin && csv.pathname.startsWith("/api/jobs/");
      el(kind, "csv").classList.toggle("hidden", !csvAvailable);
      if (csvAvailable) el(kind, "csv").href = csv.href;
      el(kind, "metadata").textContent = JSON.stringify(payload.metadata || {}, null, 2);
      const overcoordinated = payload.metadata?.coverage?.overcoordinated_atom_observations;
      const hasExcludedEnvironments = kind === "functional" && numeric(overcoordinated) && Number(overcoordinated) > 0;
      el(kind, "coverage-note").classList.toggle("hidden", !hasExcludedEnvironments);
      if (hasExcludedEnvironments) el(kind, "coverage-note").textContent = `本次统计范围内有 ${number(overcoordinated)} 次跨帧原子观测超过 H≤1、C≤4、O≤2 的局部配位规则，相关环境未计入官能团候选；原始连接仍完整统计。这是观测次数，不是不同原子的个数。`;
      el(kind, "smooth").disabled = Boolean(payload.downsampled);
      render(kind);
      if (!wasComplete) window.dispatchEvent(new CustomEvent("counts-ready", { detail: { jobId } }));
    } catch (cause) {
      if (ticket !== model.request || jobId !== model.jobId) return;
      el(kind, "status").textContent = cause.message;
      el(kind, "backfill").classList.remove("hidden");
      el(kind, "backfill").disabled = false;
    }
  }

  async function backfill(kind) {
    const model = models[kind], jobId = model.jobId;
    if (!jobId) return;
    el(kind, "backfill").disabled = true;
    el(kind, "status").textContent = "提交当前任务范围的数量补算…";
    try {
      await request(`/api/jobs/${encodeURIComponent(jobId)}/counts`, { method: "POST" });
      if (jobId === model.jobId) await load(jobId, kind, true);
    } catch (cause) {
      if (jobId !== model.jobId) return;
      el(kind, "status").textContent = cause.message;
      el(kind, "backfill").disabled = false;
    }
  }

  function movingAverage(rows, key, windowPs) {
    if (!windowPs) return rows.map((row) => numeric(row[key]) ? Number(row[key]) : null);
    let left = 0, right = 0, sum = 0, count = 0;
    return rows.map((row) => {
      const time = Number(row.time_ps);
      while (right < rows.length && Number(rows[right].time_ps) <= time + windowPs / 2) {
        if (numeric(rows[right][key])) { sum += Number(rows[right][key]); count += 1; }
        right += 1;
      }
      while (left < right && Number(rows[left].time_ps) < time - windowPs / 2) {
        if (numeric(rows[left][key])) { sum -= Number(rows[left][key]); count -= 1; }
        left += 1;
      }
      return numeric(row[key]) && count ? sum / count : null;
    });
  }

  function niceUpper(maximum) {
    const rough = Math.max(1, maximum) / 5;
    const unit = 10 ** Math.floor(Math.log10(rough));
    const ratio = rough / unit;
    const step = (ratio <= 1 ? 1 : ratio <= 2 ? 2 : ratio <= 5 ? 5 : 10) * unit;
    return Math.max(step, Math.ceil(Math.max(1, maximum) / step) * step);
  }

  function render(kind) {
    const model = models[kind], payload = model.payload;
    if (!payload || payload.status !== "completed") return;
    const from = el(kind, "from").value.trim(), to = el(kind, "to").value.trim();
    const windowPs = payload.downsampled ? 0 : Number(el(kind, "smooth").value || 0);
    if ((from && !numeric(from)) || (to && !numeric(to)) || !Number.isFinite(windowPs) || windowPs < 0 || (from && to && Number(from) > Number(to))) {
      showError(kind, "请填写有效时间范围和非负平滑窗口，起始时间不能大于结束时间。");
      return;
    }
    showError(kind);
    const allRows = (payload.rows || []).filter((row) => numeric(row.time_ps));
    const rowIndexes = new Map(allRows.map((row, index) => [row, index]));
    const rows = allRows.filter((row) => (!from || Number(row.time_ps) >= Number(from)) && (!to || Number(row.time_ps) <= Number(to)));
    const categories = (payload.categories || []).filter((c) => model.selected.has(c.key));
    model.visibleRows = rows;
    model.chart = null;
    el(kind, "coverage").textContent = `共 ${number(payload.frame_count)} 帧 · 当前范围 ${number(rows.length)} 个点${payload.downsampled ? " · 全时段极值抽样展示，平滑已关闭；完整数据见 CSV" : " · 完整逐帧展示"}${windowPs ? ` · ${number(windowPs)} ps 居中移动平均` : " · 原始线"}`;
    el(kind, "frame").max = Math.max(0, rows.length - 1);
    el(kind, "frame").value = "0";
    el(kind, "frame").disabled = !rows.length || !categories.length;
    el(kind, "svg").disabled = !rows.length || !categories.length;
    el(kind, "png").disabled = !rows.length || !categories.length;
    if (!rows.length || !categories.length) {
      el(kind, "plot").innerHTML = `<p class="counts-empty">${rows.length ? "请选择至少一条曲线。" : "当前时间范围内没有保存帧。"}</p>`;
      el(kind, "position").textContent = "";
      el(kind, "values").replaceChildren();
      return;
    }
    const chart = { rows, categories, windowPs, x0: Number(rows[0].time_ps), x1: Number(rows[rows.length - 1].time_ps), left: 110, right: 1160, top: 45, bottom: 700 };
    if (chart.x0 === chart.x1) { chart.x0 -= 0.5; chart.x1 += 0.5; }
    chart.x = (time) => chart.left + (time - chart.x0) / (chart.x1 - chart.x0) * (chart.right - chart.left);
    chart.series = categories.map((category, index) => {
      // Zooming changes the visible range, not the averaging neighborhood.
      const complete = movingAverage(allRows, category.key, windowPs);
      return { category, color: palette[index], dash: dash[index], values: rows.map((row) => complete[rowIndexes.get(row)]) };
    });
    let maximum = 0;
    for (const row of rows) for (const c of categories) if (numeric(row[c.key])) maximum = Math.max(maximum, Number(row[c.key]));
    for (const series of chart.series) for (const value of series.values) if (numeric(value)) maximum = Math.max(maximum, Number(value));
    chart.y1 = niceUpper(maximum * 1.07);
    chart.y = (value) => chart.bottom - value / chart.y1 * (chart.bottom - chart.top);
    model.chart = chart;
    const axes = [];
    for (let i = 0; i <= 5; i += 1) {
      const value = chart.y1 * i / 5, y = chart.y(value);
      axes.push(`<line x1="${chart.left}" x2="${chart.right}" y1="${y}" y2="${y}" stroke="#ededed" stroke-width="0.7"/><line x1="${chart.left}" x2="${chart.left + 8}" y1="${y}" y2="${y}" stroke="#222"/><text x="${chart.left - 13}" y="${y + 5}" text-anchor="end">${escape(number(value))}</text>`);
      const time = chart.x0 + (chart.x1 - chart.x0) * i / 5, x = chart.x(time);
      axes.push(`<line x1="${x}" x2="${x}" y1="${chart.bottom}" y2="${chart.bottom - 8}" stroke="#222"/><text x="${x}" y="${chart.bottom + 26}" text-anchor="middle">${escape(number(time))}</text>`);
    }
    const pathFor = (key, values) => {
      let path = "", connected = false;
      rows.forEach((row, index) => {
        const value = values ? values[index] : row[key];
        if (!numeric(value)) { connected = false; return; }
        path += `${connected ? "L" : "M"}${chart.x(Number(row.time_ps)).toFixed(2)},${chart.y(Number(value)).toFixed(2)} `;
        connected = true;
      });
      return path;
    };
    const lines = chart.series.map((series) => `${windowPs ? `<path d="${pathFor(series.category.key)}" fill="none" stroke="${series.color}" stroke-width="0.85" stroke-opacity="0.5"/>` : ""}<path d="${pathFor(series.category.key, series.values)}" fill="none" stroke="${series.color}" stroke-width="${windowPs ? 3 : 2.8}" stroke-dasharray="${series.dash}" stroke-linejoin="round"/>${rows.length === 1 && numeric(series.values[0]) ? `<circle cx="${chart.x(Number(rows[0].time_ps))}" cy="${chart.y(series.values[0])}" r="4" fill="${series.color}"/>` : ""}`).join("");
    const legendWidth = 245, legendHeight = 24 * categories.length + 8;
    const corners = [{ x: chart.left + 18, y: chart.top + 17 }, { x: chart.right - legendWidth - 16, y: chart.top + 17 }, { x: chart.left + 18, y: chart.bottom - legendHeight - 16 }, { x: chart.right - legendWidth - 16, y: chart.bottom - legendHeight - 16 }];
    const score = (corner) => chart.series.reduce((total, series) => total + series.values.reduce((hits, value, i) => {
      if (!numeric(value)) return hits;
      const x = chart.x(Number(rows[i].time_ps)), y = chart.y(Number(value));
      return hits + Number(x >= corner.x - 8 && x <= corner.x + legendWidth && y >= corner.y - 12 && y <= corner.y + legendHeight);
    }, 0), 0);
    const corner = corners.sort((a, b) => score(a) - score(b))[0];
    const legend = chart.series.map((series, index) => `<g transform="translate(${corner.x} ${corner.y + index * 24})"><line x1="0" x2="30" y1="0" y2="0" stroke="${series.color}" stroke-width="3" stroke-dasharray="${series.dash}"/><text x="40" y="5" font-size="13">${escape(figureLabel(series.category))}</text></g>`).join("");
    const ylabel = kind === "functional" ? (payload.metric === "molecules" ? "Number of molecules" : "Number of candidate groups") : payload.metric === "formed" ? "New connections per interval" : payload.metric === "broken" ? "Lost connections per interval" : "Number of connections";
    const clipId = `counts-${kind}-clip`;
    el(kind, "plot").innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="800" viewBox="0 0 1200 800" role="img" aria-label="${escape(kind === "functional" ? "官能团候选" : "元素对连接")}数量随时间变化折线图" style="font-family:'Times New Roman',serif;font-size:14.7px;fill:#222;background:white"><defs><clipPath id="${clipId}"><rect x="${chart.left - 1}" y="${chart.top - 1}" width="${chart.right - chart.left + 2}" height="${chart.bottom - chart.top + 2}"/></clipPath></defs><rect width="1200" height="800" fill="white"/>${axes.join("")}<path d="M${chart.left},${chart.top}V${chart.bottom}H${chart.right}" fill="none" stroke="#222" stroke-width="1.3"/><text x="635" y="768" text-anchor="middle" font-size="16">Time (ps)</text><text transform="translate(28 373) rotate(-90)" text-anchor="middle" font-size="16">${ylabel}</text><g clip-path="url(#${clipId})">${lines}</g>${legend}<g data-count-cursor data-export-exclude="true" visibility="hidden" pointer-events="none"><line y1="${chart.top}" y2="${chart.bottom}" stroke="#667565" stroke-width="1" stroke-dasharray="4 4"/>${chart.series.map((s) => `<circle data-cursor-key="${escape(s.category.key)}" r="4" fill="${s.color}" stroke="#555" stroke-width="0.8"/>`).join("")}</g><rect data-count-hover data-export-exclude="true" x="${chart.left}" y="${chart.top}" width="${chart.right - chart.left}" height="${chart.bottom - chart.top}" fill="transparent"/></svg>`;
    const svg = el(kind, "plot").querySelector("svg");
    svg.querySelector("[data-count-hover]").addEventListener("pointermove", (event) => {
      const point = svg.createSVGPoint(); point.x = event.clientX; point.y = event.clientY;
      const local = point.matrixTransform(svg.getScreenCTM().inverse());
      const time = chart.x0 + (local.x - chart.left) / (chart.right - chart.left) * (chart.x1 - chart.x0);
      let low = 0, high = rows.length - 1;
      while (low < high) { const middle = Math.floor((low + high) / 2); if (Number(rows[middle].time_ps) < time) low = middle + 1; else high = middle; }
      if (low && Math.abs(Number(rows[low - 1].time_ps) - time) < Math.abs(Number(rows[low].time_ps) - time)) low -= 1;
      el(kind, "frame").value = String(low);
      selectFrame(kind, low);
    });
    selectFrame(kind, 0);
  }

  function selectFrame(kind, index) {
    const model = models[kind], chart = model.chart;
    if (!chart || !chart.rows[index]) return;
    const row = chart.rows[index];
    const slider = el(kind, "frame");
    slider.setAttribute("aria-valuetext", `${number(row.time_ps)} ps，step ${number(row.timestep)}`);
    el(kind, "position").textContent = `第 ${number(index + 1)} / ${number(chart.rows.length)} 个显示帧 · step ${number(row.timestep)} · ${number(row.time_ps)} ps${chart.windowPs ? "；下表为原始计数，均值另列" : ""}`;
    el(kind, "values").innerHTML = `<table><thead><tr><th>类别</th><th>原始计数</th>${chart.windowPs ? '<th>窗口均值</th>' : ''}</tr></thead><tbody>${chart.series.map((s) => `<tr><td>${escape(s.category.label || s.category.key)}</td><td>${numeric(row[s.category.key]) ? number(row[s.category.key]) : '—（无前一帧）'}</td>${chart.windowPs ? `<td>${number(s.values[index])}</td>` : ''}</tr>`).join("")}</tbody></table>`;
    const cursor = el(kind, "plot").querySelector("[data-count-cursor]");
    const x = chart.x(Number(row.time_ps));
    cursor.setAttribute("visibility", "visible");
    cursor.querySelector("line").setAttribute("x1", String(x));
    cursor.querySelector("line").setAttribute("x2", String(x));
    cursor.querySelectorAll("circle").forEach((circle, i) => {
      const value = chart.series[i].values[index];
      circle.setAttribute("visibility", numeric(value) ? "visible" : "hidden");
      if (numeric(value)) { circle.setAttribute("cx", String(x)); circle.setAttribute("cy", String(chart.y(Number(value)))); }
    });
  }

  function download(blob, name) {
    const url = URL.createObjectURL(blob), link = document.createElement("a");
    link.href = url; link.download = name; document.body.append(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  async function pngAt300Dpi(blob) {
    const bytes = new Uint8Array(await blob.arrayBuffer());
    const chunk = new Uint8Array(21), view = new DataView(chunk.buffer);
    view.setUint32(0, 9); chunk.set([112, 72, 89, 115], 4);
    view.setUint32(8, 11811); view.setUint32(12, 11811); chunk[16] = 1;
    let crc = 0xffffffff;
    for (const value of chunk.slice(4, 17)) { crc ^= value; for (let bit = 0; bit < 8; bit += 1) crc = (crc >>> 1) ^ (crc & 1 ? 0xedb88320 : 0); }
    view.setUint32(17, (crc ^ 0xffffffff) >>> 0);
    const parts = [bytes.slice(0, 8)];
    let offset = 8;
    while (offset + 12 <= bytes.length) {
      const length = new DataView(bytes.buffer, bytes.byteOffset + offset, 4).getUint32(0);
      const type = String.fromCharCode(...bytes.slice(offset + 4, offset + 8));
      const end = offset + length + 12;
      if (end > bytes.length) throw new Error("PNG 数据不完整");
      if (type !== "pHYs") parts.push(bytes.slice(offset, end));
      if (type === "IHDR") parts.push(chunk);
      offset = end;
    }
    return new Blob(parts, { type: "image/png" });
  }

  async function exportFigure(kind, png) {
    const model = models[kind], original = el(kind, "plot").querySelector("svg");
    if (!original || !model.chart) return;
    const svg = original.cloneNode(true);
    svg.querySelectorAll("[data-export-exclude]").forEach((node) => node.remove());
    const text = new XMLSerializer().serializeToString(svg);
    const blob = new Blob([text], { type: "image/svg+xml;charset=utf-8" });
    const filename = `${kind}_${model.payload.metric}_${model.jobId}_${Date.now()}`;
    if (!png) { download(blob, `${filename}.svg`); return; }
    const button = el(kind, "png"); button.disabled = true;
    const url = URL.createObjectURL(blob);
    try {
      const image = new Image();
      await new Promise((resolve, reject) => { image.onload = resolve; image.onerror = () => reject(new Error("图像转换失败，可导出 SVG")); image.src = url; });
      const canvas = document.createElement("canvas"); canvas.width = 1200; canvas.height = 800;
      const context = canvas.getContext("2d"); context.fillStyle = "white"; context.fillRect(0, 0, 1200, 800); context.drawImage(image, 0, 0, 1200, 800);
      const rendered = await new Promise((resolve) => canvas.toBlob(resolve, "image/png"));
      if (!rendered) throw new Error("PNG 导出失败，可导出 SVG");
      download(await pngAt300Dpi(rendered), `${filename}_300dpi.png`);
    } catch (cause) { showError(kind, cause.message); }
    finally { URL.revokeObjectURL(url); button.disabled = false; }
  }

  return { init, reset, load };
})();
