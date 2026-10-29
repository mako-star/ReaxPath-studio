"use strict";

const CHUNK_SIZE = 8 * 1024 * 1024;
const state = {
  mode: "local",
  jobId: null,
  network: null,
  pollTimer: null,
  polling: false,
  resultsReady: false,
  speciesLoaded: false,
};

const $ = (selector) => document.querySelector(selector);
const $$ = (selector) => Array.from(document.querySelectorAll(selector));

document.addEventListener("DOMContentLoaded", () => {
  bindControls();
  checkHealth();
  restoreJobFromUrl();
});

function bindControls() {
  $$(".segment").forEach((button) => {
    button.addEventListener("click", () => setSourceMode(button.dataset.mode));
  });
  $$(".result-tab").forEach((button) => {
    button.addEventListener("click", () => setResultTab(button.dataset.tab));
  });
  $("#bonds-file").addEventListener("change", updateFileLabel);
  $("#analysis-form").addEventListener("submit", submitAnalysis);
  window.EventBrowser.init();
  window.CountCharts.init();
  window.addEventListener("counts-ready", async (event) => {
    if (event.detail.jobId !== state.jobId) return;
    try {
      const outputs = await api(`/api/jobs/${encodeURIComponent(state.jobId)}/outputs`);
      if (event.detail.jobId === state.jobId) renderOutputs(outputs.files || []);
    } catch (_) { /* The chart still offers direct CSV downloads. */ }
  });
}

async function checkHealth() {
  const target = $("#service-state");
  try {
    const health = await api("/healthz");
    target.className = "service-state online";
    target.querySelector("span:last-child").textContent = `服务在线 / 上限 ${formatBytes(health.max_file_size)}`;
  } catch (error) {
    target.className = "service-state offline";
    target.querySelector("span:last-child").textContent = "服务不可用";
  }
}

function setSourceMode(mode) {
  state.mode = mode;
  $$(".segment").forEach((button) => {
    const active = button.dataset.mode === mode;
    button.classList.toggle("active", active);
    button.setAttribute("aria-selected", String(active));
  });
  $("#local-source").classList.toggle("hidden", mode !== "local");
  $("#upload-source").classList.toggle("hidden", mode !== "upload");
  hideError();
}

function setResultTab(name) {
  $$(".result-tab").forEach((button) => {
    const active = button.dataset.tab === name;
    button.classList.toggle("active", active);
    button.setAttribute("aria-selected", String(active));
  });
  $$(".tab-panel").forEach((panel) => {
    panel.classList.toggle("hidden", panel.dataset.panel !== name);
  });
  if (state.resultsReady && name === "network") loadNetwork();
  if (state.resultsReady && name === "species") loadSpecies();
  if (state.resultsReady && (name === "bonds" || name === "functional")) window.CountCharts.load(state.jobId, name);
}

function updateFileLabel() {
  const file = $("#bonds-file").files[0];
  $("#file-name").textContent = file ? file.name : "选择 bonds 文件";
  $("#file-size").textContent = file ? formatBytes(file.size) : "8 MiB 分块";
}

async function submitAnalysis(event) {
  event.preventDefault();
  hideError();
  clearInterval(state.pollTimer);
  const runButton = $("#run-button");
  runButton.disabled = true;

  try {
    const config = readConfig();
    showJobShell(state.mode === "local" ? basename($("#bonds-path").value) : selectedFilename());
    let job;
    if (state.mode === "local") {
      const path = $("#bonds-path").value.trim();
      if (!path) throw new Error("请输入 bonds 文件的绝对路径");
      updateProgress(1, "提交任务", "验证本机文件路径", "");
      job = await createLocalJob(path, config);
    } else {
      const file = $("#bonds-file").files[0];
      if (!file) throw new Error("请选择 bonds 文件");
      const uploadId = await uploadInChunks(file);
      updateProgress(99, "提交任务", "创建流式分析任务", formatBytes(file.size));
      job = await createUploadJob(uploadId, config);
    }
    state.jobId = job.job_id;
    const url = new URL(window.location.href);
    url.searchParams.set("job", job.job_id);
    url.searchParams.delete("event");
    window.history.replaceState(null, "", url);
    $("#job-short-id").textContent = job.job_id.slice(0, 8).toUpperCase();
    runButton.textContent = "分析进行中";
    await pollJob();
  } catch (error) {
    showError(error.message || String(error));
    $("#job-status").textContent = "失败";
    $("#job-status").className = "status-badge failed";
    updateProgress(0, "任务中止", error.message || String(error), "");
    runButton.disabled = false;
    runButton.textContent = "开始分析";
  }
}

function readConfig() {
  const mapping = $("#element-mapping").value.trim();
  if (!mapping) throw new Error("请输入 LAMMPS 原子类型映射");
  const scale = Number($("#timestep-scale").value);
  const recrossing = Number($("#recrossing-window").value);
  const intervalText = $("#expected-interval").value.trim();
  if (!Number.isFinite(scale) || !(scale > 0)) throw new Error("时间步换算必须是大于 0 的有限数值");
  if (!Number.isFinite(recrossing) || !(recrossing >= 0)) throw new Error("往返窗口必须是非负有限数值");
  if (intervalText && (!Number.isSafeInteger(Number(intervalText)) || Number(intervalText) <= 0)) throw new Error("保存帧间隔必须是正整数");
  return {
    element_mapping: mapping,
    timestep_scale: scale,
    recrossing_window: recrossing,
    expected_interval: intervalText ? Number(intervalText) : null,
    generate_identifiers: $("#generate-identifiers").checked,
  };
}

function configParams(config) {
  const params = new URLSearchParams({
    element_mapping: config.element_mapping,
    timestep_scale: String(config.timestep_scale),
    recrossing_window: String(config.recrossing_window),
    generate_identifiers: String(config.generate_identifiers),
  });
  if (config.expected_interval !== null) params.set("expected_interval", config.expected_interval);
  return params;
}

async function createLocalJob(path, config) {
  const params = configParams(config);
  params.set("bonds_path", path);
  return api(`/api/jobs/local?${params}`, { method: "POST" });
}

async function createUploadJob(uploadId, config) {
  const params = configParams(config);
  params.set("upload_id", uploadId);
  return api(`/api/jobs?${params}`, { method: "POST" });
}

async function uploadInChunks(file) {
  const resumeKey = `reaxpath-upload:${file.name}:${file.size}:${file.lastModified}`;
  let uploadId = localStorage.getItem(resumeKey);
  let status = null;
  if (uploadId) {
    try {
      status = await api(`/api/uploads/${encodeURIComponent(uploadId)}`);
      if (status.total_size !== file.size || status.filename !== file.name) status = null;
    } catch (error) {
      localStorage.removeItem(resumeKey);
      uploadId = null;
    }
  }

  if (!status) {
    const params = new URLSearchParams({
      filename: file.name,
      total_size: String(file.size),
      chunk_size: String(CHUNK_SIZE),
    });
    status = await api(`/api/uploads?${params}`, { method: "POST" });
    uploadId = status.upload_id;
    localStorage.setItem(resumeKey, uploadId);
  }

  if (!status.completed) {
    const missing = status.missing_chunks || Array.from({ length: status.total_chunks }, (_, index) => index);
    for (let order = 0; order < missing.length; order += 1) {
      const index = missing[order];
      const start = index * status.chunk_size;
      const end = Math.min(start + status.chunk_size, file.size);
      const chunk = file.slice(start, end);
      status = await api(`/api/uploads/${encodeURIComponent(uploadId)}/chunks/${index}`, {
        method: "PUT",
        headers: { "Content-Type": "application/octet-stream" },
        body: chunk,
      });
      updateProgress(
        Math.min(98, status.progress * 0.98),
        "分块上传",
        `已提交 ${status.received_chunks} / ${status.total_chunks} 个分块`,
        `${formatBytes(Math.min(end, file.size))} / ${formatBytes(file.size)}`,
      );
    }
    status = await api(`/api/uploads/${encodeURIComponent(uploadId)}/complete`, { method: "POST" });
  }
  localStorage.removeItem(resumeKey);
  return uploadId;
}

async function pollJob() {
  const runButton = $("#run-button");
  const pollOnce = async () => {
    if (state.polling) return;
    state.polling = true;
    try {
      const job = await api(`/api/jobs/${encodeURIComponent(state.jobId)}`);
      renderJobStatus(job);
      if (job.status === "completed") {
        clearInterval(state.pollTimer);
        await loadResults(job);
        runButton.disabled = false;
        runButton.textContent = "再次分析";
      } else if (job.status === "failed") {
        clearInterval(state.pollTimer);
        throw new Error(job.error || "分析任务失败");
      }
    } catch (error) {
      clearInterval(state.pollTimer);
      showError(error.message || String(error));
      $("#job-status").textContent = "失败";
      $("#job-status").className = "status-badge failed";
      runButton.disabled = false;
      runButton.textContent = "开始分析";
    } finally {
      state.polling = false;
    }
  };
  await pollOnce();
  if (!runButton.disabled) return;
  state.pollTimer = setInterval(pollOnce, 1000);
}

function renderJobStatus(job) {
  const labels = { queued: "排队", running: "分析中", completed: "已完成", failed: "失败" };
  const badge = $("#job-status");
  badge.textContent = labels[job.status] || job.status;
  badge.className = `status-badge ${job.status}`;
  const summary = job.summary || {};
  const volume = summary.bytes_total
    ? `${formatBytes(summary.bytes_processed || 0)} / ${formatBytes(summary.bytes_total)}`
    : "";
  updateProgress(job.progress || 0, job.stage || "等待", job.message || "", volume);
}

async function loadResults(job) {
  const outputs = await api(`/api/jobs/${state.jobId}/outputs`);
  state.resultsReady = true;
  $("#results").classList.remove("hidden");
  renderMetrics(job.summary || {});
  renderOutputs(outputs.files || []);
  setResultTab("events");
  await window.EventBrowser.load(state.jobId);
}

async function loadNetwork() {
  if (state.network) return;
  const jobId = state.jobId;
  $("#network-canvas").textContent = "读取汇总网络…";
  try {
    const network = await api(`/api/jobs/${jobId}/outputs/reaction_network`);
    if (jobId !== state.jobId) return;
    state.network = network;
    renderNetwork(network);
  } catch (error) { $("#network-canvas").textContent = error.message; }
}

async function loadSpecies() {
  if (state.speciesLoaded) return;
  const jobId = state.jobId;
  $("#species-count-label").textContent = "读取物种目录…";
  try {
    const species = await api(`/api/jobs/${jobId}/preview/molecule_registry?limit=500`);
    if (jobId !== state.jobId) return;
    renderSpecies(species);
    state.speciesLoaded = true;
  } catch (error) { $("#species-count-label").textContent = error.message; }
}

async function restoreJobFromUrl() {
  const jobId = new URL(window.location.href).searchParams.get("job");
  if (!jobId) return;
  try {
    const job = await api(`/api/jobs/${encodeURIComponent(jobId)}`);
    showJobShell(job.filename);
    state.jobId = job.job_id;
    const mapping = job.element_mapping || {};
    $("#element-mapping").value = mapping._raw || Object.entries(mapping).map(([type, symbol]) => `${type}:${symbol}`).join(",");
    const config = job.config || {};
    $("#timestep-scale").value = config.timestep_scale ?? 1;
    $("#expected-interval").value = config.expected_interval ?? "";
    $("#recrossing-window").value = config.recrossing_window ?? 0;
    $("#generate-identifiers").checked = config.generate_identifiers !== false;
    $("#job-short-id").textContent = job.job_id.slice(0, 8).toUpperCase();
    $("#run-button").disabled = true;
    await pollJob();
  } catch (error) { showError(`无法恢复任务：${error.message}`); }
}

function renderMetrics(summary) {
  const metrics = [
    ["轨迹帧", summary.frames],
    ["物种", summary.species],
    ["检测事件", summary.events_total],
    ["聚合反应", summary.reactions],
  ];
  $("#metrics").innerHTML = metrics.map(([label, value]) => (
    `<div class="metric"><span>${escapeHtml(label)}</span><strong>${formatNumber(value)}</strong></div>`
  )).join("");
}

function renderNetwork(network) {
  const reactions = [...(network.reactions || [])]
    .sort((a, b) => b.event_count - a.event_count)
    .slice(0, 18);
  const speciesMap = new Map((network.species || []).map((item) => [item.id, item]));
  const referenced = new Set();
  reactions.forEach((reaction) => {
    reaction.reactants.forEach((item) => referenced.add(item.species_key));
    reaction.products.forEach((item) => referenced.add(item.species_key));
  });
  let species = [...referenced].map((key) => speciesMap.get(key)).filter(Boolean);
  species.sort((a, b) => b.total_observations - a.total_observations);
  // Every displayed reaction keeps all of its reactant/product participants.
  const visibleSpecies = new Set(species.map((item) => item.id));
  const canvas = $("#network-canvas");
  $("#network-caption").textContent = `${network.species?.length || 0} 个结构 / ${network.reactions?.length || 0} 个聚合反应；展示前 ${reactions.length} 个高频反应及全部参与结构。逐次证据见反应事件页。`;

  if (!reactions.length || !species.length) {
    canvas.innerHTML = '<div class="network-empty">当前参数下未识别到已接受的反应事件</div>';
    $("#reaction-list").innerHTML = "";
    return;
  }

  const center = { x: 500, y: 305 };
  const speciesPositions = new Map();
  const reactionPositions = new Map();
  species.forEach((item, index) => {
    const angle = -Math.PI / 2 + (index / species.length) * Math.PI * 2;
    speciesPositions.set(item.id, {
      x: center.x + Math.cos(angle) * 255,
      y: center.y + Math.sin(angle) * 235,
    });
  });
  reactions.forEach((item, index) => {
    const angle = -Math.PI / 2 + (index / reactions.length) * Math.PI * 2;
    reactionPositions.set(item.id, {
      x: center.x + Math.cos(angle) * 120,
      y: center.y + Math.sin(angle) * 105,
    });
  });

  const edges = [];
  reactions.forEach((reaction) => {
    const rp = reactionPositions.get(reaction.id);
    reaction.reactants.forEach((item) => {
      if (!visibleSpecies.has(item.species_key)) return;
      const sp = speciesPositions.get(item.species_key);
      edges.push(`<line class="network-edge" x1="${sp.x}" y1="${sp.y}" x2="${rp.x}" y2="${rp.y}"><title>${escapeHtml(labelFor(item))} to ${reaction.id}</title></line>`);
    });
    reaction.products.forEach((item) => {
      if (!visibleSpecies.has(item.species_key)) return;
      const sp = speciesPositions.get(item.species_key);
      edges.push(`<line class="network-edge product" x1="${rp.x}" y1="${rp.y}" x2="${sp.x}" y2="${sp.y}"><title>${reaction.id} to ${escapeHtml(labelFor(item))}</title></line>`);
    });
  });

  const speciesNodes = species.map((item) => {
    const p = speciesPositions.get(item.id);
    const label = truncate(item.label || item.formula || item.id, 18);
    return `<g><circle class="species-node" cx="${p.x}" cy="${p.y}" r="16"><title>${escapeHtml(item.label || item.formula)} | observations ${item.total_observations}</title></circle><text class="node-label" x="${p.x}" y="${p.y + 31}">${escapeHtml(label)}</text></g>`;
  }).join("");
  const reactionNodes = reactions.map((item) => {
    const p = reactionPositions.get(item.id);
    return `<g transform="translate(${p.x} ${p.y})"><rect class="reaction-node" x="-13" y="-13" width="26" height="26" transform="rotate(45)"><title>${item.id} | ${item.event_count} events</title></rect><text class="reaction-label" x="0" y="3">${item.id.slice(1)}</text></g>`;
  }).join("");

  canvas.innerHTML = `<svg viewBox="0 0 1000 620" role="img" aria-label="反应网络图"><defs><marker id="arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7 z" fill="#638d85"></path></marker></defs>${edges.join("")}${speciesNodes}${reactionNodes}</svg>`;
  renderReactionList(reactions);
}

function renderReactionList(reactions) {
  $("#reaction-list").innerHTML = reactions.slice(0, 10).map((reaction) => {
    const left = reaction.reactants.map(labelFor).join(" + ");
    const right = reaction.products.map(labelFor).join(" + ");
    return `<div class="reaction-row"><strong>${escapeHtml(reaction.id)}</strong><code title="${escapeHtml(`${left} -> ${right}`)}">${escapeHtml(`${left} -> ${right}`)}</code><span>${reaction.forward_count} / ${reaction.reverse_count}</span></div>`;
  }).join("");
}

function renderSpecies(payload) {
  const rows = payload.rows || [];
  $("#species-count-label").textContent = payload.has_more ? `显示前 ${rows.length} 行，完整目录可在输出文件中下载` : `${rows.length} 行`;
  renderTable($("#species-table"), rows, [
    ["structure_id", "结构编号"], ["formula", "分子式"], ["smiles", "SMILES"], ["inchi", "InChI"],
    ["inchikey", "InChIKey"], ["atom_count", "原子数"], ["bond_count", "键数"],
    ["status", "状态"],
  ]);
}

function renderTable(table, rows, columns, statusKey = null) {
  const head = columns.map(([, label]) => `<th>${escapeHtml(label)}</th>`).join("");
  const body = rows.map((row) => {
    const cells = columns.map(([key]) => {
      const value = row[key] ?? "";
      let className = "";
      let display = value;
      if (statusKey && key === statusKey) {
        className = value === "true" ? "accepted" : "filtered";
        display = value === "true" ? "接受" : "过滤";
      }
      return `<td class="${className}" title="${escapeHtml(value)}">${escapeHtml(display)}</td>`;
    }).join("");
    return `<tr>${cells}</tr>`;
  }).join("");
  table.innerHTML = `<thead><tr>${head}</tr></thead><tbody>${body || `<tr><td colspan="${columns.length}">无数据</td></tr>`}</tbody>`;
}

function renderOutputs(files) {
  const descriptions = {
    species_timeseries: "逐帧物种计数与分子标识",
    reaction_events: "断键、成键、方向及过滤状态",
    molecule_registry: "去重后的物种身份注册表",
    reaction_network: "物种节点与聚合反应超边",
    event_index: "完整事件证据与分页索引（SQLite）",
    run_record: "任务参数、处理统计和运行状态（JSON）",
    bond_counts: "逐帧元素对连接现存数量",
    bond_changes: "相邻帧各类连接新增数与消失数",
    functional_groups: "逐帧官能团拓扑候选个数",
    functional_group_molecules: "逐帧含各类官能团候选的分子数",
    count_metadata: "数量统计范围、来源和识别规则",
  };
  const downloadable = [...files, { name: "run_record", url: `/api/jobs/${encodeURIComponent(state.jobId)}` }];
  $("#output-files").innerHTML = downloadable.map((file) => (
    `<div class="output-file"><div><strong>${escapeHtml(file.name)}</strong><span>${escapeHtml(descriptions[file.name] || "分析输出")}</span></div><a href="${escapeHtml(file.url)}"${file.name === "run_record" ? ' download="run_record.json"' : ""}>下载</a></div>`
  )).join("");
}

function showJobShell(filename) {
  state.resultsReady = false;
  state.speciesLoaded = false;
  state.network = null;
  window.EventBrowser.reset();
  window.CountCharts.reset();
  $("#empty-state").classList.add("hidden");
  $("#job-view").classList.remove("hidden");
  $("#results").classList.add("hidden");
  $("#job-filename").textContent = filename || "bonds trajectory";
  $("#job-short-id").textContent = "PENDING";
  $("#job-status").textContent = "准备";
  $("#job-status").className = "status-badge";
}

function updateProgress(percent, stage, message, volume) {
  const bounded = Math.max(0, Math.min(100, Number(percent) || 0));
  $("#progress-stage").textContent = stage || "处理中";
  $("#progress-percent").textContent = `${Math.round(bounded)}%`;
  $("#progress-bar").style.width = `${bounded}%`;
  $(".progress-track").setAttribute("aria-valuenow", String(Math.round(bounded)));
  $("#progress-message").textContent = message || "";
  $("#progress-volume").textContent = volume || "";
}

function showError(message) {
  const target = $("#form-error");
  target.textContent = message;
  target.classList.remove("hidden");
}

function hideError() {
  $("#form-error").classList.add("hidden");
  $("#form-error").textContent = "";
}

async function api(url, options = {}) {
  const response = await fetch(url, options);
  if (!response.ok) {
    let detail = `${response.status} ${response.statusText}`;
    try {
      const payload = await response.json();
      detail = payload.detail || detail;
    } catch (error) {
      // Keep the HTTP status when the response is not JSON.
    }
    throw new Error(detail);
  }
  if (response.status === 204) return null;
  return response.json();
}

function selectedFilename() {
  return $("#bonds-file").files[0]?.name || "bonds trajectory";
}

function basename(path) {
  return path.trim().replace(/[\\/]+$/, "").split(/[\\/]/).pop() || path;
}

function labelFor(item) {
  const coefficient = Number(item.stoichiometry || 1);
  const label = item.smiles || item.formula || item.species_key;
  return coefficient > 1 ? `${coefficient} ${label}` : label;
}

function formatBytes(value) {
  const bytes = Number(value) || 0;
  if (bytes < 1024) return `${bytes} B`;
  const units = ["KiB", "MiB", "GiB", "TiB"];
  let number = bytes / 1024;
  let unit = units[0];
  for (let index = 1; index < units.length && number >= 1024; index += 1) {
    number /= 1024;
    unit = units[index];
  }
  return `${number >= 10 ? number.toFixed(1) : number.toFixed(2)} ${unit}`;
}

function formatNumber(value) {
  const number = Number(value);
  return Number.isFinite(number) ? number.toLocaleString("zh-CN") : "-";
}

function truncate(value, length) {
  const text = String(value || "");
  return text.length > length ? `${text.slice(0, length - 3)}...` : text;
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}
