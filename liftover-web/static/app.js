"use strict";

const PREVIEW_LIMIT = 1000;

const STATUS_TEXT = {
  mapped: "完全映射",
  partial: "部分映射",
  ambiguous: "歧义（多位置）",
  unmapped: "无映射",
};

const form = document.getElementById("lift-form");
const bedInput = document.getElementById("bed");
const chainInput = document.getElementById("chain");
const bedHint = document.getElementById("bed-hint");
const chainHint = document.getElementById("chain-hint");
const submitBtn = document.getElementById("submit-btn");
const errorBox = document.getElementById("error-box");
const errorList = document.getElementById("error-list");
const resultBox = document.getElementById("result-box");
const statsEl = document.getElementById("stats");
const tbody = document.getElementById("result-tbody");
const previewNote = document.getElementById("preview-note");
const downloadBtn = document.getElementById("download-btn");
const filterBtns = document.querySelectorAll(".filter");

let currentToken = null;
let currentName = null;
let currentRows = [];
let activeFilter = "all";

bedInput.addEventListener("change", () => {
  bedHint.textContent = bedInput.files[0]
    ? `${bedInput.files[0].name}（${formatSize(bedInput.files[0].size)}）`
    : "未选择文件";
});
chainInput.addEventListener("change", () => {
  chainHint.textContent = chainInput.files[0]
    ? `${chainInput.files[0].name}（${formatSize(chainInput.files[0].size)}）`
    : "未选择文件";
});

filterBtns.forEach((btn) => {
  btn.addEventListener("click", () => {
    activeFilter = btn.dataset.f;
    filterBtns.forEach((b) => b.classList.toggle("active", b === btn));
    renderRows();
  });
});

form.addEventListener("submit", async (ev) => {
  ev.preventDefault();
  errorBox.classList.add("hidden");
  resultBox.classList.add("hidden");
  submitBtn.disabled = true;
  submitBtn.textContent = "转换中…";

  const fd = new FormData();
  fd.append("bed", bedInput.files[0]);
  fd.append("chain", chainInput.files[0]);
  fd.append("side", form.querySelector('input[name="side"]:checked').value);

  try {
    const resp = await fetch("/api/liftover", { method: "POST", body: fd });
    const data = await resp.json();
    if (!resp.ok) {
      showErrors(data.errors || [`请求失败（HTTP ${resp.status}）`]);
      return;
    }
    currentToken = data.download_token;
    currentName = data.download_name;
    currentRows = data.preview.rows;
    renderStats(data.preview);
    renderRows();
    resultBox.classList.remove("hidden");
    resultBox.scrollIntoView({ behavior: "smooth", block: "start" });
  } catch (err) {
    showErrors([`网络或服务错误：${err && err.message ? err.message : err}`]);
  } finally {
    submitBtn.disabled = false;
    submitBtn.textContent = "校验并转换";
  }
});

downloadBtn.addEventListener("click", () => {
  if (!currentToken) return;
  // 直接触发浏览器下载；token 只在转换完全成功后才存在
  const a = document.createElement("a");
  a.href = `/api/download/${encodeURIComponent(currentToken)}`;
  a.download = currentName || "";
  document.body.appendChild(a);
  a.click();
  a.remove();
});

function showErrors(errors) {
  errorList.innerHTML = "";
  for (const msg of errors) {
    const li = document.createElement("li");
    li.textContent = msg;
    errorList.appendChild(li);
  }
  errorBox.classList.remove("hidden");
  errorBox.scrollIntoView({ behavior: "smooth", block: "start" });
}

function renderStats(preview) {
  const c = preview.counts;
  statsEl.innerHTML = "";
  const items = [
    ["总数", preview.total, ""],
    ["完全映射", c.mapped, "mapped"],
    ["部分映射", c.partial, "partial"],
    ["歧义", c.ambiguous, "ambiguous"],
    ["无映射", c.unmapped, "unmapped"],
  ];
  for (const [label, n, cls] of items) {
    const div = document.createElement("div");
    div.className = `stat ${cls}`;
    div.innerHTML = `<b>${n}</b>${label}`;
    statsEl.appendChild(div);
  }
}

function renderRows() {
  tbody.innerHTML = "";
  const rows = currentRows.filter((r) => activeFilter === "all" || r.status === activeFilter);
  const shown = rows.slice(0, PREVIEW_LIMIT);
  previewNote.textContent =
    currentRows.length > PREVIEW_LIMIT
      ? `预览仅显示前 ${PREVIEW_LIMIT} 行（当前筛选 ${rows.length} 行 / 共 ${currentRows.length} 行）；下载的 TSV 包含全部结果。`
      : `共 ${currentRows.length} 行源区间；多片段、歧义候选在“候选 / 片段”列内分行展开。`;

  for (const row of shown) {
    const tr = document.createElement("tr");

    const tdLine = document.createElement("td");
    tdLine.className = "mono";
    tdLine.textContent = row.line_no;

    const tdSrc = document.createElement("td");
    tdSrc.className = "mono";
    tdSrc.textContent = `${row.src_chrom}:${row.src_start}-${row.src_end}`;

    const tdStatus = document.createElement("td");
    const badge = document.createElement("span");
    badge.className = `badge ${row.status}`;
    badge.textContent = STATUS_TEXT[row.status] || row.status;
    tdStatus.appendChild(badge);

    const tdDest = document.createElement("td");
    if (row.candidates.length === 0) {
      tdDest.textContent = "—";
    } else {
      row.candidates.forEach((cand, ci) => {
        const div = document.createElement("div");
        div.className = "cand";
        if (row.candidates.length > 1) {
          const head = document.createElement("div");
          head.className = "cand-head";
          head.textContent =
            `候选 ${ci + 1}/${row.candidates.length} · chain ${cand.chain_id}`;
          div.appendChild(head);
        }
        cand.fragments.forEach((frag, fi) => {
          const fl = document.createElement("div");
          fl.className = "frag-line";
          const prefix = cand.fragments.length > 1
            ? `片段${fi + 1}(${frag.src_start}-${frag.src_end}) → `
            : "";
          fl.textContent =
            `${prefix}${cand.dest_chrom}:${frag.dest_start}-${frag.dest_end} ` +
            `(${cand.dest_strand})`;
          div.appendChild(fl);
        });
        tdDest.appendChild(div);
      });
    }

    const tdReason = document.createElement("td");
    tdReason.className = "reason";
    tdReason.textContent = row.reason || "";

    tr.append(tdLine, tdSrc, tdStatus, tdDest, tdReason);
    tbody.appendChild(tr);
  }
}

function formatSize(n) {
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(1)} KB`;
  return `${(n / (1024 * 1024)).toFixed(2)} MB`;
}
