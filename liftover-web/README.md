# 离线基因组区间迁移 Web 应用（BED Liftover）

供实验室**内网/离线环境**使用的基因组坐标区间迁移工具。用户在浏览器中
上传 BED 文件与对应的 UCSC `chain` 文件，指定 BED 区间属于链的
**target 侧**还是**query 侧**，先在页面上逐条预览转换状态，确认后再
下载 TSV 结果。

- 后端：Python 3.11+ / FastAPI / Uvicorn
- 前端：原生 HTML + CSS + JavaScript（无任何 CDN / 外网依赖）
- 运行时数据不落盘：上传内容与转换结果只保存在进程内存中；校验失败时
  **不会**生成任何可下载文件

---

## 目录结构

```
liftover-web/
├── app/                  # 后端
│   ├── bed.py            # BED 解析（零基左闭右开，错误带行号）
│   ├── chain.py          # UCSC chain 解析（严格校验，错误带行号）
│   ├── engine.py         # 区间转换引擎（缺口拆分 / 歧义 / 负链）
│   ├── report.py         # TSV 导出与预览 JSON
│   └── main.py           # FastAPI 入口
├── static/               # 前端（纯本地资源）
├── examples/             # 可直接试用的 BED / chain 示例
├── tests/                # pytest 测试（26 个用例）
├── scripts/init.sh       # 本地开发依赖初始化
├── Dockerfile
├── docker-compose.yml
├── .env.example
├── requirements.txt      # 运行依赖
└── requirements-dev.txt  # 运行 + 测试依赖
```

---

## 一、Docker / Compose 启动（推荐，完全离线）

镜像在构建时安装全部 Python 依赖，运行阶段不需要访问外网。

```bash
cd liftover-web

# 构建并后台启动
docker compose up -d --build

# 查看状态 / 日志
docker compose ps
docker compose logs -f
```

访问地址：

- **http://localhost:8080** （应用页面）
- 健康检查：`http://localhost:8080/healthz`

停止：

```bash
docker compose down
```

### 端口与配置

复制环境变量模板并按需修改：

```bash
cp .env.example .env
```

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `LIFTOVER_PORT` | `8080` | 宿主机映射端口（容器内固定监听 8080） |
| `LIFTOVER_MAX_UPLOAD_MB` | `50` | 单个上传文件大小上限（MB） |
| `LIFTOVER_JOB_TTL_SECONDS` | `3600` | 转换结果在内存中的保留时长（秒） |
| `LIFTOVER_MAX_JOBS` | `100` | 内存中同时保留的结果任务数量上限 |

改端口示例：`LIFTOVER_PORT=9000 docker compose up -d`，随后访问
`http://localhost:9000`。

不用 Compose 时：

```bash
docker build -t liftover-web:1.0.0 .
docker run -d --name liftover-web -p 8080:8080 \
  -e LIFTOVER_MAX_UPLOAD_MB=100 \
  liftover-web:1.0.0
```

---

## 二、本地 Python 直接启动（开发 / 调试）

需要 Python 3.11+。

```bash
cd liftover-web

# 1) 初始化：创建虚拟环境并安装依赖（含测试依赖）
bash scripts/init.sh

# 2) 启动
.venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8080
```

访问 **http://localhost:8080**。

运行测试：

```bash
.venv/bin/python -m pytest -q
```

---

## 三、使用方法

1. 选择 BED 文件（至少前 3 列 `chrom start end`，坐标为**零基左闭右开**，
   支持第 4 列起的额外列但不参与转换；允许 `#`、`browser`、`track` 行）。
2. 选择对应的 UCSC `chain` 文件。
3. 选择 BED 区间位于 chain 的哪一侧：
   - **target 侧**：target → query；
   - **query 侧**：query → target。
4. 点击「校验并转换」。页面显示每条源区间的状态；确认无误后点击
   「下载完整结果 TSV」。

### 状态定义

| 状态 | 含义 |
| --- | --- |
| 完全映射 `mapped` | 区间全部落在对齐块内 |
| 部分映射 `partial` | 一部分落在链内缺口/未对齐区域；命中部分被拆成**按源坐标有序的片段**，未映射部分丢弃 |
| 歧义 `ambiguous` | 同一源区间命中**多个 chain / 多个目标位置**；系统列出全部候选，**绝不自动挑选** |
| 无映射 `unmapped` | 没有任何碱基可映射；结果中保留源行号、源区间与原因 |

### 坐标与负链

- 输入与输出一律是**零基、左闭右开**（`start < end`）。
- chain 负链侧的存储坐标按 UCSC 语义（以染色体长度为枢轴）归一化后
  再做区间求交，因此**负链转换后的坐标仍然是 start < end 的零基左闭
  右开区间**；结果的 `dest_strand` 列标明该候选在目标坐标系中的方向。

### 下载文件格式（TSV）

UTF-8（带 BOM，Excel 友好），列：

```
source_line, status, source_chrom, source_start, source_end,
candidate_index, chain_id, dest_chrom, dest_strand,
fragment_index, src_start, src_end, dest_start, dest_end, reason
```

- `source_line` 是源文件**物理行号**（从 1 开始）；
- 跨缺口区间每个片段占一行；歧义区间的每个候选、每个片段各占一行；
- `unmapped` 行保留源行号与原因，候选/片段列写 `.`。

### 失败即阻止导出

- BED 坏行（列数不足、坐标非整数、`start >= end` 等）或 chain 格式错误
  （头列数/字段非法、块越界、缺少终止行、gap 越界等）都会在页面上
  **带行号列出全部问题**，接口返回 400，且**不会生成或保留任何可下载
  文件**；修正后需重新提交。

---

## 四、快速试用

```bash
# 页面上分别选择：
#   BED:   examples/example_target.bed
#   chain: examples/demo.chain
#   侧别:  target 侧
```

该示例覆盖全部四种状态（完全映射 / 跨缺口拆分的部分映射 / 歧义 /
无映射）。`examples/example_query.bed` 演示 query → target，
包含负链 query 的反向转换；`examples/example_bad.bed` 演示错误定位
与导出拦截。

---

## 五、离线部署说明

- 前端资源全部由本服务在 `/static/` 下提供，页面不含任何外部字体、
  脚本或样式表；
- Docker 构建需在有网络的机器上 `docker build`（或导出镜像），之后把
  镜像 `docker save | gzip > liftover-web.tgz` 拷贝到内网，
  `docker load` 后即可运行；
- 容器以非 root 用户运行，只暴露 8080 端口，不挂载任何数据卷；
- 上传文件与结果仅存内存，进程重启即清除，适合处理受限数据。
