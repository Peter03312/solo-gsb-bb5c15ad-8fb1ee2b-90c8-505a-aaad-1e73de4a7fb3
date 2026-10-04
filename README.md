# 离线基因组区间迁移（BED LiftOver Web）

供实验室**内网离线**使用的基因组区间坐标迁移小工具：上传 BED 文件与对应的
[UCSC chain](https://genome.ucsc.edu/goldenPath/help/chain.html) 文件，指定 BED
位于 chain 的 **target** 还是 **query** 侧，在浏览器中预览每条源区间的转换状态，
然后一键下载 TSV 结果。

- 后端：Python 3.11 + FastAPI（文件全程只在内存中处理，服务端不落盘）
- 前端：单页原生 HTML/JS（无外部 CDN，离线可用）
- 交付：Dockerfile + docker-compose + 本地依赖初始化脚本

## 坐标与语义约定

1. 输入 BED 与 chain 一律按 **零基、左闭右开 `[start, end)`** 解释；
   输出结果同样保持零基左闭右开。**负链**（query strand = `-`）转换后坐标已在
   算法内部完成反转，结果坐标仍是基因组正向的零基左闭右开区间。
2. 区间跨 chain 缺口时，拆成若干**按源坐标有序**的片段输出，并标记
   `partial`（部分映射），备注命中 bp 数；目标片段列带 `frag_index/frag_count`。
3. 一个源区间若在**多条 chain** 上都有命中（多个目标位置），标记 `ambiguous`，
   列出每个候选，备注“未自动挑选”——程序不会擅自替用户选择位置。
4. 完全无映射的区间保留源行（行号、染色体、起始终止、方向），目标列留空，
   状态为 `unmapped`，备注原因（染色体不存在 / 全部落在缺口上）。
5. 方向：chain query 正链时输出方向与 BED 第 6 列一致；query 负链时翻转；
   BED 未写 strand 默认 `+`；strand 为 `.` 时输出也为 `.`。

> chain 文件 target strand 必须为 `+`（UCSC 标准）；query 支持 `+`/`-`。

## 状态含义

| 状态 | 含义 |
| --- | --- |
| `mapped` | 全部 bp 命中（可能因对齐块边界拆成多个首尾相接的有序片段） |
| `partial` | 仅部分 bp 命中，其余落在 chain 缺口，片段按源坐标有序 |
| `ambiguous` | 在多条 chain 上命中，存在多个目标位置，候选全部列出，不自动挑选 |
| `unmapped` | 完全无映射，保留源行与原因，目标列留空 |

## 快速开始（Docker，推荐离线部署）

在项目根目录执行：

```bash
docker compose up -d --build
```

然后浏览器访问：

```
http://localhost:8000
```

停止：

```bash
docker compose down
```

### 配置

通过环境变量或项目根目录的 `.env` 配置（compose 会自动读取）：

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `HOST_PORT` | `8000` | 宿主机映射端口，被占用时可改，例如 `HOST_PORT=18000` |
| `PREVIEW_LIMIT` | `1000` | 网页预览表格最多返回的行数；下载结果始终包含全部行 |

示例：

```bash
HOST_PORT=18000 PREVIEW_LIMIT=5000 docker compose up -d --build
# => http://localhost:18000
```

上传文件大小上限 50 MiB（在 `app/liftover.py` / `app/main.py` 中调整
`MAX_BYTES` 后重建镜像）。

## 本地运行（不使用 Docker）

需要 Python 3.11+。初始化依赖（创建 `.venv` 并安装 `requirements.txt`）：

```bash
./scripts/init.sh
```

启动：

```bash
source .venv/bin/activate
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

访问：

```
http://localhost:8000
```

## 使用方法

1. 选择 BED 文件：至少包含 `chrom start end` 三列；可选第 4 列 name、
   第 6 列 strand（`+`/`-`/`.`）。以 `#`、`track`、`browser` 开头的行会跳过，
   空行跳过，**结果里的源行号是物理行号**。
2. 选择对应的 UCSC `.chain` 文件。
3. 选择 BED 属于 chain 的 target 侧（转到 query）还是 query 侧（转到 target）。
4. 点「转换并预览」：页面按状态着色显示每条源区间/片段的结果与统计。
5. 点「下载结果（TSV）」保存 `liftover_result.tsv`（UTF-8，含 BOM，Excel 可直接打开）。

### 输入校验与“失败不留下文件”

- BED 坏行（字段不足、坐标非整数/为负/`start >= end`、非法 strand 等）或
  chain 格式错误（头部字段数不对、字段非法、孤立数据行、块/gap 累计与头部不符、
  对齐块越界等）会在页面顶部**逐条列出文件与行号**。
- 只要存在任何这类错误，转换与导出都会被拒绝（HTTP 422），
  浏览器收到的是错误信息而不是文件；服务端全程内存处理，**失败不会产生任何可下载文件**。

### 输出 TSV 列

`source_line, source_chrom, source_start, source_end, source_strand,
status, frag_index, frag_count, dest_chrom, dest_start, dest_end,
dest_strand, note`

- 一个拆分区间占多行，`frag_index/frag_count` 表示“第几段/共几段”；
  `unmapped` 行目标列留空；`ambiguous` 行列出每个候选位置。

## HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET | `/` | 前端页面 |
| GET | `/api/health` | 健康检查，返回 `{"status":"ok"}` |
| POST | `/api/convert` | multipart：`bed`、`chain`、`side`，返回预览 JSON（最多 `PREVIEW_LIMIT` 行） |
| POST | `/api/export` | 同上，成功返回 `liftover_result.tsv`；输入错误返回 422 JSON，不产生文件 |

## 测试

```bash
python3 tests/test_liftover.py
# 或
pytest tests/
```

覆盖：正/负链块构造、零基左闭右开边界、负链坐标反转与方向翻转、
缺口拆分与部分映射、多链歧义不挑选、无映射保留源行与原因、
BED/chain 坏行行号报错并阻止导出、TSV 结构等。

## 目录结构

```
app/
  liftover.py        # chain/BED 解析与迁移引擎（纯标准库）
  main.py            # FastAPI 路由
  static/index.html  # 单页前端
scripts/init.sh      # 本地虚拟环境依赖初始化
tests/test_liftover.py
Dockerfile
docker-compose.yml
requirements.txt
```
