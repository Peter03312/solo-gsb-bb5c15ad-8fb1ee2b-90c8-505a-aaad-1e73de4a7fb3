#!/usr/bin/env bash
# 本地（非 Docker）依赖初始化：创建虚拟环境并安装 requirements。
set -euo pipefail

cd "$(dirname "$0")"

PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-.venv}"

"$PYTHON_BIN" -m venv "$VENV_DIR"
# shellcheck disable=SC1090
source "$VENV_DIR/bin/activate"
pip install --upgrade pip
pip install -r requirements.txt

cat <<EOF

依赖初始化完成。启动服务：

  source $VENV_DIR/bin/activate
  uvicorn app.main:app --host 0.0.0.0 --port 8000

然后访问 http://localhost:8000
EOF
